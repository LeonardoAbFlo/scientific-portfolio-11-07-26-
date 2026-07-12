#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# =========================
# === CONFIG ==============
# =========================
PROJECT_ROOT="/path/to/mosquitos/F3"                   # <- ajusta si corresponde
WORKDIR="$PROJECT_ROOT/mitogenomas"                         # FASTA por barcode (p.ej., barcode89.fasta)
INPUT_DIR="$PROJECT_ROOT/dnaapler_input_trnI"               # multi-FASTA input (registro)
OUTDIR="$PROJECT_ROOT/dnaapler_trnI"                        # salida (rotados y/o dnaapler)
MITOF_DIR="$PROJECT_ROOT/mitofinder_result_v2_dnaapler"     # raíz con subcarpetas barcodeXX/...
DB="/path/to/mosquitos/db/trnI_db/fastas/trnI_nr.fasta"  # DB custom (AA=usa dnaapler; NT=usa GFF)
THREADS=8
PREFIX="culicidae_trnI"

# Si true, crea symlink *_trnI_link.scafSeq
MAKE_TRNI_SYMLINKS=true

# Modo depuración: imprime coincidencias GFF
DEBUG_GFF=false

# ===== PARCHE: forzar orientación uniforme
# Si trnI está en strand '-', RC del genoma completo -> trnI queda en '+'
FORCE_TRNI_PLUS=true

# ===== PARCHE: evitar que el pipeline use FASTA viejo
# Exporta una copia "oficial" por barcode (recomendado para pasos posteriores)
EXPORT_REORIENTED_DIR="$PROJECT_ROOT/mitogenomas_trnI_plus"
# Sobrescribe el symlink estándar (muchos pipelines usan este nombre)
OVERRIDE_STANDARD_LINK=true
# (CUIDADO) Si true, pisa WORKDIR/barcodeXX.fasta con symlink al reorientado
OVERRIDE_WORKDIR_FASTA=false

# =========================
# === PREP ===============
# =========================
mkdir -p "$OUTDIR" "$INPUT_DIR" "$(dirname "$DB")" "$EXPORT_REORIENTED_DIR"

OUT_MFA="$INPUT_DIR/all_mitos.fasta"
COMBINED_TRNI="$OUTDIR/combined_trnI_reoriented.fasta"
SUMMARY_TSV="$OUTDIR/${PREFIX}_trnI_reorientation_summary.tsv"

# =========================
# === FUNCIONES ===========
# =========================

# Devuelve 0 si la DB parece AMINOACÍDICA; 1 si parece NUCLEOTÍDICA
is_db_aa () {
  local db="$1"
  awk '
    BEGIN{bad=0}
    /^>/ {next}
    {
      s=toupper($0); gsub(/[ \t\r\n-]/,"",s)
      # Permitimos IUPAC nucleotídico. Si aparece otra letra -> asumimos AA
      if (s ~ /[^ACGTURYSWKMBDHVN]/) { bad=1; exit }
    }
    END{ exit(bad?0:1) }
  ' "$db"
}

# Longitud de un FASTA (concatenando líneas de secuencia)
fasta_len () {
  local fasta="$1"
  awk '/^>/{next} {gsub(/\r/,""); if(length($0)>0) s+=length($0)} END{print (s?s:0)}' "$fasta"
}

# Rotar un FASTA (un contig) para que inicie en POS (1-based), conservando orientación
rotate_fasta_at_pos () {
  local fasta="$1" pos="$2"
  awk -v POS="$pos" '
    BEGIN{hdr=""; seq=""}
    /^>/ { if(!hdr) hdr=$0; next }
    { gsub(/\r/,""); if(length($0)>0) seq=seq $0 }
    END{
      if(hdr==""){ exit 2 }
      L=length(seq); if(L==0){ exit 3 }
      p=POS; if(p<1 || p>L) p=1
      left=substr(seq,1,p-1); right=substr(seq,p)
      print hdr; print right left
    }
  ' "$fasta"
}

# Reverse-complement FASTA (IUPAC básico). Mantiene el header original.
revcomp_fasta () {
  local fasta="$1"
  awk '
    BEGIN{hdr=""; seq=""}
    /^>/ { if(!hdr) hdr=$0; next }
    { gsub(/\r/,""); if(length($0)>0) seq=seq toupper($0) }
    END{
      if(hdr==""){ exit 2 }
      n=length(seq); if(n==0){ exit 3 }
      out=""
      for(i=n;i>=1;i--){
        c=substr(seq,i,1)
        if(c=="A") c="T"; else if(c=="T") c="A";
        else if(c=="C") c="G"; else if(c=="G") c="C";
        else if(c=="U") c="A";
        else if(c=="R") c="Y"; else if(c=="Y") c="R";
        else if(c=="S") c="S"; else if(c=="W") c="W";
        else if(c=="K") c="M"; else if(c=="M") c="K";
        else if(c=="B") c="V"; else if(c=="V") c="B";
        else if(c=="D") c="H"; else if(c=="H") c="D";
        else if(c=="N") c="N";
        out=out c
      }
      print hdr
      print out
    }
  ' "$fasta"
}

# Buscar trnI en un GFF (robusto):
# Imprime: start<TAB>end<TAB>strand<TAB>seqid; exit 0 si encontró; 1 si no
find_trnI_from_gff () {
  local gff="$1"
  awk -F'\t' -v DBG="$DEBUG_GFF" '
    BEGIN{ IGNORECASE=1; found=0 }
    $0 !~ /^\s*#/ && NF>=9 {
      type=$3
      attrs=$9
      al=tolower(attrs)

      if (type !~ /(tRNA|gene)/) next
      gsub(/%2D/,"-",al)

      if (al ~ /trni/ || al ~ /trna[-_ ]?ile/ || (al ~ /isoleuc/ && al ~ /trna/)) {
        if (DBG=="true") print "[debug] hit GFF: " $0 > "/dev/stderr"
        print $4 "\t" $5 "\t" $7 "\t" $1
        found=1
        exit 0
      }
    }
    END{ if(!found) exit 1 }
  ' "$gff"
}

# =========================
# === 1) Construir MFA ====
# =========================
echo "Reconstruyendo multi-FASTA de entrada en: $OUT_MFA"
: > "$OUT_MFA"

FOUND=0
dirs=( "$WORKDIR"/barcode*/ )

if (( ${#dirs[@]} > 0 )); then
  echo "[modo directorios] Trabajando en $WORKDIR (carpetas barcode*/)"
  for D in "${dirs[@]}"; do
    [[ -d "$D" ]] || continue
    SAMPLE="$(basename "$D")"
    FILE=""
    for CAND in "$D"/mitogenome.fasta "$D"/*.mito*.fasta "$D"/assembly.fasta "$D"/*.fasta "$D"/*.fa; do
      [[ -s "$CAND" ]] && { FILE="$CAND"; break; }
    done
    if [[ -z "${FILE:-}" ]]; then
      echo "[skip] $SAMPLE: no se encontró FASTA"; continue
    fi
    printf ">%s\n" "$SAMPLE" >> "$OUT_MFA"
    awk '/^>/{next} {sub(/\r$/,""); if(length($0)>0) print}' "$FILE" >> "$OUT_MFA"
    printf "\n" >> "$OUT_MFA"
    ((FOUND++)) || true
    echo "[ok]  $SAMPLE  <-  $(basename "$FILE")"
  done
else
  echo "[modo archivos] Trabajando con archivos planos en $WORKDIR"
  for F in "$WORKDIR"/barcode*.fasta "$WORKDIR"/barcode*.fa "$WORKDIR"/barcode*mtDNA_contig*.fasta "$WORKDIR"/barcode*mitogenome*.fasta; do
    [[ -e "$F" ]] || continue
    base="$(basename "$F")"
    case "$base" in *genes_AA*.fasta|*genes_NT*.fasta) continue ;; esac
    SAMPLE="${base%.*}"; SAMPLE="${SAMPLE%%_*}"
    printf ">%s\n" "$SAMPLE" >> "$OUT_MFA"
    awk '/^>/{next} {sub(/\r$/,""); if(length($0)>0) print}' "$F" >> "$OUT_MFA"
    printf "\n" >> "$OUT_MFA"
    ((FOUND++)) || true
    echo "[ok]  $SAMPLE  <-  $(basename "$F")"
  done
fi

if [[ "$FOUND" -eq 0 ]]; then
  echo "ERROR: No se agregaron entradas al multi-FASTA desde $WORKDIR"
  exit 1
fi

HDRS=$(grep -c '^>' "$OUT_MFA" || true)
if [[ "$HDRS" -eq 0 ]]; then
  echo "ERROR: El multi-FASTA $OUT_MFA no contiene cabeceras"
  exit 1
fi
ls -lh "$OUT_MFA" | awk '{print $0}'

# =========================
# === 2) RUTA SEGÚN DB ====
# =========================

if is_db_aa "$DB"; then
  # -----------------------
  # Caso A: DB AMINOACÍDICA -> dnaapler
  # -----------------------
  echo "DB detectada como AMINOACÍDICA. Ejecutando dnaapler (mode=custom)…"
  dnaapler bulk \
    -i "$OUT_MFA" \
    -o "$OUTDIR" \
    -p "$PREFIX" \
    -t "$THREADS" \
    -m custom \
    -c "$DB" \
    -f

  echo "Dnaapler finalizado. Resultados en: $OUTDIR"
  echo "  - ${OUTDIR}/${PREFIX}_reoriented.fasta"
  echo "  - ${OUTDIR}/${PREFIX}_failed_to_reorient.fasta"
  echo "  - ${OUTDIR}/${PREFIX}_bulk_reorientation_summary.tsv"

  # Split por barcode
  IN_OK="${OUTDIR}/${PREFIX}_reoriented.fasta"
  OUT_SPLIT="${OUTDIR}/barcodes_re"
  [[ -d "$OUT_SPLIT" ]] && mv "$OUT_SPLIT" "${OUT_SPLIT}.prev.$(date +%F_%H%M%S)"
  mkdir -p "$OUT_SPLIT"

  if [[ -s "$IN_OK" ]]; then
    echo "Separando reorientados por barcode en: $OUT_SPLIT"
    awk -v out="$OUT_SPLIT" '
      BEGIN{fn=""}
      /^>/{
        name=$0; sub(/^>/,"",name)
        gsub(/[ \t].*$/,"",name)
        gsub(/[^A-Za-z0-9._-]/,"_",name)
        fn=sprintf("%s/%s.fasta", out, name)
        print ">" name > fn; next
      }
      { print >> fn }
    ' "$IN_OK"
    N=$(grep -c '^>' "$IN_OK" || true)
    echo "OK: generados $N archivos en $OUT_SPLIT"
  else
    echo "ATENCIÓN: No existe archivo reorientado: $IN_OK"
  fi

else
  # -----------------------
  # Caso B: DB NUCLEOTÍDICA (trnI) -> rotación con GFF de MitoFinder
  # PARCHE: si strand '-' y FORCE_TRNI_PLUS=true:
  #   RC(genoma completo) + CUT = L - END + 1 + rotación en CUT
  #   => trnI queda en '+' y empieza en 1
  # -----------------------
  echo "DB detectada como NUCLEOTÍDICA (trnI). Omitiendo dnaapler y reorientando con GFF de MitoFinder…"

  : > "$COMBINED_TRNI"
  echo -e "barcode\tstart\tend\tstrand\tcontig_len\tgff_used\trotated_fasta\tdid_revcomp\tcut_pos_used" > "$SUMMARY_TSV"

  OUT_ROT_DIR="$OUTDIR/trnI_reoriented"
  [[ -d "$OUT_ROT_DIR" ]] && mv "$OUT_ROT_DIR" "${OUT_ROT_DIR}.prev.$(date +%F_%H%M%S)"
  mkdir -p "$OUT_ROT_DIR"

  ANY=0
  for BDIR in "$MITOF_DIR"/barcode*/ ; do
    [[ -d "$BDIR" ]] || continue
    BARCODE="$(basename "$BDIR")"

    RESDIR_A="$BDIR/${BARCODE}_MitoFinder_mitfi_Final_Results"
    RESDIR_B="$BDIR/${BARCODE}/${BARCODE}_MitoFinder_mitfi_Final_Results"
    TMPDIR_A="$BDIR/${BARCODE}_tmp"
    TMPDIR_B="$BDIR/${BARCODE}/${BARCODE}_tmp"

    RESDIR=""
    if   [[ -d "$RESDIR_A" ]]; then RESDIR="$RESDIR_A"
    elif [[ -d "$RESDIR_B" ]]; then RESDIR="$RESDIR_B"
    else
      echo "[skip] $BARCODE: no encuentro *_MitoFinder_mitfi_Final_Results"
      continue
    fi

    GFF_MAIN="$RESDIR/${BARCODE}_mtDNA_contig.gff"
    FA_MAIN="$RESDIR/${BARCODE}_mtDNA_contig.fasta"

    if [[ ! -s "$FA_MAIN" ]]; then
      echo "[skip] $BARCODE: falta FASTA principal en $RESDIR"
      continue
    fi

    # Candidatos GFF (primero final, luego raw del tmp)
    GFF_CAND=()
    [[ -s "$GFF_MAIN" ]] && GFF_CAND+=( "$GFF_MAIN" )
    [[ -d "$TMPDIR_A" && -s "$TMPDIR_A/${BARCODE}_mtDNA_contig_raw.gff" ]] && GFF_CAND+=( "$TMPDIR_A/${BARCODE}_mtDNA_contig_raw.gff" )
    [[ -d "$TMPDIR_B" && -s "$TMPDIR_B/${BARCODE}_mtDNA_contig_raw.gff" ]] && GFF_CAND+=( "$TMPDIR_B/${BARCODE}_mtDNA_contig_raw.gff" )

    if (( ${#GFF_CAND[@]} == 0 )); then
      echo "[skip] $BARCODE: no encuentro GFF (final ni raw)"
      continue
    fi

    HIT=""
    GFF_USED=""
    for G in "${GFF_CAND[@]}"; do
      if COORDS="$(find_trnI_from_gff "$G" 2>/dev/null)"; then
        HIT="$COORDS"; GFF_USED="$G"; break
      elif [[ "$DEBUG_GFF" == "true" ]]; then
        echo "[debug] no hit en: $G" >&2
      fi
    done

    if [[ -z "$HIT" ]]; then
      echo "[skip] $BARCODE: no se encontró trnI en ${GFF_CAND[*]}"
      continue
    fi

    read -r START END STRAND SEQID <<< "$HIT"
    CLEN="$(fasta_len "$FA_MAIN")"
    [[ "$CLEN" -lt 1 ]] && { echo "[skip] $BARCODE: FASTA vacío"; continue; }

    OUT_FASTA="$OUT_ROT_DIR/${BARCODE}_trnI_reoriented.fasta"
    DID_RC="false"
    CUT_POS=""

    if [[ "$FORCE_TRNI_PLUS" == "true" && "$STRAND" == "-" ]]; then
      # RC primero, luego cortar en la coordenada equivalente al 5' real
      DID_RC="true"
      CUT_POS=$(( CLEN - END + 1 ))

      tmp_rc="$(mktemp)"
      revcomp_fasta "$FA_MAIN" > "$tmp_rc"
      rotate_fasta_at_pos "$tmp_rc" "$CUT_POS" > "$OUT_FASTA"
      rm -f "$tmp_rc"
    else
      # strand '+' (o no forzar)
      CUT_POS="$START"
      rotate_fasta_at_pos "$FA_MAIN" "$CUT_POS" > "$OUT_FASTA"
    fi

    # Normaliza encabezado a >barcodeXX
    tmpf="$(mktemp)"
    awk -v b="$BARCODE" 'NR==1{print ">" b; next} {print}' "$OUT_FASTA" > "$tmpf" && mv "$tmpf" "$OUT_FASTA"

    echo "[ok] $BARCODE: trnI start=$START end=$END strand=$STRAND | did_rc=$DID_RC cut=$CUT_POS | out=$(basename "$OUT_FASTA")"
    cat "$OUT_FASTA" >> "$COMBINED_TRNI"
    echo -e "${BARCODE}\t${START}\t${END}\t${STRAND}\t${CLEN}\t${GFF_USED}\t${OUT_FASTA}\t${DID_RC}\t${CUT_POS}" >> "$SUMMARY_TSV"
    ANY=1

    # ===== PARCHE USO-CORRECTO: export + symlinks
    cp -f "$OUT_FASTA" "$EXPORT_REORIENTED_DIR/${BARCODE}.fasta"

    if [[ "${MAKE_TRNI_SYMLINKS}" == "true" ]]; then
      LINK_BASE="$BDIR"
      [[ -d "$BDIR/$BARCODE" ]] && LINK_BASE="$BDIR/$BARCODE"

      # tu link "trnI"
      ln -sf "$OUT_FASTA" "$LINK_BASE/${BARCODE}_trnI_link.scafSeq"

      # link estándar (muchos pipelines usan SOLO este)
      if [[ "${OVERRIDE_STANDARD_LINK}" == "true" ]]; then
        ln -sf "$OUT_FASTA" "$LINK_BASE/${BARCODE}_link.scafSeq"
      fi
    fi

    # opcional: pisar input del WORKDIR
    if [[ "${OVERRIDE_WORKDIR_FASTA}" == "true" ]]; then
      ln -sf "$EXPORT_REORIENTED_DIR/${BARCODE}.fasta" "$WORKDIR/${BARCODE}.fasta"
    fi
  done

  if [[ "$ANY" -eq 0 ]]; then
    echo "ERROR: No se generó ningún rotado a trnI. ¿MITOF_DIR correcto? ¿Existen GFF final o raw?"
    exit 1
  fi

  echo "Resumen TSV: $SUMMARY_TSV"
  echo "Combinado FASTA: $COMBINED_TRNI"
  echo "FASTA 'oficiales' (para que uses en el realineamiento): $EXPORT_REORIENTED_DIR/*.fasta"
fi

echo "Hecho."
