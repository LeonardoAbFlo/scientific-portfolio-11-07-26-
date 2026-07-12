#!/usr/bin/env bash
set -euo pipefail

# =========================
# === CONFIG ==============
# =========================
PROJECT_ROOT="/path/to/mosquitos/F4"
WORKDIR="$PROJECT_ROOT/mitogenomas"
INPUT_DIR="$PROJECT_ROOT/dnaapler_input_cr"
OUTDIR="$PROJECT_ROOT/dnaapler_cr_before_trnI"
MITOF_DIR="$PROJECT_ROOT/mitofinder_result_v2_dnaapler"
DB="/path/to/mosquitos/db/trnI_db/fastas/trnI_nr.fasta"
THREADS=8
PREFIX="culicidae_CR_before_trnI"

# Si true, crea symlink *_CR_before_trnI_link.scafSeq
MAKE_CR_SYMLINKS=true

# Depuración
DEBUG_GFF=false

# Fuerza orientación uniforme: si trnI está en strand '-', hace RC para dejarlo en '+'
FORCE_TRNI_PLUS=true

# Exporta copia "oficial" de cada FASTA reorientado
EXPORT_REORIENTED_DIR="$PROJECT_ROOT/mitogenomas_CR_before_trnI"

# Si true, sobreescribe el symlink estándar que muchos pipelines usan: barcodeXX_link.scafSeq
OVERRIDE_STANDARD_LINK=true

# Si true, sobreescribe WORKDIR/barcodeXX.fasta con symlink al reorientado (CUIDADO)
OVERRIDE_WORKDIR_FASTA=false

# === NUEVO ===
# Si no se encuentra rrnS, corta esta cantidad de nt antes de trnI
CR_FALLBACK_BP=1000

# Rango "razonable" para la CR inferida rrnS -> trnI; si cae fuera, usa fallback
MIN_CR_BP=50
MAX_CR_BP=3000

# =========================
# === PREP ================
# =========================
mkdir -p "$OUTDIR" "$INPUT_DIR" "$(dirname "$DB")" "$EXPORT_REORIENTED_DIR"

OUT_MFA="$INPUT_DIR/all_mitos.fasta"
COMBINED_CR="$OUTDIR/combined_CR_before_trnI_reoriented.fasta"
SUMMARY_TSV="$OUTDIR/${PREFIX}_reorientation_summary.tsv"

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
      if (s ~ /[^ACGTURYSWKMBDHVN]/) { bad=1; exit }
    }
    END{ exit(bad?0:1) }
  ' "$db"
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

# Buscar trnI en un GFF
# Imprime: start<TAB>end<TAB>strand<TAB>seqid
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
        if (DBG=="true") print "[debug] trnI hit GFF: " $0 > "/dev/stderr"
        print $4 "\t" $5 "\t" $7 "\t" $1
        found=1
        exit 0
      }
    }
    END{ if(!found) exit 1 }
  ' "$gff"
}

# Buscar rrnS / 12S en un GFF
# Imprime: start<TAB>end<TAB>strand<TAB>seqid
find_rrnS_from_gff () {
  local gff="$1"
  awk -F'\t' -v DBG="$DEBUG_GFF" '
    BEGIN{ IGNORECASE=1; found=0 }
    $0 !~ /^\s*#/ && NF>=9 {
      type=$3
      attrs=$9
      al=tolower(attrs)
      if (type !~ /(rRNA|gene)/) next
      gsub(/%2D/,"-",al)
      if (
           al ~ /rrns/ ||
           al ~ /12s/ ||
           al ~ /small[_ -]?subunit/ ||
           al ~ /s[-_ ]?rrna/ ||
           al ~ /12s ribosomal rna/ ||
           al ~ /srna/
         ) {
        if (DBG=="true") print "[debug] rrnS hit GFF: " $0 > "/dev/stderr"
        print $4 "\t" $5 "\t" $7 "\t" $1
        found=1
        exit 0
      }
    }
    END{ if(!found) exit 1 }
  ' "$gff"
}

# Longitud de un FASTA
fasta_len () {
  local fasta="$1"
  awk '/^>/{next} {gsub(/\r/,""); if(length($0)>0) s+=length($0)} END{print (s?s:0)}' "$fasta"
}

# Remapea coordenadas de una feature tras hacer reverse-complement del contig
# Entrada: L start end strand
# Salida: new_start<TAB>new_end<TAB>new_strand
remap_feature_after_rc () {
  local L="$1" START="$2" END="$3" STRAND="$4"
  local NS NE NSTR
  NS=$(( L - END + 1 ))
  NE=$(( L - START + 1 ))
  NSTR="$STRAND"
  if [[ "$STRAND" == "+" ]]; then
    NSTR="-"
  elif [[ "$STRAND" == "-" ]]; then
    NSTR="+"
  fi
  printf "%s\t%s\t%s\n" "$NS" "$NE" "$NSTR"
}

# Envuelve posición circular al rango [1, L]
wrap_pos () {
  local p="$1" L="$2"
  while (( p < 1 )); do p=$(( p + L )); done
  while (( p > L )); do p=$(( p - L )); done
  echo "$p"
}

# Longitud del intervalo circular "entre" from_end y to_start (sin incluir extremos)
# Útil para CR = rrnS_end+1 ... trnI_start-1
circular_gap_len () {
  local from_end="$1" to_start="$2" L="$3"
  local d=$(( to_start - from_end - 1 ))
  while (( d < 0 )); do d=$(( d + L )); done
  echo "$d"
}

# =========================
# === 1) Construir MFA ====
# =========================
echo "Reconstruyendo multi-FASTA de entrada en: $OUT_MFA"
: > "$OUT_MFA"

shopt -s nullglob
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
    [[ -z "${FILE:-}" ]] && { echo "[skip] $SAMPLE: no se encontró FASTA"; continue; }
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

[[ "$FOUND" -eq 0 ]] && { echo "ERROR: No se agregaron entradas al multi-FASTA desde $WORKDIR"; exit 1; }
HDRS=$(grep -c '^>' "$OUT_MFA" || true)
[[ "$HDRS" -eq 0 ]] && { echo "ERROR: El multi-FASTA $OUT_MFA no contiene cabeceras"; exit 1; }
ls -lh "$OUT_MFA" | awk '{print $0}'

# =========================
# === 2) RUTA SEGÚN DB ====
# =========================

if is_db_aa "$DB"; then
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
  echo "DB detectada como NUCLEOTÍDICA. Reorientando al INICIO DE LA CR (antes de trnI) usando GFF de MitoFinder…"

  : > "$COMBINED_CR"
  echo -e "barcode\ttrnI_start\ttrnI_end\ttrnI_strand\trrnS_start\trrnS_end\trrnS_strand\tcontig_len\tgff_trnI\tgff_rrnS\trotated_fasta\tdid_revcomp\tcut_pos_used\tcut_reason\tinferred_CR_len\tused_fallback" > "$SUMMARY_TSV"

  OUT_ROT_DIR="$OUTDIR/CR_before_trnI_reoriented"
  [[ -d "$OUT_ROT_DIR" ]] && mv "$OUT_ROT_DIR" "${OUT_ROT_DIR}.prev.$(date +%F_%H%M%S)"
  mkdir -p "$OUT_ROT_DIR"

  shopt -s nullglob
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
    [[ -s "$FA_MAIN" ]] || { echo "[skip] $BARCODE: falta FASTA principal en $RESDIR"; continue; }

    GFF_CAND=()
    [[ -s "$GFF_MAIN" ]] && GFF_CAND+=( "$GFF_MAIN" )
    [[ -d "$TMPDIR_A" && -s "$TMPDIR_A/${BARCODE}_mtDNA_contig_raw.gff" ]] && GFF_CAND+=( "$TMPDIR_A/${BARCODE}_mtDNA_contig_raw.gff" )
    [[ -d "$TMPDIR_B" && -s "$TMPDIR_B/${BARCODE}_mtDNA_contig_raw.gff" ]] && GFF_CAND+=( "$TMPDIR_B/${BARCODE}_mtDNA_contig_raw.gff" )
    (( ${#GFF_CAND[@]} == 0 )) && { echo "[skip] $BARCODE: no encuentro GFF"; continue; }

    TRNI_HIT=""
    GFF_TRNI=""
    RRNS_HIT=""
    GFF_RRNS=""

    for G in "${GFF_CAND[@]}"; do
      if [[ -z "$TRNI_HIT" ]]; then
        if COORDS="$(find_trnI_from_gff "$G" 2>/dev/null)"; then
          TRNI_HIT="$COORDS"
          GFF_TRNI="$G"
        elif [[ "$DEBUG_GFF" == "true" ]]; then
          echo "[debug] no hit trnI en: $G" >&2
        fi
      fi

      if [[ -z "$RRNS_HIT" ]]; then
        if COORDS="$(find_rrnS_from_gff "$G" 2>/dev/null)"; then
          RRNS_HIT="$COORDS"
          GFF_RRNS="$G"
        elif [[ "$DEBUG_GFF" == "true" ]]; then
          echo "[debug] no hit rrnS en: $G" >&2
        fi
      fi

      [[ -n "$TRNI_HIT" && -n "$RRNS_HIT" ]] && break
    done

    [[ -n "$TRNI_HIT" ]] || { echo "[skip] $BARCODE: no se encontró trnI en ${GFF_CAND[*]}"; continue; }

    read -r TRNI_START TRNI_END TRNI_STRAND TRNI_SEQID <<< "$TRNI_HIT"

    RRNS_START="NA"
    RRNS_END="NA"
    RRNS_STRAND="NA"
    if [[ -n "$RRNS_HIT" ]]; then
      read -r RRNS_START RRNS_END RRNS_STRAND RRNS_SEQID <<< "$RRNS_HIT"
    fi

    CLEN="$(fasta_len "$FA_MAIN")"
    [[ "$CLEN" -ge 1 ]] || { echo "[skip] $BARCODE: FASTA vacío"; continue; }

    WORK_FASTA="$FA_MAIN"
    TMP_RC=""
    DID_RC="false"

    # 1) Forzar trnI en + si corresponde
    if [[ "$FORCE_TRNI_PLUS" == "true" && "$TRNI_STRAND" == "-" ]]; then
      DID_RC="true"
      TMP_RC="$(mktemp)"
      revcomp_fasta "$FA_MAIN" > "$TMP_RC"
      WORK_FASTA="$TMP_RC"

      read -r TRNI_START TRNI_END TRNI_STRAND <<< "$(remap_feature_after_rc "$CLEN" "$TRNI_START" "$TRNI_END" "$TRNI_STRAND")"

      if [[ "$RRNS_START" != "NA" ]]; then
        read -r RRNS_START RRNS_END RRNS_STRAND <<< "$(remap_feature_after_rc "$CLEN" "$RRNS_START" "$RRNS_END" "$RRNS_STRAND")"
      fi
    fi

    # 2) Elegir punto de corte:
    #    preferido = inicio de CR = rrnS_end + 1
    #    fallback  = CR_FALLBACK_BP antes de trnI
    USED_FALLBACK="false"
    CUT_REASON=""
    INFERRED_CR_LEN="NA"

    if [[ "$RRNS_START" != "NA" ]]; then
      INFERRED_CR_LEN="$(circular_gap_len "$RRNS_END" "$TRNI_START" "$CLEN")"
      if (( INFERRED_CR_LEN >= MIN_CR_BP && INFERRED_CR_LEN <= MAX_CR_BP )); then
        CUT_POS="$(wrap_pos $((RRNS_END + 1)) "$CLEN")"
        CUT_REASON="CR_start_from_rrnS_end_plus_1"
      else
        USED_FALLBACK="true"
        CUT_POS="$(wrap_pos $((TRNI_START - CR_FALLBACK_BP)) "$CLEN")"
        CUT_REASON="fallback_${CR_FALLBACK_BP}bp_before_trnI_suspicious_CRlen_${INFERRED_CR_LEN}"
      fi
    else
      USED_FALLBACK="true"
      CUT_POS="$(wrap_pos $((TRNI_START - CR_FALLBACK_BP)) "$CLEN")"
      CUT_REASON="fallback_${CR_FALLBACK_BP}bp_before_trnI_no_rrnS"
    fi

    OUT_FASTA="$OUT_ROT_DIR/${BARCODE}_CR_before_trnI_reoriented.fasta"
    rotate_fasta_at_pos "$WORK_FASTA" "$CUT_POS" > "$OUT_FASTA"
    [[ -n "$TMP_RC" ]] && rm -f "$TMP_RC"

    # Normaliza encabezado a >barcodeXX
    tmpf="$(mktemp)"
    awk -v b="$BARCODE" 'NR==1{print ">" b; next} {print}' "$OUT_FASTA" > "$tmpf" && mv "$tmpf" "$OUT_FASTA"

    echo "[ok] $BARCODE: trnI=$TRNI_START-$TRNI_END($TRNI_STRAND) rrnS=$RRNS_START-$RRNS_END($RRNS_STRAND) | did_rc=$DID_RC | cut=$CUT_POS | CR_len=$INFERRED_CR_LEN | fallback=$USED_FALLBACK"

    cat "$OUT_FASTA" >> "$COMBINED_CR"
    echo -e "${BARCODE}\t${TRNI_START}\t${TRNI_END}\t${TRNI_STRAND}\t${RRNS_START}\t${RRNS_END}\t${RRNS_STRAND}\t${CLEN}\t${GFF_TRNI}\t${GFF_RRNS}\t${OUT_FASTA}\t${DID_RC}\t${CUT_POS}\t${CUT_REASON}\t${INFERRED_CR_LEN}\t${USED_FALLBACK}" >> "$SUMMARY_TSV"
    ANY=1

    # Exportar copia "oficial" por barcode
    cp -f "$OUT_FASTA" "$EXPORT_REORIENTED_DIR/${BARCODE}.fasta"

    # Symlinks para forzar que el pipeline use el reorientado
    if [[ "${MAKE_CR_SYMLINKS}" == "true" ]]; then
      LINK_BASE="$BDIR"
      [[ -d "$BDIR/$BARCODE" ]] && LINK_BASE="$BDIR/$BARCODE"

      ln -sf "$OUT_FASTA" "$LINK_BASE/${BARCODE}_CR_before_trnI_link.scafSeq"

      if [[ "${OVERRIDE_STANDARD_LINK}" == "true" ]]; then
        ln -sf "$OUT_FASTA" "$LINK_BASE/${BARCODE}_link.scafSeq"
      fi
    fi

    if [[ "${OVERRIDE_WORKDIR_FASTA}" == "true" ]]; then
      ln -sf "$EXPORT_REORIENTED_DIR/${BARCODE}.fasta" "$WORKDIR/${BARCODE}.fasta"
    fi
  done

  [[ "$ANY" -eq 1 ]] || { echo "ERROR: No se generó ningún rotado al inicio de la CR."; exit 1; }

  echo "Resumen TSV: $SUMMARY_TSV"
  echo "Combinado FASTA: $COMBINED_CR"
  echo "Exportados (uso recomendado para siguientes pasos): $EXPORT_REORIENTED_DIR/*.fasta"
fi

echo "Hecho."
