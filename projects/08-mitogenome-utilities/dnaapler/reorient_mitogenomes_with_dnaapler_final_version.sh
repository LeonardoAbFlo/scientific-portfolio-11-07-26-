#!/usr/bin/env bash
set -euo pipefail

# === Config ===
PROJECT_ROOT="/path/to/Desktop/mosquitos"
WORKDIR="$PROJECT_ROOT"                    # aquí están tus archivos barcodeXX_*.fasta
OUTDIR="$PROJECT_ROOT/dnaapler"
DB="$PROJECT_ROOT/db/cox1_culicidae.faa"
THREADS=8
OUT_MFA="$PROJECT_ROOT/all_mitos.fasta"
PREFIX="culicidae_cox1"

mkdir -p "$OUTDIR" "$(dirname "$DB")"

# === Checks ===
if [[ ! -s "$DB" ]]; then
  echo "ERROR: No encuentro DB de proteínas (aa) para COX1: $DB"
  exit 1
fi

echo "Reconstruyendo multi-FASTA en: $OUT_MFA"
: > "$OUT_MFA"

shopt -s nullglob

# Detecta si hay carpetas barcode*/ (modo 1) o solo archivos planos (modo 2)
dirs=( "$WORKDIR"/barcode*/ )

FOUND=0
if (( ${#dirs[@]} > 0 )); then
  echo "[modo directorios] Encontradas carpetas barcode*/"
  for D in "${dirs[@]}"; do
    [[ -d "$D" ]] || continue
    SAMPLE="$(basename "$D")"

    FILE=""
    for CAND in \
        "$D"/mitogenome.fasta \
        "$D"/*.mito*.fasta \
        "$D"/assembly.fasta \
        "$D"/*.fasta \
        "$D"/*.fa
    do
      [[ -s "$CAND" ]] && { FILE="$CAND"; break; }
    done

    if [[ -z "${FILE:-}" ]]; then
      echo "[skip] $SAMPLE: no se encontró FASTA"
      continue
    fi

    printf ">%s\n" "$SAMPLE" >> "$OUT_MFA"
    awk 'BEGIN{len=0}
         /^>/ {next}
         {sub(/\r$/,"");}
         length($0)>0 {print; len+=length($0)}
         END{ if(len==0) exit 1 }' "$FILE" >> "$OUT_MFA"
    printf "\n" >> "$OUT_MFA"

    ((FOUND++)) || true
    echo "[ok]  $SAMPLE  <-  $(basename "$FILE")"
  done

else
  echo "[modo archivos] No hay carpetas barcode*/; trabajando con archivos planos"

  # Recorre archivos de contig mito en la raíz, excluyendo *_genes_AA / *_genes_NT
  for F in "$WORKDIR"/barcode*mtDNA_contig*.fasta "$WORKDIR"/barcode*mitogenome*.fasta; do
    [[ -e "$F" ]] || continue
    base="$(basename "$F")"
    case "$base" in
      *genes_AA*.fasta|*genes_NT*.fasta) continue ;;
    esac

    # SAMPLE = prefijo hasta el primer underscore -> 'barcode81', 'barcode83', etc.
    SAMPLE="${base%.*}"
    SAMPLE="${SAMPLE%%_*}"

    printf ">%s\n" "$SAMPLE" >> "$OUT_MFA"
    awk 'BEGIN{len=0}
         /^>/ {next}
         {sub(/\r$/,"");}
         length($0)>0 {print; len+=length($0)}
         END{ if(len==0) exit 1 }' "$F" >> "$OUT_MFA"
    printf "\n" >> "$OUT_MFA"

    ((FOUND++)) || true
    echo "[ok]  $SAMPLE  <-  $(basename "$F")"
  done
fi

if [[ "$FOUND" -eq 0 ]]; then
  echo "ERROR: No se agregaron entradas al multi-FASTA."
  echo "Sugerencias:"
  echo "  - Verifica que existan archivos como barcode*mtDNA_contig*.fasta o barcode*mitogenome*.fasta"
  echo "  - Ejemplo: ls -1 $WORKDIR/barcode*mtDNA_contig*.fasta || true"
  exit 1
fi

HDRS=$(grep -c '^>' "$OUT_MFA" || true)
if [[ "$HDRS" -eq 0 ]]; then
  echo "ERROR: El multi-FASTA no contiene cabeceras '>'."
  exit 1
fi
echo "Listo: $HDRS entradas en $OUT_MFA"

echo "Ejecutando Dnaapler (threads=$THREADS, prefix=$PREFIX, mode=custom)…"
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

# === Split: volver a separar por barcode en OUTDIR/barcodes_re ===
IN_OK="${OUTDIR}/${PREFIX}_reoriented.fasta"
OUT_SPLIT="${OUTDIR}/barcodes_re"

if [[ -d "$OUT_SPLIT" ]]; then
  mv "$OUT_SPLIT" "${OUT_SPLIT}.prev.$(date +%F_%H%M%S)"
fi
mkdir -p "$OUT_SPLIT"

if [[ -s "$IN_OK" ]]; then
  echo "Separando reorientados por barcode en: $OUT_SPLIT"
  # Divide por cada cabecera; sanea el nombre para usarlo como filename
  awk -v out="$OUT_SPLIT" '
    BEGIN{fn=""}
    /^>/{
      name=$0; sub(/^>/,"",name);
      gsub(/[ \t].*$/,"",name);                 # corta descripción tras primer espacio
      gsub(/[^A-Za-z0-9._-]/,"_",name);         # sanea a filename
      fn=sprintf("%s/%s.fasta", out, name);
      print ">" name > fn;                      # crea/trunca el archivo del sample
      next
    }
    {print >> fn}
  ' "$IN_OK"

  N=$(grep -c '^>' "$IN_OK" || true)
  echo "OK: generados $N archivos en $OUT_SPLIT"
else
  echo "ATENCIÓN: No existe archivo reorientado: $IN_OK"
fi
