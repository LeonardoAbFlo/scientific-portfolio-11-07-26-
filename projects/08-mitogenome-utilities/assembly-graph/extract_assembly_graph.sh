#!/usr/bin/env bash

set -euo pipefail

ROOT="$HOME/mosquitos/flye_data"
OUT="$HOME/mosquitos/bandage_qc"

mkdir -p "$OUT/images" "$OUT/logs"

echo "Buscando assembly_graph.gfa en: $ROOT"

find "$ROOT" -mindepth 3 -maxdepth 3 -name "assembly_graph.gfa" | sort > "$OUT/flye_graphs.list"

N=$(wc -l < "$OUT/flye_graphs.list")
echo "Grafos encontrados: $N"

if [ "$N" -eq 0 ]; then
    echo "No se encontraron archivos assembly_graph.gfa"
    exit 1
fi

QC="$OUT/bandage_qc_notes.tsv"

echo -e "sample\tgroup\tbarcode\tgraph_file\tassembly_fasta\tassembly_info\tbandage_png\tcircular_like\tbranches_or_bubbles\tcoverage_consistent\tmitochondrial_component\tnotes\tdecision" > "$QC"

while read GFA; do
    GROUP=$(basename "$(dirname "$(dirname "$GFA")")")
    BARCODE=$(basename "$(dirname "$GFA")")
    SAMPLE="${GROUP}_${BARCODE}"
    SAMPLE_DIR=$(dirname "$GFA")

    FASTA="$SAMPLE_DIR/assembly.fasta"
    INFO="$SAMPLE_DIR/assembly_info.txt"
    PNG="$OUT/images/${SAMPLE}_bandage_graph.png"
    LOG="$OUT/logs/${SAMPLE}_bandage.log"

    echo "----------------------------------------"
    echo "Procesando muestra: $SAMPLE"
    echo "GFA: $GFA"

    if BandageNG image "$GFA" "$PNG" --width 2500 --height 2000 --names --lengths --depth > "$LOG" 2>&1; then
        echo "Imagen creada: $PNG"
    else
        echo "Falló con opciones avanzadas. Intentando comando simple..."
        if BandageNG image "$GFA" "$PNG" > "$LOG" 2>&1; then
            echo "Imagen creada con comando simple: $PNG"
        else
            echo "No se pudo crear imagen para $SAMPLE. Revisa: $LOG"
            PNG="FAILED"
        fi
    fi

    echo -e "${SAMPLE}\t${GROUP}\t${BARCODE}\t${GFA}\t${FASTA}\t${INFO}\t${PNG}\tNA\tNA\tNA\tNA\tNA\tNA" >> "$QC"

done < "$OUT/flye_graphs.list"

echo "----------------------------------------"
echo "Finalizado."
echo "Lista de grafos: $OUT/flye_graphs.list"
echo "Imágenes: $OUT/images"
echo "Logs: $OUT/logs"
echo "Tabla QC: $QC"
