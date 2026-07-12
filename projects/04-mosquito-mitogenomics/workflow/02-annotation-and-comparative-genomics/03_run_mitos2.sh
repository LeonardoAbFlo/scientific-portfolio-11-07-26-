#!/usr/bin/env bash
set -euo pipefail

# Run MITOS2 on rrnS-corrected mitochondrial FASTA files.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd -P)}"
ANNOTATION_ROOT="${ANNOTATION_ROOT:-$WORKSPACE_ROOT/results/annotation}"

PROJECT_ROOT="${PROJECT_ROOT:-$ANNOTATION_ROOT}"
INDIR="${INDIR:-$ANNOTATION_ROOT/mitogenomes_medaka_rrnS}"
OUTDIR="${OUTDIR:-$ANNOTATION_ROOT/medaka_mitos2_rrnS}"

LABEL="MITOS2"
GENETIC_CODE="${GENETIC_CODE:-5}"

MITOS_REFDIR="${MITOS_REFDIR:-$WORKSPACE_ROOT/resources/mitos2_refdata}"
MITOS_REFSEQVER="${MITOS_REFSEQVER:-refseq89m}"
MITOS_REF_TARBALL="${MITOS_REFSEQVER}.tar.bz2"

MITOS_REF_URLS=(
    "https://zenodo.org/record/3685310/files/${MITOS_REF_TARBALL}?download=1"
    "https://zenodo.org/records/3685310/files/${MITOS_REF_TARBALL}?download=1"
    "https://zenodo.org/record/4284483/files/${MITOS_REF_TARBALL}?download=1"
    "https://zenodo.org/records/4284483/files/${MITOS_REF_TARBALL}?download=1"
)

LOGDIR="${OUTDIR}/logs"
MANIFEST="${OUTDIR}/medaka_${LABEL}_rrnS_manifest.tsv"
FAILED="${OUTDIR}/medaka_${LABEL}_rrnS_failed.tsv"

mkdir -p "$OUTDIR" "$LOGDIR" "$MITOS_REFDIR"

have_cmd() {
    command -v "$1" >/dev/null 2>&1
}

download_file() {
    local out="$1"
    shift
    local urls=("$@")

    for url in "${urls[@]}"; do
        echo "[INFO] Downloading: $url"

        if have_cmd wget; then
            if wget -O "$out" "$url"; then
                return 0
            fi
        elif have_cmd curl; then
            if curl -L -o "$out" "$url"; then
                return 0
            fi
        else
            echo "ERROR: Neither wget nor curl is available."
            return 1
        fi
    done

    return 1
}

detect_runmitos() {
    if have_cmd runmitos.py; then
        command -v runmitos.py
        return 0
    fi

    local candidates=(
        "${CONDA_PREFIX:-}/bin/runmitos.py"
        "$HOME/.local/bin/runmitos.py"
        "$HOME/anaconda3/envs/mitos2_t/bin/runmitos.py"
        "$HOME/anaconda3/envs/mitos2/bin/runmitos.py"
    )

    local candidate
    for candidate in "${candidates[@]}"; do
        if [[ -f "$candidate" ]]; then
            echo "$candidate"
            return 0
        fi
    done

    return 1
}

if [[ -n "${CONDA_PREFIX:-}" ]]; then
    export PATH="${CONDA_PREFIX}/bin:${PATH}"
    hash -r
fi

if [[ ! -d "$INDIR" ]]; then
    echo "ERROR: Input directory does not exist: $INDIR"
    exit 1
fi

RUNMITOS="$(detect_runmitos || true)"

if [[ -z "$RUNMITOS" ]]; then
    echo "ERROR: Could not find runmitos.py"
    exit 1
fi

chmod +x "$RUNMITOS" 2>/dev/null || true

if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/python" ]]; then
    PYTHON_BIN="${CONDA_PREFIX}/bin/python"
else
    PYTHON_BIN=""
fi

if [[ ! -d "${MITOS_REFDIR}/${MITOS_REFSEQVER}" ]]; then
    cd "$MITOS_REFDIR"

    if [[ ! -s "$MITOS_REF_TARBALL" ]]; then
        if ! download_file "$MITOS_REF_TARBALL" "${MITOS_REF_URLS[@]}"; then
            echo "ERROR: Could not download ${MITOS_REF_TARBALL}"
            exit 1
        fi
    fi

    tar -xjf "$MITOS_REF_TARBALL"

    if [[ ! -d "${MITOS_REFDIR}/${MITOS_REFSEQVER}" ]]; then
        echo "ERROR: ${MITOS_REFDIR}/${MITOS_REFSEQVER} was not created"
        exit 1
    fi
fi

echo -e "genus\tsample\tinput_fasta\toutput_dir\tstatus" > "$MANIFEST"
echo -e "genus\tsample\tinput_fasta\toutput_dir\terror_log" > "$FAILED"

N_FASTA=$(find "$INDIR" -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) | wc -l)

if [[ "$N_FASTA" -eq 0 ]]; then
    echo "ERROR: No FASTA files were found in: $INDIR"
    exit 1
fi

run_mitos2_one() {
    local fasta="$1"
    local outdir="$2"
    local stdout_log="$3"
    local stderr_log="$4"

    if [[ -n "$PYTHON_BIN" ]]; then
        "$PYTHON_BIN" "$RUNMITOS" \
            -i "$fasta" \
            -c "$GENETIC_CODE" \
            -o "$outdir" \
            --linear \
            -r "$MITOS_REFSEQVER" \
            -R "$MITOS_REFDIR" \
            > "$stdout_log" 2> "$stderr_log"
    else
        "$RUNMITOS" \
            -i "$fasta" \
            -c "$GENETIC_CODE" \
            -o "$outdir" \
            --linear \
            -r "$MITOS_REFSEQVER" \
            -R "$MITOS_REFDIR" \
            > "$stdout_log" 2> "$stderr_log"
    fi
}

find "$INDIR" -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) | sort | while read -r FASTA; do

    RELPATH="${FASTA#$INDIR/}"

    GENUS="$(echo "$RELPATH" | cut -d'/' -f1)"
    SAMPLE="$(echo "$RELPATH" | cut -d'/' -f2)"

    SAMPLE_OUTDIR="${OUTDIR}/${GENUS}/${SAMPLE}"
    MITOS_OUT="${SAMPLE_OUTDIR}/${LABEL}"

    mkdir -p "$MITOS_OUT"

    BASENAME="$(basename "$FASTA")"
    PREFIX="${BASENAME%.*}"

    SAFE_GENUS="$(echo "$GENUS" | tr '/ ' '__')"
    SAFE_SAMPLE="$(echo "$SAMPLE" | tr '/ ' '__')"
    SAFE_PREFIX="$(echo "$PREFIX" | tr '/ ' '__')"

    STDOUT_LOG="${LOGDIR}/${SAFE_GENUS}_${SAFE_SAMPLE}_${SAFE_PREFIX}.${LABEL}.stdout.log"
    STDERR_LOG="${LOGDIR}/${SAFE_GENUS}_${SAFE_SAMPLE}_${SAFE_PREFIX}.${LABEL}.stderr.log"

    echo "[INFO] Processing: $GENUS / $SAMPLE"

    if find "$MITOS_OUT" -maxdepth 1 -type f \( -name "*.gff" -o -name "*.bed" -o -name "*.fas" -o -name "*.fasta" -o -name "*.gb" -o -name "*.tbl" \) | grep -q .; then
        echo -e "${GENUS}\t${SAMPLE}\t${FASTA}\t${MITOS_OUT}\tskipped_existing" >> "$MANIFEST"
        continue
    fi

    rm -rf "$MITOS_OUT"
    mkdir -p "$MITOS_OUT"

    set +e
    run_mitos2_one "$FASTA" "$MITOS_OUT" "$STDOUT_LOG" "$STDERR_LOG"
    EXIT_CODE=$?
    set -e

    if [[ "$EXIT_CODE" -eq 0 ]]; then
        echo -e "${GENUS}\t${SAMPLE}\t${FASTA}\t${MITOS_OUT}\tok" >> "$MANIFEST"
    else
        echo "[ERROR] MITOS2 failed for $SAMPLE"
        tail -n 30 "$STDERR_LOG" || true
        echo -e "${GENUS}\t${SAMPLE}\t${FASTA}\t${MITOS_OUT}\t${STDERR_LOG}" >> "$FAILED"
        echo -e "${GENUS}\t${SAMPLE}\t${FASTA}\t${MITOS_OUT}\tfailed" >> "$MANIFEST"
    fi

done

echo "Manifest: $MANIFEST"
echo "Failed: $FAILED"
echo "Logs: $LOGDIR"
echo "Output: $OUTDIR"
