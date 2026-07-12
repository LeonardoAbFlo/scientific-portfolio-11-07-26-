#!/usr/bin/env bash
set -euo pipefail

# Polish selected assemblies with Dorado and Medaka.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-$SCRIPT_DIR}"
RUNS_ROOT="${RUNS_ROOT:-$WORKSPACE_ROOT/runs}"

BASE_DIR="${BASE_DIR:-$WORKSPACE_ROOT}"
READS_ROOT="${READS_ROOT:-$RUNS_ROOT}"
READS_LAYOUT="${READS_LAYOUT:-auto}"
ASSEMBLY_ROOT="${ASSEMBLY_ROOT:-$WORKSPACE_ROOT/results/assemblies/selected}"
OUTPUT_BASE="${OUTPUT_BASE:-$WORKSPACE_ROOT/results/polishing/dorado_medaka}"
DORADO_MODELS_DIR="${DORADO_MODELS_DIR:-$WORKSPACE_ROOT/resources/models}"

THREADS="${THREADS:-$(nproc 2>/dev/null || printf '8')}"
DEVICE="${DEVICE:-auto}"
FORCE="${FORCE:-0}"
SAMPLES="${SAMPLES:-}"

DORADO_ENV="${DORADO_ENV:-dorado}"
MEDAKA_ENV="${MEDAKA_ENV:-medaka_gpu}"
SAMTOOLS_ENV="${SAMTOOLS_ENV:-$MEDAKA_ENV}"

# Use this model tag in the read group metadata for Dorado polish.
BASECALL_MODEL="${BASECALL_MODEL:-dna_r10.4.1_e8.2_400bps_sup@v5.0.0}"

# Keep this aligned with the Medaka installation in MEDAKA_ENV.
MEDAKA_MODEL="${MEDAKA_MODEL:-r1041_e82_400bps_sup_v4.2.0}"

DORADO_ALIGNER_OPTS="${DORADO_ALIGNER_OPTS:-}"
DORADO_POLISH_OPTS="${DORADO_POLISH_OPTS:-}"
MEDAKA_OPTS="${MEDAKA_OPTS:-}"
DORADO_QUALITIES="${DORADO_QUALITIES:-1}"
KEEP_RAW_BAM="${KEEP_RAW_BAM:-1}"
KEEP_RG_BAM="${KEEP_RG_BAM:-1}"

for tool_dir in \
    "$HOME/bin" \
    "$HOME/Desktop/RZ/bin" \
    "$HOME/anaconda3/condabin" \
    "$HOME/anaconda3/bin"
do
    if [[ -d "$tool_dir" ]]; then
        PATH="$tool_dir:$PATH"
    fi
done
export PATH

usage() {
    cat <<EOF
Usage:
  $(basename "$0")

Key variables:
    READS_ROOT=/path/to/runs_or_mapq60
    READS_LAYOUT=auto|workspace|legacy|flat
  THREADS=32
  DEVICE=auto|cuda:all|cuda:0|cpu
  FORCE=1
  SAMPLES="F1_barcode81 barcode90 F3/barcode92"
  OUTPUT_BASE=/path/to/output
  DORADO_MODELS_DIR=/path/to/dorado_models
  BASECALL_MODEL=dna_r10.4.1_e8.2_400bps_sup@v5.0.0
  MEDAKA_MODEL=r1041_e82_400bps_sup_v4.2.0
  DORADO_QUALITIES=0                 Skip consensus.fastq from Dorado.
  DORADO_ALIGNER_OPTS="..."
  DORADO_POLISH_OPTS="..."
  MEDAKA_OPTS="..."
  KEEP_RAW_BAM=0
  KEEP_RG_BAM=0
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

mkdir -p "$OUTPUT_BASE"
mkdir -p "$DORADO_MODELS_DIR"
RUN_LOG="$OUTPUT_BASE/polish_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$RUN_LOG") 2>&1

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

env_exists() {
    local mgr="$1"
    local env="$2"
    "$mgr" run -n "$env" true >/dev/null 2>&1
}

env_has_tools() {
    local mgr="$1"
    local env="$2"
    shift 2
    "$mgr" run -n "$env" bash -lc 'for tool in "$@"; do command -v "$tool" >/dev/null 2>&1 || exit 1; done' bash "$@" >/dev/null 2>&1
}

find_env_runner() {
    local env="$1"
    local mgr
    for mgr in micromamba mamba conda; do
        if command -v "$mgr" >/dev/null 2>&1 && env_exists "$mgr" "$env"; then
            printf '%s\n' "$mgr"
            return 0
        fi
    done
    return 1
}

run_in_env() {
    local mgr="$1"
    local env="$2"
    shift 2

    if [[ "$mgr" == "current" ]]; then
        "$@"
    elif [[ "$mgr" == "path" ]]; then
        PATH="$env/bin:$PATH" "$@"
    else
        "$mgr" run -n "$env" "$@"
    fi
}

runner_label() {
    local mgr="$1"
    local env="$2"

    if [[ "$mgr" == "current" ]]; then
        printf 'active env'
    elif [[ "$mgr" == "path" ]]; then
        printf 'PATH=%s/bin' "$env"
    else
        printf '%s run -n %s' "$mgr" "$env"
    fi
}

find_env_dir_with_tools() {
    local env="$1"
    local prefix
    local tool
    shift

    for prefix in \
        "$HOME/micromamba/envs/$env" \
        "$HOME/anaconda3/envs/$env" \
        "$HOME/.local/share/r-miniconda/envs/$env" \
        "$HOME/pymol/envs/$env"
    do
        [[ -d "$prefix" ]] || continue
        for tool in "$@"; do
            [[ -x "$prefix/bin/$tool" ]] || break
        done
        if [[ "$tool" == "${!#}" && -x "$prefix/bin/$tool" ]]; then
            printf '%s\n' "$prefix"
            return 0
        fi
    done

    return 1
}

is_current() {
    local output="$1"
    local input
    shift

    [[ -s "$output" ]] || return 1
    for input in "$@"; do
        [[ "$output" -nt "$input" ]] || return 1
    done
}

strip_fastq_ext() {
    local name="$1"
    name="${name%.fastq.gz}"
    name="${name%.fq.gz}"
    name="${name%.fastq}"
    name="${name%.fq}"
    printf '%s\n' "$name"
}

parse_flat_fastq_name() {
    local name="$1"
    local normalized="${name,,}"

    if [[ "$normalized" =~ ^f([0-9]+)b([0-9]+)(\..+)?$ ]]; then
        printf 'F%s\tbarcode%s\n' "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"
        return 0
    fi

    return 1
}

sample_selected() {
    local flye_set="$1"
    local barcode="$2"
    local sample="$3"
    local wanted

    [[ -z "$SAMPLES" ]] && return 0

    for wanted in ${SAMPLES//,/ }; do
        if [[ "$wanted" == "$sample" || "$wanted" == "$barcode" || "$wanted" == "$flye_set/$barcode" ]]; then
            return 0
        fi
    done

    return 1
}

fastq_to_fasta() {
    local fastq="$1"
    local fasta="$2"

    awk 'NR % 4 == 1 {sub(/^@/, ">"); print; next} NR % 4 == 2 {print}' "$fastq" > "$fasta"
}

validate_medaka_model() {
    local mgr="$1"
    local model="$2"
    local models
    local models_spaced
    local default_model

    [[ -z "$model" || -f "$model" ]] && return 0

    if ! models="$(run_in_env "$mgr" "$MEDAKA_ENV" medaka tools list_models 2>/dev/null)"; then
        log "WARN: could not list Medaka models; will try '$model'."
        return 0
    fi

    models_spaced="${models//,/ }"
    if [[ " $models_spaced " == *" $model "* ]]; then
        return 0
    fi

    default_model="$(sed -n 's/^Default consensus:[[:space:]]*//p' <<< "$models" | head -1)"
    if [[ -n "$default_model" ]]; then
        die "Medaka model '$model' is not available in '$MEDAKA_ENV'. Use MEDAKA_MODEL=$default_model for this installation."
    fi

    die "Medaka model '$model' is not available in '$MEDAKA_ENV'. Check medaka tools list_models."
}

run_dorado_aligner() {
    local reads="$1"
    local draft="$2"
    local bam_raw="$3"
    local log_file="$4"
    local tmp_bam="${bam_raw}.tmp.$$"
    local aligner_extra=()

    if [[ -n "$DORADO_ALIGNER_OPTS" ]]; then
        read -r -a aligner_extra <<< "$DORADO_ALIGNER_OPTS"
    fi

    if [[ "$FORCE" != "1" ]] && is_current "$bam_raw" "$reads" "$draft"; then
        log "Dorado aligner: reusing existing BAM $bam_raw"
        return 0
    fi

    log "Dorado aligner: aligning reads against draft"
    run_in_env "$DORADO_MGR" "$DORADO_ENV" dorado aligner \
        --threads "$THREADS" \
        "${aligner_extra[@]}" \
        "$draft" "$reads" \
        > "$tmp_bam" 2> "$log_file"
    mv "$tmp_bam" "$bam_raw"

    run_in_env "$SAMTOOLS_MGR" "$SAMTOOLS_ENV" samtools quickcheck "$bam_raw"
}

prepare_bam_for_polish() {
    local sample="$1"
    local bam_raw="$2"
    local bam_rg="$3"
    local bam_sorted="$4"
    local log_file="$5"
    local tmp_rg="${bam_rg}.tmp.$$"
    local tmp_sorted="${bam_sorted}.tmp.$$"

    if [[ "$FORCE" != "1" ]] && is_current "$bam_rg" "$bam_raw"; then
        log "Samtools RG: reusing existing BAM $bam_rg"
    else
        log "Samtools RG: adding read group with BASECALL_MODEL=$BASECALL_MODEL"
        run_in_env "$SAMTOOLS_MGR" "$SAMTOOLS_ENV" samtools addreplacerg \
            -r "ID:${sample}" \
            -r "SM:${sample}" \
            -r "PL:ONT" \
            -r "DS:basecall_model=${BASECALL_MODEL}" \
            -o "$tmp_rg" \
            "$bam_raw" \
            >> "$log_file" 2>&1
        mv "$tmp_rg" "$bam_rg"
    fi

    if [[ "$FORCE" != "1" ]] && is_current "$bam_sorted" "$bam_rg"; then
        log "Samtools sort: reusing existing sorted BAM $bam_sorted"
    else
        log "Samtools sort: sorting BAM"
        run_in_env "$SAMTOOLS_MGR" "$SAMTOOLS_ENV" samtools sort \
            -@ "$THREADS" \
            -o "$tmp_sorted" \
            "$bam_rg" \
            >> "$log_file" 2>&1
        mv "$tmp_sorted" "$bam_sorted"
    fi

    if [[ "$FORCE" == "1" || ! -s "${bam_sorted}.bai" || "$bam_sorted" -nt "${bam_sorted}.bai" ]]; then
        log "Samtools index: indexing BAM"
        run_in_env "$SAMTOOLS_MGR" "$SAMTOOLS_ENV" samtools index "$bam_sorted" >> "$log_file" 2>&1
    fi

    run_in_env "$SAMTOOLS_MGR" "$SAMTOOLS_ENV" samtools quickcheck "$bam_sorted"
    if ! run_in_env "$SAMTOOLS_MGR" "$SAMTOOLS_ENV" samtools view -H "$bam_sorted" | grep -q '^@RG'; then
        die "Missing @RG header in $bam_sorted"
    fi

    if [[ "$KEEP_RAW_BAM" == "0" ]]; then
        rm -f "$bam_raw"
    fi
    if [[ "$KEEP_RG_BAM" == "0" ]]; then
        rm -f "$bam_rg"
    fi
}

run_dorado_polish() {
    local bam_sorted="$1"
    local draft="$2"
    local dorado_dir="$3"
    local dorado_draft_fasta="$4"
    local log_file="$5"
    local dorado_extra=()
    local quality_args=()
    local consensus_fasta="$dorado_dir/consensus.fasta"
    local consensus_fastq="$dorado_dir/consensus.fastq"

    if [[ -n "$DORADO_POLISH_OPTS" ]]; then
        read -r -a dorado_extra <<< "$DORADO_POLISH_OPTS"
    fi
    if [[ "$DORADO_QUALITIES" == "1" ]]; then
        quality_args=(--qualities)
    fi

    if [[ "$FORCE" != "1" ]] && is_current "$dorado_draft_fasta" "$bam_sorted" "$draft"; then
        log "Dorado polish: reusing existing draft for Medaka $dorado_draft_fasta"
        return 0
    fi

    if [[ -d "$dorado_dir" && "$FORCE" == "1" ]]; then
        local backup_dir="${dorado_dir}.previous.$(date +%Y%m%d_%H%M%S)"
        log "Dorado polish: FORCE=1, moving previous output to $backup_dir"
        mv "$dorado_dir" "$backup_dir"
    fi
    mkdir -p "$dorado_dir"

    if [[ "$FORCE" == "1" || ! -s "$consensus_fasta" && ! -s "$consensus_fastq" ]]; then
        log "Dorado polish: running dorado polish"
        run_in_env "$DORADO_MGR" "$DORADO_ENV" dorado polish \
            --threads "$THREADS" \
            --device "$DEVICE" \
            --models-directory "$DORADO_MODELS_DIR" \
            "${quality_args[@]}" \
            --output-dir "$dorado_dir" \
            "${dorado_extra[@]}" \
            "$bam_sorted" "$draft" \
            > "$log_file" 2>&1
    else
        log "Dorado polish: reusing existing output in $dorado_dir"
    fi

    if [[ -s "$consensus_fasta" ]]; then
        cp "$consensus_fasta" "$dorado_draft_fasta"
    elif [[ -s "$consensus_fastq" ]]; then
        log "Dorado polish: converting consensus.fastq to FASTA for Medaka"
        fastq_to_fasta "$consensus_fastq" "$dorado_draft_fasta"
    else
        log "Files generated by Dorado:"
        find "$dorado_dir" -maxdepth 2 -type f -printf '  %p\n' >&2
        die "Dorado polish did not produce consensus.fasta or consensus.fastq in $dorado_dir"
    fi

    [[ -s "$dorado_draft_fasta" ]] || die "Could not prepare Dorado FASTA for Medaka: $dorado_draft_fasta"
}

run_medaka() {
    local reads="$1"
    local draft="$2"
    local medaka_dir="$3"
    local final_fasta="$4"
    local log_file="$5"
    local model_args=()
    local medaka_extra=()

    if [[ -n "$MEDAKA_MODEL" ]]; then
        model_args=(-m "$MEDAKA_MODEL")
    fi
    if [[ -n "$MEDAKA_OPTS" ]]; then
        read -r -a medaka_extra <<< "$MEDAKA_OPTS"
    fi

    if [[ "$FORCE" != "1" ]] && is_current "$final_fasta" "$reads" "$draft"; then
        log "Medaka: reusing existing consensus $final_fasta"
        return 0
    fi

    if [[ "$FORCE" != "1" ]] && is_current "$medaka_dir/consensus.fasta" "$reads" "$draft"; then
        log "Medaka: copying existing consensus from $medaka_dir/consensus.fasta"
        cp "$medaka_dir/consensus.fasta" "$final_fasta"
        return 0
    fi

    if [[ -e "$medaka_dir" ]]; then
        if [[ "$FORCE" == "1" ]]; then
            local backup_dir="${medaka_dir}.previous.$(date +%Y%m%d_%H%M%S)"
            log "Medaka: FORCE=1, moving previous output to $backup_dir"
            mv "$medaka_dir" "$backup_dir"
        else
            die "Found $medaka_dir but no current consensus. Review it or rerun with FORCE=1."
        fi
    fi

    log "Medaka: polishing the Dorado consensus"
    run_in_env "$MEDAKA_MGR" "$MEDAKA_ENV" medaka_consensus \
        -i "$reads" \
        -d "$draft" \
        -o "$medaka_dir" \
        -t "$THREADS" \
        "${model_args[@]}" \
        "${medaka_extra[@]}" \
        > "$log_file" 2>&1

    [[ -s "$medaka_dir/consensus.fasta" ]] || die "Medaka did not produce consensus.fasta in $medaka_dir"
    cp "$medaka_dir/consensus.fasta" "$final_fasta"
}

detect_reads_layout() {
    local reads_dir
    local reads

    case "$READS_LAYOUT" in
        auto)
            for reads_dir in "$READS_ROOT"/*/filtered_reads; do
                if [[ -d "$reads_dir" ]]; then
                    printf 'workspace\n'
                    return 0
                fi
            done

            for reads_dir in "$READS_ROOT"/F*_p1_filtered_fastq; do
                if [[ -d "$reads_dir" ]]; then
                    printf 'legacy\n'
                    return 0
                fi
            done

            for reads in "$READS_ROOT"/*.fastq.gz "$READS_ROOT"/*.fq.gz "$READS_ROOT"/*.fastq "$READS_ROOT"/*.fq; do
                if [[ -f "$reads" ]]; then
                    printf 'flat\n'
                    return 0
                fi
            done
            ;;
        workspace|legacy|flat)
            printf '%s\n' "$READS_LAYOUT"
            return 0
            ;;
        *)
            die "READS_LAYOUT must be auto, workspace, legacy, or flat; got '$READS_LAYOUT'."
            ;;
    esac

    die "No reads found in READS_ROOT=$READS_ROOT"
}

process_sample() {
    local flye_set="$1"
    local barcode="$2"
    local reads="$3"
    local sample="${flye_set}_${barcode}"
    local assembly="$ASSEMBLY_ROOT/$flye_set/$barcode/assembly.fasta"
    local out_dir
    local aln_dir
    local dorado_dir
    local medaka_dir
    local bam_raw
    local bam_rg
    local bam_sorted
    local align_log
    local samtools_log
    local dorado_log
    local medaka_log
    local dorado_draft
    local final_fasta

    if ! sample_selected "$flye_set" "$barcode" "$sample"; then
        return 0
    fi

    if [[ ! -s "$assembly" ]]; then
        log "Skipping $sample: assembly not found at $assembly"
        skipped=$((skipped + 1))
        return 0
    fi

    log "Processing $sample"
    out_dir="$OUTPUT_BASE/$flye_set/$barcode"
    aln_dir="$out_dir/aln"
    dorado_dir="$out_dir/dorado"
    medaka_dir="$out_dir/medaka"
    mkdir -p "$out_dir" "$aln_dir"

    bam_raw="$aln_dir/${sample}.raw.bam"
    bam_rg="$aln_dir/${sample}.rg.bam"
    bam_sorted="$aln_dir/${sample}.sorted.bam"
    align_log="$aln_dir/dorado_aligner.log"
    samtools_log="$aln_dir/samtools.log"
    dorado_log="$out_dir/dorado_polish.log"
    medaka_log="$out_dir/medaka.log"
    dorado_draft="$out_dir/${sample}.dorado_polished.fasta"
    final_fasta="$out_dir/${sample}.dorado.medaka.consensus.fasta"

    run_dorado_aligner "$reads" "$assembly" "$bam_raw" "$align_log"
    prepare_bam_for_polish "$sample" "$bam_raw" "$bam_rg" "$bam_sorted" "$samtools_log"
    run_dorado_polish "$bam_sorted" "$assembly" "$dorado_dir" "$dorado_draft" "$dorado_log"
    run_medaka "$reads" "$dorado_draft" "$medaka_dir" "$final_fasta" "$medaka_log"

    printf '%s\t%s\t%s\t%s\t%s\n' "$sample" "$reads" "$assembly" "$dorado_draft" "$final_fasta" >> "$MANIFEST"
    processed=$((processed + 1))
    log "Ready: $sample -> $final_fasta"
}

DORADO_MGR="$(find_env_runner "$DORADO_ENV" || true)"
if [[ -z "$DORADO_MGR" ]] && command -v dorado >/dev/null 2>&1; then
    DORADO_MGR="current"
    DORADO_ENV="current"
fi
if [[ -z "$DORADO_MGR" ]]; then
    DORADO_ENV_DIR="$(find_env_dir_with_tools "$DORADO_ENV" dorado || true)"
    if [[ -n "$DORADO_ENV_DIR" ]]; then
        DORADO_MGR="path"
        DORADO_ENV="$DORADO_ENV_DIR"
    fi
fi
[[ -n "$DORADO_MGR" ]] || die "Could not find environment '$DORADO_ENV' with micromamba, mamba, or conda."
if [[ "$DORADO_MGR" != "current" ]]; then
    if [[ "$DORADO_MGR" == "path" ]]; then
        [[ -x "$DORADO_ENV/bin/dorado" ]] || die "Could not find dorado in $DORADO_ENV/bin."
    else
        env_has_tools "$DORADO_MGR" "$DORADO_ENV" dorado || die "Environment '$DORADO_ENV' does not provide dorado."
    fi
fi

MEDAKA_MGR="$(find_env_runner "$MEDAKA_ENV" || true)"
if [[ -z "$MEDAKA_MGR" ]] && command -v medaka_consensus >/dev/null 2>&1; then
    MEDAKA_MGR="current"
    MEDAKA_ENV="current"
fi
if [[ -z "$MEDAKA_MGR" ]]; then
    MEDAKA_ENV_DIR="$(find_env_dir_with_tools "$MEDAKA_ENV" medaka_consensus || true)"
    if [[ -n "$MEDAKA_ENV_DIR" ]]; then
        MEDAKA_MGR="path"
        MEDAKA_ENV="$MEDAKA_ENV_DIR"
    fi
fi
[[ -n "$MEDAKA_MGR" ]] || die "Could not find environment '$MEDAKA_ENV' with micromamba, mamba, or conda."
if [[ "$MEDAKA_MGR" != "current" ]]; then
    if [[ "$MEDAKA_MGR" == "path" ]]; then
        [[ -x "$MEDAKA_ENV/bin/medaka_consensus" ]] || die "Could not find medaka_consensus in $MEDAKA_ENV/bin."
    else
        env_has_tools "$MEDAKA_MGR" "$MEDAKA_ENV" medaka_consensus || die "Environment '$MEDAKA_ENV' does not provide medaka_consensus."
    fi
fi

SAMTOOLS_MGR="$(find_env_runner "$SAMTOOLS_ENV" || true)"
if [[ -z "$SAMTOOLS_MGR" ]] && command -v samtools >/dev/null 2>&1; then
    SAMTOOLS_MGR="current"
    SAMTOOLS_ENV="current"
fi
if [[ -z "$SAMTOOLS_MGR" ]]; then
    SAMTOOLS_ENV_DIR="$(find_env_dir_with_tools "$SAMTOOLS_ENV" samtools || true)"
    if [[ -n "$SAMTOOLS_ENV_DIR" ]]; then
        SAMTOOLS_MGR="path"
        SAMTOOLS_ENV="$SAMTOOLS_ENV_DIR"
    fi
fi
[[ -n "$SAMTOOLS_MGR" ]] || die "Could not find environment '$SAMTOOLS_ENV' with micromamba, mamba, or conda."
if [[ "$SAMTOOLS_MGR" != "current" ]]; then
    if [[ "$SAMTOOLS_MGR" == "path" ]]; then
        [[ -x "$SAMTOOLS_ENV/bin/samtools" ]] || die "Could not find samtools in $SAMTOOLS_ENV/bin."
    else
        env_has_tools "$SAMTOOLS_MGR" "$SAMTOOLS_ENV" samtools || die "Environment '$SAMTOOLS_ENV' does not provide samtools."
    fi
fi

validate_medaka_model "$MEDAKA_MGR" "$MEDAKA_MODEL"

log "Using $(runner_label "$DORADO_MGR" "$DORADO_ENV") for Dorado."
log "Using $(runner_label "$SAMTOOLS_MGR" "$SAMTOOLS_ENV") for samtools."
log "Using $(runner_label "$MEDAKA_MGR" "$MEDAKA_ENV") for Medaka."
log "Threads: $THREADS"
log "Dorado device: $DEVICE"
log "Dorado models dir: $DORADO_MODELS_DIR"
log "Basecall model for RG: $BASECALL_MODEL"
log "Medaka model: $MEDAKA_MODEL"
log "Output: $OUTPUT_BASE"
if [[ -n "$SAMPLES" ]]; then
    log "Sample filter: $SAMPLES"
fi

MANIFEST="$OUTPUT_BASE/polishing_manifest.tsv"
printf 'sample\treads\tassembly\tdorado_draft\tmedaka_consensus\n' > "$MANIFEST"

shopt -s nullglob
processed=0
skipped=0

reads_layout="$(detect_reads_layout)"
log "Reads layout: $reads_layout"

if [[ "$reads_layout" == "workspace" ]]; then
    for reads_dir in "$READS_ROOT"/*/filtered_reads; do
        [[ -d "$reads_dir" ]] || continue
        flye_set="$(basename "$(dirname "$reads_dir")")"

        for reads in "$reads_dir"/*.fastq.gz "$reads_dir"/*.fq.gz "$reads_dir"/*.fastq "$reads_dir"/*.fq; do
            [[ -f "$reads" ]] || continue
            barcode="$(strip_fastq_ext "$(basename "$reads")")"
            process_sample "$flye_set" "$barcode" "$reads"
        done
    done
elif [[ "$reads_layout" == "legacy" ]]; then
    for reads_dir in "$READS_ROOT"/F*_p1_filtered_fastq; do
        [[ -d "$reads_dir" ]] || continue
        flye_set="$(basename "$reads_dir" _p1_filtered_fastq)"

        for reads in "$reads_dir"/*.fastq.gz "$reads_dir"/*.fq.gz "$reads_dir"/*.fastq "$reads_dir"/*.fq; do
            [[ -f "$reads" ]] || continue
            barcode="$(strip_fastq_ext "$(basename "$reads")")"
            process_sample "$flye_set" "$barcode" "$reads"
        done
    done
else
    for reads in "$READS_ROOT"/*.fastq.gz "$READS_ROOT"/*.fq.gz "$READS_ROOT"/*.fastq "$READS_ROOT"/*.fq; do
        [[ -f "$reads" ]] || continue
        flat_name="$(strip_fastq_ext "$(basename "$reads")")"

        if ! IFS=$'\t' read -r flye_set barcode < <(parse_flat_fastq_name "$flat_name"); then
            log "Skipping $(basename "$reads"): name does not match fNbNN[.tag].fastq.gz"
            skipped=$((skipped + 1))
            continue
        fi

        process_sample "$flye_set" "$barcode" "$reads"
    done
fi

echo
echo "========================================"
echo "[DONE] Polishing complete"
echo "Processed: $processed"
echo "Skipped: $skipped"
echo "Output: $OUTPUT_BASE"
echo "Manifest: $MANIFEST"
echo "Log: $RUN_LOG"
echo "========================================"
