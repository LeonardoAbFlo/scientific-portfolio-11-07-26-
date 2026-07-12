#!/usr/bin/env bash
set -Eeuo pipefail

FORCE_FIELD_DIR="${FORCE_FIELD_DIR:-$PWD/force-fields}"
REPOSITORY="${FORCE_FIELD_REPOSITORY:-https://github.com/intbio/gromacs_ff.git}"
FORCE_FIELD_NAME="${FORCE_FIELD_NAME:-charmm36-mar2019.ff}"

command -v git >/dev/null 2>&1 || { echo "git is required" >&2; exit 1; }
mkdir -p "$FORCE_FIELD_DIR"

checkout="$FORCE_FIELD_DIR/gromacs_ff"
if [[ ! -d "$checkout/.git" ]]; then
  git clone --depth 1 "$REPOSITORY" "$checkout"
else
  git -C "$checkout" pull --ff-only
fi

source_force_field="$checkout/$FORCE_FIELD_NAME"
[[ -d "$source_force_field" ]] || { echo "Force field not found: $source_force_field" >&2; exit 1; }
destination="$FORCE_FIELD_DIR/$FORCE_FIELD_NAME"
if [[ ! -d "$destination" ]]; then
  cp -R "$source_force_field" "$destination"
fi

printf 'Force field ready: %s\n' "$destination"
printf 'For later stages, run: export GMXLIB=%q\n' "$FORCE_FIELD_DIR"

