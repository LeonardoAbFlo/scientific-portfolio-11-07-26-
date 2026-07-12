#!/usr/bin/env bash
set -euo pipefail

GNINA_VERSION="${GNINA_VERSION:-1.3.2}"
OUTDIR="tools/gnina"
mkdir -p "$OUTDIR"

CANDIDATES=(
  "https://github.com/gnina/gnina/releases/download/v${GNINA_VERSION}/gnina.${GNINA_VERSION}"
  "https://github.com/gnina/gnina/releases/download/v${GNINA_VERSION}/gnina"
  "https://github.com/gnina/gnina/releases/download/${GNINA_VERSION}/gnina.${GNINA_VERSION}"
)

cd "$OUTDIR"

if [[ -x "gnina" ]]; then
  echo "[OK] GNINA already present: $OUTDIR/gnina"
  exit 0
fi

echo "[INFO] Downloading GNINA version: ${GNINA_VERSION}"
for url in "${CANDIDATES[@]}"; do
  echo "  - trying: $url"
  if command -v curl >/dev/null 2>&1; then
    if curl -fsSL "$url" -o gnina_download; then
      break
    fi
  else
    if wget -q "$url" -O gnina_download; then
      break
    fi
  fi
done

if [[ ! -f gnina_download ]]; then
  echo "[ERR] Could not download GNINA. Set GNINA_VERSION or download manually."
  exit 1
fi

# sanity checks
file gnina_download || true
chmod +x gnina_download
mv gnina_download gnina

echo "[OK] GNINA ready: $OUTDIR/gnina"
