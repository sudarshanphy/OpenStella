#!/usr/bin/env bash
set -euo pipefail

MODEL="${1:?give model name, e.g. s15model_comb04}"
INPUT_DIR="${2:-.}"

HYD="${INPUT_DIR}/${MODEL}.hyd"
ABN="${INPUT_DIR}/${MODEL}.abn"
DAT="${INPUT_DIR}/${MODEL}.dat"

for f in "$HYD" "$ABN" "$DAT"; do
    if [ ! -f "$f" ]; then
        echo "Missing file: $f"
        exit 1
    fi
done

cp "$HYD" "modmake/${MODEL}.hyd"
cp "$ABN" "modmake/${MODEL}.abn"
cp "$DAT" "strad/run/${MODEL}.dat"
cp "$DAT" "run/strad/${MODEL}.dat"

echo "Copied:"
echo "  $HYD -> modmake/${MODEL}.hyd"
echo "  $ABN -> modmake/${MODEL}.abn"
echo "  $DAT -> strad/run/${MODEL}.dat"
echo "  $DAT -> run/strad/${MODEL}.dat"

NZON_HYD=$(awk 'NR==1 {print $2}' "$HYD")
NZON_ABN=$(wc -l < "$ABN")

echo
echo "Nzon from hyd header = $NZON_HYD"
echo "lines in abn         = $NZON_ABN"

if [ "$NZON_HYD" != "$NZON_ABN" ]; then
    echo "WARNING: hyd Nzon and abn line count do not match!"
else
    echo "hyd/abn zone counts match."
fi
