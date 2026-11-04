#!/usr/bin/env bash
set -euo pipefail

# input CSV needs columns: date, med_lat, med_lon

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_csv> <output_csv>"
    exit 1
fi

INPUT_PATH="$1"
OUTPUT_PATH="$2"
SCRIPT_PATH="data/extract_ghrsst.py"

# Resolve absolute directory paths
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd)
SCRIPT_NAME=$(basename "$SCRIPT_PATH")

INPUT_DIR=$(cd "$(dirname "$INPUT_PATH")" && pwd)
INPUT_NAME=$(basename "$INPUT_PATH")

# Pack files using absolute directory references
tar -C "$SCRIPT_DIR" -cf - "$SCRIPT_NAME" -C "$INPUT_DIR" "$INPUT_NAME" | ssh mike@delphinus '
  export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"
  TMP_DIR=$(mktemp -d)
  tar -xf - -C "$TMP_DIR"
  uv run "$TMP_DIR/'"$SCRIPT_NAME"'" -i "$TMP_DIR/'"$INPUT_NAME"'"
  rm -rf "$TMP_DIR"
' > "$OUTPUT_PATH"