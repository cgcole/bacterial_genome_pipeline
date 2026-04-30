#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "Error: No active conda environment detected. Activate one first." >&2
    exit 1
fi

DB_DIR="$CONDA_PREFIX/plassembler_db"

if [[ ! -d "$DB_DIR" ]]; then
    echo "Downloading Plassembler DB to $DB_DIR ..."
    plassembler download -d "$DB_DIR"
else
    echo "Plassembler DB already exists at $DB_DIR, skipping download."
fi