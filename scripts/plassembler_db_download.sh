#!/usr/bin/env bash
# Must be run AFTER initiating conda environment

mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"

cat << 'EOF' > "$CONDA_PREFIX/etc/conda/activate.d/plassembler_db.sh"
#!/bin/bash
if [ ! -d "$CONDA_PREFIX/plassembler_db" ]; then
    echo "Downloading Plassembler DB into $CONDA_PREFIX/plassembler_db ..."
    plassembler download -d "$CONDA_PREFIX/plassembler_db"
fi
EOF

chmod +x "$CONDA_PREFIX/etc/conda/activate.d/plassembler_db.sh"