#!/bin/bash
# Downloads the Kraken2 Standard-16 database (about 12 GB download, 16 GB extracted)
# to the kraken_db path in the HoMi configs.
# HoMi downloads core_nt (about 240 GB) if one of the required files is missing,
# so this script stops with an error if the extracted database is not complete.
# Usage, from benchmarking/: bash download_kraken_standard16.sh
set -euo pipefail

DB_DIR=synthetic/data/kraken2_db_2
DB_NAME=k2_standard_16gb_20250402
REQUIRED_FILES="hash.k2d opts.k2d taxo.k2d seqid2taxid.map database150mers.kmer_distrib"

mkdir -p "$DB_DIR"
wget -c -O "$DB_DIR/$DB_NAME.tar.gz" "https://genome-idx.s3.amazonaws.com/kraken/$DB_NAME.tar.gz"
tar -xzf "$DB_DIR/$DB_NAME.tar.gz" -C "$DB_DIR"
rm "$DB_DIR/$DB_NAME.tar.gz"

for file in $REQUIRED_FILES; do
    if [[ ! -f "$DB_DIR/$file" ]]; then
        echo "ERROR: $DB_DIR/$file is missing." >&2
        exit 1
    fi
done
echo "Kraken2 Standard-16 database is complete in $DB_DIR"
