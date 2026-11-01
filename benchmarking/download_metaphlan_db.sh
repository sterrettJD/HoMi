#!/bin/bash
# Downloads the MetaPhlAn database and its prebuilt Bowtie2 index (about 26 GB)
# to the metaphlan_bowtie_db path in the HoMi configs.
# MetaPhlAn 4.0.3 in HoMi builds the Bowtie2 index if it is missing. The build
# needs more than 28 GB of memory, so this script downloads the prebuilt index.
# HUMAnN 3.9 accepts only MetaPhlAn profiles made with a vJun23 database.
# The HoMi configs set metaphlan_index_name to the same DB_NAME.
# Usage, from benchmarking/: bash download_metaphlan_db.sh
set -euo pipefail

DB_DIR=synthetic/data/metaphlan_db
DB_NAME=mpa_vJun23_CHOCOPhlAnSGB_202307
BASE_URL=http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases
REQUIRED_FILES="$DB_NAME.pkl
$DB_NAME.1.bt2l $DB_NAME.2.bt2l $DB_NAME.3.bt2l $DB_NAME.4.bt2l
$DB_NAME.rev.1.bt2l $DB_NAME.rev.2.bt2l"

mkdir -p "$DB_DIR"
wget -c -P "$DB_DIR" "$BASE_URL/$DB_NAME.tar" "$BASE_URL/$DB_NAME.md5" \
    "$BASE_URL/bowtie2_indexes/${DB_NAME}_bt2.tar" "$BASE_URL/bowtie2_indexes/${DB_NAME}_bt2.md5"

# The md5 files contain relative file names, so run md5sum in DB_DIR
(cd "$DB_DIR" && md5sum -c "$DB_NAME.md5" "${DB_NAME}_bt2.md5")

tar -xf "$DB_DIR/$DB_NAME.tar" -C "$DB_DIR"
tar -xf "$DB_DIR/${DB_NAME}_bt2.tar" -C "$DB_DIR"
rm "$DB_DIR/$DB_NAME.tar" "$DB_DIR/${DB_NAME}_bt2.tar"

for file in $REQUIRED_FILES; do
    if [[ ! -f "$DB_DIR/$file" ]]; then
        echo "ERROR: $DB_DIR/$file is missing." >&2
        exit 1
    fi
done
echo "MetaPhlAn database $DB_NAME is complete in $DB_DIR"
