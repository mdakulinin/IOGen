#!/bin/bash

declare -A BLAST_DBS
BLAST_DBS["CARD"]="/home/galanova/students/misha/kulakova_2025/MetaGeneMark/CARD_databases/prot/CARD_protein_db"
#BLAST_DBS["CATALOG"]="catalog.09.05.2021"
#BLAST_DBS["OTHER"]="modified_VFDB_setB_full_pro"

# Папка с входными файлами .fasta
INPUT_DIR="/home/galanova/students/misha/kulakova_2025/MetaGeneMark/aa_fasta_results"

# Главная папка для результатов
MAIN_OUTPUT_DIR="blastp_all_results"
mkdir -p "$MAIN_OUTPUT_DIR"

# Перебираем все базы
for DB_NAME in "${!BLAST_DBS[@]}"; do
    DB_PATH="${BLAST_DBS[$DB_NAME]}"
    OUTPUT_DIR="$MAIN_OUTPUT_DIR/${DB_NAME}_results"
    mkdir -p "$OUTPUT_DIR"

    echo "Поиск по базе $DB_NAME: $DB_PATH"

    # Перебираем .fasta-файлы
    for FASTA_FILE in "$INPUT_DIR"/*.fasta; do
        FILENAME=$(basename "$FASTA_FILE")
        OUTPUT_FILE="$OUTPUT_DIR/${FILENAME%.fasta}_blast.txt"

        echo "  Обрабатываю $FILENAME → $OUTPUT_FILE"
        blastp -query "$FASTA_FILE" -db "$DB_PATH" -out "$OUTPUT_FILE" -outfmt 6 -evalue 1e-10 -num_threads 32
    done
done

echo "Все завершено! Результаты сохранены в: $MAIN_OUTPUT_DIR"
