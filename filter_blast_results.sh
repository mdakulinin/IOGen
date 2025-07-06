#!/bin/bash

# Папка с подкаталогами, содержащими BLAST-результаты по разным базам
INPUT_DIR="/home/galanova/students/misha/kulakova_2025/MetaGeneMark/blastp_all_results"

# Папка для сохранения отфильтрованных результатов
OUTPUT_DIR="/home/galanova/students/misha/kulakova_2025/MetaGeneMark/filtered_blast_results"
mkdir -p "$OUTPUT_DIR"

# Получаем список всех уникальных имен образцов
SAMPLE_NAMES=$(find "$INPUT_DIR" -type f -name "*.txt" | xargs -n1 basename | sed 's/_orf_aa_blast.txt//' | sort | uniq)

# Перебираем каждый образец
for SAMPLE in $SAMPLE_NAMES; do
    MERGED_FILE=$(mktemp)

    echo "Обработка образца: $SAMPLE"

    # Объединяем соответствующие файлы из всех подкаталогов
    for DB_DIR in "$INPUT_DIR"/*; do
        FILE="$DB_DIR/${SAMPLE}_orf_aa_blast.txt"
        if [[ -f "$FILE" ]]; then
            cat "$FILE" >> "$MERGED_FILE"
        fi
    done

    # Фильтрация: длина выравнивания > 100, % идентичности > 90
    FILTERED_FILE="$OUTPUT_DIR/${SAMPLE}_orf_aa_blast_filtered.txt"
    awk '$4 > 100 && $3 > 90' "$MERGED_FILE" > "$FILTERED_FILE"

    # Статистика
    INPUT_LINES=$(wc -l < "$MERGED_FILE")
    OUTPUT_LINES=$(wc -l < "$FILTERED_FILE")
    echo "Файл: $SAMPLE | Совпадений до/после фильтрации: $INPUT_LINES / $OUTPUT_LINES"

    # Удаляем временный файл
    rm "$MERGED_FILE"
done

echo "Готово! Все отфильтрованные файлы сохранены в: $OUTPUT_DIR"
