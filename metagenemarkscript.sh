#!/bin/bash

MAKE_NT=false
MAKE_AA=true


if [ $MAKE_NT = true]; then
    # Создаем папку для результатов (если её нет)
    mkdir -p nt_fasta_results
    
    # Обрабатываем каждый GFF файл
    for gff_file in *.gff; do
        # Получаем имя образца (удаляем расширение .gff)
        sample_name="${gff_file%.gff}"
        
        echo "Обрабатываю файл: $gff_file -> ${sample_name}_nt.fasta"
        
        # Запускаем преобразование
        /home/galanova/students/misha/MetaGeneMark/nt_from_gff.pl < "$gff_file" > "nt_fasta_results/${sample_name}_nt.fasta"
        
        # Проверяем успешность выполнения
        if [ $? -eq 0 ]; then
            echo "Успешно создан: nt_fasta_results/${sample_name}_nt.fasta"
        else
            echo "Ошибка при обработке файла $gff_file"
        fi
    done
    echo "Все файлы обработаны. Результаты в папке nt_fasta_results/"
fi


if [ $MAKE_AA = true ]; then
    # Пути к файлам
    INPUT_DIR="/home/galanova/students/misha/kulakova_2025/MetaGeneMark/MetaGeneMark_gff"
    OUTPUT_DIR="aa_fasta_results"
    SCRIPT="/home/galanova/students/misha/MetaGeneMark/aa_from_gff.pl"  
    
    # Создаем выходную директорию
    mkdir -p "$OUTPUT_DIR"
    
    # Обрабатываем каждый GFF файл
    for gff_file in "$INPUT_DIR"/*.gff; do
        # Получаем имя файла без пути и расширения
        sample_name=$(basename "$gff_file" .gff)
        
        echo "Обрабатываю файл: $gff_file"
        echo "Сохраняю результат в: $OUTPUT_DIR/${sample_name}_aa.fasta"
        
        # Запускаем преобразование
        $SCRIPT < "$gff_file" > "$OUTPUT_DIR/${sample_name}_aa.fasta"
        
        # Проверяем успешность выполнения
        if [ $? -eq 0 ]; then
            echo "Успешно обработан: $sample_name"
        else
            echo "Ошибка при обработке файла: $gff_file"
        fi
    done
    
    echo "Все файлы обработаны. Результаты сохранены в $OUTPUT_DIR/"
fi