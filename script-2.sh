#!/bin/bash
eval "$(/home/galanova/anaconda3/bin/conda shell.bash hook)"

# Константы
MAKE_FASTQC=false
MAKE_TRIMMING=false
MAKE_FASTQC_TRIM=false
MAKE_BWA=false
MAKE_BOWTIE2=false
MAKE_TABLE=false
MAKE_RAR_CURV=false
MAKE_RAR_CURV2=false
MAKE_KRAKEN=false
MAKE_BRACKEN=false
MAKE_SPADES=false
MAKE_MEGAHIT=false
MAKE_QUAST=false
MAKE_METAGENEMARK=true

# Пути
FASTQC_DIR=fastqc_results
MULTIQC_DIR=multiqc_results
FASTQ_DIR=/home/galanova/projects/kulakova_deti/metagenome_may2025

TRIM_DIR=trimming

FASTQC_TRIM_DIR=fastqc_trim_results
MULTIQC_TRIM_DIR=multiqc_trim_results


if [ $MAKE_FASTQC = true ]; then
    conda activate multiqc
    
    # Создание папок, если не существуют
    mkdir -p "$FASTQC_DIR"
    mkdir -p "$MULTIQC_DIR"

    # Запуск FastQC
    fastqc -t 16 -o "$FASTQC_DIR" "$FASTQ_DIR"/*/*/*.fq.gz

    # Запуск MultiQC в директории fastqc
    multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR"
fi

if [ "$MAKE_TRIMMING" = true ]; then
    conda activate multiqc
    
    mkdir -p "$TRIM_DIR"
    
    for r1 in /home/galanova/projects/kulakova_deti/metagenome_may2025/*-metagenome/raw_data_*/*1.fq.gz; do
        r2=${r1/%1.fq.gz/2.fq.gz}
        sample_name=$(basename "$r1" | cut -d'_' -f1)
        base_output="$TRIM_DIR/${sample_name}_"
        
        echo "Processing $sample_name (without adapter trimming)"
        
        trimmomatic PE -phred33 \
            "$r1" "$r2" \
            "${base_output}1_paired.fq.gz" "${base_output}1_unpaired.fq.gz" \
            "${base_output}2_paired.fq.gz" "${base_output}2_unpaired.fq.gz" \
            LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50 \
            -threads 64
    done
fi
    
if [ $MAKE_FASTQC_TRIM = true ]; then
    conda activate multiqc
    
    # Создание папок, если не существуют
    mkdir -p "$FASTQC_TRIM_DIR"
    mkdir -p "$MULTIQC_TRIM_DIR"

    # Запуск FastQC
    fastqc -t 16 -o "$FASTQC_TRIM_DIR" "$TRIM_DIR"/*_paired.fq.gz

    # Запуск MultiQC в директории fastqc
    multiqc "$FASTQC_TRIM_DIR" -o "$MULTIQC_TRIM_DIR"
fi
    
if [ $MAKE_BWA = true ]; then
    conda activate hg38_filter_env
    GENOME_INDEX="/home/galanova/human_genome/hg38.fa"
    READ_DIR="trimming"
    OUTDIR="filtered_samples"
    mkdir -p "$OUTDIR"

    # Сюда запишем процент человеческих ридов по всем образцам
    REJECTED_SAMPLES="human_mapped_percent.txt"
    echo -e "sample\tmapped_percent" > "$REJECTED_SAMPLES"

    for r1 in "$READ_DIR"/*1_paired.fq.gz; do
        r2="${r1/_1_paired.fq.gz/_2_paired.fq.gz}"
        sample=$(basename "$r1" _1_paired.fq.gz)

        echo "Обработка образца: $sample"

        # Выравнивание
        bwa mem -t 64 "$GENOME_INDEX" "$r1" "$r2" | \
        samtools view -@ 64 -b -o "$sample.bam" -

        # Подсчеты
        total_reads=$(samtools view -c "$sample.bam")
        mapped_reads=$(samtools view -c -F 4 "$sample.bam")
        percent=$(echo "$mapped_reads $total_reads" | awk '{printf "%.2f", ($1/$2)*100}')

        echo "Выровнено: $percent%"

        # Лог
        echo -e "$sample\t$percent" >> "$REJECTED_SAMPLES"

        # Извлекаем только невыравненные пары ридов
        samtools view -@ 64 -b -f 12 "$sample.bam" | \
        samtools sort -n -@ 64 -o "$sample.unmapped.bam" -

        # Преобразуем обратно в fastq
        bedtools bamtofastq \
            -i "$sample.unmapped.bam" \
            -fq "$OUTDIR/${sample}_R1_unmapped.fq" \
            -fq2 "$OUTDIR/${sample}_R2_unmapped.fq"

        # Сжимаем
        gzip "$OUTDIR/${sample}_R1_unmapped.fq"
        gzip "$OUTDIR/${sample}_R2_unmapped.fq"

        # Удаляем временные файлы
        rm "$sample.bam" "$sample.unmapped.bam"

        echo "$sample — завершено"
    done
fi


if [ "$MAKE_BOWTIE2" = true ]; then
    conda activate hg38_filter_env
    GENOME_INDEX="/home/galanova/human_genome/GRCh38_index"
    READ_DIR="trimming"
    OUTDIR="filtered_samples_bowtie2"
    mkdir -p "$OUTDIR"

    REJECTED_SAMPLES="human_mapped_percent_bowtie2.txt"
    echo -e "sample\tmapped_percent" > "$REJECTED_SAMPLES"

    for r1 in "$READ_DIR"/*1_paired.fq.gz; do
        r2="${r1/_1_paired.fq.gz/_2_paired.fq.gz}"
        sample=$(basename "$r1" _1_paired.fq.gz)

        echo "Обработка образца: $sample"

        # Выравнивание с сохранением SAM
        bowtie2 -x "$GENOME_INDEX" \
            -1 "$r1" -2 "$r2" \
            -S "$sample.sam" \
            -p 32 --very-sensitive

        # Конвертация SAM в BAM
        samtools view -@ 32 -b -o "$sample.bam" "$sample.sam"
        rm "$sample.sam"

        # Подсчёты
        total_reads=$(samtools view -c "$sample.bam")
        mapped_reads=$(samtools view -c -F 4 "$sample.bam")
        percent=$(echo "$mapped_reads $total_reads" | awk '{printf "%.2f", ($1/$2)*100}')
        echo "Выровнено: $percent%"

        echo -e "$sample\t$percent" >> "$REJECTED_SAMPLES"

        # Извлекаем невыравненные ПАРЫ ридов
        samtools view -@ 32 -b -f 12 "$sample.bam" | \
        samtools sort -n -@ 32 -o "$sample.unmapped.bam" -

        # Обратно в FASTQ
        bedtools bamtofastq \
            -i "$sample.unmapped.bam" \
            -fq "$OUTDIR/${sample}_R1_unmapped.fq" \
            -fq2 "$OUTDIR/${sample}_R2_unmapped.fq"

        gzip "$OUTDIR/${sample}_R1_unmapped.fq"
        gzip "$OUTDIR/${sample}_R2_unmapped.fq"

        # Удаление временных файлов
        rm "$sample.bam" "$sample.unmapped.bam"

        echo "$sample — завершено"
    done
fi


# Сводная таблица
if [ "$MAKE_TABLE" = true ]; then
    echo -e "sample\traw_reads\ttrimmed_reads\tfiltered_reads" > reads_summary.tsv

    for r1 in "$FASTQ_DIR"/*/*/*1.fq.gz; do
        sample=$(basename "$r1" | cut -d'_' -f1)

        # Найти все соответствующие файлы
        trimmed_r1="$TRIM_DIR/${sample}_1_paired.fq.gz"
        filtered_r1="$OUTDIR/${sample}_R1_unmapped.fq.gz"

        # Подсчет числа ридов в сыром файле
        raw_count=$(zcat "$r1" | awk 'END {print NR/4}')

        # Подсчет числа ридов после тримминга
        if [ -f "$trimmed_r1" ]; then
            trim_count=$(zcat "$trimmed_r1" | awk 'END {print NR/4}')
        else
            trim_count=0
        fi

        # Подсчет числа ридов после фильтрации человеческого генома
        if [ -f "$filtered_r1" ]; then
            filt_count=$(zcat "$filtered_r1" | awk 'END {print NR/4}')
        else
            filt_count=0
        fi

        # Запись в таблицу
        echo -e "$sample\t$raw_count\t$trim_count\t$filt_count" >> reads_summary.tsv
    done

    echo "Сводная таблица сохранена в reads_summary.tsv"
fi


if [ $MAKE_RAR_CURV = true ]; then

    KRAKEN2=/home/galanova/kraken2/KRAKEN2/kraken2
    DB=/home/kovtunas/kraken2_pluspf
    READ_DIR="trimming"
    OUTDIR="kraken2_rarefaction"
    mkdir -p "$OUTDIR"

    # Параметры
    CONFIDENCES=(0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.25 0.5 0.75 0.9 1)
    READ_COUNTS=(50000 100000 250000 500000 1000000 1500000 2000000 2500000 3000000 3500000)
    THREADS=16

    # Главный цикл
    for r1 in "$READ_DIR"/*1_paired.fq.gz; do
        r2="${r1/_1_paired.fq.gz/_2_paired.fq.gz}"
        sample=$(basename "$r1" _1_paired.fq.gz)

        for count in "${READ_COUNTS[@]}"; do
            echo "==> [$sample] Сабсэмплируем $count ридов"

            sub_r1="$OUTDIR/${sample}_R1_${count}.fq"
            sub_r2="$OUTDIR/${sample}_R2_${count}.fq"

            zcat "$r1" | head -n $((count * 4)) > "$sub_r1"
            zcat "$r2" | head -n $((count * 4)) > "$sub_r2"

            for conf in "${CONFIDENCES[@]}"; do
                echo "    Kraken2 c confidence=$conf"

                $KRAKEN2 --paired \
                  --db "$DB" \
                  --threads "$THREADS" \
                  --confidence "$conf" \
                  --report "$OUTDIR/${sample}_${count}_c${conf}.report" \
                  --output /dev/null \
                  "$sub_r1" "$sub_r2"
            done
        done
    done

fi


if [ $MAKE_KRAKEN = true ]; then
    KRAKEN2=/home/galanova/kraken2/KRAKEN2/kraken2
    DB=/home/kovtunas/kraken2_pluspf
    
    READ_DIR="/home/galanova/students/misha/kulakova_2025/filtered_samples_bowtie2"
    OUTDIR="kraken2_results_025"
    mkdir -p "$OUTDIR"
    
    THREADS=32
    CONF=0.25
    
    for r1 in "$READ_DIR"/*_R1_unmapped.fq.gz; do
        r2="${r1/_R1_/_R2_}"
        sample=$(basename "$r1" _R1_unmapped.fq.gz)

        echo "==> Обработка образца: $sample"

        $KRAKEN2 --db "$DB" \
            --threads "$THREADS" \
            --confidence "$CONF" \
            --report "$OUTDIR/${sample}.report" \
            --output "$OUTDIR/${sample}.output" \
            --use-names \
            --report-zero-counts \
            --paired "$r1" "$r2"
    done

    echo "Kraken2 завершен для всех образцов"
fi



if [ "$MAKE_RAR_CURV2" = true ]; then
    # Пути
    KRAKEN2=/home/galanova/kraken2/KRAKEN2/kraken2
    DB=/home/kovtunas/kraken2_pluspf
    READ_DIR="filtered_samples"
    OUTDIR="kraken2_rarefaction"
    TMP_FASTQ="/dev/shm/kraken2_rarefaction_tmp"
    mkdir -p "$OUTDIR" "$TMP_FASTQ"
    
    # Параметры
    CONFIDENCES=(0.001 0.005 0.01 0.05 0.1 0.25 0.5 0.75 1)
    READ_COUNTS=(50000 100000 250000 500000 1000000 1500000 2000000)
    THREADS=32

    echo "===> Сабсэмплирование FASTQ и запуск Kraken2..."

    sample_counter=0

    for r1 in "$READ_DIR"/*R1_unmapped.fq.gz; do
        if [ "$sample_counter" = 0 ]; then
            ((sample_counter++))
            continue
        fi
        r2="${r1/_R1_unmapped.fq.gz/_R2_unmapped.fq.gz}"
        sample=$(basename "$r1" _R1_unmapped.fq.gz)

        ((sample_counter++))
        if [ "$sample_counter" -gt 4 ]; then
            break
        fi

        echo ">>> [$sample]"

        for count in "${READ_COUNTS[@]}"; do
            echo "  └─ сабсэмплируем $count ридов"

            sub_r1="$TMP_FASTQ/${sample}_R1_${count}.fq"
            sub_r2="$TMP_FASTQ/${sample}_R2_${count}.fq"

            zcat "$r1" | head -n $((count * 4)) > "$sub_r1"
            zcat "$r2" | head -n $((count * 4)) > "$sub_r2"

            for conf in "${CONFIDENCES[@]}"; do
                echo "      → Kraken2: conf=$conf"

                $KRAKEN2 --db "$DB" \
                    --threads "$THREADS" \
                    --confidence "$conf" \
                    --use-names \
                    --quick \
                    --report-zero-counts \
                    --paired "$sub_r1" "$sub_r2" \
                    --report "$OUTDIR/${sample}_R1_unmapped.fq.gz_${count}_c${conf}.report" \
                    --output /dev/null
            done

            # Удаление текущего сабсэмпла
            rm -f "$sub_r1" "$sub_r2"
        done
    done

    echo "Завершено. Результаты сохранены в $OUTDIR"
fi

if [ $MAKE_BRACKEN = true ]; then

    BRACKEN_DB=/home/kovtunas/kraken2_pluspf
    REPORT_DIR="kraken2_results_025"
    OUTDIR="bracken_output"
    mkdir -p "${OUTDIR}_S_025"
    mkdir -p "${OUTDIR}_G_025"
    mkdir -p "${OUTDIR}_F_025"
    mkdir -p "${OUTDIR}_P_025"
    
    # Параметры
    READ_LEN=150           # Длина ридов
    #LEVEL="S"              # Уровень: S (species), G (genus), P (phylum), etc.
    
    for r1 in "$REPORT_DIR"/*.report; do
    	sample=$(basename "$r1" .report)
    
    	echo "===> Bracken: Обработка отчета $sample на уровне $LEVEL"
    
        bracken -d "$BRACKEN_DB" \
                -i "$r1" \
                -o "${OUTDIR}_S_025/${sample}_bracken_S.txt" \
                -r "$READ_LEN" \
                -l S

        bracken -d "$BRACKEN_DB" \
                -i "$r1" \
                -o "${OUTDIR}_G_025/${sample}_bracken_G.txt" \
                -r "$READ_LEN" \
                -l G

        bracken -d "$BRACKEN_DB" \
                -i "$r1" \
                -o "${OUTDIR}_F_025/${sample}_bracken_F.txt" \
                -r "$READ_LEN" \
                -l F

        bracken -d "$BRACKEN_DB" \
                -i "$r1" \
                -o "${OUTDIR}_P_025/${sample}_bracken_P.txt" \
                -r "$READ_LEN" \
                -l P

    
    done

fi


if [ $MAKE_SPADES = true ]; then
    INPUT_DIR=filtered_samples_bowtie2
    OUTPUT_BASE="spades_output"
    mkdir -p "$OUTPUT_BASE"
    for R1 in "$INPUT_DIR"/*_R1_unmapped.fq.gz; do
    # Получаем соответствующий файл R2
    R2="${R1/_R1_/_R2_}"
    # Проверяем, существует ли второй файл
    if [[ -f "$R2" ]]; then
        # Извлекаем префикс, например: ZVU939
        SAMPLE=$(basename "$R1" | cut -d'-' -f1)

        # Создаем папку для вывода
        OUTDIR="${OUTPUT_BASE}/${SAMPLE}"

        # Запускаем SPAdes
        spades.py --meta --only-assembler -1 "$R1" -2 "$R2" -o "$OUTDIR" -t 32

        echo "Готово: $SAMPLE"
    else
        echo "Файл не найден: $R2"
    fi
    done

fi

if [ $MAKE_MEGAHIT = true ]; then
    INPUT_DIR=filtered_samples_bowtie2
    OUTPUT_BASE="megahit_output"
    mkdir -p "$OUTPUT_BASE"
    for R1 in "$INPUT_DIR"/*_R1_unmapped.fq.gz; do
        # Получаем соответствующий файл R2
        R2="${R1/_R1_/_R2_}"

        # Проверяем, существует ли второй файл
        if [[ -f "$R2" ]]; then
            # Извлекаем префикс, например: ZVU939
            SAMPLE=$(basename "$R1" | cut -d'-' -f1)

            # Создаем папку для вывода
            OUTDIR="${OUTPUT_BASE}/${SAMPLE}"

            # Запускаем MEGAHIT
            megahit -1 "$R1" -2 "$R2" -o "$OUTDIR" -t 32 --presets meta-sensitive

            echo "Готово: $SAMPLE"
        else
            echo "Файл не найден: $R2"
        fi
    done
fi


if [ $MAKE_QUAST = true ]; then
    SPADES_DIR="megahit_output"
    QUAST_OUTDIR="quast_megahit_output"
    mkdir -p "$QUAST_OUTDIR"
    
    for ASSEMBLY_DIR in "$SPADES_DIR"/*; do
        if [[ -d "$ASSEMBLY_DIR" ]]; then
            SAMPLE=$(basename "$ASSEMBLY_DIR")
            CONTIGS="${ASSEMBLY_DIR}/final.contigs.fa"
    
            if [[ -f "$CONTIGS" ]]; then
                OUTDIR="${QUAST_OUTDIR}/${SAMPLE}"
                mkdir -p "$OUTDIR"
    
                quast.py "$CONTIGS" -o "$OUTDIR" -t 16 --min-contig 500
    
                echo "QUAST завершён для: $SAMPLE"
            else
                echo "contigs.fasta не найден в $ASSEMBLY_DIR"
            fi
        fi
    done
fi

if [ $MAKE_METAGENEMARK = true ]; then
    /mnt/storage/kovtunas/_Ira_/pipeline/MetaGeneMark/gmhmmp -a -d -f G -m /mnt/storage/kovtunas/_Ira_/pipeline/MetaGeneMark/MetaGeneMark_v1.mod /home/galanova/students/misha/kulakova_2025/spades_output/AEE519/contigs.fasta -o RESULT.gff

fi







