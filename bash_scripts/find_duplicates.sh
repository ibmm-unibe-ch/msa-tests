#!/bin/bash

INPUT_DIR=default_porter_all
DB_PATH=/scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db
OUTPUT_DIR=data/filter_results

while getopts "i:d:o:" opt; do
	case $opt in
		i) INPUT_DIR="$OPTARG" ;;
		d) DB_PATH="$OPTARG" ;;
		o) OUTPUT_DIR="$OPTARG" ;;
		\?) echo "Usage: $0 [-i input_dir] [-d db_path] [-o output_dir]" >&2; exit 1 ;;
	esac
done

touch combined.a3m
for file in $INPUT_DIR/*/*_conf.a3m     # list directories in the form "/tmp/dirname/"
do
    tail -n +3 $file >> combined.a3m
done

# alignment-mode 3 --> https://mmseqs.com/latest/userguide.pdf page 73
mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_99 mmseqs_tmp --min-seq-id 0.99 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_95 mmseqs_tmp --min-seq-id 0.95 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_90 mmseqs_tmp --min-seq-id 0.9 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_80 mmseqs_tmp --min-seq-id 0.8 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_70 mmseqs_tmp --min-seq-id 0.7 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_60 mmseqs_tmp --min-seq-id 0.6 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_50 mmseqs_tmp --min-seq-id 0.5 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1


mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_99_default mmseqs_tmp --min-seq-id 0.99 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_95_default mmseqs_tmp --min-seq-id 0.95 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_90_default mmseqs_tmp --min-seq-id 0.9 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_80_default mmseqs_tmp --min-seq-id 0.8 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_70_default mmseqs_tmp --min-seq-id 0.7 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_60_default mmseqs_tmp --min-seq-id 0.6 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m $DB_PATH $OUTPUT_DIR/filtered_out_50_default mmseqs_tmp --min-seq-id 0.5 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
