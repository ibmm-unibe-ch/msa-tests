#!/bin/bash

touch combined.a3m
for file in default_porter_all/*/*_conf.a3m     # list directories in the form "/tmp/dirname/"
do
    tail -n +3 $file >> combined.a3m
done

# alignment-mode 3 --> https://mmseqs.com/latest/userguide.pdf page 73
mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_99 mmseqs_tmp --min-seq-id 0.99 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_95 mmseqs_tmp --min-seq-id 0.95 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_90 mmseqs_tmp --min-seq-id 0.9 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_80 mmseqs_tmp --min-seq-id 0.8 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_70 mmseqs_tmp --min-seq-id 0.7 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_60 mmseqs_tmp --min-seq-id 0.6 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_50 mmseqs_tmp --min-seq-id 0.5 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1


mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_99_default mmseqs_tmp --min-seq-id 0.99 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_95_default mmseqs_tmp --min-seq-id 0.95 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_90_default mmseqs_tmp --min-seq-id 0.9 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_80_default mmseqs_tmp --min-seq-id 0.8 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_70_default mmseqs_tmp --min-seq-id 0.7 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_60_default mmseqs_tmp --min-seq-id 0.6 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
mmseqs easy-search combined.a3m /scratch/alphafold_database/mmseqs_databases/uniprot/uniprot_db filtered_out_50_default mmseqs_tmp --min-seq-id 0.5 --threads 32 --alignment-mode 3 --format-output "query,qseq" --max-accept 1
