#!/bin/bash

export CUDA_VISIBLE_DEVICES=0
SEEDS=1
MODELS=5
RECYCLES=3
BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
PARTIAL_T=20
RF_DESIGNS=10
SAMPLING_TEMP=0.5
SEED=6217
CONF_A=/data/jgut/msa-tests/porter/1miqB1qs8B/1qs8B_conf.pdb
PARENT_PATH=/data/jgut/msa-tests/test
HHFILTER_SIMILARITY=99
BLAST_SIMILARITY=100

function get_pdbs() {
	# $FULL_A $CONF_A $LENGTH $OFFSET_A $FULL_FASTA_A
    if true; then #[ ! -f $1 ]; then
        IDENTIFIER=$(basename $1 _full.pdb)
		PAR_PATH="$(dirname "$1")"
		CHAIN=${IDENTIFIER:4:9999}
        PDB_ID=${IDENTIFIER:0:4}
		pdb_fetch ${PDB_ID} | pdb_selmodel -1 | pdb_selchain -$CHAIN | pdb_rplchain -$CHAIN:A | pdb_delhetatm | pdb_reres -1 | pdb_tidy | grep ^ATOM >$1
		pdb_tofasta $1>$5
		wget -q -O $PAR_PATH/$PDB_ID.fasta https://www.rcsb.org/fasta/entry/${PDB_ID}
		PDB_OFFSET=$(python find_pdb_offset.py --realfasta $PAR_PATH/$PDB_ID.fasta --foundfasta $5 --realchain $CHAIN --realstart $4 --length $3)
	fi
	echo "PDB_OFFSET"
	echo $PDB_OFFSET
    if true; then #[ ! -f $2 ]; then
        pdb_selres -$(($PDB_OFFSET+1)):$(($3+$PDB_OFFSET)) $1  >$2
    fi
}

function fold_alpha() {
    colabfold_batch --overwrite-existing-results --random-seed $SEED --num-seeds $SEEDS --num-models $MODELS --num-recycle $RECYCLES --num-relax 0 $1 $2
    BEST=$(find $2 -name $BESTNAME | tail -1)
    cp $BEST $2/best.pdb
}

function prot_MPNN() {
	rm -r ${2}_folder
	micromamba run -n proteinmpnn python ~/GitHub/msa-diffusion/ProteinMPNN/protein_mpnn_run.py --num_seq_per_target 128 --sampling_temp $SAMPLING_TEMP --pdb_path $1 --pdb_path_chains A --out_folder ${2}_folder --seed $SEED --batch_size 1 
    CURR=$(find ${2}_folder/seqs | tail -1)
    cp $CURR $2
}

function score() {
	if [ -d $1 ]; then
		REFERENCE=$1/best.pdb
	else
		REFERENCE=$1
	fi
	if [ -d $2 ]; then
		MODEL=$2/best.pdb
	else
		MODEL=$2
	fi
    podman run --rm -v $(pwd):$(pwd) openstructure:latest compare-structures --model $MODEL --reference $REFERENCE --output $3 --lddt --local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
}

function rf_diffusion() {
	pa_dir="$(dirname "$3")"
	rm -r $pa_dir
	micromamba run -n SE3nv /data/jgut/template-analysis/RFdiffusion/scripts/run_inference.py "contigmap.contigs=[${5}-${5}]" "contigmap.provide_seq=[${5}-${5}]" inference.output_prefix=$3 inference.input_pdb=$1 inference.num_designs=$RF_DESIGNS diffuser.partial_T=$PARTIAL_T
	python select_lowest.py --parentpath $3 --confa $1 --confb $2 --output ${3}_best.pdb
	cp ${3}_best.pdb $4
}

function filter_diff() {
    hhfilter -id $HHFILTER_SIMILARITY -i $1 -o $2 
}

function filter_unk() {
    blastp -remote -db nr -query $1 -out blast_filter_results.csv -outfmt "10 qseqid sseqid ppos positive" -max_hsps 1 -num_alignments 1
	# queryname, ID, hit percentage, hit_LENGTH
	# first,gb|AAP36446.1|,100.00,110
	# second,sp|P01317.2|,100.00,60
	# no-hit is missing
	python remove_blastp.py --query $1 --scores blast_filter_results.csv --output $2 --similarity $BLAST_SIMILARITY
}

PROT_MPNN_A=$PARENT_PATH/conf
PROT_MPNN_A3M_A=$PARENT_PATH/a_conf.a3m
SCORE_AA=$PARENT_PATH/score.json
prot_MPNN $CONF_A $PROT_MPNN_A3M_A
fold_alpha $PROT_MPNN_A3M_A $PROT_MPNN_A
score $CONF_A $PROT_MPNN_A $SCORE_AA
