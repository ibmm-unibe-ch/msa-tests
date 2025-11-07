#!/bin/bash

SEEDS=1
MODELS=5
RECYCLES=3
BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
SAMPLING_TEMP=1
SEED=6217
PARENT_PATH=/data/jgut/msa-tests/single_protein_test
HHFILTER_SIMILARITY=99
BLAST_SIMILARITY=100
MAXIT_PATH=/data/jgut/template-analysis/maxit-v10.200-prod-src/bin/maxit
FASPR_PATH=/data/jgut/template-analysis/FASPR/FASPR

function get_pdbs() {
	# $FULL_A $Alphafold_output $A3M_output
	CHAIN="A"
	IDENTIFIER=$(basename $1 .pdb)
	CHAIN=${IDENTIFIER:4:9999}
	PDB_ID=${IDENTIFIER:0:4}
	pdb_fetch ${PDB_ID} | pdb_selmodel -1 | pdb_selchain -$CHAIN | pdb_rplchain -$CHAIN:A | pdb_delhetatm | pdb_delinsertion | pdb_reres -1 | pdb_tidy | grep ^ATOM >$1
	wget -O $2 https://www.rcsb.org/fasta/entry/${PDB_ID}
	colabfold_batch --overwrite-existing-results --random-seed $SEED --num-seeds 1 --num-models 1 --num-recycle 0 --templates --custom-template-path $PARENT_PATH/temp_path --num-relax 0 $3 $PARENT_PATH/out
	mv $PARENT_PATH/out/${PDB_ID}${CHAIN}_full_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_6217.pdb $2
	rm -r $PARENT_PATH/temp_path
}

function cut_pdb() {
	# 1 --> whole file
	# 2 --> output file
	# 3 --> start
	# 4 --> length
	if (( $4 )) ; then
    	pdb_selres -$(($3+1)):$4 $1 >$2
 	else
    	pdb_selres -$(($3+1)): $1 >$2
 	fi		
}

function fold_alpha() {
	colabfold_batch --overwrite-existing-results --random-seed $SEED --num-seeds $SEEDS --num-models $MODELS --num-recycle $RECYCLES --num-relax 0 $1 $2
	BEST=$(find $2 -name $BESTNAME | tail -1)
	cp $BEST $2/best.pdb
}

function prot_MPNN() {
	rm -r ${2}_folder
	micromamba run -n RF2 python ~/GitHub/msa-diffusion/ProteinMPNN/protein_mpnn_run.py --num_seq_per_target 128 --sampling_temp $SAMPLING_TEMP --pdb_path $1 --pdb_path_chains A --out_folder ${2}_folder --seed $SEED --batch_size 1 
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
	docker run --rm -v $(pwd):$(pwd) -v /scratch/alphafold_database/DEShaw_simulations:/scratch/alphafold_database/DEShaw_simulations registry.scicore.unibas.ch/schwede/openstructure:latest compare-structures --model $MODEL --reference $REFERENCE --output $3 --residue-number-alignment --lddt --local-lddt --bb-lddt --bb-local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
}

function score_both(){
	PAR_PATH=$1
	IDENTIFIER1=$(basename $2 .pdb)
	IDENTIFIER2=$(basename $3 .pdb)
	score $1 $2 ${PAR_PATH}/score_${IDENTIFIER1}.json
	score $1 $3 ${PAR_PATH}/score_${IDENTIFIER2}.json
}

function first_and_rest(){
	touch $3
	head -n 2 $1 > $3
	if [ -d $2 ]; then
		for file in $2/*.a3m ; do                           
    		A3M_file=$file
		done
	else
		A3M_file=$2
	fi
	tail -n +4 $A3M_file >> $3
}

function get_secstruc(){ #$1 --> .pdb #$2 output_sec_struc.json
		pdb_tofasta $1 >$1.fasta
		python utils/get_secstrucs.py --input_pdb $1 --output_json $2
}

while IFS=, read -r GROUP ID_A ID_B START END
do
	CURR_PATH=$PARENT_PATH/$GROUP/${ID_A}${ID_B}
	echo The current path is: $CURR_PATH
	mkdir -p $CURR_PATH
    STRUC_A=$CURR_PATH/${ID_A}.pdb
	FASTA_PATH_A=$CURR_PATH/$ID_A.fasta
    PROT_MPNN_A=$CURR_PATH/${ID_A}_protmpnn
    PROT_MPNN_A3M_A=$CURR_PATH/${ID_A}_protmpnn/prot_mpnn.a3m
    single_seq_A=$CURR_PATH/$ID_A.a3m
    single_out_A=$CURR_PATH/${ID_A}_single
    ALPHAFOLD_A=$CURR_PATH/${ID_A}_full
    STRUC_B=$CURR_PATH/${ID_B}.pdb
    FASTA_PATH_B=$CURR_PATH/$ID_B.fasta
    PROT_MPNN_B=$CURR_PATH/${ID_B}_protmpnn
    PROT_MPNN_A3M_B=$CURR_PATH/${ID_B}_protmpnn/prot_mpnn.a3m
    single_seq_B=$CURR_PATH/$ID_B.a3m
    single_out_B=$CURR_PATH/${ID_B}_single
    ALPHAFOLD_B=$CURR_PATH/${ID_B}_full
    MIXED_A_B=$CURR_PATH/${ID_A}_main_${ID_B}_rest
    MIXED_A_B_MPNN=$CURR_PATH/${ID_A}_main_${ID_B}_protmpnn
	MIXED_B_A=$CURR_PATH/${ID_B}_main_${ID_A}_rest
    MIXED_B_A_MPNN=$CURR_PATH/${ID_B}_main_${ID_A}_protmpnn
	MIXED_A3M_A_B=$CURR_PATH/${ID_A}_main_${ID_B}_rest.a3m
    MIXED_A3M_A_B_MPNN=$CURR_PATH/${ID_A}_main_${ID_B}_protmpnn.a3m
	MIXED_A3M_B_A=$CURR_PATH/${ID_B}_main_${ID_A}_rest.a3m
    MIXED_A3M_B_A_MPNN=$CURR_PATH/${ID_B}_main_${ID_A}_protmpnn.a3m
    get_pdbs $STRUC_A $FASTA_PATH_A || (echo "Problem with $FULL_A"; continue)
	get_pdbs $STRUC_B $FASTA_PATH_B || (echo "Problem with $FULL_B"; continue)
	cut_pdb $STRUC_A process_single_tmp $START $END
	mv process_single_tmp $STRUC_A
	get_secstruc $STRUC_A ${STRUC_A%????}_sec_struc.json
	get_secstruc $STRUC_B ${STRUC_B%????}_sec_struc.json
	pdb_tofasta $STRUC_A > $FASTA_PATH_A
	lower=$(tail -n +2 $FASTA_PATH_A | tr -d '\n')
	higher=$(head -1 $FASTA_PATH_A)
	echo -e $higher"\n"$lower > $FASTA_PATH_A
    prot_MPNN $STRUC_A $PROT_MPNN_A3M_A
	prot_MPNN $STRUC_B $PROT_MPNN_A3M_B
	echo "Classical $ID_A"
	cp $FASTA_PATH_A $single_seq_A
    fold_alpha $single_seq_A $single_out_A
	score_both $single_out_A $STRUC_A $STRUC_B
    fold_alpha $FASTA_PATH_A $ALPHAFOLD_A
    score_both $ALPHAFOLD_A $STRUC_A $STRUC_B
    fold_alpha $PROT_MPNN_A3M_A $PROT_MPNN_A
	score_both $PROT_MPNN_A $STRUC_A $STRUC_B
    echo "Classical $ID_B"
	cp $FASTA_PATH_B $single_seq_B
    fold_alpha $single_seq_B $single_out_B  
	score_both $single_out_B $STRUC_A $STRUC_B
    fold_alpha $FASTA_PATH_B $ALPHAFOLD_B
    score_both $ALPHAFOLD_B $STRUC_A $STRUC_B
    fold_alpha $PROT_MPNN_A3M_B $PROT_MPNN_B
	score_both $PROT_MPNN_B $STRUC_A $STRUC_B
	echo "Mixed $ID_A"
    first_and_rest $single_seq_A $ALPHAFOLD_B $MIXED_A3M_A_B
    first_and_rest $single_seq_A $PROT_MPNN_A3M_B $MIXED_A3M_A_B_MPNN
    fold_alpha $MIXED_A3M_A_B $MIXED_A_B 
    fold_alpha $MIXED_A3M_A_B_MPNN $MIXED_A_B_MPNN
	get_secstruc $MIXED_A_B_MPNN/best.pdb ${MIXED_A_B_MPNN}best_sec_struc.json
    score_both $MIXED_A_B $STRUC_A $STRUC_B
    score_both $MIXED_A_B_MPNN $STRUC_A $STRUC_B
    echo "Mixed $ID_B"
    first_and_rest $single_seq_B $ALPHAFOLD_A $MIXED_A3M_B_A
    first_and_rest $single_seq_B $PROT_MPNN_A3M_A $MIXED_A3M_B_A_MPNN
    fold_alpha $MIXED_A3M_B_A $MIXED_B_A 
    fold_alpha $MIXED_A3M_B_A_MPNN $MIXED_B_A_MPNN
    score_both $MIXED_B_A $STRUC_A $STRUC_B
    score_both $MIXED_B_A_MPNN $STRUC_A $STRUC_B
	echo "Done with $CURR_PATH"
done < single_proteins.csv