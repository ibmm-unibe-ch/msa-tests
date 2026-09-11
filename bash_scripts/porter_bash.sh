#!/bin/bash

export CUDA_VISIBLE_DEVICES=0
SEEDS=1
MODELS=5
RECYCLES=3
BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
PARTIAL_T=20
RF_DESIGNS=10
SAMPLING_TEMP=1
SEED=6217
PARENT_PATH=/data/jgut/msa-tests/aaa_porter_all_models/porter_all_models
SEQ_IDENTITY=0.3
RCSBROOT=/data/jgut/template-analysis/maxit-v10.200-prod-src; export RCSBROOT
MAXIT_PATH=/data/jgut/template-analysis/maxit-v10.200-prod-src/bin/maxit
FASPR_PATH=/data/jgut/template-analysis/FASPR/FASPR
ROSETTAFOLD_PATH=/home/jgut/tools/RoseTTAFold2/network/predict.py
ROSETTAFOLD_WEIGHTS=/home/jgut/tools/RoseTTAFold2/network/RF2_jan24.pt

function get_pdbs() {
	# $FULL_A $Alphafold_output $A3M_output 
	# if [ ! -f $1 ]; then
	IDENTIFIER=$(basename $1 _full.pdb)
	PAR_PATH="$(dirname "$1")"
	CHAIN=${IDENTIFIER:4:9999}
	PDB_ID=${IDENTIFIER:0:4}
	rm -r $PARENT_PATH/temp_path
	mkdir -p $PARENT_PATH/temp_path
	pdb_fetch ${PDB_ID} | pdb_selmodel -1 | pdb_selchain -$CHAIN | pdb_rplchain -$CHAIN:A | pdb_delhetatm | pdb_delinsertion | pdb_reres -1 | pdb_tidy | grep ^ATOM >$PARENT_PATH/temp_path/temp.pdb
	#relax $PARENT_PATH/temp_path/temp.pdb
	#mv $PARENT_PATH/temp_path/temp.pdb_relaxed.pdb $PARENT_PATH/temp_path/temp.pdb 
	cp $PARENT_PATH/temp_path/temp.pdb $1
	#$MAXIT_PATH -input temp_path/temp.pdb -output temp_path/1inp.cif -o 1
	wget -O $1.fasta https://www.rcsb.org/fasta/entry/${PDB_ID}
	python utils/sel_chain.py --fastapath $1.fasta --chain $CHAIN --outputpath $3
	cp $3 $1.fasta
	colabfold_batch --overwrite-existing-results --random-seed $SEED --num-seeds 1 --num-models 1 --num-recycle 0 --templates --custom-template-path $PARENT_PATH/temp_path --num-relax 0 $3 $PARENT_PATH/out
	mv $PARENT_PATH/out/${PDB_ID}${CHAIN}_full_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_6217.pdb $2
	rm -r $PARENT_PATH/temp_path
}

function cut_pdb() {
	# 1 --> whole file
	# 2 --> output file
	# 3 --> start
	# 4 --> length
	pdb_selres -$(($3+1)):$(($3+$4)) $1 >$2
}

function fold_alpha2() {
	rm -r $2
	colabfold_batch --overwrite-existing-results --random-seed $SEED --num-seeds $SEEDS --num-models $MODELS --num-recycle $RECYCLES --num-relax 0 $1 $2
	BEST=$(find $2 -name $BESTNAME | tail -1)
	cp $BEST $2/best.pdb
}

function msa_to_input() {
	PARENT_DIR="$(dirname "$2")"
	mkdir -p $PARENT_DIR
	echo "{">$2
	echo '"name": "Fold",'>>$2
	echo $"\"modelSeeds\": [$SEED],">>$2
	echo '"dialect": "alphafold3",'>>$2
	echo '"version": 1,'>>$2
	echo '"sequences": [{'>>$2
	echo '"protein": {'>>$2
	SECOND="$(sed -n '2p' $1)"
	echo '"id": "A",'>>$2
	echo $"\"sequence\": \"$SECOND\",">>$2
	MSA="$(sed '$!s/$/\\n/' $1 | tr -d '\n')"
	echo $"\"unpairedMsa\": \"$MSA\",">>$2
	echo '"pairedMsa": "",'>>$2
	echo '"templates": []'>>$2
	echo '}}]'>>$2
	echo '}'>>$2
}

function fold_alpha3() {
	msa_to_input $1 ${1%????}/af_input/fold_input.json
	docker context use default
	docker run --runtime=nvidia --gpus device=$CUDA_VISIBLE_DEVICES --rm --volume ${1%????}/af_input:/root/af_input --volume $2:/root/af_output --volume /scratch/alphafold_database/alphafold3_weights:/root/models --volume /scratch/alphafold_database:/root/public_databases alphafold3 python run_alphafold.py --json_path=/root/af_input/fold_input.json --model_dir=/root/models --output_dir=/root/af_output
	declare -a numbers=("0" "1" "2" "3" "4") 
	for number in "${numbers[@]}"
	do
		mv $2/fold/seed-${SEED}_sample-${number}/model.cif $2/sample-${number}.cif
		$MAXIT_PATH -i $2/sample-${number}.cif -output $2/sample-it-${number}.pdb -o 2
	done
}

function prot_MPNN() {
	rm -r ${2}_folder
	micromamba run -n RF2 python ~/GitHub/msa-diffusion/ProteinMPNN/protein_mpnn_run.py --num_seq_per_target 128 --sampling_temp $SAMPLING_TEMP --pdb_path $1 --pdb_path_chains A --out_folder ${2}_folder --seed $SEED --batch_size 1 
	CURR=$(find ${2}_folder/seqs | tail -1)
	cp $CURR $2
}

function align() {
	pa_dir="$(dirname "$3")"
	MODEL_SEQ_PATH=${1}.fasta
	REF_SEQ_PATH=${2}.fasta
	rm $MODEL_SEQ_PATH
	rm $REF_SEQ_PATH
	pdb_tofasta $1>>$MODEL_SEQ_PATH
	pdb_tofasta $2>>$REF_SEQ_PATH
	OFFSET=$(blastp -subject $REF_SEQ_PATH -query $MODEL_SEQ_PATH -outfmt "10 qstart")
	re='^[0-9]+$'
	if ! [[ $OFFSET =~ $re ]] ; then
		return 0
	else
		pdb_reres -$OFFSET $2>$pa_dir/temp
		mv $pa_dir/temp $2
		return 1
	fi
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
	#align $MODEL $REFERENCE
	#ALIGNED_STATE=$?
	#if [ "$ALIGNED_STATE" -eq "1" ]; then
	#	docker run --rm -v $(pwd):$(pwd) registry.scicore.unibas.ch/schwede/openstructure:latest compare-structures --model $MODEL --reference $REFERENCE --output $3 --lddt --local-lddt --bb-lddt --bb-local-lddt --tm-score --rigid-scores --lddt-no-stereochecks -rna
	#else
	docker run --rm -v $(pwd):$(pwd) registry.scicore.unibas.ch/schwede/openstructure:latest compare-structures --model $MODEL --reference $REFERENCE --output $3 --lddt --local-lddt --bb-lddt --bb-local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
	#fi
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
	python utils/remove_blastp.py --query $1 --scores blast_filter_results.csv --output $2 --similarity $BLAST_SIMILARITY
}

function fold_rosetta() {
	PAR_PATH_RF="$(dirname "$1")"
	INPUT_PATH=$PAR_PATH_RF/CURR_ROSETTA_input.a3m
	first_chars=$(head -c 3 $1)
	if [[ $first_chars == \#* ]]; then tail -n +2 $1 > $INPUT_PATH; else cp $1 $INPUT_PATH; fi
	micromamba run -n RF2 python $ROSETTAFOLD_PATH -inputs $INPUT_PATH -prefix $PAR_PATH_RF/CURR_ROSETTA -model $ROSETTAFOLD_WEIGHTS -nseqs 128 -nseqs_full 1
	mv $PAR_PATH_RF/CURR_ROSETTA_00_pred.pdb $2
	rm $PAR_PATH_RF/CURR_ROSETTA_*
}

function af_cluster() {
	# $Protein_name $AF_A3M $CLUSTER_OUTPUT_FOLDER
	python utils/AF-cluster.py $1 -i $2 -o $3
}

readarray -t array < leftover_singles.csv
for a in "${array[@]}"; do
	echo $a |
	 while IFS=, read -r ID_A ID_B DATE_A DATE_B LENGTH_A LENGTH_B TYPE OFFSET_A OFFSET_B
	 do
	 	CURR_PATH=$PARENT_PATH/${ID_A}${ID_B}
	 	echo The current path is: $CURR_PATH
	 	mkdir -p $CURR_PATH 
	 	FULL_A=$CURR_PATH/${ID_A}_full.pdb
	 	FULL_B=$CURR_PATH/${ID_B}_full.pdb
	 	FULL_A3M_A=$CURR_PATH/${ID_A}_full.a3m
	 	FULL_A3M_B=$CURR_PATH/${ID_B}_full.a3m
	 	CONF_A=$CURR_PATH/${ID_A}_conf.pdb
	 	CONF_B=$CURR_PATH/${ID_B}_conf.pdb
	 	CONF_R=$CURR_PATH/rfdiff_conf.pdb
	 	CONF_R_FOLDER=$CURR_PATH/rfdiff/rfdiff
	 	PROT_MPNN_A3M_A=$CURR_PATH/${ID_A}_conf.a3m
	 	PROT_MPNN_A3M_B=$CURR_PATH/${ID_B}_conf.a3m
	 	PROT_MPNN_A3M_R=$CURR_PATH/rfdiff_conf.a3m
	 	AF_A=$CURR_PATH/${ID_A}_AF.pdb
	 	AF_B=$CURR_PATH/${ID_B}_AF.pdb
	 	ALPHAFOLD_A=$CURR_PATH/${ID_A}_alphafold
	 	ALPHAFOLD_B=$CURR_PATH/${ID_B}_alphafold
	 	ALPHAFOLD_CONF_A=$CURR_PATH/${ID_A}_alphafold_conf
		mkdir -p $ALPHAFOLD_CONF_A
		AF2_SCORE_DIR_A=$CURR_PATH/${ID_A}_score_dir_af2_nat
		mkdir -p $AF2_SCORE_DIR_A
	 	ROSETTAFOLD_CONF_A=$CURR_PATH/${ID_A}_rosettafold_conf.pdb
	 	ALPHAFOLD_CONF_B=$CURR_PATH/${ID_B}_alphafold_conf
		mkdir -p $ALPHAFOLD_CONF_B
		AF2_SCORE_DIR_B=$CURR_PATH/${ID_B}_score_dir_af2_nat
	 	mkdir -p $AF2_SCORE_DIR_B
		ROSETTAFOLD_CONF_B=$CURR_PATH/${ID_B}_rosettafold_conf.pdb
	 	PROT_MPNN_A=$CURR_PATH/${ID_A}_prot
	 	PROT_MPNN_B=$CURR_PATH/${ID_B}_prot
	 	PROT_MPNN_CONF_A=$CURR_PATH/${ID_A}_prot_conf.pdb
	 	PROT_MPNN_CONF_B=$CURR_PATH/${ID_B}_prot_conf.pdb
	 	PROT_MPNN_CONF_A_RF=$CURR_PATH/${ID_A}_prot_conf_rf.pdb
	 	PROT_MPNN_CONF_B_RF=$CURR_PATH/${ID_B}_prot_conf_rf.pdb
	 	PROT_MPNN_R=$CURR_PATH/rfdiff_prot
	 	SCORE_REF_AB=$CURR_PATH/score_ref_AB.json
	 	SCORE_REF_AR=$CURR_PATH/score_ref_AR.json
	 	SCORE_REF_BR=$CURR_PATH/score_ref_BR.json
	 	SCORE_AA=$CURR_PATH/score_AA.json
	 	SCORE_AA_RF=$CURR_PATH/score_AA_RF.json
	 	SCORE_AB=$CURR_PATH/score_AB.json
	 	SCORE_BA=$CURR_PATH/score_BA.json
	 	SCORE_BB=$CURR_PATH/score_BB.json
	 	SCORE_BB_RF=$CURR_PATH/score_BB_RF.json
	 	SCORE_APROT=$CURR_PATH/score_Aprot.json
	 	SCORE_BPROT=$CURR_PATH/score_Bprot.json
	 	SCORE_RPROT=$CURR_PATH/score_Rprot.json
	 	SCORE_APROT_RF=$CURR_PATH/score_Aprot_rf.json
	 	SCORE_BPROT_RF=$CURR_PATH/score_Bprot_rf.json
	 	SCORE_RPROT_RF=$CURR_PATH/score_Rprot_rf.json
 		# All models
	 	PROT_MPNN_CONF_A_DIR=$CURR_PATH/${ID_A}_conf_dir
	 	mkdir -p $PROT_MPNN_CONF_A_DIR
	 	PROT_MPNN_CONF_B_DIR=$CURR_PATH/${ID_B}_conf_dir
	 	mkdir -p $PROT_MPNN_CONF_B_DIR
	 	PROT_MPNN_CONF_A_DIR_PACKED=$CURR_PATH/${ID_A}_conf_dir_packed
	 	mkdir -p $PROT_MPNN_CONF_A_DIR_PACKED
	 	PROT_MPNN_CONF_B_DIR_PACKED=$CURR_PATH/${ID_B}_conf_dir_packed
	 	mkdir -p $PROT_MPNN_CONF_B_DIR_PACKED
	 	PROT_MPNN_A_DIR_PACKED=$CURR_PATH/${ID_A}_dir_packed
	 	mkdir -p $PROT_MPNN_A_DIR_PACKED
	 	PROT_MPNN_B_DIR_PACKED=$CURR_PATH/${ID_B}_dir_packed
	 	mkdir -p $PROT_MPNN_B_DIR_PACKED
	 	SCORE_APROT_DIR=$CURR_PATH/${ID_A}_prot_dir
	 	mkdir -p $SCORE_APROT_DIR
	 	SCORE_BPROT_DIR=$CURR_PATH/${ID_B}_prot_dir
	 	mkdir -p $SCORE_BPROT_DIR
	 	SCORE_APROT_DIR_PACKED=$CURR_PATH/${ID_A}_prot_dir_packed
	 	mkdir -p $SCORE_APROT_DIR_PACKED
	 	SCORE_BPROT_DIR_PACKED=$CURR_PATH/${ID_B}_prot_dir_packed
	 	mkdir -p $SCORE_BPROT_DIR_PACKED
	 
	 	#AF3
	 	PROT_MPNN_CONF_A_DIR_AF3=$CURR_PATH/${ID_A}_conf_dir_af3
	 	mkdir -p $PROT_MPNN_CONF_A_DIR_AF3
	 	PROT_MPNN_CONF_B_DIR_AF3=$CURR_PATH/${ID_B}_conf_dir_af3
	 	mkdir -p $PROT_MPNN_CONF_B_DIR_AF3
	 	PROT_MPNN_A_DIR_AF3=$CURR_PATH/${ID_A}_prot_dir_af3
	 	mkdir -p $PROT_MPNN_A_DIR_AF3
	 	PROT_MPNN_B_DIR_AF3=$CURR_PATH/${ID_B}_prot_dir_af3
	 	mkdir -p $PROT_MPNN_B_DIR_AF3
	 	SCORE_APROT_DIR_AF3=$CURR_PATH/${ID_A}_prot_dir_af3
	 	mkdir -p $SCORE_APROT_DIR_AF3
	 	SCORE_BPROT_DIR_AF3=$CURR_PATH/${ID_B}_prot_dir_af3
	 	mkdir -p $SCORE_BPROT_DIR_AF3
	 
		DIFF_SEQ_A3M_A=$CURR_PATH/${ID_A}_diff.a3m
		DIFF_SEQ_A3M_B=$CURR_PATH/${ID_B}_diff.a3m
		DIFF_SEQ_ALPHAFOLD_A=$CURR_PATH/${ID_A}_diff_alphafold
		mkdir -p $DIFF_SEQ_ALPHAFOLD_A
		DIFF_SEQ_ALPHAFOLD_B=$CURR_PATH/${ID_B}_diff_alphafold
		mkdir -p $DIFF_SEQ_ALPHAFOLD_B
		PROT_MPNN_CONF_A_DIFF=$CURR_PATH/${ID_A}_prot_conf_diff/
		mkdir -p $PROT_MPNN_CONF_A_DIFF 
		PROT_MPNN_CONF_B_DIFF=$CURR_PATH/${ID_B}_prot_conf_diff/
		mkdir -p $PROT_MPNN_CONF_B_DIFF
		SCORE_A_DIFF=$CURR_PATH/score_${ID_A}_diff/
		mkdir -p $SCORE_A_DIFF
		SCORE_B_DIFF=$CURR_PATH/score_${ID_B}_diff/
		mkdir -p $SCORE_B_DIFF
		get_pdbs $FULL_A $AF_A $FULL_A3M_A || (echo "Problem with $FULL_A"; continue)
		get_pdbs $FULL_B $AF_B $FULL_A3M_B || (echo "Problem with $FULL_B"; continue)
		#Make targets
		cut_pdb $AF_A $CONF_A $OFFSET_A $LENGTH_A
		cut_pdb $AF_B $CONF_B $OFFSET_B $LENGTH_B
		#ProteinMPNN
		prot_MPNN $AF_A $PROT_MPNN_A3M_A
		prot_MPNN $AF_B $PROT_MPNN_A3M_B
		echo "Folding $ID_A"
		fold_alpha2 $FULL_A.fasta $ALPHAFOLD_A
		for file in ${ALPHAFOLD_A}/*.a3m ; do
			A3M_file_A=$file
		done
		fold_rosetta $A3M_file_A ${ALPHAFOLD_A}/"Rosetta.pdb"
		echo "Folding $ID_B"
		fold_alpha2 $FULL_B.fasta $ALPHAFOLD_B
		for file in ${ALPHAFOLD_B}/*.a3m ; do
			A3M_file_B=$file
		done
		fold_rosetta $A3M_file_B ${ALPHAFOLD_B}/"Rosetta.pdb"
	 	echo "Folding Prot $ID_A"
		fold_alpha2 $PROT_MPNN_A3M_A $PROT_MPNN_A
	 	fold_alpha3 $PROT_MPNN_A3M_A $PROT_MPNN_A_DIR_AF3
	 	fold_rosetta $PROT_MPNN_A3M_A ${PROT_MPNN_A}/Rosetta.pdb
	 	echo "Folding Prot $ID_B"
		fold_alpha2 $PROT_MPNN_A3M_B $PROT_MPNN_B
	 	fold_alpha3 $PROT_MPNN_A3M_B $PROT_MPNN_B_DIR_AF3
	 	fold_rosetta $PROT_MPNN_A3M_B ${PROT_MPNN_B}/Rosetta.pdb
		cut_pdb ${ALPHAFOLD_A}/Rosetta.pdb $ROSETTAFOLD_CONF_A $OFFSET_A $LENGTH_A
		cut_pdb ${ALPHAFOLD_B}/Rosetta.pdb $ROSETTAFOLD_CONF_B $OFFSET_B $LENGTH_B
		cut_pdb ${PROT_MPNN_A}/Rosetta.pdb $PROT_MPNN_CONF_A_RF $OFFSET_A $LENGTH_A
		cut_pdb ${PROT_MPNN_B}/Rosetta.pdb $PROT_MPNN_CONF_B_RF $OFFSET_B $LENGTH_B
		score $ROSETTAFOLD_CONF_A $FULL_A $SCORE_AA_RF
		score $ROSETTAFOLD_CONF_B $FULL_B $SCORE_BB_RF
		score $PROT_MPNN_CONF_A_RF $FULL_A $SCORE_APROT_RF
		score $PROT_MPNN_CONF_B_RF $FULL_B $SCORE_BPROT_RF
		ALIGNMENT_A=${FULL_A%????}_align.fasta
		pdb_tofasta $FULL_A >$FULL_A.fasta
		cat $FULL_A.fasta $FULL_A3M_A > temp_A.fasta
		clustalo -i temp_A.fasta --outfmt=a2m --wrap=99999 > $ALIGNMENT_A
		python utils/get_secstrucs.py --input_pdb $FULL_A --output_json ${FULL_A%????}_sec_struc.json --alignment_path $ALIGNMENT_A --cut_start $OFFSET_A --cut_length $LENGTH_A
		ALIGNMENT_B=${FULL_B%????}_align.fasta
		pdb_tofasta $FULL_B >$FULL_B.fasta
		cat $FULL_B.fasta $FULL_A3M_B > temp_B.fasta
		clustalo -i temp_B.fasta --outfmt=a2m --wrap=99999 > $ALIGNMENT_B
		python utils/get_secstrucs.py --input_pdb $FULL_B --output_json ${FULL_B%????}_sec_struc.json --alignment_path $ALIGNMENT_B --cut_start $OFFSET_B --cut_length $LENGTH_B
		python utils/get_secstrucs.py --input_pdb $PROT_MPNN_CONF_A_RF --output_json ${PROT_MPNN_CONF_A_RF%????}_sec_struc.json
		python utils/get_secstrucs.py --input_pdb $PROT_MPNN_CONF_B_RF --output_json ${PROT_MPNN_CONF_B_RF%????}_sec_struc.json
	 	# AF2, Natural models
		for FOUND_FILE in $ALPHAFOLD_A/*rank_00[1-5]*.pdb; do
			FOUND_FILENAME=$(basename $FOUND_FILE)
			cut_pdb $FOUND_FILE $ALPHAFOLD_CONF_A/$FOUND_FILENAME $OFFSET_A $LENGTH_A
			score $ALPHAFOLD_CONF_A/$FOUND_FILENAME $FULL_A $AF2_SCORE_DIR_A/${FOUND_FILENAME%????}_score.json
 			python utils/get_secstrucs.py --input_pdb $ALPHAFOLD_CONF_A/$FOUND_FILENAME --output_json $AF2_SCORE_DIR_A/${FOUND_FILENAME%????}_sec_struc.json
		done
		for FOUND_FILE in $ALPHAFOLD_B/*rank_00[1-5]*.pdb; do
			FOUND_FILENAME=$(basename $FOUND_FILE)
			cut_pdb $FOUND_FILE $ALPHAFOLD_CONF_B/$FOUND_FILENAME $OFFSET_B $LENGTH_B
			score $ALPHAFOLD_CONF_B/$FOUND_FILENAME $FULL_B $AF2_SCORE_DIR_B/${FOUND_FILENAME%????}_score.json
 			python utils/get_secstrucs.py --input_pdb $ALPHAFOLD_CONF_B/$FOUND_FILENAME --output_json $AF2_SCORE_DIR_B/${FOUND_FILENAME%????}_sec_struc.json
		done
	 	# AF2, ProteinMPNN models
		for FOUND_FILE in $PROT_MPNN_A/*rank_00[1-5]*.pdb; do
			FOUND_FILENAME=$(basename $FOUND_FILE)
			cut_pdb $FOUND_FILE $PROT_MPNN_CONF_A_DIR/$FOUND_FILENAME $OFFSET_A $LENGTH_A
			score $PROT_MPNN_CONF_A_DIR/$FOUND_FILENAME $FULL_A $SCORE_APROT_DIR/${FOUND_FILENAME%????}_score.json
 			python utils/get_secstrucs.py --input_pdb $PROT_MPNN_CONF_A_DIR/$FOUND_FILENAME --output_json $SCORE_APROT_DIR/${FOUND_FILENAME%????}_sec_struc.json
	 	done
		for FOUND_FILE in $PROT_MPNN_B/*rank_00[1-5]*.pdb; do
			FOUND_FILENAME=$(basename $FOUND_FILE)
			cut_pdb $FOUND_FILE $PROT_MPNN_CONF_B_DIR/$FOUND_FILENAME $OFFSET_B $LENGTH_B
			score $PROT_MPNN_CONF_B_DIR/$FOUND_FILENAME $FULL_B $SCORE_BPROT_DIR/${FOUND_FILENAME%????}_score.json
 			python utils/get_secstrucs.py --input_pdb $PROT_MPNN_CONF_B_DIR/$FOUND_FILENAME --output_json $SCORE_BPROT_DIR/${FOUND_FILENAME%????}_sec_struc.json
		done
		
		# Do all possible models
		declare -a numbers=("0" "1" "2" "3" "4") 
		for number in "${numbers[@]}"; do
			FOUND_FILE=$PROT_MPNN_A_DIR_AF3/sample-it-${number}.pdb
			FOUND_FILENAME=$(basename $FOUND_FILE)
			cut_pdb $FOUND_FILE $PROT_MPNN_CONF_A_DIR_AF3/$FOUND_FILENAME $OFFSET_A $LENGTH_A
			score $PROT_MPNN_CONF_A_DIR_AF3/$FOUND_FILENAME $FULL_A $SCORE_APROT_DIR_AF3/${FOUND_FILENAME%????}_score.json
 			mkdir -p ${SCORE_APROT_DIR_AF3}_sec_struc
			python utils/get_secstrucs.py --input_pdb $FOUND_FILE --output_json ${SCORE_APROT_DIR_AF3}_sec_struc/${FOUND_FILENAME%????}_sec_struc.json
		done
		for number in "${numbers[@]}"; do
			FOUND_FILE=$PROT_MPNN_B_DIR_AF3/sample-it-${number}.pdb
			FOUND_FILENAME=$(basename $FOUND_FILE)
			cut_pdb $FOUND_FILE $PROT_MPNN_CONF_B_DIR_AF3/$FOUND_FILENAME $OFFSET_B $LENGTH_B
			score $PROT_MPNN_CONF_B_DIR_AF3/$FOUND_FILENAME $FULL_B $SCORE_BPROT_DIR_AF3/${FOUND_FILENAME%????}_score.json
 			mkdir -p ${SCORE_BPROT_DIR_AF3}_sec_struc
			python utils/get_secstrucs.py --input_pdb $FOUND_FILE --output_json ${SCORE_BPROT_DIR_AF3}_sec_struc/${FOUND_FILENAME%????}_sec_struc.json
		done
		#Filter out similar sequences
		echo $ID_A
		python utils/filter_similar.py --input $PROT_MPNN_A3M_A --output $DIFF_SEQ_A3M_A --seqid $SEQ_IDENTITY 
		python utils/filter_similar.py --input $PROT_MPNN_A3M_B --output $DIFF_SEQ_A3M_B --seqid $SEQ_IDENTITY 
		fold_alpha2 $DIFF_SEQ_A3M_A $DIFF_SEQ_ALPHAFOLD_A
		fold_alpha2 $DIFF_SEQ_A3M_B $DIFF_SEQ_ALPHAFOLD_B
		# AF2, ProteinMPNN models
		for FOUND_FILE in $DIFF_SEQ_ALPHAFOLD_A/*rank_00[1-5]*.pdb; do
			FOUND_FILENAME=$(basename $FOUND_FILE)
			cut_pdb $FOUND_FILE $PROT_MPNN_CONF_A_DIFF/$FOUND_FILENAME $OFFSET_A $LENGTH_A
			score $PROT_MPNN_CONF_A_DIFF/$FOUND_FILENAME $FULL_A $SCORE_A_DIFF/${FOUND_FILENAME%????}_score.json
	 	done
		echo $ID_B
		for FOUND_FILE in $DIFF_SEQ_ALPHAFOLD_B/*rank_00[1-5]*.pdb; do
			FOUND_FILENAME=$(basename $FOUND_FILE)
			cut_pdb $FOUND_FILE $PROT_MPNN_CONF_B_DIFF/$FOUND_FILENAME $OFFSET_B $LENGTH_B
			score $PROT_MPNN_CONF_B_DIFF/$FOUND_FILENAME $FULL_B $SCORE_B_DIFF/${FOUND_FILENAME%????}_score.json
		done
		AF_CLUSTERS_A=$CURR_PATH/${ID_A}_clusters
		AF_CLUSTERS_B=$CURR_PATH/${ID_B}_clusters
		mkdir -p $AF_CLUSTERS_A
		mkdir -p $AF_CLUSTERS_B
		af_cluster ${ID_A} ${ALPHAFOLD_A}/1inp.a3m $AF_CLUSTERS_A
		af_cluster ${ID_B} ${ALPHAFOLD_B}/1inp.a3m $AF_CLUSTERS_B
		for CLUSTER_FILE in $AF_CLUSTERS_A/*.a3m; do
			CLUSTER_DIR=${CLUSTER_FILE%????}
			CUT_DIR=${CLUSTER_DIR}_cut
			mkdir -p $CLUSTER_DIR
			mkdir -p $CUT_DIR
			fold_alpha2 $CLUSTER_FILE $CLUSTER_DIR
			for FOUND_FILE in $CLUSTER_DIR/*rank_00[1-5]*.pdb; do
				FOUND_FILENAME=$(basename $FOUND_FILE)
				cut_pdb $FOUND_FILE $CUT_DIR/$FOUND_FILENAME $OFFSET_A $LENGTH_A
				score $CUT_DIR/$FOUND_FILENAME $FULL_A $CUT_DIR/${FOUND_FILENAME%????}_score_A.json
				score $CUT_DIR/$FOUND_FILENAME $FULL_B $CUT_DIR/${FOUND_FILENAME%????}_score_B.json
			done 
		done
		for CLUSTER_FILE in $AF_CLUSTERS_B/*.a3m; do
			CLUSTER_DIR=${CLUSTER_FILE%????}
			CUT_DIR=${CLUSTER_DIR}_cut
			mkdir -p $CLUSTER_DIR
			mkdir -p $CUT_DIR
			fold_alpha2 $CLUSTER_FILE $CLUSTER_DIR
			for FOUND_FILE in $CLUSTER_DIR/*rank_00[1-5]*.pdb; do
				FOUND_FILENAME=$(basename $FOUND_FILE)
				cut_pdb $FOUND_FILE $CUT_DIR/$FOUND_FILENAME $OFFSET_B $LENGTH_B
				score $CUT_DIR/$FOUND_FILENAME $FULL_A $CUT_DIR/${FOUND_FILENAME%????}_score_A.json
				score $CUT_DIR/$FOUND_FILENAME $FULL_B $CUT_DIR/${FOUND_FILENAME%????}_score_B.json
			done 
		done
		echo "Done with $CURR_PATH"
	done
done