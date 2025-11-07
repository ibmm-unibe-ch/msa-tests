#!/bin/bash

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

function pack_sidechains(){
	micromamba run -n attnpacker python utils/attnpack.py --inputpdb $1 --outputpath ${2}_temp.pdb --numcores 8
	tail -n+2 ${2}_temp.pdb| head -n -1 > $2
}

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
	cp $PARENT_PATH/temp_path/temp.pdb $1
	$MAXIT_PATH -input temp_path/temp.pdb -output temp_path/1inp.cif -o 1
	wget -O $1.fasta https://www.rcsb.org/fasta/entry/${PDB_ID}
	python sel_chain.py --fastapath $1.fasta --chain $CHAIN --outputpath $3
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
	echo "CURR"
	echo $CURR
	cp $CURR $2
}

function carbon() {
	pdb_tofasta $1>$2
	micromamba run -n carbonara python carbonar.py --num_sequences 127 --imprint_ratio 0.5 $1 $2
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

function relax(){
    micromamba run -n relaxation pdbfixer $1 --output=${1}_relaxed.pdb --add-atoms=all --keep-heterogens=none
	micromamba run -n relaxation python amber_relax.py --input ${1}_relaxed.pdb --output ${1}_relaxed.pdb
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
	align $MODEL $REFERENCE
	ALIGNED_STATE=$?
	if [ "$ALIGNED_STATE" -eq "1" ]; then
		docker run --rm -v $(pwd):$(pwd) registry.scicore.unibas.ch/schwede/openstructure:latest compare-structures --model $MODEL --reference $REFERENCE --output $3 --lddt --local-lddt --bb-lddt --bb-local-lddt --tm-score --rigid-scores --lddt-no-stereochecks -rna
	else
		docker run --rm -v $(pwd):$(pwd) registry.scicore.unibas.ch/schwede/openstructure:latest compare-structures --model $MODEL --reference $REFERENCE --output $3 --lddt --local-lddt --bb-lddt --bb-local-lddt --tm-score --rigid-scores --lddt-no-stereochecks
	fi
}

function rf_diffusion() {
	pa_dir="$(dirname "$3")"
	rm -r $pa_dir
	micromamba run -n SE3nv /data/jgut/template-analysis/RFdiffusion/scripts/run_inference.py "contigmap.contigs=[${5}-${5}]" "contigmap.provide_seq=[${5}-${5}]" inference.output_prefix=$3 inference.input_pdb=$1 inference.num_designs=$RF_DESIGNS diffuser.partial_T=$PARTIAL_T
	python select_lowest.py --parentpath $3 --confa $1 --confb $2 --output ${3}_best.pdb
	pdb_tofasta $1 > ${3}_rec.fasta
	$FASPR_PATH -i ${3}_best.pdb -o ${3}_best_packed.pdb -s ${3}_rec.fasta
	cp ${3}_best_packed.pdb $4
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

function fold_rosetta() {
	PAR_PATH_RF="$(dirname "$1")"
	INPUT_PATH=$PAR_PATH_RF/CURR_ROSETTA_input.a3m
	first_chars=$(head -c 3 $1)
	if [[ $first_chars == \#* ]]; then tail -n +2 $1 > $INPUT_PATH; else cp $1 $INPUT_PATH; fi
	micromamba run -n RF2 python /home/jgut/tools/RoseTTAFold2/network/predict.py -inputs $INPUT_PATH -prefix $PAR_PATH_RF/CURR_ROSETTA -model /home/jgut/tools/RoseTTAFold2/network/RF2_jan24.pt -nseqs 128 -nseqs_full 1
	mv $PAR_PATH_RF/CURR_ROSETTA_00_pred.pdb $2
	rm $PAR_PATH_RF/CURR_ROSETTA_*
}

function clean_pdb() {
	pdb_delhetatm $1 | pdb_delinsertion | pdb_reres -1 | pdb_tidy | grep ^ATOM | grep -E "ALA|ARG|ASN|ASP|CYS|GLU|GLN|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|SEC|PYL|HCY" > $2
}

INPUT_STRUC=$1
OUTPUT_DIR=$2

mkdir -p $2
PROT_MPNN_A3M=${OUTPUT_DIR}/prot_mpnn.a3m
clean_pdb $INPUT_STRUC temp
prot_MPNN temp $PROT_MPNN_A3M
fold_alpha2 $PROT_MPNN_A3M $OUTPUT_DIR
