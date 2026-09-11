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

INPUT_A3M=$1
OUTPUT_DIR=$2
mkdir -p $OUTPUT_DIR
tail -n +4 $INPUT_A3M > ${INPUT_A3M%????}_cut.a3m
colabfold_batch --overwrite-existing-results --random-seed $SEED --num-seeds $SEEDS --num-models $MODELS --num-recycle $RECYCLES --num-relax 0 ${INPUT_A3M%????}_cut.a3m $OUTPUT_DIR
BEST=$(find $OUTPUT_DIR -name $BESTNAME | tail -1)
cp $BEST $OUTPUT_DIR/best.pdb