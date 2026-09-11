#!/bin/bash

SEEDS=1
MODELS=5
RECYCLES=3
BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
SEED=6217

INPUT_A3M=$1
OUTPUT_DIR=$2
mkdir -p $OUTPUT_DIR
tail -n +4 $INPUT_A3M > ${INPUT_A3M%????}_cut.a3m
colabfold_batch --overwrite-existing-results --random-seed $SEED --num-seeds $SEEDS --num-models $MODELS --num-recycle $RECYCLES --num-relax 0 ${INPUT_A3M%????}_cut.a3m $OUTPUT_DIR
BEST=$(find $OUTPUT_DIR -name $BESTNAME | tail -1)
cp $BEST $OUTPUT_DIR/best.pdb
