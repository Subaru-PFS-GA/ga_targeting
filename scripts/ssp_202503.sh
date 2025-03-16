#!/bin/bash

set -e

PREFIX=SSP
VERSION=001

FIELD=draco
STAGES="0 1"

# FIELD=ursaminor
# STAGES="0 1 2 3"

NVISITS=6
NFRAMES=4
EXP_TIME=1800
OBS_TIME="2025-03-25T10:00:00"
PROPOSAL_ID=S25A-OT02

EXTRA_OPTIONS="--skip-notebooks"

FIELD_DIR=$PFS_TARGETING_DATA/data/targeting/dSph/${FIELD}

PMAP_DIR=${FIELD_DIR}/pmap/${FIELD}_nb
SAMPLE_DIR=${FIELD_DIR}/sample/$PREFIX/${FIELD}_${VERSION}
IMPORT_DIR=${FIELD_DIR}/import/$PREFIX/${FIELD}_${VERSION}
EXPORT_DIR=${FIELD_DIR}/export/$PREFIX/${FIELD}_${VERSION}

if [ ! -d "$PMAP_DIR" ]; then
    ga-pmap --dsph ${FIELD} \
        --config ./configs/pmap/SSP/dSph/${FIELD}.py \
        --out $PMAP_DIR \
        ${EXTRA_OPTIONS}
fi

if [ ! -d "$SAMPLE_DIR" ]; then
    ga-sample --dsph ${FIELD} \
        --config ./configs/sample/SSP/dSph/${FIELD}.py \
        --out $SAMPLE_DIR \
        ${EXTRA_OPTIONS}
fi

# NOTE: input filies to ga-import are listed in the config file

if [ ! -d "$IMPORT_DIR" ]; then
    ga-import --dsph ${FIELD} \
        --config ./configs/netflow/SSP/dSph/_common.py ./configs/netflow/SSP/dSph/${FIELD}.py \
        --exp-time 1800 \
        --out $IMPORT_DIR \
        ${EXTRA_OPTIONS}
fi

indir=$IMPORT_DIR
for stage in $STAGES; do
    outdir=${FIELD_DIR}/netflow/${PREFIX}/${FIELD}_${NVISITS}_${stage}_${VERSION}
    if [ ! -d "$outdir" ]; then
        ga-netflow --dsph ${FIELD} \
            --config ./configs/netflow/SSP/dSph/_common.py ./configs/netflow/SSP/dSph/${FIELD}.py \
            --stage ${stage} \
            --nvisits ${NVISITS} \
            --exp-time ${EXP_TIME} \
            --obs-time "${OBS_TIME}" \
            --in $indir \
            --out $outdir \
            ${EXTRA_OPTIONS}
    fi
    indir=$outdir
done

if [ ! -d "$EXPORT_DIR" ]; then
    ga-export --in ${FIELD_DIR}/netflow/${PREFIX}/${FIELD}_${NVISITS}_?_${VERSION} \
        --out ${EXPORT_DIR} \
        --input-catalog-id 10092 \
        --proposal-id $PROPOSAL_ID \
        --nframes $NFRAMES \
        ${EXTRA_OPTIONS}
fi