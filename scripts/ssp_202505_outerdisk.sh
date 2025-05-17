#!/bin/bash

set -e

PREFIX=SSP
VERSION="006"

NVISITS=6
NFRAMES=4
EXP_TIME=1800
OBS_TIME="2025-05-28T12:00:00"
OBS_RUN="2025-05"
PROPOSAL_ID="S25A-OT02"
INPUT_CATALOG_ID="10092"

# EXTRA_OPTIONS="--debug"
EXTRA_OPTIONS=""

for FIELD in l90_b28 l90_b29 l90_bm28 l90_bm29; do
# for FIELD in l90_b28; do
# for FIELD in l90_b29; do
# for FIELD in l90_bm28; do
# for FIELD in l90_bm29; do

    FIELD_DIR=$PFS_TARGETING_DATA/data/targeting/MW/outerdisk_${FIELD}_${PREFIX}
    IMPORT_DIR=${FIELD_DIR}/import/outerdisk_${FIELD}_${PREFIX}_${VERSION}
    NETFLOW_DIR=${FIELD_DIR}/netflow/outerdisk_${FIELD}_${NVISITS}_${PREFIX}_${VERSION}
    EXPORT_DIR=${FIELD_DIR}/export/outerdisk_${FIELD}_${NVISITS}_${PREFIX}_${VERSION}

    # rm -Rf "${IMPORT_DIR}"
    # rm -Rf "${NETFLOW_DIR}"
    rm -Rf "${EXPORT_DIR}"

    # if [ ! -d "$IMPORT_DIR" ]; then
    #     ga-import \
    #         --config \
    #             ./configs/netflow/${PREFIX}/MW/outerdisk_common.py \
    #             ./configs/netflow/SSP/MW/outerdisk_${FIELD}.py \
    #         --exp-time ${EXP_TIME} \
    #         --out "${IMPORT_DIR}" \
    #         ${EXTRA_OPTIONS}
    # fi

    # if [ ! -d "$NETFLOW_DIR" ]; then
    #     ga-netflow \
    #         --config \
    #             ./configs/netflow/SSP/MW/outerdisk_common.py \
    #             ./configs/netflow/SSP/MW/outerdisk_${FIELD}.py \
    #         --nvisits ${NVISITS} \
    #         --exp-time ${EXP_TIME} \
    #         --obs-time ${OBS_TIME} \
    #         --in "${IMPORT_DIR}" \
    #         --out "${NETFLOW_DIR}" \
    #         ${EXTRA_OPTIONS}
    # fi

    if [ ! -d "$EXPORT_DIR" ]; then
        ga-export \
            --in "${NETFLOW_DIR}" \
            --out "$EXPORT_DIR" \
            --input-catalog-id $INPUT_CATALOG_ID \
            --proposal-id "$PROPOSAL_ID" \
            --nframes $NFRAMES \
            --obs-run "$OBS_RUN" \
            ${EXTRA_OPTIONS}
    fi

done