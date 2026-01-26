#!/bin/bash

set -e

PREFIX=SSP
VERSION="001"

NVISITS=6
NFRAMES=2
EXP_TIME=1800
OBS_TIME="2026-01-17T06:00:00"
OBS_RUN="2026-01"
PROPOSAL_ID="S25B-OT02"
INPUT_CATALOG_ID="10092"

# EXTRA_OPTIONS="--debug"
EXTRA_OPTIONS=""

# for FIELD in l180_bm39; do
# for FIELD in l181_bm60; do
# for FIELD in l180_bm31; do
for FIELD in l180_bm39 l181_bm60 l180_bm31; do

    FIELD_DIR=$PFS_TARGETING_DATA/data/targeting/MW/outerhalo_${FIELD}_${PREFIX}
    IMPORT_DIR=${FIELD_DIR}/import/outerhalo_${FIELD}_${PREFIX}_${VERSION}
    NETFLOW_DIR=${FIELD_DIR}/netflow/outerhalo_${FIELD}_${NVISITS}_${PREFIX}_${VERSION}
    EXPORT_DIR=${FIELD_DIR}/export/outerhalo_${FIELD}_${NVISITS}_${PREFIX}_${VERSION}

    # rm -Rf "${IMPORT_DIR}"
    # if [ ! -d "$IMPORT_DIR" ]; then
    #    ga-import \
    #        --config \
    #             ./configs/netflow/${PREFIX}/MW/outerhalo_common.py \
    #             ./configs/netflow/SSP/MW/outerhalo_${FIELD}.py \
    #         --exp-time ${EXP_TIME} \
    #         --out "${IMPORT_DIR}" \
    #         ${EXTRA_OPTIONS}
    # fi

    # rm -Rf "${NETFLOW_DIR}"
    # if [ ! -d "$NETFLOW_DIR" ]; then
    #     ga-netflow \
    #         --config \
    #             ./configs/netflow/SSP/MW/outerhalo_common.py \
    #             ./configs/netflow/SSP/MW/outerhalo_${FIELD}.py \
    #         --nvisits ${NVISITS} \
    #         --exp-time ${EXP_TIME} \
    #         --obs-time ${OBS_TIME} \
    #         --in "${IMPORT_DIR}" \
    #         --out "${NETFLOW_DIR}" \
    #         ${EXTRA_OPTIONS}
    # fi

    # rm -Rf "${EXPORT_DIR}"
    # if [ ! -d "$EXPORT_DIR" ]; then
    #     ga-export \
    #         --in "${NETFLOW_DIR}" \
    #         --out "$EXPORT_DIR" \
    #         --input-catalog-id $INPUT_CATALOG_ID \
    #         --proposal-id "$PROPOSAL_ID" \
    #         --nframes $NFRAMES \
    #         --obs-run "$OBS_RUN" \
    #         ${EXTRA_OPTIONS}
    # fi

    # cp -r "${EXPORT_DIR}/runs/${OBS_RUN}/targets/GA" "/home/dobos/project/Subaru-PFS/spt_ssp_observation/runs/${OBS_RUN}/targets/"

    cat "${EXPORT_DIR}/runs/${OBS_RUN}/targets/GA/ppcList.ecsv" \
        | sed '/^#/d' | sed '/^ppc_code/d' \
        >> "/home/dobos/project/Subaru-PFS/spt_ssp_observation/runs/${OBS_RUN}/targets/GA/ppcList.ecsv"

done