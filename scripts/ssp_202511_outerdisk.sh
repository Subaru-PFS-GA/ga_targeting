#!/bin/bash

set -e

PREFIX=SSP
VERSION="004"

NVISITS=6
NFRAMES=2
EXP_TIME=1800
OBS_TIME="2025-11-18T13:00:00"
OBS_RUN="2025-11"
PROPOSAL_ID="S25B-OT02"
INPUT_CATALOG_ID="10092"

# EXTRA_OPTIONS="--debug"
EXTRA_OPTIONS=""

# for FIELD in l180_b16; do
# for FIELD in l180_b28; do
for FIELD in l180_b16 l180_b17 l180_b18 l180_b20 l180_b21 l180_b22 l180_b24 l180_b25 l180_b27 l180_b28 l180_b29; do

    FIELD_DIR=$PFS_TARGETING_DATA/data/targeting/MW/outerdisk_${FIELD}_${PREFIX}
    IMPORT_DIR=${FIELD_DIR}/import/outerdisk_${FIELD}_${PREFIX}_${VERSION}
    NETFLOW_DIR=${FIELD_DIR}/netflow/outerdisk_${FIELD}_${NVISITS}_${PREFIX}_${VERSION}
    EXPORT_DIR=${FIELD_DIR}/export/outerdisk_${FIELD}_${NVISITS}_${PREFIX}_${VERSION}

    # rm -Rf "${IMPORT_DIR}"
    # if [ ! -d "$IMPORT_DIR" ]; then
    #     ga-import \
    #         --config \
    #             ./configs/netflow/${PREFIX}/MW/outerdisk_common.py \
    #             ./configs/netflow/SSP/MW/outerdisk_${FIELD}.py \
    #         --exp-time ${EXP_TIME} \
    #         --out "${IMPORT_DIR}" \
    #         ${EXTRA_OPTIONS}
    # fi

    # rm -Rf "${NETFLOW_DIR}"
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