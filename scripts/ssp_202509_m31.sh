#!/bin/bash

set -e

PREFIX=SSP
VERSION="001"

NVISITS=1
NFRAMES=2
EXP_TIME=18000
OBS_TIME="2025-09-14T09:00:00"      # GMT
OBS_RUN="2025-09"
PROPOSAL_ID="S25B-OT02"
INPUT_CATALOG_ID="10092"

EXTRA_OPTIONS="--debug"
# EXTRA_OPTIONS=""

for SECTOR in 'PFS_1'; do

    SECTOR_DIR="$PFS_TARGETING_DATA/data/targeting/m31/m31_${SECTOR}_${PREFIX}"
    IMPORT_DIR="${SECTOR_DIR}/import/m31_${SECTOR}_${PREFIX}_${VERSION}"
    NETFLOW_DIR="${SECTOR_DIR}/netflow/m31_${SECTOR}_${NVISITS}_${PREFIX}_${VERSION}"
    EXPORT_DIR="${SECTOR_DIR}/export/m31_${SECTOR}_${NVISITS}_${PREFIX}_${VERSION}"

    # CONFIG_FILE="./configs/netflow/${PREFIX}/m31/${SECTOR}.py"
    CONFIG_FILE="./configs/netflow/${PREFIX}/m31/m31.py"

    mkdir -p "${SECTOR_DIR}"

    # rm -Rf "${IMPORT_DIR}"
    # rm -Rf "${NETFLOW_DIR}"
    rm -Rf "${EXPORT_DIR}"

    # if [ ! -d "$IMPORT_DIR" ]; then
    #     ga-import \
    #         --m31 ${SECTOR} \
    #         --config ./configs/netflow/${PREFIX}/m31/_common.py ${CONFIG_FILE} \
    #         --exp-time ${EXP_TIME} \
    #         --out "${IMPORT_DIR}" \
    #         ${EXTRA_OPTIONS}
    # fi

    # if [ ! -d "$NETFLOW_DIR" ]; then
    #     ga-netflow \
    #         --m31 ${SECTOR} \
    #         --config ./configs/netflow/${PREFIX}/m31/_common.py ${CONFIG_FILE} \
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