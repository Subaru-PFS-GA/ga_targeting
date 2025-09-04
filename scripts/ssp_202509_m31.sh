#!/bin/bash

set -e

PREFIX=SSP
VERSION="001"

NVISITS=10
NFRAMES=2
EXP_TIME=18000
OBS_TIME="2025-09-14T09:00:00"      # GMT
OBS_RUN="2025-09"
PROPOSAL_ID="S25B-OT02"
INPUT_CATALOG_ID="10092"

# EXTRA_OPTIONS="--debug --skip-notebooks"
# EXTRA_OPTIONS="--skip-notebooks"

for SECTOR in 'm31_E0'; do

    SECTOR_DIR="$PFS_TARGETING_DATA/data/targeting/m31/${SECTOR}_${PREFIX}"
    PMAP_DIR="${SECTOR_DIR}/pmap/${SECTOR}_${PREFIX}_${VERSION}"
    SAMPLE_DIR="$PFS_TARGETING_DATA/data/targeting/m31/m31_all_${PREFIX}/${SECTOR_DIR}/sample/m31_all_${PREFIX}"
    IMPORT_DIR="${SECTOR_DIR}/import/${SECTOR}_${PREFIX}_${VERSION}"
    NETFLOW_DIR="${SECTOR_DIR}/netflow/${SECTOR}_${NVISITS}_${PREFIX}_${VERSION}"
    EXPORT_DIR="${SECTOR_DIR}/export/${SECTOR}_${NVISITS}_${PREFIX}_${VERSION}"

    # Line 106 of ${CONFIG_FILE} needs to be changed to point to the correct directory.
    # CONFIG_FILE="./configs/netflow/${PREFIX}/m31/${SECTOR}.py"    
    M31_CONFIG_FILE="./configs/netflow/${PREFIX}/m31/m31.py"
    SECTOR_CONFIG_FILE="./configs/netflow/${PREFIX}/m31/${SECTOR}.py"

    mkdir -p "${SECTOR_DIR}"

    # rm -Rf "${PMAP_DIR}"
    # if [ ! -d "$PMAP_DIR" ]; then
    #     ga-pmap \
    #         --m31 ${SECTOR} \
    #         --config ./configs/pmap/${PREFIX}/m31/m31.py \
    #         --out "${PMAP_DIR}" \
    #         ${EXTRA_OPTIONS}
    # fi

    # rm -Rf "${SAMPLE_DIR}"
    # if [ ! -d "$SAMPLE_DIR" ]; then
    #     ga-sample \
    #         --m31 ${SECTOR} \
    #         --config ./configs/sample/${PREFIX}/m31/m31.py \
    #         --out "${SAMPLE_DIR}" \
    #         ${EXTRA_OPTIONS}
    # fi

    # rm -Rf "${IMPORT_DIR}"
    # if [ ! -d "$IMPORT_DIR" ]; then
    #     ga-import \
    #         --m31 ${SECTOR} \
    #         --config ./configs/netflow/${PREFIX}/m31/_common.py ${M31_CONFIG_FILE} ${SECTOR_CONFIG_FILE} \
    #         --exp-time ${EXP_TIME} \
    #         --out "${IMPORT_DIR}" \
    #         ${EXTRA_OPTIONS}
    # fi

    # rm -Rf "${NETFLOW_DIR}"
    # if [ ! -d "$NETFLOW_DIR" ]; then
    #     ga-netflow \
    #         --m31 ${SECTOR} \
    #         --config ./configs/netflow/${PREFIX}/m31/_common.py ${M31_CONFIG_FILE} ${SECTOR_CONFIG_FILE} \
    #         --nvisits ${NVISITS} \
    #         --exp-time ${EXP_TIME} \
    #         --obs-time ${OBS_TIME} \
    #         --in "${IMPORT_DIR}" \
    #         --out "${NETFLOW_DIR}" \
    #         ${EXTRA_OPTIONS}
    # fi

    rm -Rf "${EXPORT_DIR}"
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