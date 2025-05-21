#!/bin/bash

set -e

PREFIX=SSP
VERSION="004"

NVISITS=1
NFRAMES=2
EXP_TIME=300
OBS_TIME="2025-06-01T15:00:00"          # UTC
OBS_RUN="2025-05"
PROPOSAL_ID="S25A-OT02"
INPUT_CATALOG_ID="10092"

SSP_OBS_PATH="/home/dobos/project/Subaru-PFS/spt_ssp_observation"

# EXTRA_OPTIONS="--debug"
EXTRA_OPTIONS=""

FIELDS="ra288_dec-11 ra288_dec-17 ra288_dec-22 ra336_dec-12"
# FIELDS="ra288_dec-22"

# for FIELD in $FIELDS; do

#     FIELD_DIR=$PFS_TARGETING_DATA/data/targeting/CC/crosscalib_${FIELD}
#     IMPORT_DIR=${FIELD_DIR}/import/crosscalib_${FIELD}_${PREFIX}_${VERSION}
#     NETFLOW_DIR=${FIELD_DIR}/netflow/crosscalib_${FIELD}_${NVISITS}_${PREFIX}_${VERSION}
#     EXPORT_DIR=${FIELD_DIR}/export/crosscalib_${FIELD}_${NVISITS}_${PREFIX}_${VERSION}

#     rm -Rf "${IMPORT_DIR}"
#     rm -Rf "${NETFLOW_DIR}"
#     rm -Rf "${EXPORT_DIR}"

#     if [ ! -d "$IMPORT_DIR" ]; then
#         ga-import \
#             --config \
#                 ./configs/netflow/${PREFIX}/CC/crosscalib_common.py \
#                 ./configs/netflow/SSP/CC/crosscalib_${FIELD}.py \
#             --exp-time ${EXP_TIME} \
#             --out "${IMPORT_DIR}" \
#             ${EXTRA_OPTIONS}
#     fi

#     if [ ! -d "$NETFLOW_DIR" ]; then
#         ga-netflow \
#             --config \
#                 ./configs/netflow/SSP/CC/crosscalib_common.py \
#                 ./configs/netflow/SSP/CC/crosscalib_${FIELD}.py \
#             --nvisits ${NVISITS} \
#             --exp-time ${EXP_TIME} \
#             --obs-time ${OBS_TIME} \
#             --in "${IMPORT_DIR}" \
#             --out "${NETFLOW_DIR}" \
#             ${EXTRA_OPTIONS}
#      fi

#     if [ ! -d "$EXPORT_DIR" ]; then
#         ga-export \
#             --in "${NETFLOW_DIR}" \
#             --out "$EXPORT_DIR" \
#             --input-catalog-id $INPUT_CATALOG_ID \
#             --proposal-id "$PROPOSAL_ID" \
#             --nframes $NFRAMES \
#             --obs-run "$OBS_RUN" \
#             ${EXTRA_OPTIONS}
#     fi

# done


# Clean up spt_ssp_observation

# rm ${SSP_OBS_PATH}/runs/${OBS_RUN}/pfs_designs/GA/*.fits
# rm ${SSP_OBS_PATH}/runs/${OBS_RUN}/targets/GA/fluxstd/*.ecsv
# rm ${SSP_OBS_PATH}/runs/${OBS_RUN}/targets/GA/science/*.ecsv
# rm ${SSP_OBS_PATH}/runs/${OBS_RUN}/targets/GA/sky/*.ecsv
# rm ${SSP_OBS_PATH}/runs/${OBS_RUN}/targets/GA/ppcList.ecsv


# Copy target lists, design files and concatenate the ppcList.ecsv file
for FIELD in $FIELDS; do

    FIELD_DIR=$PFS_TARGETING_DATA/data/targeting/CC/crosscalib_${FIELD}
    EXPORT_DIR=${FIELD_DIR}/export/crosscalib_${FIELD}_${NVISITS}_${PREFIX}_${VERSION}

    # Copy from the export dir to the spt_ssp_observation dir
    echo "Copying ${FIELD} files to ${SSP_OBS_PATH}"
    cp -R ${EXPORT_DIR}/runs ${SSP_OBS_PATH}

    # Concatenate the ppcList.ecsv files
    if [ ! -f ${EXPORT_DIR}/runs/${OBS_RUN}/targets/GA/ppcList.ecsv ]; then
        echo "File ${EXPORT_DIR}/runs/${OBS_RUN}/targets/GA/ppcList.ecsv does not exist"
    else
        echo "Concatenating ${EXPORT_DIR}/runs/${OBS_RUN}/targets/GA/ppcList.ecsv to ${SSP_OBS_PATH}/runs/${OBS_RUN}/targets/GA/ppcList.ecsv"
        cat ${EXPORT_DIR}/runs/${OBS_RUN}/targets/GA/ppcList.ecsv >> ${SSP_OBS_PATH}/runs/${OBS_RUN}/targets/GA/ppcList.ecsv
    fi
done