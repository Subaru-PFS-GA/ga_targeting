#!/bin/bash

set -e

PREFIX=SSP
VERSION=001

FIELD=bootes
STAGES="0" # 1 2 3"        # TODO: update stages as needed, depending on number of pointings

NVISITS=6
NFRAMES=2
EXP_TIME=1800
OBS_TIME="2026-03-16T12:00:00"          # UTC, 2am in Hawaii
OBS_RUN="2026-03"
PROPOSAL_ID="S25B-OT02"
INPUT_CATALOG_ID="10092"

# EXTRA_OPTIONS="--skip-notebooks --debug"
EXTRA_OPTIONS=""

FIELD_DIR=$PFS_TARGETING_DATA/data/targeting/dSph/${FIELD}

PMAP_DIR=${FIELD_DIR}/pmap/${PREFIX}/${FIELD}_${PREFIX}_${VERSION}
SAMPLE_DIR=${FIELD_DIR}/sample/${PREFIX}/${FIELD}_${PREFIX}_${VERSION}
IMPORT_DIR=${FIELD_DIR}/import/${PREFIX}/${FIELD}_${PREFIX}_${VERSION}
EXPORT_DIR=${FIELD_DIR}/export/${PREFIX}/${FIELD}_${PREFIX}_${VERSION}

# NOTE: the NB cut for generating the pmap is different from the NB
#       cut used for sampling from the observations. This is because
#       the NB model magnitudes are off. Refer to 
#       python/pfs/ga/targeting/targets/dsph/sextans.py:189

# rm -r "$PMAP_DIR"
# if [ ! -d "$PMAP_DIR" ]; then
#     ga-pmap --dsph ${FIELD} \
#         --config ./configs/pmap/${PREFIX}/dSph/${FIELD}.py \
#         --out $PMAP_DIR \
#         ${EXTRA_OPTIONS}
# fi

# NOTE: pmap file is listed in the config file
#       make sure the file path matches the version number!
#       also make sure the NB cut is the correct one, use the
#       revised NB cut for Sextans here.

# rm -r "$SAMPLE_DIR"
# if [ ! -d "$SAMPLE_DIR" ]; then
#     ga-sample --dsph ${FIELD} \
#         --config ./configs/sample/${PREFIX}/dSph/${FIELD}.py \
#         --out $SAMPLE_DIR \
#         --obs-time "${OBS_TIME}" \
#         ${EXTRA_OPTIONS}
# fi

# NOTE: input files to ga-import are listed in the netflow config file
#       make sure the file path matches the version number!

rm -r "$IMPORT_DIR"
if [ ! -d "$IMPORT_DIR" ]; then
    ga-import --dsph ${FIELD} \
        --config \
            ./configs/netflow/${PREFIX}/dSph/_common.py \
            ./configs/netflow/${PREFIX}/dSph/${FIELD}.py \
        --exp-time ${EXP_TIME} \
        --out ${IMPORT_DIR} \
        ${EXTRA_OPTIONS}
fi

# indir=$IMPORT_DIR
# for stage in $STAGES; do
#     outdir=${FIELD_DIR}/netflow/${PREFIX}/${FIELD}_${NVISITS}_${stage}_${VERSION}
#     if [ ! -d "$outdir" ]; then
#         ga-netflow --dsph ${FIELD} \
#             --config \
#                 ./configs/netflow/${PREFIX}/dSph/_common.py \
#                 ./configs/netflow/${PREFIX}/dSph/${FIELD}.py \
#             --stage ${stage} \
#             --nvisits ${NVISITS} \
#             --exp-time ${EXP_TIME} \
#             --obs-time "${OBS_TIME}" \
#             --in $indir \
#             --out $outdir \
#             ${EXTRA_OPTIONS}
#     fi
#     indir=$outdir
# done

# rm -r "$EXPORT_DIR"
# if [ ! -d "$EXPORT_DIR" ]; then
#     ga-export \
#         --in ${FIELD_DIR}/netflow/${PREFIX}/${FIELD}_${NVISITS}_?_${VERSION} \
#         --out ${EXPORT_DIR} \
#         --input-catalog-id $INPUT_CATALOG_ID \
#         --proposal-id $PROPOSAL_ID \
#         --nframes $NFRAMES \
#         --obs-run "$OBS_RUN" \
#         ${EXTRA_OPTIONS}
# fi


# cp -r "${EXPORT_DIR}/runs/${OBS_RUN}/targets/GA" "/home/dobos/project/Subaru-PFS/spt_ssp_observation/runs/${OBS_RUN}/targets/"

# cat "${EXPORT_DIR}/runs/${OBS_RUN}/targets/GA/ppcList.ecsv" \
#     | sed '/^#/d' | sed '/^ppc_code/d' \
#     >> "/home/dobos/project/Subaru-PFS/spt_ssp_observation/runs/${OBS_RUN}/targets/GA/ppcList.ecsv"