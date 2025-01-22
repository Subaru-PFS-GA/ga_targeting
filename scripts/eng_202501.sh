# Process design files for the engineering run starting in January 2025

# Outer disk fields

function rmdirs() {
    out=$1

    if [ -d $out ]; then
        echo Removing $out
        rm -r $out
    fi
}

function run_netflow() {
    extra="$1"
    config="$2"
    out="$3"
    obs_time="$4"

    # echo $config
    # echo $out

    # Delete the output directory
    # if [ -d $out ]; then
    #     rm -r $out
    # fi

    # If the output directory doesn't exit, run netflow
    if [ ! -d $out ]; then
        ga-netflow \
            $extra \
            --config $config \
            --out $out \
            --obs-time $obs_time \
            --nvisits 1 \
            --no-nan-values --debug
    fi

    # If the output directory exists, zip up its contents
    if [ -d $out ]; then

        pushd $out > /dev/null

        # Get name of current directory
        dirname=${PWD##*/}

        # Zip up the output but exclude large files
        zip -r $dirname.zip * -x *_targets_*.feather *_summary.feather

        popd > /dev/null  
    fi
}

VERSION="003"

# Outer disk fields

# # for field in l180_b22 l180_b25; do
# for field in l180_b22; do
#     # for br in bright faint; do
#     for br in faint; do

#         extra=""
#         configs="./configs/netflow/ENG/MW/outerdisk_common.py ./configs/netflow/ENG/MW/outerdisk_${field}_${br}.py"
#         outdir="$PFS_TARGETING_DATA/data/targeting/MW/outerdisk_${field}_ENG/netflow/outerdisk_${field}_${br}_$VERSION"
#         obs_time="2025-01-24T10:00:00"
        
#         # rmdirs "$outdir"
#         run_netflow "$extra" "$configs" "$outdir" "$obs_time"
#     done
# done

##########################################################################################

# M31

configs="./configs/netflow/ENG/M31/m31.py"
outdir="$PFS_TARGETING_DATA/data/targeting/M31/M31_ENG/netflow/M31_ENG_$VERSION"
obs_time="2025-01-25T06:00:00"

# # rmdirs "$outdir"
run_netflow "--m31 m31" "$configs" "$outdir" "$obs_time"

##########################################################################################

# Bootes I

# configs="./configs/netflow/ENG/dSph/bootes.py"
# outdir="$PFS_TARGETING_DATA/data/targeting/dSph/bootes_ENG/netflow/bootes_ENG_$VERSION"
# obs_time="2025-01-24T14:00:00"

# rmdirs "$outdir"
# run_netflow "--dsph booi" "$configs" "$outdir" "$obs_time"