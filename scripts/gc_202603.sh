set -e

VERSION="001"
IMF="chab"

# for gc in NGC1904_M79; do
# for gc in NGC6341_M92; do
# for gc in NGC7078_M15; do
# for gc in NGC7089_M2; do
# for gc in NGC7099_M30; do
for gc in NGC1904_M79 NGC6341_M92 NGC7078_M15 NGC7089_M2 NGC7099_M30; do
    # Take name before the underscore to match the GC class names
    NGCNAME="${gc%%_*}"
    MNAME="${gc##*_}"
    SIMDIR="/datascope/subaru/data/cmdfit/run/$gc/sim/bin_${IMF}_250k_$VERSION"
    OUTDIR="/datascope/subaru/data/cmdfit/run/$gc/pmap/bin_${IMF}_250k_$VERSION"

    # # Create the pmap
    # rm -rf "$OUTDIR"
    # if [ ! -d "$OUTDIR" ]; then
    #     ga-pmap --gc $NGCNAME \
    #         --config ./configs/pmap/PI/GC/$NGCNAME.py \
    #         --sim-path $SIMDIR \
    #         --obs-path /datascope/subaru/data/targeting/GC/${NGCNAME}${MNAME}_07deg_ext.csv \
    #         --out $OUTDIR
    # fi

    # Compress contents of output directory into a zip file and copy to .
    ZIPFILE="${gc}_pmap.zip"
    zip -r "$ZIPFILE" -j "$OUTDIR"
done