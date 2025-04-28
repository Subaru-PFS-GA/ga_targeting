# 5 Executing netflow

These steps assume that the configuration file is already created and all input target lists are available. For dSphs and M31, see Chapter 2 on how to generate the input target lists by selecting targets and assigning priorities.

## Import step

    $ ga-import --dsph ursaminor --config ./configs/netflow/SSP/dSph/_common.py ./configs/netflow/SSP/dSph/ursaminor.py --exp-time 1800 --out /datascope/subaru/data/targeting/dSph/ursaminor/import/SSP/ursaminor

## Netflow steps

Netflow can be executed in a single run or in stages. When executed in stages, fiber assignments are calculated for a subset of pointings. This is done in order to ensure early completion of the observations, see Chapter 1 and 4.

When executing netflow for all pointings, without stages, the command-line will be similar to the following.

    $ ga-netflow --dsph umi --config ./configs/netflow/SSP/dSph/_common.py ./configs/netflow/SSP/dSph/ursaminor.py --nvisits 6 --exp-time 1800 --obs-time 2025-03-25T10:00:00 --in /datascope/subaru/data/targeting/dSph/ursaminor/import/SSP/ursaminor --out /datascope/subaru/data/targeting/dSph/ursaminor/netflow/SSP/ursaminor_6_all

When executing netflow in stages, the first command (stage 0) must specify the output of `ga-import` as the input directory. Note the `--stage 0` argument.

    $ ga-netflow --dsph ursaminor --config ./configs/netflow/SSP/dSph/_common.py ./configs/netflow/SSP/dSph/ursaminor.py --stage 0 --nvisits 6 --exp-time 1800 --obs-time 2025-03-25T10:00:00 --in /datascope/subaru/data/targeting/dSph/ursaminor/import/SSP/ursaminor --out /datascope/subaru/data/targeting/dSph/ursaminor/netflow/SSP/ursaminor_6_0

All subsequent steps must specify the output directory of the previous step as the input. Note the different `--in` argument and `--out` directory.

    $ ga-netflow --dsph ursaminor --config ./configs/netflow/SSP/dSph/_common.py ./configs/netflow/SSP/dSph/ursaminor.py --stage 1 --nvisits 6 --exp-time 1800 --obs-time 2025-03-25T10:00:00 --in /datascope/subaru/data/targeting/dSph/ursaminor/netflow/SSP/ursaminor_6_0 --out /datascope/subaru/data/targeting/dSph/ursaminor/netflow/SSP/ursaminor_6_1

Repeat the above as many stages are defined.

## Export step

The export step converts the output of the netflow step into the format required by the observatory. The output of `ga-export` can directly be copied into the repository `spt_ssp_observation`.

The command-line of the export step will be very similar to the following:

    $ ga-export --in /datascope/subaru/data/targeting/dSph/ursaminor/netflow/SSP/ursaminor_6_0 /datascope/subaru/data/targeting/dSph/ursaminor/netflow/SSP/ursaminor_6_1 /datascope/subaru/data/targeting/dSph/ursaminor/netflow/SSP/ursaminor_6_2 /datascope/subaru/data/targeting/dSph/ursaminor/netflow/SSP/ursaminor_6_3 --out /datascope/subaru/data/targeting/dSph/ursaminor/export/SSP/ursaminor --input-catalog-id 10092 --proposal-id S25A-OT02 --nframes 4

