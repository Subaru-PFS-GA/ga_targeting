# 2 Sample compilation

This chapter applies to dSph and M31 field only. For these field, the selection of potential targets is implemented as part of the targeting tool. For other fields, selections by color, magnitude, etc. must be implemented elsewhere and the input target lists are to be provided instead of the raw catalogs.

This chapter describes how to create membership probability maps from CMD simulations and how to use them to assign priorities to stars within dSphs.

## Generating pmaps

Membership probability maps are generated from simulated CMDs. First, the color and magnitude selections are applied to the simulation, then a histogram of stellar counts is calculated for all simulated stellar populations (dSph + foreground). Finally, the histograms are normalized in each histogram cell to obtain the membership probability for each population.

The input to `ga-pmap` is a configuration file that contains all settings of the pmap generation, including input file locations. The config files are available under `./configs/pmap/`.

An example command of running `ga-pmap` is the following.

    $ ga-pmap --dsph ursaminor --config ./configs/pmap/SSP/dSph/ursaminor.py --out /datascope/subaru/data/targeting/dSph/ursaminor/pmap/ursaminor_nb

The argument `--dsph` selects the field we want to process. This tells the script which field configuration to be used. The color and magnitude selections are implemented in code under the namespace `pfs.ga.targeting.targets.dsph`. For M31 field, the argument `--m31` can be used and the selections are implemented under `pfs.ga.targeting.targets.m31`.

The `--out` sets the output directory. The directory must not exist and if it does an error will be raised to avoid accidental overwriting the output of previous runs.

### The pmap configuration file

The pmap configuration files are saved in `./configs/pmap` in the form of Python files defining a single dictionary stored in the variable `config`. The dictionary keys are the following:

* *cut_nb*: Set to `True` to make cuts based on the narrow-band magnitude.
* *keep_blue*: When `True`, ignore that color cut that excludes stars bluer than the blue edge of the halo.
* *extents*: List of lists, four elements all together, that define the extents of the probability map in terms of the broadband color HSC g - i and HSC g magnitude.
* *bins*: Number of histogram bins in color and magnitude.
* *population_names*
* *population_weights*: The prior weight of each population _before_ applying the selections.
* *merge_list*: List of lists that define which populations to be merged as 'foreground' and as 'members'.
* *sim_path*: Location of the simulated CMDs from CMDFit.
* *obs_path*: Path to the observed data. Used to generate the evaluation plots only.

## Generating the sample

The script `ga-sample` generates an input target list from a full catalog of observations by applying the selection criteria defined in code under the namespace `pfs.ga.targeting.targets.dsph` and `pfs.ga.targeting.targets.m31`. It also assigns the target priorities based on pmaps.

An example command-line of running `ga-sample` is the following:

    $ ga-sample --dsph ursaminor --config ./configs/sample/SSP/dSph/ursaminor.py --out /datascope/subaru/data/targeting/dSph/ursaminor/sample/SSP/ursaminor_nb

The script generates a file named `hsc_{target}_priorities.feather` which contains the input target lists for netflow including target priorities. It might also generate the file `gaia_{target}.feather` which contains the GAIA catalog matches to the HSC data.

### The sample configuration file

The sample configuration files are saved in `./configs/sample` in the form of Python files defining a single dictionary stored in the variable `config`. The dictionary keys are the following:

* *obs_path*: Path to the HSC catalog file.
* *pmap_path*: Path to the pmap file generated using `ga-pmap`
* *isochrones_path*: Path to the isochrones file.
* *cut_nb*: Set to `True` to make cuts based on the narrow-band magnitude. Must be consistent with the pmap.
* *keep_blue*: When `True`, ignore that color cut that excludes stars bluer than the blue edge of the halo. Must be consitent with the pmap.
* *lp_member_limit*: Probability limit (in log p units) at which stars should be considered members.
* *gaia_crossmatch*: Cross-match HSC data with GAIA.
* *gaia_crossmatch_radius*: Cross-match radius in arc sec.