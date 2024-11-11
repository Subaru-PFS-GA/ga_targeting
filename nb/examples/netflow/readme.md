# Netflow example

Run notebooks in this order:

1. `pmap.ipynb`

Generate the probability map from a CMD simulation.

2. `selection.ipynb`

Apply color-cuts, assign priorities, calculate required exposure time and generate the input target list from the observations.

3. `targets.ipynb` 

Plot the input target lists including the photometric data, the flux standards and the sky positions. This notebook loads the input target list generated in the previous notebook. File names are looked up in the configuration file.

4. `calibration.ipynb`

Plot the calibration targets (sky and flux standards) for each pointing, as well as the histogram of visible targets per cobra and visible targets per cobra group. These are usuful to detect problems when the number of available calibration targets is too low and results in an unfeasible netflow solution.

5. Run the command-line script to solve the netflow problem

To run netflow on well-known fields (i.e. those that have their own classes defined in the library, like dSphs), use the command-line script:

    $ ./ga-netflow --dsph fornax --config ./configs/netflow/dSph/fornax.py --out /datascope/subaru/user/dobos/netflow/run/fornax_1 --nvisit 1 --exp-time 10800 --time-limit 300

This script will copy the input target lists and intermediate results of the processing into the output directory. All subsequent notebooks require these files.

6. `solution.ipynb`

This a very preliminary notebook. I'll add more routines here to verify the netflow solution and identify potential problems such as the constraints causing infeasibility, very high final cost, etc.

7. `assignments.ipynb`

Plot the assigned targets and various statistics.

8. `cobra_groups.ipynb`

Generate plots of the calibration targets and the cobra groups. These plots are
useful to verify that the calibration targets are uniformly distributed in the field of view or along the spectrograph slits.

TODO: generate plots for each visit and each calibration target type