# 4 Netflow configuration reference

TBW

## 4.1 Target lists

### Target list settings

* *prefix*: Must be one of `sci`, `cal` or `sky`. It defines the kind of targets which is the same of all targets within the list.
* *frame*: Reference frame of the coordinates. Defaults to `ICRS` but supports anything that astropy understands.
* *epoch*: Overrides the epoch of the coordinates. By default, the epoch is taken from the `epoch` column of the input file. If neither the column nor the `epoch` configuration option is found, the default epoch of J2016.0 is assumed.
* *catid*: Overrides the `catid` column of the target list. The value is provided by the observatory.

### Column map

Allows for mapping between the column names of the input target list and the column names required by the targeting tools. It must be a dictionary with the column names of the file as keys and the internal column names as values. This parameter is passed to Pandas's `DataFrame.rename` function.

### Value maps

Allow mapping values from those in the input target list to anything else. It is a dictionary of dictionaries. The keys on the first level specify the columns whereas keys on the second level define the input values and values on the second level define target values. This feature can be used to map values without modifying the input target lists. This configuration parameter is passed to Panda's `DataFrame.replace` function.

### Generated extra columns

A list of columns can be generated on the fly from configuration, without materializing them in the input target lists files. This is convenient to add columns required by the targeting tools but are not available as part of the target lists. For example, the ˙obcode` column can be easily generated this way.

Extra columns can have a constant value, can specify a string pattern that references other, already existing columns or can specify a lambda function along with its arguments.

TBW

## 4.2 Photometry

### Specifying per-filter magnitudes and fluxes

### Specifying per-band magnitudes and fluxes

### Photometric limits

## 4.3 Defining pointings

## 4.4 Defining stages

TBW: Write about the requirement of early completion and what this means for targeting: cannot split visits between pointings because they'd be observed months apart.

## 4.5 Netflow options

TBW

## 4.6 Gurobi options

TBW

These options are passed to the Gurobi optimizer transparently.

## 4.7 Special configurations for dwarf spheroidal galaxies

You need to update the following variables and functions:

* add new ID prefix in `python/pfs/ga/targeting/targets/ids.py`, make sure this is the one used in the python class for the object
* `pointings` in `__init__`
* `get_pmap_config` to change the color and magnitude limits
* `get_text_observation_reader` if the data file format is different from the default
* `get_filter_map`, to match the i filter (old or new) used for HSC photometry
* `get_nb_selection_mask`
* `get_selection_mask`
* `assign_priorities`