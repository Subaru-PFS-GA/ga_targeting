# 3 Input target lists

Science targets, as well as flux standards and sky positions must be provided as input target lists. The configuration file offers a way of specifying the kind of targets, such as science or calibration, see Chapter 4.

## Target list preparation

Input targets lists must contain all potential targets and no stars that should not be targeted. To generate input target lists for dSph and M31, use the corresponding command-line tool `ga-sample`. For all other fields, generate the target list with custom code.

Do not cross-match catalogs but provide them as they are. The targeting tools will cross-match input target lists, so it is not necessary to deal with cross-matching when creating the input target lists, unless selection based on multiple catalogs is necessary. When a match is found, collisions in parameters are resolved as follows:

* always the highest priority is used for a star
* always the longest exposure time is used for a star
* coordinates, etc. are taken as found in the first input target list the star appears in

Provide `priority` as an integer number. Smaller numbers generally mean higher priority, but priorities and their associated netflow costs are defined in a configuration file which can be freely modified, see Chapter 4.

## Target list file formats

The targeting tools currently support the file formats `csv` and `feather` but can easily extended to any format that Pandas supports. In addition, `fits` files with a single BINTABLE can also be provided. Data types are automatically converted to the data type required by netflow and the PfsDesign files.

## Target list columns

TBW

Column names can be anything because name mapping can be defined in the configuration file. Also, some mandatory columns can be created based on a definition in the configuration file. Here we list the name of the columns that the targeting tools assume. If column names are different, they must be mapped to these names in the configuration. 

Column names are case-sensitive!

Mandatory target list columns:

* *targetid*: Unique identifier of the target. It must correspond to some identifier in the original source catalog. It must be unique within the input target list but not unique overall. The targeting tools internally use a target index that is unique across input target lists.
* *RA*: Right ascension, preferably in ICRS as J2016.0
* *Dec*: Declination, preferably in ICRS as J2016.0
* *proposalid*: Proposal ID, this is provided by the observatory for each observation run. The column can more easily be defined in the config file instead of actually adding it to the input target list file.
* *catid*: Catalog ID, provided by the observatory. This is unique for an SSP sub-project, such as Galactic Archeology. The column can more easily be defined in the config file instead of actually adding it to the input target list file. 
* *obcode*: Must be a globally unique string identifying the target. This column should be generated from the config file instead of actually adding it to the input target list file. 
* *priority*: Priority of the target.
* *exp_time*: Minimum exposure time of the target.

Optional target list columns:

* *tract*: HSC catalog tract. Optional.
* *patch*: HSC catalog patch. Optional.
* *epoch*: Epoch of the coordinates. If not provided, J2016.0 is assumed. Internally, all coordinates are propagated to this epoch when the input epoch is different and proper motions are available.
* *pmra*: Proper motion in the RA direction, mas/yr.
* *err_pmra*: Error of the proper motion in the RA direction, mas/yr.
* *pmdec*: Proper motion in the Dec direction, mas/yr.
* *err_pmdec*: Error of the proper motion in the Dec direction, mas/yr.
* *parallax*: Parallax in mas.
* *err_parallax*: Error of the parallax, in mas.
* *rv*: Radial velocity in km/s.
* *err_rv*: Error of the radial velocity in km/s.

### Magnitudes and fluxes

In addition to the column listed, fluxes, magnitudes and their errors in various filters can be provided. Providing `ps1`, `gaia` or `hsc` fluxes is mandatory for flux standards. The targeting tools handle three different types of magnitudes: PSF, fiber (aperture) and total.

There are two different ways of providing the fluxes and magnitudes. In the first mode, called `filters`, each flux and magnitude column corresponds to a specific filter of an instrument. In the second mode, referred to as `bands` the magnitude and flux columns correspond to a band and an additional, string-valued column contains the name of the filter. The way of specifying the photometry for a target list is described in Chapter 4.

The targeting tools can automatically convert between fluxes and magnitudes assuming AB magnitudes. If any of the PSF, fiber or total magnitudes are missing they will be substituted from the other types without any numerical conversion.

## Specifying priorities

TBW

## Specifying the exposure time

In fact, moving fibers between visits from one target to another is possible if more than one bright target can be reached by the same cobra and the total exposure time is enough to cover the minimum exposure time of both targets.