# PFS target selection for GA

TBW

## 0.1 Installation

The targeting tool is currently distributed as Python source. It is intended to be used as a command-line tool and not as a library, though parts of it, such as plotting, can be useful for other purposes.

The targeting tool depends on the `ga_isochones` library, the `gurobipy` linear program optimization library and several standard Python libraries. See Chapter 1 for more details on installation.

## 0.2 Terminology

A *field* is an area of a sky, for example, the field around a dwarf spheroidal galaxy. A field can be covered by several *pointings* which can overlap. The footprint of a pointing is determined by the coverage of the Prime Focus Instrument (PFI). A pointing is *visited* multiple times to accumulate enough *exposure time*. The telescope boresight position and the instrument position angle should be the same for each visit of a pointing but the fiber configuration can be different. A *visit* refers to one or more *exposures* which are reduced by the 2DRP into 1-dimensional spectra.

Pointings can be organized into *stages* and stages can have *priorities*. This information is used by the observation scheduler (Yasuda, 2025) to ensure *early completion* of GA pointings.

A *source target lists* is a data file containing all potential targets. The targeting tool will automatically *assign fibers* to a subset of the potential targets based on *target priority*. Target priority is an integer number defined by the user for each target. By default, smaller integers indicate higher priorities. A *minimum exposure time* must be defined for every target which constraints how many visits will be allocated to a target.

*Flux standards* are F stars used for flux calibration and should be selected to be slightly brighter than typical targets and to cover the field of view as uniformly as possible. *Sky fibers* are fibers positioned on coordinates with no known objects. The sky fibers must cover the field of view, as well as the spectrograph slits as uniformly as possible. Flux standards and sky positions are supplied to the targeting tool as files containing data tables.

`Netflow` (Fabricius et al.) is an implementation of the Network Flow algorithm for fiber assignment. The GA targeting tool uses a slightly updated version of the algorithm with a significantly optimized codebase.

`Gurobi` is a linear program solver used by `Netflow` to optimize the network flow problem. A good solution of the network flow problem is a good fiber assignment. Some solution might be better than others based on the various *costs* associated with observing or not observing a particular targets. Finding the best solution (or solutions) is often not possible due to computational time constraints. A Netflow problem is *infeasible* if no fiber assignment exists which can satisfy all constraints - for example, it is possible that low number of flux standard renders a Netflow problem infeasible. In this case some constraints, such as the uniformity of sky coverage, must be relaxed or the number of flux standards increased, if possible.

A *configuration file* is a JSON file summarizing all settings of a targeting run which includes the list of pointing, targets lists and all configuration settings of the various components.

A *feather file* a binary file format to store data tables, such as Pandas data frames. The targeting tool uses this format internally to write out intermediate data.

A *PfsDesign* file is a FITS file that contains the final fiber assignment and is an input to the telescope operation during the observing runs.

A *proposal ID* or `proposalId` is an identifier assigned to the entire SSP observing run by the observatory.

A *catalog ID* or `catId` is the same for all GA targets and is assigned by the observatory. Flux standards and sky positions have their of `catId`s.

An *observation code* or `obCode` is a unique identifier assigned to the targets by the observer or SSP working group.

## 0.1 Overview

Target selection consist of the following steps.

1. Prepare input target lists
2. Create configuration file
3. Import targets lists
4. Run netflow
5. Extract and validate results.

### 0.1.1 Prepare input target lists

Target lists must be provided in a format that either pandas can read or in a standard FITS file with a BINTABLE extension.

Although a target list can be further filtered by simple criteria, such as magnitude cuts, it should only contain targets that are potentially good GA targets.

If possible, each target list should contain sources selected from a single catalog. There is no need to cross-match the catalogs unless required by the target selection or priority assignment. For example, to combine NIR and optical bands from different instruments. The targeting tool will cross-match the input target lists, identify and correctly handle objects that appear in multiple input target lists.

When the same object is found in several input target lists, the highest priority and longest exposure time will be used. If coordinates slightly differ, the coordinates from the first target list will be used in which the object appears.

Input target lists should contain the coordinates (preferably in the ICRS frame) at a known epoch (preferably at J2016.0), proper motion components and parallaxes, as well as fluxes and magnitudes at all available bands. If not given at the reference epoch of J2016.0, the targeting tool will propagate the position of the targets to J2016.0 automatically. If only the flux or magnitude is given for a certain band, the targeting tool will automatically calculate the other assuming AB magnitudes. The PfsDesign files require fluxes.

There are two different ways of specifying fluxes and magnitudes: one data file column per filter or one data file column per band. The first approach assumes that a data file column contains magnitudes measured in the same filter while the second approach uses three data file columns to store a magnitude: two columns store the numerical value and its error while a third column stores the name of the filter.

In addition to the aforementioned columns, target lists must contain a column with the minimum required exposure time and can contain a column with the priority of the target, but priorities can be defined later, in the configuration file.

Because they're large, we store the input target lists on the file system somewhere, never in the code repository (git).

### 0.1.2 Create a configuration file

The targeting tool uses a configuration file to define the pointing centers, position angles, exposure times, input target lists, object class priorities, Netflow settings, Gurobi settings etc.

### 0.1.3 Import targets

### 0.1.4 Run netflow

### 0.1.5 Evaluate results

TBW: write about some generic features of the command-line scripts:

* output must be an empty directory
* command-line and arguments are saved to the output directory
* log files
* evaluation notebooks