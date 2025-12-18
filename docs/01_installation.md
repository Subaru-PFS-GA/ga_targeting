# 1 Installation

## Prerequisites

The targeting toolkit depends on the following libraries:

* *ga_isochrones*: Isochrone interpolation, for isochrones based selections.
* *ga_cmdfit*: Weak dependency, for membership-based selections.

In addition, the following PFS libraries are necessary:

* *pfs_utils*: Coordinate transformation, focal plane positions.
* *pfs_instdata*: Instrument calibration and configuration.
* *ics_cobraOps*: Cobra operations.
* *ics_cobraCharmer*: Cobra operations.
* *ics_fpsActor*: Cobra trajectory simulations.
* *ets_fiberalloc*: "Original" netflow implementation, might be a dependency of other PFS libraries.
* *spt_operational_database*: Operational database, dependency of other PFS libraries.

Either install these libraries in the conda environment, or clone them from github and set their paths in the environment configuration file (see below).

## Cloning the library

Please make sure that you set up github access with an SSH key instead of relying on password authentication. Certain git operations might require entering the password several times. Refer to the following github page for details:

* https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent

Create a directory named `Subaru-PFS-GA` and clone the necessary libraries:

    $ mkdir ~/Subaru-PFS-GA
    $ cd ~/Subaru-PFS-GA
    $ git clone git@github.com:Subaru-PFS-GA/ga_common.git
    $ git clone git@github.com:Subaru-PFS-GA/ga_isochrones.git
    $ git clone git@github.com:Subaru-PFS-GA/ga_targeting.git

Similarly, clone the PFS libraries, if not already installed in the environment:

    $ mkdir ~/Subaru-PFS
    $ cd ~/Subaru-PFS
    $ git clone ...

Where the following libraries must be cloned from https://github.com/Subaru-PFS:

* pfs_utils
* pfs_instdata
* ics_cobraOps
* ics_cobraCharmer
* ics_fpsActor
* ics_utils
* spt_operational_database
* (ets_fiberalloc)

Most of these are necessary to calculate fiber collisions. `ets_fiberalloc` is not strictly a dependency as its functionality is re-implemented in the GA targeting library.

## Install a new conda environment

Assumed that you already have a working Anaconda installation, run the following `conda` command to initialize a new environment.

    $ cd ~/Subaru-PFS-GA/ga_targeting
    $ conda env create -n targeting -f ./setup/targeting.yaml

All necessary python dependecies are listed in `./setup/targeting.yaml`, including Gurobi, but Gurobi requires additional manual setup.

## Installing Gurobi

Gurobi requires an academic licence to run. Follow the instructions here:

* https://support.gurobi.com/hc/en-us/articles/360040541251-How-do-I-obtain-a-free-academic-license

Basically, there are two types of free academic licenses: a long-term license that's associated with a user and a host computer and a short-term license which is associated with a user but can be moved from computer to computer. Either type of license should work but for server installation the former might be better as it requires less attention. On the other hand, Gurobi issues a single free license per registered academic user, so when portability is a goal, the second type of license should be chosen.

Once the license file is obtained, copy it to `~/gurobi.lic`.

    $ conda install -c gurobi gurobi

The license file looks something like this:

```gurobi.lic
LICENSEID=0000000
TYPE=ACADEMIC
VERSION=12
HOSTNAME=xxx.xxx.edu
HOSTID=xxxxxxxx
SOCKETS=2
USERNAME=dobos
EXPIRATION=2026-03-11
KEY=XXXXXXXX
CKEY=XXXXXXXX
```

## Environment configuration file

This environment configuration must be done independently of setting up the Anaconda environment.

Every terminal session must be initialized to work with the targeting tools. This is achieved by creating an environment configuration file and sourcing the init script that will load the environment configuration file and set the necessary environmental variables. The init script also creates a `.env` file in the root of `ga_targeting` that can be used to initialize the environment for development with `vscode`.

The environment configuration files are located under `./configs/envs/`. The initialization script loads `default.sh` by default. To create a new `default.sh`, make a copy of `example.sh` and rename it to `default.sh` then modify its contents according to the actual system and file locations.

* You can ignore settings irrelevant for your use-case, such as variables starting with `PFSSPEC_` and `CMDFIT_`.
* Make sure you set the correct path to the Anaconda installation.
* Set the name of the conda environment.
* Add more paths to `PYTHONPATH` if necessary.

Here is a minimal example environment file that works with the directory names used above:

```./configs/envs/default.sh
SUBARU_PFS_ROOT="~/Subaru-PFS"
SUBARU_PFS_GA_ROOT="$HOME/Subaru-PFS-GA"

export PFSSPEC_ROOT="$SUBARU_PFS_GA_ROOT/ga_pfsspec-all"
export PFSSPEC_DATA="$HOME/data/pfsspec"

export CMDFIT_ROOT="$SUBARU_PFS_GA_ROOT/ga_cmdfit"
export CMDFIT_DATA="$HOME/data/cmdfit"

export PFS_UTILS="$SUBARU_PFS_ROOT/pfs_utils"
export PFS_INSTDATA="$SUBARU_PFS_ROOT/pfs_instdata"
export PFS_ICS_COBRAOPS="$SUBARU_PFS_ROOT/ics_cobraOps"
export PFS_ICS_COBRA_CHARMER="$SUBARU_PFS_ROOT/ics_cobraCharmer"
export PFS_ICS_FPSACTOR="$SUBARU_PFS_ROOT/ics_fpsActor"
export PFS_ICS_UTILS="$SUBARU_PFS_ROOT/ics_utils"
export PFS_ETS_FIBER_ASSIGNER="$SUBARU_PFS_ROOT/ets_fiberalloc"
export PFS_SPT_OPDB="$SUBARU_PFS_ROOT/spt_operational_database"
export PFS_GA_ISOCHRONES="$SUBARU_PFS_GA_ROOT/ga_isochrones"
export PFS_GA_COMMON="$SUBARU_PFS_GA_ROOT/ga_common"

export PFS_TARGETING_DEBUGPORT=5678
export PFS_TARGETING_ROOT="$SUBARU_PFS_GA_ROOT/ga_targeting"
export PFS_TARGETING_DATA="$HOME/data/targeting"
export PFS_TARGETING_TEMP="$HOME/tmp/targeting"
export PFS_TARGETING_CONDAPATH="$HOME/miniconda3"
export PFS_TARGETING_CONDAENV="targeting"

export PFS_TARGETING_MODULES="pfs_utils:$PFS_UTILS:python
pfs_instdata:$PFS_INSTDATA:python
ics_cobraOps:$PFS_ICS_COBRAOPS:python
ics_cobraCharmer:$PFS_ICS_COBRA_CHARMER:python
ics_fpsActor:$PFS_ICS_FPSACTOR:python
ics_utils:$PFS_ICS_UTILS:python
ets_fiberalloc:$PFS_ETS_FIBER_ASSIGNER:.
spt_operational_database:$PFS_SPT_OPDB:python
datamodel:$PFS_DATAMODEL:python
ga_isochrones:$PFS_GA_ISOCHRONES:python
ga_common:$PFS_GA_COMMON:python
ga_targeting:$PFS_TARGETING_ROOT:python
ga_pfsspec:$PFSSPEC_ROOT:python"
```

Note, that `PFS_TARGETING_CONDAPATH` must point to your Anaconda installation. Find it by running

    $ echo $CONDA_EXE

The path should not contain the `bin/conda` part of the path.

## Initialize the environment from a terminal

From a bash terminal, enter the project root directory and source the init script.

    $ cd ~/Subaru-PFS-GA/ga_targeting
    $ source bin/init

Note, that the environment can only be initialized when the current working directory is the project root.

The init script will automatically activate the conda environment defined in the environment configuration file. The init script also registers a set of bash functions that can be called similarly to executables.

The init script registers the following bash functions that can be executed from the command-line:

* `ga-import`: Imports all input target lists defined in a netflow configuration file. It also cross-matches the input target lists. Its outputs are the imported, cross-matched target list and each imported input target list converted to feather format.
* `ga-netflow`: Executed the netflow algorithm and saves intermediate and final fiber assignment results, as well as generates the `PfsDesign` files.
* `ga-export`: Converts fiber assignments into a file format defined by the observatory as an input for design file generation.

In addition to these, the targeting library implements target list preparation based on membership probability. These tools are only used for dSph and M31 and do not apply to other fields. For these purposes, two additional commands are available: 

* `ga-pmap`: Generate membership probability maps from simulated CMDs of dSph
* `ga-sample`: Generate a sample from a raw catalog based on magnitude, color and probability-based selections implemented for dSphs as part of the targeting library. Its output is an input target list that contains potential targets only.

## Testing the installation

Test the installation by running

    $ cd ~/Subaru-PFS-GA/ga_targeting
    $ ga-netflow --help

## Updating the installation

In order to bump the version of the libraries, pull the latest changes from the remote git respositories.

To update the targeting library, run

    $ cd $PFS_TARGETING_ROOT
    $ git fetch
    $ git pull origin master

In certain cases you might need to pull a specific branch, or a tag, so replace `master` with the desired branch or tag name.

To update the PFS libraries, you typically need to pull a specific release verions. For example, you can update `pfs_utils` by running

    $ cd $PFS_UTILS/..
    $ git fetch
    $ git checkout tags/w.2025.41

