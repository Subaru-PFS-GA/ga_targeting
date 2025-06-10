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
* spt_operational_database
* (ets_fiberalloc)

Most of these are necessary to calculate fiber collisions.

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
export PFS_UTILS="$HOME/Subaru-PFS/pfs_utils/python"
export PFS_INSTDATA="$HOME/Subaru-PFS/pfs_instdata/python"
export PFS_ICS_COBRAOPS="$HOME/Subaru-PFS/ics_cobraOps/python"
export PFS_ICS_COBRA_CHARMER="$HOME/Subaru-PFS/ics_cobraCharmer/python"
export PFS_ICS_FPSACTOR="$HOME/Subaru-PFS/ics_fpsActor/python"
export PFS_ETS_FIBER_ASSIGNER="$HOME/Subaru-PFS/ets_fiberalloc"
export PFS_SPT_OPDB="$HOME/Subaru-PFS/spt_operational_database/python"

export PFS_ISOCHRONES="$HOME/Subaru-PFS-GA/pfs_isochrones/python"

export PFS_TARGETING_ROOT="$HOME/Subaru-PFS-GA/ga_targeting"
export PFS_TARGETING_DATA="$HOME/data/targeting"
export PFS_TARGETING_TEMP="$HOME/tmp/targeting"
export PFS_TARGETING_CONDAPATH="$HOME/miniconda3"
export PFS_TARGETING_CONDAENV="targeting"
export PFS_TARGETING_PYTHONPATH="$PFS_UTILS:$PFS_INSTDATA:$PFS_ICS_COBRAOPS:$PFS_ICS_COBRA_CHARMER:$PFS_ETS_FIBER_ASSIGNER:$PFS_SPT_OPDB:$PFS_ISOCHRONES"
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