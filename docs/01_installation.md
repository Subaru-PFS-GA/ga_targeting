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

The libraries depends on, among others, the following standard Python packages:

* numpy
* pandas
* matplotlib

## Installing Gurobi

Gurobi requires an academic licence to run. Follow the instructions here:

* https://support.gurobi.com/hc/en-us/articles/360040541251-How-do-I-obtain-a-free-academic-license

Once the license file is obtained and copied to `~/gurobi.lic` install the `gurobi` package with conda:

    $ conda install -c gurobi gurobi

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
    $ git clone 

## Environment configuration file

Every terminal session must be initialized to work with the targeting tools. This is achieved by creating an environment configuration file and sourcing the init script that will load the environment configuration file and set the necessary environmental variables. The init script also creates a `.env` file in the root of `ga_targeting` that can be used to initialize the environment for development with `vscode`.

The environment configuration files are located under `./configs/envs/`. The initialization script loads `default.sh` by default. To create a new `default.sh`, make a copy of `example.sh` and rename it to `default.sh` then modify its contents according to the actual system and file locations.

* You can ignore settings irrelevant for your use-case, such as variables starting with `PFSSPEC_` and `CMDFIT_`.
* Make sure you set the correct path to the Anaconda installation.
* Set the name of the conda environment.
* Add more paths to `PYTHONPATH` if necessary.

## Initialize the environment from a terminal

From a bash terminal, enter the project root directory and source the init script.

    $ cd ~/Subaru-PFS-GA/ga_targeting
    $ source bin/init

Note, that the environment can only be initialized when the current working directory is the project root.

The init script will automatically activate the conda environment defined in the environment configuration file. The init script also registers a set of bash functions that can be called similarly to executables.

The init script registers the following bash functions that can be executed from the command-line:

* `ga-pmap`: Generate membership probability maps from simulated CMDs of dSph
* `ga-sample`: Generate a sample from a raw catalog based on magnitude, color and probability-based selections implemented for dSphs as part of the targeting library. Its output is an input target list that contains potential targets only.
* `ga-import`: Imports all input target lists defined in a netflow configuration file. It also cross-matches the input target lists. Its outputs are the imported, cross-matched target list and each imported input target list converted to feather format.
* `ga-netflow`: Executed the netflow algorithm and saves intermediate and final fiber assignment results, as well as generates the `PfsDesign` files.
* `ga-export`: Converts fiber assignments into a file format defined by the observatory as an input for design file generation.

Some of the above tools are only used for dSph and M31 and do not apply to other fields.

## Testing the installation

Test the installation by running

    $ ga-netflow --help