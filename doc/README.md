---
output: html_document
---
# CANDOCK Modules

## Abstract

These scripts are designed to run various parts of candock in a modular fashion. For example: they can be used to generate the fragments for docking **without** actually doing any docking. These modular can be run independently of each other, and dependencies between modules are taken care of automatically ( IE binding site identification will take place before docking of fragments ). The list of modules is given below:

1. *find_centroids* Identifies the binding site of a protein.
2. *prep_fragments* Determines which bonds to cut for all ligands.
3. *make_fragments* Produces a PDB file for each of the fragments.
4. *dock_fragments* Docks the given fragments to a protein.
5. *link_fragments* Links the docked fragments together to form the original ligands.
6. *extract_result* Extracts all the important parts and makes a PyMOL session.
7. *design_ligands* Designs new ligands 

## General Usage

The modules can be used several different ways and all of these ways are controlled by the same set of variables.

### Submission Script Variables

#### **MCANDOCK_PATH** (*@CMAKE_INSTALL_PREFIX@*)
Path containing the candock directories.

#### **MCANDOCK_MOD_PATH** (*@CMAKE_INSTALL_PREFIX@*)
Path containing all of the the scripts! It's recommended that one exports this variable in their `.bash_profile` file. Otherwise, the script will attempt to determine this automatically, which may not work in all cases. You have been warned.

#### **MCANDOCK_VER** (*v@CANDOCK_MAJOR_VERSION@.@CANDOCK_MINOR_VERSION@.@CANDOCK_TWEAK_VERSION@*)
Version of candock to use. 

#### **MCANDOCK_NCPU** (*All availible cores*)
Number of CPUs to use. 

### Important Program Variables

#### **CANDOCK_receptor** (*receptor.pdb*)
File containing the protein that one wishes to dock to.

#### **CANDOCK_ligands** (*ligands.mol2*)
File containing all the ligands that one wishes to use.

#### **CANDOCK_centroid** (*site.cen*)
File to read the binding site from. If this file is not present or is empty, then the binding site will be determined automatically and saved to this file.

#### **CANDOCK_prep** (*prepared_ligands.pdb*)
File containing all the ligands that have been tagged with which fragments they consist of. If this file does not exist or is empty, then the ligands in `$ligands` will be tagged as such and saved to this file.

#### **CANDOCK_top_seeds_dir** (*top_seeds*)
Directory containing the docked fragments. If it does not exist, then the fragments given in `$prep` will be docked to the binding site given in `$bsite`.

#### **CANDOCK_iterative** (*0*)
Flag controlling rigid vs flexible docking. If set to `0`, then rigid docking will occur. If set to any value greater than `0`, then flexible docking will occur.

#### **CANDOCK_top_percent** (*0.05*)
Number (as a percentage) of each of the highest scored seeds to use. Default is `0.02` which means 2%.

#### **CANDOCK_max_iter** (*100*)
Maximum number of iterations to preform while linking.

#### **CANDOCK_max_possible_conf** (*20*)
Maximum of clustered confirmations to link.

### How to run CANDOCK Jobs

There are three ways to use the modules. Each way has advantages and disadvantages and the correct

#### Function Invocation

Each module exists as both a Bash script and a Bash function. This method is the quickest and dirtiest way of using the modules. It should not be used for long jobs and is best for quickly checking the modules work properly with a given set of variables. To use, start by using the following to load the functions in to the Shell environment `source $MCANDOCK_MOD_PATH/load_variables.sh`. Now you can invoke a module by simply typing the module name as you would a program name. For example, one could simply type `bsite` to do binding site identification. To change a variable, simply say `varname=varvalue`. To change the receptor, for example, use the following: `CANDOCK_receptor=1aaq.pdb`.

#### Script Invocation

This method is the recommended method for invoking a module when using a Lab or local machine. Do not use it for jobs that are to be run on the cluster. To invoke a module, simply add `$MCANDOCK_MOD_PATH` to your `$PATH` and type the module name. For example, just type `prep_fragments.sh` and your off and running! To change a variable, you **must** export it first. Do so like the following example for ligands: `export CANDOCK_ligands=new_drug.mol2`.

#### PBS Submission

If you're using the modules on an RCAC cluster, you **must** use this method if your jobs is to run for more than a few minutes. It is the most convoluted, but the most powerful. To start, jobs must be submitted using `qsub` and the full path to each module must be given. Variables are given to `qsub` as comma separated equalities using the `-v` option. For example: `qsub -v var1=value1,var2=value2`. Do not use spaces to separate values. Alternately, one can export variables and use `qsub -V`. If using a different `$MCANDOCK_MOD_PATH` than the default, ensure that this variable is exported properly as PBS blocks the scripts ability to determine it's location. For an interactive job, export `$MCANDOCK_MOD_PATH` before using the script.

So, the above was really confusing ~and badly written~ - so now there's an easier way! (tm). The `submit_candock_module.sh` command simplifies things by running qsub for you. To run `dock_fragments`, for example, use `submit_candock_module.sh dock_fragments`. You can pass `qsub` arguments as well, thus `submit_candock_module.sh dock_fragments -v var1=value1 -l myarg=whatever` is valid. Note that `-V` is passed automatically, so make sure your environment is setup properly!

## Examples

### Create a seeds data base for FLT3

The following will create a seeds database on the standby queue with 20 cores.

```bash
cd $RCAC_SCRATCH

mkdir flt3
cd flt3

dlrcsb.pl 4xuf > 4xuf.pdb
grep '^ATOM' 4xuf.pdb | grep ' A ' > receptor.pdb

cp /depot/gchopra/data/seeds_data_base/cando-ligands-3d.mol2 ./ligands.mol2

submit_candock_module.sh dock_fragments
```

### Dock molecules to a single receptor

The following will docking a small number of ligands to a single protein. It will have 200 hours to complete and will run in the gchopra queue.

```bash
cd $RCAC_SCRATCH

mkdir working_dir
cd working_dir

cp /path/to/my/protein.pdb ./receptor.pdb
cp /path/to/my/drugs.mol2 ./ligands.mol2

submit_candock_module.sh link_fragments
```

### Dock molecules to several receptors

```bash
cd $RCAC_SCRATCH

mkdir working_dir
cd working_dir

mkdir structures compounds seeds_database docking

cd structures

# Do this for all pdb codes you want to dock against
dlrcsb.pl 1PDB > 1pdb.pdb
dlrcsb.pl 2PDB > 2pdb.pdb
# .....
basename -s .pdb -a *.pdb > ../all.lst
# You can also place your own binding sites here, named
# 1pdb.cen, 2pdb.cen, etc

cd ../compounds
cp /path/to/directory/with/your/ligands.mol2 .

# If there's a few ligands ( <100ish )
prep_fragments.sh
# A lot of fagmenets
submit_candock_module.sh prep_fragments

cd ../seeds_database
dock_multiple.sh dock_fragments
# wait for all jobs to finish

cd ../docking
mkdir -p top_0.02/all_conf
cd top_0.02/all_conf
dock_multiple.sh link_fragments
```
### All Variables
