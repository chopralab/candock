---
output: html_document
---
# CANDOCK Modules

## Abstract

These scripts are designed to run various parts of candock in a modular fashion. For example: they can be used to generate the fragments for docking **without** actually doing any docking. These modular can be run independently of each other, and dependencies between modules are taken care of automatically ( IE binding site identification will take place before docking of fragments ). The list of modules is given below:

1. *bsite* Identifies the binding site of a protein.
2. *prep_fragments* Determines which bonds to cut for all ligands.
3. *make_fragments* Produces a PDB file for each of the fragments.
4. *dock_fragments* Docks the given fragments to a protein.
5. *link_fragments* Links the docked fragments together to form the original ligands.
6. *extract_result* Extracts all the important parts and makes a PyMOL session.

## General Usage

The modules can be used several different ways and all of these ways are controlled by the same set of variables.

### Important Variables

###### Please note that these will be renamed in a future version to avoid variable name clashes and take advantage of the new command line processing availible in 0.2.0. It will also make the transition to the next version easier ( nicknamed candock_oo ).

#### **MCANDOCK_MOD_PATH** (*/depot/gchopra/apps/scripts/jonathan_fine/modules*)
Path containing all of the the scripts, including this README! It's recommended that one exports this variable in their `.bash_profile` file. Otherwise, the script will attempt to determine this automatically, which may not work in all cases. You have been warned.

#### **MCANDOCK_PATH** (*/depot/gchopra/apps/bin/candock*)
Path containing the candock directories.

#### **MCANDOCK_VER** (*current*)
Version of candock to use. 

#### **ncpu** (*8*)
Number of CPUs to use. 

#### **receptor** (*receptor.pdb*)
File containing the protein that one wishes to dock to.

#### **ligands** (*ligands.mol2*)
File containing all the ligands that one wishes to use.

#### **bsite** (*site.cen*)
File to read the binding site from. If this file is not present or is empty, then the binding site will be determined automatically and saved to this file.

#### **prep** (*prepared_ligands.pdb*)
File containing all the ligands that have been tagged with which fragments they consist of. If this file does not exist or is empty, then the ligands in `$ligands` will be tagged as such and saved to this file.

#### **top_seeds_dir** (*top_seeds*)
Directory containing the docked fragments. If it does not exist, then the fragments given in `$prep` will be docked to the binding site given in `$bsite`.

#### **iterative** (*1*)
Flag controlling rigid vs flexible docking. If set to `0`, then rigid docking will occur. If set to any value greater than `0`, then flexible docking will occur.

#### **top_percent** (*0.02*)
Number (as a percentage) of each of the highest scored seeds to use. Default is `0.02` which means 2%.

#### **max_iter** (*10*)
Maximum number of iterations to preform while linking.

### Calling the Modules

There are three ways to use the modules. Each way has advantages and disadvantages and the correct

#### Function Invocation

Each module exists as both a Bash script and a Bash function. This method is the quickest and dirtiest way of using the modules. It should not be used for long jobs and is best for quickly checking the modules work properly with a given set of variables. To use, start by using the following to load the functions in to the Shell environment `source $CANDOCK_MOD_PATH/load_variables.sh`. Now you can invoke a module by simply typing the module name as you would a program name. For example, one could simply type `bsite` to do binding site identification. To change a variable, simply say `varname=varvalue`. To change the receptor, for example, use the following: `receptor=1aaq.pdb`.

#### Script Invocation

This method is the recommended method for invoking a module when using a Lab or local machine. Do not use it for jobs that are to be run on the cluster. To invoke a module, simply add `$CANDOCK_MOD_PATH` to your `$PATH` and type the module name. For example, just type `prep_fragments.sh` and your off and running! To change a variable, you **must** export it first. Do so like the following example for ligands: `export ligands=new_drug.mol2`.

#### PBS Submission

If you're using the modules on an RCAC cluster, you **must** use this method if your jobs is to run for more than a few minutes. It is the most convoluted, but the most powerful. To start, jobs must be submitted using `qsub` and the full path to each module must be given. Variables are given to `qsub` as comma separated equalities using the `-v` option. For example: `qsub -v var1=value1,var2=value2`. Do not use spaces to separate values. Alternately, one can export variables and use `qsub -V`. If using a different `$CANDOCK_MOD_PATH` than the default, ensure that this variable is exported properly as PBS blocks the scripts ability to determine it's location. For an interactive job, export `$CANDOCK_MOD_PATH` before using the script.

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
