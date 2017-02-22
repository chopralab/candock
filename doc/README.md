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

#### **MCANDOCK_MOD_PATH** (*/depot/gchopra/apps/scripts/jonathan_fine/modules*)
Path containing all of the the scripts! It's recommended that one exports this variable in their `.bash_profile` file. Otherwise, the script will attempt to determine this automatically, which may not work in all cases. You have been warned.

#### **MCANDOCK_PATH** (*/depot/gchopra/apps/bin/candock*)
Path containing the candock directories.

#### **MCANDOCK_VER** (*current*)
Version of candock to use. 

#### **MCANDOCK_NCPU** (*All availible cores*)
Number of CPUs to use. 

### Program Variables

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

### All Variables

Automated Design Options:

|Option Name | Default | Description |
|------------|--------| ----------------------------- |
| seeds_till_bad | -1 | Number of times a seed must be present in the top_seeds for antitargets until it is removed from the good list |
| change_terminal_atom | ~None~ | Change non-hydrogen atoms that terminate chains to given atoms. Multiple atoms can be given. |
| fragment_bag | fragment_bag.mol2 (Implicit) | Additional fragments to be added to seeds.pdb |
| antitarget_dir | atargets (Implicit) | Directory containing PDB files. These are docked against and labeled as antitargets |
| antitarget_linking | true | Shoutd the ligands be linked for antitargets |
| add_single_atoms | ~None~ | Change hydrogens to given atoms. Multiple atoms can be given. |
| fragment_mol | fragment_mol.mol2 (Implicit) | Additional fragments to be added to seeds.pdb without rotatable bonds beingcut. |
| target_linking | true | Should the ligands be linked for target |
| seeds_to_avoid | 50 | Number of seeds from seeds.pdb to be considered for removal from determined from seeds_to_add |
| target_dir | targets (Implicit) | Directory containing PDB files. These are docked against and labeled as targets.  |
| seeds_till_good | -1 | Number of times a seed must be present in the top_seeds for targets until it is considered for addition |
| seeds_to_add | 50 | Number of seeds from seeds.pdb to be considered for addition to the ligands in prepared_ligands.pdb |
| force_seed | ~None~ | Force addition of a certain seed from seeds.pdb. Multiple seeds can be given |

Forcefield and Minimization Options:

|Option Name | Default | Description |
|------------|--------| ----------------------------- |
| gaff_dat | data/gaff.dat | Gaff DAT forcefield input file |
| max_iter_final | 100 | Maximum iterations for final minimization |
| update_freq | 10 | Update non-bond frequency |
| gaff_xml | data/gaff.xml | Gaff XML forcefield and ligand topologyoutput file |
| pos_tol | 0.00000000001 | Minimization position tolerance in Angstroms - only for KB |
| water_xml | data/tip3p.xml | Water XML parameters (and topology) input file |
| amber_xml | data/amber10.xml | Receptor XML parameters (and topology) input file |
| max_iter | 100 | Maximum iterations for minimization during linking |
| fftype | kb | Forcefield to use 'kb' (knowledge-based) or 'phy' (physics-based) |
| mini_tol | 0.0001 | Minimization tolerance |

Ligand Fragmention Options:

|Option Name | Default | Description |
|------------|--------| ----------------------------- |
| seeds_pdb | seeds.pdb | File to save full seeds into. |
| max_num_ligands | 10 | Maximum number of ligands to read in one chunk |
| prep | prepared_ligands.pdb | Prepared small molecule(s) are outputted to this filename |
| seeds | seeds.txt | Read unique seeds from this file, if itexists, and append new unique seeds if found |

Scoring Function Arguments:

|Option Name | Default | Description |
|------------|--------| ----------------------------- |
| dist | data/csd_complete_distance_distributions.txt | Select one of the interatomic distance distribution file(s) provided with thisprogram |
| func | radial | Function for calculating scores 'radial' or 'normalized_frequency' |
| comp | reduced | Atom types used in calculating reference state 'reduced' or 'complete'('reduced' includes only those atom types present in the specified receptorand small molecule, whereas 'complete' includes all atom types) |
| potential_file | potentials.txt | Output file for potentials and derivatives |
| cutoff | 6 | Cutoff length [4-15]. |
| ref | mean | Normalization method for the reference state ('mean' is averaged over all atomtype pairs, whereas 'cumulative' is a summation for atom type pairs) |
| scale | 10.0 | Scale non-bonded forces and energy for knowledge-based potential [0.0-1000.0] |
| obj_dir | obj | Output directory for objective functionand derivatives |
| step | 0.01 | Step for spline generation of non-bonded knowledge-based potential [0.0-1.0] |

Starting Input Files:

|Option Name | Default | Description |
|------------|--------| ----------------------------- |
| receptor | receptor.pdb | Receptor filename |
| ligand | ligands.mol2 | Ligand filename |

Fragment Linking Options:

|Option Name | Default | Description |
|------------|--------| ----------------------------- |
| spin | 60 | Spin degrees to rotate ligand. Allowed values are 5, 10, 15, 20, 30, 60, 90 |
| link_iter | 1000 | Maximum iterations for linking procedure |
| lower_tol_seed_dist | 2.0 | Lower tolerance on seed distance for getting initial conformations of dockedfragments |
| cuda | 1 (Implicit) | (=false)            Enable cuda iterative linker during linking |
| max_num_possibles | 200000 | Maximum number of possibles conformations considered for clustering |
| clash_coeff | 0.75 | Clash coefficient for determining whether two atoms clash by eq. dist12 s< C * (vdw1 + vdw2) |
| tol_seed_dist | 2.0 | Tolerance on seed distance in-between linking |
| max_allow_energy | 0.0 | Maximum allowed energy for seed conformations |
| max_clique_size | 3 | Maximum clique size for initial partialconformations generation |
| iterative | 1 (Implicit) | (=false)       Enable iterative minimization during linking |
| top_percent | 0.05 | Top percent of each docked seed to extend to full molecule |
| docked_dir | docked | Docked ligands output directory |
| max_possible_conf | 20 | Maximum number of possible conformations to link (-1 means unlimited) |
| docked_clus_rad | 2.0 | Cluster radius between docked ligand conformations |
| upper_tol_seed_dist | 2.0 | Upper tolerance on seed distance for getting initial conformations of dockedfragments |

Probis (binding site indentification) Options:

|Option Name | Default | Description |
|------------|--------| ----------------------------- |
| probis_min_pts | 10 | The minimum number of points (for predicted ligands) required to form a cluster |
| num_bsites | 3 | Maximum number of predicted (or given) binding sites to consider for docking |
| lig_clus_file | ligand_clusters.pdb | Ligand clusters found by ProBiS are outputted to this file |
| json | probis.json | Json-formatted ProBiS alignments outputfile |
| names | bslibdb/data/names | Directory with ligand names |
| bslib | bslibdb/bslib.txt | Read binding sites library from this file |
| probis_min_z_score | 2.5 | Minimium z-score of ligands to be considered in clustering |
| neighb | false | Allow only ligands that are in the similar regions according to REMARKs |
| probis_clus_rad | 3.0 | Cluster radius for predicted ligands byprobis |
| centroid | site.cen | Filename for reading and writing centroids |
| z_scores_file | z_scores.pdb | Binding site z-scores are outputted to this file |
| nosql | probis.nosql | NoSql-formatted ProBiS alignments output file |
| centro_clus_rad | 3.0 | Cluster radius for centroid centers |
| jsonwl | probis_with_ligands.json | Json-formatted ProBiS alignments with transposed ligands output file |
| bio | bslibdb/data/bio | Directory with ProBiS-ligands bio database |

Fragment Docking Options:

|Option Name | Default | Description |
|------------|--------| ----------------------------- |
| excluded | 0.8 | Excluded radius |
| num_univec | 256 | Number of unit vectors evenly distributed on a sphere for conformation generation |
| conf_spin | 10 | Spin degrees for conformation generation |
| gridpdb_hcp | gridpdb_hcp.pdb | Grid pdb hcp file for output |
| max_frag_radius | 16.0 | Maximum fragment radius for creating the initial rotamers |
| top_seeds_dir | top_seeds | Directory for saving top docked seeds |
| interatomic | 8.0 | Maximum interatomic distance |
| clusterfile | clustered_seeds.txt | Clustered representative docked-seed conformations output file |
| topseedsfile | top_seeds.pdb | Top seeds output file |
| clus_rad | 2.0 | Cluster radius for docked seeds |
| grid | 0.375 | Grid spacing |

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
