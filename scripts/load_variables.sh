#!/usr/bin/env bash

# Load default variables

: ${MCANDOCK_PATH:=/depot/gchopra/apps/bin/candock}
: ${MCANDOCK_VER:=v0.2.1}

# Note: Untested with use in load_variables.sh as this script should be sourced
# Using other scripts!
: ${MCANDOCK_MOD_PATH:="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"}
: ${MCANDOCK_SUP_PATH:=/depot/gchopra/apps/bin/candock/}

: ${receptor:=receptor.pdb}
: ${ligands:=ligands.mol2}

: ${iterative:=1}

: ${ncpu:=8}

: ${bsite:=site.cen}

: ${prep:=prepared_ligands.pdb}
: ${seeds:=seeds.txt}

: ${top_seeds_dir:=top_seeds}

: ${top_percent:=0.02}
: ${max_iter:=10}
: ${CANDOCK_max_possible_conf:=-1}

: ${frag_separator:=$MCANDOCK_MOD_PATH/splitfragments.pl}

: ${extractor:=$MCANDOCK_MOD_PATH/extract_candock_to_scores_and_models.pl}
: ${poses:=1}
: ${pymol:=/depot/gchopra/apps/pymol/pymol}

# Check if were running in a PBS environment

if [[ ! -z $PBS_O_WORKDIR ]]
then
	: ${working_dir:=$PBS_O_WORKDIR}
else
	: ${working_dir:=$PWD}
fi

function bsite() {
	cd $working_dir
	echo "Running bsite in $PWD"

	if [[ -s $bsite ]]
	then
		echo "Existing binding site found in $bsite"
		return
	fi

	if [[ ! -s $receptor ]]
	then
		echo "No receptor found in $receptor"
		exit 1
	fi

	ln -s $MCANDOCK_SUP_PATH/bslibdb  bslibdb

	module load gcc/5.2.0

	SECONDS=0

	$MCANDOCK_PATH/$MCANDOCK_VER/test_bsite --ncpu $ncpu \
		--receptor $receptor \
		--centroid $bsite

	echo "bsite took $SECONDS seconds to complete"

	module unload gcc

	unlink bslibdb

	chmod g+r $bsite
	chmod g+r gridpdb_hcp.pdb
	chmod g+r probis.json
	chmod g+r probis.nosql
	chmod g+r probis_with_ligands.json
	chmod g+r "`basename $receptor .pdb`.srf"
	chmod g+r z_scores.pdb
}

function prep_fragments() {
	cd $working_dir
	echo "Running fragment in $PWD"

	if [[ -s $prep ]]
	then
		echo "Existing fragmentation found in $prep"
		return
	fi

	if [[ ! -s $ligands ]]
	then
		echo "No ligands found in $ligands"
		exit 1
	fi

	module load gcc/5.2.0

	SECONDS=0

	$MCANDOCK_PATH/$MCANDOCK_VER/test_fragmenting --ncpu $ncpu \
		--ligand $ligands \
		--prep $prep

	echo "Fragmenting took $SECONDS seconds to complete"

	module unload gcc

	chmod g+r $prep
	chmod g+r seeds.txt
}

function make_fragments() {
	cd $working_dir
	echo "Running make_fragments in $PWD"

	prep_fragments

	perl $frag_separator $prep > separated_fragments.pdb

	chmod g+r separated_fragments.pdb	
}

function dock_fragments() {
	cd $working_dir
	echo "Running dock in $PWD"

	if [[ ! -s $receptor ]]
	then
		echo "No receptor found in $receptor"
		exit 1
	fi

	if [[ -d $top_seeds_dir ]]
	then
		echo "Seeds database found at $top_seeds_dir"
		return
	fi

	bsite
	prep_fragments

	ln -s $MCANDOCK_SUP_PATH/obj   obj
	ln -s $MCANDOCK_SUP_PATH/data  data

	module load gcc/5.2.0

	SECONDS=0

	$MCANDOCK_PATH/$MCANDOCK_VER/test_clq --ncpu $ncpu \
		--receptor $receptor \
		--ligand $prep \
		--prep   $prep \
		--centroid $bsite \
		--top_seeds_dir $top_seeds_dir

	echo "Docking of fragments took $SECONDS seconds to complete"

	module unload gcc

	unlink obj
	unlink data

	chmod g+r -R $top_seeds_dir
}

function link_fragments() {
	cd $working_dir
	echo "Running link in $PWD"

	if [[ ! -s $receptor ]]
	then
		echo "No receptor found in $receptor"
		exit 1
	fi

	dock_fragments

	if [[ -s $seeds ]]
	then
		numseeds=$(wc -l $seeds | awk '{print $1}')
		numdockd=$(ls $top_seeds_dir/*/top_seeds.pdb | wc -w)
		if [[ "$numseeds" -ne "$numdockd" ]]
		then
			echo "Warning: inconsistant seeds.txt($numseeds) and seeds data_base($numdockd)"
		fi
	else
		echo "Warning: no seeds.txt file, cannot check consistancy"
	fi

	ln -s $MCANDOCK_SUP_PATH/obj   obj
	ln -s $MCANDOCK_SUP_PATH/data  data

	module load gcc/5.2.0

	SECONDS=0

	$MCANDOCK_PATH/$MCANDOCK_VER/test_link --ncpu $ncpu \
		--receptor $receptor \
		--prep $prep \
		--ligand $prep \
		--top_seeds_dir $top_seeds_dir \
		--top_percent $top_percent \
		--max_iter $max_iter \
		$([[ "$iterative" -gt 0 ]] && echo --iterative)

	echo "Linking of fragments took $SECONDS seconds to complete"

	module unload gcc

	unlink obj
	unlink data

	chmod -R g+r docked/
}

function extract_result() {
	cd $working_dir
	echo "Running extract in $PWD"
	
	mkdir extract
	cd extract

	perl $extractor -n $poses ../docked/*.pdb > ../scores.csv

	$pymol -c create_pse.py
	mv all.pse ..

	cd ..
	rm -fr extract

	chmod g+r scores.csv
	chmod g+r all.pse
}

