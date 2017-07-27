#!/usr/bin/env bash

: ${MCANDOCK_MOD_PATH:=$( echo $( cd $( dirname ${BASH_SOURCE[0]} ) && pwd ) )}

: ${MCANDOCK_poses:=1}

echo "Running extract in $PWD"

mkdir extract
cd extract

$MCANDOCK_MOD_PATH/extract_candock_to_scores_and_models.pl $([[ "$MCANDOCK_pymol" -eq 0 ]] && echo -f) -n $MCANDOCK_poses ../docked/*.pdb > ../scores.csv

if [[ ! -z $MCANDOCK_pymol ]]
then
    $MCANDOCK_pymol -c create_pse.py
    mv all.pse ..
fi

cd ..
rm -fr extract
