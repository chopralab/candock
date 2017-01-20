#!/usr/bin/env bash

: ${MCANDOCK_MOD_PATH:=$( cd $( dirname ${BASH_SOURCE[0]} ) && pwd )}

if [[ ! -z $PBS_ENVIRONMENT ]]
then
    echo "DO NOT 'qsub' this script :-)"
    exit 1
fi

command=$1

working_dir=$2
prot_list=$working_dir/all.lst

shift 2

export CANDOCK_prep=$working_dir/compounds/prepared_ligands.pdb
export CANDOCK_seeds=$working_dir/compounds/seeds.txt

: ${limit:=1000}
: ${start:=1}

count=0
current=1

for i in `cat $prot_list`
do

    if [[ "$count" -ge "$limit" ]]
    then
        break
    fi

    if [[ "$current" -lt "$start" ]]
    then
        current=$((current+1))
        continue
    fi

    export CANDOCK_receptor=$working_dir/structures/$i.pdb
    export CANDOCK_centroid=$working_dir/structures/$i.cen
    export CANDOCK_top_seeds_dir=$working_dir/seeds_database/$i/top_seeds

    if [[ ! -d $i ]]
    then
        mkdir $i
    else
        echo "Warning: previous run for $i, using previous results"
        continue
    fi

    cd $i

    $MCANDOCK_MOD_PATH/submit_candock_module.sh $command -N $(basename $working_dir)_$i $([[ ! -z "$@" ]] && echo "$@")

    count=$((count+1))

    cd ..
done
