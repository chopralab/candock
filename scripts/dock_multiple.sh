#!/usr/bin/env bash

: ${MCANDOCK_MOD_PATH:=$( cd $( dirname ${BASH_SOURCE[0]} ) && pwd )}

if [[ ! -z $PBS_ENVIRONMENT ]]
then
    echo "DO NOT 'qsub' this script :-)"
    exit 1
fi

__command=$1

__working_dir=$2
__prot_list=$__working_dir/all.lst

shift 2

export CANDOCK_prep=$__working_dir/compounds/prepared_ligands.pdb
export CANDOCK_seeds=$__working_dir/compounds/seeds.txt
export CANDOCK_seeds_pdb=$__working_dir/compounds/seeds.pdb

: ${MCANDOCK_relaunch:=0}

: ${MCANDOCK_limit:=1000}
: ${MCANDOCK_start:=1}

__count=0
__current=1

for i in `cat $__prot_list`
do

    if [[ "$__count" -ge "$MCANDOCK_limit" ]]
    then
        break
    fi

    if [[ "$__current" -lt "$MCANDOCK_start" ]]
    then
        __current=$((__current+1))
        continue
    fi

    export CANDOCK_receptor=$working_dir/structures/$i.pdb
    export CANDOCK_centroid=$working_dir/structures/$i.cen
    export CANDOCK_top_seeds_dir=$working_dir/seeds_database/$i/top_seeds

    if [[ ! -d $i ]]
    then
        mkdir $i
    elif [[ "$MCANDOCK_relaunch" -eq "0" ]]
        echo "Warning: previous run for $i, using previous results"
        continue
    fi

    $MCANDOCK_MOD_PATH/submit_candock_module.sh $__command -N $(basename $__working_dir)_$i $([[ ! -z "$@" ]] && echo "$@") -o $i/$(basename $working_dir)_${i}_output.log -e $i/$(basename $working_dir)_${i}_errors.log

    __count=$((__count+1))

done
