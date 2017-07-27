#!/usr/bin/env bash

: ${MCANDOCK_MOD_PATH:=$( cd $( dirname ${BASH_SOURCE[0]} ) && pwd )}

if [[ $PBS_ENVIRONMENT == PBS_BATCH ]]
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

: ${MCANDOCK_limit:=1000}
: ${MCANDOCK_start:=1}

__count=0
__current=1

while read i || [[ -n $i ]]; do

    if [[ "$__count" -ge "$MCANDOCK_limit" ]]
    then
        break
    fi

    if [[ "$__current" -lt "$MCANDOCK_start" ]]
    then
        __current=$((__current+1))
        continue
    fi

    export CANDOCK_receptor=$__working_dir/structures/$i.pdb
    export CANDOCK_centroid=$__working_dir/structures/$i.cen
    export CANDOCK_top_seeds_dir=$__working_dir/seeds_database/$i/top_seeds

    if [[ ! -d $i ]]
    then
        mkdir $i
    elif [[ -z $MCANDOCK_relaunch ]]
    then
        echo "Warning: previous run for $i, using previous results"
        continue
    fi

    if [[ -z $MCANDOCK_test ]]
    then
        $MCANDOCK_MOD_PATH/submit_candock_module.sh $__command -N $(basename $__working_dir)_$i \
                                                               -o $i/$(basename $__working_dir)_${i}_output.log \
                                                               -e $i/$(basename $__working_dir)_${i}_errors.log \
                                                                  $([[ ! -z "$@" ]] && echo "$@") \

    else
        $MCANDOCK_MOD_PATH/${__command}.sh
    fi
    __count=$((__count+1))

done <$__prot_list
