#!/usr/bin/env bash

: ${MCANDOCK_MOD_PATH:=$( cd $( dirname ${BASH_SOURCE[0]} ) && pwd )}

if [[ $PBS_ENVIRONMENT == PBS_BATCH ]]
then
    echo "DO NOT 'qsub' this script :-)"
    exit 1
fi

export parent_working_dir=$1

export prot_list=$parent_working_dir/all.lst

__start=$2
__step=$3
__end=$4

shift 4

for i in `seq $__start $__step $__end`
do

    export MCANDOCK_start=$i
    export MCANDOCK_limit=$step

    if [[ -z $MCANDOCK_test ]]
    then
        $MCANDOCK_MOD_PATH/submit_candock_module.sh dock_many_proteins -d . -N $(basename $parent_working_dir)_$i $([[ ! -z "$@" ]] && echo "$@")
    else
        $MCANDOCK_MOD_PATH/dock_many_proteins.sh
    fi
done
