#!/usr/bin/env bash

: ${MCANDOCK_MOD_PATH:=$( cd $( dirname ${BASH_SOURCE[0]} ) && pwd )}

export parent_working_dir=$1

export prot_list=$parent_working_dir/all.lst

start=$2
step=$3
end=$4

shift 4

for i in `seq $start $step $end`
do

    export p_start=$i
    export p_limit=$step

    if [[ -z $MCANDOCK_test ]]
    then
        $MCANDOCK_MOD_PATH/submit_candock_module.sh dock_many_proteins -d . -N $(basename $parent_working_dir)_$i $([[ ! -z "$@" ]] && echo "$@")
    else
        $MCANDOCK_MOD_PATH/dock_many_proteins.sh
    fi
done

