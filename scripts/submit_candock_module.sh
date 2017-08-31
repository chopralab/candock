#!/usr/bin/env bash

: ${MCANDOCK_MOD_PATH:=$( cd $( dirname ${BASH_SOURCE[0]} ) && pwd )}


# If NOT undefined PBS Variable
if [[ $PBS_ENVIRONMENT == PBS_BATCH ]]
then
    echo "DO NOT 'qsub' this script :-)"
    exit 1
fi

export module_to_run=$1

shift 1

__variable_list=`env | grep CANDOCK | sed -e ':a;N;$!ba;s/\n/,/g'`

qsub -v module_to_run,$__variable_list $MCANDOCK_MOD_PATH/${module_to_run}.sh -d . $([[ ! -z "$@" ]] && echo "$@")

