#!/usr/bin/env bash

: ${MCANDOCK_MOD_PATH:=$( cd $( dirname ${BASH_SOURCE[0]} ) && pwd )}


# If NOT undefined PBS Variable
if [[ ! -z $PBS_ENVIRONMENT ]]
then
    echo "DO NOT 'qsub' this script :-)"
    exit 1
fi

if [[ -s settings.sh ]]
then
    source settings.sh
fi

export module_to_run=$1

shift 1

qsub -V $MCANDOCK_MOD_PATH/${module_to_run}.sh $([[ ! -z "$@" ]] && echo "$@")
