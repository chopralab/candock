#!/usr/bin/env bash
#PBS -l walltime=3:50:00
#PBS -q standby
#PBS -l nodes=1:ppn=8
#PBS -l naccesspolicy=shared

: ${MCANDOCK_MOD_PATH:=$([[ "$PBS_ENVIRONMENT" -eq "PBS_BATCH" ]] \
	&& echo /depot/gchopra/apps/bin/candock/modules/v0.2.0 \
	|| echo $( cd $( dirname ${BASH_SOURCE[0]} ) && pwd ) )}

: ${module_to_run:=$PBS_JOBNAME}

echo "Using modules in $MCANDOCK_MOD_PATH"

if [[ ! -s $MCANDOCK_MOD_PATH/load_variables.sh ]]
then
	echo "Modules not found in \$MCANDOCK_MOD_PATH!"
	echo "Please set this variable!"
	exit 1
fi

source $MCANDOCK_MOD_PATH/load_variables.sh

if [[ ! -z $PBS_ENVIRONMENT && "$PBS_ENVIRONMENT" -eq "PBS_BATCH" ]]
then
	$(basename $module_to_run .sh)
else
	$(basename $(basename $0 .sh) )
fi

echo "All caclulations for $(basename $0) are complete"

