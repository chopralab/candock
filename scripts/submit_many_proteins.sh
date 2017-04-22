#!/usr/bin/env bash

export parent_working_dir=$1

export prot_list=$parent_working_dir/all.lst

start=$2
step=$3
end=$4

shift 4

if [[ -s settings.sh ]]
then
	echo "Reading settings.sh"
	source settings.sh
fi

for i in `seq $start $step $end`
do

    export p_start=$i
    export p_limit=$step
    
    submit_candock_module.sh dock_many_proteins -N $(basename $parent_working_dir)_$i $([[ ! -z "$@" ]] && echo "$@")
done

