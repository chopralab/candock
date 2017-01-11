#!/usr/bin/env bash

command=$1

working_dir=$2
prot_list=$working_dir/all.lst

shift 2

export prep=$working_dir/compounds/prepared_ligands.pdb
export seeds=$working_dir/compounds/seeds.txt

if [[ -s settings.sh ]]
then
	source settings.sh
fi

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

	export receptor=$working_dir/structures/$i.pdb
	export bsite=$working_dir/structures/$i.cen
	export top_seeds_dir=$working_dir/seeds_database/$i/top_seeds

	if [[ ! -d $i ]]
	then
		mkdir $i
	else
		echo "Warning: previous run for $i, using previous results"
		continue
	fi

	cd $i

	submit_candock_module.sh $command -N $(basename $working_dir)_$i $([[ ! -z "$@" ]] && echo "$@")

	count=$((count+1))

	cd ..
done
