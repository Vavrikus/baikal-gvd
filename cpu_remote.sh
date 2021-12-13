#!/bin/bash

while getopts "S:a" opt; do
    case $opt in
        S)  set -f # disable glob
            IFS=',' # split on comma characters
            array=($OPTARG) ;; # use the split+glob operator
        a)	array=("aries" "perseus" "orion" "lyra"
        	"dragon" "cassiopeia" "gidra" "hercules" "ursa" "val") ;;
    esac
done

for i in "${array[@]}"; do
	echo "Connecting to ${i}"

	ssh ${i} '
		awk '"'"'{u=$2+$4; t=$2+$4+$5; if (NR==1){u1=u; t1=t;} else print ($2+$4-u1) * 100 / (t-t1) "%"; }'"'"' \
		<(grep '"'"'cpu '"'"' /proc/stat) <(sleep 1;grep '"'"'cpu '"'"' /proc/stat)
'
done