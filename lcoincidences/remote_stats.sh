#!/bin/bash
cd $(dirname $0)

[ -z $2 ] && echo "2 info missing" && exit
[ -z $3 ] && echo "3 info missing" && exit

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

	ssh ${i} "
        ./storage/lcoincidences/getstats.sh $2 $3
"
done