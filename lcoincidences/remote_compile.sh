#!/bin/bash
cd $(dirname $0)

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
        home_dir=$(pwd | cat)
        if [ "$home_dir" = "/home/vavrik" ]
        then
            source /opt/root/bin/thisroot.sh
        fi
        ./storage/lcoincidences/compile.sh
'
done