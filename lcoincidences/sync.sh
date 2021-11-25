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
        cd storage/lcoincidences/data || (mkdir -p storage/lcoincidences/data && echo "Created directories.")
'
    rsync compile.sh data_pdf.C getstats.sh pseudo_exp.C cos_theta.root skyfit.C ${i}":~/storage/lcoincidences"
    rsync ../horizon.C ../EventLoop.cpp ../EventLoop.h ../Instrumentor.h ../threading.h ../profilling.h ../transformations.h ${i}":~/storage"

    if [ "${i}" != "val" ]
    then
        rsync prototype_bashrc.txt ${i}":~/.bashrc"
    fi
done