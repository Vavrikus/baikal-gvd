#!/bin/bash
cd $(dirname $0)

while getopts "S:al" opt; do
    case $opt in
        S)  set -f # disable glob
            IFS=',' # split on comma characters
            array=($OPTARG) ;; # use the split+glob operator
        a)	array=("aries" "perseus" "orion" "lyra"
        	"dragon" "cassiopeia" "gidra" "hercules" "ursa" "val") ;;
        l)  cd ~/work/Baikal-GVD/lcoincidences/data
            nSign=' ls | grep "nSign" | head -1 | grep -o -P "(?<=nSign).*(?=_)" '
            cat data_nSign$(eval $nSign)* > fdata_nSign$(eval $nSign)_all.txt && rm data_nSign$(eval $nSign)*

            tStat=' ls | grep "tStat" | head -1 | grep -o -P "(?<=tStat).*(?=_)" '
            cat data_tStat$(eval $tStat)* > fdata_tStat$(eval $tStat)_all.txt && rm data_tStat$(eval $tStat)*

            cd ~/work/Baikal-GVD/lcoincidences/data/local || mkdir ~/work/Baikal-GVD/lcoincidences/data/local
            cd ~/work/Baikal-GVD/lcoincidences/data/local 2> /dev/null

            mv ../fdata* .
            ;;
    esac
done

for i in "${array[@]}"; do
	echo "Connecting to ${i}"

	ssh ${i} '
        cd ~/storage/lcoincidences/data
        nSign='"'"' ls | grep "nSign" | head -1 | grep -o -P "(?<=nSign).*(?=_)" '"'"'
        cat data_nSign$(eval $nSign)* > fdata_nSign$(eval $nSign)_all.txt && rm data_nSign$(eval $nSign)*

        tStat='"'"' ls | grep "tStat" | head -1 | grep -o -P "(?<=tStat).*(?=_)" '"'"'
        cat data_tStat$(eval $tStat)* > fdata_tStat$(eval $tStat)_all.txt && rm data_tStat$(eval $tStat)*
'
    cd ~/work/Baikal-GVD/lcoincidences/data/${i} || mkdir ~/work/Baikal-GVD/lcoincidences/data/${i}
    cd ~/work/Baikal-GVD/lcoincidences/data/${i} 2> /dev/null

    rsync --progress ${i}:~/storage/lcoincidences/data/fdata* .

    ssh ${i} '
        cd ~/storage/lcoincidences/data
        rm fdata*
'
done

#test -e filename && echo exists || echo does not