#!/bin/bash
cd $(dirname $0)

[ -z $1 ] && echo "1 info missing" && exit
[ -z $2 ] && echo "2 info missing" && exit

seq -w 0 $(($1-1)) | parallel -j+0 --progress ./PE.exe $2 {} $3
# seq -w 0 $1 | parallel -j+0 -Sval, --progress ./PE.exe $2 {}  
# seq -w 0 2 | parallel -j+0 -Sval,gidra,: --progress 'ls && echo {}'
# usage: ./getstats.sh (number of runs) (declination) (iterate right ascension)
# number of runs is number of simulations divided by 10000