#!/bin/bash
cd $(dirname $0)
cd data

nSign=' find ./*/ | grep "nSign" | head -1 | grep -Po "(?<=/).*(?=txt)"txt | grep -Po "(?<=/).*(?=txt)"txt '
cat `find -wholename "./*/$(eval $nSign)" -not -path "./merged/*"` >> merged/$(eval $nSign)
rm `find -wholename "./*/$(eval $nSign)" -not -path "./merged/*"`

tStat=' find ./*/ | grep "tStat" | head -1 | grep -Po "(?<=/).*(?=txt)"txt | grep -Po "(?<=/).*(?=txt)"txt '
cat `find -wholename "./*/$(eval $tStat)" -not -path "./merged/*"` >> merged/$(eval $tStat)
rm `find -wholename "./*/$(eval $tStat)" -not -path "./merged/*"`