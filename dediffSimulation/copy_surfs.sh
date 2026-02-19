#!/bin/bash

from=${1}
ois=${2}
opref=${3}
oTp=${4}

nis=${5}
npref=${6}
nTp=${7}
surfDir=prepSurf_perturb

for set in 1 2 3 4 5 6 7 8 #9 10 11 12 13 14 15 16
do

olddir=${from}/Set_${set}_in${ois}_${opref}_Tp${oTp}/${surfDir}
newdir=Set_${set}_in${nis}_${npref}_Tp${nTp}


mkdir -p ${newdir}
mkdir -p ${newdir}/${surfDir}

echo ${olddir}

echo "to"

echo ${newdir}/${surfDir}

cp -r ${olddir}/* ${newdir}/${surfDir}/

done

