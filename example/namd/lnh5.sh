#!/bin/bash

npools=2 # No. of processes in each part of ephmat calculation.
nkparts=7 # No. of k-points in each part of ephmat calculation.
prefix='graphene'
pertpath="../perturbo/ephmat"

ndigit=$(echo -n ${nkparts} | wc -c)
fullpath=$(cd $pertpath; pwd)

mkdir -p h5files
cd h5files

for ((ipart=1; ipart<=$nkparts; ipart++))
do
    for ((ipool=1; ipool<=$npools; ipool++))
    do
        jpart=$[$ipart * $npools - $npools + $ipool]
        ifile=${prefix}_ephmat_p${ipool}.h5
        jfile=${prefix}_ephmat_p${jpart}.h5
        pnum=`printf "%0${ndigit}d" ${ipart}`
        ln -sf ${fullpath}/${pnum}P/$ifile $jfile
    done
done

cd ..
