#!/bin/bash

#This script collects the data generated from PHonon calculation and collects it to a 'save' directory for use in qe2pert.x

#Change the prefix below
PREFIX='graphene'

#ph-collect.sh should be in the work directory of PHonon calculation

echo `date`
echo `pwd`

echo 'PREFIX: ' $PREFIX
echo "Creating a save dir..."
mkdir -p save/${PREFIX}.phsave

PH0_DIR="tmp/_ph0"

echo "Copying prefix.phsave..."
cp ${PH0_DIR}/${PREFIX}.phsave/* save/${PREFIX}.phsave/

echo "Copying dyn files..."
cp ./${PREFIX}.dyn* save/

echo "Copying the dvscf file for the first q-point..."
cp ${PH0_DIR}/${PREFIX}.dvscf1 save/${PREFIX}.dvscf_q1

echo "Copy the dvscf for q-points > 1..."
for q_folder in ${PH0_DIR}/${PREFIX}.q_*; do
   echo $q_folder;
   NQ=`echo $q_folder | awk -F_ '{print $NF}'`;
   cp ${PH0_DIR}/${PREFIX}.q_${NQ}/${PREFIX}.dvscf1 save/${PREFIX}.dvscf_q${NQ}
done

echo "Done!"
echo `date`
