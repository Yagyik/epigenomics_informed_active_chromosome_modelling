#!/bin/bash

# targSuffix=${1}

pref=${1}
inSuffix=${2}
Tp=${3}
# fhull=${4}

# suffix=_${targSuffix}_${inSuffix}_

nset=8

Tsim=0.0001
# sgacdel=0.0001
sigdel=0.001
eps_scale=8
# eps_diffu=4.0
# lmn=0.1
sfdf=0.001
# sgactau=2
# actscale=0.1

if [ ${pref} == ell ]
then
ella=9
ellb=6
ellc=2.5
fi

if [ ${pref} == sph ]
then
ella=5.5
ellb=5.5
ellc=5.5
fi

for lmn in 0.1 #0.001 #0.5
do

for sgactau in 2 #0.2 1
do

eps_diffu=`echo "scale=4; ${eps_scale} * 0.5" | bc`

for sgacdel in 0.001 0.005 0.01 0.02 #0.0001 0.001 0.002 0.005
do

for epsij_mult in 0.6 #0.75 #0.5 0.4 0.6 #0.6 0.75 #0.5 0.6
do

epsij_scale=`echo "scale=4; ${eps_scale} * ${epsij_mult} " | bc`
echo ${eps_scale} ${epsij_scale} ${eps_diffu} ${sigdel} ${sgacdel} ${Tsim} ${lmn} ${sfdf} ${sgactau}

for actscale in 0.0 0.08 0.2 0.3 #0.0 0.04 0.08 0.12 0.16 0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.9 1.0 #1.0 ### baseline upper bound is 0.2, multiply divide to get diff numbers for now
do

inDir=46ch_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_as${actscale}

dirp1=/Data3/AnaCode_chromosome_modelling/prepSurf_perturb_merlin_geom${pref}/
dirp2=46ch_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_as${actscale} #_as${actscale}/

outdir=averaged_quantities_${inSuffix}_${pref}_Tp${Tp}
# outdir=averaged_quantities_LMN
mkdir -p ${outdir}
mkdir -p ${outdir}/IPD
mkdir -p ${outdir}/MSD
mkdir -p ${outdir}/HiC
mkdir -p ${outdir}/Geoms


# python3 avg_msd_single_cond.py ${dirp1} ${dirp2} ${nset} 80 ${inSuffix} ${pref} ${Tp} ${outdir}/MSD/
# python3 avg_ipd_pcc_single_cond.py ${dirp1} ${dirp2} ${nset} 80 ${inSuffix} ${pref} ${Tp} ${outdir}/IPD/
python3 avg_geoms_ori_single_cond.py ${dirp1} ${dirp2} ${nset} ${inSuffix} ${pref} ${Tp} ${outdir}/Geoms/ ${actscale}
# python3 avg_HiC_single_cond.py ${dirp1} ${dirp2} ${nset} 0.5 0.5 ${inSuffix} ${pref} ${Tp} ${outdir}/HiC/
sleep 2
done
done
done
done
done
