#!/bin/bash


intMatfile=testWhole_46c-IntMat_3_1.dat




#### assign all the values


pref=${1}
inSuffix=${2}
Tp=${3}
ref=${4}

suffix=_${targSuffix}_${inSuffix}_

### consider looping over following
Tsim=0.0001
# sgacdel=0.0001
sigdel=0.001
eps_scale=8
sfdf=0.001

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

for lmn in 0.1 #0.5
do

for eps_mem_mult in 0.5
do

for epsij_mult in 0.4 #0.5 0.6 0.75
do


for sgactau in 2  #1 2 #0.2 1
do


for sgacdel in 0.001 0.02 #0.005 0.01 0.02 #0.0001 0.001 0.002 0.005
do


for actscale in 0.0 0.04 0.2 0.4 0.7 1.0 #0.04 0.08 0.12 0.16 0.2 0.25 0.3 0.4 0.7 1.0 #1 #0.3 0.45 0.6 #0.2 # 0.45 0.7 0.95 #0.1 1
do

eps_diffu=`echo "scale=4; ${eps_scale} * ${eps_mem_mult}" | bc`
epsij_scale=`echo "scale=4; ${eps_scale} * ${epsij_mult}" | bc`



dataDir=/Data1/chromosome_modelling/prepSurf_perturb_merlin/


strpt2=_in${inSuffix}_${pref}_Tp${Tp}/46ch_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_as${actscale}/


outdir=prepSurf_perturb_merlin_geom${pref}_in${inSuffix}_Tp${Tp}/
echo ${outdir}
mkdir -p ${outdir}/


outdir=prepSurf_perturb_merlin_geom${pref}_in${inSuffix}_Tp${Tp}/46ch_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_as${actscale}/

echo ${outdir}
mkdir -p ${outdir}/



python3 analyse_IPD_tbs_cross_sim.py ${dataDir} ${strpt2} Conf-Constraints reference_constraints_IPDHiC_${ref}.dat ${outdir}/ ${targSuffix} ${inSuffix} ${pref}




done #as
done #sgacdel
done #sgactau
done #epsij_mult
done #eps_mem_mult
done #lmn

