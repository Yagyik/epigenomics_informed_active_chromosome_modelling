#!/bin/bash

mainpath=/Data1/chromosome_modelling/prepSurf_perturb_merlin/

pref=${1}
insuffix=${2}
Tp=${3}
outDir=phase_IM_merlin_in${insuffix}_pref${pref}_Tp${Tp}/
mkdir -p ${outDir}
for Tsim in 0.0001
do

for sgacdel in 0.001 0.02 #0.01 0.015 0.02 #0.0001 0.001 0.002 0.005 #0.005 0.01 0.02
do


for sigdel in 0.001
do

for eps_scale in 8 #8 16
do

eps_diffu=`echo "scale=4; ${eps_scale} * 0.5" | bc`


for epsij_mult in 0.5 0.75 #0.5 0.6
do

epsij_scale=`echo "scale=4; ${eps_scale} * ${epsij_mult} " | bc`
echo ${eps_scale} ${epsij_scale} ${eps_diffu} ${sigdel} ${sgacdel} ${Tsim} ${lmn} ${sfdf} ${sgactau}


for lmn in 0.1 #0.25 0.5 #20 # max eps_scale, do some factor so (0,1)
do

for sfdf in 0.001 # 0.02
do

for sgactau in 2 #1 # 20
do

rm ${outDir}/avgIM_v_va_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}.dat
echo ${outDir}/avgIM_v_va_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}.dat
rm ${outDir}/thermo_v_va_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_set*.dat

for actscale in 0.0 0.04 0.08 0.12 0.14 0.16 0.2 0.22 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.7 0.75 0.8 0.9 1.0 #1.0 ### baseline upper bound is 0.2, multiply divide to get diff numbers for now
do

inDir=46ch_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_as${actscale}
# outDir=prepSurf_perturb
echo ${inDir}
mkdir -p ${outDir}/${inDir}

for set in 1 2 3 4 5 6 7 8
do

thermofile=${mainpath}/Set_${set}_in${insuffix}_${pref}_Tp${Tp}/${inDir}/Conf-thermo.dat

tail -10000 ${thermofile} > cutthermo.dat

awk -v va=${actscale} '{sum1 +=$2} END {print va" "sum1/NR}' cutthermo.dat >> ${outDir}/thermo_v_va_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_set${set}.dat



done #{set}

echo ${mainpath} ${inDir} ${pref} ${insuffix} ${Tp} 40 4 ${outDir}/${inDir} 50 ${actscale} 8

python3 logBin_avgIM_nostop.py ${mainpath} ${inDir} ${pref} ${insuffix} ${Tp} 40 4 ${outDir}/${inDir} 50 ${actscale} 8



done
## now average over columns of our thermo files

for f in `ls ${outDir}/thermo_v_va_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_set*.dat`
do

awk '{print $2}' $f > $f.En.dat
# awk '{print $3}' $f > $f.KE.dat

done

paste ${outDir}/thermo_v_va_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_set*.dat.En.dat > allEn.dat
# paste ${outpath}/thermo_v_f_Tp${Tp}_taup${taup}_set*.dat.KE.dat > allKE.dat
rm ${outDir}/thermo_v_va_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_set*.dat.En.dat
# rm ${outpath}/thermo_v_f_Tp${Tp}_taup${taup}_set*.dat.KE.dat


# sum over columns

paste -d " " <(awk '{print $1}' ${outDir}/thermo_v_va_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_set1.dat) <(awk '{sum=0; for(i=1;i<=NF;i++) sum+=$i; print sum/NF}' allEn.dat) >> ${outDir}/avgIM_v_va_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}.dat

echo ${outDir}/avgIM_v_va_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}.dat

# rm ${outDir}/thermo_v_va_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_set*.dat


done



done

done

done

done

done

done

done #Tsim
