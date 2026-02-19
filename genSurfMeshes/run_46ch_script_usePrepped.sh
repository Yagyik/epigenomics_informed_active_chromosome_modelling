#!/bin/bash


outDir=$1
inSuffix=$2
targSuffix=$3

epsblob=$4
epsij_scale=$5



### other energy scales constrained
kspring=`echo "scale=4; ${epsblob} / 8" | bc`
epshertz=`echo "scale=4; ${epsblob} * 4.0" | bc`

### epsij scale not constrained (upper bound at epsblob)


# epsij_scale=`echo "scale=4; ${epsblob} / 2" | bc`



##### surface spring -- tauinverse hardcoded to maintain stress-geometry relationship
surfself=$6  ## `echo "scale=4; ${epsblob} * 0.5" | bc`
surftauinv=`echo "scale=4; ${surfself} * 1" | bc`

### modulatable
surfcross=`echo "scale=4; ${surfself} * 0.125" | bc`
# surfDiffu=`echo "scale=4; ${surftauinv} * 0.01" | bc`


sigdeli=$7
sigdelf=$8
sigactdeli=$9
sigactdelf=${10}
Temi=${11}
Temf=${12}

lmni=${13}
lmnf=${14}
lmn_scalei=`echo "scale=4; ${epsij_scale} * ${lmni}" | bc` ### upperbound at eps_blob
lmn_scalef=`echo "scale=4; ${epsij_scale} * ${lmnf}" | bc` ### upperbound at eps_blob
# lmn_scalei=`echo "scale=4; ${lmn}" | bc` ### hard-coded -- no upper bound (lower bound 2 because creates a repulsive hump at large-ish r)
# lmn_scalef=`echo "scale=4; ${lmn}" | bc` ### hard-coded -- no upper bound (lower bound 2 because creates a repulsive hump at large-ish r)


surfDiffu=${15}
sgactau=${16}
actscalei=${17}
actscalef=${18}

echo ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18}
echo ${sgactau}


pref=${19}

start=${20}
start_anneal=${21}
simlength=${22}
Tp=${23}




echo ${start} ${start_anneal} ${simlength} ${Tp}
## tau_theta and surfspring combine to produce fluctuations -- high tau_theta decreases responsivity to stress -- high spring also decreases responsivity to stress

echo ${epsblob} ${epsij_scale} ${surfself}
echo ${sigdeli} ${sigdelf} ${sigactdeli} ${sigactdeli}
echo ${Temi} ${Temf} ${lmni} ${lmnf}
echo ${surfDiffu} ${sgactau} ${surfcross} ${surftauinv}
echo ${pref}

## generate seed for targSuffix, sgacdel, sgactau, actscale
seed=`echo "scale=2; ${targSuffix} * 1000+ ${actscalei} * 100 + ${sigactdeli} * 100 + 5 * ${sgactau}" | bc`
echo ${seed}

if [ ${pref} == ell ]
then
echo ${pref} ell
ea=9
eb=6
ec=2.5
# surfDir=Set_${targSuffix}_in2_${pref}/prepSurf_perturb/
fi

if [ ${pref} == sph ]
then
echo ${pref} sph
ea=5.5
eb=5.5
ec=5.5
# surfDir=Set_${targSuffix}_insph1_${pref}/prepSurf_perturb/
fi

# if [ ${pref} == ell ]
# then
# ella=9
# ellb=6
# ellc=2.5
# fi
#
# if [ ${pref} == sph ]
# then
# ella=5.5
# ellb=5.5
# ellc=5.5
# fi





## conditional -- default, with GD annealing a parameter
maindir=Set_${targSuffix}_in${inSuffix}_${pref}_Tp${Tp}/
startsdir=${inSuffix}_${pref}_Tp${Tp}_starts/
surfcyc=1000

### old -- w/o annealing parameter -- and different charac surfcyc
if [ ${Tp} == 0 ]
then
maindir=Set_${targSuffix}_in${inSuffix}_${pref}/
startsdir=${inSuffix}_${pref}_starts/
surfcyc=500
fi

echo ${pref} ${ea} ${eb} ${ec}
mkdir -p ${maindir}
mkdir -p ${maindir}/${outDir}
echo ${maindir}/${outDir}

mkdir -p ParaFiles
rm ParaFiles/Para_test_46c_whole_${targSuffix}.par

echo "Nchrom		    46" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "nPatchTot		1058" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "Temi				${Temi}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "Temf				${Temf}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "dt				0.002" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "gamma			10" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "rotgamma		10" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "chromAct		0" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "actscalei		${actscalei}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "actscalef		${actscalef}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "rho				0.1" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "pressure		0" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "init_read		${start}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "startGeom		0" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "seed			${seed}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "Label			Conf" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "TSfact			2" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "eqRun			${start_anneal}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "totRun			${simlength}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "thermo			200" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "dump			1000" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "ellai			${ea}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "ellbi			${eb}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "ellci			${ec}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "ellaf			9" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "ellcf			2.5" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "dellwin			5000" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "eps_blob		${epsblob}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "springscale		${kspring}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "espij_scale   ${epsij_scale}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "eps_hertz		${epshertz}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "nSurfBasicPoints		2500" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "surfNeighfact		1.75" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "nSurfCycPoints		${surfcyc}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "refineRounds		5" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "delta_Ar_th		0.1" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "surfPhi			0.35" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "surfSelf		${surfself}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "surfCrossi		${surfcross}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "surfCrossf		${surfcross}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "surfDiffu		${surfDiffu}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "lamin_scalei		${lmn_scalei}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "lamin_scalef		${lmn_scalef}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "taup			100" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "sigdeli          ${sigdeli}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "sigdelf          ${sigdelf}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "nGrid           1000" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "tau_theta     ${sgactau}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "sigactdeli       ${sigactdeli}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "sigactdeli       ${sigactdelf}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "surftauinv      ${surftauinv}" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "surfrestore     1.0" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "hic_cut     0.5" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par
echo "hic_decay     0.25" >> ParaFiles/Para_test_46c_whole_${targSuffix}.par


cp ParaFiles/Para_test_46c_whole_${targSuffix}.par ${maindir}/${outDir}/


# startsdir=${inSuffix}_${pref}_Tp${Tp}_starts/
# startsdir=${inSuffix}_${pref}_starts/
surfDir=prepSurf_perturb

cp ${maindir}/${surfDir}/Conf-finalfixedSurfFile_ea${ea}_eb${eb}_ec${ec}_5.dat ${maindir}/${outDir}/
cp ${maindir}/${surfDir}/Conf-SurfCyc_ea${ea}_eb${eb}_ec${ec}_5.dat ${maindir}/${outDir}/


# nohup ./prepSurf_perturb.exe ParaFiles/Para_test_46c_whole_${targSuffix}.par ${startsdir}/testWhole_46c-SimStart_${inSuffix}_${targSuffix}.dat ${startsdir}/testWhole_46c-IntMat_${inSuffix}_${targSuffix}.dat Set_${targSuffix}_in${inSuffix}_${pref}/${outDir}/ > Set_${targSuffix}_in${inSuffix}_${pref}/${outDir}/prepsurf_log.dat &

./mdChrom_batch_gr.exe ParaFiles/Para_test_46c_whole_${targSuffix}.par ${startsdir}/testWhole_46c-SimStart_${inSuffix}_${targSuffix}.dat ${startsdir}/testWhole_46c-IntMat_${inSuffix}_${targSuffix}.dat ${maindir}/${outDir}/

# nohup ./mdChrom_batch_gr.exe ParaFiles/Para_test_46c_whole_${targSuffix}.par ${startsdir}/testWhole_46c-SimStart_${inSuffix}_${targSuffix}.dat ${startsdir}/testWhole_46c-IntMat_${inSuffix}_${targSuffix}.dat ${maindir}/${outDir}/ > ${maindir}/${outDir}/log.dat &
