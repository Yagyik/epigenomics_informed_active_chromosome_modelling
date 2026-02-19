#!/bin/bash



# dataDir=/Data1/chromosome_modelling/stability_tests/46chrom_reorg_batch_4_1/
# dataDir=/Data1/chromosome_modelling/stability_tests/46chrom_reorg_batch_4_10/

# dataDir=../AnaCode_chrom_model/46chrom_batch_geom1/
# intMatfile=testWhole_46c-IntMat_4.dat
intMatfile=testWhole_46c-IntMat_3_1.dat
# midfix=${1}
# suffix=${2} ## _smpatch_ _smlam_ or _



# dataDir=/Data1/chromosome_modelling/stability_tests/46chrom_bat4/46chrom_reorg_batch_4_${suffix}/
# #
# dataDir=/Data1/chromosome_modelling/stability_tests/46chrom_${midfix}_batch_${suffix}/
# outdir=${midfix}_batch_${suffix}
# dataDir=../46chrom_reorg_batch_4_${suffix}/ ## suffix = _reorg_

# dataDir=/Data1/chromosome_modelling/stability_tests/46chrom_sph_batch_${suffix}/
# outdir=sph_batch_1_${suffix}

#### assign all the values

targSuffix=${1}
pref=${2}
inSuffix=${3}

fhull=${4}
Tp=${5}
ref=${6}

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


mainpath=/Data3/AnaCode_chromosome_modelling/prepSurf_perturb_merlin_geom${pref}


for lmn in 0.1 #0.5
do

for eps_mem_mult in 0.5
do

for epsij_mult in 0.4 0.5 0.6 0.75
do


for sgactau in 2  #1 2 #0.2 1
do


for sgacdel in 0.001 0.005 0.01 0.02 #0.0001 0.001 0.002 0.005
do


for actscale in 0.0 0.04 0.08 0.12 0.16 0.2 0.25 0.3 0.4 0.7 1.0 #1 #0.3 0.45 0.6 #0.2 # 0.45 0.7 0.95 #0.1 1
do

eps_diffu=`echo "scale=4; ${eps_scale} * ${eps_mem_mult}" | bc`
epsij_scale=`echo "scale=4; ${eps_scale} * ${epsij_mult}" | bc`


# dataDir=/home/yagyik/Dropbox/chromosome_modelling/prepSurf_perturb_reprog/Set_${targSuffix}_in${inSuffix}_${pref}_Tp${Tp}/46ch_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_as${actscale}/

# dataDir=/home/goswam_y/Dropbox/chromosome_modelling/prepSurf_perturb_reprog/Set_${targSuffix}_in${inSuffix}_${pref}_Tp${Tp}/46ch_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_as${actscale}/

dataDir=/Data1/chromosome_modelling/prepSurf_perturb_merlin/Set_${targSuffix}_in${inSuffix}_${pref}_Tp${Tp}/46ch_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_as${actscale}/

outdir=${mainpath}
echo ${outdir}
mkdir -p ${outdir}/

outdir=${mainpath}/OutPS${targSuffix}_in${inSuffix}_Tp${Tp}/
mkdir -p ${outdir}/


outdir=${mainpath}/OutPS${targSuffix}_in${inSuffix}_Tp${Tp}/46ch_T${Tsim}_sga${sgacdel}_sgd${sigdel}_epb${eps_scale}_epa${epsij_scale}_epm${eps_diffu}_lmn${lmn}_sfdf${sfdf}_sgat${sgactau}_as${actscale}/









echo ${outdir}
mkdir -p ${outdir}/

# python3 analyse_IPD_tbs.py ${dataDir} Conf-Constraints reference_constraints_IPDHiC_${ref}.dat ${outdir}/ ${targSuffix} ${inSuffix} ${pref}


mkdir -p SimParamFiles
mkdir -p AnaParamFiles

rm SimParamFiles/simparams${suffix}geom.json
rm AnaParamFiles/simparams${suffix}geom.json

rm SimParamFiles/simparams${suffix}msd.json
rm AnaParamFiles/simparams${suffix}msd.json




echo ${outdir}
mkdir -p ${outdir}/
mkdir -p ${outdir}/HiC_dir
mkdir -p ${outdir}/msd_dir
mkdir -p ${outdir}/velAC_dir
mkdir -p ${outdir}/geom_dir
mkdir -p ${outdir}/packing_dir
mkdir -p ${outdir}/IMperco_dir
mkdir -p ${outdir}/Intensity_dir
mkdir -p ${outdir}/Chperco_dir

# rm ${outdir}/geom_dir/orientations*


echo "{" >> SimParamFiles/simparams${suffix}msd.json
echo "    \"simlength\": [100000,10000,2500000]," >> SimParamFiles/simparams${suffix}msd.json
echo "    \"dt\": 0.001," >> SimParamFiles/simparams${suffix}msd.json
echo "    \"nsurf\": 3500," >> SimParamFiles/simparams${suffix}msd.json
echo "    \"nchrom\": 46," >> SimParamFiles/simparams${suffix}msd.json
echo "    \"npatch\": 23," >> SimParamFiles/simparams${suffix}msd.json
echo "    \"confdir\": \"${dataDir}\"," >> SimParamFiles/simparams${suffix}msd.json
echo "    \"infileprefix\": \"Conf-Anadump_\"," >> SimParamFiles/simparams${suffix}msd.json
echo "    \"intmatfile\": \"${intMatfile}\"" >> SimParamFiles/simparams${suffix}msd.json
echo "}"  >> SimParamFiles/simparams${suffix}msd.json



echo "{" >> AnaParamFiles/simparams${suffix}msd.json
echo "    \"outprefix\": \"${outdir}/\"," >> AnaParamFiles/simparams${suffix}msd.json
echo "    \"intermingling\": [0,\"HiC_dir/\",1,0.5,0.5]," >> AnaParamFiles/simparams${suffix}msd.json
# echo "    \"msd\": [1,\"msd_dir\"]," >> AnaParamFiles/simparams${suffix}msd.json
echo "    \"msd\": [1,\"msd_dir/\"]," >> AnaParamFiles/simparams${suffix}msd.json
echo "    \"gca_greenkubo\": [0,\"velAC_dir/\"]," >> AnaParamFiles/simparams${suffix}msd.json
echo "    \"geom_props\": [0,\"geom_dir/\",${fhull}]," >> AnaParamFiles/simparams${suffix}msd.json
echo "    \"packing\": [0,\"packing_dir/\"]," >> AnaParamFiles/simparams${suffix}msd.json
echo "    \"IMpercolation\": [0,\"IMperco_dir/\",0.1,${ella},${ellb},${ellc}]," >> AnaParamFiles/simparams${suffix}msd.json
echo "    \"IntensityMap\": [0,\"Intensity_dir/\",0.2,${ella},${ellb},${ellc},8]," >> AnaParamFiles/simparams${suffix}msd.json
echo "    \"Channelpercolation\": [0,\"Chperco_dir/\",0.2,${ella},${ellb},${ellc},8,1.0]," >> AnaParamFiles/simparams${suffix}msd.json
echo "    \"tau\": 0.035" >> AnaParamFiles/simparams${suffix}msd.json
echo "}" >> AnaParamFiles/simparams${suffix}msd.json


echo "{" >> SimParamFiles/simparams${suffix}geom.json
echo "    \"simlength\": [10000,200000,4900000]," >> SimParamFiles/simparams${suffix}geom.json
echo "    \"dt\": 0.001," >> SimParamFiles/simparams${suffix}geom.json
echo "    \"nsurf\": 3500," >> SimParamFiles/simparams${suffix}geom.json
echo "    \"nchrom\": 46," >> SimParamFiles/simparams${suffix}geom.json
echo "    \"npatch\": 23," >> SimParamFiles/simparams${suffix}geom.json
echo "    \"confdir\": \"${dataDir}\"," >> SimParamFiles/simparams${suffix}geom.json
echo "    \"infileprefix\": \"Conf-Anadump_\"," >> SimParamFiles/simparams${suffix}geom.json
echo "    \"intmatfile\": \"${intMatfile}\"" >> SimParamFiles/simparams${suffix}geom.json
echo "}"  >> SimParamFiles/simparams${suffix}geom.json

echo "{" >> AnaParamFiles/simparams${suffix}geom.json
echo "    \"outprefix\": \"${outdir}/\"," >> AnaParamFiles/simparams${suffix}geom.json
echo "    \"intermingling\": [0,\"HiC_dir/\",1,0.5,0.5]," >> AnaParamFiles/simparams${suffix}geom.json
# echo "    \"msd\": [1,\"msd_dir\"]," >> AnaParamFiles/simparams${suffix}geom.json
echo "    \"msd\": [0,\"msd_dir/\"]," >> AnaParamFiles/simparams${suffix}geom.json
echo "    \"gca_greenkubo\": [0,\"velAC_dir/\"]," >> AnaParamFiles/simparams${suffix}geom.json
echo "    \"geom_props\": [0,\"geom_dir/\",${fhull}]," >> AnaParamFiles/simparams${suffix}geom.json
echo "    \"packing\": [0,\"packing_dir/\",0.1,20000]," >> AnaParamFiles/simparams${suffix}geom.json
echo "    \"IMpercolation\": [0,\"IMperco_dir/\",0.1,${ella},${ellb},${ellc}]," >> AnaParamFiles/simparams${suffix}geom.json
echo "    \"IntensityMap\": [1,\"Intensity_dir/\",0.4,${ella},${ellb},${ellc},8]," >> AnaParamFiles/simparams${suffix}geom.json
echo "    \"Channelpercolation\": [0,\"Chperco_dir/\",0.4,${ella},${ellb},${ellc},8,1.0]," >> AnaParamFiles/simparams${suffix}geom.json
echo "    \"tau\": 0.035" >> AnaParamFiles/simparams${suffix}geom.json
echo "}" >> AnaParamFiles/simparams${suffix}geom.json

#
# python3 run_conf_analysis.py AnaParamFiles/simparams${suffix}msd.json SimParamFiles/simparams${suffix}msd.json
python3 run_conf_analysis.py AnaParamFiles/simparams${suffix}geom.json SimParamFiles/simparams${suffix}geom.json



echo ${outdir}

done
done
done
done
done
done
