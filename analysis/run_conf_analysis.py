import numpy as np
import scipy
import matplotlib.pyplot as plt
import sys
import subprocess
import pickle as pkl


from ChromUtils.Structures import *
from ChromUtils.utils import *
from ChromUtils.conf_reader import *
from ChromUtils.ChromLib.intermingling import *
from ChromUtils.ChromLib.geometrics import *
from ChromUtils.ChromLib.percolation import *
from ChromUtils.TSeriesLib.MSD import *
from ChromUtils.TSeriesLib.Coactivity_GK import *


def main():
    global anaparams
    global simparams
    anaparams_file = sys.argv[1]
    simparams_file = sys.argv[2]

    anaparams = parse_file(anaparams_file)

    outprefix = anaparams["outprefix"]

    
    simparams = parse_file(simparams_file)

    
    start = simparams["simlength"][0]
    step = simparams["simlength"][1]
    stop = simparams["simlength"][2]
    nsurf = simparams["nsurf"]
    nchrom = simparams["nchrom"]
    npatch = simparams["npatch"]
    nato = nsurf+nchrom+nchrom*2*npatch ## 1600 + 8 + 8x2x4 OR 2000 + 46 + 46x2x23
    

    indir = simparams["confdir"]
    infileprefix = simparams["infileprefix"]

    

    nfiles = 1 + int((stop - start)/step)
    print(nfiles," files")

    configs = np.empty(nfiles,dtype=object)

    for i in range(nfiles):
        t = start + i*step
        filename = indir + "/" + infileprefix + str(t) + ".dat"

        configs[i] = read_config(filename,t,nsurf,nchrom,npatch)

##    if anaparams["intermingling"][0] == 1:
    intMat = readIntMat(simparams["intmatfile"],configs[0])

    print("using int mat\n",intMat.shape)
##    print(intMat[intMat>0])

##    print(configs)
##    print(configs[0])
##    print(configs[0].chromarr[0].patchlist)

##    pf_arr = np.empty(nfiles,dtype=float64)
##    geom_arr = np.empty(nfiles,dtype=float64)
##    orient_arr = np.empty(nfiles,dtype=float64)
##    disparr = np.empty(nfiles,dtype=float64)
##    area_arr = np.empty(nfiles,dtype=float64)

    pf_series = np.array([],np.float64)
    geom_series = np.array([],np.float64)
    disp_series = np.array([],np.float64)
    AC_series = np.array([],np.float64)
    HiC_series = np.array([],np.float64)
    decomp_series = np.array([],np.float64)
    allHiC_series = np.array([],np.float64)
    
    area_series = np.array([],np.float64)
    intMat_series = np.array([],np.float64)
    ov_series = np.array([],np.float64)
    ov_tseries = np.array([],np.float64)
    
    for i in range(nfiles):

        t = start + i*step

        if anaparams["packing"][0] == 1:

            res = calcPackingFrac(configs[i],anaparams["packing"][2],anaparams["packing"][3],outprefix,anaparams["packing"][1])
            pf = np.array([x for x in res],np.float64)

            pf_series = stackDim(pf,pf_series)

        if anaparams["geom_props"][0] == 1:
##            print(anaparams["geom_props"])

            geom_arr = calcChromGeometries(configs[i],bool(anaparams["geom_props"][2]))
            ## R1,R2,R3,Rg,kappasq,asphere,cent,decomp,orient for config i

            geom_series = stackDim(geom_arr,geom_series)

            


        if anaparams["msd"][0] == 1:

            disp = calcMSD(configs[0],configs[i])
##            print(disp)
##            print(disp.shape)
            disp_series = stackDim(disp,disp_series)
            ov_mat,ov_ratio = calcOverlaps(configs[i])
            ov_series = stackDim(ov_mat,ov_series)
            ov_tseries = stackDim(np.array([t,ov_ratio]),ov_tseries)
            if i==0 or i == nfiles-1:              
                shx = ov_mat.shape[0]
                shy = ov_mat.shape[1]
                ov_plot = np.reshape(ov_mat,(1,shx,shy))
                plotOvMat(ov_plot,ov_tseries,i,anaparams["msd"],outprefix)

        if anaparams["gca_greenkubo"][0] == 1:

            AC = ACTwoTime(configs[0],configs[i],intMat)
            AC_series = stackDim(AC,AC_series)


        if anaparams["intermingling"][0] == 1:
            
            if anaparams["intermingling"][2] == 1:
                HiCMat = calcHiCPatches(configs[i],intMat,anaparams["intermingling"][3:])
##            if anaparams["intermingling"][1] == 2:
##                HiCMat = calcHiCHulls(configs[i],intMat,
##                                        anaparams["intermingling"][2],
##                                        anaparams["intermingling"][3])

            HiC_series = stackDim(HiCMat,HiC_series)

            decompMat,allHiCMat = correlDecompHiC(configs[i],intMat,anaparams["intermingling"][3:])
            

            decomp_series = stackDim(decompMat,decomp_series)
            allHiC_series = stackDim(allHiCMat,allHiC_series)

##        if anaparams["area_fluctuations"] == 1:
##
##            props = areaProps(configs[i])
##            area_arr[i] = np.array([prop for prop in props])
        if anaparams["IntensityMap"][0] == 1:
            print("calc intensity\n",i,t)
            intensityMat = intensityMaps(configs[i],anaparams["IntensityMap"],makefig=False)

            intMat_series = stackDim(intensityMat,intMat_series)

        if anaparams["Channelpercolation"][0] == 1:
            print("find channels\n")
            channelPercolation(configs[i],anaparams["Channelpercolation"])

        if anaparams["IMpercolation"][0] == 1:
            print("finding rigid structures\n")
            connectedGraph(configs[i],intMat,anaparams["IMpercolation"])
             


    if anaparams["msd"][0] == 1:
##        print(disp_series)
        msd,chrmsd,patchmsd,lmnmsd = splitMSD(configs[0],disp_series)
        timeax = anaparams["tau"]*simparams["dt"]*np.arange(start-start,stop+step-start,step)
##        print(timeax)

        plotMSD(msd,chrmsd,patchmsd,lmnmsd,timeax,anaparams["msd"],outprefix)

        plotOvMat(ov_series,ov_tseries,"all",anaparams["msd"],outprefix)

        
    if anaparams["gca_greenkubo"][0] == 1:
        timeax = anaparams["tau"]*simparams["dt"]*np.arange(start-start,stop+step-start,step)
        plotAC(AC_series,timeax,anaparams["gca_greenkubo"],outprefix)
        

    if anaparams["intermingling"][0] == 1:
        avgHiC = compileHiC(HiC_series)

        plotHiC(avgHiC,anaparams["intermingling"],outprefix)

        pickledecomp = outprefix + anaparams["intermingling"][1] + "/pickleDecomp"
        fileObj = open(pickledecomp,'wb')

        pkl.dump(decomp_series,fileObj)
        fileObj.close()

        pickleAllHiC = outprefix + anaparams["intermingling"][1] + "/pickleAllHiC"
        fileObj = open(pickleAllHiC,'wb')

        print(allHiC_series.shape)
        pkl.dump(allHiC_series,fileObj)
        fileObj.close()

    if anaparams["geom_props"][0] == 1:

        plotGeom(configs[0],geom_series,anaparams["geom_props"],outprefix,bool(anaparams["geom_props"][2]))

    if anaparams["IntensityMap"][0] == 1:

        avgIntensity = compileIntensity(intMat_series)
        plotIntensity(avgIntensity,anaparams["IntensityMap"],outprefix)

        
    
if __name__ == "__main__":
    main()


    
