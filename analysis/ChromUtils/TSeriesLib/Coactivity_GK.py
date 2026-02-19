import numpy as np
import matplotlib.pyplot as plt


""" stuff to calculate a gene-coactivity measure,
looking at relative velocity between ij pairs at different timepoints
then performing an integration over time, to get an integrated auto-correlation
type function, i.e., a Green-Kubo relation, measuring if implied gene expression
activity is similar on the basis of relative motion"""


def calcRelVel(initconf,finconf,intMat):
    """ calculate relative velocities for chosen patches (or all patches)
        then find average and integrate

        """

    nato = initconf.nchrom+initconf.nchrom*initconf.npatch
    spam = finconf.nchrom+finconf.nchrom*finconf.npatch
    
    if nato != spam:
        print("config mismatch at ",initconf.timestep,finconf.timestep)
        exit(0)
    disparr = np.zeros([nato,3],np.float64)
    relvellist = []
    filterlist = []
    parcount = 0
    for i in range(initconf.nchrom-1):
        patchlisti = initconf.chromarr[i].patchlist
        itype = initconf.chromarr[i].chrtype

        for j in range(i+1,initconf.nchrom):
            jtype = initconf.chromarr[j].chrtype
            if itype == jtype:
                continue
            patchlistj = initconf.chromarr[j].patchlist

            for ii in range(1,patchlisti[0]):              
                v10x = initconf.chromarr[i].patchvx[ii]
                v10y = initconf.chromarr[i].patchvy[ii]
                v10z = initconf.chromarr[i].patchvz[ii]
                point1 = initconf.chromarr[i].gcount[ii]
                r1 = initconf.chromarr[i].patchrad[ii]

                v11x = finconf.chromarr[i].patchvx[ii]
                v11y = finconf.chromarr[i].patchvy[ii]
                v11z = finconf.chromarr[i].patchvz[ii]
                
                for jj in range(1,patchlistj[0]):

                    v20x = initconf.chromarr[j].patchvx[jj]
                    v20y = initconf.chromarr[j].patchvy[jj]
                    v20z = initconf.chromarr[j].patchvz[jj]
                    point2 = initconf.chromarr[j].gcount[jj]
                    r2 = initconf.chromarr[j].patchrad[jj]

                    v21x = finconf.chromarr[j].patchvx[jj]
                    v21y = finconf.chromarr[j].patchvy[jj]
                    v21z = finconf.chromarr[j].patchvz[jj]

                    dv0x = v20x - v10x
                    dv0y = v20y - v10y
                    dv0z = v20z - v10z

                    dv1x = v21x - v11x
                    dv1y = v21y - v11y
                    dv1z = v21z - v11z

                    ### append 
                    relvellist.append(dv0x*dv1x + dv0y*dv1y + dv0z*dv1z)
                    ### append 1-0 to a filterlist
                    if intMat[point1,point2] == 0 and point1 != point2: ## only off-diagonals
                        filterlist.append(0)
                    else:
                        filterlist.append(1)

                    ### later can append interaction weight list, or separation-based filter list

    return relvellist,filterlist


def ACTwoTime(initconf,finconf,intMat):

    relvellist,filterlist = calcRelVel(initconf,finconf,intMat)

    ### now do an averaging over relvellist as per needs
    relvelarr = np.array(relvellist)
    mean1 = np.mean(relvelarr)
    mean2 = np.sum(relvelarr * np.array(filterlist))/np.sum(filterlist)

    ### returns a value, parameterised on the time of second conf -- this is stackDim'd
    return np.array([mean1,mean2])


def plotAC(AC_series,timeax,parlist,outprefix):

    ### subroutine that writes the AC time series to a file and plots it
    plt.clf()
    plt.plot(timeax[1:],AC_series[1:,0],lw=3,color="black",label="all")
    plt.plot(timeax[1:],AC_series[1:,1],lw=3,color="red",label="IntMatFiltered")
    xaxstr="time [s]"
    yaxstr="C$_{v_{12}}$(t) [$\mu$ m$^2$ / s$^2$]"
    plt.xscale("log")
##    plt.yscale("log")
    plt.xlabel(xaxstr,fontsize=25)
    plt.ylabel(yaxstr,fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
##    plt.ylim(0.000001,100)
    plt.xlim(0.5*timeax[1],5*timeax[-1])
    plt.legend(loc="lower right", prop={'size': 25})
    plt.tight_layout()


    plotname = outprefix + parlist[1] + "/allVelAC.png"
    plt.savefig(plotname)

    filename=outprefix + parlist[1] + '/allVelAC.dat'
    print(filename)
    f=open(filename,"w+") # put proper output file name with temp and fblamb here
    ##f.write("%8d %4.4f\n" % (0,1.0))
    for i in range(1,timeax.shape[0]):
    ##    if denomlogMSD[i]>0 and i > 0:
        f.write("%4.4f %4.8f %4.8f\n" % (timeax[i],AC_series[i,0],AC_series[i,1]))
    f.close()

    ### the integral needs to be reported -- maybe the time series' can be read in an avg script and put in a diff. file vs va etc.
    ### needs to be 

                    

                    
    
