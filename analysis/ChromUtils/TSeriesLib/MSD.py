import numpy as np
import matplotlib.pyplot as plt


def calcOverlaps(currconf):

    npar = currconf.nchrom
    ovmat = np.zeros([npar,npar],np.float64)
    tot_ov = 0
    counts = 0
    for i in range(npar-1):
        pos1x = currconf.chromarr[i].chrcentx
        pos1y = currconf.chromarr[i].chrcenty
        pos1z = currconf.chromarr[i].chrcentz

        rad1 = currconf.chromarr[i].chrrad

        for j in range(i+1,npar):
            
            pos2x = currconf.chromarr[j].chrcentx
            pos2y = currconf.chromarr[j].chrcenty
            pos2z = currconf.chromarr[j].chrcentz

            

            rad2 = currconf.chromarr[j].chrrad

##            print(pos1x,pos1y,pos1z,pos2x,pos2y,pos2z,rad1,rad2)

            sep = np.sqrt((pos2x-pos1x)**2 + (pos2y-pos1y)**2 + (pos2z-pos1z)**2 )

            overlap = sep - rad1 - rad2

            if overlap < 0:
                x_int = (sep**2 - rad1**2 + rad2**2)/(2*sep)
                d1 = x_int
                d2 = sep - x_int
                h1 = (rad1 - d1) #/rad1
                h2 = (rad2 - d2) #/rad2
                ov1 = h1/rad1
                ov2 = h2/rad2

                ovmat[i,j] = abs(ov1)
                ovmat[j,i] = abs(ov2)

                tot_ov += abs(ov1) + abs(ov2)
            counts +=1

    

    # print(tot_ov,counts)
    return ovmat,tot_ov/counts


def plotOvMat(ovmat_series,ov_tseries,step,parlist,outprefix):
    
    ovmat = np.mean(ovmat_series,axis=0)
    npar = ovmat.shape[0]
    fig,axs=plt.subplots(nrows=1, ncols=1,sharex=True,sharey=True,figsize=(15,10),gridspec_kw={'width_ratios': [0.8]})
    cmax = np.max(ovmat)
    cmin = np.min(ovmat)
    # print(cmin,cmax,100)
    xaxis = list(range(0,npar,2))
    yaxis = list(range(0,npar,2))
    ##    cbticklabels=list(np.arange(1000,50000,5000))
    cbticklabels = list(np.linspace(cmin,cmax,10))
    ##    print(cmin,cmax)
    im1 = axs.imshow(ovmat,vmin=cmin,vmax=cmax,aspect='auto', cmap='Blues')
    cb = fig.colorbar(im1)
    cb.ax.set_yticklabels([f'{x:.2f}' for x in cbticklabels],size=20)
    cb.set_label("Overlap", labelpad=-1, size=25)
    ##    print(xaxis,yaxis)
    axs.set_xticks(xaxis)
    axs.set_yticks(yaxis)
    axs.set_xticklabels(xaxis, fontsize=20)
    axs.set_yticklabels(yaxis, fontsize=20)
    plotfile = outprefix + parlist[1] + "/Overlaps_MD_"+str(step)+".png"
    plt.savefig(plotfile)

    if step == "all":
        outfilename = outprefix + parlist[1] + "/overlap_t_series.dat"
        f=open(outfilename,"w+") # put proper output file name with temp and fblamb here
        ##f.write("%8d %4.4f\n" % (0,1.0))
        for i in range(1,ov_tseries.shape[0]):
        ##    if denomlogMSD[i]>0 and i > 0:
            f.write("%4.4f %4.8f\n" % (ov_tseries[i,0],ov_tseries[i,1]))
        f.close()
        
        

    

        
        

def calcMSD(initconf,finconf):
    """ subroutine to calculate the displacement wrt initial position
    takes two configs as input
    """
    nato = initconf.nchrom+initconf.nchrom*initconf.npatch
    spam = finconf.nchrom+finconf.nchrom*finconf.npatch
    
    if nato != spam:
        print("config mismatch at ",initconf.timestep,finconf.timestep)
        exit(0)
    disparr = np.zeros([nato,3],np.float64)
    parcount = 0
    for i in range(initconf.nchrom):

        for j in range(initconf.chromarr[i].patchlist[0]+1):
            if j==0:
                oldposx = initconf.chromarr[i].chrcentx
                oldposy = initconf.chromarr[i].chrcenty
                oldposz = initconf.chromarr[i].chrcentz
                currposx = finconf.chromarr[i].chrcentx
                currposy = finconf.chromarr[i].chrcenty
                currposz = finconf.chromarr[i].chrcentz
                
            else:
                oldposx = initconf.chromarr[i].patchx[j]
                oldposy = initconf.chromarr[i].patchy[j]
                oldposz = initconf.chromarr[i].patchz[j]
                currposx = finconf.chromarr[i].patchx[j]
                currposy = finconf.chromarr[i].patchy[j]
                currposz = finconf.chromarr[i].patchz[j]

            disparr[parcount,:] = np.array([currposx - oldposx,currposy - oldposy, currposz - oldposz])
##            print(disparr[parcount])
            parcount+=1


    return disparr

        
        
    



def splitMSD(config,disparr_series):
    """ subroutine to split the displacement array into bead types
    return type-specific arrays for further analysis
    """

    nconfig = disparr_series.shape[0]
    nato = disparr_series.shape[1]

    msd = np.zeros([nconfig],np.float64)
    msdchrcent = np.zeros([nconfig],np.float64)
    msdpatch = np.zeros([nconfig],np.float64)
    msdlamin = np.zeros([nconfig],np.float64)

    spam = config.nchrom+config.nchrom*config.npatch

    if nato != spam:
        print("config mismatch at splitMSD")
        exit(0)

    nchrom = config.nchrom
    npatch = config.npatch

    typearr = np.zeros([nato],np.float64)
    radarr = np.zeros([nato],np.float64)

    parcount=0
    for i in range(config.nchrom):
        patchlist = config.chromarr[i].patchlist
        

        for j in range(patchlist[0] + 1):
            if j==0:
                typearr[parcount] = 0
                radarr[parcount] = config.chromarr[i].chrrad

            elif j > 0 and j < config.chromarr[i].patchlist[0]:
                typearr[parcount] = 1
                radarr[parcount] = config.chromarr[i].patchrad[j]

            else:
                typearr[parcount] = 2
                radarr[parcount] = config.chromarr[i].patchrad[j]
            parcount+=1

##    print(typearr)
    for t in range(nconfig):
        disparr = disparr_series[t,:,:]
        sqdisp = np.square(disparr)
        # print(sqdisp.shape)
##        print(sqdisp)
        sumsqdisp = np.sum(sqdisp,axis=1)
##        print(sumsqdisp[typearr==0].shape,
##              sumsqdisp[typearr==1].shape,
##              sumsqdisp[typearr==2].shape)
        

        msd[t] = np.mean(sumsqdisp,axis=0) #/sumsqdisp.shape[0]
        msdchrcent[t] = np.mean(sumsqdisp[typearr==0],axis=0) #/sumsqdisp[typearr==0].shape[0]
        msdpatch[t] = np.mean(sumsqdisp[typearr==1],axis=0) #/sumsqdisp[typearr==1].shape[0]
        msdlamin[t] = np.mean(sumsqdisp[typearr==2],axis=0) #/sumsqdisp[typearr==2].shape[0]
        
##    print(msd)
##    print(msdchrcent)
##    print(msdpatch)
##    print(msdlamin)
    return msd,msdchrcent,msdpatch,msdlamin
                

def plotMSD(msd,chrmsd,patchmsd,lmnmsd,timeax,parlist,outprefix):

    plt.clf()
    plt.plot(timeax[1:],msd[1:],lw=3,color="black",label="all")
    plt.plot(timeax[1:],chrmsd[1:],lw=3,color="red",label="chrcent")
    plt.plot(timeax[1:],patchmsd[1:],lw=3,color="blue",label="patch")
    plt.plot(timeax[1:],lmnmsd[1:],lw=3,color="purple",label="lmn")
    xaxstr="time [s]"
    yaxstr="MSD [$\mu$ m$^2$]"
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(xaxstr,fontsize=25)
    plt.ylabel(yaxstr,fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(0.0001,100)
    plt.xlim(0.5*timeax[1],5*timeax[-1])
    plt.legend(loc="lower right", prop={'size': 25})
    plt.tight_layout()


    plotname = outprefix + parlist[1] + "/allMSD.png"
    plt.savefig(plotname)

    filename=outprefix + parlist[1] + '/allMSD.dat'
    # print(filename)
    f=open(filename,"w+") # put proper output file name with temp and fblamb here
    ##f.write("%8d %4.4f\n" % (0,1.0))
    for i in range(1,timeax.shape[0]):
    ##    if denomlogMSD[i]>0 and i > 0:
        f.write("%4.4f %4.8f %4.8f %4.8f %4.8f\n" % (timeax[i],msd[i],chrmsd[i],patchmsd[i],lmnmsd[i]))
    f.close()


    plt.clf()
##    plt.plot(timeax[1:],msd[1:],lw=3,color="black",label="all")
    plt.plot(timeax[1:],chrmsd[1:],lw=3,color="red",label="chrcent")
##    plt.plot(timeax[1:],patchmsd[1:],lw=3,color="blue",label="patch")
##    plt.plot(timeax[1:],lmnmsd[1:],lw=3,color="purple",label="lmn")
    xaxstr="time [s]"
    yaxstr="MSD"
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(xaxstr,fontsize=25)
    plt.ylabel(yaxstr,fontsize=25)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylim(0.0001,10)
    plt.legend(loc="lower right",prop={'size': 15})
    plt.tight_layout()

    plotname = outprefix + parlist[1] + "/centMSD.png"
    plt.savefig(plotname)


    plt.clf()
##    plt.plot(timeax[1:],msd[1:],lw=3,color="black",label="all")
##    plt.plot(timeax[1:],chrmsd[1:],lw=3,color="red",label="chrcent")
    plt.plot(timeax[1:],patchmsd[1:],lw=3,color="blue",label="patch")
##    plt.plot(timeax[1:],lmnmsd[1:],lw=3,color="purple",label="lmn")
    xaxstr="time [s]"
    yaxstr="MSD"
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(xaxstr,fontsize=20)
    plt.ylabel(yaxstr,fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylim(0.0001,50)
    plt.legend(loc="lower right", prop={'size': 15})
    plt.tight_layout()

    plotname = outprefix + parlist[1] + "/patchMSD.png"
    plt.savefig(plotname)

##    plt.show()
    
    
        

    

    

    

    
