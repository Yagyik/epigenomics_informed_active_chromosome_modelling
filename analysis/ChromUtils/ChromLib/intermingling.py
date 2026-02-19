import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial.distance import cdist


def getDecompMat(config): #,intMat,paramslist):
    """ compute the decompaction of patches from surface of chrom centre for each
        we skip patchpointer info, instead just looking at patch distances.
        bookkeeping of which chrom pairs to look at etc are taken care of later
        """
    nchrom = config.nchrom # 46
    npatch = config.npatch # 23
    nato = nchrom+nchrom*npatch
    
    decompMat = np.zeros([nchrom,npatch+1],np.float64)
    for i in range(nchrom):
        patchlisti = config.chromarr[i].patchlist
        chrom = config.chromarr[i]
        for j in range(1,patchlisti[0]+1):
            dx = chrom.patchx[j] - chrom.chrcentx
            dy = chrom.patchy[j] - chrom.chrcenty
            dz = chrom.patchz[j] - chrom.chrcentz
            patchdist = np.sqrt(dx*dx + dy*dy + dz*dz)

##            ppoint = chrom.patchpointer[i][j] ### takes values from 1-23 for chromtype, 1107 for none, 1108 for lmn
##            if ppoint <= npatch:
##                decompMat[i,ppoint] = patchdist/chrom.chrrad
##            elif ppoint == 1107:
##                decompMat[i,
            decompMat[i,j] = patchdist/(chrom.chrrad + chrom.patchrad[j])
        
##        print("chr p dist rad",chrom.chrid,j,patchdist,chrom.chrrad)

    return decompMat




def getIMMat(config, intMat, paramslist):
    """
    Computes pairwise distances between domains, applies sigmoid activation, 
    and stores intermingling scores in a structured format.

    Parameters:
    - snapshot: Object containing chromosome data, where snapshot.chromosomes is a list of chromosome objects
    - interaction_matrix: dict { (chrom1, domain1, chrom2, domain2): strength } - Encoded interactions
    - threshold: float - Threshold for classifying interactions
    - domain_radius: float - Characteristic radius of a domain

    Returns:
    - intermingling_matrix: np.ndarray of shape (46, 46, 23, 23, 2) with scores
    - df_interactions: Pandas DataFrame for categorized interactions
    """

    sepcut = paramslist[0]
    threshold = paramslist[1]
    nchrom = config.nchrom
    npatch = config.npatch
    nato = nchrom+nchrom*npatch
    hicAllMat = np.zeros([nchrom,nchrom,npatch+1,npatch+1,2],np.float64)


    interaction_list = []

    # Extract domain positions from snapshot
    chromosome_data = {}
    # for chrom in config.chromarr:
    #     for domain_idx in chrom.patchlist:
    #         chromosome_data[(chrom.id, domain_idx)] = chrom.patchx[ii]


    for i in range(nchrom-1):
        for ii in range(1,npatch+1):
            pos1x = config.chromarr[i].patchx[ii]
            pos1y = config.chromarr[i].patchy[ii]
            pos1z = config.chromarr[i].patchz[ii]
            r1 = config.chromarr[i].patchrad[ii]
            chromosome_data[(i, ii)] = {
                    'pos': np.array([pos1x, pos1y, pos1z]),
                    'radius': r1
                }
    
    for (chrom1, domain1), data1 in chromosome_data.items():
        for (chrom2, domain2), data2 in chromosome_data.items():
            if chrom1 == chrom2 and domain1 == domain2:
                continue  # Skip self-interactions

            pos1 = data1['pos']
            pos2 = data2['pos']
            radius1 = data1['radius']
            radius2 = data2['radius']
            domain_radius = sepcut*(radius1 + radius2)

            distance = np.linalg.norm(pos1 - pos2)

            # get global indices
            gid1 = config.chromarr[chrom1].patchlist[domain1]
            gid2 = config.chromarr[chrom2].patchlist[domain2]
            interaction_strength = intMat[gid1,gid2]

            # Classify interaction
            interaction_type = "positive" if interaction_strength >= threshold else "negative"

            # Compute sigmoid activation: 1 if close, 0 if beyond 2 * domain_radius
            sigmoid_value = 1 / (1 + np.exp((distance - 2 * domain_radius) * 10))

            # Store in intermingling matrix
            hicAllMat[chrom1, chrom2, domain1, domain2, 0] = interaction_strength
            hicAllMat[chrom1, chrom2, domain1, domain2, 1] = sigmoid_value

            # Store in DataFrame-friendly list
            interaction_list.append([chrom1, domain1, chrom2, domain2, distance, interaction_strength, sigmoid_value, interaction_type])

    # Create DataFrame for easy access and filtering
    # df_interactions = pd.DataFrame(interaction_list, columns=["Chrom1", "Domain1", "Chrom2", "Domain2", 
                                                            #   "Distance", "InteractionStrength", 
                                                            #   "InterminglingScore", "Type"])

    return hicAllMat



def getIMMat_old(config,intMat,paramslist):
    """ compute tanh activated intermingling degree between specific pairs of chromosomes
        Return a matrix of intermingling degrees between these pairs nchrom x nchrom rather than
        nchrom x npatch
        """

    sepcut = paramslist[0]
    decay = paramslist[1]
    nchrom = config.nchrom
    npatch = config.npatch
    nato = nchrom+nchrom*npatch
    hicAllMat = np.zeros([nchrom,nchrom,npatch+1,npatch+1],np.float64)

    for i in range(nchrom-1):
        patchlisti = config.chromarr[i].patchlist
        itype = config.chromarr[i].chrtype
            
        for j in range(i+1,nchrom):
            jtype = config.chromarr[j].chrtype
            patchlistj = config.chromarr[j].patchlist
            if itype == jtype:
                continue
            
            for ii in range(1,patchlisti[0]+1):

                pos1x = config.chromarr[i].patchx[ii]
                pos1y = config.chromarr[i].patchy[ii]
                pos1z = config.chromarr[i].patchz[ii]
                point1 = config.chromarr[i].gcount[ii]
                r1 = config.chromarr[i].patchrad[ii]
                
                
                for jj in range(1,patchlistj[0]+1):

                    pos2x = config.chromarr[j].patchx[jj]
                    pos2y = config.chromarr[j].patchy[jj]
                    pos2z = config.chromarr[j].patchz[jj]
                    point2 = config.chromarr[j].gcount[jj]
                    r2 = config.chromarr[j].patchrad[jj]

                    dx = pos2x - pos1x
                    dy = pos2y - pos1y
                    dz = pos2z - pos1z

                    sep = np.sqrt(dx*dx + dy*dy + dz*dz)

                    ### assign this intermingling to the respective homologues
                    hicAllMat[i,j,ii,jj] = (0.5 - 0.5*np.tanh((sep - sepcut)/(decay*sepcut)))

##                    hicmat[itype,jtype] += (0.5 - 0.5*np.tanh((sep - sepcut)/(decay*sepcut)))
##                    hicmat[jtype,itype] += (0.5 - 0.5*np.tanh((sep - sepcut)/(decay*sepcut)))
    return hicAllMat


def correlDecompHiC(config,intMat,paramslist):
    """ makes calls to the decomp mat routine and to the HiCAll routine
        gets us two big matrices -- here we filter and look for pairwise patterns
        """
    decompMat = getDecompMat(config) ### nchrom x (npatch+1)
    hicAllMat = getIMMat(config,intMat,paramslist) ### nchrom x nchrom x (npatch+1) x (npatch+1)
    print(hicAllMat.shape)

    print(np.unique(hicAllMat))

    return decompMat,hicAllMat
    
def calcHiCPatches(config,intMat,paramslist):
    """ compute contacts from patch positions in xyz
    """

    sepcut = paramslist[0]
    decay = paramslist[1]
    nchrom = config.nchrom
    npatch = config.npatch
    nato = nchrom+nchrom*npatch
    pararr = np.linspace(0,nato-1,nato,dtype=np.int64)
    hicmat = np.zeros([int(0.5*nchrom),int(0.5*nchrom)],np.float64)
    parcount = 0
    for i in range(nchrom-1):
        patchlisti = config.chromarr[i].patchlist
        itype = config.chromarr[i].chrtype
            
        for j in range(i+1,nchrom):
            jtype = config.chromarr[j].chrtype
            if itype == jtype:
                continue
            patchlistj = config.chromarr[j].patchlist

            for ii in range(1,patchlisti[0]):

                pos1x = config.chromarr[i].patchx[ii]
                pos1y = config.chromarr[i].patchy[ii]
                pos1z = config.chromarr[i].patchz[ii]
                point1 = config.chromarr[i].gcount[ii]
                r1 = config.chromarr[i].patchrad[ii]
                
                for jj in range(1,patchlistj[0]):

                    pos2x = config.chromarr[j].patchx[jj]
                    pos2y = config.chromarr[j].patchy[jj]
                    pos2z = config.chromarr[j].patchz[jj]
                    point2 = config.chromarr[j].gcount[jj]
                    r2 = config.chromarr[j].patchrad[jj]

                    dx = pos2x - pos1x
                    dy = pos2y - pos1y
                    dz = pos2z - pos1z
##                    print(i,j,itype,jtype,ii,jj,point1,point2)

##                    if intMat[point1,point2] == 0:
##                        continue

                    sep = np.sqrt(dx*dx + dy*dy + dz*dz)

##                    hicmat[itype,jtype] += 1.0/(r1**2+r2**2)*(0.5 - 0.5*np.tanh((sep - sepcut)/(decay*sepcut)))
##                    hicmat[jtype,itype] += 1.0/(r1**2+r2**2)*(0.5 - 0.5*np.tanh((sep - sepcut)/(decay*sepcut)))
                    hicmat[itype,jtype] += (0.5 - 0.5*np.tanh((sep - sepcut)/(decay*sepcut)))
                    hicmat[jtype,itype] += (0.5 - 0.5*np.tanh((sep - sepcut)/(decay*sepcut)))

                    

    return hicmat        
                    
#### include subroutines to filter in other contacts and to 
    


##
##def calcHiCHulls(config):
##    """ compute contacts from hull intersections
##    """
##
##
##
##
##
def compileHiC(HiC_series):

    nmat = HiC_series.shape[0]

    avgHiC = np.mean(HiC_series,axis=0)

    return avgHiC
    



def plotHiC(HiCMat,parlist,outprefix):
    fig,axs=plt.subplots(nrows=1, ncols=1,sharex=True,sharey=True,figsize=(15,10),gridspec_kw={'width_ratios': [0.8]})
    ##
    ##cmax = np.max(np.log(HiCMat))
    ##cmin = 1
    cmax = np.max(HiCMat)
    cmin = np.min(HiCMat)
    print(cmin,cmax,100)
    xaxis = list(range(0,23,1))
    yaxis = list(range(0,23,1))
##    cbticklabels=list(np.arange(1000,50000,5000))
    cbticklabels = list(np.linspace(cmin,10,10))
##    print(cmin,cmax)
    im1 = axs.imshow(HiCMat[:-1,:-1],vmin=cmin,vmax=10,aspect='auto', cmap='Blues')
    cb = fig.colorbar(im1)
    cb.ax.set_yticklabels([f'{x:.2f}' for x in cbticklabels],size=20)
    cb.set_label("Contact frequency", labelpad=-1, size=25)
##    print(xaxis,yaxis)
    axs.set_xticks(xaxis)
    axs.set_yticks(yaxis)
    axs.set_xticklabels(xaxis, fontsize=20)
    axs.set_yticklabels(yaxis, fontsize=20)
    outdir = outprefix + parlist[1]

    sepcut = parlist[3]
    decay = parlist[4]
    plotname=outdir + '/HiCMat_cl_'+str(sepcut)+'_'+str(decay)+'.png'
##    plt.show()
    plt.savefig(plotname)


    filename=outdir + '/HiCMat_cl_'+str(sepcut)+'_'+str(decay)+'.dat'

    print(filename)
    f=open(filename,"w+") # put proper output file name with temp and fblamb here

    for i in range(HiCMat.shape[0]-1):
        for j in range(HiCMat.shape[1]-1):
            f.write("%d %d %4.4f\n" % (i+1,j+1,HiCMat[i,j]))
    f.close()















    
