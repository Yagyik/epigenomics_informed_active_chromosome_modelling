import numpy as np
import matplotlib.pyplot as plt

import ChromUtils.utils as utils
from ChromUtils.Structures import *
from ChromUtils.ChromLib.hullprops import *

def calcGeometry(chrom,ella,ellb,ellc,fromHull: bool = False):
    """ calculates geometric properties based on the supplied chrom points
    """

    patchlist = chrom.patchlist   
    if fromHull==True:
        hull,collect = computeHullSingleChrom(chrom,toplot=False)
##        denom = collect.shape[0]
##        denom = patchlist[0]+1
        denom = len(hull.vertices)
##        print("DENOM",denom,len(hull.vertices))
        xc = collect[hull.vertices]
        
        print("denoms",patchlist[0]+1,xc.shape[0],collect.shape[0])
##        denom = xc.shape[0]
        ##### calculate gyration tensor and geometric properties
        cent = np.mean(xc,axis=0)
##        cent = np.mean(collect,axis=0)
        print(np.mean(xc,axis=0),np.mean(collect,axis=0))
##        print("chrom centre from hull",chrom.chrid,cent,collect.shape[0],pt.shape[0])
    if fromHull==False:
        npatch = patchlist[0]
        collect = np.zeros([npatch+1,5],np.float64)
        
        
##        denom = collect.shape[0]
        denom = patchlist[0] + 1
        
        for j in range(patchlist[0]+1):
            if j==0:
                radspam = chrom.chrrad
                currposx = chrom.chrcentx
                currposy = chrom.chrcenty
                currposz = chrom.chrcentz
                tc = chrom.chrtype
                
            else:
    ##            print(j,patchrad[chrid][j],patchx[chrid][j],patchy[chrid][j],patchz[chrid][j])
                radspam = chrom.patchrad[j]
                currposx = chrom.patchx[j]
                currposy = chrom.patchy[j]
                currposz = chrom.patchz[j]
                tc = chrom.chrtype
            collect[j,:] = np.array([currposx,currposy,currposz,radspam,tc])
        cent = np.mean(collect,axis=0)
        xc = np.copy(collect)
##        print("chrom centre from points",chrom.chrid,cent,collect.shape[0],pt.shape[0])
    ### calc geometric properties for the collection
    gyr_tensor = np.zeros([3,3],np.float64)
    for i in range(3):
##            print(i,chrid,cent[i])
        for j in range(3):
##                print(i,j,chrid,cent[j])
            for pt in xc:
                gyr_tensor[i,j] += (1.0/denom)*(pt[i] - cent[i])*(pt[j]-cent[j]) ## treated as assembly
##                gyr_tensor[i,j] += (pt[i] - cent[i])*(pt[j]-cent[j]) ### treated as one object
     

##    print(gyr_tensor)
    eigval,eigvec = np.linalg.eigh(gyr_tensor)
    D = np.diag(eigval)
    P_inv = np.linalg.inv(eigvec)
    gyr_tensor_diag = np.dot(P_inv,np.dot(gyr_tensor,eigvec))
##    print(gyr_tensor_diag)
##    print(D,"\n\n")

    lambsq = eigval
    print("old lamb",np.sqrt(lambsq))
##    print(collect.shape[0],np.sum(lambsq))
##    Rg = np.sqrt((1.0/collect.shape[0])*np.sum(lambsq)) ### treated as assembly
    Rg = np.sqrt(np.sum(lambsq)/denom) ## treated as one object
    kappasq = 1.5*(lambsq[0]**2 + lambsq[1]**2 + lambsq[2]**2)/(np.sum(lambsq)**2) - 0.5
    asphere = 1.5*lambsq[2] - 0.5*Rg**2


    ### centrality and decompaction
    cent = 1 - np.sqrt((chrom.chrcentx/ella)**2 + (chrom.chrcenty/ellb)**2 + (chrom.chrcentz/ellc)**2 )
##    print("chrom centre",chrid,np.mean(chromcollect,axis=0),centchrom,cent,chromcollect.shape[0],ptcollect.shape[0],xc.shape[0])
    decomp = 0
    for j in range(1,patchlist[0]+1):
        dx = chrom.patchx[j] - chrom.chrcentx
        dy = chrom.patchy[j] - chrom.chrcenty
        dz = chrom.patchz[j] - chrom.chrcentz
        patchdist = np.sqrt(dx*dx + dy*dy + dz*dz)
##        print("chr p dist rad",chrom.chrid,j,patchdist,chrom.chrrad)
        decomp += patchdist/chrom.chrrad
    decomp /= patchlist[0]

    orient = (180.0/np.pi)*abs(np.arccos(np.dot(eigvec[2],np.array([1,0,0]))/np.linalg.norm(eigvec[2])))

    ### factor of sqrt(3) goes into lambda -> ea,eb,ec conversion based on ref.
    ## Influence of Solvent Quality on Conformations of Crowded Polymers - Wyatt J. Davis and Alan R. Denton
    return chrom.timestep,np.sqrt(3*lambsq[2]),np.sqrt(3*lambsq[1]),np.sqrt(3*lambsq[0]),Rg,kappasq,asphere,cent,decomp,orient
    


def calcChromGeometries(config,fromHull: bool = False):
    """ subroutine to aggregate all chromosome geometric properties
    save as an array
    """
##    print("printing from hull?? ", fromHull)
    geom_arr = np.array([],np.float64)

    for i in range(config.nchrom):
        res = calcGeometry(config.chromarr[i],config.ella,config.ellb,config.ellc,fromHull)
        
##        geom_arr = np.array([x for x in res]) ## R1,R2,R3,Rg,kappasq,asphere,cent,decomp,orient

        geom_arr = utils.stackDim(np.array([x for x in res]),geom_arr)

    return geom_arr


def plotGeom(config,geom_mat,parlist,outprefix,fromHull: bool = False):

    refquants = referenceQuants()

    lengthlist = refquants.lengthlist
    ecfraclist = refquants.ecfraclist

    explengthlist = refquants.explengthlist
    expecfraclist = refquants.expecfraclist
    expidlist = refquants.expidlist
    expRglist = refquants.expRglist
    expkappasqlist = refquants.expkappasqlist

    expr1list = refquants.expr1list
    expr2list = refquants.expr2list
    expr3list = refquants.expr3list


    explengtharr = np.array(explengthlist)
    expr1arr = np.array(expr1list)
    expr2arr = np.array(expr2list)
    expr3arr = np.array(expr3list)
    expecfracarr = np.array(expecfraclist)
    expkappasqarr = np.array(expkappasqlist)

    

    r1mat = geom_mat[:,:,1]
    r2mat = geom_mat[:,:,2]
    r3mat = geom_mat[:,:,3]

    print(r1mat.shape,geom_mat.shape)

    lengtharr = np.array(lengthlist)
    ecfracarr = np.array(ecfraclist)
    kappasqmat = geom_mat[:,:,5]
    centmat = geom_mat[:,:,7]
    decompmat = geom_mat[:,:,8]
    orientmat = geom_mat[:,:,9]

    outdir = outprefix + parlist[1]
    if fromHull== True:
        geomfilename = outdir+'geometries_hull.dat'
        orientsfilename = outdir+'orientations_hull.dat'
    else:
        geomfilename = outdir+'geometries_points.dat'
        orientsfilename = outdir+'orientations_points.dat'
    print(geomfilename)
    geomfile= open(geomfilename,"w")

    orientsfile = open(orientsfilename,"w")
    

    ### store geom info with mean and stdev
    for i in range(lengtharr.shape[0]):
        chrid = config.chromarr[i].chrid
##        print(chrid,lengthlist[chrid],ecfraclist[chrid],np.mean(r1mat,axis=0),np.std(r1mat,axis=0),
##                     np.mean(r2mat,axis=0),np.std(r2mat,axis=0),np.mean(r3mat,axis=0),np.std(r3mat,axis=0),
##                     np.mean(kappasqmat,axis=0),np.std(kappasqmat,axis=0),
##                     np.mean(1-centmat,axis=0),np.std(1-centmat,axis=0),
##                     np.mean(decompmat,axis=0),np.std(decompmat,axis=0),
##                     np.mean(orientmat,axis=0),np.std(orientmat,axis=0))
        geomfile.write("%d %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n"
                   %(chrid,lengthlist[chrid],ecfraclist[chrid],np.mean(r1mat,axis=0)[i],np.std(r1mat,axis=0)[i],
                     np.mean(r2mat,axis=0)[i],np.std(r2mat,axis=0)[i],np.mean(r3mat,axis=0)[i],np.std(r3mat,axis=0)[i],
                     np.mean(kappasqmat,axis=0)[i],np.std(kappasqmat,axis=0)[i],
                     np.mean(1-centmat,axis=0)[i],np.std(1-centmat,axis=0)[i],
                     np.mean(decompmat,axis=0)[i],np.std(decompmat,axis=0)[i]))
##                     np.mean(orientmat,axis=0)[i],np.std(orientmat,axis=0)[i]))
        for j in range(geom_mat.shape[0]):
            orient = orientmat[j,i]
            if orient < 90.0 and orient > 0.0:
                orientsfile.write("%d %4.4f %4.4f %4.4f\n" %(chrid,lengthlist[chrid],ecfraclist[chrid],orient))
            elif orient > 90.0:
                orientsfile.write("%d %4.4f %4.4f %4.4f\n" %(chrid,lengthlist[chrid],ecfraclist[chrid],180 - orient))
            elif orient < 0.0:
                orientsfile.write("%d %4.4f %4.4f %4.4f\n" %(chrid,lengthlist[chrid],ecfraclist[chrid],-orient))

    geomfile.close()
    orientsfile.close()



    ##### store orientation list
    
    
    
    
    ### make pairs


##    r1mean = np.mean(r1mat,axis=0)
##    r1std = np.std(r1mat,axis=0)
##
##    print(r1mat)
##    print(r1mean.shape)
##    print(r1std.shape)
##
##    plt.errorbar(lengtharr,r1mean,yerr=r1std,fmt='o')
##    plt.scatter(explengthlist,expr1list,marker='o',s=50,color='red')
##    plt.show()

    
    pairs=[[lengtharr,r1mat,'Chrom length','R1','R1_v_L_',explengtharr,expr1arr],
       [lengtharr,r2mat,'Chrom length','R2','R2_v_L_',explengtharr,expr2arr],
       [lengtharr,r3mat,'Chrom length','R3','R3_v_L_',explengtharr,expr3arr],
       [ecfracarr,kappasqmat,'EC fraction','shape anisotropy','shape_aniso_v_ec_frac_',expecfracarr,expkappasqarr],
        [lengtharr,kappasqmat,'Chrom length','shape anisotropy','shape_aniso_v_L_',explengtharr,expkappasqarr],
       [lengtharr,centmat,'Chrom length','centrality','centrality_v_L_'],
       [centmat,decompmat,'centrality','decompaction','decomp_v_centrality_'],
       [ecfracarr,centmat,'EC fraction','centrality','centrality_v_ec_frac_'],
        [lengtharr,decompmat,'Chrom length','decompaction','decomp_v_L_']]
    
    correls=[]
    rmses=[]

    for pair in pairs:
##        print(pair[0].shape,pair[2])
##        print(pair[1].shape,pair[3])
        if len(pair[0].shape)==2:
            xdata = np.mean(pair[0],axis=0)
            xerr = np.std(pair[0],axis=0)
        elif len(pair[0].shape)==1:
            xdata = pair[0]
            xerr = np.zeros([pair[0].shape[0]],np.float64)

        if len(pair[1].shape)==2:
            ydata = np.mean(pair[1],axis=0)
            yerr = np.std(pair[1],axis=0)
        elif len(pair[1].shape)==1:
            ydata = pair[1]
            yerr = np.zeros([pair[1].shape[0]],np.float64)

        print(pair[0].shape,len(pair[0].shape))
            
        ydata = np.mean(pair[1],axis=0)
        
        yerr = np.std(pair[1],axis=0)
        xaxstr = pair[2]
        yaxstr = pair[3]
        plotnamestr = pair[4]

        plt.clf()
        plt.figure(figsize=(10,8))


        plt.errorbar(xdata,ydata,xerr=xerr,yerr=yerr,fmt='o',markersize=5,color='black')
        cov = np.cov(xdata,ydata)
        lamb,vec = np.linalg.eig(cov)
        pearson = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
        plotslope = np.sign(pearson)*abs(vec[0,1]/vec[0,0])
        plotlinex = np.linspace(min(xdata)-0.5*min(xdata),1.5*max(xdata),100)
        yshift = np.mean(ydata) - plotslope*np.mean(xdata)
        
        plt.plot(plotlinex,plotlinex*plotslope+yshift,lw=3,color='black')
    ##    plt.scatter(xdata,yshift+xdata*plotslope,s=60,color='black')
        plt.xlabel(xaxstr,fontsize=30)
        plt.ylabel(yaxstr,fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)

        rmse = np.sqrt(np.mean((ydata - (yshift+xdata*plotslope))**2))

        plt.figtext(0.20, 0.89, f'sPCC = {pearson:.2f}',size=18)
        plt.figtext(0.20, 0.83, f'sRMSE = {rmse:.2f}',size=18)
        plt.figtext(0.20, 0.77, f'sScSlope = {np.mean(xdata)/np.mean(ydata)*plotslope:.2f}',size=18)

        correls.append(pearson)
        rmses.append(rmse)

        #### now overlay experimental data points

        if len(pair) > 5:

            xdata = pair[5]
            ydata = pair[6]
            plt.scatter(xdata,ydata,s=80,color='red')
            cov = np.cov(xdata,ydata)
            lamb,vec = np.linalg.eig(cov)
            pearson = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
            plotslope = np.sign(pearson)*abs(vec[0,1]/vec[0,0])
            plotlinex = np.linspace(min(xdata)-0.5*min(xdata),1.5*max(xdata),100)
            yshift = np.mean(ydata) - plotslope*np.mean(xdata)
            
            plt.plot(plotlinex,plotlinex*plotslope+yshift,lw=3,color='red')
            rmse = np.sqrt(np.mean((ydata - (yshift+xdata*plotslope))**2))
        ##    print(ydata,yshift+xdata*plotslope)
        ##    print(ydata - (yshift+xdata*plotslope),(ydata - (yshift+xdata*plotslope))**2)
        ##    print(cov)
        ##    print(vec)
        ##    print(lamb,pearson,plotslope,rmse)
            plt.figtext(0.70, 0.89, f'ePCC = {pearson:.2f}',size=18)
            plt.figtext(0.70, 0.83, f'eRMSE = {rmse:.2f}',size=18)
            plt.figtext(0.70, 0.77, f'eScSlope = {np.mean(xdata)/np.mean(ydata)*plotslope:.2f}',size=18)

        

        plt.tight_layout()
        plt.rcParams["figure.figsize"] = (25,20)
        if fromHull == True:
            plotname = outdir + plotnamestr + 'hull.png'
        else:
            plotname = outdir + plotnamestr + 'points.png'
        print(plotname)
        plt.savefig(plotname)
    ##    plt.show()
        ## append the entries (total 4)
        
    

    

    
        





    
