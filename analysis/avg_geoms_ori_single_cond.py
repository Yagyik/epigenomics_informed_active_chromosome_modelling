import numpy as np
import matplotlib.pyplot as plt
import time
import math
import sys
from scipy.optimize import curve_fit

plt.rcParams["figure.figsize"] = (15, 10)

def linear_fit(x,a,b):

    return a*x + b

strp1=sys.argv[1]
strp2=sys.argv[2]
nset=int(sys.argv[3])
insuffix=sys.argv[4]
pref=sys.argv[5]
Tp=sys.argv[6]
outdir=sys.argv[7]
fact=float(sys.argv[8])
logbase=4



lengthlist=[24.8 ,24.3 ,19.8 ,19.1 ,18.1 ,17.1 ,15.9 ,14.6, 14.1, 13.6, 13.5,
            13.4, 11.5, 10.7, 10.3, 9, 8.1, 7.8, 5.9, 6.3, 4.8, 5.1, 15.5,
            24.8,24.3 ,19.8 ,19.1 ,18.1 ,17.1 ,15.9 ,14.6, 14.1, 13.6, 13.5,
            13.4, 11.5, 10.7, 10.3, 9, 8.1, 7.8, 5.9, 6.3, 4.8, 5.1, 15.5]

ecfraclist=[0.60,0.52,0.63,0.63,0.59,0.64,0.51,0.66,0.62,0.596,0.65,
    0.44,0.49,0.62,0.72,0.67,0.73,0.67,0.61,0.54,0.646,0.42,0.5,
    0.60,0.52,0.63,0.63,0.59,0.64,0.51,0.66,0.62,0.596,0.65,
    0.44,0.49,0.62,0.72,0.67,0.73,0.67,0.61,0.54,0.646,0.42,0.5]


explengthlist=[24.3,19.1,14.6,11.5,10.7,6.3]
expecfraclist=[0.52,0.63,0.66,0.49,0.62,0.54]
expidlist=[2,4,8,13,14,20]

expRglist=[3.10793471445247, 2.82976842052184,
                2.8714241285949 , 2.35809288836717,
                2.43914060391052, 2.78441863921879]

expkappasqlist=[0.331649822034014, 0.244570725050393,
                     0.19792293250514, 0.238979172798409,
                     0.227756935517807, 0.224796702887166]

expr1list = [2.61,2.29,2.26,1.90,1.96,2.23]
expr2list = [1.47,1.34,1.43,1.16,1.18,1.33]
expr3list = [0.82,0.97,1.03,0.78,0.85,0.99]

explengtharr = np.array(explengthlist)
expr1arr = np.array(expr1list)
expr2arr = np.array(expr2list)
expr3arr = np.array(expr3list)
expecfracarr = np.array(expecfraclist)
expkappasqarr = np.array(expkappasqlist)


### read, open geometries files 
##geomfile.write("%d %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n"
##                   %(chrid,lengthlist[chrid],ecfraclist[chrid],np.mean(r1mat,axis=0)[i],np.std(r1mat,axis=0)[i],
##                     np.mean(r2mat,axis=0)[i],np.std(r2mat,axis=0)[i],np.mean(r3mat,axis=0)[i],np.std(r3mat,axis=0)[i],
##                     np.mean(kappasqmat,axis=0)[i],np.std(kappasqmat,axis=0)[i],
##                     np.mean(1-centmat,axis=0)[i],np.std(1-centmat,axis=0)[i],
##                     np.mean(decompmat,axis=0)[i],np.std(decompmat,axis=0)[i]))


oporientlist=[]

filename=strp1+'/OutPS1_in'+insuffix+'_Tp'+Tp+'/'+strp2+'/geom_dir/geometries_hull.dat'
infile = open(filename,"r");
lines=infile.readlines()
infile.close()
flength = len(lines)

lengtharr = np.array(lengthlist)
ecfracarr = np.array(ecfraclist)
nchrom=len(lengthlist)
r1mat=np.zeros([nset,nchrom],np.float64)
r3mat=np.zeros([nset,nchrom],np.float64)
r2mat=np.zeros([nset,nchrom],np.float64)
kappasqmat=np.zeros([nset,nchrom],np.float64)
centmat=np.zeros([nset,nchrom],np.float64)
decompmat=np.zeros([nset,nchrom],np.float64)

r1stmat=np.zeros([nset,nchrom],np.float64)
r3stmat=np.zeros([nset,nchrom],np.float64)
r2stmat=np.zeros([nset,nchrom],np.float64)
kappasqstmat=np.zeros([nset,nchrom],np.float64)
centstmat=np.zeros([nset,nchrom],np.float64)
decompstmat=np.zeros([nset,nchrom],np.float64)



for i in range(1,nset+1):
    ## get dir
    filename=strp1+'/OutPS'+str(i)+'_in'+insuffix+'_Tp'+Tp+'/'+strp2+'/geom_dir/geometries_hull.dat'
    infile = open(filename,"r")
    lines=infile.readlines()
    infile.close()
    linecount=0

    for line in lines:
        tokens=line.split()

        chrid= int(tokens[0])
        if chrid != linecount:
            print("some pathology",chrid,linecount)
            exit(0)
            
        chrlen = float(tokens[1])
        if chrlen != lengthlist[chrid]:
            print("more pathology!",chrlen,lengthlist[chrid],lengthlist[linecount])
            exit(0)
        ecfrac = float(tokens[2])
        r1mat[i-1,chrid] = float(tokens[3])
        r1stmat[i-1,chrid] = float(tokens[4])
        r2mat[i-1,chrid] = float(tokens[5])
        r2stmat[i-1,chrid] = float(tokens[6])
        r3mat[i-1,chrid] = float(tokens[7])
        r3stmat[i-1,chrid] = float(tokens[8])
        kappasqmat[i-1,chrid] = float(tokens[9])
        kappasqstmat[i-1,chrid] = float(tokens[10])
        centmat[i-1,chrid] = 1-float(tokens[11])
        
##                centmat[i-1,chrid,0] = 1-float(tokens[11])
##                centmat[i-1,chrid,1] = 1-float(tokens[12])
        centstmat[i-1,chrid] = float(tokens[12])
        decompmat[i-1,chrid] = float(tokens[13])
        decompstmat[i-1,chrid] = float(tokens[14])

        linecount +=1



#### write decomp vs centrality to file
print(len(strp2),len(strp2[:-12]),strp2[:-12])
decompopfile = outdir+'/decompaction_'+strp2[:66]+'.dat'

mean_decomp_v_l = np.mean(decompmat,axis=0)
std_decomp_v_l = np.mean(decompstmat,axis=0)

#### do straight line fit
popt,pocv = curve_fit(linear_fit,lengtharr,mean_decomp_v_l,p0=[-1,2],nan_policy='omit',check_finite=False)

if fact==0: ## new file
    f=open(decompopfile,"w+")
else:
    f=open(decompopfile,"a")
print(fact,popt,np.mean(mean_decomp_v_l),np.mean(std_decomp_v_l))
f.write("%4.2f %4.12f %4.12f %4.12f %4.12f\n" % (fact,popt[0],popt[1],np.mean(mean_decomp_v_l),np.mean(std_decomp_v_l)))
f.close()

### read orientations
for i in range(1,nset+1):
    filename=strp1+'/OutPS'+str(i)+'_in'+insuffix+'_Tp'+Tp+'/'+strp2+'/geom_dir/orientations_hull.dat'
    infile = open(filename,"r");
    lines=infile.readlines()
    infile.close()
    linecount=0
    print("orientations lines",len(lines))

    for line in lines:
        tokens=line.split()

##        if round(float(tokens[0])) != float(tokens[0]): ## not an int

##                "%d %4.4f %4.4f %4.4f\n" %(chrid,lengthlist[chrid],ecfraclist[chrid],orient))
        try:
            chrid = int(tokens[0])
        except:
            print(tokens,"skipping")
            continue
            
        chrlen = float(tokens[1])
        ecfrac = float(tokens[2])
        oporientlist.append(float(tokens[3]))
        linecount += 1
        
#### write average orientations        
orientsopfile = outdir+'/orientations_'+strp2+'.dat'
outfile = open(orientsopfile,"w")
for i in range(len(oporientlist)):
    outfile.write("%f \n" %(oporientlist[i]))
outfile.close()

maxr = np.max(r1mat)
    
pairs=[[lengtharr,np.zeros([lengtharr.shape[0]]),r1mat,r1stmat,
        'Chrom length [10 MBp]','R1','R1_v_L_',np.array(explengthlist),np.array(expr1list),(0.5,maxr)],
       [lengtharr,np.zeros([lengtharr.shape[0]]),r2mat,r2stmat,
        'Chrom length [10 MBp]','R2','R2_v_L_',np.array(explengthlist),np.array(expr2list),(0.5,maxr)],
       [lengtharr,np.zeros([lengtharr.shape[0]]),r3mat,r3stmat,
        'Chrom length [10 MBp]','R3','R3_v_L_',np.array(explengthlist),np.array(expr3list),(0.5,maxr)],
       [ecfracarr,np.zeros([lengtharr.shape[0]]),kappasqmat,kappasqstmat,
        'EC fraction','shape anisotropy','shape_aniso_v_ec_frac_',expecfracarr,expkappasqarr],
       [lengtharr,np.zeros([lengtharr.shape[0]]),kappasqmat,kappasqstmat,
        'Chrom length [10 MBp]','shape anisotropy','shape_aniso_v_L_',np.array(explengthlist),expkappasqarr],
       [lengtharr,np.zeros([lengtharr.shape[0]]),centmat,centstmat,
        'Chrom length [10 MBp]','centrality','centrality_v_L_'],
       [centmat,centstmat,decompmat,decompstmat,
        'centrality','decompaction','decomp_v_centrality_'],
       [ecfracarr,np.zeros([lengtharr.shape[0]]),centmat,centstmat,
        'EC fraction','centrality','centrality_v_ec_frac_'],
       [lengtharr,np.zeros([lengtharr.shape[0]]),decompmat,decompstmat,
        'Chrom length [10 MBp]','decompaction','decomp_v_L_']]
### now plot the average
for pair in pairs:
##        print(pair[0].shape,pair[2])
##        print(pair[1].shape,pair[3])
    if len(pair[0].shape)==2:
        xdata = np.mean(pair[0],axis=0)              
        xerr = np.sqrt(np.mean(pair[1]**2,axis=0))
        
    elif len(pair[0].shape)==1:
        xdata = pair[0]
        xerr = pair[1] #np.zeros([pair[0].shape[0]],np.float64)

    if len(pair[2].shape)==2:
        ydata = np.mean(pair[2],axis=0)
        yerr = np.sqrt(np.mean(pair[3]**2,axis=0))
        print(yerr.shape,pair[3].shape)
##        yerr = np.std(pair[1],axis=0)
    elif len(pair[2].shape)==1:
        ydata = pair[2]
        yerr = pair[3] #np.zeros([pair[2].shape[0]],np.float64)

    if pair[5] == 'decompaction' or pair[5] == 'shape anisotropy':
        yerr = yerr/nset
        print("yerr",pair[5],yerr,ydata,strp2)
    print(pair[0].shape,len(pair[0].shape))
        
##    ydata = np.mean(pair[1],axis=0)
##    
##    yerr = np.std(pair[1],axis=0)
    xaxstr = pair[4]
    yaxstr = pair[5]
    plotnamestr = outdir+pair[6]+strp2   #+'.dat'

    plt.clf()
    plt.figure(figsize=[10,8])
    plt.errorbar(xdata,ydata,xerr=xerr,yerr=yerr,fmt='o',markersize=5,color='black')
    spx = (xdata-np.mean(xdata))/np.std(xdata)
    spy = (ydata-np.mean(ydata))/np.std(ydata)
    cov = np.cov(spx,spy)
    lamb,vec = np.linalg.eig(cov)
    pearson = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
    plotslope = np.sign(pearson)*abs(vec[0,1]/vec[0,0])*np.std(ydata)/np.std(xdata)
    plotlinex = np.linspace(min(xdata)-0.5*min(xdata),1.5*max(xdata),100)
    yshift = np.mean(ydata) - plotslope*np.mean(xdata)
    
    plt.plot(plotlinex,plotlinex*plotslope+yshift,lw=3,color='black')
    plt.xlabel(xaxstr,fontsize=30)
    plt.ylabel(yaxstr,fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)

    # Make plot borders (spines) thicker
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(3)

    # Make tick marks more prominent
    ax.tick_params(axis='both', which='both', length=10, width=2, direction='in', top=True, right=True)

    rmse = np.sqrt(np.mean((ydata - (yshift+xdata*plotslope))**2))

    # Choose figure text position based on slope sign
    # Define possible annotation positions: (xpos, ypos1, ypos2, ypos3)
    annotation_positions = {
        "top_left":    (0.2, 0.89, 0.83, 0.77),
        "top_right":   (0.70, 0.89, 0.83, 0.77),
        "bottom_left": (0.2, 0.25, 0.19, 0.13),
        "bottom_right":(0.70, 0.25, 0.19, 0.13)
    }

    # Helper to check if data is in a region
    def data_in_region(x, y, region):
        ax = plt.gca()
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        xspan = xlim[1] - xlim[0]
        yspan = ylim[1] - ylim[0]
        # Define region bounds in axes fraction
        if region == "top_left":
            xbounds = (xlim[0], xlim[0] + 0.3 * xspan)
            ybounds = (ylim[1] - 0.3 * yspan, ylim[1])
        elif region == "top_right":
            xbounds = (xlim[1] - 0.3 * xspan, xlim[1])
            ybounds = (ylim[1] - 0.3 * yspan, ylim[1])
        elif region == "bottom_left":
            xbounds = (xlim[0], xlim[0] + 0.3 * xspan)
            ybounds = (ylim[0], ylim[0] + 0.3 * yspan)
        elif region == "bottom_right":
            xbounds = (xlim[1] - 0.3 * xspan, xlim[1])
            ybounds = (ylim[0], ylim[0] + 0.3 * yspan)
        else:
            return False
        # Check if any data point falls in this region
        return np.any((x >= xbounds[0]) & (x <= xbounds[1]) & (y >= ybounds[0]) & (y <= ybounds[1]))

    # Choose annotation positions for model and experimental data
    # Try to avoid overlap with data
    model_x, model_y = xdata, ydata
    exp_x, exp_y = None, None
    if len(pair) > 7:
        exp_x, exp_y = pair[7], pair[8]

    # Find available regions for model and experimental annotations
    regions = ["top_left", "top_right", "bottom_left", "bottom_right"]
    model_region = None
    exp_region = None

    # Check which regions are free of model data
    free_regions_model = [r for r in regions if not data_in_region(model_x, model_y, r)]
    # For experimental data, if present
    free_regions_exp = regions.copy()
    if exp_x is not None and exp_y is not None:
        free_regions_exp = [r for r in regions if not data_in_region(exp_x, exp_y, r)]

    
    # Assign regions, prefer top corners for model, bottom for exp, but avoid overlap
    if free_regions_model:
        if "top_left" in free_regions_model:
            model_region = "top_left"
        elif "top_right" in free_regions_model:
            model_region = "top_right"
        else:
            model_region = free_regions_model[0]
    else:
        model_region = "top_left"  # fallback

    if exp_x is not None and exp_y is not None:
        # Avoid putting exp in same region as model
        possible_exp_regions = [r for r in free_regions_exp if r != model_region]
        if possible_exp_regions:
            if "bottom_right" in possible_exp_regions:
                exp_region = "bottom_right"
            elif "bottom_left" in possible_exp_regions:
                exp_region = "bottom_left"
            else:
                exp_region = possible_exp_regions[0]
        else:
            # If all regions have data, just pick one not used by model
            exp_region = [r for r in regions if r != model_region][0]
    if pair[5] == 'decompaction':
        model_region = "top_right"  # Force model data to top right for decompaction plots
    elif pair[5] == 'shape anisotropy':
        model_region = "top_right"
        exp_region = "bottom_left"  # Force experimental data to bottom left for shape anisotropy plots

    # Place model annotation
    xpos, ypos1, ypos2, ypos3 = annotation_positions[model_region]
    plt.figtext(xpos, ypos1, f'sPCC = {pearson:.2f}', size=20)
    plt.figtext(xpos, ypos2, f'sRMSE = {rmse:.2f}', size=20)
    plt.figtext(xpos, ypos3, f'sScSlope = {np.mean(model_x)/np.mean(model_y)*plotslope:.2f}', size=20)

    # Place experimental annotation if present
    if exp_x is not None and exp_y is not None:
        plt.scatter(exp_x, exp_y, s=80, color='red')
        cov = np.cov(exp_x, exp_y)
        lamb, vec = np.linalg.eig(cov)
        pearson_exp = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
        plotslope_exp = np.sign(pearson_exp)*abs(vec[0,1]/vec[0,0])
        plotlinex = np.linspace(min(exp_x)-0.5*min(exp_x), 1.5*max(exp_x), 100)
        yshift = np.mean(exp_y) - plotslope_exp*np.mean(exp_x)
        plt.plot(plotlinex, plotlinex*plotslope_exp + yshift, lw=3, color='red')
        rmse_exp = np.sqrt(np.mean((exp_y - (yshift + exp_x*plotslope_exp))**2))

        # Use chosen region for experimental annotation
        xpos_e, ypos1_e, ypos2_e, ypos3_e = annotation_positions[exp_region]
        plt.figtext(xpos_e, ypos1_e, f'ePCC = {pearson_exp:.2f}', size=20, color='red')
        plt.figtext(xpos_e, ypos2_e, f'eRMSE = {rmse_exp:.2f}', size=20, color='red')
        plt.figtext(xpos_e, ypos3_e, f'eScSlope = {np.mean(exp_x)/np.mean(exp_y)*plotslope_exp:.2f}', size=20, color='red')

        print(len(pair), pair[-1])

        if len(pair) > 9:
            plt.ylim(pair[-1][0], pair[-1][1])

    plt.tight_layout()
    plotname = plotnamestr + '.png'
    print(plotname)
    plt.savefig(plotname)

    plotname = plotnamestr + '.eps'
    print(plotname)
    plt.savefig(plotname,format='eps',dpi=300)


