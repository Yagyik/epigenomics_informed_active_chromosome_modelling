import numpy as np
import scipy
import matplotlib.pyplot as plt
import sys
import subprocess



##indir = '/Data1/chromosome_modelling/stability_tests/46chrom_stable_batch_2_'+dirsuffix+'/'
##indir = '/Data1/chromosome_modelling/stability_tests/46chrom_bat4/46chrom_reorg_batch_4_'+dirsuffix+'/'
##indir = '../46chrom_reorg_batch_4_'+dirsuffix+'/'
##indir='46chrom_geometries_'+dirsuffix+'/'

### first read reference constraints, in particular, c0vec and flag vec



#### redo the path specification and run
##dirsuffix=sys.argv[1]
##outdir=sys.argv[2]
##
####indir = '/Data1/chromosome_modelling/stability_tests/46chrom_reorg_batch_7_'+dirsuffix+'/'
##
####indir = '/Data1/chromosome_modelling/stability_tests/46chrom_sph_batch_sph1_'+dirsuffix+'/'
##
##conf_prefix = 'Conf-Constraints'
##
##reffilename = 'reference_constraints.dat' #indir+conf_prefix+'.dat'

indir=sys.argv[1]
conf_prefix=sys.argv[2]
reffilename=sys.argv[3]
outdir=sys.argv[4]

print("new version")
print(indir)
print(conf_prefix)
print(reffilename)
print(sys.argv)
##fi=int(sys.argv[2])
nsurf=3500
nchrom=46
npatch=23

start=0

vec_tseries = np.array([],np.float64)
tseries = np.array([],np.float64)
nconstraints=1 + int((0.5*(23+1)*(23)) + 23)


opt_vec_tseries = np.array([],np.float64)
opt_tseries = np.array([],np.float64)
print(reffilename)
reffile =open(reffilename,"r")
lines=reffile.readlines()
reffile.close()

for line in lines[:10]:
    tokens = line.split()
##    timeindex = float(tokens[0])
    res = [eval(i) for i in tokens[1:nconstraints]]
##    print(tokens[0])
    if float(tokens[0]) == 1.1:
        c0vec = np.array(res)
        print(len(tokens),len(res),nconstraints)

    elif float(tokens[0]) == 1.2:
        calcflag = np.array(res)

    elif float(tokens[0]) == start:
        opt_vec_tseries = np.array(res)
        opt_tseries = np.array([float(tokens[0])+1])
    elif float(tokens[0]) > start :
        opt_vec_tseries = np.vstack((opt_vec_tseries,np.array(res)))
        opt_tseries = np.append(tseries,float(tokens[0])+1)

## convert reference constraints to matrix

print(calcflag)

##calcflag = np.ones([len(res)],np.int64)

orig_ind = np.array([])
for i in range(calcflag.shape[0]):
    if calcflag[i] > 0:
        if orig_ind.size ==0:
            orig_ind = np.array(i)
        else:
            orig_ind = np.append(orig_ind,i)
##print(orig_ind,orig_ind.shape)

filename=indir+conf_prefix+'.dat'
print(filename)
infile=open(filename,"r")
lines=infile.readlines()
infile.close()

for line in lines[:-1]:
    tokens = line.split()
##    timeindex = int(tokens[0])
    if len(tokens) > nconstraints+2:
        print(len(tokens),nconstraints)
        continue
    res = [eval(i) for i in tokens[1:nconstraints]]

    if float(tokens[0]) == 1.1:
        c0vec = np.array(res)
        print("tokens tokens",len(tokens),len(res),nconstraints)

    elif float(tokens[0]) ==start:
        vec_tseries = np.array(res)
        tseries = np.array([float(tokens[0])+1])
    elif float(tokens[0]) > start: # and int(float(tokens[0]))%1000 ==0:
        vec_tseries = np.vstack((vec_tseries,np.array(res)))
        tseries = np.append(tseries,float(tokens[0])+1)


#### we have the full time series of vectors, now analyse it
print(vec_tseries.shape)

### get a pcc time series, linf time series and overlap time series
pearson_tseries= np.zeros([vec_tseries.shape[0]],np.float64)
linf_tseries= np.zeros([vec_tseries.shape[0]],np.float64)
linf_id_tseries= np.zeros([vec_tseries.shape[0]],np.float64)
Q_tseries= np.zeros([vec_tseries.shape[0]],np.float64)
stst=0

for i in range(stst,vec_tseries.shape[0]):

    pcc_xdata = calcflag
    pcc_ydata = np.zeros_like(pcc_xdata)
    pcc_ydata[vec_tseries[i,:] > c0vec - 0.1] = 1
    pcc_ydata[calcflag==0] = 0

##    print(pcc_xdata)
##    print(pcc_ydata)
    
    cov = np.cov(pcc_xdata,pcc_ydata)
##    print(cov)
    lamb,vec = np.linalg.eig(cov)
    pearson = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
    pearson_tseries[i] = pearson


    #### now compute the overlap function
    
    q_xdata = vec_tseries[stst,calcflag > 0]
    q_ydata = vec_tseries[i,calcflag > 0]
    q_ncons = q_xdata.shape[0]
##    print(q_ncons,vec_tseries.shape)
    

    Q_tseries[i] = 1.0 - np.where(abs(q_ydata - q_xdata) > 0.1)[0].shape[0]/q_ncons
    
    
    linf_xdata = c0vec[calcflag > 0]
    linf_ydata = vec_tseries[i,calcflag > 0]


    flag = linf_ydata <= linf_xdata #+ 0.05
##    print(flag)

    xdata = linf_xdata * flag
    ydata = linf_ydata * flag
##
##    xdata = xdata[flag]
##    ydata = ydata[flag]

    linf_tseries[i] = np.max(abs(ydata - xdata))

##    xdata = vec_tseries[stst,calcflag > 0] # - c0vec[calcflag > 0]
##    ydata = vec_tseries[i,calcflag > 0] # - c0vec[calcflag > 0]

##    xdata = vec_tseries[stst,calcflag > 0] - c0vec[calcflag > 0]
##    ydata = vec_tseries[i,calcflag > 0] - c0vec[calcflag > 0]

    
##    diffs = [abs(min(ydata[j] - xdata[j],0)) for j in range(xdata.shape[0])]
##    print(diffs)
##    print(xdata)
##    print(ydata)
##    linf_tseries[i] = np.max(np.array(diffs))
##
##    linf_id_tseries[i] = orig_ind[np.argmax(abs(ydata - xdata))]

    

##    print("calcflag",np.sum(calcflag),xdata.shape,ydata.shape)

##    xdata = c0vec
##    ydata = vec_tseries[i,:]

    

##print(pearson_tseries)
##print(linf_tseries)
##print(Q_tseries)


##print(pearson_tseries)
##print(linf_tseries)
##print(Q_tseries)



plt.clf()
plt.plot(0.00035*tseries[stst:],pearson_tseries[stst:],color='red',lw=4,label="Pearson correlation",linestyle='--', marker='o')
plt.xlabel('time [s]',fontsize=20)
plt.ylabel('PCC$_{IPD}(t)$',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(0.5,1)
plt.xscale('log')
plt.tight_layout()

plotname=outdir+"/PCC.png"
plt.savefig(plotname)
##plt.show()


filename=outdir+'/PCC.dat'
print(filename,stst,vec_tseries.shape[0])
f=open(filename,"w+") # put proper output file name with temp and fblamb here
##f.write("%8d %4.4f\n" % (0,1.0))
for i in range(stst,vec_tseries.shape[0]):
##    if denomlogMSD[i]>0 and i > 0:
    f.write("%4.4f %4.8f %4.8f %4.8f\n" % (tseries[i],pearson_tseries[i],Q_tseries[i],linf_tseries[i]))
f.close()



plt.clf()
plt.plot(0.00035*tseries[stst:],linf_tseries[stst:],color='black',label="L$_{\infty}$",linestyle='--', marker='o')
plt.xlabel('time [s]',fontsize=20)
plt.ylabel('L$_{\infty}(t)$',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(0,0.3)
plt.xscale('log')
plt.tight_layout()

plotname=outdir+"/Linf.png"
plt.savefig(plotname)

plt.clf()
plt.plot(0.00035*tseries[stst:],linf_id_tseries[stst:],color='green',lw=4,label="L$_{\infty}$ ID",linestyle='--', marker='o')
plt.xlabel('time [s]',fontsize=20)
plt.ylabel('ID L$_{\infty}(t)$',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(0,300)
plt.xscale('log')
plt.tight_layout()

plotname=outdir+"/Linf_id.png"
plt.savefig(plotname)
##plt.show()

print("plotting Q(t)")
plt.clf()
plt.plot(0.00035*tseries[stst:],Q_tseries[stst:],color='blue',lw=4,label="overlap function Q(t)")
plt.plot(0.00035*tseries,pearson_tseries[stst:],color='red',lw=4,label="PCC",linestyle='--', marker='o')
plt.plot(0.00035*tseries,linf_tseries[stst:],color='black',lw=4,label="L$_{\infty}$")
plt.legend(loc="center left", prop={'size': 15})
plt.xlabel('time [s]',fontsize=20)
plt.ylabel('Q$_{IPD}(t)$',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(0,1.05)
##plt.xlim(0.000035,100)
plt.tight_layout()
plt.xscale('log')
##plotname=indir+"Decay.png"
##outdir='46chrom_batch_geom'+dirsuffix
subprocess.run(["mkdir","-p",outdir])
plotname=outdir+"/Q_t.png"
print(plotname)
plt.savefig(plotname)
##plt.show()
print("done plotting Q(t)")
    
    
##### now output histos

#### need an index mapper

maxind = (23+1)*22 - int(22*(22-1)/2) + 23 - 22   ##(p->nChrom+1)*ii - (int)(ii*(ii-1)/2) + jj - ii
print(maxind)
ind_map = np.array([],np.float64)

for i in range(24*24-1):
    indi = int( i/24)
    indj = i % 24

    if indj >= indi:
        
        if i ==0:
            ind_map = np.array([indi,indj])
        if i > 0:
            ind_map = np.vstack((ind_map,np.array([indi,indj])))

##        print(i,indi,indj)

print(maxind,ind_map.shape)

c0mat = np.zeros([24,24],np.float64)
cflagmat = np.zeros([24,24],np.float64)
cvmat = np.zeros([24,24],np.float64)
for i in range(len(ind_map)):
    c0mat[ind_map[i][0],ind_map[i][1]] = c0vec[i]
    c0mat[ind_map[i][1],ind_map[i][0]] = c0vec[i]
    cflagmat[ind_map[i][0],ind_map[i][1]] = calcflag[i]
    cflagmat[ind_map[i][1],ind_map[i][0]] = calcflag[i]
##    print(i,ind_map[i][0],ind_map[i][1],calcflag[i])

    
filename = 'ipd_init_mat.dat'
outfile = open(filename,"w")
for i in range(c0mat.shape[0]-1):
    for j in range(c0mat.shape[1]-1):
        outfile.write("%d %d %f\n" %(i,j,c0mat[i][j]))


outfile.close()
for i in range(len(ind_map)):
    cvmat[ind_map[i][0],ind_map[i][1]] = vec_tseries[vec_tseries.shape[0]-1,i]
    cvmat[ind_map[i][1],ind_map[i][0]] = vec_tseries[vec_tseries.shape[0]-1,i]
plt.clf()
filename = outdir+'/ipd_quench_mat.dat'
outfile = open(filename,"w")
for i in range(c0mat.shape[0]-1):
    for j in range(c0mat.shape[1]-1):
        
        outfile.write("%d %d %2.4f %2.4f %d\n" %(i,j,c0mat[i][j],cvmat[i][j],cflagmat[i][j]))
plt.imshow(cvmat[:23,:23])
##plt.show()
outfile.close()
##print(ind_map)

blacklist = np.unique(np.where(ind_map>=22)[0][:-2])
selflist = np.where(ind_map[:,0]==ind_map[:,1])[0]
smalllist = np.unique(np.where((ind_map >=18) & (ind_map < 23))[0])
##print("filt blacklist",blacklist)


plt.clf()
xdata = c0vec[calcflag > 0]
fin_ind = np.arange(vec_tseries[0,:].shape[0])
fin_ind_filt = fin_ind[calcflag>0]
ydata_init = vec_tseries[0,calcflag > 0]
y_bl_init = vec_tseries[0,blacklist]
x_bl = c0vec[blacklist]
##print(fin_ind,fin_ind_filt)
##print(xdata)
##print(ydata_init)
##print(abs(ydata_init - xdata))

##print(np.column_stack((xdata,ydata_init,abs(ydata_init - xdata))))
##ind_map_filt = ind_map[calcflag>0]
##for i in range(xdata.shape[0]):
####    print(np.column_stack((xdata,ydata_fin,abs(ydata_fin - xdata),
####                       ind_map[calcflag>0,0],ind_map[calcflag>0,1],
####                       fin_ind_filt)))
##    print("%4.3f %4.3f %4.3f %d %d %d" 
##          %(xdata[i],ydata_init[i],abs(ydata_init[i] - xdata[i]),
##          ind_map_filt[i,0],ind_map_filt[i,1],fin_ind_filt[i]))
    

linf_init = np.max(abs(ydata_init - xdata))
linf_id_init = orig_ind[np.argmax(abs(ydata_init - xdata))]
id_init = ind_map[linf_id_init,:]

print(linf_init,linf_id_init,id_init)

counts,bins=np.histogram(abs(ydata_init - xdata),bins=10)
plt.hist(bins[:-1],bins,weights=counts,color='red',label='all',histtype='step',linewidth=4)
counts,bins=np.histogram(abs(vec_tseries[0,blacklist] - c0vec[blacklist]),bins=10)
plt.hist(bins[:-1],bins,weights=counts,color='blue',label='chrX',histtype='step',linewidth=4)
counts,bins=np.histogram(abs(vec_tseries[0,selflist] - c0vec[selflist]),bins=10)
plt.hist(bins[:-1],bins,weights=counts,color='green',label='homo',histtype='step',linewidth=4)
counts,bins=np.histogram(abs(vec_tseries[0,smalllist] - c0vec[smalllist]),bins=10)
plt.hist(bins[:-1],bins,weights=counts,color='violet',label='small',histtype='step',linewidth=4)
plt.legend(loc="upper right", prop={'size': 15})
plt.figtext(0.20, 0.88, "$L_{\infty}$ : "+str(linf_init)+" : "+str(linf_id_init) ,size=20)
plt.figtext(0.20, 0.78, "id: "+str(id_init[0])+" , " + str(id_init[1]),size=20)
plotname=outdir+"/L_histo_init_"+sys.argv[5]+"_inp"+sys.argv[6]+"_"+sys.argv[7]+".png"
plt.savefig(plotname)



plt.clf()


ydata_fin = vec_tseries[-1,calcflag > 0]
y_bl_fin = vec_tseries[-1,blacklist]
##print(xdata)
##print(ydata_fin)
##print(abs(ydata_fin - xdata))
print(xdata.shape,ind_map.shape,fin_ind_filt.shape)
##print(ind_map[calcflag>0,:])
ind_map_filt = ind_map[calcflag>0]
##for i in range(xdata.shape[0]):
##    print(np.column_stack((xdata,ydata_fin,abs(ydata_fin - xdata),
##                       ind_map[calcflag>0,0],ind_map[calcflag>0,1],
##                       fin_ind_filt)))
##    print("%4.3f %4.3f %4.3f %4.3f %4.3f %d %d %d" 
##          %(xdata[i],ydata_init[i],abs(ydata_init[i] - xdata[i]),
##            ydata_fin[i],abs(ydata_fin[i] - xdata[i]),
##          ind_map_filt[i,0],ind_map_filt[i,1],fin_ind_filt[i]))


linf_fin = np.max(abs(ydata_fin - xdata))
linf_id_fin = orig_ind[np.argmax(abs(ydata_fin - xdata))]
id_fin = ind_map[linf_id_fin,:]



print(linf_fin,linf_id_fin,id_fin)


##print(abs(vec_tseries[-1,blacklist] - c0vec[blacklist]))
##print(blacklist)
##print(c0vec[blacklist])
##print(vec_tseries[-1,blacklist])
counts,bins=np.histogram(abs(ydata_fin - xdata),bins=10)
plt.hist(bins[:-1],bins,weights=counts,color='red',label='all',histtype='step',linewidth=4)
counts,bins=np.histogram(abs(vec_tseries[-1,blacklist] - c0vec[blacklist]),bins=10)
plt.hist(bins[:-1],bins,weights=counts,color='blue',label='chrX',histtype='step',linewidth=4)
counts,bins=np.histogram(abs(vec_tseries[-1,selflist] - c0vec[selflist]),bins=10)
plt.hist(bins[:-1],bins,weights=counts,color='green',label='homo',histtype='step',linewidth=4)
counts,bins=np.histogram(abs(vec_tseries[-1,smalllist] - c0vec[smalllist]),bins=10)
plt.hist(bins[:-1],bins,weights=counts,color='violet',label='small',histtype='step',linewidth=4)
plt.legend(loc="upper right", prop={'size': 15})
plt.figtext(0.20, 0.88, "$L_{\infty}$ : "+str(linf_fin)+" : "+str(linf_id_fin) ,size=20)
plt.figtext(0.20, 0.78, "id: "+str(id_fin[0])+" , " + str(id_fin[1]),size=20)
plotname=outdir+"/L_histo_fin_"+sys.argv[5]+"_inp"+sys.argv[6]+"_"+sys.argv[7]+".png"
plt.savefig(plotname)   
