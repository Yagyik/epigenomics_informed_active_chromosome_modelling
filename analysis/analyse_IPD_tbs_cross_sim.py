import numpy as np
import scipy
import matplotlib.pyplot as plt
import sys
import subprocess


indir=sys.argv[1]
strpt2=sys.argv[2]
conf_prefix=sys.argv[3]
reffilename=sys.argv[4]
outdir=sys.argv[5]

print("new version")
print(indir)
print(conf_prefix)
print(reffilename)
print(sys.argv)
##fi=int(sys.argv[2])
nsurf=3500
nchrom=46
npatch=23
nruns = 8
nconstraints=1 + int((0.5*(23+1)*(23)) + 23)



start=0
stop=2000000
step=20000
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

orig_ind = np.array([])
for i in range(calcflag.shape[0]):
    if calcflag[i] > 0:
        if orig_ind.size ==0:
            orig_ind = np.array(i)
        else:
            orig_ind = np.append(orig_ind,i)

##### done with reference file -- reading main now


vec_tseries = np.array([],np.float64)
fullset = np.array([],np.float64)

countlines = 0     
for ss in range(1,nruns+1):

    filename=indir+"Set_"+str(ss)+strpt2+"/"+conf_prefix+".dat"
    print(filename)
    infile=open(filename,"r")
    lines=infile.readlines()
    infile.close()
    t = start

    for line in lines[:-1]:
        tokens = line.split()

        
        if len(tokens) > nconstraints+2:
            print(len(tokens),nconstraints)
            continue
        
        res = [eval(i) for i in tokens[1:nconstraints]]

        if float(tokens[0]) == 1.1:
            c0vec = np.array(res)
            print("tokens tokens",len(tokens),len(res),nconstraints)

        timeindex = int(tokens[0])
##        print("time",timeindex)
        if timeindex >= stop:
            break
        if timeindex % step != 0:
            continue

        elif float(tokens[0]) == start and ss == 1:
            vec_tseries = np.array(res)
            tseries = np.array([float(tokens[0])+1])
            fullset=np.array([countlines,ss,timeindex])
            countlines +=1

            
        elif float(tokens[0]) > start: # and int(float(tokens[0]))%1000 ==0:

            if timeindex % step == 0:
                print(ss,countlines,timeindex)
                vec_tseries = np.vstack((vec_tseries,np.array(res)))
                tseries = np.append(tseries,float(tokens[0])+1)
                fullset=np.vstack((fullset,np.array([countlines,ss,timeindex])))
                countlines +=1

        
    
final_select = np.arange(fullset.shape[0])
    

print(final_select) 
##select = np.arange(fullset.shape[0])
##print(select)
selection = fullset[final_select]
print(selection)
print(vec_tseries.shape,selection.shape)

length = selection.shape[0]
self_sim=np.zeros([length,length],np.float64)
pcc_xdata = np.copy(calcflag)
for ii in range(length):
    i = selection[ii,0]
    veci = vec_tseries[i]
##    normveci = np.linalg.norm(veci)
    
    pcc_y1data = np.zeros_like(pcc_xdata)
    pcc_y1data[veci > c0vec - 0.1] = 1
    pcc_y1data[calcflag==0] = 0
    
    print("going through rows",ii)
    for jj in range(length):
        j = selection[jj,0]
        vecj = vec_tseries[j]
##        normvecj = np.linalg.norm(vecj)

##        cov = np.cov(veci,vecj)
##        pearson = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
##        self_sim[ii,jj] = pearson

        pcc_y2data = np.zeros_like(pcc_xdata)
        pcc_y2data[vecj > c0vec - 0.1] = 1
        pcc_y2data[calcflag==0] = 0
        
        cov = np.cov(pcc_y1data,pcc_y2data)
    ##    print(cov)
        lamb,vec = np.linalg.eig(cov)
        pearson = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
        self_sim[ii,jj] = pearson


        
        



print("done calculating correlation matrix\n")
fig,axs=plt.subplots(nrows=1, ncols=1,sharex=True,sharey=True,figsize=(15,10),gridspec_kw={'width_ratios': [0.8]})
##
im1 = axs.imshow(self_sim,vmin=0.0,vmax=1.0,aspect='auto', cmap='viridis')
axs.set_title("Batched PCC$_{IPD}(t)$",fontsize=40)
        
cb = fig.colorbar(im1)
cb.ax.tick_params(labelsize=25)
xaxis = list(range(0,length,int(0.1*length)))
yaxis = list(range(0,length,int(0.1*length)))
axs.set_xticks(xaxis)
axs.set_yticks(yaxis)
axs.set_xticklabels(xaxis, fontsize=20)
axs.set_yticklabels(yaxis, fontsize=20)
##axs.set_xlabel("Sequence id",fontsize=20)
plt.tight_layout()
plotname = outdir + 'batch_IPD_constraint_self.png'
plt.savefig(plotname)        
