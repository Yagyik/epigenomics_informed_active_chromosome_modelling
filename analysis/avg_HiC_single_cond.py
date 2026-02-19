import numpy as np
import matplotlib.pyplot as plt
import time
import math
import sys
import os

aslist=['0.2','0.45','0.7','0.95']
sgadlist=['0.008','0.01','0.015','0.02']

strp1=sys.argv[1]
strp2=sys.argv[2]
nset = int(sys.argv[3])
sepcut = sys.argv[4]
decay = sys.argv[5]
insuffix=sys.argv[6]
pref=sys.argv[7]
Tp=sys.argv[8]
outdir=sys.argv[9]


filename=strp1+'/OutPS1_in'+insuffix+'_Tp'+Tp+'/'+strp2+'/HiC_dir/HiCMat_cl_'+str(sepcut)+'_'+str(decay)+'.dat'
infile = open(filename,"r");
lines=infile.readlines()
infile.close()
flength = len(lines)
maxx=int(lines[flength-1].split()[0])
maxy=int(lines[flength-1].split()[1])
print(maxx,maxy)

avg_HiC = np.zeros([maxx,maxy],np.float64)
plt.clf()
fig,axs=plt.subplots(nrows=1, ncols=1,sharex=True,sharey=True,figsize=(15,10),gridspec_kw={'width_ratios': [0.8]})
for i in range(1,nset+1):
    ## get dir
    filename=strp1+'/OutPS'+str(i)+'_in'+insuffix+'_Tp'+Tp+'/'+strp2+'/HiC_dir/HiCMat_cl_'+str(sepcut)+'_'+str(decay)+'.dat'
    print("now opening",filename)
    if os.path.isfile(filename):
        infile = open(filename,"r");
        lines=infile.readlines()
        infile.close()
    else:
        continue
                 
                
    linecount=0
                 
    otime=-1

    for line in lines:
        tokens = line.split()

        chri = int(tokens[0])-1
        chrj = int(tokens[1])-1
        avg_HiC[chri,chrj] += float(tokens[2])

avg_HiC /= nset

cmax = np.max(avg_HiC)
cmin = np.min(avg_HiC)
print(cmin,cmax,100)
xaxis = list(range(0,maxx,1))
yaxis = list(range(0,maxy,1))
##    cbticklabels=list(np.arange(1000,50000,5000))
cbticklabels = list(np.linspace(cmin,1.2*cmax,10))
##    print(cmin,cmax)
im1 = axs.imshow(avg_HiC[:,:],vmin=cmin,vmax=1.2*cmax,aspect='auto', cmap='Blues')
cb = fig.colorbar(im1)
cb.ax.set_yticklabels([f'{x:.2f}' for x in cbticklabels],size=20)
cb.set_label("Contact frequency", labelpad=-1, size=25)
##    print(xaxis,yaxis)
axs.set_xticks(xaxis)
axs.set_yticks(yaxis)
axs.set_xticklabels(xaxis, fontsize=20)
axs.set_yticklabels(yaxis, fontsize=20)



plotname = outdir + '/HiCMat_'+strp2+'.png'
##    plt.show()
plt.savefig(plotname)


filename = outdir + '/HiCMat_'+strp2+'.dat'

print(filename)
f=open(filename,"w+") # put proper output file name with temp and fblamb here

for i in range(avg_HiC.shape[0]):
    for j in range(avg_HiC.shape[1]):
        f.write("%d %d %4.4f\n" % (i+1,j+1,avg_HiC[i,j]))
f.close()

        



            
