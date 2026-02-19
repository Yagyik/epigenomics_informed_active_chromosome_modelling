import numpy as np
import matplotlib.pyplot as plt
import time
import math
import sys
import os

##aslist=['0.2','0.45','0.7','0.95']
##sgadlist=['0.008','0.01','0.015','0.02']
##colorlist=['black','red','green','blue']

##
##aslist=['0.3','0.38','0.46','0.54','0.6','0.68']
##sgadlist=['0.001','0.01']
##colorlist=['black','red','green','blue','violet','turquoise'] 

strp1=sys.argv[1]
strp2=sys.argv[2]
nset = int(sys.argv[3])
nbins = int(sys.argv[4])
insuffix=sys.argv[5]
pref=sys.argv[6]
Tp=sys.argv[7]
outdir=sys.argv[8]
logbase=4

filename=strp1+'/OutPS1_in'+insuffix+'_Tp'+Tp+'/'+strp2+'/msd_dir/allMSD.dat'
if os.path.isfile(filename):
    infile = open(filename,"r");
    lines=infile.readlines()
    infile.close()
else:
    print("file path error",filename)
    exit(0)
flength = len(lines)

avg_allmsd = np.zeros([flength],np.float64)
avg_chrmsd = np.zeros([flength],np.float64)
avg_patchmsd = np.zeros([flength],np.float64)
avg_lmnmsd = np.zeros([flength],np.float64)
avg_time = np.zeros([flength],np.float64)
denom = np.zeros([flength],np.float64)
for i in range(1,nset+1):
    ## get dir
    filename=strp1+'/OutPS'+str(i)+'_in'+insuffix+'_Tp'+Tp+'/'+strp2+'/msd_dir/allMSD.dat'
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
        if linecount < flength and float(tokens[0]) < 4000000:
            avg_time[linecount] += float(tokens[0])
            avg_allmsd[linecount] += float(tokens[1])
            avg_chrmsd[linecount] += float(tokens[2])
            avg_patchmsd[linecount] += float(tokens[3])
            avg_lmnmsd[linecount] += float(tokens[4])
            denom[linecount] += 1.0
            otime = avg_time[linecount]
            linecount+=1


avg_time /= denom
avg_allmsd /= denom
avg_chrmsd /= denom
avg_patchmsd /= denom
avg_lmnmsd /= denom
print(avg_time,"time")
print(avg_allmsd)
start = avg_time[0]
stop = avg_time[-1]
skip = avg_time[1] - avg_time[0]

log_time = np.zeros([nbins],np.float64)
log_allmsd = np.zeros([nbins],np.float64)
log_chrmsd = np.zeros([nbins],np.float64)
log_patchmsd = np.zeros([nbins],np.float64)
log_lmnmsd = np.zeros([nbins],np.float64)
denom_log_msd = np.zeros([nbins],np.float64)


maxlogt = math.log((stop-start),logbase)
minlogt = math.log(skip,logbase)

dlogt = maxlogt/nbins
print(start,stop,skip,flength,maxlogt,minlogt,dlogt)
for i in range(flength):
    
    if i>0:
        skip = avg_time[i] - avg_time[i-1]
        log10t = math.log(i*skip,logbase)
        spamindex = int(log10t/dlogt)
##        print(i,i*skip,log10t,spamindex)
    else:
        spamindex = 0

    if spamindex < nbins:
        log_time[spamindex] += start + math.pow(logbase,((i+0.5)*dlogt))
        denom_log_msd[spamindex] += 1.0
        log_allmsd[spamindex] += avg_allmsd[i]
        log_chrmsd[spamindex] += avg_chrmsd[i]
        log_patchmsd[spamindex] += avg_patchmsd[i]
        log_lmnmsd[spamindex] += avg_lmnmsd[i]

##log_time /= denom_log_pcc
##log_pcc /= denom_log_pcc

print(denom_log_msd)


##filename=sys.argv[8]+'/MSD/MSD'+sys.argv[5]+'_lmn'+sys.argv[6]+'_sgat'+sys.argv[7]+'_sgad'+sgad+'_as'+a+'.dat'
filename=outdir+'/MSD_'+strp2+'.dat'
if os.path.isfile(filename):
    os.remove(filename)
print(filename)
f=open(filename,"w+") # put proper output file name with temp and fblamb here
##f.write("%8d %4.4f\n" % (0,1.0))
for i in range(nbins):
    if denom_log_msd[i]>0 and i > 0:
        f.write("%4.4f %4.4f %4.4f %4.4f %4.4f\n" % (start + math.pow(logbase,((i+0.5)*dlogt)),log_allmsd[i]/denom_log_msd[i],
                                                   log_chrmsd[i]/denom_log_msd[i],
                                                   log_patchmsd[i]/denom_log_msd[i],
                                                   log_lmnmsd[i]/denom_log_msd[i]))
f.close()

##
##plt.plot(0.00035*avg_time,avg_chrmsd,color=colorlist[ai],label="$v_a$="+aslist[ai],linestyle='--', linewidth=3+len(aslist)-ai) #marker='o')
##
##        
##
##plt.figtext(0.20, 0.88, "Geom : "+sys.argv[5][-3:],size=20)
##plt.figtext(0.20, 0.78, "$\Delta\sigma_{a}$ = "+sgad,size=20)
####    plt.figtext(0.20, 0.77, f'sScSlope = {np.mean(xdata)/np.mean(ydata)*plotslope:.2f}',size=25)    
##plt.xlabel('time [s]',fontsize=20)
##plt.ylabel('MSD Chrom centre [$\mu$m$^2$]',fontsize=20)
##plt.xticks(fontsize=18)
##plt.yticks(fontsize=18)
####        plt.ylim(0,1)
##plt.xscale('log')
##plt.yscale('log')
##plt.legend(loc="lower right", prop={'size': 15})
##plt.tight_layout()
##plt.savefig(sys.argv[8]+'/MSD/MSDChr'+sys.argv[5]+'_lmn'+sys.argv[6]+'_sgat'+sys.argv[7]+'_sgad'+sgad+'.png')
##plt.show()


##
##plt.clf()
##plt.plot(,linf_tseries[stst:],color='black',label="L$_{\infty}$",linestyle='--', marker='o')
##plt.xlabel('time [s]',fontsize=20)
##plt.ylabel('L$_{\infty}(t)$',fontsize=20)
##plt.xticks(fontsize=18)
##plt.yticks(fontsize=18)
##plt.ylim(0,1)
##plt.xscale('log')
##plt.tight_layout()
##plt.show()


