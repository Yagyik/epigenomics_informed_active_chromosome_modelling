import numpy as np
import matplotlib.pyplot as plt
import time
import math
import sys
import os

##
##aslist=['0.3','0.38','0.46','0.54','0.6','0.68']
##sgadlist=['0.001','0.01']
##colorlist=['black','red','green','blue','violet','turquoise']

##filestr=sys.argv[1]
strp1=sys.argv[1]
strp2=sys.argv[2]
nset = int(sys.argv[3])
nbins = int(sys.argv[4])
insuffix=sys.argv[5]
pref=sys.argv[6]
Tp=sys.argv[7]
outdir=sys.argv[8]
logbase=4


filename=strp1+'/OutPS8_in'+insuffix+'_Tp'+Tp+'/'+strp2+'/PCC.dat'
if os.path.isfile(filename):
    infile = open(filename,"r");
    lines=infile.readlines()
    infile.close()
else:
    print("file path error",filename)
    exit(0)
flength = len(lines)

avg_pcc = np.zeros([flength],np.float64)
avg_qt = np.zeros([flength],np.float64)
avg_linf = np.zeros([flength],np.float64)
avg_time = np.zeros([flength],np.float64)
denom = np.zeros([flength],np.float64)
for i in range(1,nset+1):
    ## get dir  
    filename=strp1+'/OutPS'+str(i)+'_in'+insuffix+'_Tp'+Tp+'/'+strp2+'/PCC.dat'
    print(i,filename)
    infile = open(filename,"r");
    lines=infile.readlines()
    infile.close()
    linecount=0
    otime=-1


    for line in lines:
        tokens = line.split()
        tt = float(tokens[0])
        if linecount < flength and int(tt)%500 and tt < 4000000 and otime < tt:
            avg_time[linecount] += tt #/0.00035
            avg_pcc[linecount] += float(tokens[1])
            avg_qt[linecount] += float(tokens[2])
            avg_linf[linecount] += float(tokens[3])
            denom[linecount] += 1.0
            otime = tt
            linecount+=1


avg_time /= denom
avg_pcc /= denom
avg_qt /= denom
avg_linf /= denom
##print(avg_time)
##print(avg_pcc)
##print(denom)
filename=outdir+'/PCC'+strp2+'.dat'
print(filename)
if os.path.isfile(filename):
    os.remove(filename)
f=open(filename,"w+") # put proper output file name with temp and fblamb here

##f.write("%8d %4.4f\n" % (0,1.0))
for i in range(len(denom)):
    if denom[i]>0 and i > 0:
        f.write("%4.8f %4.8f %4.8f %4.8f\n" % (avg_time[i]*0.001*0.035,avg_pcc[i],
                                             avg_qt[i],avg_linf[i]))
f.close()



start = avg_time[0] #/(0.035*0.001)
stop = avg_time[-1] #/(0.035*0.001)
skip = (avg_time[1] - avg_time[0]) #/(0.035*0.001)

log_time = np.zeros([nbins],np.float64)
log_pcc = np.zeros([nbins],np.float64)
log_qt = np.zeros([nbins],np.float64)
log_linf = np.zeros([nbins],np.float64)
denom_log_pcc = np.zeros([nbins],np.float64)


maxlogt = math.log((stop-start),logbase)
minlogt = math.log(skip,logbase)

dlogt = (maxlogt - minlogt)/nbins
print(start,stop,skip,flength,maxlogt,minlogt,dlogt)
for i in range(flength):
    
    if i>0:
        skip = (avg_time[i] - avg_time[i-1]) #/(0.035*0.001)
        if skip <= 0.0:
            continue
##        print(skip,logbase)
        log10t = math.log(i*skip,logbase)
        spamindex = int(log10t/dlogt)
##        print(i,start+i*skip,log10t,spamindex)
    else:
        spamindex = 0

    if spamindex < nbins:
        log_time[spamindex] += start + math.pow(logbase,((i+0.5)*dlogt))
        denom_log_pcc[spamindex] += 1.0
        log_pcc[spamindex] += avg_pcc[i]
        log_qt[spamindex] += avg_qt[i]
        log_linf[spamindex] += avg_linf[i]

##log_time /= denom_log_pcc
##log_pcc /= denom_log_pcc

print(denom_log_pcc)


filename=outdir+'/PCC_logBin_'+strp2+'.dat'
if os.path.isfile(filename):
    os.remove(filename)
print(filename)
f=open(filename,"w+") # put proper output file name with temp and fblamb here
##f.write("%8d %4.4f\n" % (0,1.0))
for i in range(nbins):
    if denom_log_pcc[i]>0 and i > 0:
        f.write("%4.8f %4.8f %4.8f %4.8f\n" % ((start + math.pow(logbase,((i+0.5)*dlogt)))*(0.035*0.001),log_pcc[i]/denom_log_pcc[i],
                                             log_qt[i]/denom_log_pcc[i],log_linf[i]/denom_log_pcc[i]))
f.close()



        ##
        





