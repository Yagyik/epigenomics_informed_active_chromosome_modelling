import numpy
import matplotlib.pyplot as plt
import time
import math
import sys
from scipy.optimize import curve_fit


## program to read a thermo file and log bin it,
## provided a start, stop, skip, Tp, taup, f, nbins


mainpath=sys.argv[1]
indir=sys.argv[2]
pref=sys.argv[3]
insuffix=sys.argv[4]
Tp=sys.argv[5]
nbins=int(sys.argv[6])
logbase=float(sys.argv[7])
outdir=sys.argv[8]
skip=int(sys.argv[9])
fact=float(sys.argv[10])

print("mainpath,indir,pref,insuffix,Tp,nbins,logbase,outdir,skip,fact",mainpath,indir,pref,insuffix,Tp,nbins,logbase,outdir,skip,fact)
nset=int(sys.argv[11])
##avgdir=sys.argv[11]
### first we read the original thermo, keep 

##flength = int ((stop-start)/skip)
##print(flength,nbins)
##time=numpy.zeros([flength],numpy.float64)
##thermo=numpy.zeros([flength],numpy.float64)
##avgtime=numpy.zeros([flength],numpy.float64)
##avgthermo=numpy.zeros([flength],numpy.float64)

def stretch_exp(x,U0,deltaU,tau,beta):
    return U0 + deltaU*(numpy.exp(-(x/tau)**beta))



avgtime=numpy.array([])
avgthermo=numpy.array([])
avgdenom=numpy.array([])
maxflength = 0
timelist=[]
thermolist=[]
thermosqlist=[]
for ss in range(1,nset):
##    if fact == 0.0:
##        filename=sys.argv[1]+str(ss)+"/Tpar"+sys.argv[3]+"/taup0.01/Conf_N1000_rho1.2_T0.001_taup0.01_fact"+sys.argv[5]+"-thermo.dat"
##    else:
##        filename=sys.argv[1]+str(ss)+"/Tpar"+sys.argv[3]+"/taup"+sys.argv[4]+"/Conf_N1000_rho1.2_T0.001_taup"+sys.argv[4]+"_fact"+sys.argv[5]+"-thermo.dat"
##    Set_5_insph3_sph_Tp0.2
    filename=mainpath+"Set_"+str(ss)+"_in"+insuffix+"_"+pref+"_Tp"+Tp+"/"+indir+"/Conf-thermo.dat"
    print(filename)
    infile=open(filename,"r")
    lines=infile.readlines()
    infile.close()
    
##    time=numpy.zeros([flength],numpy.float64)
##    thermo=numpy.zeros([flength],numpy.float64)
    time=numpy.array([])
    thermo=numpy.array([])
    denom=numpy.array([])

    ## now store the file sequentially, trying to ensure that we don't backtrack over our files
    otime=-1
    linecount=0
    for line in lines:
##        if linecount >= flength:
##            break
        tokens=line.split()
##        print(line)
        ## do the store only if we are in multiples of skip
        if int(tokens[0]) % skip == 0: # and int(tokens[0]) >= start and int(tokens[0]) <= stop:
            ## check if the current time is greater than the old time
            if int(tokens[0]) >= otime:
                time = numpy.append(time,float(tokens[0]))
                thermo = numpy.append(thermo,float(tokens[1]))
                denom = numpy.append(denom,1)
                
##                time[linecount] = float(tokens[0])
##                thermo[linecount] = float(tokens[1])
                otime = float(tokens[0]) # time[linecount] ## otime remains the last read entry to the array
                linecount += 1

        if linecount %10000 == 0:
            print(line)

    ### curve-fit for this data set -- no good, too noisy
##    lb=thermo[0]-50
##    ub=thermo[0]+50
##    popt,pocv = curve_fit(stretch_exp,time,thermo,p0=[thermo[0],20,10000,0.5],bounds=([lb,-4000,1,0.25],[ub,6000,numpy.inf,0.75]))
##    U0 = popt[0]
##    deltaU = popt[1]
##    tau = popt[2]
##    beta = popt[3]
##    print(ss,fact,tau,beta,U0,deltaU)

    ## update the average thermo array
    flength = time.size
    if flength > maxflength:
        maxflength = flength

    ## get array of arrays
    if ss==1:
        setsize=2000
        collectthermo=thermo[-setsize:].reshape((1,setsize))
        print(collectthermo.shape)
        
    else:
        collectthermo=numpy.vstack((collectthermo,thermo[-setsize:].reshape((1,setsize))))
        print(collectthermo.shape)
            

    avglength = avgtime.size
    timelist.append(time)
    thermolist.append(thermo)
    thermosqlist.append((thermo)*(thermo))
    
    for i in range(flength):
        if i >= avglength:
            avgtime = numpy.append(avgtime,time[i])
            avgthermo = numpy.append(avgthermo,thermo[i])
            avgdenom = numpy.append(avgdenom,denom[i])
        else:
            avgtime[i] += time[i]
            avgthermo[i] += thermo[i]
            avgdenom[i] += denom[i]
        if i%int(0.05*flength) ==0:
            print(i,time[i],thermo[i],flength)
    
    ## sanity check
    print(maxflength,avgtime.size,"should be equal")

print([x.shape[0] for x in thermolist])
nruns=collectthermo.shape[0]
minlength = min([x.shape[0] for x in thermolist])
spamarr = numpy.zeros([nruns,minlength],numpy.float64)
spamsqarr = numpy.zeros([nruns,minlength],numpy.float64)
for i in range(len(thermolist)):
    spamsqarr[i,:] = thermosqlist[i][:minlength]
    spamarr[i,:] = thermolist[i][:minlength]


spamsqmean = numpy.mean(spamsqarr,axis=0)
spammean = numpy.mean(spamarr,axis=0)

print("sq mean",spamsqmean)
print("mean",spammean)
print("varr",spamsqmean - spammean*spammean)

thermostdarr=numpy.sqrt(spamsqmean - spammean*spammean)

print("std",thermostdarr)
         
#### check collectthermo size
print("\n\n")
print("shape of all things",collectthermo.shape)
avgrange=setsize
print("mean",numpy.mean(collectthermo[:,-avgrange:])) ### mean both axes or ,axis=0)) ## mean over runs
print("std",numpy.std(collectthermo[:,-avgrange:])) ## mean both axes or ,axis=0))
spammean=numpy.mean(collectthermo[:,-avgrange:])
spamstdev=numpy.std(collectthermo[:,-avgrange:])
filename=outdir+"/avgStDThermo_in"+insuffix+"_pref"+pref+"_v_f.dat"
if fact==0: ## new file
    f=open(filename,"w+")
else:
    f=open(filename,"a")
f.write("%4.2f %4.12f %4.12f\n" % (fact,spammean,spamstdev))
f.close()
print("\n\n")

for i in range(avgtime.size):
    avgtime[i] /= avgdenom[i]
    avgthermo[i] /= avgdenom[i]


##print(avgtime,avgtime.size)
##print(avgthermo,avgthermo.size)
print(avgthermo[0],avgtime[:-1])
th0=avgthermo[0]
lb=avgthermo[0]-20
ub=avgthermo[0]+20
maxtime=avgtime[avgtime.size-1]
## now fit avg thermo to stretched exponential form
##popt,pocv = curve_fit(stretch_exp,avgtime,avgthermo,p0=[avgthermo[0],20,10000,0.5],bounds=([lb,-4000,1,0.45],[ub,6000,numpy.inf,0.55]))





logthermo = numpy.zeros([nbins],numpy.float64)
denomlogthermo = numpy.zeros([nbins],numpy.float64)
logtime=numpy.zeros([nbins],numpy.float64)

collectlogthermo = numpy.zeros([nruns,nbins],numpy.float64)
collectdenomlogthermo = numpy.zeros([nruns,nbins],numpy.float64)
start = avgtime[0]
stop = avgtime[-1]

maxlogt = math.log((stop-start),logbase)

minlogt = math.log(skip,logbase)

##dlogt = (maxlogt - minlogt)/nbins

dlogt = maxlogt/nbins
##print(minlogt,maxlogt,dlogt)
for i in range(flength):
    
    if i>0:
        log10t = math.log(i*skip,logbase)
        spamindex = int(log10t/dlogt)
        
##        print(i,i*skip,maxlogt,minlogt,dlogt,log10t,spamindex)
    else:
        spamindex = 0

    if spamindex < nbins:
        denomlogthermo[spamindex] += 1.0
        logthermo[spamindex] += avgthermo[i]



for j in range(nruns):
    print(j,maxflength,flength)
    print("list shape",j,thermolist[j].shape)
    for i in range(maxflength):
        if i>0:
            log10t = math.log(i*skip,logbase)
            spamindex = int(log10t/dlogt)
        
##        print(i,i*skip,maxlogt,minlogt,dlogt,log10t,spamindex)
        else:
            spamindex = 0
        
        if i < thermolist[j].shape[0] and spamindex < nbins:
            collectlogthermo[j,spamindex] += thermolist[j][i]
            collectdenomlogthermo[j,spamindex] += 1.0

collectlogthermostd = numpy.zeros([nbins],numpy.float64)
collectdenomlogthermostd = numpy.zeros([nbins],numpy.float64)

for i in range(maxflength):
    if i>0:
        log10t = math.log(i*skip,logbase)
        spamindex = int(log10t/dlogt)
    
##        print(i,i*skip,maxlogt,minlogt,dlogt,log10t,spamindex)
    else:
        spamindex = 0
    
    if i < thermostdarr.shape[0] and spamindex < nbins:
        collectlogthermostd[spamindex] += thermostdarr[i]
        collectdenomlogthermostd[spamindex] += 1.0
print("\n\n")
print("log shapes",collectlogthermo.shape)
print("log mean",numpy.mean(collectlogthermo/collectdenomlogthermo,axis=0,where=collectdenomlogthermo>0))
print("log meaner",numpy.std(collectlogthermo/collectdenomlogthermo,axis=0,where=collectdenomlogthermo>0))
print("log stder",numpy.divide(collectlogthermostd,collectdenomlogthermostd))
print("log denom",collectdenomlogthermo[-1,:])

##spam = numpy.mean(collectlogthermosq,axis=0,where=collectdenomlogthermo>0) - numpy.mean(collectlogthermo,axis=0,where=collectdenomlogthermo>0)*numpy.mean(collectlogthermo,axis=0,where=collectdenomlogthermo>0)
##spamspam = numpy.sqrt(numpy.divide(spam,collectdenomlogthermo,where=collectdenomlogthermo>0),where=collectdenomlogthermo>0)
####spamspamspam = numpy.mean(spamspam,axis=0,where=collectdenomlogthermo>0)
##print("all", spamspam)
##print("log bang average",numpy.average(collectlogthermo/collectdenomlogthermo,weights=(collectdenomlogthermo > 0),axis=0))
print("\n\n")
meanforf=numpy.mean(collectlogthermo/collectdenomlogthermo,axis=0,where=collectdenomlogthermo>0)
##stdforf=numpy.std(collectlogthermo/collectdenomlogthermo,axis=0,where=collectdenomlogthermo>0)
stdforf=numpy.divide(collectlogthermostd,collectdenomlogthermostd)
filename=outdir+"/logBinAvgStDThermoNS_v_time_in"+insuffix+"_pref"+pref+"_Tp"+Tp+".dat"
print(filename)
f=open(filename,"w+") # put proper output file name with temp and fblamb here
##f.write("%8d %4.4f\n" % (0,1.0))
for i in range(nbins):
    if math.isnan(meanforf[i]) == False and math.isnan(stdforf[i])==False:
        f.write("%8d %4.4f %4.4f\n" % (start + math.pow(logbase,((i+0.5)*dlogt)),meanforf[i],stdforf[i]))
f.close()
##print(denomlogthermo)

### curve fit the log-binned data




##filename=outdir+"/logBinAvgThermoNS_v_time_Tp"+sys.argv[3]+"_taup"+sys.argv[4]+"_f"+sys.argv[5]+".dat"
##print(filename)
##f=open(filename,"w+") # put proper output file name with temp and fblamb here
####f.write("%8d %4.4f\n" % (0,1.0))
##for i in range(nbins):
##    if denomlogthermo[i]>0:
##        f.write("%8d %4.4f\n" % (start + math.pow(logbase,((i+0.5)*dlogt)),logthermo[i]/denomlogthermo[i]))
##f.close()

##fitcurve = stretch_exp(avgtime,popt[0],popt[1],popt[2],popt[3])
filename=outdir+"/avgThermoNS_v_time_in"+insuffix+"_pref"+pref+"_Tp"+Tp+".dat"
print(filename)
f=open(filename,"w+") # put proper output file name with temp and fblamb here
##f.write("%8d %4.4f\n" % (0,1.0))
for i in range(flength):
    f.write("%8d %4.4f \n" % (avgtime[i],avgthermo[i]))
f.close()


