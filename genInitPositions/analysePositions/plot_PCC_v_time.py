import numpy as np
import matplotlib.pyplot as plt
import sys
##import matplotlib.animation as animation



indir = sys.argv[1]+'/'
##checkpoint=int(sys.argv[2]) ## config to check
expectpoint=float(sys.argv[2]) ## expected values should be 1.1
figtext=sys.argv[3]
##tskip=int(sys.argv[3])
print(expectpoint)
nchrom=23
nconstraints = int(0.5*(nchrom+1)*(nchrom)) + nchrom

nruns = 16


for i in range(1,nruns+1):
    pearsonarr = np.array([],np.float64)
    tseriesarr = np.array([],np.float64)

    matfile = indir + str(i) + "/testWhole_46c-constraints.dat"
    f = open(matfile,"r")

    lines = f.readlines()
    f.close()
    ydatalist=[]
    print(len(lines))
##    print(lines[-805],lines[-1])

    for line in lines[:10]:
        tokens=line.split()
        print(float(tokens[0]),expectpoint,len(tokens))
        if float(tokens[0])==expectpoint:
            print("expected constraint values")
            expectvec = np.array([float(x) for x in tokens[1:-2]])
        elif float(tokens[0])==1.2:
            calcflag = np.array([float(x) for x in tokens[1:-2]])

        
    

    for line in lines[:-2]:
##        print(line)
        tokens=line.split()
##        print(len(tokens),np.sum(calcflag))
        if float(tokens[0])==0:
            firstvec = np.array([float(x) for x in tokens[1:-2]])
            checkvec = np.copy(firstvec)
            ydatalist.append(firstvec)
        elif float(tokens[0])==expectpoint:
            print("expected constraint values")
            expectvec = np.array([float(x) for x in tokens[1:-2]])
##            linfindex = np.argmax(abs(firstvec-expectvec))
##            specialpoint=[expectvec[linfindex],firstvec[linfindex]]
##            specialpointlist.append(specialpoint)
    ##        print(expectvec)
        elif float(tokens[0])==1.2:
            calcflag = np.array([float(x) for x in tokens[1:-2]])
        if float(tokens[0]) != expectpoint and int(float(tokens[0]))%100==0: ## last digit of time index != 1
            checkvec = np.array([float(x) for x in tokens[1:-2]])


            
##            pxdata=expectvec[calcflag > 0]
##            pydata=checkvec[calcflag > 0]
##
##            
##
##
##            flag = pydata <= pxdata
####            print(flag)
##
##            xdata = pxdata * flag
##            ydata = pydata * flag
##
##            xdata = xdata[flag]
##            ydata = ydata[flag]

            pcc_xdata = calcflag
            pcc_ydata = np.zeros_like(pcc_xdata)
            pcc_ydata[checkvec[:] > expectvec[:] - 0.1] = 1
            pcc_ydata[calcflag==0] = 0

            xdata = np.copy(pcc_xdata)
            ydata = np.copy(pcc_ydata)

            

##            print(xdata)
##            print(ydata)
            cov = np.cov(xdata,ydata)
            pearson = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
            pearsonarr = np.append(pearsonarr,pearson)
            tseriesarr = np.append(tseriesarr,int(tokens[0])+1)

    
    plt.plot(tseriesarr,pearsonarr,lw=2,label="r"+str(i))
    
    

plt.xlabel("Iterations",fontsize=20)
plt.ylabel("PCC",fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim(1000,5000000)
plt.ylim(0.1,0.99) #,0.2)
plt.xscale('log')
##plt.legend(loc="lower right", prop={'size': 20})
plt.tight_layout()
##plt.savefig("/home/goswam_y/Dropbox/chromosome_modelling/PCC_v_time.png")
plt.text(0.2, 0.8, figtext.upper(), transform=plt.gca().transAxes, fontsize=20, va='top', ha='left', color='red')
plt.savefig("PCC_v_time_"+figtext+".png")

# plt.show()
        
        

    
    
