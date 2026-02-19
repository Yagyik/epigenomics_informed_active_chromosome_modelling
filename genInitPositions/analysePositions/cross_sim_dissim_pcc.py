import numpy as np
import matplotlib.pyplot as plt
import sys


##### script to demonstrate the dissimilarity of configurations based on raw positions
##### as a bonus would be good to show the self-similiarity of IPD entries for same configs

dirpath1 = sys.argv[1]
dirpath2 = sys.argv[2]
constraintfile = "/testWhole_46c-constraints.dat"
configprefix = "/testWhole_46c-dump_"


#### we iterate over the runs and dump number -- gathering configurations 3N vec
#### store the index and timepoint -- we find the constraint matrix for the corresponding
#### we will calculate the cosine similarity matrix for constraints + expected
#### we will also calculate the cosine similarity matrix for positions (3N vec)

nruns=1 #16
config_tseries1=np.array([],np.float64)
config_tseries2=np.array([],np.float64)
countlines=0
pselect = 1
fullset1=np.array([],np.float64)
fullset2=np.array([],np.float64)
targetselection=50
for runindex in range(1,nruns+1):

##    for timeindex in range(10000000,19990000,10000): ## 99 files
    for timeindex in range(100000,499000,1000): ## 99 files
        filename1 = dirpath1+"/"+str(runindex)+configprefix+str(timeindex)+".dat"
        filename2 = dirpath2+"/"+str(runindex)+configprefix+str(timeindex)+".dat"
        print(filename1,filename2)
        configfile1=open(filename1,"r")
        lines1=configfile1.readlines()
        configfile1.close()

        configfile2=open(filename2,"r")
        lines2=configfile2.readlines()
        configfile2.close()
##        print(len(lines),lines[2002])
        res1 = []
        res2 = []
        

        for line in lines1[2002:]: ### chromosome positions start here
##            print(line)
            tokens=line.split()

            ### flat append xyz to get 3N vec in order
            x = float(tokens[4])
            res1.append(x)
            y = float(tokens[5])
            res1.append(y)
            z = float(tokens[6])
            res1.append(z)
        ### got full config
##        print(len(res),res)
            

        if countlines==0:
            config_tseries1=np.array(res1)
            ### store the runindex and the timeindex
##            selection.append([runindex,timeindex])
            fullset1=np.array([countlines,runindex,timeindex])
            ### update the probability of next selection
##            pselect = np.random.uniform(0,(targetselection - len(selection))/targetselection)
            print(fullset1,pselect)
        elif countlines > 0:
            #### insert code here to make a selection from the full list
            #### need to store the runindex and the corresponding timeindex
            config_tseries1=np.vstack((config_tseries1,np.array(res1)))
            fullset1=np.vstack((fullset1,np.array([countlines,runindex,timeindex])))


        for line in lines2[2002:]: ### chromosome positions start here
##            print(line)
            tokens=line.split()

            ### flat append xyz to get 3N vec in order
            x = float(tokens[4])
            res2.append(x)
            y = float(tokens[5])
            res2.append(y)
            z = float(tokens[6])
            res2.append(z)
        ### got full config
##        print(len(res),res)
            

        if countlines==0:
            config_tseries2=np.array(res2)
            ### store the runindex and the timeindex
##            selection.append([runindex,timeindex])
            fullset2=np.array([countlines,runindex,timeindex])
            ### update the probability of next selection
##            pselect = np.random.uniform(0,(targetselection - len(selection))/targetselection)
            print(fullset2,pselect)
        elif countlines > 0:
            #### insert code here to make a selection from the full list
            #### need to store the runindex and the corresponding timeindex
            config_tseries2=np.vstack((config_tseries2,np.array(res2)))
            fullset2=np.vstack((fullset2,np.array([countlines,runindex,timeindex])))


        countlines+=1

##print(config_tseries.shape,pselect,len(selection))

#### generate a selection
print(fullset1.shape[0],pselect)
##select = np.random.choice(np.arange(fullset.shape[0]),size=targetselection,replace=False)
select1 = np.arange(fullset1.shape[0])
print(select1)
selection1 = fullset1[select1]
print(selection1)

##select2 = np.arange(fullset2.shape[0])
##print(select2)
selection2 = fullset2[select1]
print(selection2)


##length=config_tseries.shape[0]
length = selection1.shape[0]
self_sim=np.zeros([length,length],np.float64)
for ii in range(length):
    i = selection1[ii,0]
    veci = config_tseries1[i]
    normveci = np.linalg.norm(veci)

    
    for jj in range(length):
        j = selection2[jj,0]
        vecj = config_tseries2[j]
        normvecj = np.linalg.norm(vecj)
        cov = np.cov(veci,vecj)
        pearson = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
        self_sim[ii,jj] = pearson
##        self_sim[ii,jj] = np.dot(veci,vecj)/(normveci*normvecj)


## now show the self_sim matrix
fig,axs=plt.subplots(nrows=1, ncols=1,sharex=True,sharey=True,figsize=(15,10),gridspec_kw={'width_ratios': [0.8]})
##
im1 = axs.imshow(self_sim,vmin=-1.0,vmax=1.0,aspect='auto', cmap='viridis')
axs.set_title("Raw position cross similarity (PCC)",fontsize=40)

##cbticklabels= ['{:.1f}'.format(x) for x in np.arange(-1,1.1,0.2)]
##print(cbticklabels)       
##cb = fig.colorbar(im1)
##cb.ax.set_yticklabels(cbticklabels,size=20)
##cb.set_label("PCC", labelpad=-1, size=25)
##print(xaxis,yaxis)
cb = fig.colorbar(im1)
cb.ax.tick_params(labelsize=25)
xaxis = list(range(0,length,int(0.1*length)))
yaxis = list(range(0,length,int(0.1*length)))
axs.set_xticks(xaxis)
axs.set_yticks(yaxis)
axs.set_xticklabels(xaxis, fontsize=20)
axs.set_yticklabels(yaxis, fontsize=20)
plt.tight_layout()
#plt.show()
plt.savefig('/home/goswam_y/Dropbox/chromosome_modelling/AnaCode_chrom_model_OOP/GD_GD_pos_cross.png')



#### from the selection of (runindex,timeindex) we store the corresponding constraint matrices
##print(selection)
vec_tseries1 = np.array([],np.float64)
vec_tseries2 = np.array([],np.float64)
nselected=0
nconstraints=int((0.5*(46+1)*(46)) + 46)
newfullset1=np.array([],np.float64)
newfullset2=np.array([],np.float64)


countlines=0
for runindex in range(1,nruns+1):
    filename1 = dirpath1+"/"+str(runindex)+constraintfile
    print(filename1,countlines)
    cfile1=open(filename1,"r")
    lines1=cfile1.readlines()
    cfile1.close()


    filename2 = dirpath2+"/"+str(runindex)+constraintfile
    print(filename2,countlines)
    cfile2=open(filename2,"r")
    lines2=cfile2.readlines()
    cfile2.close()
##        print(len(lines),lines[2002])

    minlines = min(len(lines1),len(lines2))
    print(len(lines1),len(lines2),minlines)
    for i in range(minlines):
        line1 = lines1[i]   
        tokens1=line1.split()
        timeindex1 = int(tokens1[0])
        res1 = [eval(j) for j in tokens1[1:nconstraints]]

        line2 = lines2[i]   
        tokens2=line2.split()
        timeindex2 = int(tokens2[0])
        res2 = [eval(j) for j in tokens2[1:nconstraints]]


##        print(i,len(res1),len(res2),tokens1[0],tokens2[0].timeindex1,timeindex2)

        ### get this out of the way
        if int(tokens1[0])==1:
            expectvec=np.array(res1)


        ### store all like earlier but keeping in mind the time interval for the configs
        if int(tokens1[0])>=100000 and int(tokens1[0])<499000 and int(tokens1[0])%1000==0:
            print(tokens1[0],tokens2[0])
            if countlines==0:
                vec_tseries1=np.array(res1)
                newfullset1 = np.array([countlines,runindex,int(tokens1[0])])

                vec_tseries2=np.array(res2)
                newfullset2 = np.array([countlines,runindex,int(tokens2[0])])
            elif countlines>0:
                vec_tseries1=np.vstack((vec_tseries1,np.array(res1)))
                newfullset1 = np.vstack((newfullset1,np.array([countlines,runindex,int(tokens1[0])])))

                vec_tseries2=np.vstack((vec_tseries2,np.array(res2)))
                newfullset2 = np.vstack((newfullset2,np.array([countlines,runindex,int(tokens2[0])])))
                
            countlines+=1

               


print(fullset1.shape,newfullset1.shape)
                                    
newselection1=newfullset1[select1]
newselection2=newfullset2[select1]
print(np.array_equal(selection1,newselection1))
print(selection1)
print(newselection1)

##length=vec_tseries.shape[0]
length = newselection1.shape[0]
self_sim=np.zeros([length,length],np.float64)
for ii in range(length):
    i=newselection1[ii,0]
    veci = vec_tseries1[i]
    normveci = np.linalg.norm(veci)
    for jj in range(length):
        j=newselection2[jj,0]
        vecj = vec_tseries2[j]
        normvecj = np.linalg.norm(vecj)
        cov = np.cov(veci,vecj)
        pearson = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
        self_sim[ii,jj] = pearson
##        self_sim[ii,jj] = np.dot(veci,vecj)/(normveci*normvecj)



#### now show the self_sim matrix
fig,axs=plt.subplots(nrows=1, ncols=1,sharex=True,sharey=True,figsize=(15,10),gridspec_kw={'width_ratios': [0.8]})
##
im1 = axs.imshow(self_sim,vmin=-1.0,vmax=1.0,aspect='auto', cmap='viridis')
axs.set_title("Constraint cross similarity (PCC)",fontsize=40)
##        
##cbticklabels= ['{:.1f}'.format(x) for x in np.arange(-1,1.1,0.2)]
###list(np.arange(-1,1,0.2))      
##cb = fig.colorbar(im1)
##cb.ax.set_yticklabels(cbticklabels,size=20)
##cb.set_label("PCC", labelpad=-1, size=25)
cb = fig.colorbar(im1)
cb.ax.tick_params(labelsize=25)
##print(xaxis,yaxis)
xaxis = list(range(0,length,int(0.1*length)))
yaxis = list(range(0,length,int(0.1*length)))
axs.set_xticks(xaxis)
axs.set_yticks(yaxis)
axs.set_xticklabels(xaxis, fontsize=20)
axs.set_yticklabels(yaxis, fontsize=20)
plt.tight_layout()
#plt.show()
plt.savefig('/home/goswam_y/Dropbox/chromosome_modelling/AnaCode_chrom_model_OOP/GD_GD_constraint_cross.png')

                
                
                

            



            
