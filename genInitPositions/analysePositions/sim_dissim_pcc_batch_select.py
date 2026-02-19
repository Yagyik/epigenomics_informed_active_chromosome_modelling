import numpy as np
import matplotlib.pyplot as plt
import sys


##### script to demonstrate the dissimilarity of configurations based on raw positions
##### as a bonus would be good to show the self-similiarity of IPD entries for same configs

dirpath = sys.argv[1]
opprefix=sys.argv[2]
blocksize=sys.argv[3] # use 500
constraintfile = "/testWhole_46c-constraints.dat"
configprefix = "/testWhole_46c-dump_"


#### we iterate over the runs and dump number -- gathering configurations 3N vec
#### store the index and timepoint -- we find the constraint matrix for the corresponding
#### we will calculate the cosine similarity matrix for constraints + expected
#### we will also calculate the cosine similarity matrix for positions (3N vec)

nruns= 16
config_tseries=np.array([],np.float64)
countlines=0
pselect = 1
fullset=np.array([],np.float64)
targetselection=60
for runindex in range(1,nruns+1):

    for timeindex in range(1000000,1990000,100000): ## 99 files
        filename = dirpath+"/"+str(runindex)+configprefix+str(timeindex)+".dat"
        print(filename)
        configfile=open(filename,"r")
        lines=configfile.readlines()
        configfile.close()
##        print(len(lines),lines[2002])
        res = []

        for line in lines[2002:]: ### chromosome positions start here
##            print(line)
            tokens=line.split()

            ### flat append xyz to get 3N vec in order
            x = float(tokens[4])
            res.append(x)
            y = float(tokens[5])
            res.append(y)
            z = float(tokens[6])
            res.append(z)
        ### got full config
##        print(len(res),res)
            

        if countlines==0:
            config_tseries=np.array(res)
            ### store the runindex and the timeindex
##            selection.append([runindex,timeindex])
            fullset=np.array([countlines,runindex,timeindex])
            ### update the probability of next selection
##            pselect = np.random.uniform(0,(targetselection - len(selection))/targetselection)
            print(fullset,pselect)
        elif countlines > 0:
            #### insert code here to make a selection from the full list
            #### need to store the runindex and the corresponding timeindex
            config_tseries=np.vstack((config_tseries,np.array(res)))
            fullset=np.vstack((fullset,np.array([countlines,runindex,timeindex])))
##            if rand < pselect:
##                
##                ### store this config
##                config_tseries=np.vstack((config_tseries,np.array(res)))
##                ### store the runindex and the timeindex
##                selection.append([runindex,timeindex])
                ### update the probability of next selection
##                probmin = min(0.5,(targetselection - len(selection))/targetselection)
##                probmax = max(0.5,(targetselection - len(selection))/targetselection)
##                pselect = np.random.uniform(probmin,probmax)
##                print(len(selection),pselect)

        countlines+=1

##print(config_tseries.shape,pselect,len(selection))

#### generate a selection
print(fullset.shape[0],pselect)
##select = np.random.choice(np.arange(fullset.shape[0]),size=targetselection,replace=False)
##final_select = np.array([],np.int64)
##for i in range(select.shape[0]-2):
##    linei = select[i]
##    print(linei)
##    final_select = np.append(final_select,np.array([linei,linei+1,linei+2]))
##    print(final_select)
##    if i == 0:
##        print(fullset[linei])
##        final_select = fullset[linei,:]
##        print(final_select)
##        final_select = np.vstack((final_select,fullset[linei+1,:]))
##        final_select = np.vstack((final_select,fullset[linei+2,:]))
##        print(fullset[linei,:])
##    else:
##        final_select = np.vstack((final_select,fullset[linei,:]))
##        final_select = np.vstack((final_select,fullset[linei+1,:]))
##        final_select = np.vstack((final_select,fullset[linei+2,:]))



##### FINAL VERSION


##final_select = np.array([],np.int64)
##while final_select.shape[0] < targetselection:
##    linei = np.random.choice(np.arange(fullset.shape[0]),size=1,replace=False)[0]
##    print("linei",linei,final_select.shape[0])
##    flag_check = 0
##    for i in range(final_select.shape[0]):
####        print(linei,final_select[i])
##        if abs(linei - final_select[i]) <= 5:
##            flag_check = 1
##            break
##    if flag_check == 0:
##        print("appending",np.array([linei,linei+1,linei+2]))
##        final_select = np.append(final_select,np.array([linei,linei+1,linei+2]))
##        print("curr shape",final_select.shape[0])


#### END FINAL VERSION

final_select = np.arange(fullset.shape[0])
    

print(final_select) 
##select = np.arange(fullset.shape[0])
##print(select)
selection = fullset[final_select]
print(selection)


##length=config_tseries.shape[0]
length = selection.shape[0]
self_sim=np.zeros([length,length],np.float64)
for ii in range(length):
    i = selection[ii,0]
    veci = config_tseries[i]
    normveci = np.linalg.norm(veci)
    for jj in range(length):
        j = selection[jj,0]
        vecj = config_tseries[j]
        normvecj = np.linalg.norm(vecj)
        cov = np.cov(veci,vecj)
        pearson = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
        self_sim[ii,jj] = pearson
##        self_sim[ii,jj] = np.dot(veci,vecj)/(normveci*normvecj)


## now show the self_sim matrix
fig,axs=plt.subplots(nrows=1, ncols=1,sharex=True,sharey=True,figsize=(15,10),gridspec_kw={'width_ratios': [0.8]})
##
im1 = axs.imshow(self_sim,vmin=0.0,vmax=1.0,aspect='auto', cmap='viridis')
axs.set_title("Batched - Raw position self similarity (PCC)",fontsize=40)
        
cb = fig.colorbar(im1)
cb.ax.tick_params(labelsize=25)
xaxis = list(range(0,length,int(0.1*length)))
yaxis = list(range(0,length,int(0.1*length)))
axs.set_xticks(xaxis)
axs.set_yticks(yaxis)
axs.set_xticklabels(xaxis, fontsize=20)
axs.set_yticklabels(yaxis, fontsize=20)
plt.tight_layout()
##plotname = '/home/goswam_y/Dropbox/chromosome_modelling/AnaCode_chrom_model/'+opprefix+'_batch_GD_pos_self.png'
plotname = opprefix+'_batch_GD_pos_self.png'
plt.savefig(plotname)
##plt.show()


#### from the selection of (runindex,timeindex) we store the corresponding constraint matrices
##print(selection)
vec_tseries = np.array([],np.float64)
nselected=0
nconstraints=int((0.5*(46+1)*(46)) + 46)
newfullset=np.array([],np.float64)


countlines=0
##nruns = 3
for runindex in range(1,nruns+1):
    filename = dirpath+"/"+str(runindex)+constraintfile
    
    cfile=open(filename,"r")
    lines=cfile.readlines()
##    print(lines)

    
    cfile.close()


    for line in lines[:10]:
        tokens = line.split()
    ##    timeindex = float(tokens[0])
        res = [eval(i) for i in tokens[1:nconstraints]]
    ##    print(tokens[0])
        if float(tokens[0]) == 1.1:
            expectvec = np.array(res)
            print(len(tokens),len(res),nconstraints)

        elif float(tokens[0]) == 1.2:
            calcflag = np.array(res)

        
    last = int(lines[-1].split()[0])
    print(filename,len(lines),last,int(len(lines)/float(blocksize)))
    loop = np.arange(1000000,last,100000)
    print("loop",loop,loop.shape,int(len(lines)/float(blocksize)))
##    raw = np.linspace(0,last,int(len(lines)/float(sys.argv[3])))
##    raw = np.linspace(0,last+int(sys.argv[3]),int(sys.argv[3]))
    raw = np.arange(0,last+int(blocksize),int(blocksize))
    print("raw",raw,raw.shape)
    
    orig_index = np.where(np.in1d(raw,loop))[0]
    print("orig_index",orig_index)
    print("loop index",loop)
##    orig_index = loop
    print("set",raw[orig_index])
##        print(len(lines),lines[2002])
##    for i in orig_index:
    i=0
    for line in lines[2+orig_index[0]:]:
        if i >= len(loop):
            break
##        line = lines[i]
        tokens=line.split()
##        print("line",i,tokens[0],loop[i])
        timeindex = int(tokens[0])
        res = [eval(i) for i in tokens[1:nconstraints]]

        ### get this out of the way
        if float(tokens[0]) == 1.1:
            expectvec=np.array(res)
        elif float(tokens[0]) == 1.2:
            calcflag = np.array(res)

        

        ### store all like earlier but keeping in mind the time interval for the configs
##        if int(tokens[0])>=5000000 and int(tokens[0])<9990000 and int(tokens[0])%500000==0:
        if int(tokens[0]) == loop[i]:
            print("adding",i,tokens[0],loop[i])
            if countlines==0:
                vec_tseries=np.array(res)
                newfullset = np.array([countlines,runindex,int(tokens[0])])
            elif countlines>0:
                vec_tseries=np.vstack((vec_tseries,np.array(res)))
                newfullset = np.vstack((newfullset,np.array([countlines,runindex,int(tokens[0])])))
            countlines+=1
            i+=1


print(fullset.shape,newfullset.shape)
                                    
newselection=newfullset[final_select]
print(np.array_equal(selection,newselection))
print(selection)
print(newselection)

##length=vec_tseries.shape[0]
length = newselection.shape[0]
self_sim=np.zeros([length,length],np.float64)

for ii in range(length):
    i=newselection[ii,0]
##    veci = vec_tseries[i]
    veci = np.zeros_like(calcflag)
    veci[vec_tseries[i,:] > expectvec - 0.1] = 1
    veci[calcflag==0] = 0
    normveci = np.linalg.norm(veci)
    for jj in range(length):
        j=newselection[jj,0]
##        vecj = vec_tseries[j]
        vecj = np.zeros_like(calcflag)
        vecj[vec_tseries[j,:] > expectvec - 0.1] = 1
        vecj[calcflag==0] = 0
        normvecj = np.linalg.norm(vecj)

        
        cov = np.cov(veci,vecj)
        pearson = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
        self_sim[ii,jj] = pearson
##        self_sim[ii,jj] = np.dot(veci,vecj)/(normveci*normvecj)



#### now show the self_sim matrix
fig,axs=plt.subplots(nrows=1, ncols=1,sharex=True,sharey=True,figsize=(15,10),gridspec_kw={'width_ratios': [0.8]})
##
im1 = axs.imshow(self_sim,vmin=0.0,vmax=1.0,aspect='auto', cmap='viridis')
axs.set_title("Batched - constraint self similarity (PCC)",fontsize=40)
        
cb = fig.colorbar(im1)
cb.ax.tick_params(labelsize=25)
xaxis = list(range(0,length,int(0.1*length)))
yaxis = list(range(0,length,int(0.1*length)))
axs.set_xticks(xaxis)
axs.set_yticks(yaxis)
axs.set_xticklabels(xaxis, fontsize=20)
axs.set_yticklabels(yaxis, fontsize=20)
plt.tight_layout()
##plotname = '/home/goswam_y/Dropbox/chromosome_modelling/AnaCode_chrom_model/'+opprefix+'_batch_GD_constraint_self.png'
plotname = opprefix+'_batch_GD_constraint_self.png'
plt.savefig(plotname)
##plt.show()
        
                
                
                

            



            
