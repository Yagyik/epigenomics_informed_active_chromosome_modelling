import numpy
import matplotlib.pyplot as plt
import sys
from sklearn.metrics.pairwise import cosine_similarity



matfile = sys.argv[1] + "/testWhole_4c-constraintMat.dat"
f = open(matfile,"r")

lines= f.readlines()

x=numpy.array([1,2,3,4,5])
y=numpy.array([1,2,3,4,5])
z=numpy.zeros([5,5],numpy.float64)
z2=numpy.array([[0.65,0.4,0.35,0.3,0.3],
                [0.4,0.5,0.2,0.15,0.25],
                [0.35,0.2,0.25,0.15,0.2],
                [0.3,0.15,0.15,0.2,0.15],
                [0.3,0.25,0.2,0.15,0]])
for line in lines:
    tokens=line.split()
    ci = int(tokens[0])
    cj = int(tokens[1])
    z[ci,cj] = float(tokens[2])*250

# read array now plot
print(z)
xmesh,ymesh = numpy.meshgrid(x,y)

##fig = plt.figure()
##ax1 = plt.contourf(xmesh,ymesh,z)
##plt.subplot(1,2,1)


##fig,axs=plt.subplots(nrows=1, ncols=2,sharex=True,sharey=True,figsize=(15,5),gridspec_kw={'width_ratios': [0.8, 1]})
##
##im1 = axs[0].imshow(z2,vmin=0,vmax=0.8,aspect='auto', cmap='viridis')
##axs[0].set_title("Input",fontsize=40)
##
##im2 = axs[1].imshow(z,vmin=0, vmax=0.8,aspect='auto', cmap='viridis')
##axs[1].set_title("Output",fontsize=40)
z1_flat= z.flatten()
z2_flat = z2.flatten()

z3_flat = numpy.correlate(z1_flat,z2_flat,"same")

z3 = numpy.reshape(z3_flat,(5,5))
cos_sim = cosine_similarity(z,z2)


print(cos_sim)
fig,axs=plt.subplots(nrows=1, ncols=1,sharex=True,sharey=True,figsize=(15,5),gridspec_kw={'width_ratios': [0.8]})

##im1 = axs.imshow(abs(z2-z)/z2,vmin=0,vmax=0.4,aspect='auto', cmap='viridis')
##axs.set_title("L1 error",fontsize=40)

im1 = axs.imshow(cos_sim,vmin=0.8,vmax=1.0,aspect='auto', cmap='viridis')
axs.set_title("cosine similarity",fontsize=40)

##plt.subplot(1,2,2)
##plt.colorbar()
fig.colorbar(im1)
plt.show()
    
