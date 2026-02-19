import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx
import sys
import subprocess
import itertools
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import plotly.graph_objects as go

from deform_grid_utils import *
from placement_utils import *
from maxHole_utils import *

    

###### execution

conffile = sys.argv[1]
nsurf = int(sys.argv[2])
ncyc= int(sys.argv[3])
nGrid = int(sys.argv[4])
infile = open(conffile,"r")
lines = infile.readlines()
ella = float(sys.argv[5])
ellb = float(sys.argv[6])
ellc = float(sys.argv[7])
sfneighfact = float(sys.argv[8])
outfilename=sys.argv[9]
##conffile = 'test_surfs/Conf-fixedSurfFile_ea5.5_eb5.5_ec5.5.dat'
##nsurf = 2500
##infile = open(conffile,"r")
##lines = infile.readlines()
##ella = 5.5
##ellb = 5.5
##ellc = 5.5



surfx = np.zeros([nsurf],np.float64)
surfy = np.zeros([nsurf],np.float64)
surfz = np.zeros([nsurf],np.float64)
surfid = np.zeros([nsurf],np.int64)
surfrad = 0.0
sfid = 0
count = 0
for line in lines:
    count +=1
    tokens = line.split()
    sfid = int(tokens[0])
    surfid[sfid] = sfid
    surfrad = float(tokens[1])
    surfx[sfid] = float(tokens[2])
    surfy[sfid] = float(tokens[3])
    surfz[sfid] = float(tokens[4])


surfxp=surfx.reshape((1,surfx.shape[0]))
surfyp=surfy.reshape((1,surfy.shape[0]))
surfzp=surfz.reshape((1,surfz.shape[0]))
surfpoints=np.hstack((surfxp.T,surfyp.T,surfzp.T))
print(surfpoints.shape)

oldsurfpoints = np.copy(surfpoints)

## define a grid for current geometry and points
old_grid_indices,old_gridcent_surfpoints = grid(oldsurfpoints,ella,ellb,ellc,
                                        ella,ellb,ellc,surfrad,nGrid)
print("made basic grid",old_grid_indices.shape)
## create a modification schedule -- starting with base geometry and different mods
if ella > ellc:
    final = [1.5,1.5,1.5]
if ella == ellc:
    final = [1.5,1,0.5]

if ncyc == 2000:   
    modifications = [[[1.5,1,1],200],
                     [[1,1.5,1],200],
                     [[1,1,1.5],200],
                     [[1.5,1.5,1],200],
                     [[1.5,1,1.5],200],
                     [[1,1.5,1.5],200],
                     [[1,1,1],200],
                     [final,600]]
if ncyc == 1000:
    modifications = [[[1.5,1,1],100],
                     [[1,1.5,1],100],
                     [[1,1,1.5],100],
                     [[1.5,1.5,1],100],
                     [[1.5,1,1.5],100],
                     [[1,1.5,1.5],100],
                     [[1,1,1],100],
                     [final,300]]
if ncyc == 500:
    modifications = [[[1.5,1,1],50],
                     [[1,1.5,1],50],
                     [[1,1,1.5],50],
                     [[1.5,1.5,1],50],
                     [[1.5,1,1.5],50],
                     [[1,1.5,1.5],50],
                     [[1,1,1],50],
                     [final,150]]

totcyc = np.sum(np.array([mod[1] for mod in modifications[:]]))
if totcyc != ncyc:
    print("budget vs tot points mismatch -- provided",ncyc," but budget comes to ", totcyc)
    exit(0)
## each mod has an associated budget of points to add


accumcyc = 0
accumsurfpoints = np.copy(oldsurfpoints)
grid_indices = np.copy(old_grid_indices)
gridcent_surfpoints = np.copy(old_gridcent_surfpoints)
accumtris = np.array([],np.float64)
accumindices = np.array([],np.float64)
## for each modification

for m in range(len(modifications)):
    mod = modifications[m][0]
    neighscale = max(mod[0],mod[1])
    neighscale = max(neighscale,mod[2])

    ## deform the ellipsoid -- generate new positions
    print("deforming main set",mod,old_grid_indices.shape)
    modsurfpoints = deform_from_grid(oldsurfpoints,old_grid_indices,
                                     old_gridcent_surfpoints,surfrad,ella,ellb,ellc,mod,
                                     nGrid)
    
    
    ## find gaps in those positions, add points till budget filled
    pt_budget = modifications[m][1]
    print("extracting cycles")
##    coms,tris,indices = get_new_points_tris(modsurfpoints,nsurf,surfrad,
##                                            ella,ellb,ellc,mod)
    print("deforming full set",mod)
    modaccum = deform_from_grid(accumsurfpoints,grid_indices,
                                gridcent_surfpoints,surfrad,ella,ellb,ellc,mod,
                                nGrid)
    

##    sorted_coms,sorted_indices,sorted_tris = get_points_tris_maxhole(modsurfpoints,neighscale*sfneighfact,
##                                                nsurf,surfrad,pt_budget)

    
    sorted_coms,sorted_indices,sorted_tris = get_points_tris_maxhole(modaccum,neighscale*sfneighfact,
                                                modaccum.shape[0],surfrad,pt_budget)
    
    
    #### filter these by overlaps with existing points

    ### first modify all points for current deformation -- takes undeformed accumpoints as input
    
    ## compare sorted_coms with these modified points to filter
##    print("extracting top tris",pt_budget)
##    sorted_coms,sorted_tris,sorted_indices = filter_tris(coms,tris,indices,
##                                                         modaccum,surfrad,pt_budget)
##    print(sorted_tris[-1],"smallest of the extracted triangles")
              
    ## new array with extra sizing -- takes deformed accum and new points -- returns deformed
    
    print("old shape",accumsurfpoints.shape,totcyc,pt_budget)
    accumsurfpoints = np.vstack((modaccum,sorted_coms))
    print("new shape",accumsurfpoints.shape,totcyc,pt_budget)

    #### write this to file
    name = conffile[:-4] +"_"+str(mod[0]*ella)+"_"+str(mod[1]*ellb)+"_"+str(mod[2]*ellc)+"_hole_ovito.dat"

    outfile = open(name,"w")
    outfile.write("%d\n" %(accumsurfpoints.shape[0]))
    outfile.write("spam\n")
    for i in range(nsurf):
        outfile.write("%d %d %f %f %f %f %f\n" %(i,1,0.2,surfrad,accumsurfpoints[i,0],accumsurfpoints[i,1],accumsurfpoints[i,2]))

    for i in range(nsurf,accumsurfpoints.shape[0]):
        outfile.write("%d %d %f %f %f %f %f\n" %(i,2,0.5,surfrad,accumsurfpoints[i,0],accumsurfpoints[i,1],accumsurfpoints[i,2]))

    outfile.close()

    ### regrid with true geom -- grid always uses baseline geom as reference
    grid_indices,gridcent_surfpoints = grid(accumsurfpoints,ella,ellb,ellc,
                                                 mod[0]*ella,mod[1]*ellb,mod[2]*ellc,surfrad,nGrid)

    print("grid for all",grid_indices.shape)

    ### un-deform for next round
    print("deformed",accumsurfpoints[0,:])
    accumsurfpoints = deform_from_grid(accumsurfpoints,grid_indices,
                                      gridcent_surfpoints,surfrad,ella,ellb,ellc,[1,1,1],
                                      nGrid)
    print("undeformed",accumsurfpoints[0,:])
    print("reversed deformation",accumsurfpoints.shape,mod)

    #### write this to file
    name = conffile[:-4] +"_"+str(mod[0]*ella)+"_"+str(mod[1]*ellb)+"_"+str(mod[2]*ellc)+"_hole_undeform_ovito.dat"

    outfile = open(name,"w")
    outfile.write("%d\n" %(accumsurfpoints.shape[0]))
    outfile.write("spam\n")
    for i in range(nsurf):
        outfile.write("%d %d %f %f %f %f %f\n" %(i,1,0.2,surfrad,accumsurfpoints[i,0],accumsurfpoints[i,1],accumsurfpoints[i,2]))

    for i in range(nsurf,accumsurfpoints.shape[0]):
        outfile.write("%d %d %f %f %f %f %f\n" %(i,2,0.5,surfrad,accumsurfpoints[i,0],accumsurfpoints[i,1],accumsurfpoints[i,2]))

    outfile.close()

    
    
    accumtris = np.append(accumtris,sorted_tris)
    accumcyc +=  pt_budget
    if m==0:
        accumindices = np.copy(sorted_indices)
        
    else:
        accumindices = np.vstack((accumindices,sorted_indices))


    name = conffile[:-4] +"_"+str(mod[0]*ella)+"_"+str(mod[1]*ellb)+"_"+str(mod[2]*ellc)+"_hole_final_ovito.dat"

    outfile = open(name,"w")
    outfile.write("%d\n" %(accumsurfpoints.shape[0]))
    outfile.write("spam\n")
    for i in range(nsurf):
        outfile.write("%d %d %f %f %f %f %f\n" %(i,1,0.2,surfrad,accumsurfpoints[i,0],accumsurfpoints[i,1],accumsurfpoints[i,2]))

    for i in range(nsurf,accumsurfpoints.shape[0]):
        outfile.write("%d %d %f %f %f %f %f\n" %(i,2,0.5,surfrad,accumsurfpoints[i,0],accumsurfpoints[i,1],accumsurfpoints[i,2]))

    outfile.close()
        
    
    





##print(sorted_tri)
##print(sorted_com)
##print(sorted_indices)
##### print out the first ncyc
outfile=open(outfilename,"w")
print(outfilename,ncyc,nsurf,accumsurfpoints.shape,accumindices.shape)
for i in range(ncyc):
    j = i+nsurf
    print(i,nsurf,j)
    outfile.write("%d %f %f %f %f %d %d %d\n" %(i,surfrad,accumsurfpoints[j][0],accumsurfpoints[j][1],accumsurfpoints[j][2],
                                             accumindices[i][0],accumindices[i][1],accumindices[i][2]))

outfile.close()


#### write this to file
name = conffile[:-4] +"_"+str(mod[0]*ella)+"_"+str(mod[1]*ellb)+"_"+str(mod[2]*ellc)+"_hole_final_ovito.dat"

outfile = open(name,"w")
outfile.write("%d\n" %(accumsurfpoints.shape[0]))
outfile.write("spam\n")
for i in range(nsurf):
    outfile.write("%d %d %f %f %f %f %f\n" %(i,1,0.2,surfrad,accumsurfpoints[i,0],accumsurfpoints[i,1],accumsurfpoints[i,2]))

for i in range(nsurf,accumsurfpoints.shape[0]):
    outfile.write("%d %d %f %f %f %f %f\n" %(i,2,0.5,surfrad,accumsurfpoints[i,0],accumsurfpoints[i,1],accumsurfpoints[i,2]))

outfile.close()

##tri_indices = np.where(tri_area_arr > np.pi*surfrad*surfrad)[0]
##print(tri_indices)
##print(np.where(tri_area_arr > np.pi*surfrad*surfrad))
##print(com_arr[tri_indices])
##
##
##filt_tri = tri_area_arr[tri_indices]
##filt_com = com_arr[tri_indices]
##sorted_tri_indices = filt_tri.argsort()
##sorted_arr1 = filt_tri[sorted_tri_indices[::-1]]
##sorted_arr2 = filt_com[sorted_tri_indices[::-1]]
##
##print(sorted_tri_indices,len(sorted_tri_indices))


        

