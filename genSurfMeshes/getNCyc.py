import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx
import sys
import subprocess
import itertools
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import plotly.graph_objects as go



def gendistmat(surfpoints,nsurf):
    distmat = np.zeros([nsurf,nsurf],np.float64)

    for i in range(nsurf-1):
        for j in range(i+1,nsurf):
            fpt = surfpoints[i,:]
            spt = surfpoints[j,:]

    ##            print(fpt,spt)

            distmat[i,j] = np.linalg.norm(fpt-spt)

            distmat[j,1] = np.linalg.norm(fpt-spt)
            ### make sure we dont consider zero entries (diagonals)
    dist1d = np.ravel(distmat)
    return distmat

def genneighlist(surfpoints,nsurf,surfrad):

    neighlist = np.zeros([nsurf,nsurf+1],np.int64)

    for i in range(nsurf-1):
        for j in range(i+1,nsurf):
            fpt = surfpoints[i,:]
            spt = surfpoints[j,:]

    ##            print(fpt,spt)
            dist = np.linalg.norm(fpt - spt)

            if dist < 4.0*surfrad:
                neighlist[i,0] +=1
##                print(neighlist[i,0])
                neighlist[i,int(neighlist[i,0])] = j

                neighlist[j,0] +=1
                neighlist[j,int(neighlist[j,0])] = i

##    print(neighlist)

    return neighlist

def genadjlist(surfpoints,nsurf,distcut):

    neighlist = np.zeros([nsurf,nsurf+1],np.int64)

    for i in range(nsurf-1):
        for j in range(i+1,nsurf):
            fpt = surfpoints[i,:]
            spt = surfpoints[j,:]

    ##            print(fpt,spt)
            dist = np.linalg.norm(fpt - spt)

            if dist < distcut:
                neighlist[i,0] +=1
##                print(neighlist[i,0])
                neighlist[i,int(neighlist[i,0])] = j

                neighlist[j,0] +=1
                neighlist[j,int(neighlist[j,0])] = i

##    print(neighlist)

    return neighlist



###### execution

conffile = sys.argv[1]
nsurf = int(sys.argv[2])
ncyc= int(sys.argv[3])
infile = open(conffile,"r")
lines = infile.readlines()
ella = float(sys.argv[4])
ellb = float(sys.argv[5])
ellc = float(sys.argv[6])
sfneighfact = float(sys.argv[7])
outfilename=sys.argv[8]
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


distmat = gendistmat(surfpoints,nsurf)
print("distmat",distmat.shape)

neighlist = genneighlist(surfpoints,nsurf,surfrad)

print("neighlist",neighlist.shape)

### get characteristic distances
distlist = []
for i in range(nsurf):
    si = surfpoints[i]
    
    for j in range(1,neighlist[i,0]+1):
        sj = surfpoints[neighlist[i,j]]
#         print(i,neighlist[i,0],neighlist[i,j])
        
        if neighlist[i,j]>i:
            distlist.append(np.linalg.norm(sj-si))
            
# print(len(distlist))

distarr = np.array(distlist)

##counts,bins = np.histogram(distarr,bins=40)
##print(bins[0],bins[-1],bins.shape,surfrad)
##
##shells = np.zeros([bins.shape[0]-1],np.float64)
##for i in range(bins.shape[0]-1):
##    shells[i] = 4.0*np.pi/3*(bins[i+1]**3 - bins[i]**3)
###     print(bins[i],bins[i+1],shells[i])
##plt.plot(bins[:-1],counts/(nsurf*shells))

print("radius and factor",surfrad,2*sfneighfact*surfrad)
adj_graph = genadjlist(surfpoints,nsurf,2*sfneighfact*surfrad)

print("adj graph",adj_graph.shape,np.min(adj_graph[:][0]),np.mean(adj_graph[:][0]))

### now make a graph of connected surf points

# G = nx.DiGraph()
G = nx.Graph()
G.add_nodes_from([(i,{"pos":[surfpoints[i]]}) for i in range(nsurf)])
##print(G)
# print(G.nodes())
print("node 0 position",G.nodes[0]["pos"])
for i in range(nsurf):
    n_adj = adj_graph[i,0]
    adj_list = adj_graph[i,1:n_adj+1]
    
    for j in adj_list:
        if j > i:
            G.add_edge(i,j)

##print(G.edges())


DG = G.to_directed()
cycles = nx.simple_cycles(DG,3)
listcyc = sorted(cycles)
print("list cyc",len(listcyc))
##print(G)

##print(len(listcyc))
##real_cyc=[]
####for cyc in listcyc:
####    if len(cyc) > 2 and cyc[-1] > cyc[-2]:
#####         print(cyc)
####        real_cyc.append([cyc])
##
##
##print(len(listcyc))



real_cyc=[]
for cyc in listcyc:
    if len(cyc) > 2 and cyc[-1] > cyc[-2]:
        
        ### check that these 3 don't share a common neighbour
#         print(cyc,cyc[0])
#         print(adj_graph[cyc[0],:])
        row1 = adj_graph[cyc[0],1:adj_graph[cyc[0],0]+1]
#         print(row1)
        row2 = adj_graph[cyc[1],1:adj_graph[cyc[1],0]+1]
#         print(row2)
        row3 = adj_graph[cyc[2],1:adj_graph[cyc[2],0]+1]
#         print(row3)
        
        
        ### make sure now that no common neighbours are present (that are not each other)
        common1 = np.intersect1d(row1,row2)
        common2 = np.intersect1d(row1,row3)
        
#         print(common1)
#         print(common2)
        
#         print(np.intersect1d(common1,common2))
        
        totcommon = np.intersect1d(common1,common2)
        
        if len(totcommon) <=3:
            ### check positions
#             print("Yaay",cyc)
#             print(surfpoints[cyc[0]])
#             print(surfpoints[cyc[1]])
#             print(surfpoints[cyc[2]])
            
#             print(totcommon)
            
#             print(surfpoints[totcommon])
            real_cyc.append([cyc])
            
            v1 = surfpoints[cyc[1]] - surfpoints[cyc[0]]
            v2 = surfpoints[cyc[2]] - surfpoints[cyc[0]]
            
            tri_arr = 0.5*np.linalg.norm(np.cross(v1,v2))
            
#             print(tri_arr,surfrad,np.pi*surfrad*surfrad)
            
print("real cyc",len(real_cyc))
coms = []
tri_areas = []
tri_indices = []
for cyc in real_cyc:
    pts = surfpoints[cyc][0]
#     print(cyc,pts.shape)
#     print(pts)
    com = np.mean(pts,axis=0)
#     print("com",com,surfpoints[0])
    pti = [surfpoints[i] for i in cyc]
    
#     for pt in pti[0]:
#         print(pt)
    coms.append(com)
    
    tri_area = 0.5*np.linalg.norm(np.cross(pti[0][1]-pti[0][0],pti[0][2]-pti[0][0]))
    tri_areas.append(tri_area)
    tri_indices.append(np.array(cyc[0]))
    dists = [np.linalg.norm(pt - com) for pt in pti[0]]
#     print("dists",dists,pti[0].shape)

com_arr = np.array(coms)
tri_area_arr = np.array(tri_areas)
tri_indices_arr = np.array(tri_indices)
print("com tri",com_arr.shape,tri_area_arr.shape)

#### now sort triangle areas, drop the first ncyc

sorted_tri_indices = tri_area_arr.argsort()

print(sorted_tri_indices)

sorted_tri = tri_area_arr[sorted_tri_indices[::-1]]
sorted_com = com_arr[sorted_tri_indices[::-1]]
sorted_indices = tri_indices_arr[sorted_tri_indices[::-1]]

print(sorted_tri)
print(sorted_com)
print(sorted_indices)
##### print out the first ncyc
outfile=open(outfilename,"w")

for i in range(ncyc):
    outfile.write("%d %f %f %f %f %d %d %d\n" %(i,surfrad,sorted_com[i][0],sorted_com[i][1],sorted_com[i][2],
                                             sorted_indices[i][0],sorted_indices[i][1],sorted_indices[i][2]))

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


        

