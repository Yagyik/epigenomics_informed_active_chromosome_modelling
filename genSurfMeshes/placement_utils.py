import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx
import itertools


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

            if dist < 8.0*surfrad:
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


def proximity_array(exist_surfpoints,candidate_points, surfrad):

    ## check that generated candidate points are not too close to existing points
    ## avoids redundancy -- make a proximity array size of candidate_points -- 0 where no proximity
    outlen = candidate_points.shape[0]
    proximity_arr = np.zeros([outlen],np.int64)


    for c in range(outlen):
        cand_pt = candidate_points[c]
        dist_arr = np.array([np.linalg.norm(ept - cand_pt) for ept in exist_surfpoints])
        proximity_arr[c] = np.any(dist_arr < 0.5*surfrad)

    return proximity_arr
    
