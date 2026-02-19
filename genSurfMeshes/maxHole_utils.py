import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import sys
import subprocess
import itertools
from scipy.spatial import ConvexHull, convex_hull_plot_2d

from placement_utils import *


def get_new_points_tris(surfpoints,nsurf,surfrad,sfneighfact,ella,ellb,ellc,mod):

    new_ella = mod[0]*ella
    new_ellb = mod[1]*ellb
    new_ellc = mod[2]*ellc

    neighscale = max(mod[0],mod[1])
    neighscale = max(neighscale,mod[2])

    distmat = gendistmat(surfpoints,nsurf)
    print("distmat",distmat.shape)

    neighlist = genneighlist(surfpoints,nsurf,surfrad*neighscale)

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


    print("radius and factor",surfrad,2*sfneighfact*surfrad)
    adj_graph = genadjlist(surfpoints,nsurf,2*sfneighfact*surfrad)

    print("adj graph",adj_graph.shape,np.min(adj_graph[:][0]),np.mean(adj_graph[:][0]))

    ### now make a graph of connected surf points

    G = nx.Graph()
    G.add_nodes_from([(i,{"pos":[surfpoints[i]]}) for i in range(nsurf)])

    print("node 0 position",G.nodes[0]["pos"])
    for i in range(nsurf):
        n_adj = adj_graph[i,0]
        adj_list = adj_graph[i,1:n_adj+1]
        
        for j in adj_list:
            if j > i:
                G.add_edge(i,j)



    DG = G.to_directed()
    cycles = nx.simple_cycles(DG,3)
    listcyc = sorted(cycles)
    print("list cyc",len(listcyc))




    real_cyc=[]
    for cyc in listcyc:
        if len(cyc) > 2 and cyc[-1] > cyc[-2]:
            
            ### check that these 3 don't share a common neighbour
            row1 = adj_graph[cyc[0],1:adj_graph[cyc[0],0]+1]
            row2 = adj_graph[cyc[1],1:adj_graph[cyc[1],0]+1]
            row3 = adj_graph[cyc[2],1:adj_graph[cyc[2],0]+1]

            
            
            ### make sure now that no common neighbours are present (that are not each other)
            common1 = np.intersect1d(row1,row2)
            common2 = np.intersect1d(row1,row3)
            
            totcommon = np.intersect1d(common1,common2)
            
            if len(totcommon) <=3:
                ### check positions
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

    

##    outfile.write("%d %f %f %f %f %d %d %d\n" %(i,surfrad,sorted_com[i][0],sorted_com[i][1],sorted_com[i][2],
##                                             sorted_indices[i][0],sorted_indices[i][1],sorted_indices[i][2]))


    return com_arr,tri_area_arr,tri_indices_arr




def filter_tris(coms_arr,tri_area_arr,tri_indices_arr,
                exist_surfpoints,surfrad,npts):
    ## filter out those that are too close to existing points
    proximity_arr = proximity_array(exist_surfpoints,coms_arr,surfrad)
    
    filt_com_arr = coms_arr[proximity_arr == 0]
    filt_tri_area_arr = tri_area_arr[proximity_arr == 0]
    filt_tri_indices_arr = tri_indices_arr[proximity_arr == 0]
    print("filtered com tri",coms_arr.shape,tri_area_arr.shape)

    #### now sort triangle areas, drop the first ncyc 

    sorted_tri_indices = tri_area_arr.argsort()

    print(sorted_tri_indices)

    sorted_tri = tri_area_arr[sorted_tri_indices[::-1]]
    sorted_com = coms_arr[sorted_tri_indices[::-1]]
    sorted_indices = tri_indices_arr[sorted_tri_indices[::-1]]

    return sorted_com[:npts],sorted_tri[:npts],sorted_indices[:npts]



def checkOthers(pi,pj,neighlist,surfpoints,rad):

    midpt = 0.5*(surfpoints[pi] + surfpoints[pj])
    flag = 0

    neighsi = neighlist[pi,1:neighlist[pi,0]+1]
    neighsj = neighlist[pj,1:neighlist[pj,0]+1]

    allneighs= list(neighsi) + list(neighsj)
##    print(allneighs,pi,pj)

    
    
    for ni in allneighs:

        if ni == pi or ni == pj:
            continue

        dist = np.linalg.norm(surfpoints[ni] - midpt)

        if dist < 0.6*rad:
##            print("non-simplex neighbour here",ni,surfpoints[ni],midpt,dist)
            flag = 1
            break
##    print(pi,pj,midpt,flag)
    return flag


def findArea(vertices,rad):
##    print(vertices.shape)
    ab = vertices[1,:] - vertices[0,:]
    ac = vertices[2,:] - vertices[0,:]
##    print(ab.shape,ac.shape)
    area = 0.5*np.linalg.norm(np.cross(ab,ac))

    pmab = vertices[0,:] - vertices[1,:]
    pmac = vertices[2,:] - vertices[1,:]
    pmarea = 0.5*np.linalg.norm(np.cross(pmab,pmac))   
    centroid = np.mean(vertices,axis=0)
    dcent = np.linalg.norm(vertices - centroid,axis=1)
    
    return area

def tris_from_pts(vertices,surfpoints,surfrad):

    if len(vertices) == 3:
        return vertices

    
    elif len(vertices) > 3:

        trialtris = list(itertools.combinations(vertices,3))

                
        areas = []
        tris = []
        for trialtri in trialtris:
            trtr = np.array(trialtri)
##            print(trtr,surfpoints[trtr,:])
            spamarea = findArea(surfpoints[trtr,:],surfrad)
            areas.append(spamarea)
            tris.append(trtr)

        ## find max
        maxarr_tri = tris[np.argmax(np.array(areas))]

        return maxarr_tri,areas[np.argmax(np.array(areas))]
            



def maxHole(surfpoints,neighscale,nsurf,surfrad,toreplace,meansep):

    print("replacing!",toreplace,surfpoints.shape)


    surfhull=ConvexHull(surfpoints)

    neighlist = genneighlist(surfpoints,nsurf,surfrad*neighscale)
    
    tris =[]
    j=0 ### reset index because previous loop
    for s in surfhull.simplices:
        tris.extend([(s[0],s[1],s[2])])

    utri = list(set(tris))  #utri is the list of unique triangles


    #### checking triangles that have toreplace
    autri = np.asarray(utri)
##    print("participating tris",np.where(autri == toreplace))
##    print("old triangles")
##    print(autri[np.where(autri == toreplace)[0]])
    

    distarr = np.array([],np.float64)
    pointlist = []

    for tri in utri:
        ### each triangle -- indices composing
        ltri = list(tri)
##        print(ltri)

        ## list of edges
        combs = list(itertools.combinations(ltri,2))

##        print(list(combs))
        spamdistarr = np.array([],np.float64)

        ## find positions of edge ends and get distances
        for comb in combs:
            pari= comb[0]
            posi = surfpoints[pari,:]

            parj= comb[1]
            posj = surfpoints[parj,:]

            spam = np.sqrt(np.sum((posj - posi)**2))

            dist = np.linalg.norm(posj - posi)
            spamdistarr = np.append(spamdistarr,dist)

            

        ### get the participants of the largest edge, provided it is > 1.5*rad
        ### can definitely squeeze in one particle between these two --
        ### provided that the resultant transverse crowding is not too high!
        if np.max(spamdistarr) > 1.5*surfrad:
##            print(comb,spamdistarr)
            (pi,pj) = combs[np.argmax(spamdistarr)]
##            print(surfpoints[pi],surfpoints[pj])

            ### checkOthers returns a flag
            othersFlag = checkOthers(pi,pj,neighlist,surfpoints,surfrad) ### check if midpoint is close to any others (not part of simplices)
            if othersFlag == 0:
                
                pointlist.append((pi,pj))
                distarr = np.append(distarr,np.linalg.norm(surfpoints[pi] - surfpoints[pj]))
            
    ### now we have array of all the filtered distances and their points
##    print(pointlist)
##    print(distarr)

    print("num candidate points",len(pointlist),distarr.shape)
    if len(pointlist) <= 0:
        print("ran out of holes")
        return np.array([-nsurf,-nsurf,-nsurf])
    ### for each of these pairs of points, complete the quadrilateral and find the mean
    rndpt = np.random.choice(np.arange(len(pointlist)),1)[0]
    print(rndpt,pointlist[0])
    pi = pointlist[rndpt][0]
    pj = pointlist[rndpt][1]
    meanpos = 0.5*(surfpoints[pi,:]+surfpoints[pj,:])
    poslist = []
    avgdistlist = []
    arealist = []
    sharedlist = []
    while len(avgdistlist) <= 0 and meansep > 0.01:
        accept_counter = 0
        reject_counter = 0
        full_counter = 0
        for pts in pointlist:
##            print("pts",pts,len(pts),full_counter)
            full_counter += 1
            pi = pts[0]
            pj = pts[1]

            ### find triangles (indices thereof) that include these two points
            pitris = np.where(autri == pi)
            pjtris = np.where(autri == pj)


            ### get the set of unique points from within the two sets of triangles -- will have at least one overlap each (pi in pjtris and pj in pitris)
            piuniq = np.unique(autri[pitris[0]])
            pjuniq = np.unique(autri[pjtris[0]])

            sharedpts = []
            fullpts = list(piuniq) + list(pjuniq)

            fullpts_trim = np.unique(np.array(fullpts))
        #     print(fullpts,fullpts_trim)

            ## shared points will store the set of points that participate in triangles containing pi and triangles containing pj -- 3rd vertices in the triangle including pi,pj
            for i in fullpts_trim:
                if np.isin(i,piuniq) and np.isin(i,pjuniq):
                    sharedpts.append(i)

##            print("quadrangle",sharedpts)
            while len(sharedpts) < 4:
                print("triangle subset i",autri[pitris[0]])
                print("triangle subset j",autri[pjtris[0]])
                print("points shared by two triangles -- should number 4 (incl pi,pj)",sharedpts,len(sharedpts),pi,pj)
                for i in fullpts_trim: ### emergency add more points to make something happen
                    if np.isin(i,piuniq) or np.isin(i,pjuniq):
                        sharedpts.append(i)
            
            ### for this, get the mean pos, dist of mean from each.
            ### check if dist of mean from each > 0.8 sigma for all
            meanpos = np.mean(surfpoints[sharedpts],axis=0)
            distlist = [np.linalg.norm(x) for x in surfpoints[sharedpts] - np.mean(surfpoints[sharedpts],axis=0)]
            
##            if np.all(np.array(distlist)>meansep*surfrad):
            if len(np.where(np.array(distlist)>= meansep*surfrad)[0]) >= 2:

                ### do a sanity check to ensure this pos not anywhere near our particles before appending
                testpos = meanpos.reshape(1,3)
                normarr = np.linalg.norm(surfpoints-testpos,axis=1)
    ##            print("checking",normarr.shape,testpos,surfpoints[0],surfpoints[0]-testpos,normarr[0])
    ##            print(normarr)
                if np.any(normarr<0.5*surfrad):             
##                    print("far too close",np.where(normarr<surfrad),testpos,surfpoints[np.where(normarr<surfrad)])
                    continue
                poslist.append(meanpos)
                avgdistlist.append(np.mean(np.array(distlist)))
                sharedlist.append(tuple(sharedpts))
                trialtris = list(itertools.combinations(sharedpts,3))
            
                
                spamarea = 0
                counttris = 0
                for trialtri in trialtris:
                    trtr = np.array(trialtri)
                    if np.isin(pi,trtr) and np.isin(pj,trtr):

                        
                        ### sum their two areas
                        spamarea += findArea(surfpoints[trtr],surfrad)
                        counttris +=1
                
        #         print(counttris,spamarea)
                if counttris < 2:
                    print("WHOA WHOA!!!!")
                    print(pi,pj)
                    print(sharedpts)
                    print(trialtris)
                    exit(0)
                arealist.append(spamarea)
                accept_counter += 1
            else:
                reject_counter += 1
                

##        if meansep < 0.5*surfrad:
##            rndpt = np.random.choice(np.arange(len(pointlist)),1)
##            pi = pointlist[rndpt][0]
##            pj = pointlist[rndpt][1]
##            meanpos = 0.5*(surfpoints[pi,:]+surfpoints[pj,:])
##            print(meanpos)
##            poslist.append(meanpos)
##            distlist = [np.linalg.norm(meanpos - surfpoints[pi,:]),np.linalg.norm(meanpost - surfpoints[pj,:])]
##            avgdistlist.append(np.mean(np.array(distlist)))
            
            
        if len(avgdistlist) <= 1:
            meansep = 0.95*meansep
            print("final candidate hole number",len(avgdistlist),meansep)
        else:
            print("final candidate hole number",len(avgdistlist),meansep,np.min(avgdistlist),np.mean(avgdistlist),np.max(avgdistlist))
    print("filter state",meansep*surfrad,len(pointlist),len(poslist),accept_counter,reject_counter,full_counter)

    if len(avgdistlist) <= 0:
        print("ran out of holes -- avg dist based",meansep,avgdistlist)
        rndpt = np.random.choice(np.arange(len(pointlist)),1)[0]
        print(rndpt,pointlist[0])
        pi = pointlist[rndpt][0]
        pj = pointlist[rndpt][1]
        meanpos = 0.5*(surfpoints[pi,:]+surfpoints[pj,:])
        print(meanpos)
        
        distlist = [np.linalg.norm(meanpos - surfpoints[pi,:]),np.linalg.norm(meanpost - surfpoints[pj,:])]
        avgdistlist.append(np.mean(np.array(distlist)))

        pitris = np.where(autri == pi)
        pjtris = np.where(autri == pj)


        ### get the set of unique points from within the two sets of triangles -- will have at least one overlap each (pi in pjtris and pj in pitris)
        piuniq = np.unique(autri[pitris[0]])
        pjuniq = np.unique(autri[pjtris[0]])

        sharedpts = []
        fullpts = list(piuniq) + list(pjuniq)

        fullpts_trim = np.unique(np.array(fullpts))
    #     print(fullpts,fullpts_trim)

        ## shared points will store the set of points that participate in triangles containing pi and triangles containing pj -- 3rd vertices in the triangle including pi,pj
        for i in fullpts_trim:
            if np.isin(i,piuniq) or np.isin(i,pjuniq):
                sharedpts.append(i)

        trialtris = list(itertools.combinations(sharedpts,3))
    
        
        spamarea = 0
        counttris = 0
        for trialtri in trialtris:
            trtr = np.array(trialtri)
            if np.isin(pi,trtr) and np.isin(pj,trtr):
                ### sum their two areas
                spamarea += findArea(surfpoints[trtr],surfrad)
                counttris +=1
        
#         print(counttris,spamarea)
        if counttris < 2:
            print("WHOA WHOA!!!!")
            print(pi,pj)
            print(sharedpts)
            print(trialtris)
            exit(0)
        
        return meanpos,tuple(sharepts),spamarea,np.mean(np.array(distlist))

    max1 = np.argmax(np.array(avgdistlist))
    max2 = np.argmax(np.array(arealist))
    print("max dist index, area index",max1,max2)
    print("corr posns",poslist[max1],poslist[max2])
    print("dists",avgdistlist[max1],avgdistlist[max2])
    print("areas (returning basis)",arealist[max1],arealist[max2])
##    if avgdistlist[max2] <= 0.25*surfrad:
##        print("biggest hole not big enough")
##        return np.array([-nsurf,-nsurf,-nsurf])

    print("returning ",poslist[max2])


    
    return poslist[max2],sharedlist[max2],arealist[max2],meansep




def get_points_tris_maxhole(surfpoints,neighscale,nsurf,surfrad,npts):
    out_surfpoints = np.zeros([npts,3],np.float64)
    out_vertices = np.zeros([npts,3],np.float64)
    out_areas = np.zeros([npts],np.float64)

    inp_surfpoints = np.copy(surfpoints)
    meansep = 1.5
    for i in range(npts):
        ## get new point
        newpos,main_neighs,area,meansep = maxHole(inp_surfpoints,neighscale,nsurf,surfrad,nsurf+i,meansep)

        print("newpos",newpos,main_neighs,area,meansep)

        inp_surfpoints = np.vstack((inp_surfpoints,np.array([newpos[0],newpos[1],newpos[2]])))

        ### its 3-4 nearest neighbours -- get largest triangle to get,
        out_surfpoints[i,:] = np.array([newpos[0],newpos[1],newpos[2]])
        print(out_surfpoints[:i+1,:])

        vertices,area = tris_from_pts(main_neighs,surfpoints,surfrad)
        print("triangle",vertices,area)

        out_vertices[i,:] = np.copy(vertices)
        out_areas[i] = area

    return out_surfpoints,out_vertices,out_areas

        
