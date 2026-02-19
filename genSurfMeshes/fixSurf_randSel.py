import matplotlib.pyplot as plt
import numpy as np
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

            if dist < 6.0*surfrad:
                neighlist[i,0] +=1
##                print(neighlist[i,0])
                neighlist[i,int(neighlist[i,0])] = j

                neighlist[j,0] +=1
                neighlist[j,int(neighlist[j,0])] = i

##    print(neighlist)

    return neighlist
                          


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

        if dist < 1.15*rad:
##            print("non-simplex neighbour here",ni,surfpoints[ni],midpt,dist)
            flag = 1
            break
##    print(pi,pj,midpt,flag)
    return flag
        

def findcrowded(surfpoints,nsurf,rad,ella,ellb,ellc):

    distmat = np.zeros([nsurf,nsurf],np.float64)

    ### first pop out beads away from surface
    offell = np.zeros([nsurf],np.float64)
    for i in range(nsurf):
        offell[i] = (surfpoints[i,0]/ella)**2 + (surfpoints[i,1]/ellb)**2 + (surfpoints[i,2]/ellc)**2 - 1.0

    popouts = np.where(offell > 0.5*rad)[0]

    print("list of popouts",popouts,len(popouts))

    if len(popouts) > 0:
        spam = np.random.choice(np.arange(0,len(popouts)),1,replace=True)
        print("returning a pop out -- %d\n",spam,popouts[spam],len(popouts))
        return popouts[spam]
        

    for i in range(nsurf-1):
        for j in range(i+1,nsurf):
            fpt = surfpoints[i,:]
            spt = surfpoints[j,:]

    ##            print(fpt,spt)

            distmat[i,j] = np.linalg.norm(fpt-spt)

            distmat[j,i] = np.linalg.norm(fpt-spt)
            ### make sure we dont consider zero entries (diagonals)
    dist1d = np.ravel(distmat)

##    outmat = np.zeros([nsurf],np.float64)
##    outmat = np.array([(surfpoints[i,0]/ella)**2 + (surfpoints[i,1]/ellb)**2 + (surfpoints[i,2]/ellc)**2 for i in range(nsurf)])
##    print(np.max(outmat))
##    plt.hist(dist1d[dist1d>0], bins=np.arange(0,0.5,0.01), 
##         color='blue', alpha=0.5, density=True, label='area histo')
##
##    plt.show()
##    print(distmat)
##    print("\n\nmin non-zero entries",np.where(distmat == np.min(distmat[distmat>0])),np.min(distmat[distmat>0])/rad,np.min(distmat[distmat>0]))

##    mindist = np.min(distmat[distmat>0])
##    mini = np.where(distmat == np.min(distmat[distmat>0]))[0][0]
##    minj = np.where(distmat == np.min(distmat[distmat>0]))[1][0]

    mindist = distmat[(distmat>0) &(distmat < 1.4*rad)]
    print("\n\n points too close",mindist,len(mindist))
##    print(np.logical_and(distmat > 0, distmat < 1.6*rad))
##    print(np.where((distmat > 0) &(distmat < 1.4*rad)))
    minilist = np.where((distmat > 0) &(distmat < 1.4*rad))[0]
    minjlist = np.where((distmat > 0) &(distmat < 1.4*rad))[1]
    print(minilist.shape,minjlist.shape)

    fsrl = 1
    ssrl = 1
    iters = 0

    while max(fsrl,ssrl)< 4: ## do until we get a point crowded by at least 4 others

        if iters > 500:
            fsrl = 0
            ssrl = 0
            break
        spam = np.random.choice(np.arange(0,minilist.shape[0]),1,replace=True)
        mini = minilist[spam][0]
        minj = minjlist[spam][0]

        print(iters," -- ",mini,minj,distmat[mini,minj],distmat[minj,mini])
        
    ##    mini = np.where(distmat == np.min(distmat[distmat>0]))[0][0]
    ##    minj = np.where(distmat == np.min(distmat[distmat>0]))[1][0]
    ##    (mini,minj) = np.argmin(distmat)
        print(iters," -- indices of min participants",mini,minj)
        print("locations",surfpoints[mini],surfpoints[minj])

        frow = distmat[mini,:]
        srow = distmat[minj,:]

        fsubrow = frow[(frow > 0) & (frow < 1.4*rad)] #np.where(frow > 0)[0]
        ssubrow = srow[(srow > 0) & (srow < 1.4*rad)] #np.where(srow > 0)[0]
        print(iters," -- rows with other indices",fsubrow,ssubrow,len(fsubrow),len(ssubrow),rad)
        print(iters," -- row indices (compare triangles)",np.where((frow > 0) &(frow < 1.4*rad)),np.where((srow > 0) &(srow < 1.4*rad)))

    ##    fsrl = len(np.where(frow < 1.6*rad))
    ##    ssrl = len(np.where(srow < 1.6*rad))

        fsrl = len(frow[(frow >0) & (frow < 1.4*rad)])
        ssrl = len(srow[(srow >0) & (srow < 1.4*rad)])

        
        iters +=1

        
    print("final close neighbours ",fsrl,ssrl,mini,minj," -- took many iters ", iters)
        

    if fsrl <=2 and ssrl <= 2:
        return -1
##    return min(fsrl,ssrl)

    if fsrl >=ssrl: ### mini more crowded than minj
        return mini
    else:
        return minj
    
    

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


def maxHole(surfpoints,neighlist,nsurf,surfrad,toreplace):

    print("replacing!",toreplace,surfpoints[toreplace])


    surfhull=ConvexHull(surfpoints)
    
    tris =[]
    j=0 ### reset index because previous loop
    for s in surfhull.simplices:
        tris.extend([(s[0],s[1],s[2])])

    utri = list(set(tris))  #utri is the list of unique triangles


    #### checking triangles that have toreplace
    autri = np.asarray(utri)
##    print("participating tris",np.where(autri == toreplace))
    print("old triangles")
    print(autri[np.where(autri == toreplace)[0]])
    

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

            

        ### get the participants of the largest edge, provided it is > 1.8*rad
        ### can definitely squeeze in one particle between these two --
        ### provided that the resultant transverse crowding is not too high!
        if np.max(spamdistarr) > 1.8*surfrad:
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
    
    poslist = []
    avgdistlist = []
    arealist = []
    for pts in pointlist:
        pi = pts[0]
        pj = pts[1]

        ### find triangles (indices thereof) that include these two points
        pitris = np.where(autri == pi)
        pjtris = np.where(autri == pj)
    #     print(pitris,pjtris)
        

        # print(np.unique(autri[pitris[0]]))

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
        if(len(sharedpts) <4):
            print("triangle subset i",autri[pitris[0]])
            print("triangle subset j",autri[pjtris[0]])
            print("points shared by two triangles -- should number 4 (incl pi,pj)",sharedpts,len(sharedpts),pi,pj)
        
        ### for this, get the mean pos, dist of mean from each.
        ### check if dist of mean from each > 0.8 sigma for all
        
    #     print(surfpoints[sharedpts])
    #     print(np.mean(surfpoints[sharedpts],axis=0))
    #     print([np.linalg.norm(x) for x in surfpoints[sharedpts] - np.mean(surfpoints[sharedpts],axis=0)])
    #     print([np.linalg.norm(x) for x in surfpoints[sharedpts] - midpt])
        
        meanpos = np.mean(surfpoints[sharedpts],axis=0)
        distlist = [np.linalg.norm(x) for x in surfpoints[sharedpts] - np.mean(surfpoints[sharedpts],axis=0)]
        
        if np.all(np.array(distlist)>1.2*surfrad):

            ### do a sanity check to ensure this pos not anywhere near our particles before appending
            testpos = meanpos.reshape(1,3)
            normarr = np.linalg.norm(surfpoints-testpos,axis=1)
            print("checking",normarr.shape,testpos,surfpoints[0],surfpoints[0]-testpos,normarr[0])
##            print(normarr)
            if np.any(normarr<surfrad):
##                print(testpos,meanpos)
##                print(surfpoints[0]-testpos,surfpoints[0]-meanpos)
##                
                print(np.where(normarr<surfrad),testpos,surfpoints[np.where(normarr<surfrad)])
                continue
            poslist.append(meanpos)
            avgdistlist.append(np.mean(np.array(distlist)))
            
            trialtris = list(itertools.combinations(sharedpts,3))
            
    #         print("as list",trialtris)
    #         print("as arr",np.array([trialtris]))
    #         print(pi,pj)
            
            spamarea = 0
            counttris = 0
            for trialtri in trialtris:
                trtr = np.array(trialtri)
                if np.isin(pi,trtr) and np.isin(pj,trtr):
##                    print("yaay",trtr)
##                    print("pts yay",surfpoints[trtr])
##                    print("area yaay",findArea(surfpoints[trtr],surfrad))
                    
                    ### sum there two areas
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
                
##   print(len(arealist),meanpos,np.mean(np.array(distlist)),spamarea)
##        else:
##            print("too close for ",pi,pj,sharedpts,distlist,1.2*surfrad)
##    print(len(poslist))
    print("final candidate hole number",len(avgdistlist))
##    print(len(arealist))

    # print(avgdistlist)

    if len(avgdistlist) <= 0:
        print("ran out of holes -- avg dist based")
        return np.array([-nsurf,-nsurf,-nsurf])

    max1 = np.argmax(np.array(avgdistlist))
    max2 = np.argmax(np.array(arealist))
    print("max dist index, area index",max1,max2)
    print("corr posns",poslist[max1],poslist[max2])
    print("dists",avgdistlist[max1],avgdistlist[max2])
    print("areas (returning basis)",arealist[max1],arealist[max2])
    if avgdistlist[max2] <= 1.0*surfrad:
        print("biggest hole not big enough")
        return np.array([-nsurf,-nsurf,-nsurf])

    print("returning ",poslist[max2]) 
    return poslist[max2]

    

    


    


conffile = sys.argv[1]
nsurf = int(sys.argv[2])
infile = open(conffile,"r")
lines = infile.readlines()
ella = float(sys.argv[4])
ellb = float(sys.argv[5])
ellc = float(sys.argv[6])

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




neighlist = genneighlist(surfpoints,nsurf,surfrad)

surfhull=ConvexHull(surfpoints)
    
tri =[]
j=0 ### reset index because previous loop
for s in surfhull.simplices:
    tri.extend([(s[0],s[1],s[2])])

utri = list(set(tri))  #utri is the list of unique triangles

##print(utri[0])
I, J, K = np.asarray(utri).T

x, y, z = surfhull.points.T
verts = surfhull.points
tri_vertices= verts[np.asarray(utri)]

Xe = []
Ye = []
Ze = []
for T in tri_vertices:
    Xe += [T[k%3][0] for k in range(4)]+[ None]
    Ye += [T[k%3][1] for k in range(4)]+[ None]
    Ze += [T[k%3][2] for k in range(4)]+[ None]



lines= go.Scatter3d(
                x=Xe,
                y=Ye,
                z=Ze,
                mode='lines',
                name='',
                line=dict(color= 'rgba(50,50,50,0.5)', width=1.5),opacity=0.4) 
fig=go.Figure(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers", marker_size=15,marker_color='blue',opacity=0.05))
##fig=go.Figure(go.Scatter3d(x=x, y=y, z=z, mode="markers", marker_size=15,marker_color='red',opacity=0.05))

##fig.add_scatter3d(x=x, y=y, z=z, mode="markers", marker_size=10,marker_color='red',opacity=0.15)
fig.add_mesh3d(x=x, y=y, z=z, i=I, j=J, k=K, color='red',opacity=0.15) #'rgba(30, 218, 245, 0.35)')
fig.add_trace(lines)
fig.update_layout(width=800, height=800, showlegend=False)
fig.update_layout(
    scene = dict(
        xaxis = dict(nticks=4, range=[-1.5*ella,1.5*ella],),
        yaxis = dict(nticks=4, range=[-1.5*ella,1.5*ella],),
        zaxis = dict(nticks=4, range=[-1.0*ellb,1.0*ellb],),),
    width=700,
    margin=dict(r=20, l=10, b=10, t=10))

#fig.show()




##print(findcrowded(surfpoints,nsurf,surfrad))

###### first find the worst crowder
##toreplace = findcrowded(surfpoints,nsurf,surfrad)
##
##### now find the biggest hole
##
##
##newpos = maxHole(surfpoints,nsurf,surfrad)
##print(newpos,toreplace,surfpoints[toreplace])
##surfpoints[toreplace] = newpos


##### looped segment
##nCrowd = 1
changelist=[]
candpositions = []
while len(changelist) < 40: ## as long as crowders exist
    #### first find the worst crowder
    toreplace = findcrowded(surfpoints,nsurf,surfrad,ella,ellb,ellc)
    nCrowd = toreplace
    
    

    ### now find the biggest hole

##    surfhull=ConvexHull(surfpoints)
##        
##    tri =[]
##    j=0 ### reset index because previous loop
##    for s in surfhull.simplices:
##        tri.extend([(s[0],s[1],s[2])])
##
##    utri = list(set(tri))  #utri is the list of unique triangles
##    I, J, K = np.asarray(utri).T
##
##    x, y, z = surfhull.points.T
##    verts = surfhull.points
##
##    tri_vertices= verts[np.asarray(utri)]
    
    newpos = maxHole(surfpoints,neighlist,nsurf,surfrad,toreplace)

    
    
    if toreplace == -1 or np.sum(newpos)==-3*nsurf:
        break

    print("old ",toreplace,surfpoints[toreplace])
    surfpoints[toreplace] = newpos
    print("new ",toreplace,surfpoints[toreplace])
    changelist.append(toreplace)
    candpositions.append(newpos)
    print("changed so far... ",changelist)
    print("new positions placed at...",candpositions)

    print("ellipsoidicity",[(x[0]/ella)**2 + (x[1]/ellb)**2 + (x[2]/ellc)**2 -1 for x in candpositions])

    print("old neigh",toreplace,neighlist[toreplace,:neighlist[toreplace,0]+1])

    neighlist = genneighlist(surfpoints,nsurf,surfrad)

    print("new neigh",toreplace,neighlist[toreplace,:neighlist[toreplace,0]+1])

    print("changes made!!!!",len(changelist),"\n\n")

##    for cand in changelist:
##        print("new neighs",cand,neighlist[cand,:neighlist[cand,0]+1])
    

#### find triangle areas
##print(utri)
##print(tri_vertices)

##areas = np.zeros([len(utri)],np.float64)
##areas = np.array([findArea(x,surfrad) for x in tri_vertices[:100]])
##print(areas.shape)
##plt.hist(areas, bins=np.arange(0,0.5,0.01), 
##         color='blue', alpha=0.5, density=True, label='area histo')
##
##plt.show()
print("changed these many",len(changelist))
print("changed indices",changelist)
print("positions",candpositions)

### make new hull for modified surf points


surfhull=ConvexHull(surfpoints)
    
tri =[]
j=0 ### reset index because previous loop
for s in surfhull.simplices:
    tri.extend([(s[0],s[1],s[2])])

utri = list(set(tri))  #utri is the list of unique triangles

##print(utri[0])
I, J, K = np.asarray(utri).T

x, y, z = surfhull.points.T
verts = surfhull.points
tri_vertices= verts[np.asarray(utri)]

Xe = []
Ye = []
Ze = []
for T in tri_vertices:
    Xe += [T[k%3][0] for k in range(4)]+[ None]
    Ye += [T[k%3][1] for k in range(4)]+[ None]
    Ze += [T[k%3][2] for k in range(4)]+[ None]



lines= go.Scatter3d(
                x=Xe,
                y=Ye,
                z=Ze,
                mode='lines',
                name='',
                line=dict(color= 'rgba(50,50,50,0.5)', width=1.5),opacity=0.4) 
fig=go.Figure(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers", marker_size=15,marker_color='blue',opacity=0.05))
##fig=go.Figure(go.Scatter3d(x=x, y=y, z=z, mode="markers", marker_size=15,marker_color='red',opacity=0.05))

##fig.add_scatter3d(x=x, y=y, z=z, mode="markers", marker_size=10,marker_color='red',opacity=0.15)
fig.add_mesh3d(x=x, y=y, z=z, i=I, j=J, k=K, color='red',opacity=0.15) #'rgba(30, 218, 245, 0.35)')
fig.add_trace(lines)
fig.update_layout(width=800, height=800, showlegend=False)
fig.update_layout(
    scene = dict(
        xaxis = dict(nticks=4, range=[-1.5*ella,1.5*ella],),
        yaxis = dict(nticks=4, range=[-1.5*ella,1.5*ella],),
        zaxis = dict(nticks=4, range=[-1.0*ellb,1.0*ellb],),),
    width=700,
    margin=dict(r=20, l=10, b=10, t=10))

#fig.show()

outfile = open(sys.argv[3],"w")

for i in range(nsurf):
    outfile.write("%d %f %f %f %f\n" %(i,surfrad,surfpoints[i,0],surfpoints[i,1],surfpoints[i,2]))
outfile.close()


name = sys.argv[3][:-4] + '_old_ovito.dat'

outfile = open(name,"w")
outfile.write("%d\n" %(2500))
outfile.write("spam\n")
for i in range(nsurf):
    if np.isin(i,np.array(changelist)):
        outfile.write("%d %d %f %f %f %f %f\n" %(i,2,0.2,surfrad,oldsurfpoints[i,0],oldsurfpoints[i,1],oldsurfpoints[i,2]))
    else:
        outfile.write("%d %d %f %f %f %f %f\n" %(i,1,0.5,surfrad,oldsurfpoints[i,0],oldsurfpoints[i,1],oldsurfpoints[i,2]))
outfile.close()



name = sys.argv[3][:-4] + '_ovito.dat'

outfile = open(name,"w")
outfile.write("%d\n" %(2500))
outfile.write("spam\n")
for i in range(nsurf):
    if np.isin(i,np.array(changelist)):
        outfile.write("%d %d %f %f %f %f %f\n" %(i,2,0.2,surfrad,surfpoints[i,0],surfpoints[i,1],surfpoints[i,2]))
    else:
        outfile.write("%d %d %f %f %f %f %f\n" %(i,1,0.5,surfrad,surfpoints[i,0],surfpoints[i,1],surfpoints[i,2]))
outfile.close()

