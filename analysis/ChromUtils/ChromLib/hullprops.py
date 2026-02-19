import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import plotly.graph_objects as go
import plotly.express as px


from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl

def get_colors(cmap_name, n=256):
    cmap = plt.get_cmap(cmap_name)
    colors = cmap(np.linspace(0, 1, n))
    
#     print(colors)
    spam = []
    spam.extend(colors)
#     print(spam)
    np.random.shuffle(spam)
#     print("hey",spam)
    spam[0]=[0,0,0,1]
    
    new_cmap = LinearSegmentedColormap.from_list('new_cmap', spam)
    return new_cmap,spam


def computeHullMembrane(config,realsurf:int=2500):
    """ this subroutine computes the convex hull of the membrane
    and can also compute the hulls of all chromosomes in it
    the entire config is supplied, including membrane geometry etc.
    """

    surfx = config.surfx[:realsurf]
    surfy = config.surfy[:realsurf]
    surfz = config.surfz[:realsurf]

    ella = config.ella
    ellb = config.ellb
    ellc = config.ellc
    surfrad = config.surfrad
    
    surfxp=surfx.reshape((1,surfx.shape[0]))
    surfyp=surfy.reshape((1,surfy.shape[0]))
    surfzp=surfz.reshape((1,surfz.shape[0]))
    surfpoints=np.hstack((surfxp.T,surfyp.T,surfzp.T))

    ### add a shift to get a different, smaller mesh
    spamnorm = np.sqrt(ella**2+ellb**2+ellc**2)
    shift = [ella*surfrad,ellb*surfrad,ellc*surfrad]/spamnorm
##    print(shift,surfrad)
    shiftsurfpoints = surfpoints[0,:]
    shiftsurfpoints = shiftsurfpoints.reshape((1,3))

    for surfpt in surfpoints:
        newpoint = surfpt - np.sign(surfpt)*np.array(shift)
        shiftsurfpoints = np.vstack((shiftsurfpoints,newpoint))
        
        
##    print(surfpoints.shape,shiftsurfpoints.shape)
    surfhull=ConvexHull(shiftsurfpoints)


##    return surfhull,shiftsurfpoints
    
    tri =[]
    j=0 ### reset index because previous loop
    for s in surfhull.simplices:
        tri.extend([(s[0],s[1],s[2])])

    utri = list(set(tri))  #utri is the list of unique triangles
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
    fig=go.Figure(go.Scatter3d(x=x, y=y, z=z, mode="markers", marker_size=0.05,marker_color='red',opacity=0.05))
##    fig.add_scatter3d(x=x, y=y, z=z, mode="markers", marker_size=10,marker_color='red',opacity=0.15)
    fig.add_mesh3d(x=x, y=y, z=z, i=I, j=J, k=K, color='red',opacity=0.05) #'rgba(30, 218, 245, 0.35)')
    fig.add_trace(lines)
    fig.update_layout(width=800, height=800, showlegend=False)
    fig.update_layout(
        scene = dict(
            xaxis = dict(nticks=4, range=[-1.5*ella,1.5*ella],),
            yaxis = dict(nticks=4, range=[-1.5*ella,1.5*ella],),
            zaxis = dict(nticks=4, range=[-1.0*ellb,1.0*ellb],),),
        width=700,
        margin=dict(r=20, l=10, b=10, t=10))

    return fig,surfhull,shiftsurfpoints






def visuHull(hull,inputpoints,
             colorlist: chr = ['rgba(150, 00, 55, 0.95)',
                                'rgba(105, 165, 0, 0.75)',
                                'rgba(5, 100, 150, 0.85)']):
    """ visualises the hull for a collection of points
    requires supplying the hull and the points themselves (a Nx5 array which also includes type info for coloring
    the colorlist should include a color for each type value
    """

##    fig=go.Figure(go.Scatter3d(x=[np.mean(inputpoints,axis=0)[0]],
##                               y=[np.mean(inputpoints,axis=0)[1]],
##                               z=[np.mean(inputpoints,axis=0)[2]], mode="markers", marker_size=5))
    
    x, y, z = hull.points.T
    verts = hull.points

##    print("IP",inputpoints.shape[0],inputpoints)
    xc = inputpoints[hull.vertices]
##    print(xc)

    ### identify simplices and faces thereof
    tri =[]
    j=0 ### reset index because previous loop
    for s in hull.simplices:
        tri.extend([(s[0],s[1],s[2])])

        
    utri = list(set(tri))  #utri is the list of unique triangles
    I, J, K = np.asarray(utri).T

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
                    line=dict(color= 'rgba(50,50,50,0.65)', width=2.5))

    
    mindc = np.min(inputpoints[inputpoints[:,4] > 0,5])
    maxdc = np.max(inputpoints[:,5])

        
    scatter_df = pd.DataFrame({'xc' : inputpoints[:,0],
                               'yc' : inputpoints[:,1],
                               'zc' : inputpoints[:,2],
                               'rc' : 50*inputpoints[:,3],
                               'tc' : inputpoints[:,4],
                               'decompaction' : inputpoints[:,5]})

    print(scatter_df)
    
    fig = px.scatter_3d(scatter_df,x='xc',y='yc',z='zc',size='rc',size_max=150,
                        color='decompaction',
                        color_continuous_scale='inferno',
                        range_color=[0,3],
                        opacity=0.7)

    fig.add_mesh3d(x=x, y=y, z=z, i=I, j=J, k=K,color='rgba(5, 100, 150, 0.35)',opacity=0.35)
    fig.add_trace(lines)
    fig.update_layout(width=600, height=600, showlegend=False)
    extentlow = abs(np.mean(x) - np.min(x))
    extenthigh = abs(np.max(x) - np.mean(x))
    fig.update_layout(
    scene = dict(
        xaxis = dict(nticks=4, range=[np.mean(x)-1.5*extentlow,np.mean(x)+1.5*extenthigh],),
        yaxis = dict(nticks=4, range=[np.mean(y)-1.5*extentlow,np.mean(y)+1.5*extenthigh],),
        zaxis = dict(nticks=4, range=[np.mean(z)-1.5*extentlow,np.mean(z)+1.5*extenthigh],),),
    width=700,
    margin=dict(r=20, l=10, b=10, t=10))
    fig.show()


    return fig



def computeHullSingleChrom(chrom, toplot: bool = False):
    """ this subroutine is supplied a chromosome and computes the hull for it
    point augmentation is a step in this process, chrom class has basic information
    chromHull class (an instance of which is created here) has augmented points,
    hull (and associated attributes) and geometric properties that are calculated here
    """
    chrid = chrom.chrid
    npatch = chrom.npatch
    chrrad = chrom.chrrad
    chrcentx = chrom.chrcentx
    chrcenty = chrom.chrcenty
    chrcentz = chrom.chrcentz
    patchrad = chrom.patchrad
    patchx = chrom.patchx
    patchy = chrom.patchy
    patchz = chrom.patchz
    patchlist = chrom.patchlist
    

    
    ptcollect = np.zeros([7*npatch+7,6],np.float64)
    chromcollect = np.zeros([npatch+1,6],np.float64)
##    ptcollecty = np.zeros([2*7*npatch+7],np.float64)
##    ptcollectz = np.zeros([2*7*npatch+7],np.float64)
    counter=0
    tp = 0
##    print("chrom and npatches",chrid,patchlist[0])



    for j in range(patchlist[0]+1):
        spampos = np.zeros([7,3],np.float64)
        if j==0:
            tp = 0
            radspam = chrrad
            currposx = chrcentx
            currposy = chrcenty
            currposz = chrcentz
            dc = 0
            
        else:
##            print(j,patchrad[chrid][j],patchx[chrid][j],patchy[chrid][j],patchz[chrid][j])
            radspam = patchrad[j]
            currposx = patchx[j]
            currposy = patchy[j]
            currposz = patchz[j]
            if j < patchlist[0]:
                tp = 1
            else:
                tp = 2

            dx = patchx[j] - chrcentx
            dy = patchy[j] - chrcenty
            dz = patchz[j] - chrcentz
            patchdist = np.sqrt(dx*dx + dy*dy + dz*dz)

            dc = patchdist/chrrad
            
        chromcollect[j,:] = np.array([currposx,currposy,currposz,radspam,tp,dc])
        
        
        spampos[0,:] = np.array([0,0,0])
        spampos[1,:] = np.array([-radspam,0,0])
        spampos[2,:] = np.array([radspam,0,0])
        spampos[3,:] = np.array([0,-radspam,0])
        spampos[4,:] = np.array([0,radspam,0])
        spampos[5,:] = np.array([0,0,-radspam])
        spampos[6,:] = np.array([0,0,radspam])
        thischromcounter = chrid*(npatch*7+7) + j

        for i in range(7):
            ptcollect[counter+i,0] = currposx + spampos[i,0]
            ptcollect[counter+i,1] = currposy + spampos[i,1]
            ptcollect[counter+i,2] = currposz + spampos[i,2]
            if i == 0:
                ptcollect[counter+i,3] = radspam
            else:
                ptcollect[counter+i,3] = 0.0005*radspam
            ptcollect[counter+i,4] = tp
            ptcollect[counter+i,5] = dc
##            allcollect[thischromcounter+i,0] = currposx + spampos[i,0]
##            allcollect[thischromcounter+i,1] = currposy + spampos[i,1]
##            allcollect[thischromcounter+i,2] = currposz + spampos[i,2]

        counter+=7
##        print(counter)
##        thischromcounter+=7

        ## isolated a chromosome -- find and plot it's hull

    ### now find and plot hull
    hull = ConvexHull(ptcollect[:,:3])
    
    xc = ptcollect[hull.vertices]

    ### identify simplices and faces thereof
    tri =[]
    j=0 ### reset index because previous loop
    for s in hull.simplices:
        tri.extend([(s[0],s[1],s[2])])
        
    utri = list(set(tri))  #utri is the list of unique triangles
    I, J, K = np.asarray(utri).T

    
    x, y, z = hull.points.T
    verts = hull.points

    if toplot==True:
        visuHull(hull,ptcollect)

    return hull,ptcollect[:,:3]

def computeHullSet(config,chromset,colorset=["red",
##                                             "yellow",
##                                             "darkviolet",
                                             "cyan"]
                   ):
    " compute the hull for a subset of chromosomes"

##    fig=go.Figure(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers", marker_size=1))

    #### get surfhull fig and hull etc
    ella = config.ella
    ellb = config.ellb
    ellc = config.ellc
    
    fig,surfhull,shiftsurfpoints = computeHullMembrane(config)

    
    for ii in range(len(chromset)):
        
        i = chromset[ii]
        chrid = config.chromarr[i].chrid
        npatch = config.chromarr[i].npatch
        chrrad = config.chromarr[i].chrrad
        chrcentx = config.chromarr[i].chrcentx
        chrcenty = config.chromarr[i].chrcenty
        chrcentz = config.chromarr[i].chrcentz
        patchrad = config.chromarr[i].patchrad
        patchx = config.chromarr[i].patchx
        patchy = config.chromarr[i].patchy
        patchz = config.chromarr[i].patchz
        patchlist = config.chromarr[i].patchlist

        ptcollect = np.zeros([7*npatch+7,6],np.float64)
        chromcollect = np.zeros([npatch+1,6],np.float64)
    ##    ptcollecty = np.zeros([2*7*npatch+7],np.float64)
    ##    ptcollectz = np.zeros([2*7*npatch+7],np.float64)
        counter=0
        tp = 0
    ##    print("chrom and npatches",chrid,patchlist[0])



        for j in range(patchlist[0]+1):
            spampos = np.zeros([7,3],np.float64)
            if j==0:
                tp = 0
                radspam = chrrad
                currposx = chrcentx
                currposy = chrcenty
                currposz = chrcentz
                dc = 0
                
            else:
    ##            print(j,patchrad[chrid][j],patchx[chrid][j],patchy[chrid][j],patchz[chrid][j])
                radspam = patchrad[j]
                currposx = patchx[j]
                currposy = patchy[j]
                currposz = patchz[j]
                if j < patchlist[0]:
                    tp = 1
                else:
                    tp = 2

                dx = patchx[j] - chrcentx
                dy = patchy[j] - chrcenty
                dz = patchz[j] - chrcentz
                patchdist = np.sqrt(dx*dx + dy*dy + dz*dz)

                dc = patchdist/chrrad
                
            chromcollect[j,:] = np.array([currposx,currposy,currposz,radspam,tp,dc])
            
            
            spampos[0,:] = np.array([0,0,0])
            spampos[1,:] = np.array([-radspam,0,0])
            spampos[2,:] = np.array([radspam,0,0])
            spampos[3,:] = np.array([0,-radspam,0])
            spampos[4,:] = np.array([0,radspam,0])
            spampos[5,:] = np.array([0,0,-radspam])
            spampos[6,:] = np.array([0,0,radspam])
            thischromcounter = chrid*(npatch*7+7) + j

            for k in range(7):
                ptcollect[counter+k,0] = currposx + spampos[k,0]
                ptcollect[counter+k,1] = currposy + spampos[k,1]
                ptcollect[counter+k,2] = currposz + spampos[k,2]
                if k == 0:
                    ptcollect[counter+k,3] = radspam
                else:
                    ptcollect[counter+k,3] = 0.0005*radspam
                ptcollect[counter+k,4] = tp
                ptcollect[counter+k,5] = dc
##            thischromcounter+=7
            counter +=7

        hull = ConvexHull(ptcollect[:,:3])
        x, y, z = hull.points.T
        xc = ptcollect[hull.vertices]
        verts = hull.points



        ### identify simplices and faces thereof
        tri =[]
        j=0 ### reset index because previous loop
        for s in hull.simplices:
            tri.extend([(s[0],s[1],s[2])])

        utri = list(set(tri))  #utri is the list of unique triangles
        I, J, K = np.asarray(utri).T
        tri_vertices= verts[np.asarray(utri)]
        Xe = []
        Ye = []
        Ze = []
        for T in tri_vertices:
            Xe += [T[k%3][0] for k in range(4)]+[ None]
            Ye += [T[k%3][1] for k in range(4)]+[ None]
            Ze += [T[k%3][2] for k in range(4)]+[ None]

        ## hard code colour selection
        
##        chrtype = chrid %(int(0.5*nchrom))
##        print("checking type for meshcolor",chrid,chrtype)
        meshcolor=colorset[ii % len(colorset)]

        #define the trace consisting in all triangle edges
        lines= go.Scatter3d(
                    x=Xe,
                    y=Ye,
                    z=Ze,
                    mode='lines',
                    name='',
                    line=dict(color= 'rgba(50,50,50,0.65)', width=2.5))

        #### make dataframe and add patches

        scatter_df = pd.DataFrame({'xc' : ptcollect[:,0],
                               'yc' : ptcollect[:,1],
                               'zc' : ptcollect[:,2],
                               'rc' : 50*ptcollect[:,3],
                               'tc' : ptcollect[:,4],
                               'decompaction' : ptcollect[:,5]})

##        print(scatter_df)
        
##        fig.add_trace(px.scatter_3d(scatter_df,x='xc',y='yc',z='zc',size='rc',size_max=150,
##                            color='decompaction',
##                            color_continuous_scale='inferno',
##                            range_color=[0,6],
##                            opacity=0.7))
        fig.add_scatter3d(x=ptcollect[:,0],y=ptcollect[:,1],z=ptcollect[:,2],mode='markers',
                          marker=dict(size=ptcollect[:,3],color=ptcollect[:,5],colorscale='inferno'))

    
##        fig.add_scatter3d(x=x, y=y, z=z, mode="markers", marker_size=2,marker_color='rgba(00, 218, 245, 0.25)')
##        for j in range(npatch+1):
##            xc, yc, zc, rc, tc = chromcollect[j,:]
##            if tc ==0:
##                colorstring='rgba(150, 00, 55, 0.75)'
##            if tc ==1:
##                colorstring='rgba(255,165,0, 0.75)'
##            if tc ==2:
##                colorstring='rgba(5, 100, 20, 0.95)'
##            fig.add_scatter3d(x=[xc], y=[yc], z=[zc], mode="markers", marker_size=10*rc,marker_color=colorstring)
    ##    fig.add_scatter3d(x=xc, y=yc, z=zc, mode="markers", marker_size=50)
    ##    fig.add_mesh3d(x=x, y=y, z=z, i=I, j=J, k=K, color=meshcolor,opacity=0.35) #'rgba(30, 218, 245, 0.35)')
        fig.add_mesh3d(x=x, y=y, z=z, i=I, j=J, k=K,color=meshcolor,opacity=0.25)
        fig.add_trace(lines)
        fig.update_layout(width=600, height=600, showlegend=False)
        fig.update_layout(
        scene = dict(
            xaxis = dict(nticks=4, range=[-1.5*ella,1.5*ella],),
            yaxis = dict(nticks=4, range=[-1.5*ella,1.5*ella],),
            zaxis = dict(nticks=4, range=[-1.0*ellb,1.0*ellb],),),
        width=700,
        margin=dict(r=20, l=10, b=10, t=10))
##        fig.show()

    return fig


def computeHullAll(config):

    """ subroutine to calculate the hulls of all chromosomes and of the membrane
        use calls to individual chromosome hulls and to membrane hull
    """

    nchrom = config.nchrom

    hull_arr = np.empty(nchrom,dtype=object)
    pt_arr = np.array([],np.float64)

    for i in range(nchrom):
        hullstore,ptstore = convexHullSingleChrom(config.chromarr[i])
        hull_arr[i] = hullstore
        pt_arr = stackDim(ptstore,pt_arr)


    surffig,surfhull,surfpoints = computeHullMembrane(config)

    return hull_arr,pt_arr,surfhull,surfpoints











def calcPackingFrac(config,proberad,nsamples,outprefix,pfdir):
    """ compute the packing fraction for the assembly by perturbatively analysing the change in hull size for each chrom
    given a probe insertion.
    """


    surfx = config.surfx
    surfy = config.surfy
    surfz = config.surfz

    
    ella = config.ella
    ellb = config.ellb
    ellc = config.ellc

    surfrad = config.surfrad
    
    surfxp=surfx.reshape((1,surfx.shape[0]))
    surfyp=surfy.reshape((1,surfy.shape[0]))
    surfzp=surfz.reshape((1,surfz.shape[0]))
    surfpoints=np.hstack((surfxp.T,surfyp.T,surfzp.T))

    surffig,surfhull,shiftsurfpoints = computeHullMembrane(config)

    full_vol=4.0*np.pi*ella*ellb*ellc/3.0
    volnumefree=0
    volnumeocc=0
    voldenom=0
    outsidecount=0
    membranecount=0
    pffile = outprefix + pfdir +'packing_frac_46c_'+str(config.timestep)+'.dat'
    outfile = open(pffile,"w")
    for t in range(nsamples):
        ## generate a random particle in 3d -- add centre and make a mesh of points around it
        centx = np.random.uniform(-1.2*ella,1.2*ella)
        centy = np.random.uniform(-1.2*ellb,1.2*ellb)
        centz = np.random.uniform(-1.2*ellc,1.2*ellc)

        if centx**2/ella**2 + centy**2/ellb**2 + centz**2/ellc**2 >=1.0:
            ### outside our ellipsoid
            outsidecount+=1
            continue

        u_ = np.linspace(0, 2 * np.pi, 100)
        v_ = np.linspace(0, np.pi, 100)

        ### code block below to add points

        xset = np.random.normal(0,1,100)
        yset = np.random.normal(0,1,100)
        zset = np.random.normal(0,1,100)

        normset = np.sqrt(xset**2 + yset**2 + zset**2)

        
        x = centx + proberad*xset/normset
        y = centy + proberad*yset/normset
        z = centz + proberad*zset/normset

        xp = x.reshape((1,x.shape[0]))
        yp = y.reshape((1,y.shape[0]))
        zp = z.reshape((1,z.shape[0]))

        probepts = np.hstack((xp.T,yp.T,zp.T))

        ## check if probe points within hull of nucleus

        overlapmembrane=0
        
        newsurfpoints = np.append(surfpoints,probepts,axis=0)
        newsurfhull = ConvexHull(newsurfpoints)

        ### compare the hull vertices, before and after
        ### if new point does not alter the membrane hull then is interior
        if len(surfhull.vertices) == len(newsurfhull.vertices):
            ## overlapping with membrane
            overlapmembrane=1
            break

        if overlapmembrane==0:
            membranecount+=1
        ### if overlapmembrane=0 and not overlapping with a chromosome (checked later) then update numerator-free vol
        ### else update numerator-occupied
                
        ## check if within hull of any of the chromosomes
        
    ##    for pt in probepts:
        ### append probepoints to each chrome in turn
        ### check if hull is modified by inclusion ##
        ### -- if yes then point is exterior
        ### -- if no then point is interior
        ### if any of probepts in interior to any chrom then update numerator-occupied
        ### if none of the probepts is in the interior of any of the chrom then update numerator-free vol
        overlapchrom=0


        for chrid in range(config.nchrom):
            oldchromhull,ptcollect = computeHullSingleChrom(config.chromarr[chrid])

            newptcollect = np.append(ptcollect,probepts[0,:].reshape((1,3)),axis=0)
            newchromhull = ConvexHull(newptcollect)

            if len(oldchromhull.vertices) == len(newchromhull.vertices):
                ### inside a chrom
                overlapchrom=1
                break

        ## now check counts
        if overlapmembrane==0 and overlapchrom==0:
            ### inside free volume
            volnumefree +=1.0
        elif overlapmembrane==0 and overlapchrom==1:
            volnumeocc +=1.0
        elif overlapmembrane==1 and overlapchrom==0:
            volnumeocc +=1.0
        else:
            volnumeocc +=1.0
            
        ### remember to always update denominator
        voldenom +=1.0

        if t%int(0.01*nsamples)==0:
    ##        print("old vertex count",chrid,len(oldchromhull.vertices))
    ##        print("new vertex count",chrid,len(newchromhull.vertices))
            print(volnumefree/voldenom,volnumeocc/voldenom,membranecount,t-outsidecount,t)
            outfile.write("%d %d %d %4.4f %4.4f \n" %(t,t-outsidecount,membranecount,volnumefree/voldenom,volnumeocc/voldenom))

##    print(volnumefree/voldenom,volnumeocc/voldenom,outsidecount)
    outfile.close()


    return volnumefree/voldenom,volnumeocc/voldenom,outsidecount



    
