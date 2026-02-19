import numpy as np
import math
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import plotly.express as px

from scipy import ndimage
from itertools import combinations

import ChromUtils.utils
from ChromUtils.Structures import *
from ChromUtils.ChromLib.hullprops import *





def connectedGraph(config,intmat,arglist):

    nchrom = config.nchrom
    nato = config.nchrom+config.nchrom*config.npatch

    dist_thresh = float(arglist[2])
    baseline_ella = float(arglist[3])
    baseline_ellb = float(arglist[4])
    baseline_ellc = float(arglist[5])

    G = nx.empty_graph(n=nato)
##    attrs = {0: {"attr1": 20, "attr2": "nothing"}, 1: {"attr2": 3}}
##    nx_latt.set_node_attributes(G, attrs)

    parcount=0
    for i in range(config.nchrom):
        patchlist = config.chromarr[i].patchlist
        
        

        for j in range(patchlist[0] + 1):
            if j==0:
                attr = {parcount: {"chromind": i,
                                   "patchind": j,
                                   "radius": config.chromarr[i].chrrad,
                                   "pos":(config.chromarr[i].chrcentx,
                                          config.chromarr[i].chrcenty,
                                          config.chromarr[i].chrcentz),
                                   "type": 0 ,
                                   "isBound": 0, 
                                   "isIntMing": 0,
                                   }
                        }
                nx.set_node_attributes(G,attr)
                                   
##                typearr[parcount] = 0
##                radarr[parcount] = config.chromarr[i].chrrad

            elif j > 0 and j < config.chromarr[i].patchlist[0]:

                attr = {parcount: {"chromind": i,
                                   "patchind": j,
                                   "radius": config.chromarr[i].patchrad[j],
                                   "pos":(config.chromarr[i].patchx[j],
                                          config.chromarr[i].patchy[j],
                                          config.chromarr[i].patchz[j]),
                                   "type": 1 ,
                                   "isBound": 0, 
                                   "isIntMing": 0,
                                   }
                        }
                nx.set_node_attributes(G,attr)
                
##                typearr[parcount] = 1
##                radarr[parcount] = config.chromarr[i].patchrad[j]

            else:
                attr = {parcount: {"chromind": i,
                                   "patchind": j,
                                   "radius": config.chromarr[i].patchrad[j],
                                   "pos":(config.chromarr[i].patchx[j],
                                          config.chromarr[i].patchy[j],
                                          config.chromarr[i].patchz[j]),
                                   "type": 2 ,
                                   "isBound": 0, 
                                   "isIntMing": 0,
                                   }
                        }
                nx.set_node_attributes(G,attr)
##                typearr[parcount] = 2
##                radarr[parcount] = config.chromarr[i].patchrad[j]
            parcount+=1


    ### selecting subsets of nodes
##    selected_nodes = [n for n,v in G.nodes(data=True) if v['since'] == 'December 2008']

    #### add edges between each chrom centre and its patches
##    [x for x,y in g.nodes(data=True) if y['d']==1 and y['a'] == 9]
    chromnodes = [x for x,y in G.nodes(data=True) if y['type'] == 0]
##    chromnodes = list(range(config.nchrom))
##    print("chromnodes",chromnodes)
    ### all nodes have to have a dict with key 'type' for this to not give err
    for chrom in chromnodes:
        chid = G.nodes[chrom]['chromind']
        chromnode = G.nodes(data=True)[chrom]
##        print(chrom,chid,chromnode)

        patches = [x for x,y in G.nodes(data=True) if y['chromind']== chid and y['type'] == 1]
##        print(chid,patches)
        for patchid in patches:
            G.add_edge(chrom,patchid)
            



    ### now add edges between patches that are intemingling not on same chrom

    patchnodes = [x for x,y in G.nodes(data=True) if y['type'] == 1]

    patchcombs = combinations(patchnodes,2)

    for pairs in patchcombs:
        p1 = pairs[0]
        p2 = pairs[1]

        ## check that not on same chromosome
        if G.nodes[p1]['chromind'] != G.nodes[p2]['chromind']:

            ### get refsep
            refsep = G.nodes[p1]['radius'] + G.nodes[p2]['radius']
            

            ### check that distance between them is small
            px1 = G.nodes[p1]['pos'][0]
            py1 = G.nodes[p1]['pos'][1]
            pz1 = G.nodes[p1]['pos'][2]

            px2 = G.nodes[p2]['pos'][0]
            py2 = G.nodes[p2]['pos'][1]
            pz2 = G.nodes[p2]['pos'][2]

            dist = np.sqrt((px1 - px2)**2 + (py1 - py2)**2 + (pz1 - pz2)**2)

            

            if dist < dist_thresh*refsep:
##                print(p1,p2,G.nodes[p1]['chromind'],G.nodes[p2]['chromind'],dist,dist_thresh*refsep)
                G.add_edge(p1,p2)
                ### reset their intemingling states to 1
                G.nodes[p1]['isIntMing'] = 1
                G.nodes[p2]['isIntMing'] = 1

                ch1 = [x for x,y in G.nodes(data=True) if y['chromind'] == G.nodes[p1]['chromind'] and y['type'] == 0]
                ch2 = [x for x,y in G.nodes(data=True) if y['chromind'] == G.nodes[p2]['chromind'] and y['type'] == 0]

                print("chroms",ch1,ch2,"patches",p1,p2)
                ch1 = ch1[0]
                ch2 = ch2[0]

                G.nodes[ch1]['isIntMing'] = 1
                G.nodes[ch2]['isIntMing'] = 1
                ### if the corresponding intmat entry for the node indices is non-zero -- set intermingling to 2
                if intmat[p1,p2] > 0.0:
                    G.nodes[p1]['isIntMing'] = 2
                    G.nodes[p2]['isIntMing'] = 2

                    ### if 2 then set chrom int index to 2
                    ch1 = [x for x,y in G.nodes(data=True) if y['chromind'] == G.nodes[p1]['chromind'] and y['type'] == 0]
                    ch2 = [x for x,y in G.nodes(data=True) if y['chromind'] == G.nodes[p2]['chromind'] and y['type'] == 0]

                    print("chroms",ch1,ch2,"patches",p1,p2)
                    ch1 = ch1[0]
                    ch2 = ch2[0]
##                    print("chrom types",G.nodes[ch1]['chromind'],G.nodes[ch2]['chromind'])

                    G.nodes[ch1]['isIntMing'] = 2
                    G.nodes[ch2]['isIntMing'] = 2
    
    ### then add a isBound attribute change to the lamin patches based on their proximity to one or more surf points
    lamnodes = [x for x,y in G.nodes(data=True) if y['type'] == 2]
##    print(lamnodes)
    for n in lamnodes:
         ### get refsep
        refsep = G.nodes[n]['radius'] + config.surfrad
##        print(n,G.nodes[n]['type'])
        ### position lamin node
        px1 = G.nodes[n]['pos'][0]
        py1 = G.nodes[n]['pos'][1]
        pz1 = G.nodes[n]['pos'][2]
        done =0
        for i in range(config.nsurf):
            px2 = config.surfx[i]
            py2 = config.surfy[i]
            pz2 = config.surfz[2]

            dist = np.sqrt((px1 - px2)**2 + (py1 - py2)**2 + (pz1 - pz2)**2)

            if dist < dist_thresh*refsep:
                G.nodes[n]['isBound'] = 2
                G.nodes[n]['isIntMing'] = 2

                
                ch1 = [x for x,y in G.nodes(data=True) if y['chromind'] == G.nodes[n]['chromind'] and y['type'] == 0]

                ch1 = ch1[0]

                G.nodes[ch1]['isBound'] = 2
                G.nodes[ch1]['isIntMing'] = 2
                print("lamin edge",n,ch1)

                G.add_edge(n,ch1)
                done = 1
                break
        if done ==1:
            done =0



    ### plot the full graph of 2s


    # Get the maximum number of edges adjacent to a single node
    edge_max = max([G.degree(i) for i in range(nato)])
    print(G.edges())
    ## colour based on degree
##    colors = [plt.cm.plasma(G.degree(i)/edge_max) for i in range(nato)]

##    ## define colors based on node centrality
##
##    centlist = [ (G.nodes[i]['pos'][0]/baseline_ella)**2 +
##                 (G.nodes[i]['pos'][1]/baseline_ellb)**2 +
##                 (G.nodes[i]['pos'][2]/baseline_ellc)**2 - 1 for i in range(n)]
##
##    colors = [ plt.cm.plasma(abs(cent)) for cent in centlist]
    ## plotly 3D network plot
    pos = nx.get_node_attributes(G, 'pos')
##    print(pos)
    
    # Get number of nodes
    n = G.number_of_nodes()
    
    fig=go.Figure(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers", marker_size=1))
    # Loop on the pos dictionary to extract the x,y,z coordinates of each node
    for key, value in pos.items():
##        print(G.nodes[key],G.nodes[key]['pos'])
        if G.nodes[key]['isIntMing'] ==0 and G.nodes[key]['isBound']==0 and G.nodes[key]['type'] > 0:
            continue
        xi = value[0]
        yi = value[1]
        zi = value[2]

##        cent = (G.nodes[key]['pos'][0]/baseline_ella)**2 + (G.nodes[key]['pos'][1]/baseline_ellb)**2 + (G.nodes[key]['pos'][2]/baseline_ellc)**2 - 1
##        colindex = int(-9*cent) #int(10*(cent - min(cent))/(max(cent) - min(cent)))
##        print(px.colors.sequential.Plasma)
##        print(px.colors.named_colorscales(viridis))
##        print(px.colors.sequential.Rainbow)
##        print(cent,int(9*cent))
        cent = G.nodes[key]['isIntMing']
        centmax = 2
        colindex = int(9*cent/centmax)
        
##        print(colindex,len(px.colors.sequential.Rainbow))
        if colindex >= 9:
            colindex = 8
        if colindex < 0:
            colindex = 0
##        color = px.colors.sequential.Rainbow[colindex]
        color = px.colors.sequential.Viridis[colindex]
              
##        color = px.colors.sequential.Plasma(abs(cent))

##            print(key,xi,yi,zi)
        
        # Scatter plot
##            ax.scatter(xi, yi, zi, c=np.atleast_2d(np.array([colors[key]])), s=10+10*G.degree(key), edgecolors='k', alpha=0.7)
        fig.add_scatter3d(x=[xi], y=[yi], z=[zi], mode="markers", marker_size=25*G.nodes[key]['radius'],marker_color=color,opacity=0.75)
##        fig.add_scatter3d(x=[xi], y=[yi], z=[zi], mode="markers",
##                          marker_size=5*G.nodes[key]['radius'],
##                          marker_color='viridis',
##                          color_continuous_scale=px.colors.sequential.Rainbow,
##                          opacity=0.5)
        
                 
    
    # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
    # Those two points are the extrema of the line to be plotted
    for i,j in enumerate(G.edges()):
        print(i,j,pos[j[0]][0],pos[j[1]][0])
        x = np.array((pos[j[0]][0], pos[j[1]][0]))
        y = np.array((pos[j[0]][1], pos[j[1]][1]))
        z = np.array((pos[j[0]][2], pos[j[1]][2]))

##        if G.nodes[j[0]]['isIntMing'] == 0 and G.nodes[j[1]]['isIntMing'] == 0:
##            continue
        print(x,y,z)
    
    # Plot the connecting lines
##            ax.plot(x, y, z, c='black', alpha=0.5)
        if G.nodes[j[0]]['isIntMing'] + G.nodes[j[1]]['isIntMing'] >= 3:
            linecol = 'red'
            linewdt = 6
        elif G.nodes[j[0]]['isIntMing'] + G.nodes[j[1]]['isIntMing'] == 2:
            linecol = 'orange'
            linewdt = 4
        elif G.nodes[j[0]]['isIntMing'] + G.nodes[j[1]]['isIntMing'] < 2 and G.nodes[j[0]]['isIntMing'] + G.nodes[j[1]]['isIntMing'] > 0:
            linecol = 'blue'
            linewdt = 2
        elif G.nodes[j[0]]['isIntMing'] + G.nodes[j[1]]['isIntMing'] == 0:
            linecol = 'black'
            linewdt = 1
##        line = go.Scatter3d(x=x,
##                            y=y,
##                            z=z,
##                            mode='lines',
##                            name='',
##                            line=dict(color= 'rgba(50,50,50,0.65)', width=2.5))
        line = go.Scatter3d(x=x,
                            y=y,
                            z=z,
                            mode='lines',
                            name='',
                            line=dict(color= linecol, width=linewdt),
                            opacity=0.5)
        fig.add_trace(line)
    fig.update_layout(width=600, height=600, showlegend=False)
    fig.update_layout(
    scene = dict(
        xaxis = dict(nticks=4, range=[-1.2*baseline_ella,1.2*baseline_ella],),
        yaxis = dict(nticks=4, range=[-1.2*baseline_ella,1.2*baseline_ella],),
        zaxis = dict(nticks=4, range=[-1.0*baseline_ellb,1.0*baseline_ellb],),),
    width=700,
    margin=dict(r=20, l=10, b=10, t=10))
    fig.show()
##    # 3D network plot
##    with plt.style.context(('ggplot')):
##        
##        fig = plt.figure(figsize=(10,7))
##        ax = Axes3D(fig)
##        
##        # Loop on the pos dictionary to extract the x,y,z coordinates of each node
##        for key, value in pos.items():
##            xi = value[0]
##            yi = value[1]
##            zi = value[2]
##            
##            # Scatter plot
##            ax.scatter(xi, yi, zi, c=colors[key], s=20+20*G.nodes[key]['radius'], edgecolors='k', alpha=0.7)
##        
##        # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
##        # Those two points are the extrema of the line to be plotted
##        for i,j in enumerate(G.edges()):
##
##            x = np.array((pos[j[0]][0], pos[j[1]][0]))
##            y = np.array((pos[j[0]][1], pos[j[1]][1]))
##            z = np.array((pos[j[0]][2], pos[j[1]][2]))
##        
##        # Plot the connecting lines
##            ax.plot(x, y, z, c='black', alpha=0.5)
##    
##    # Set the initial view
##    ax.view_init(30, angle)
##
##    # Hide the axes
##    ax.set_axis_off()
##    save = False
##    if save is not False:
##        plt.savefig(outdir+ outfile+".png")
##        plt.close('all')
##    else:
##        plt.show()

    ### calculate shortest path between all pairs of isBounds passing through 1s and passing through 2s
    ### write largest distances to file (like in orientations/geometrics)

def intensityMaps(config,arglist,makefig=False):

    ### this subroutine adds a gaussian core potential of appropriate epsilon and lengthscale according to radius
    ### we then get a 3D lattice intensity map of the nucleus
    ### make a node for every lattice point
    ### connect edge if intensity below threshold
    ### mark edge nodes if within hull but close
    ### get largest path of low intensity nodes
    ### write to file
    ### save image of the intensity map in 3D



    ### to make intensity matrix
    outdir = arglist[1]
    lattice_size = float(arglist[2])
    baseline_ella = float(arglist[3])
    baseline_ellb = float(arglist[4])
    baseline_ellc = float(arglist[5])
    eps_blob = float(arglist[6])

    
    
    
    nx_latt = int(2*baseline_ella/lattice_size)
    ny_latt = int(2*baseline_ellb/lattice_size)
    nz_latt = int(2*baseline_ellc/lattice_size)


    intensity_lattice = np.zeros([nz_latt,ny_latt,nx_latt],np.float64)
    pos_lattice = np.zeros([nz_latt,ny_latt,nx_latt],np.float64)
    
    for i in range(config.nchrom):
        # print("chrom",i,nx_latt,ny_latt,nz_latt)

        patchlist = config.chromarr[i].patchlist

        for j in range(patchlist[0]+1):

            if j ==0:
                centx = config.chromarr[i].chrcentx
                centy = config.chromarr[i].chrcenty
                centz = config.chromarr[i].chrcentz
                rad = config.chromarr[i].chrrad
                alphij = 1.0/(1.04*rad)

                lattx = math.floor((centx + baseline_ella)/(2*baseline_ella))
                latty = math.floor((centy + baseline_ellb)/(2*baseline_ellb))
                lattz = math.floor((centz + baseline_ellc)/(2*baseline_ellc))

                ## add gaussian core centred at centx, centy, centz --
                for iz in range(nz_latt):
                    for iy in range(ny_latt):
                        for ix in range(nx_latt):
                            ptx = -baseline_ella + (ix+0.5)*lattice_size
                            pty = -baseline_ellb + (iy+0.5)*lattice_size
                            ptz = -baseline_ellc + (iz+0.5)*lattice_size
                
                            rijsq = (ptx - centx)**2 + (pty - centy)**2 + (ptz - centz)**2
                            if rijsq < 0.05:
                                print(i,eps_blob,rad,-alphij,rijsq,eps_blob*rad*np.exp(-alphij*rijsq))
                            intensity_lattice[iz,iy,ix] += eps_blob*rad*np.exp(-alphij*rijsq) ### eps_blob of chr cent is 4 times that of patches

            elif j > 0 and j < patchlist[0]:
                centx = config.chromarr[i].patchx[j]
                centy = config.chromarr[i].patchy[j]
                centz = config.chromarr[i].patchz[j]
                rad = config.chromarr[i].patchrad[j]
                alphij = 1.0/(1.04*rad)

                lattx = math.floor((centx + baseline_ella)/(2*baseline_ella))
                latty = math.floor((centy + baseline_ellb)/(2*baseline_ellb))
                lattz = math.floor((centz + baseline_ellc)/(2*baseline_ellc))

                ## add gaussian core centred at centx, centy, centz --
                for iz in range(nz_latt):
                    for iy in range(ny_latt):
                        for ix in range(nx_latt):
                            ptx = -baseline_ella + (ix+0.5)*lattice_size
                            pty = -baseline_ellb + (iy+0.5)*lattice_size
                            ptz = -baseline_ellc + (iz+0.5)*lattice_size
                
                            rijsq = (ptx - centx)**2 + (pty - centy)**2 + (ptz - centz)**2

                            intensity_lattice[iz,iy,ix] += eps_blob*rad*np.exp(-alphij*rijsq) ### eps_blob of patch is gaussian core peak height
                            
            elif j ==patchlist[0]:
                centx = config.chromarr[i].patchx[j]
                centy = config.chromarr[i].patchy[j]
                centz = config.chromarr[i].patchz[j]
                rad = config.chromarr[i].patchrad[j]
                alphij = 1.0/(1.04*rad)

                lattx = math.floor((centx + baseline_ella)/(2*baseline_ella))
                latty = math.floor((centy + baseline_ellb)/(2*baseline_ellb))
                lattz = math.floor((centz + baseline_ellc)/(2*baseline_ellc))

                ## add gaussian core centred at centx, centy, centz --
                for iz in range(nz_latt):
                    for iy in range(ny_latt):
                        for ix in range(nx_latt):
                            ptx = -baseline_ella + (ix+0.5)*lattice_size
                            pty = -baseline_ellb + (iy+0.5)*lattice_size
                            ptz = -baseline_ellc + (iz+0.5)*lattice_size
                
                            rijsq = (ptx - centx)**2 + (pty - centy)**2 + (ptz - centz)**2
                            

                            intensity_lattice[iz,iy,ix] += eps_blob*rad*np.exp(-alphij*rijsq) ### eps_blob of lamin patch is kept same as centre
                

                ## NOTE: next update for efficiency update full intensity_lattice based on distance from centx,centy,centz lattice points

    ### alternately can specify points where stuff is located and use scipy ndimage gaussian filter
##    pts = (l * np.random.rand(3, 15)).astype(np.int)
##    vol[tuple(indices for indices in pts)] = 1
##    vol = ndimage.gaussian_filter(vol, 4)
##    vol /= vol.max()
    #### now we have the full intensity matrix
    # print(intensity_lattice)

    # print(np.max(intensity_lattice.flatten()))

    ### do a Poisson sampling of the intensity matrix

    sampled_intensity_lattice = np.zeros_like(intensity_lattice)
    for iz in range(nz_latt):
        for iy in range(ny_latt):
            for ix in range(nx_latt):
                if intensity_lattice[iz,iy,ix] > 10:
                    print(ix,iy,iz,intensity_lattice[iz,iy,ix])
                sampled_intensity_lattice[iz,iy,ix] = np.random.poisson(lam=intensity_lattice[iz,iy,ix],size=1)

    

    #### we need to find a way to visualise it
    Z, Y, X = np.mgrid[:nz_latt, :ny_latt, :nx_latt]
    # print(sampled_intensity_lattice[sampled_intensity_lattice > 0])
    
##    vol = np.zeros((l, l, l))
    vol = sampled_intensity_lattice

    
##    
    if makefig:
        print("making figure")
        fig = go.Figure(data=go.Volume(
            x=X.flatten(), y=Y.flatten(), z=Z.flatten(),
            value=vol.flatten()/np.max(vol.flatten()),
            isomin=0.0,
            isomax=1.0,
            opacity=0.1,
            surface_count=25,
            ))
        fig.update_layout(scene_xaxis_showticklabels=False,
                        scene_yaxis_showticklabels=False,
                        scene_zaxis_showticklabels=False)
        fig.show()

    ### save the figure


    return intensity_lattice

def compileIntensity(intMat_series):

    nmat = intMat_series.shape[0]

    avgIntMat = np.sum(intMat_series,axis=0)

    return avgIntMat

def plotIntensity(intMat,arglist,outprefix):
    fig,axs=plt.subplots(nrows=1, ncols=1,sharex=True,sharey=True,figsize=(15,10),gridspec_kw={'width_ratios': [0.8]})


    nz_latt = intMat.shape[0]
    ny_latt = intMat.shape[1]
    nx_latt = intMat.shape[2]
    lattice_size = float(arglist[2])
    baseline_ella = float(arglist[3])
    baseline_ellb = float(arglist[4])
    baseline_ellc = float(arglist[5])

    centMat = np.zeros_like(intMat)

    for iz in range(nz_latt):
        for iy in range(ny_latt):
            for ix in range(nx_latt):
                ptx = -baseline_ella + (ix+0.5)*lattice_size
                pty = -baseline_ellb + (iy+0.5)*lattice_size
                ptz = -baseline_ellc + (iz+0.5)*lattice_size

                centMat[iz,iy,ix] = (ptx/baseline_ella)**2 + (pty/baseline_ellb)**2 + (ptz/baseline_ellc)**2 - 1

    print(np.min(centMat),np.max(centMat))
    flatMat = intMat[centMat<0.05].flatten()
    counts,bins = np.histogram(flatMat,bins=40)
    plt.hist(bins[:-1],bins,weights=counts,color='blue',label='rawIntensity',histtype='step',linewidth=4)

    plt.xlabel("Raw Intensity", fontsize=20)  
    plt.ylabel("Frequency", fontsize=20)
    plt.xticks(fontsize=18)  
    plt.yticks(fontsize=18)
    labels=["Intensity distribution"]
    plt.legend(labels,fontsize=20)
    outdir = outprefix + arglist[1]

    plotname=outdir + '/Intensity_distribution_latt'+str(arglist[2])+'.png'
    print(plotname,counts,bins)
##    plt.show()
    plt.savefig(plotname)

    filename = outdir + '/Intensity_distribution_latt'+str(arglist[2])+'.dat'
    spfile = open(filename,"w")
    for i in range(len(list(bins))-2):
        spfile.write("%4.4f %4.4f\n" %(0.5*(bins[i]+bins[i+1]),counts[i]))
    spfile.close()
##    plt.show()

    # Calculate the normalized distance from the ellipsoid boundary
    norm_centMat = np.abs(centMat)

    # Flatten the matrices for easier processing
    flat_intMat = intMat.flatten()
    flat_norm_centMat = norm_centMat.flatten()

    # Define bins for the normalized distance
    bins = np.linspace(0, np.max(flat_norm_centMat), 50)

    # Initialize an array to store intensity distributions
    intensity_distribution = []

    # Calculate the intensity distribution for each bin
    for i in range(len(bins) - 1):
        mask = (flat_norm_centMat >= bins[i]) & (flat_norm_centMat < bins[i + 1])
        intensities_in_bin = flat_intMat[mask]
        if len(intensities_in_bin) > 0:
            avg_intensity = np.mean(intensities_in_bin)
        else:
            avg_intensity = 0
        intensity_distribution.append(avg_intensity)

    # Write the distribution to a .dat file
    outdir = outprefix + arglist[1]
    filename = outdir + '/Intensity_vs_Distance_distribution.dat'
    with open(filename, "w") as spfile:
        for i in range(len(bins) - 1):
            mid_bin = 0.5 * (bins[i] + bins[i + 1])
            spfile.write(f"{mid_bin} {intensity_distribution[i]:.4f}\n")

    print(f"Intensity vs Distance distribution written to {filename}")

    # Generate a 2D maxz projection of intMat
    maxz_projection = np.max(intMat, axis=0)

    # Create a 2D plot
    plt.figure(figsize=(10, 8))
    plt.imshow(maxz_projection, cmap='viridis', origin='lower', extent=[-baseline_ella, baseline_ella, -baseline_ellb, baseline_ellb])
    plt.colorbar(label='Intensity')
    plt.title('MaxZ Projection of Intensity Matrix', fontsize=16)
    plt.xlabel('X-axis (units)', fontsize=14)
    plt.ylabel('Y-axis (units)', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # Save the plot
    outdir = outprefix + arglist[1]
    plotname = outdir + '/MaxZ_Projection_Intensity.png'
    plt.savefig(plotname)
    print(f"MaxZ projection plot saved to {plotname}")
    plt.close()

    # Write the maxz_projection matrix to a .dat file
    filename = outdir + '/MaxZ_Projection_Intensity.dat'
    with open(filename, "w") as spfile:
        for iy in range(maxz_projection.shape[0]):
            for ix in range(maxz_projection.shape[1]):
                x_coord = -baseline_ella + (ix + 0.5) * lattice_size
                y_coord = -baseline_ellb + (iy + 0.5) * lattice_size
                intensity_value = maxz_projection[iy, ix]
                spfile.write(f"{x_coord:.4f} {y_coord:.4f} {intensity_value:.4f}\n")
    print(f"MaxZ projection data written to {filename}")
    


def channelPercolation(config,arglist):

    
    outdir=arglist[1]
##    outfile=arglist[1]
    lattice_size = float(arglist[2])
    baseline_ella = float(arglist[3])
    baseline_ellb = float(arglist[4])
    baseline_ellc = float(arglist[5])
    eps_blob = float(arglist[6])
    hc_fact = float(arglist[7])

    intensity_lattice = intensityMaps(config,arglist)

    flat = intensity_lattice.flatten()

##    hc_thresh = np.percentile(flat,75) #,hc_ptile) ## using 75th percentile in debug
    hc_thresh = np.mean(flat) + hc_fact*np.std(flat)
    ec_thresh = np.mean(flat) - hc_fact*np.std(flat)

    print("hc ec intensity thresholds -- median",hc_thresh,ec_thresh,np.median(flat))

    ec_thresh = 0.6*hc_thresh

    #### we then need to represent as a graph and classify nodes based on intensity
    G = nx.Graph()

    nz_latt = intensity_lattice.shape[0]
    ny_latt = intensity_lattice.shape[1]
    nx_latt = intensity_lattice.shape[2]
    nodecount = 0
    for iz in range(nz_latt):
        for iy in range(ny_latt):
            for ix in range(nx_latt):
                ellfact =  ((-baseline_ella + (ix+0.5)*lattice_size)/baseline_ella)**2 
                ellfact += ((-baseline_ellb + (iy+0.5)*lattice_size)/baseline_ellb)**2 
                ellfact += ((-baseline_ellc + (iz+0.5)*lattice_size)/baseline_ellc)**2
##                if ellfact < 1:
##                    print(ix,iy,iz,ellfact,intensity_lattice[iz,iy,ix])
                ## figure out if intensity below HC thresh
                if intensity_lattice[iz,iy,ix] < ec_thresh and ellfact < 0.9:
                    G.add_node(nodecount)
                    attr = {nodecount: {"pos": (-baseline_ella + (ix+0.5)*lattice_size,
                                                -baseline_ellb + (iy+0.5)*lattice_size,
                                                -baseline_ellc + (iz+0.5)*lattice_size)
                                        }
                           }
                    nx.set_node_attributes(G,attr)
                    nodecount +=1


    print("created graph with nodes",nodecount)
    #### then see if we can make a low-intensity percolating cluster by adding relevant/chosen edges
    pairs = combinations(G.nodes,2)

    for pair in pairs:
        n1 = pair[0]
        n2 = pair[1]

        x1 = G.nodes[n1]['pos'][0]
        x2 = G.nodes[n2]['pos'][0]

        y1 = G.nodes[n1]['pos'][1]
        y2 = G.nodes[n2]['pos'][1]

        z1 = G.nodes[n1]['pos'][2]
        z2 = G.nodes[n2]['pos'][2]


        dist = np.sqrt((x1 - x2)**2 + (y1-y2)**2 + (z1-z2)**2)

        if dist < 0.05+np.sqrt(3.0)*lattice_size:
            G.add_edge(n1,n2,weight=dist)

    print("done adding edges")
##    print(G.edges())
    ### visualise the cluster of low intensity edges
    # Get node positions
    pos = nx.get_node_attributes(G, 'pos')
##    print(pos)
    
    # Get number of nodes
    n = G.number_of_nodes()

    # Get the maximum number of edges adjacent to a single node
##    edge_max = max([G.degree(i) for i in range(n)])
##
##    # Define color range proportional to number of edges adjacent to a single node
##    colors = [plt.cm.plasma(G.degree(i)/edge_max) for i in range(n)]

    ## define colors based on node centrality
    print(baseline_ella,baseline_ellb,baseline_ellc)
    print(G.nodes[0]['pos'][0],
          G.nodes[0]['pos'][1],
          G.nodes[0]['pos'][2],
          (G.nodes[0]['pos'][0]/baseline_ella)**2 +
          (G.nodes[0]['pos'][1]/baseline_ellb)**2 +
          (G.nodes[0]['pos'][2]/baseline_ellc)**2 - 1)

    print(G.nodes[10]['pos'][0],
          G.nodes[10]['pos'][1],
          G.nodes[10]['pos'][2],
          (G.nodes[10]['pos'][0]/baseline_ella)**2 +
          (G.nodes[10]['pos'][1]/baseline_ellb)**2 +
          (G.nodes[10]['pos'][2]/baseline_ellc)**2 - 1)

    
    centlist = [ (G.nodes[i]['pos'][0]/baseline_ella)**2 +
                 (G.nodes[i]['pos'][1]/baseline_ellb)**2 +
                 (G.nodes[i]['pos'][2]/baseline_ellc)**2 - 1 for i in range(n)]

##    print(np.array(centlist)[centlist > 0])
    colors = [ plt.cm.plasma(abs(cent)) for cent in centlist]


##    nx.draw_networkx(G,pos)
##    plt.savefig("3d_channel_net.png")
##    plt.show()
    # 3D network plot
    fig=go.Figure(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers", marker_size=1))
    with plt.style.context(('ggplot')):
        
##        fig = plt.figure(figsize=(10,7))
####        ax = Axes3D(fig)
##        ax = fig.add_subplot(projection='3d')
        
        # Loop on the pos dictionary to extract the x,y,z coordinates of each node
        for key, value in pos.items():
            xi = value[0]
            yi = value[1]
            zi = value[2]

##            print(key,xi,yi,zi)
            
            # Scatter plot
##            ax.scatter(xi, yi, zi, c=np.atleast_2d(np.array([colors[key]])), s=10+10*G.degree(key), edgecolors='k', alpha=0.7)
            fig.add_scatter3d(x=[xi], y=[yi], z=[zi], mode="markers", marker_size=10,marker_color='red',opacity=0.5)
        
        # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
        # Those two points are the extrema of the line to be plotted
        for i,j in enumerate(G.edges()):
##            print(i,j,pos[j[0]][0],pos[j[1]][0])
            x = np.array((pos[j[0]][0], pos[j[1]][0]))
            y = np.array((pos[j[0]][1], pos[j[1]][1]))
            z = np.array((pos[j[0]][2], pos[j[1]][2]))
        
        # Plot the connecting lines
##            ax.plot(x, y, z, c='black', alpha=0.5)
            line = go.Scatter3d(x=x,
                                y=y,
                                z=z,
                                mode='lines',
                                name='',
                                line=dict(color= 'rgba(50,50,50,0.65)', width=4.5))
            fig.add_trace(line)
    
    # Set the initial view
##    angle=120
##    ax.view_init(30, angle)
##    print("done making plot also")
##    # Hide the axes
##    ax.set_axis_off()
##    ax.set_xlim(-1.2*baseline_ella,1.2*baseline_ella)
##    ax.set_ylim(-1.2*baseline_ella,1.2*baseline_ella)
##    ax.set_zlim(-1.2*baseline_ella,1.2*baseline_ella)
##    save = False
##    if save is not False:
##        plt.savefig(outdir+ "/3d_channel_net.png")
##        plt.close('all')
##    else:
##        plt.tight_layout()
##        plt.savefig("3d_channel_net.png")
##        plt.show()


    #### do a plotly 3d graph
    fig.update_layout(width=600, height=600, showlegend=False)
    fig.update_layout(
    scene = dict(
        xaxis = dict(nticks=4, range=[-1.2*baseline_ella,1.2*baseline_ella],),
        yaxis = dict(nticks=4, range=[-1.2*baseline_ella,1.2*baseline_ella],),
        zaxis = dict(nticks=4, range=[-1.0*baseline_ellb,1.0*baseline_ellb],),),
    width=700,
    margin=dict(r=20, l=10, b=10, t=10))
    fig.show()

    

    
##    return

    ### how to determine if this is percolating from boundary to boundary
    ### we stuck with visualisation now


    ### but use node attributes -- largest path length compared to geometry of ellipsoid
##    print(centlist)
    edge_nodes = [i for i in range(len(centlist)) if centlist[i] > -0.15]
    print(len(edge_nodes))
##    print(G.nodes[edge_nodes])
##    print(centlist)
    
    distlist = np.zeros([len(edge_nodes)],np.float64)
    pairlist = np.zeros([len(edge_nodes)],np.int64)
    for i in range(len(edge_nodes)):
        spamarr = np.zeros([len(edge_nodes)],np.float64)
        for j in range(len(edge_nodes)):
            try:
                spamarr[j] = nx.shortest_path_length(G,i,j,weight='weight')
            except:
                spamarr[j] = -1
            
                                            
                       
##        spamdists,spampaths = nx.multi_source_dijkstra(G, , target=edge_nodes[i], cutoff=None, weight='weight')
##        print(spamdists,edge_nodes[i])
##        if spamdists == 0:
##            continue
        print(i,np.max(spamarr),np.min(spamarr))
        distlist[i] = np.max(spamarr)
        pairlist[i] = edge_nodes[np.argmax(spamarr)]

    maxdist = np.max(distlist)
    m1 = np.argmax(distlist)
    m2 = pairlist[m1]

    print(maxdist,m1,m2,nx_latt)


    ### add code to put this in a file (append)
        

        






    
