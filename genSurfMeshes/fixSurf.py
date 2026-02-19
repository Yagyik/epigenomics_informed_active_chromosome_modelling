import matplotlib.pyplot as plt
import numpy as np
import sys
import subprocess
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import plotly.graph_objects as go



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
    
##    if
    if np.any(dcent > 2*rad):
        print("vertices\n")
        print(vertices,area,pmarea,np.pi*rad*rad)
        print("centroid",centroid)
        print(vertices - centroid)
        print("norm of vertex centroid dist",np.linalg.norm(vertices - centroid,axis=1))
        print("distances",dcent,rad)

        newvertices = np.copy(vertices)
        change = np.argmax(dcent)
        newvertices[change] = np.copy(centroid)
        print("changing",change,np.where(dcent > rad))
        print(newvertices)
        print(vertices)

        pmab = newvertices[0,:] - newvertices[1,:]
        pmac = newvertices[2,:] - newvertices[1,:]
        pmarea = 0.5*np.linalg.norm(np.cross(pmab,pmac))
        ncentroid = np.mean(newvertices,axis=0)
        ndcent = np.linalg.norm(newvertices - ncentroid,axis=1)
        print(pmarea,ncentroid,ndcent)
    return area

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


surfhull=ConvexHull(surfpoints)


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

#### find triangle areas
print(utri)
print(tri_vertices)

##areas = np.zeros([len(utri)],np.float64)
areas = np.array([findArea(x,surfrad) for x in tri_vertices[:100]])
print(areas.shape)
plt.hist(areas, bins=np.arange(0,0.5,0.01), 
         color='blue', alpha=0.5, density=True, label='area histo')

plt.show()

Xe = []
Ye = []
Ze = []
for T in tri_vertices:
    Xe += [T[k%3][0] for k in range(4)]+[ None]
    Ye += [T[k%3][1] for k in range(4)]+[ None]
    Ze += [T[k%3][2] for k in range(4)]+[ None]



##lines= go.Scatter3d(
##                x=Xe,
##                y=Ye,
##                z=Ze,
##                mode='lines',
##                name='',
##                line=dict(color= 'rgba(50,50,50,0.5)', width=1.5),opacity=0.4) 
##fig=go.Figure(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers", marker_size=15,marker_color='blue',opacity=0.05))
####fig=go.Figure(go.Scatter3d(x=x, y=y, z=z, mode="markers", marker_size=15,marker_color='red',opacity=0.05))
##
####fig.add_scatter3d(x=x, y=y, z=z, mode="markers", marker_size=10,marker_color='red',opacity=0.15)
##fig.add_mesh3d(x=x, y=y, z=z, i=I, j=J, k=K, color='red',opacity=0.15) #'rgba(30, 218, 245, 0.35)')
##fig.add_trace(lines)
##fig.update_layout(width=800, height=800, showlegend=False)
##fig.update_layout(
##    scene = dict(
##        xaxis = dict(nticks=4, range=[-1.5*ella,1.5*ella],),
##        yaxis = dict(nticks=4, range=[-1.5*ella,1.5*ella],),
##        zaxis = dict(nticks=4, range=[-1.0*ellb,1.0*ellb],),),
##    width=700,
##    margin=dict(r=20, l=10, b=10, t=10))
##
##fig.show()

outfile = open(sys.argv[3],"w")

for i in range(nsurf):
    outfile.write("%d %f %f %f %f\n" %(i,surfrad,surfx[i],surfy[i],surfz[i]))
outfile.close()

