import numpy as np




def grid(inp_surfpoints,
         ella,ellb,ellc,
         mod_ella,mod_ellb,mod_ellc,
         surfrad, nGrid):
    #### given points and geometry, place each in grid
    #### returns grid indices and grid-frame positions of each

    ### true geometry -- ella,ellb,ellc
    ### current geometry -- mod_ella,mod_ellb,mod_ellc
    boxLx = 2*(ella+2)
    boxLy = 2*(ellb+2)
    boxLz = 2*(ellc+2)
	
    grid_indices = np.zeros_like(inp_surfpoints,np.int64)
    gridcent_surfpoints = np.zeros_like(inp_surfpoints,np.float64)
    print(inp_surfpoints.shape,grid_indices.shape)

    for i in range(grid_indices.shape[0]):
        grid_indices[i,0] = int(nGrid*(inp_surfpoints[i,0]*ella/mod_ella + 0.5*boxLx)/boxLx) ## scaled to original geometry
        grid_indices[i,1] = int(nGrid*(inp_surfpoints[i,1]*ellb/mod_ellb + 0.5*boxLy)/boxLy)
        grid_indices[i,2] = int(nGrid*(inp_surfpoints[i,2]*ellc/mod_ellc + 0.5*boxLz)/boxLz)
    
        ## find positions wrt centre of box
        ## first find box centre
        

        boxcentx = -0.5*boxLx + (grid_indices[i,0]+0.5)*boxLx/nGrid
        boxcenty = -0.5*boxLy + (grid_indices[i,1]+0.5)*boxLy/nGrid
        boxcentz = -0.5*boxLz + (grid_indices[i,2]+0.5)*boxLz/nGrid

        


        ## find distance from box centre

        gridcent_surfpoints[i,0] = (inp_surfpoints[i,0]*ella/mod_ella - boxcentx)*nGrid/boxLx ## distance from centre as fraction of box length
        gridcent_surfpoints[i,1] = (inp_surfpoints[i,1]*ellb/mod_ellb - boxcenty)*nGrid/boxLy # distance from centre as fraction of box length
        gridcent_surfpoints[i,2] = (inp_surfpoints[i,2]*ellc/mod_ellc - boxcentz)*nGrid/boxLz # distance from centre as fraction of box length

    boxcentx = -0.5*boxLx + (grid_indices[0,0]+0.5)*boxLx/nGrid
    boxcenty = -0.5*boxLy + (grid_indices[0,1]+0.5)*boxLy/nGrid
    boxcentz = -0.5*boxLz + (grid_indices[0,2]+0.5)*boxLz/nGrid
    
    print("inp grid",inp_surfpoints[0,:],[boxcentx,boxcenty,boxcentz],mod_ella,mod_ellb,mod_ellc)
    print("proper check",inp_surfpoints[0,:],inp_surfpoints[0,0]*ella/mod_ella,inp_surfpoints[0,1]*ellb/mod_ellb,inp_surfpoints[0,2]*ellc/mod_ellc)
    return grid_indices,gridcent_surfpoints

    


def deform_from_grid(inp_surfpoints,grid_indices,gridcent_surfpoints,
                     surfrad,ella,ellb,ellc,mod,nGrid):
    #### given points with grid-frame locations and deformation multipliers
    #### get new positions and return
    print(mod,ella,ellb,ellc)
##    mod_ella = mod[0]*ella
##    mod_ellb = mod[1]*ellb
##    mod_ellc = mod[2]*ellc
    boxLx = 2*mod[0]*(ella+2)
    boxLy = 2*mod[1]*(ellb+2)
    boxLz = 2*mod[2]*(ellc+2)

    new_surfpoints = np.zeros_like(inp_surfpoints,np.float64)
    

    #printf("scaling surf points to %f %f %f %d\n",p->ella,p->ellb,p->ellc,p->nGrid);
    for i in range(inp_surfpoints.shape[0]):

        

        new_surfpoints[i,0] = -0.5*boxLx + (grid_indices[i,0]+0.5)*boxLx/nGrid + gridcent_surfpoints[i,0]*boxLx/nGrid
        new_surfpoints[i,1] = -0.5*boxLy + (grid_indices[i,1]+0.5)*boxLy/nGrid + gridcent_surfpoints[i,1]*boxLy/nGrid
        new_surfpoints[i,2] = -0.5*boxLz + (grid_indices[i,2]+0.5)*boxLz/nGrid + gridcent_surfpoints[i,2]*boxLz/nGrid

##        if mod[0] == mod[1] and mod[0] == mod[2]:
##            print("new",new_surfpoints[i,0],new_surfpoints[i,1],new_surfpoints[i,2])

    print("old vs new",inp_surfpoints[-1,:],new_surfpoints[-1,:],mod)  

    return new_surfpoints
