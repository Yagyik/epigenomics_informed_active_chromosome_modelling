import numpy as np


def stackDim(inputarr,resarr = np.array([],np.float64)):
    """ subroutine to stack dimensions on numpy arrays
    """

    size = inputarr.shape
    ressize = resarr.shape[0]
    
    stackto = ressize + 1
    targshape = list(size)
    targshape.insert(0,1)
    spam = np.zeros(targshape,np.float64)
    spam[0] = np.copy(inputarr)

##    spam = np.reshape(inputarr,targshape)
    if ressize == 0:
        resarr = spam
    else:
        resarr = np.vstack((resarr,spam))

    return resarr



def __init__():
    all = ['stackDim']
   
    
