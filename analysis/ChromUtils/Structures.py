import numpy as np




class Configuration(object):

    def __init__(self,timestep,nsurf,nchrom,npatch,
                 surfrad,ella,ellb,ellc,
                 surfx,surfy,surfz,surfid,
                 chrcentx,chrcenty,chrcentz,chrrad,
                 patchlist,patchfundx,patchfundy,patchfundz,patchfundrad,
                 patchx,patchy,patchz,patchvx,patchvy,patchvz,patchrad,patchpointer):
##        self.surfx = np.zeros([nsurf],np.float64)
##        self.surfy = np.zeros([nsurf],np.float64)
##        self.surfz = np.zeros([nsurf],np.float64)
##        self.surfid = np.zeros([nsurf],np.int64)
        self.nchrom = nchrom
        self.nsurf = nsurf
        self.npatch = npatch
        self.timestep = timestep
        self.surfx = surfx
        self.surfy = surfy
        self.surfz = surfz
        self.surfid = surfid
        self.surfrad = surfrad
        self.ella = ella
        self.ellb = ellb
        self.ellc = ellc
        


##        chromarr = np.empty(nchrom,dtype = object)
        self.chromarr = np.array([],dtype = object)
        basecount = 0

        for i in range(nchrom):
##            chromarr[i] = Chromosome(timestep,i,nchrom,npatch,chrcentx[i],chrcenty[i],chrcentz[i],chrrad[i],
##                                       patchlist[i],patchfundx[i],patchfundy[i],patchfundz[i],patchfundrad[i],
##                                       patchx[i],patchy[i],patchz[i],patchrad[i],patchpointer[i],basecount)
            spam = Chromosome(timestep,i,nchrom,npatch,chrcentx[i],chrcenty[i],chrcentz[i],chrrad[i],
                              patchlist[i],patchfundx[i],patchfundy[i],patchfundz[i],patchfundrad[i],
                              patchx[i],patchy[i],patchz[i],patchvx[i],patchvy[i],patchvz[i],
                              patchrad[i],patchpointer[i],basecount)
##            print(i,spam)
            self.chromarr = np.append(self.chromarr,spam)
##            print("pre ",basecount,patchpointer[i].shape[0])
            basecount += patchpointer[i].shape[0]
##            print("post ",basecount,patchpointer[i].shape[0])
##        print(self.chromarr.shape)                              


##class Chromosome(Configuration):
class Chromosome(object):


    def __init__(self,timestep,chrid,nchrom,npatch,chrcentx,chrcenty,chrcentz,chrrad,
                 patchlist,patchfundx,patchfundy,patchfundz,patchfundrad,
                 patchx,patchy,patchz,patchvx,patchvy,patchvz,patchrad,patchpointer,basecount):
        self.chrid = chrid
        self.chrtype = chrid %(int(0.5*nchrom))
        self.timestep = timestep
        self.npatch = npatch
        self.timestep = timestep
        self.chrcentx = chrcentx
        self.chrcenty = chrcenty
        self.chrcentz = chrcentz

        self.chrrad = chrrad

        self.patchlist = patchlist
        self.patchfundx = patchfundx
        self.patchfundy = patchfundy
        self.patchfundz = patchfundz
        self.patchfundrad = patchfundrad

        
        self.patchx = patchx
        self.patchy = patchy
        self.patchz = patchz
        self.patchvx = patchvx
        self.patchvy = patchvy
        self.patchvz = patchvz
        self.patchrad = patchrad
        self.patchpointer = patchpointer
##        print(self.chrid,self.chrtype,basecount,self.patchpointer.shape[0],np.max(self.patchpointer))

        self.gcount = np.linspace(basecount,basecount+patchpointer.shape[0]-1,patchpointer.shape[0],dtype=np.int64)
##        print(self.gcount,self.gcount.shape)

##        super().__init__(*args, **kwargs)
        

##        self.patchlist = np.zeros([npatch+1],np.int64)
##        self.patchfundx = np.zeros([npatch+1],np.float64)
##        self.patchfundy = np.zeros([npatch+1],np.float64)
##        self.patchfundz = np.zeros([npatch+1],np.float64)
##        self.patchfundrad = np.zeros([npatch+1],np.float64)
##
##        
##        self.patchx = np.zeros([npatch+1],np.float64)
##        self.patchy = np.zeros([npatch+1],np.float64)
##        self.patchz = np.zeros([npatch+1],np.float64)
##        self.patchrad = np.zeros([npatch+1],np.float64)
##        self.patchpointer = np.zeros([npatch+1],np.int64)

class referenceQuants(object):

    def __init__(self):
        self.lengthlist=[24.8 ,24.3 ,19.8 ,19.1 ,18.1 ,17.1 ,15.9 ,14.6, 14.1, 13.6, 13.5,
            13.4, 11.5, 10.7, 10.3, 9, 8.1, 7.8, 5.9, 6.3, 4.8, 5.1, 15.5,
            24.8,24.3 ,19.8 ,19.1 ,18.1 ,17.1 ,15.9 ,14.6, 14.1, 13.6, 13.5,
            13.4, 11.5, 10.7, 10.3, 9, 8.1, 7.8, 5.9, 6.3, 4.8, 5.1, 15.5]

        self.ecfraclist=[0.60,0.52,0.63,0.63,0.59,0.64,0.51,0.66,0.62,0.596,0.65,
            0.44,0.49,0.62,0.72,0.67,0.73,0.67,0.61,0.54,0.646,0.42,0.5,
            0.60,0.52,0.63,0.63,0.59,0.64,0.51,0.66,0.62,0.596,0.65,
            0.44,0.49,0.62,0.72,0.67,0.73,0.67,0.61,0.54,0.646,0.42,0.5]

        
        self.explengthlist=[24.3,19.1,14.6,11.5,10.7,6.3]
        self.expecfraclist=[0.52,0.63,0.66,0.49,0.62,0.54]
        self.expidlist=[2,4,8,13,14,20]

        self.expRglist=[3.10793471445247, 2.82976842052184,
                        2.8714241285949 , 2.35809288836717,
                        2.43914060391052, 2.78441863921879]

        self.expkappasqlist=[0.331649822034014, 0.244570725050393,
                             0.19792293250514, 0.238979172798409,
                             0.227756935517807, 0.224796702887166]

        self.expr1list = [2.61,2.29,2.26,1.90,1.96,2.23]
        self.expr2list = [1.47,1.34,1.43,1.16,1.18,1.33]
        self.expr3list = [0.82,0.97,1.03,0.78,0.85,0.99]

class Constraints(object):

    def __init__(self,constraint_arr):
        self.constraint_arr = constraint_arr

       


def __init__():
    all = ['Configuration','Chromosome','referenceQuants','Constraints']