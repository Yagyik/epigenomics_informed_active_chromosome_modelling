import numpy as np
import json
from ChromUtils.Structures import *


def parse_file(config_filename):
    with open(config_filename,'r') as config_file:
        params = json.load(config_file)
    return params


def read_config(conffile,time,nsurf,nchrom,npatch):
    """ code to read a configuration file and store the entries organised as per class config
        class config has sim params and a class chrom
        class chrom stores chrom centre, patch centre and radial information

        """

    
    infile = open(conffile,"r")
    lines = infile.readlines()
    infile.close()
##    print("first line >>>",lines[0])
    ##    print(lines)

    
    nato=nsurf+nchrom+nchrom*2*npatch
    
    surfx = np.zeros([nsurf],np.float64)
    surfy = np.zeros([nsurf],np.float64)
    surfz = np.zeros([nsurf],np.float64)
    surfid = np.zeros([nsurf],np.int64)
    sfid = 0
    surfrad = 0.0

    chrcentx = np.zeros([nchrom],np.float64)
    chrcenty = np.zeros([nchrom],np.float64)
    chrcentz = np.zeros([nchrom],np.float64)
    chrrad = np.zeros([nchrom],np.float64)
    patchlist = np.zeros([nchrom,npatch+1],np.int64)
    patchfundx = np.zeros([nchrom,npatch+1],np.float64)
    patchfundy = np.zeros([nchrom,npatch+1],np.float64)
    patchfundz = np.zeros([nchrom,npatch+1],np.float64)
    patchfundrad = np.zeros([nchrom,npatch+1],np.float64)
    patchx = np.zeros([nchrom,npatch+1],np.float64)
    patchy = np.zeros([nchrom,npatch+1],np.float64)
    patchz = np.zeros([nchrom,npatch+1],np.float64)
    patchvx = np.zeros([nchrom,npatch+1],np.float64)
    patchvy = np.zeros([nchrom,npatch+1],np.float64)
    patchvz = np.zeros([nchrom,npatch+1],np.float64)
    patchrad = np.zeros([nchrom,npatch+1],np.float64)
    patchpointer = np.zeros([nchrom,npatch+1],np.int64)
    count=0
    for line in lines:
        count+=1
        # use count to decide what to do with info on line
        if count==4: # check number of atoms
            tokens=line.split()
            if int(tokens[0]) != nato:
                print("check number of atoms\n")

        if count == 6: # check box length/radius
            tokens=line.split()
            ella = float(tokens[1]) # - float(tokens[0])
        if count == 7:
            tokens=line.split()
            ellb = float(tokens[1])
        if count == 8:
            tokens=line.split()
            ellc = float(tokens[1])

    ##        print(count,majorax,minorax)

        if count >=10 and count <10+nsurf:
            tokens=line.split()
            sfid = int(tokens[0])
            surfid[sfid] = sfid
            surfrad = float(tokens[2])
            surfx[sfid] = float(tokens[4])
            surfy[sfid] = float(tokens[5])
            surfz[sfid] = float(tokens[6])

        if count>=10+nsurf and count < 10+nato:
            tokens=line.split()
    ##        print(line)
            chrid = int(tokens[1])-1
            chrtype = int(tokens[3]) -1 # 1-4 -> 0-3
            checktype = chrid % int(0.5*nchrom) # 0-7 -> 0-3,0-3
    ##        print(chrid,chrtype,checktype)
    ##        print(int(tokens[0]),int(tokens[1]),int(tokens[3]),float(tokens[8])**2 + float(tokens[9])**2 + float(tokens[10])**2)
    ##        if int(tokens[1]) == int(tokens[3]): ## is a chrom
            if checktype == chrtype:
##                print("chrom line >>>",line)
                chrcentx[chrid] = float(tokens[5])
                chrcenty[chrid] = float(tokens[6])
                chrcentz[chrid] = float(tokens[7])
                chrrad[chrid] = float(tokens[4])
                patchlist[chrid][0] = int(tokens[2])
            elif checktype !=chrtype and float(tokens[8])**2 + float(tokens[9])**2 + float(tokens[10])**2 == 0.0: ## is a patch anchor
    ##            print("patch anchor >>>",int(tokens[1]),int(tokens[2]),float(tokens[8]),float(tokens[9]),float(tokens[10]))
                patchid = int(tokens[2])
                patchlist[chrid][patchid] = int(tokens[0])
                patchfundx[chrid][patchid] = float(tokens[5])
                patchfundy[chrid][patchid] = float(tokens[6])
                patchfundz[chrid][patchid] = float(tokens[7])
                patchfundrad[chrid][patchid] = float(tokens[4])
            elif checktype !=chrtype and float(tokens[8])**2 + float(tokens[9])**2 + float(tokens[10])**2 > 0.0: ## is a patch extrusion
    ##            print("patch line >>> ",line)
                patchid = int(tokens[2])
                patchlist[chrid][patchid] = int(tokens[0])
                patchx[chrid][patchid] = float(tokens[5])
                patchy[chrid][patchid] = float(tokens[6])
                patchz[chrid][patchid] = float(tokens[7])
                patchvx[chrid][patchid] = float(tokens[8])
                patchvy[chrid][patchid] = float(tokens[9])
                patchvz[chrid][patchid] = float(tokens[10])
                patchrad[chrid][patchid] = float(tokens[4])
                patchpointer[chrid][patchid] = int(tokens[3])


    config = Configuration(time,nsurf,nchrom,npatch,
                           surfrad,ella,ellb,ellc,
                           surfx,surfy,surfz,surfid,
                           chrcentx,chrcenty,chrcentz,chrrad,
                           patchlist,patchfundx,patchfundy,patchfundz,patchfundrad,
                           patchx,patchy,patchz,patchvz,patchvy,patchvz,
                           patchrad,patchpointer)

    return config



def readIntMat(intmatfile,config):
    infile = open(intmatfile,"r")
    lines = infile.readlines()
    infile.close()

    nchrom = config.nchrom
    npatch = config.npatch

    nato = nchrom+nchrom*npatch

    intMat = np.zeros([nato,nato],np.float64)
    linecount = 0
    for line in lines:
        tokens = line.split()
        par1id = int(tokens[0])
        par2id = int(tokens[1])
   
        intMat[par1id,par2id] = float(tokens[6])
        

    return intMat



def __init__():
    all = ['parse_file','read_config','readIntMat']






    
##def readConstraints(conffile):
    
