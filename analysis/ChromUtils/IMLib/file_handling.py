import json
import os
import numpy as np
from ChromUtils.conf_reader import read_config
from ChromUtils.conf_reader import readIntMat

def parse_file(filepath):
    with open(filepath, 'r') as file:
        return json.load(file)


def makepaths(anapath, confpath, vastr):
    pathsplit = confpath.split("/")
    confdir = ""
    print("old confpath", pathsplit)
    pathsplit[:] = [x for x in pathsplit if x]
    
#     print("new", pathsplit)
    for p in pathsplit[:-1]:
        confdir += "/" + p
    maindir = pathsplit[-1]
    maindirsplit = maindir.split("_")
    mainpath = ""
    for sp in maindirsplit[:-1]:
        mainpath += sp + "_"

    mainpath += vastr
    confdir += "/" + mainpath
    
    #### now doing ana
    
    pathsplit = anapath.split("/")
    print("old anapath", pathsplit)
    pathsplit[:] = [x for x in pathsplit if x]
    
    # print("new", pathsplit)
    anadir = ""
    for p in pathsplit[:-1]:
        anadir += p + "/"
    
#     print(pathsplit)
    maindir = ""
    index = -1
    while maindir == "":
        maindir = pathsplit[index]
        index += -1
    print(maindir)
    maindirsplit = maindir.split("_")
    mainpath = ""
    for sp in maindirsplit[:-1]:
        mainpath += sp + "_"
    
    mainpath += vastr
    anadir += "/" + mainpath
    print(confdir)
    print(anadir)
    return confdir, anadir

def make_allpaths(anapath,confpath,combined_str):
    confpathsplit = confpath.split("/")
    confdir = ""
    print("old confpath", confpathsplit)
    confpathsplit[:] = [x for x in confpathsplit if x]
    
#     print("new", pathsplit)
    for p in confpathsplit[:-1]:
        confdir += "/" + p
    confmaindir = confpathsplit[-1]
    print("confmaindir",confmaindir)


    anapathsplit = anapath.split("/")
    print("old anapath", anapathsplit)
    anapathsplit[:] = [x for x in anapathsplit if x]
    
    # print("new", pathsplit)
    anadir = ""
    for p in anapathsplit[:-1]:
        anadir += p + "/"
    
    print(anapathsplit)
    anamaindir = ""
    index = -1
    while anamaindir == "":
        anamaindir = anapathsplit[index]
        index += -1
    print("anamaindir",anamaindir)
    
    confmaindir_parts = confmaindir.split('_')
    anamaindir_parts = anamaindir.split('_')

    all_paths = []
    for comb in combined_str[:]:
        # parts = comb.split('_')
        # print(parts)
        for key, value in comb.items():
            for i, part in enumerate(confmaindir_parts):
                if key in part and all(c.isdigit() or c == '.' for c in part.replace(key, '')):
                    confmaindir_parts[i] = f'{key}{value}'
            for i, part in enumerate(anamaindir_parts):
                if key in part and all(c.isdigit() or c == '.' for c in part.replace(key, '')):
                    anamaindir_parts[i] = f'{key}{value}'

        new_confmaindir = '_'.join(confmaindir_parts)
        new_anamaindir = '_'.join(anamaindir_parts)

        # new_confmaindir = f"{confmaindir.split('_')[0]}_{dsa_val}_{eij_val}_{as_val}_{sa_val}_{tl_val}"
        # new_anamaindir = f"{anamaindir.split('_')[0]}_{dsa_val}_{eij_val}_{as_val}_{sa_val}_{tl_val}"

        new_confpath = "/" + "/".join(confpathsplit[:-1] + [new_confmaindir])
        # if anapathsplit[0] != "." and dirprefix is None:
        #     new_anapath = "/" + "/".join(anapathsplit[:-1] + [new_anamaindir])
        # else:
        #     new_anapath = dirprefix.join(anapathsplit[:-1] + [new_anamaindir])
        new_anapath = "/" + "/".join(anapathsplit[:-1] + [new_anamaindir])

        print("new conf",new_confpath)
        print("new ana",new_anapath)

        all_paths.append((new_confpath, new_anapath))

    return all_paths

    
    
    # return all_paths


def get_conf(confpath, vastr, nfiles, simparams, start, step, stop):
    """
    Reads configuration files and returns the interaction matrix and reference configuration.

    Parameters:
    - confpath: str - Path to the configuration files
    - vastr: str - Variation string
    - nfiles: int - Number of files to read
    - simparams: dict - Simulation parameters
    - start: int - Start time
    - step: int - Time step
    - stop: int - Stop time

    Returns:
    - intMat: np.ndarray - Interaction matrix
    - refconfig: object - Reference configuration
    """
    nsurf = simparams["nsurf"]
    nchrom = simparams["nchrom"]
    npatch = simparams["npatch"]
    infileprefix = simparams["infileprefix"]
    nato = nsurf + nchrom + nchrom * 2 * npatch  # 1600 + 8 + 8x2x4 OR 2000 + 46 + 46x2x23

    print(nfiles, "files")
    configs = np.empty(nfiles, dtype=object)

    for i in range(nfiles):
        t = start + i * step
        filename = os.path.join(confpath, f"{infileprefix}{t}.dat")
        configs[i] = read_config(filename, t, nsurf, nchrom, npatch)

    if not os.path.exists(simparams["intmatfile"]):
        raise FileNotFoundError(f"Interaction matrix file {simparams['intmatfile']} not found.")
    
    intMat = readIntMat(simparams["intmatfile"], configs[0])
    refconfig = configs[0]
    
    return intMat, refconfig