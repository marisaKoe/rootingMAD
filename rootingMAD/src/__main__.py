'''
Created on 08.08.2019

@author: marisa
'''

import glob, os, subprocess
from dendropy import Tree


def MAD_rooting_treesample(method):
    '''
    calls MAD for rooting
    all files are in one folder
    the output is written in the same folder with the ending .rooted
    :param method: the name of the method to get the right trees
    '''
    ##pmi path
    #path = glob.glob("yourFiles/*.nwk")
    ##ml
    #path = glob.glob("yourFiles/*.boottrees")
    ##mb
    path = glob.glob("yourFiles/*.nwk")
    for filename in path:
        p = subprocess.Popen('mad '+filename+' -n',shell=True)
        os.waitpid(p.pid,0)

def MAD_rooting_singleTree(path):
    '''
    calls MAD for rooting
    all files are in one folder
    the output is written in the same folder with the ending .rooted
    :param method: the name of the method to get the right trees
    '''
    ###pmi = *.nwk | iqtree = *.treefile
    path = glob.glob(path+"*.nwk")
    for filename in path:
        p = subprocess.Popen('mad '+filename+' -n',shell=True)
        os.waitpid(p.pid,0)


def convert_singleTreeMB(path):
    '''
    converts nexus (sumtrees) to newick for rooting
    saves the files in the same folder with ending .nwk
    calls MAD for rooting all files in the folder
    output is written in the same folder with the ending .rooted
    :param path:
    '''
    path = glob.glob(path+"*.tre")
    ##converts the file to newick and saves the trees
    for filename in path:
        t = Tree.get(path=filename, schema="nexus")#,rooting="default-unrooted")#,edge_length_type=float,suppress_edge_lengths=False)
        newfilename = filename.split(".")[0]
        newfilename = newfilename+".nwk"
        newT = t.as_string(schema="newick",suppress_rooting=True)
        newT = newT[:-6]+";"
        with open(newfilename,"w") as fout:
            fout.write(newT)
        
    
     

if __name__ == '__main__':
    ###########################single trees##########################
    
    ###pmi fastme trees (single trees)
    #method = "PMI_multipleData"

    ##mrBayes MCC trees
    #convert_singleTreeMB(outpath)
    #MAD_rooting_singleTree(outpath)
    
    
    ####################multiple trees########################
    
    ##only do for the optimal methods, since they are used for the HGT transfer
    ##pmi bootstrapping with Noise
    #methodBS = "pmiMultidata"
    #MAD_rooting_treesample(methodBS)
    ####iqtree
    #methodBS = "ML_ngramsNW"
    #MAD_rooting_treesample(methodBS)
    ###MrBayes
    methodBS = "MB_NW"
    MAD_rooting_treesample(methodBS)
    
    
    
    
    
