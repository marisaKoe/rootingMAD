'''
Created on 08.08.2019

@author: marisa
'''

import glob, os, subprocess
#import edit_negativeBL

def MAD_rooting_treesample(method):
    '''
    calls MAD for rooting
    all files are in one folder
    the output is written in the same folder with the ending .rooted
    :param method: the name of the method to get the right trees
    '''
    ##pmi path
    #path = glob.glob("/home/marisa/Dropbox/EVOLAEMP/projects/Project-Borrowing-hgt/rootingMAD/NELex/"+method+"/*.nwk")
    path = glob.glob("/home/marisakoe/Dropbox/EVOLAEMP/projects/Project-ConceptTrees-CharacterBased/NELex/ML_iqtree/NgramsNW/bootstrapReplicates/*.boottrees")
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
    path = glob.glob(path+"*.treefile")
    for filename in path:
        p = subprocess.Popen('mad '+filename+' -n',shell=True)
        os.waitpid(p.pid,0)

if __name__ == '__main__':
    
    ###pmi fastme trees (single trees)
    #method = "PMI_multipleData"
    #inpath = "/home/marisakoe/Dropbox/EVOLAEMP/projects/Project-ConceptTrees-DistanceMethods/NELex/PMI_based_methods/PMI_multipleData/phylipTrees/fastme/*.nwk"
    #outpath = "/home/marisakoe/Dropbox/EVOLAEMP/projects/Project-ConceptTrees-DistanceMethods/NELex/PMI_based_methods/PMI_multipleData_rootedTrees/"
    #edit_negativeBL.edit_single_tree(inpath, outpath)
    ####iqtree
    outpath = "/home/marisakoe/Dropbox/EVOLAEMP/projects/Project-ConceptTrees-CharacterBased/NELex/ML_iqtree/NgramsNW/iqTrees/"
    MAD_rooting_singleTree(outpath)
    
    ##only do for the optimal methods, since they are used for the HGT transfer
    ##pmi bootstrapping with Noise
    #methodBS = "pmiMultidata"
    #pathBS = "/home/marisa/Dropbox/EVOLAEMP/projects/Project-BootstrappingWithNoise/nelex/"+methodBS+"/trees/*.nwk"
    #edit_negativeBL.edit_treesample(methodBS, pathBS)
    #MAD_rooting_treesample(methodBS)
    ####iqtree
    methodBS = "ML_ngramsNW"
    MAD_rooting_treesample(methodBS)
    
    
    
    
    
    