'''
Created on 08.08.2019

@author: marisa
'''


import glob, random,re
from dendropy import Tree, TreeList

def edit_treesample(method,pathBS):
    '''
    read tree files with all replicates in one file
    choose 100 random trees for the hgt detection with TRex
    set the negative branch length to 0.0
    save the trees in the input folder
    :param method: the method for which this is done
    :param pathBS: the path for the bootstrapping trees
    '''
    ##regular expresion to set negative branch length to 0.0
    regex = re.compile(r'-\d\.\d+(?:\.\d+)?')
    list_files = glob.glob(pathBS)
    #list_files = ["nelexAsjp_Berg::N+fastmeTree.nwk"]
    for f in list_files:
        ##get concept for distance based trees !!!!check file name of others
        concept = f.split("/")[-1].split("+")[0]
        ##open file, read lines and take 100 random trees from the sample
        with open(f,'r') as file1:
            ##read lines
            data = file1.readlines()
        tree_sample = random.sample(data, 100)
        ##save files with 100 tree samples and without negative branch length
        path = "/home/marisa/Dropbox/EVOLAEMP/projects/Project-Borrowing-hgt/inputTrees/"+method+"/"+concept+"+treesample.nwk"
        outfile = open(path,"w")
            
        for l in tree_sample:
            #result = [float(d) for d in re.findall(regex, l)]
            l = re.sub(regex, str(0.0), l)
            outfile.write(l)
        outfile.write("\n")
        outfile.close()


#####################edit branch length for single tree#################

def edit_single_tree(inpath,outpath):
    '''
    read the tree file
    set the negative branch length to 0.0 using regex
    save the tree in a folder
    :param inpath:the input path for the tree
    :param path:the output path for the tree
    '''
    ##regular expresion to set negative branch length to 0.0
    regex = re.compile(r'-\d\.\d+(?:\.\d+)?')
    list_files = glob.glob(inpath)
    #list_files = ["nelexAsjp_Berg::N+fastmeTree.nwk"]
    for f in list_files:
        print f
        ##get concept for distance based trees !!!!check file name of others
        concept = f.split("/")[-1].split("+")[0]
        print concept
        ##open file, read lines and take 100 random trees from the sample
        with open(f,'r') as file1:
            ##read lines
            data = file1.readlines()
        ##save files with 100 tree samples and without negative branch length
        path = outpath+concept+".nwk"
        outfile = open(path,"w")
            
        for l in data:
            #result = [float(d) for d in re.findall(regex, l)]
            l = re.sub(regex, str(0.0), l)
            outfile.write(l)
        outfile.write("\n")
        outfile.close()


###################add constant to branch length#######################


def edit_BL_constant(method,pathBS):
    '''
    read tree files with all replicates in one file
    choose 100 random trees for the hgt detection with TRex
    set the negative branch length to 0.0
    save the trees in the input folder
    :param method: the method for which this is done
    :param pathBS: the path for the bootstrapping trees
    '''
    list_files = glob.glob(pathBS)
    #list_files = ["nelexAsjp_Berg::N+fastmeTree.nwk"]
    for f in list_files:
        ##get concept for distance based trees !!!!check file name of others
        concept = f.split("/")[-1].split("+")[0]
        ##open file, read lines and take 100 random trees from the sample
        with open(f,'r') as file1:
            ##read lines
            data = file1.readlines()
        tree_sample = random.sample(data, 100)

        ##initialize tree list
        treelist = TreeList()
        ##for each tree in the 100 sample, read the tree in dendropy and append it to a TreeList object
        for t in tree_sample:
            t1 = Tree.get(data=t,schema="newick")
            treelist.append(t1)
         
        ##edit the trees by adding a constant (1.0) to account for negative branche length (do this for all methods?)
        add_constant_BL(treelist, method, concept)

def add_constant_BL(treelist, method,concept):
    '''
    add a constant (1.0) to the branch length to account of negative branch length.
    read the bootstrap replicates from the distance based-methods and edit all negative branch lenght in all trees.
    save the trees to the folder.
    '''
    #treelist = TreeList.get(path="/home/marisakoe/Dropbox/EVOLAEMP/projects/Project-Borrowing-hgt/hgt-test/nelexAsjp_Berg::N+allTrees.nwk", schema="newick")
    #print len(treelist)
    treelist_corrected = TreeList()
    count = 0
    for t in treelist:
        count +=1
        for eg in t.postorder_edge_iter():
            if not str(eg.length) == "None":
                eg.length += 1.0
        treelist_corrected.append(t)
    
    #print count
    #print len(treelist_corrected)
    #treelist_corrected.write(path="/home/marisakoe/Dropbox/EVOLAEMP/projects/Project-BootstrappingWithNoise/nelex/"+method+"/Trees_corrected/"+concept, schema="newick")
    treelist_corrected.write(path="/home/marisa/Dropbox/EVOLAEMP/projects/Project-Borrowing-hgt/inputTrees/"+method+"/"+concept+"+treesample.nwk", schema="newick")


if __name__ == '__main__':
    pass