'''
Created on 08.08.2019

@author: marisa

TODO MrBayes:
+ check location for new nexus files
+ read them in dendropy and create a newick file for rooting and save them
'''


import glob, random,re, os
from dendropy import Tree, TreeList
from collections import defaultdict

def edit_treesample_nwk(method,pathBS):
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
        
        
def read_treesample_mb(pathBS):
    '''
    this method reads the posterior sample from MrBayes
    read the treefiles from the runs in mrbayes
    get the number of trees, draw 100 random numbers
    get the 100 trees matching the random numbers (number equals line count !!without the nexus head)
    save the nexus head
    create a new file out of nexus head and 100 trees
    read the new treesample in dendropy and convert nexus to newick, save the sample
    
    :param method:
    :param pathBS:
    '''
    file_list = glob.glob(pathBS)
    concept_list = set()
    for f in file_list:
        concept =  f.split("/")[-1].split(".")[0]
        concept_list.add(concept)
    #list(concept_list)
    ##create a dict with key=concept value=list of files to run the analysis for each concept
    file_dict = defaultdict()
    for concept in concept_list:
        file_dict[concept]=list()
        for i,f in enumerate(file_list):
            #print f.startswith(concept)
            if f.startswith("input/"+concept) == True:
                file_name = file_list[i]
                file_dict[concept].append(file_name)
    ##for concept and list of files, get the tree offset and number of trees for one file (both have the same number of trees)
    ##compute the random sample
    for concept, list_files in file_dict.items():
        print list_files
        numTrees, tree_offset, remaining_trees = read_single_file(list_files[0])
        print numTrees, tree_offset
        tree_sample = sorted(random.sample(xrange(tree_offset,numTrees*2), 5))
        print tree_sample
        ##for each value in the tree sample
        values1File = []
        values2File = []
        for v in tree_sample:
            ##if value is smaller or equal the number of trees in the file, take the trees from the first file
            if v <=numTrees:
                values1File.append(v)
                #print "in if "+ str(v)
                #print values1File
            ##else substract number of trees from the value to get the correct line in the sample
            else:
                v1 = v-numTrees
                values2File.append(v1)
                #print "in else " +str(v)+" "+str(v1)
                #print values2File
        pathToFile = create_new_nexus(concept,list_files, values1File, values2File)
        create_newick(pathToFile)
        

#############helper method mb##################
def read_single_file(f):
    '''
    read one file from one mr bayes run. get the number of generations/500+1 to get the number of trees in the file
    the number of trees are the same for all runs
    0.25*number of trees = the remaining 75% of the trees after burn-in (burn-in at 25% = 25% of the trees need to be disjected)
    '''
    ##read the tree file
    stdin,stdout = os.popen2("tail -n 2 "+f)
    stdin.close()
    line = stdout.readlines(); stdout.close()
    #print line[0].split()
    numGen = line[0].split()[1].split(".")[1]
    ##compute the number of trees from the generations
    numTrees = int(numGen)/500+1
    #print numTrees
    ##get the tree offset (number of trees after 25% tree from the burn-in are disjected)
    tree_offset = int(0.25*numTrees)
    ##get the remaining trees to draw the replicates from
    remaining_trees = numTrees - tree_offset
    ##return the tree offset
    return numTrees, tree_offset, remaining_trees

def create_new_nexus(concept,nexusfileList, value1File, value2File):
    '''
    create the new nexus file with the treesample
    the tree sample are 100 random posterior trees from both runs of MrBayes exluding the 25%burnin
    the trees are saved in a new nexus file
    :param concept:the name of the concept
    :param nexusfileList:the list with all nexus files from MrBayes
    :param value1File:the values for the trees for the first file
    :param value2File:the values for the trees for the second file
    :return pathToFile: the path to the new nexus file
    '''
    ##do everything for the first run
    with open(nexusfileList[0],"r") as f1:
        data1 = f1.readlines()
    ###for each line in the file, get the lines, which are not trees
    begining = []
    trees1 = []
    for l in data1:
        ##get the lines which are not trees (the begining of the nexus file including all white spaces, punctuaiton and newlines)
        if not l.startswith("   tree") and not l.startswith("end;"):
            begining.append(l)
        ##get all the lines which are trees including white spaces, punctuation and newlines
        elif l.startswith("   tree"):
            trees1.append(l)
            
    ###get the trees from the second file
    with open(nexusfileList[1],"r") as f2:
        data2 = f2.readlines()

    trees2 = []
    for l in data2:
        ##get all the lines which are trees including white spaces, punctuation and newlines
        if l.startswith("   tree"):
            trees2.append(l)
            
    ##create the tree sample for this file
    treesample=[]
    ##append the trees from the first file
    for v1 in value1File:
        treesample.append(trees1[v1])
    ##append the trees from the second file
    for v2 in value2File:
        treesample.append(trees2[v2])
    ##create the path for the outputfile and write it
    path = concept+"+treesample.nex"
    with open(path,"w") as outfile:
        for l in begining:
            outfile.write(l)
        for t in treesample:
            outfile.write(t)
        outfile.write("end;")
    ##return the path to the new nexus file
    return path
    

def create_newick(pathToFile):
    '''
    read the nexus file with the 100 sample trees for rooting
    create newick file using dendropy and save the newick file
    :param pathToFile:
    '''
    pass


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
        #outfile.write("\n")
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
    pathBS = "input/*.t"
    read_treesample_mb(pathBS)
    
    
    
    