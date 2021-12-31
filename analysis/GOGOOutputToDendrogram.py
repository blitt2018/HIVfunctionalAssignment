#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 15:05:11 2020

@author: benlitterer


INPUT: 
    MC30File - a file containing GOGO similarities for all pairs of GO terms in the MC30 output
    B62File - a file containing GOGO similarities for all pairs of GO terms in the B62 output
    
This file takes output from the GOGO algorithm and uses it to create dendrograms showing the similarity 
between all of the GO annotations output by the new (in this case MC30) matrix

an example folder containing GOGO Output is in /shareable/GOGOOutputs001/

If needed, you can use the first portion of this code to just get similarity matrices and then do another analysis (e.g. clustering)
"""

import sys 
from goatools import obo_parser

#this is the path to the gene ontology basic version
go_obo = "/home/benlitterer/Academic/Research/ProjectMatrices"
go = obo_parser.GODag(go_obo + "/go-basic.obo")


#takes aspectKey as input, which is the type of heirarachy (Biological Process, Molecular Function, Cellular Componenet) we are interested in 
#gets the similarity between annotations in the MC30 (new matrix) file, and also outputs the index labels for the rows/columns of this matrix as well as the 
#BLOSUM62 annotations 
def getSimMatrix(aspectKey): 
    MC30File = open(sys.argv[1])
    B62File = open(sys.argv[2])
    B62GOSet = set([line.strip("\n").split(" ")[0] for line in B62File])
    simDict = {}
    for line in MC30File: 
        lineList = line.strip("\n").split(" ")
        if len(lineList) == 4: 
            go1 = lineList[0]
            go2 = lineList[1]
            aspect = lineList[2]
            
            #NOTE: the distance between GO annotations is just 1 - similarity 
            #we can play with this distance to get a more desirable distribution
            #currently we have tons of annotations very far away and a few very close to eachother
            diff = 1- float(lineList[3])

            if aspect == aspectKey:
                if go1 not in simDict: 
                    simDict[go1] = {}
                simDict[go1][go2] = diff
    keys = list(simDict.keys())
    if len(keys) != 0: 
        fRowKeys = list(simDict[keys[0]].keys())
        simMatrix = []
        
        #we want the rows of our matrix to have the same order as the column 
        for key in fRowKeys: 
            simMatrix.append(list(simDict[key].values()))
        
        return (simMatrix, fRowKeys, B62GOSet)
    return (None, None, None)

simDictBP, BPKeys, BPB62Set = getSimMatrix("BPO")
simDictMF, MFKeys, MFB62Set = getSimMatrix("MFO")
simDictCC, CCKeys, CCB62Set= getSimMatrix("CCO")

#some reference code for this
#from https://stackoverflow.com/questions/41416498/dendrogram-or-other-plot-from-distance-matrix
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

#takes the similarity matrix/dictionary, the column/row labels for the matrix, the BLOSUM62 annotations, the title, and the name of the GO hierarchy aspect
#plots a dendrogram showing the distance between all annotations for a particular GO aspect
def plotDendo(aspectSimDict, aspectKeys, aspectB62Set, title, aspectName): 
    
    """
    heatmap, but need smaller amount of gos...
    plt.imshow(mat, cmap="Reds",interpolation="nearest")
    plt.xticks(range(len(BPKeys)), BPKeys, size="small", rotation=90)
    plt.yticks(range(len(BPKeys)), BPKeys, size="small", rotation=0)
    plt.show()
    """
    
    #create the dendrogram 
    dists = squareform(np.array(aspectSimDict))
    linkage_matrix = linkage(dists, "average")
    dendrogram(linkage_matrix, labels=aspectKeys, color_threshold=0.85, leaf_rotation=70, distance_sort="ascending")
    plt.title(title)
    plt.xlabel("GO annotations")
    plt.ylabel("Distance (GOGO algo)")
    plt.xticks(fontsize=5)
    
    #change color of individual labels this way https://stackoverflow.com/questions/42878462/matplotlib-changing-the-color-of-a-single-x-axis-tick-label
    #we want to find which labels of the dendrogram are annotations also found using BLOSUM62
    #the indeces of these labels are stored in B62indeces 
    dendoLabels = [label.get_text() for label in plt.gca().get_xticklabels()]
    B62GOList = list(aspectB62Set)
    B62indeces = [dendoLabels.index(annot) for annot in B62GOList if annot in dendoLabels]
    
    #the labels that actually belong to this aspect from B62
    trueB62labels  = set([annot for annot in aspectB62Set if go[annot].namespace == aspectName])
    totalLabels = set(dendoLabels)
    print(trueB62labels - totalLabels)
    [plt.gca().get_xticklabels()[index].set_color("red") for index in B62indeces]
    plt.show()
    
    #NOTE: some code can be added to save the dendrograms here, or they can be manually saved when they display with matplotlib 
    
    """
    verifying that simMatrix is working and has same indeces 
    print(fRowKeys)
    print(simMatrix[:3])
    print(simDictBP["GO:0015074"]["GO:0039520"])
    print(simDictBP["GO:0039520"]["GO:0044826"])
    """

#plot dendrograms for each apsect of the hierarchy
plotDendo(simDictBP, BPKeys, BPB62Set, "Envelope protein: Biological Process (MC30)", "biological_process")
plotDendo(simDictMF, MFKeys, MFB62Set, "Spike Protein: Molecular Function (MC30)", "molecular_function")
plotDendo(simDictCC, CCKeys, CCB62Set, "Envelope protein: Cellular Component (MC30)","cellular_component")