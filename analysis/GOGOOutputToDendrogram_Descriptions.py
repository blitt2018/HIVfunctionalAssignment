#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 12:30:16 2020

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

cwd = "/home/benlitterer/Academic/Research/ProjectMatrices/referance" 
goannots = {}

#this is the path to the gene ontology basic version
go_obo = "/home/blitt/Academic/Research/ProjectMatrices"
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
            #diff = float(lineList[3])
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

#plotting imports 
import numpy as np

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

#takes the similarity matrix/dictionary, the column/row labels for the matrix, the BLOSUM62 annotations, the title, and the name of the GO hierarchy aspect
#plots a dendrogram showing the distance between all annotations for a particular GO aspect
def plotDendo(aspectSimDict, aspectKeys, aspectB62Set, title, aspectName): 
    #from https://stackoverflow.com/questions/41416498/dendrogram-or-other-plot-from-distance-matrix
    #print names of gos that were missed by MC30
    fig = plt.figure(figsize=(18, 12))
    
    #get the annotations from this 
    trueB62GOs  = set([annot for annot in aspectB62Set if go[annot].namespace == aspectName])
    missedGOs = trueB62GOs - set(aspectKeys)
    print([go[annot].name for annot in missedGOs])
    
    #create the dendrogram
    dists = squareform(np.array(aspectSimDict))
    linkage_matrix = linkage(dists, "average")
    
    #the only difference between GOGOOutputToDendrogram.py and this file is that we get the name of annotations rather than GO:XXXXXX 
    aspectLabels = [go[key].name for key in aspectKeys]
    dendrogram(linkage_matrix, labels=aspectLabels, color_threshold=0, above_threshold_color="black", leaf_rotation=0, distance_sort="ascending", orientation="left")
    plt.title(title, fontsize=15)
    plt.xlabel("Distance (GOGOSim, Avg. Linkage)", fontsize=15)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=17)
    
    #change color of individual labels this way https://stackoverflow.com/questions/42878462/matplotlib-changing-the-color-of-a-single-x-axis-tick-label
    #we want to find which labels of the dendrogram are annotations also found using BLOSUM62
    #the indeces of these labels are stored in B62indeces 
    dendoLabels = [label.get_text() for label in plt.gca().get_yticklabels()]
    B62LabelList = [go[annot].name for annot in list(aspectB62Set)]
    B62indeces = [dendoLabels.index(annot) for annot in B62LabelList if annot in dendoLabels]
    [plt.gca().get_yticklabels()[index].set_color("red") for index in B62indeces]
    plt.tight_layout()
    plt.show()
    
    #NOTE: some code can be added to save the dendrograms here, or they can be manually saved when they display with matplotlib 
    
    """
    print(len(B62LabelList))
    print(len(dendoLabels))
    print(len(set(B62LabelList).intersection(set(dendoLabels))))
    """

#get dendrogram for each aspect 
plotDendo(simDictBP, BPKeys, BPB62Set, "Envelope protein: Biological Process (MC30)", "biological_process")
plotDendo(simDictMF, MFKeys, MFB62Set, "Envelope protein: Molecular Function (MC30)", "molecular_function")
plotDendo(simDictCC, CCKeys, CCB62Set, "Envelope protein: Cellular Component (MC30)", "cellular_component")