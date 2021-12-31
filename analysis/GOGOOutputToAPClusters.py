#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 15:14:14 2020

@author: benlitterer

INPUT: 
    MC30File: file containing the GOGO similarity for each pair of MC30 output gos 

Takes in a file with the similarities between annotatations for a protein and outputs clusters of these annotations 
"""

import sys 
import re
import pandas as pd
from goatools import obo_parser

cwd = "/home/benlitterer/Academic/Research/ProjectMatrices/referance" 
goannots = {}

#this is the path to the gene ontology basic version
go_obo = "/home/blitt/Academic/Research/ProjectMatrices"
go = obo_parser.GODag(go_obo + "/go-basic.obo")

#this gets the similarity of GO annotations for a given aspect of the GO hierarchy 
def getSimMatrix(aspectKey): 
    MC30File = open(sys.argv[1])
    simDict = {}
    for line in MC30File: 
        lineList = line.strip("\n").split(" ")
        if len(lineList) == 4: 
            go1 = lineList[0]
            go2 = lineList[1]
            aspect = lineList[2]
            
            #this should probably be named similarity, since it is actually how close two go terms are semantically 
            diff = float(lineList[3])
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
        
        return (simMatrix, fRowKeys)
    return (None, None, None)

#plotting imports 
from sklearn.cluster import AffinityPropagation 


#get clusters based off of the similarity dictionary/matrix
def getClusters(aspectSimDict, aspectKeys, namespace, fileEnding): 
    B62AspectAnnotNames = set([go[annot].name for annot in B62GOSet if go[annot].namespace == namespace])
    aspectNames = set([go[annot].name for annot in aspectKeys])
    print(fileEnding)
    print("missing: " + str(B62AspectAnnotNames - aspectNames))
    
    #scikit learn syntax 
    af = AffinityPropagation(affinity='precomputed', verbose=True, random_state=None)
    af.fit(aspectSimDict)
    #debugging
    #print(af.labels_)
    
    #these are the annotations chosen to represent given clusters 
    centerAnnots = [aspectKeys[index] for index in af.cluster_centers_indices_]
    centerNames = [go[annot].name for annot in centerAnnots]
    #debugging
    #print(len(aspectKeys))
    
    
    clustLabels = af.labels_
    numClusters = len(set(clustLabels))
    clustArr = []
    
    #output clusters to a file seperated by dashed lines 
    #we capitalize the representative annotation for a given cluster
    #we add asterisks if the annotation is also found using BLOSUM62
    for clustNum in range(0,numClusters): 
        exemplar = centerNames[clustNum]
        #thisClust = [exemplar]
        clustArr.append(exemplar.upper())
        for index in range(0, len(clustLabels)): 
            label = clustLabels[index]
            if label == clustNum:
                #print(str(clustNum) + " == " + str(label))
                currAnnot = aspectKeys[index]
                currName = go[currAnnot].name
                if currName != exemplar: 
                    if currName in B62AspectAnnotNames: 
                        currName = "*" + currName + "*"
                    clustArr.append(currName)
                        
        clustArr.append("------------------")
    
    outStr = "\n".join(clustArr) + "\n"
    #debugging
    #print("exemplars: " + str(centerNames))
    #print(outFileStem + fileEnding + ".tsv")
    
    #we automatically create an output file name 
    outFileName = re.split("(MC30|B62|Protsub)", sys.argv[1], flags=re.IGNORECASE)[0] + fileEnding + ".tsv"
    outFile = open(outFileName, "w")
    outFile.write(outStr)

#run the functions above to get clusters for each aspect of the hierarchy
#instead of using one argument, it might be best to use two, this would allow us to not assume that we can replace MC30 with B62 to get the 
#GOGO output for the BLOSUM62 files 
if ".txt" in sys.argv[1] and "MC30" in sys.argv[1]: 
    print(sys.argv[1])
    simDictBP, BPKeys = getSimMatrix("BPO")
    simDictMF, MFKeys = getSimMatrix("MFO")
    simDictCC, CCKeys = getSimMatrix("CCO")
    
    B62FilePath = sys.argv[1].replace("MC30", "B62")
    B62File = open(B62FilePath)
    B62GOSet = set([line.strip("\n").split(" ")[0] for line in B62File])
    
    #outFileStem = sys.argv[2]
        
    getClusters(simDictBP, BPKeys, "biological_process", "BP")
    getClusters(simDictMF, MFKeys, "molecular_function", "MF")
    getClusters(simDictCC, CCKeys, "cellular_component", "CC")
    print("")