#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 12:27:10 2021

@author: blitt
"""

from goatools import obo_parser
import sys
import os 
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

#need to parse the file go-basic.obo
go_obo = "/home/blitt/Academic/Research/ProjectHIV"
go = obo_parser.GODag(go_obo + "/go-basic.obo")

'''
    This function finds the common ancestors in the GO 
    tree of the list of terms in the input.
'''
def common_parent_go_ids(terms, go):
    # Find candidates from first
    rec = go[terms[0]]
    candidates = rec.get_all_parents()
    candidates.update({terms[0]})
    
    # Find intersection with second to nth term
    for term in terms[1:]:
        rec = go[term]
        parents = rec.get_all_parents()
        parents.update({term})
        
        # Find the intersection with the candidates, and update.
        candidates.intersection_update(parents)
    return candidates

def deepest_common_ancestor(terms, go):
    '''
        This function gets the nearest common ancestor 
        using the above function.
        Only returns single most specific - assumes unique exists.
    '''
    # Take the element at maximum depth. 
    try: 
        dca = max(common_parent_go_ids(terms, go), key=lambda t: go[t].depth)
        return dca
    
    #a good portion of the time it is difficult to find a common parent 
    #partially due to not being in the same ontology 
    except ValueError: 
        pass 
        
def min_branch_length(go_id1, go_id2, go):
    '''
        Finds the minimum branch length between two terms in the GO DAG.
    '''
    try: 
        # First get the deepest common ancestor
        dca = deepest_common_ancestor([go_id1, go_id2], go)
        
        # Then get the distance from the DCA to each term
        dca_depth = go[dca].depth
        d1 = go[go_id1].depth - dca_depth
        d2 = go[go_id2].depth - dca_depth
        
        # Return the total distance - i.e., to the deepest common ancestor and back.
        return d1 + d2
    except: 
        return None
#takes a GO ID and returns its closest id in the reference set 

def getClosestAnnot(inAnnot, referenceAnnots): 
    minDist = 4000 
    closestAnnot = None
    
    #loop through reference annots and keep the closest 
    for refAnnot in referenceAnnots: 
        currDist = min_branch_length(inAnnot, refAnnot, go)
        
        #make sure currDist is not none
        if currDist and currDist < minDist: 
            minDist = currDist
            closestAnnot = refAnnot
    
    #if we can't find any annots related to the current one we can't do the analysis 
    if minDist == 4000: 
        minDist = None
    else: 
        if go[inAnnot].depth < go[closestAnnot].depth: 
            minDist = -minDist
    
    #for sanity checking        
    #print("min distance: " + str(minDist))
    #print("min distance annot: " + str(closestAnnot))
    return minDist

#takes a sequence of GO IDs as the first argument and returns a dictionary with the closest ID in the reference set (second argument) of GO IDs for each one 

def getClosestAnnots(inAnnots, referenceAnnots): 
    closestAnnots = [getClosestAnnot(annot, referenceAnnots) for annot in inAnnots if getClosestAnnot(annot, referenceAnnots)]
    closestAnnots = [item for item in closestAnnots if item != None]
    return closestAnnots

#get the GO IDs from the input file and return a list of them
def parseInFile(inFilePath): 
    inFile = open(inFilePath)
    
    #the index (starts at zero) of the column that contains the GO IDs we are interested in extracting 
    GO_COLUMN = 1
    
    #strip the newline just incase GO is the last column before the next line
    return [line.split("\t")[GO_COLUMN].strip("\n") for line in inFile]

#a folder containing the annotations for the reference matrix (BLOSUM62 in this case)
referenceFolder = sys.argv[1]

#a folder containing the annotations for the new matrix (MC30 in this case)
newFolder = sys.argv[2]

referenceFiles = []
newFiles = []

referenceName = "B62"
newName = "MC30"

#get the names of the files we need in lists 
for fileName in os.listdir(referenceFolder): 
    if ".txt" in fileName and referenceName in fileName: 
        referenceFiles.append(fileName)
    
for fileName in os.listdir(newFolder): 
    if ".txt" in fileName and newName in fileName: 
        newFiles.append(fileName)

closestAnnots = []
#go through each protein and get the distance to the closest ancestor 
for referenceFile in referenceFiles: 
    newFile = [item for item in newFiles if item.replace(newName, referenceName) == referenceFile][0]
    
    #get lists of annotations for this protein 
    refAnnots = parseInFile(referenceFolder.rstrip("/") + "/" + referenceFile)
    newAnnots = parseInFile(newFolder.rstrip("/") + "/" + newFile)
    
    #get only annotations that are in the new set but NOT the reference set 
    newAnnots = set(newAnnots) - set(refAnnots)
    #refAnnots = set(refAnnots) - set(newAnnots)
    print(len(newAnnots))
    print(len(refAnnots))
    
    currClosestAnnots = getClosestAnnots(newAnnots, refAnnots)
    print(len(closestAnnots))
    closestAnnots += currClosestAnnots
    
posClosest = [item for item in closestAnnots if item > 0] 
print("percentage greater than zero" + str(len(posClosest)/len(closestAnnots)))
freqs = plt.hist(closestAnnots, bins=np.arange(min(closestAnnots), max(closestAnnots), 1), color="#5ab4ac")
plt.axvline(0, color="black")
plt.xlim([min(closestAnnots)-1, max(closestAnnots) +1])
plt.ylim([0, max(freqs[0])+1])
plt.fill_between([min(closestAnnots)-1,0], [0],[max(freqs[0])+1], color="#efefefff")  


plt.savefig("INSERT/YOUR/PATH", dpi=1000)
plt.savefig("INSERT/YOUR/PATH", dpi=1000)

