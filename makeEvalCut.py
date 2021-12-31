#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 19:39:06 2020

@author: benlitterer
"""

import sys 
import os 
import shutil 

evalCut = float(sys.argv[1])
inDir = sys.argv[2]
outDir = sys.argv[3]

if os.path.exists(outDir): 
    shutil.rmtree(outDir)
    
#make the new output directory 
os.mkdir(outDir)

#a dict where the keys are filenames and the values are the prop of rows KEPT for that file 
propDict = {}

for fileName in os.listdir(inDir): 
    currFile = open(inDir.rstrip("/") + "/" + fileName)
    currOutFile = open(outDir.rstrip("/") + "/" + fileName, "w")
    
    outStr = ""
    lessThanCount = 0 
    totalCount = 0 
    for line in currFile: 
        if float(line.split("\t")[2]) <= evalCut: 
            outStr += line
            lessThanCount += 1
        totalCount += 1
    if totalCount != 0: 
        propDict[fileName] = lessThanCount/totalCount
    else: 
        propDict[fileName] = 1
        
    currOutFile.write(outStr)
    currFile.close()
    currOutFile.close()

analysisFile = open(outDir.rstrip("/") + "/" + "ANALYSIS.txt", "w")
analysisStr = "\n".join([str(key) + ": " + str(val) for key,val in propDict.items()])
analysisStr += "\nmean: " + str(sum(float(val) for val in propDict.values())/float(len(propDict.keys())))
analysisFile.write(analysisStr)
analysisFile.close()
