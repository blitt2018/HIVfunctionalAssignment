#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 21:30:45 2020

@author: benlitterer
"""

import sys
import os 
import shutil 

#this is the file containing cross references from nr to uniprot
refFile = sys.argv[1]

#this is the directory that our (raw) blast output is stored in 
originalOutputDir = sys.argv[2]

#this is the new directory that will hold the same output but with cross referenced ids instead of nr (if we have them)
newOutputDir = sys.argv[3]

def updateOutputDict(crossRefDict, origDir, outputDict): 
    
    #for each file in inDir
    for fileName in os.listdir(origDir): 
        if "ANALYTICS" not in fileName: 
            #account for fact that the inDir could have trialing slash or not
            filePath = origDir.rstrip("/") + "/" + fileName
            file = open(filePath, "r")
            
            fileLines = file.readlines()
            
            #for testing mostly 
            for i in range(1,len(fileLines)):
                
                #line in the original output file
                line = fileLines[i]
                
                lineList = line.strip("\n").split("\t")
                
                if len(lineList) >= 2: 
                    #the uniprot id is in the second column of the blast output 
                    currID = lineList[1].split(".")[0]
                    
                    if currID in crossRefDict: 
                        goSet = crossRefDict[currID]
                        
                        if currID not in outputDict[fileName]: 
                            outputDict[fileName][currID] = goSet
                        else: 
                            outputDict[fileName][currID] = outputDict[fileName][currID].union(goSet)
            file.close()
    return outputDict 
        
#initialize dict that will hold files as keys and a uniID:goSet dictionary as values   
outputDict = {filename:{} for filename in os.listdir(originalOutputDir) if "ANALYTICS" not in filename}

crossRefDict = {}
for line in open(refFile): 
    lineList = line.split("\t")
    if len(lineList) == 2: 
        uniAcc = lineList[0]
        goTerm = lineList[1].strip("\n").replace(" ","").split(";")
        
        """
        #possible exceptions or weird things that could happen with formatting? 
        if "." in uniAcc: 
            wOutDot = uniAcc.split(".")[0]
            if wOutDot not in  crossRefDict: 
                crossRefDict[wOutDot] = set([goTerm])
            else: 
                crossRefDict[wOutDot].add(goTerm)
            
        if "_" in uniAcc: 
            wOutUS = uniAcc.split("_")[1]
            if wOutUS not in  crossRefDict: 
                crossRefDict[wOutUS] = set([goTerm])
            else: 
                crossRefDict[wOutUS].add(goTerm)
        """
        
        if uniAcc not in  crossRefDict: 
            crossRefDict[uniAcc] = set(goTerm)
            
outputDict = updateOutputDict(crossRefDict, originalOutputDir, outputDict) 


#remove directory if it exists, even if it has files 
if os.path.exists(newOutputDir): 
    shutil.rmtree(newOutputDir)
    
#make the new output directory 
os.mkdir(newOutputDir)

#write each inner dictionary to the file corresponding to its key 
for key, value in outputDict.items(): 
    outfile = open(newOutputDir.rstrip("/") + "/" +  key.split(".")[0] + "UNI.tsv", "w")
    
    outStr = ""
    #inner dictionary with uniID:{GO1, GO2, GO3} format 
    for innerKey, innerVal in value.items(): 
        outStr += innerKey + "\t" + ";".join(list(innerVal)) + "\n"
    outfile.write(outStr)
  


outAnalysisStr = ""
for key, value in outputDict.items():
    outfile = open(newOutputDir.rstrip("/") + "/" +  key.split(".")[0] + "UNI.tsv")
    
    totalLen = len(outfile.readlines())
    newLen = len(list(value.keys()))
    
    propLabeled = 0 
    if totalLen != 0: 
        propLabeled = newLen/float(totalLen)
    outAnalysisStr += key + ":\t" + str(propLabeled) + "\n"
    
analyticsOutfile = open(newOutputDir.rstrip("/") + "/" + "ANALYTICS.txt", "w")
analyticsOutfile.write(outAnalysisStr)

