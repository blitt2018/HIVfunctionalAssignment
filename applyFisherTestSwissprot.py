#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  6 23:22:56 2020

@author: benlitterer
"""


import sys 
import os 
import glob
import scipy.stats as stats 
from statsmodels.stats.multitest import multipletests 

goAnnotFile = open(sys.argv[1])

#inDir should be the head file of the protein go label hierarchy
inDir = sys.argv[2]

supportDict = {}
totalCount = 0 

for line in goAnnotFile: 
    lineList = line.split("\t")
    
    if len(lineList) == 2: 
        goList = set(lineList[1].strip("\n").replace(" ","").split(";"))
        for go in goList: 
            if go not in supportDict: 
                supportDict[go] = 1 
            else: 
                supportDict[go] += 1
            totalCount += 1
            
#a nested dictionary where the first keys are files, and the vals are dicts with the supporting prots for each go in that file
sampleSupportDict = {}

#count up all of the unique proteins that support a particular go for each file in the prot hier
#(repeat blast output results do not add to the counts)
for outFileName in glob.glob(inDir.rstrip("/") + "/*"):
    if "ANALYSIS" not in outFileName: 
        outFile = open(outFileName)
        
        #a dict to hold the support for gos in this file
        goSupport = {}
        
        #the total amount of unique output proteins in this file
        outTotalCount = 0
        
        for line in outFile: 
            lineList = line.rstrip("\n").split("\t")
            prot = lineList[0]
            if len(lineList) > 1: 
                goList = lineList[1].split(";")
                
                for go in goList: 
                    if go not in goSupport: 
                        outTotalCount += 1
                        goSupport[go] = {prot}
                    else: 
                        #if this prot is new then add to the total amount of annotations 
                        if prot not in goSupport[go]: 
                            outTotalCount +=1
                        goSupport[go].add(prot)
            #connect this file to this support dict 
        sampleSupportDict[outFileName] = [goSupport, outTotalCount]
               

#key: file, val: list of gos with their exact test results + more info 
outputDict = {}

for fileName, infoList in sampleSupportDict.items(): 
    
    #a list of gos + their test results etc... for each go in this file
    outList = []
    
    #dict with support for gos in this file
    goSupportDict = infoList[0]
    
    #integer with the amount of total support proteins in this file
    outTotalCount = infoList[1]
    
    #for go, support set in dict
    for go, protSet in goSupportDict.items(): 
        
        if go in supportDict: 
            oddsratio, pvalue = stats.fisher_exact([[len(protSet),supportDict[go]],[outTotalCount, totalCount]])
            outList.append([outFileName, go, len(protSet), outTotalCount, supportDict[go], totalCount, oddsratio, pvalue])
            
    outputDict[fileName] = outList

for fileName, outList in outputDict.items():
    #print("####### " + fileName + " #############################################")
    if len(outList) > 0: 
        bfout = multipletests([stats[7] for stats in outList], method="bonferroni")
        correctedPvals = list(bfout[1])
        print(correctedPvals[:10])
        
        correctedOutList = []
        
        for i in range(len(outList)): 
            currList = outList[i]
            pval = correctedPvals[i]
            
            #if this is a significant go and it as an odds ratio of greater than one
            if pval < .05 and currList[6] > 1: 
                newList = currList + [correctedPvals[i]]
                correctedOutList.append(newList)
        
        correctedOutList.sort(key = lambda subList: subList[-1])
        
        outFile = open(fileName.rstrip("/") + "GOLIST.txt", "w")
        outFile.write("\n".join(["\t".join(map(str,line)) for line in correctedOutList]))


