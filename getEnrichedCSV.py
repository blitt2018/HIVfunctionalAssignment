#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 21:50:46 2020

@author: benlitterer
"""

import sys 
import glob 
import pandas as pd 

inDirName = sys.argv[1]
outFileName = sys.argv[2]

lenDict = {}

for fileName in glob.glob(inDirName.rstrip("/") + "/*"):
    if "ANALYTICS" not in fileName: 
        wOutMatrix = fileName.replace("B62", "").replace("MC30","")
        protName = wOutMatrix.split("/")[-1].split("Blast")[0]
        
        if protName not in lenDict: 
            lenDict[protName] = [0, 0]
            
        if "B62" in fileName: 
            lenDict[protName][0] = len(open(fileName).readlines())
        
        if "MC30" in fileName: 
            lenDict[protName][1] = len(open(fileName).readlines())

#reshape lenDict
lenDict = {"prot":[key for key in lenDict.keys()],"B62":[val[0] for val in lenDict.values()],"MC30":[val[1] for val in lenDict.values()]}

df = pd.DataFrame(data=lenDict)
df.to_csv(outFileName)

