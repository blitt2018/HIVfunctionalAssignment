#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 09:13:43 2020

@author: benlitterer
"""

#NOTE: this is mean to be run in parallel where we list all the enrichment files in a directory and pipe that into the gnu parallel command 
#we then provide one output directory, and each input file is given a modified name and the corresponding output file is put in that directory 

#Example of output is in GOGOInputs001

#this file takes in the enrichment file tsv input (one file for each protein, each annotation has a row with stats on the enrichment for that protein)
#we then output a file that has every pair of GO annotations contained in the enrichment file 
#we assume that every GO annotation in the enrichment file passes the .05 corrected p-value cutoff. If not, the file needs to be filtered 

from goatools import obo_parser
import sys

#the folder containing the obo file from the gene ontology website representing the GO hierarchy 
go_obo = "/home/benlitterer/Academic/Research/ProjectHIV"

#tried the standard and basic versions... doesn't appear to be much of a difference
go = obo_parser.GODag(go_obo + "/go-basic.obo")

#the enriched file (should be in a folder containing many of these files)
enrichedFile = open(sys.argv[1])

#get GOs from the file 
goList = [line.split("\t")[1] for line in enrichedFile if len(line.split("\t")) > 1]
outText = ""

#create pairs of GO annots and output them 
for go1 in goList: 
    for go2 in goList: 
        
        #NOTE: we may lose some annotations here if they aren't in the go_obo file (nearly all are)
        if go1 in go and go2 in go: 
            
            #namespace just means which one of the three hierarchies (Biological Process, Molecular Function, Cellular Component)
            if go[go1].namespace == go[go2].namespace: 
                outText += go1 + " " + go2 + "\n"

        #just to manually count the number of missing pairs 
        else: 
            print("no")

outDirName = sys.argv[2]

#NOTE: we are automatically creating an output file name here... could be changed if needed 
outFile = open(outDirName.rstrip("/") + "/" + sys.argv[1].split("/")[-1].split("Blast")[0] + "GOGOInput.txt", "w")
outFile.write(outText)