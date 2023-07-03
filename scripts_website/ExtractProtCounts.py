#!/usr/bin/env python3

import os
import sys

def ExtractProtCounts(filename):
    """Extract the total number of proteins within a given phylogenetic group's proteomes"""

    #Open file
    try:
        infile = open('../generated_data/data_website/SP_types_counts_phylogenetic_groups/'+filename, "r")

            #Protein count will always 2nd tab-separated element on 2nd line
        infile.readline()
        
        #2nd line
        data_line = infile.readline()

        #Extract protein- and proteome count
        prot_count = data_line.split("\t")[1]
        proteome_count = data_line.split("\t")[2]

        infile.close()
    
    except IOError as err:
        #print(filename)
        #print("Error: ", str(err))
        prot_count = 0
        proteome_count = 0
        #sys.exit(1)

    return prot_count, proteome_count