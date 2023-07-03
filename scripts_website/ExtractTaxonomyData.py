#!/usr/bin/env python3

import os
import sys
from ExtractProtCounts import *

def ExtractTaxonomyData():
    """Extract taxonomical information about each group"""
    """Returns a list of all Tax IDs with corresponding Tax names"""
    """Returns a dictionary that links TAX ID to names"""
    
    #Save all taxonomic ranks in list
    files_tax_ranks = os.listdir("../generated_data/data_website/TaxonomyData")

    #Initialize 
    tax_data = dict()
    all_taxa = list()
    all_tax_ids = list()

    #Append empty object to initialize webpage without a pre-selected group
    all_taxa.append("")


    
    for i in range(len(files_tax_ranks)):
        
        #Open file with information about taxonomical 
        #groups at the given rank
        infile = open("../generated_data/data_website/TaxonomyData/" + files_tax_ranks[i], "r")
        
        #Each line in file corresponds to the tax ID and 
        #name of a taxonomical group on that rank
        for line in infile:
            taxa = line.split("\t")             #Split line into Tax ID and name(s)
            try:
                #Extract proteome count of particular group
                proteome_count = ExtractProtCounts(taxa[0]+"_SP_types_counts.tab")[1]

                #Filter away groups without proteome data
                if int(proteome_count) > 0:
                    
                    #Add information to dict about taxonomical groups
                    #This adds; tax ID as key, 
                                #list of the tax group scientific name and tax rank as value
                    #taxa[0] is tax id, files_tax_ranks[i][:-4] is rank
                    tax_data[taxa[0]] = [files_tax_ranks[i][:-4]]
                    #Add several names, if the group has more than one
                    if len(taxa) > 2:
                        for j in range(len(taxa)-2):
                            tax_data[taxa[0]].append(taxa[j+1])
                    #Add name
                    elif len(taxa) == 2:
                        tax_data[taxa[0]].append(taxa[-1][:-1])

                    #Append all tax IDs to list
                    all_taxa.append(taxa[0])
                    all_tax_ids.append(taxa[0])

                    #Append all tax group names to list, except the ones without rank
                    if files_tax_ranks[i] != "no rank.txt":
                        if len(taxa) > 2:
                            for j in range(len(taxa)-2):
                                all_taxa.append(taxa[j+1])

                        elif len(taxa) == 2:
                            all_taxa.append(taxa[-1][:-1])
            except IOError as err:
                continue
                

        infile.close()

    all_taxa = tuple(all_taxa)

    
    #all_taxa is a tuple of all possible taxonomic IDs and names
    #tax_data is a dictionary linking tax IDs to tax name(s) for each group
    return all_taxa, tax_data


ExtractTaxonomyData()