#!/usr/bin/env python3

import os
import sys
import statistics


def CountSPs(filename):
    """Count region lengths for the different SP types within a phylogenetic group"""

    #Open file with SP region length counts
    try:
        infile = open("../generated_data/data_website/SP_regions_counts_phylogenetic_groups/" + filename, "r")
    except IOError as err:
        print("Error: ", str(err))
        sys.exit(1)

    #Initialize
    n_region_sp = []
    h_region_sp = []
    c_region_sp = []
    n_region_lipo = []
    h_region_lipo = []
    n_region_tat = []
    h_region_tat = []
    c_region_tat = []
    n_region_tatlipo = []
    h_region_tatlipo = []
    n_region_pilin = []

    #Control check; remove wrongly predicted signal peptides
    count_wrong = 0
    count_right = 0

    infile.readline()

    #Go through all SPs, one at a time
    for line in infile:
        #Collect Sec SPI region lengths
        if line.split()[-2] == "SP":
            #Filter out wrongly predicted SPs from analyses
            if line.split()[3] != "0" and line.split()[2] != "0" and line.split()[1] != "0":
                #Append region length observations to separate lists
                n_region_sp.append(int(line.split()[1]))
                h_region_sp.append(int(line.split()[2]))
                c_region_sp.append(int(line.split()[3]))
                count_right += 1
        
        #Collect Sec SPII region lengths
        elif line.split()[-2] == "LIPO":
            #Filter out wrongly predicted SPs from analyses
            if line.split()[3] == "0" and line.split()[2] != "0" and line.split()[1] != "0":
                #Append region length observations to separate lists
                n_region_lipo.append(int(line.split()[1]))
                h_region_lipo.append(int(line.split()[2]))
                count_right += 1
        
        #Collect Tat SPI region lengths
        elif line.split()[-2] == "TAT":
            #Filter out wrongly predicted SPs from analyses
            if line.split()[3] != "0" and line.split()[2] != "0" and line.split()[1] != "0":
                #Append region length observations to separate lists
                n_region_tat.append(int(line.split()[1]))
                h_region_tat.append(int(line.split()[2]))
                c_region_tat.append(int(line.split()[3]))
                count_right += 1
        
        #Collect Tat SPII region lengths
        elif line.split()[-2] == "TATLIPO":
            #Filter out wrongly predicted SPs from analyses
            if line.split()[3] == "0" and line.split()[2] != "0" and line.split()[1] != "0":
                #Append region length observations to separate lists
                n_region_tatlipo.append(int(line.split()[1]))
                h_region_tatlipo.append(int(line.split()[2]))
                count_right += 1

        #Collect Sec SPIII region lengths
        elif line.split()[-2] == "PILIN":
            #Filter out wrongly predicted SPs from analyses
            if line.split()[3] == "0" and line.split()[2] == "0":
                #Append region length observations to separate lists
                n_region_pilin.append(int(line.split()[1]))
                count_right += 1

    infile.close()

    #Return all region counts
    return n_region_sp, h_region_sp, c_region_sp, n_region_lipo, h_region_lipo, n_region_tat, h_region_tat, c_region_tat, n_region_tatlipo, h_region_tatlipo, n_region_pilin

