#!/usr/bin/env python3

import streamlit as st
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import math
import sys
import pickle
import logomaker
from PIL import Image
import gzip
import io
import statistics

##Functions required to run app
from CountSPs import *
from MakeHistograms import *
from ExtractProtCounts import *
from VisualizeProtCounts import *
from VisualizeSPCounts import *
from VisualizeSPFrequencies import *
from TreeSearchDict import *
from ExtractTaxonomyData import *
from MakeLogos import *
from ExtractEukaryoticGroups import *
from VisualizeGCContent import *

#Use this when implementing download of SPs
#from DownloadSPs import *


#Remember to be in main directory with all subdirectories

#######################Collect data from functions###############################

#Get dict of taxonomy data for tree-structured search
tree_search_file_subgroups = TreeSearchDict("all_domains")[0]
tree_search_file_parents = TreeSearchDict("all_domains")[1]

#Extract list of all taxa for search field and dict og taxonomy data connections
all_taxa = ExtractTaxonomyData()[0]
tax_data = ExtractTaxonomyData()[1]

#Insert images
seq_logo_color_schemes = Image.open('../website_additions/seq_logo_color_schemes.PNG')
dtu_logo = Image.open('../website_additions/DTU_logo.png')

#Extract list of eukaryotic taxonomies
euk_tax_groups = ExtractEukaryoticGroups()


#############################################################Streamlit design############################################################################

st.title("SPEx (Signal Peptide Explorer)")

st.subheader("Signal Peptides Across the Tree of Life")

option = st.selectbox(
     "Select the phylogenetic group you would like to view analyses for",
     all_taxa)
with st.expander("Information about data and distribution of phylogenetic groups"):
     st.write("""
         **The data**

         The analyses presented on this website are based on data from 323 archaeal organisms, 
         7920 bacterial organisms and 1527 eukaryotic organisms, giving a total of 9770 organisms
         spanning a broad phylogenetic spectrum. All organisms included in the dataset have
         been selected based on having reference proteomes as defined by 
         [UniProt](https://www.uniprot.org/proteomes/?query=*&fil=reference%3Ayes), from which the 
         proteomes have also been downloaded. All signal peptides have been 
         predicted from the proteomes using the signal peptide prediction program 
         [SignalP 6.0](https://services.healthtech.dtu.dk/service.php?SignalP-6.0), which the analyses
         presented on this website are based on. 


         **Distribution of phylogenetic groups**

         Signal peptide analyses can be viewed for phylogenetic groups
         on any taxonomical level. The taxonomical distribution of species is based on
         [NCBI's Taxonomy Database](https://www.ncbi.nlm.nih.gov/taxonomy).

         **Searching for phylogenetic groups**

         In the search box above, you can search for a phylogenetic group or organism by either searching for a name or for the 
         taxonomical ID defining that particular group.
     """)

#Extract taxonomical information about the group selected using selectbox
#Extract Tax ID, scientific name of group and rank
for key, value in tax_data.items():
    if option in value[1:]:
        name = option
        tax_id = key
        clean_rank = value[0]
        
        #Define scientific name for tree search
        for tax_name in value:
            if tax_name in tree_search_file_subgroups.keys():
                scient_name = tax_name
    
    elif option == key:
        tax_id = option
        name = value[1]
        clean_rank = value[0]
        
        #Define scientific name for tree search
        for tax_name in value:
            if tax_name in tree_search_file_subgroups.keys():
                scient_name = tax_name

#Initialize 
list_check = []
counter_par = 0
counter = 0

#Run this only if input to selectbox has been given
if option != "":
    #Make Tree-like structure in sidebar if possible
    #Make dict where each Tax ID is key with its' children as values for Tree-like structure in sidebar
    try:
        st.sidebar.header("Tree-structured search for subgroups")
        
        #Make checkbox linking to parent group of a given group

        #Placeholder used to modify checkbox once clicked
        placeholder = st.sidebar.empty()
        counter_par += 1
        
        #If we are at domain-level, no parent group can be displayed
        if tree_search_file_parents[tax_id] == []:
            placeholder.write("")
            run_again = False
        else:
        #Create checkbox - once clicked move to parent group
            button_ = placeholder.checkbox("Go to parent group: " + tree_search_file_parents[tax_id][2], key = counter_par)
            if button_:
                clean_rank = tree_search_file_parents[tax_id][0]
                scient_name = tree_search_file_parents[tax_id][2]
                tax_id = tree_search_file_parents[tax_id][1]

                #Initialize
                run_again = True

                while run_again:
                    #Re-initialize
                    run_again = False
                    counter_par += 1

                    #Keep looping to next parent group until domain level is reached
                    if tree_search_file_parents[tax_id] == []:
                        placeholder.write("")
                        run_again = False
                    else:
                        button_ = placeholder.checkbox("Go to parent group: " + tree_search_file_parents[tax_id][2], key = counter_par)
                        if button_:
                            clean_rank = tree_search_file_parents[tax_id][0]
                            scient_name = tree_search_file_parents[tax_id][2]
                            tax_id = tree_search_file_parents[tax_id][1]
                            run_again = True

        #Make radio to display subgroups that one can go into
        #Display "name (rank)" in radio - first extract this
        children_names = tree_search_file_subgroups[scient_name][::2]
        children_ranks = tree_search_file_subgroups[scient_name][1::2]

        children = []

        #Children contains list of "name (rank)" options for subgroups
        for i in range(len(children_names)):
            child_name = children_names[i] + " "+children_ranks[i]
            children.append(child_name)

        #Initialize first choice to be the pre-selected group itself
        children[0] = children[0] + "("+ scient_name + ")"

        #Sort groups in sidebar to be shown alphabetically, 
        #except the pre-selected group, which is positioned first
        pre_group = [children[0]]
        sub_groups = children[1:]
        sub_groups.sort(key=lambda s: s.lower())

        children_sorted = pre_group + sub_groups

        #Stop when the "lowest" branch is reached
        if len(children) == 2 and "(strain)" in children[1]:
            stop = True
        #Keep showing subgroups for the selected group
        else:
            counter += 1
            option = st.sidebar.radio(
            "Click on one of the subgroups belonging to "+ scient_name + ", if you would like to look into it.", children_sorted,
            key = counter)
            with st.sidebar.expander("Further information"):
                st.write("""
                Here you can dig further into one of the subgroups belonging to """ + scient_name + """. When clicking 
                on one of these groups, the interface will change and you will be presented with analyses of this subgroup.
                You can continue exploring subgroups of subgroups, until you reach the lowest level. 
             """)

            stop = False

        #Re-define taxonomical information for newly selected group
        while option != children[0] and not stop:
            option = option.split(" (")[0]

            for key, value in tax_data.items():
                if option in value[1:]:
                    name = option
                    tax_id = key
                    clean_rank = value[0]
                    
                    #Define scientific name for tree search
                    for tax_name in value:
                        if tax_name in tree_search_file_subgroups.keys():
                            scient_name = tax_name
                
                elif option == key:
                    tax_id = option
                    name = value[1]
                    clean_rank = value[0]
                    
                    #Define scientific name for tree search
                    for tax_name in value:
                        if tax_name in tree_search_file_subgroups.keys():
                            scient_name = tax_name
            
            children_names = tree_search_file_subgroups[scient_name][::2]
            children_ranks = tree_search_file_subgroups[scient_name][1::2]

            children = []
            for i in range(len(children_names)):
                child_name = children_names[i] + " "+children_ranks[i]
                children.append(child_name)
            if len(children) == 1: 
                stop = True
            if len(children) == 2 and "(strain)" in children[1]:
                stop = True
            else:
                children[0] = children[0] + " ("+ scient_name + ")"

                #sort groups in sidebar to come alphabetically, except the pre-selected group, which is positioned first
                pre_group = [children[0]]
                sub_groups = children[1:]
                sub_groups.sort(key=lambda s: s.lower())

                children_sorted = pre_group + sub_groups

                counter += 1
                st.sidebar.subheader("Currently viewing "+ scient_name)
                option = st.sidebar.radio("Please click on a subgroup, to see the analyses of this group.", children_sorted,
                key = counter)

                #Split string to extract name, rank and ID of selected group
                option_redo = option.split(" (")[0]

                stop = True
                for key, value in tax_data.items():
                    if option_redo in value[1:]:
                        stop = False
                        name = option_redo                          #Save name (selected) of group
                        tax_id = key                                #Save tax ID of the group
                        clean_rank = value[0]                       #Save rank of the group
                        #Define scientific name for tree search
                        for tax_name in value:
                            if tax_name in tree_search_file_subgroups.keys():
                                scient_name = tax_name              #Save scientific name
                    elif option_redo == key:
                        stop = False
                        tax_id = option_redo                        #Save tax ID of the group
                        name = value[1]                             #Save name (selected) of group
                        clean_rank = value[0]                       #Save rank of the group
                        #Define scientific name for tree search
                        for tax_name in value:
                            if tax_name in tree_search_file_subgroups.keys():
                                scient_name = tax_name              #Save scientific name


        #Display analyses of group
        st.header('You are viewing analyses for the ' + clean_rank + ' '+ scient_name + ' ('+ tax_id + ')')
        
        st.subheader("Proteome and protein counts")
        VisualizeProtCounts(tax_id)

        st.subheader("Counts of each signal peptide type")
        VisualizeSPCounts(tax_id)

        st.subheader("Fraction of proteins tagged with each signal peptide type")
        VisualizeSPFrequencies(tax_id)

        st.subheader("GC-content")
        VisualizeGCContent(tax_id)
        with st.expander("Information about data and calculation of GC-content"):
            st.write("""
                Genomes have been downloaded from the 
                [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/home) in
                EMBL format. The GC-content within each group have been calculated
                based on this genomic data. For some of the organisms, it was not 
                possible to download the genomes and thereby not possible to calculate GC-content.
            """)

        #Specialize design if group belongs to eukarya (only Sec SPI)
        if tax_id in euk_tax_groups:
            
            with st.expander("Description of the signal peptide types"):
             st.write("""
                 **Sec SPI**: Signal peptides transported by the Sec translocon and cleaved by Signal Peptidase I.
                 This is the most frequent type, and the only one present in eukaryotes.

                 **Sec SPII**: Lipoprotein signal peptides transported by the Sec translocon and cleaved by Signal Peptidase II.

                 **Tat SPI**: Tat signal peptides transported by the Tat translocon and cleaved by Signal Peptidase I.

                 **Tat SPII**: Tat lipoprotein signal peptides transported by the Tat translocon and cleaved by Signal Peptidase II.

                 **Sec SPIII**: Pilin and pilin-like signal peptides transported by the Sec translocon and cleaved by Signal Peptidase III.

                 Regarding eukaryotic groups and organisms, it should be noted that the Tat translocon is present in chloroplasts and some mitochondria, 
                 however this is not captured on this website, and only Sec SPI analyses will be shown for eukaryotic groups and organisms.
             """)

            st.subheader("Showing analyses for Sec SPI")
        
            #Sec SPI analyses
            select_analyses = st.radio("Which analyses would you like to take a look at?", ("Full signal peptide length distribution", "Length distributions of each signal peptide region", "Sequence logos", "Download protein sequences tagged with signal peptide"))
            if select_analyses == "Full signal peptide length distribution":
                st.subheader("Histogram of signal peptide length distribution")
                MakeHistograms(tax_id+"_SP_regions_counts.tab", "SP", group = scient_name, see_regions = False, rank = clean_rank)
            elif select_analyses == "Length distributions of each signal peptide region":
                st.subheader("Histograms of region length distributions")
                MakeHistograms(tax_id+"_SP_regions_counts.tab", "SP", group = scient_name, see_regions = True, rank = clean_rank)
            elif select_analyses == "Sequence logos":
                st.subheader("Sequence logos aligned at region borders")
                
                #Create menu for dynamic settings (color and position definement)
                option1,option2 = st.columns(2)
                option1_inf,option2_inf = st.columns(2)
                color_scheme = option1.radio("Select a color scheme for the sequence logos", ("Default", "Hydrophobicity", "Charge", "Chemistry"))
                with option1_inf.expander("Information about color schemes"):
                    st.image(seq_logo_color_schemes)
                    st.write("""
                        Default color scheme divides amino acids into negatively charged (red), 
                        positively charged (blue), hydrophobic (black) and others (green). This color scheme
                        is useful when examining important signal peptides properties such as
                        the positive charge of n-regions and hydrophobicity of h-regions. 
                        """)

                keep_zero = option2.radio("Do you want to include 0 as the position before alignment site?", ("Yes", "No"))
                with option2_inf.expander("Information about position assignments"):
                    st.write("""
                        Choosing either of these options do not influence the logos themselves, 
                        except changing whether the position before alignment site should be 
                        written as position 0 or position -1. 
                        """)
                
                #Define logo positions to show
                min_pos_n = ExtractAlignmentPosition(tax_id, "SP", "N")[0]
                len_seq_n = ExtractAlignmentPosition(tax_id, "SP", "N")[1]
                max_pos_n = int(len_seq_n) - int(min_pos_n)
                st.write("**Select the range of positions in the logo aligned at the n-h region border**")
                low,high = st.columns(2)

                #Potential error-handling for short sequences
                if int(min_pos_n) >= 5:
                    lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -5, step = 1)
                else:
                    lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -int(min_pos_n), step = 1)
                
                if int(max_pos_n) >= 7:
                    upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = 7, step = 1)
                else: 
                    upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = max_pos_n, step = 1)
                
                #Display logo
                MakeLogo(tax_id,"SP","n", lower_n, upper_n, scient_name, color_scheme, keep_zero)
                
                min_pos_h = ExtractAlignmentPosition(tax_id, "SP", "H")[0]
                len_seq_h = ExtractAlignmentPosition(tax_id, "SP", "H")[1]
                max_pos_h = int(len_seq_h) - int(min_pos_h)
                st.write("**Select the range of positions in the logo aligned at the h-c region border**")
                low,high = st.columns(2)
                
                if int(min_pos_h) >= 15:
                    lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -15, step = 1)
                else:
                    lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -int(min_pos_h), step = 1)
                
                if int(max_pos_h) >= 5:
                    upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = 5, step = 1)
                else:
                    upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = max_pos_h, step = 1)

                MakeLogo(tax_id,"SP","h", lower_h, upper_h, scient_name, color_scheme, keep_zero)

                min_pos_c = ExtractAlignmentPosition(tax_id, "SP", "C")[0]
                len_seq_c = ExtractAlignmentPosition(tax_id, "SP", "C")[1]
                max_pos_c = int(len_seq_c) - int(min_pos_c)
                st.write("**Select the range of positions in the logo aligned at cleavage site**")
                low,high = st.columns(2)
                if int(min_pos_c) >= 7:
                    lower_c = low.number_input('Lower value', min_value = -int(min_pos_c), max_value = -1, value = -7, step = 1)
                else:
                    lower_c = low.number_input('Lower value', min_value = -int(min_pos_c), max_value = -1, value = -int(min_pos_c), step = 1)
                
                if int(max_pos_c) >= 5:
                    upper_c = high.number_input('Upper value', min_value = 1, max_value = max_pos_c, value = 5, step = 1)
                else:
                    upper_c = high.number_input('Upper value', min_value = 1, max_value = max_pos_c, value = max_pos_c, step = 1)

                MakeLogo(tax_id,"SP","c", lower_c, upper_c, scient_name, color_scheme, keep_zero)

                #Explanations about sequence logo creation
                with st.expander("Creation of sequence logos"):
                    st.write("""
                        The sequence logos presented on this website are based on Shannon Entropy,
                        as described in "Sequence logos: a new way to display consensus sequences" 
                        by T. D. Schneider and R. M. Stephens in 1990.
                        """)
                with st.expander("Calculation of errorbars"):
                    st.write("""
                        The magnitude of the errorbars are based on the small sample correction
                        factor. The small sample correction factor for sequence logos with 
                        amino acids is calculated as
                        """)
                    st.latex(r'''e_n = \frac{1}{\ln{2}} \cdot \frac{20-1}{2\cdot n},''')
                    st.write("""
                        with *n* being the number of amino acids at a particular aligned position. 
                        The height of the errorbars correspond to twice the small sample correction
                        factor. 
                        """)
            #elif select_analyses == "Download protein sequences tagged with signal peptide":
            #    st.download_button('Download protein sequences tagged with Sec SPI', prepareDownloadSPs(group_number = str(tax_id), sp_type = "Sec SPI"))

        #Prokaryotic groups
        else:
            select_sp = st.radio("Please select the signal peptide type you wish to explore", ('Sec SPI', "Sec SPII", "Tat SPI", "Tat SPII", "Sec SPIII"))
            with st.expander("Description of the signal peptide types"):
             st.write("""
                 **Sec SPI**: Signal peptides transported by the Sec translocon and cleaved by Signal Peptidase I.
                 This is the most frequent type, and the only one present in eukaryotes.

                 **Sec SPII**: Lipoprotein signal peptides transported by the Sec translocon and cleaved by Signal Peptidase II.

                 **Tat SPI**: Tat signal peptides transported by the Tat translocon and cleaved by Signal Peptidase I.

                 **Tat SPII**: Tat lipoprotein signal peptides transported by the Tat translocon and cleaved by Signal Peptidase II.

                 **Sec SPIII**: Pilin and pilin-like signal peptides transported by the Sec translocon and cleaved by Signal Peptidase III.

                 Regarding eukaryotic groups and organisms, it should be noted that the Tat translocon is present in chloroplasts and some mitochondria, 
                 however this is not captured on this website, and only Sec SPI analyses will be shown for eukaryotic groups and organisms.
             """)

            st.subheader("You have chosen to view analyses of " + select_sp)
            
            #Sec SPI analyses
            if select_sp == "Sec SPI":
                select_analyses = st.radio("Which analyses would you like to take a look at?", ("Full signal peptide length distribution", "Length distributions of each signal peptide region", "Sequence logos", "Download protein sequences tagged with signal peptide"))
                if select_analyses == "Full signal peptide length distribution":
                    st.subheader("Histogram of signal peptide length distribution")
                    MakeHistograms(tax_id+"_SP_regions_counts.tab", "SP", group = scient_name, see_regions = False, rank = clean_rank)
                elif select_analyses == "Length distributions of each signal peptide region":
                    st.subheader("Histograms of region length distributions")
                    MakeHistograms(tax_id+"_SP_regions_counts.tab", "SP", group = scient_name, see_regions = True, rank = clean_rank)
                elif select_analyses == "Sequence logos":
                    st.subheader("Sequence logos aligned at region borders")
                    
                    option1,option2 = st.columns(2)
                    option1_inf,option2_inf = st.columns(2)
                    color_scheme = option1.radio("Select a color scheme for the sequence logos", ("Default", "Hydrophobicity", "Charge", "Chemistry"))
                    with option1_inf.expander("Information about color schemes"):
                        st.image(seq_logo_color_schemes)
                        st.write("""
                            Default color scheme divides amino acids into negatively charged (red), 
                            positively charged (blue), hydrophobic (black) and others (green). This color scheme
                            is useful when examining important signal peptides properties such as
                            the positive charge of n-regions and hydrophobicity of h-regions. 
                            """)

                    keep_zero = option2.radio("Do you want to include 0 as the position before alignment site?", ("Yes", "No"))
                    with option2_inf.expander("Information about position assignments"):
                        st.write("""
                            Choosing either of these options do not influence the logos themselves, 
                            except changing whether the position before alignment site should be 
                            written as position 0 or position -1. 
                            """)
                    #st.caption("Alignment at N-H region border")
                    min_pos_n = ExtractAlignmentPosition(tax_id, "SP", "N")[0]
                    len_seq_n = ExtractAlignmentPosition(tax_id, "SP", "N")[1]
                    max_pos_n = int(len_seq_n) - int(min_pos_n)
                    st.write("**Select the range of positions in the logo aligned at the n-h region border**")
                    low,high = st.columns(2)
                    

                    if int(min_pos_n) >= 5:
                        lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -5, step = 1)
                    else:
                        lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -int(min_pos_n), step = 1)
                
                    if int(max_pos_n) >= 7:
                        upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = 7, step = 1)
                    else: 
                        upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = max_pos_n, step = 1)
                    
                    MakeLogo(tax_id,"SP","n", lower_n, upper_n, scient_name, color_scheme, keep_zero)
                    #except ValueError:
                    #    st.warning("Your selected value for lower part of the range is too small. Please try with a greater value.")
                    min_pos_h = ExtractAlignmentPosition(tax_id, "SP", "H")[0]
                    len_seq_h = ExtractAlignmentPosition(tax_id, "SP", "H")[1]
                    max_pos_h = int(len_seq_h) - int(min_pos_h)
                    st.write("**Select the range of positions in the logo aligned at the h-c region border**")
                    low,high = st.columns(2)
                    if int(min_pos_h) >= 15:
                        lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -15, step = 1)
                    else:
                        lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -int(min_pos_h), step = 1)
                
                    if int(max_pos_h) >= 5:
                        upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = 5, step = 1)
                    else:
                        upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = max_pos_h, step = 1)

                    MakeLogo(tax_id,"SP","h", lower_h, upper_h, scient_name, color_scheme, keep_zero)

                    min_pos_c = ExtractAlignmentPosition(tax_id, "SP", "C")[0]
                    len_seq_c = ExtractAlignmentPosition(tax_id, "SP", "C")[1]
                    max_pos_c = int(len_seq_c) - int(min_pos_c)
                    st.write("**Select the range of positions in the logo aligned at cleavage site**")
                    low,high = st.columns(2)
                    

                    if int(min_pos_c) >= 7:
                        lower_c = low.number_input('Lower value', min_value = -int(min_pos_c), max_value = -1, value = -7, step = 1)
                    else:
                        lower_c = low.number_input('Lower value', min_value = -int(min_pos_c), max_value = -1, value = -int(min_pos_c), step = 1)
                
                    if int(max_pos_c) >= 5:
                        upper_c = high.number_input('Upper value', min_value = 1, max_value = max_pos_c, value = 5, step = 1)
                    else:
                        upper_c = high.number_input('Upper value', min_value = 1, max_value = max_pos_c, value = max_pos_c, step = 1)

                    MakeLogo(tax_id,"SP","c", lower_c, upper_c, scient_name, color_scheme, keep_zero)

                    with st.expander("Creation of sequence logos"):
                        st.write("""
                            The sequence logos presented on this website are based on Shannon Entropy,
                            as described in "Sequence logos: a new way to display consensus sequences" 
                            by T. D. Schneider and R. M. Stephens in 1990.
                            """)
                    with st.expander("Calculation of errorbars"):
                        st.write("""
                            The magnitude of the errorbars are based on the small sample correction
                            factor. The small sample correction factor for sequence logos with 
                            amino acids is calculated as
                            """)
                        st.latex(r'''e_n = \frac{1}{\ln{2}} \cdot \frac{20-1}{2\cdot n}''')
                        st.write("""
                            The height of the errorbars correspond to twice the small sample correction
                            factor. 
                            """)
                #elif select_analyses == "Download protein sequences tagged with signal peptide":
                #    st.download_button('Download protein sequences tagged with Sec SPI', prepareDownloadSPs(group_number = str(tax_id), sp_type = select_sp))


            #Sec SPII analyses
            elif select_sp == "Sec SPII":
                select_analyses = st.radio("Which analyses would you like to take a look at?", ("Full signal peptide length distribution", "Length distributions of each signal peptide region", "Sequence logos", "Download protein sequences tagged with signal peptide"))
                if select_analyses == "Full signal peptide length distribution":
                    st.subheader("Histogram of signal peptide length distribution")
                    MakeHistograms(tax_id+"_SP_regions_counts.tab", "LIPO", group = scient_name, see_regions = False, rank = clean_rank)
                elif select_analyses == "Length distributions of each signal peptide region":
                    st.subheader("Histograms of region length distributions")
                    MakeHistograms(tax_id+"_SP_regions_counts.tab", "LIPO", group = scient_name, see_regions = True, rank = clean_rank)
                elif select_analyses == "Sequence logos":
                    st.subheader("Sequence logos aligned at region borders")
                    option1,option2 = st.columns(2)
                    option1_inf,option2_inf = st.columns(2)
                    color_scheme = option1.radio("Select a color scheme for the sequence logos", ("Default", "Hydrophobicity", "Charge", "Chemistry"))
                    with option1_inf.expander("Information about color schemes"):
                        st.image(seq_logo_color_schemes)
                        st.write("""
                            Default color scheme divides amino acids into negatively charged (red), 
                            positively charged (blue), hydrophobic (black) and others (green). This color scheme
                            is useful when examining important signal peptides properties such as
                            the positive charge of n-regions and hydrophobicity of h-regions. 
                            """)

                    keep_zero = option2.radio("Do you want to include 0 as the position before alignment site?", ("Yes", "No"))
                    with option2_inf.expander("Information about position assignments"):
                        st.write("""
                            Choosing either of these options do not influence the logos themselves, 
                            except changing whether the position before alignment site should be 
                            written as position 0 or position -1. 
                            """)
                    #st.caption("Alignment at N-H region border")
                    min_pos_n = ExtractAlignmentPosition(tax_id, "LIPO", "N")[0]
                    len_seq_n = ExtractAlignmentPosition(tax_id, "LIPO", "N")[1]
                    max_pos_n = int(len_seq_n) - int(min_pos_n)
                    st.write("**Select the range of positions in the logo aligned at the n-h region border**")
                    low,high = st.columns(2)
                    
                    if int(min_pos_n) >= 5:
                        lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -5, step = 1)
                    else: 
                        lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -int(min_pos_n), step = 1)
                    
                    if int(max_pos_n) >= 7:
                        upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = 7, step = 1)
                    else: 
                        upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = max_pos_n, step = 1)
                    
                    MakeLogo(tax_id,"LIPO","n", lower_n, upper_n, scient_name, color_scheme, keep_zero)
                    min_pos_h = ExtractAlignmentPosition(tax_id, "LIPO", "H")[0]
                    len_seq_h = ExtractAlignmentPosition(tax_id, "LIPO", "H")[1]
                    max_pos_h = int(len_seq_h) - int(min_pos_h)
                    st.write("**Select the range of positions in the logo aligned at the h-c region border**")
                    low,high = st.columns(2)
                    if int(min_pos_h) >= 7:
                        lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -7, step = 1)
                    else:
                        lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -int(min_pos_h), step = 1)
                    
                    if int(max_pos_h) >= 5:
                        upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = 5, step = 1)
                    else:
                        upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = max_pos_h, step = 1)
                    
                    MakeLogo(tax_id,"LIPO","h", lower_h, upper_h, scient_name, color_scheme, keep_zero)

                    with st.expander("Creation of sequence logos"):
                        st.write("""
                            The sequence logos presented on this website are based on Shannon Entropy,
                            as described in "Sequence logos: a new way to display consensus sequences" 
                            by T. D. Schneider and R. M. Stephens in 1990.
                            """)
                    with st.expander("Calculation of errorbars"):
                        st.write("""
                            The magnitude of the errorbars are based on the small sample correction
                            factor. The small sample correction factor for sequence logos with 
                            amino acids is calculated as
                            """)
                        st.latex(r'''e_n = \frac{1}{\ln{2}} \cdot \frac{20-1}{2\cdot n}''')
                        st.write("""
                            The height of the errorbars correspond to twice the small sample correction
                            factor. 
                            """)

                #elif select_analyses == "Download protein sequences tagged with signal peptide":
                #    st.download_button('Download protein sequences tagged with Sec SPII', prepareDownloadSPs(group_number = str(tax_id), sp_type = select_sp))

            
            #Sec SPIII analyses
            elif select_sp == "Sec SPIII":
                select_analyses = st.radio("Which analyses would you like to take a look at?", ("Full signal peptide length distribution", "Sequence logo", "Download protein sequences tagged with signal peptide"))
                with st.expander("Further information"):
                    st.write("""
                        The Sec SPIII signal peptide only contains one region, which is why histograms of region length distributions cannot be viewed for this signal peptide type.
                        """)

                if select_analyses == "Full signal peptide length distribution":
                    st.subheader("Histogram of signal peptide length distribution")
                    MakeHistograms(tax_id+"_SP_regions_counts.tab", "PILIN", group = scient_name, see_regions = False, rank = clean_rank)
                elif select_analyses == "Sequence logo":
                    st.subheader("Sequence logos aligned at region borders")
                    option1,option2 = st.columns(2)
                    option1_inf,option2_inf = st.columns(2)
                    color_scheme = option1.radio("Select a color scheme for the sequence logos", ("Default", "Hydrophobicity", "Charge", "Chemistry"))
                    with option1_inf.expander("Information about color schemes"):
                        st.image(seq_logo_color_schemes)
                        st.write("""
                            Default color scheme divides amino acids into negatively charged (red), 
                            positively charged (blue), hydrophobic (black) and others (green). This color scheme
                            is useful when examining important signal peptides properties such as
                            the positive charge of n-regions and hydrophobicity of h-regions. 
                            """)

                    keep_zero = option2.radio("Do you want to include 0 as the position before alignment site?", ("Yes", "No"))
                    with option2_inf.expander("Information about position assignments"):
                        st.write("""
                            Choosing either of these options do not influence the logos themselves, 
                            except changing whether the position before alignment site should be 
                            written as position 0 or position -1. 
                            """)
                    #st.caption("Alignment at N-H region border")
                    min_pos_n = ExtractAlignmentPosition(tax_id, "PILIN", "N")[0]
                    len_seq_n = ExtractAlignmentPosition(tax_id, "PILIN", "N")[1]
                    max_pos_n = int(len_seq_n) - int(min_pos_n)
                    st.write("**Select the range of positions in the logo aligned at cleavage site**")
                    low,high = st.columns(2)
                    
                    if int(min_pos_n) >= 7:
                        lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -7, step = 1)
                    else: 
                        lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -int(min_pos_n), step = 1)

                    if int(max_pos_n) >= 12:
                        upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = 12, step = 1)
                    else: 
                        upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = max_pos_n, step = 1)
                    
                    MakeLogo(tax_id,"PILIN","n", lower_n, upper_n, scient_name, color_scheme, keep_zero)
                    
                    with st.expander("Creation of sequence logos"):
                        st.write("""
                            The sequence logos presented on this website are based on Shannon Entropy,
                            as described in "Sequence logos: a new way to display consensus sequences" 
                            by T. D. Schneider and R. M. Stephens in 1990.
                            """)
                    with st.expander("Calculation of errorbars"):
                        st.write("""
                            The magnitude of the errorbars are based on the small sample correction
                            factor. The small sample correction factor for sequence logos with 
                            amino acids is calculated as
                            """)
                        st.latex(r'''e_n = \frac{1}{\ln{2}} \cdot \frac{20-1}{2\cdot n}''')
                        st.write("""
                            The height of the errorbars correspond to twice the small sample correction
                            factor. 
                            """)
                #elif select_analyses == "Download protein sequences tagged with signal peptide":
                #    st.download_button('Download protein sequences tagged with Sec SPIII', prepareDownloadSPs(group_number = str(tax_id), sp_type = select_sp))


            #Tat SPI analyses
            elif select_sp == "Tat SPI":
                select_analyses = st.radio("Which analyses would you like to take a look at?", ("Full signal peptide length distribution", "Length distributions of each signal peptide region", "Sequence logos", "Download protein sequences tagged with signal peptide"))
                if select_analyses == "Full signal peptide length distribution":
                    st.subheader("Histogram of signal peptide length distribution")
                    MakeHistograms(tax_id+"_SP_regions_counts.tab", "TAT", group = scient_name, see_regions = False, rank = clean_rank)
                elif select_analyses == "Length distributions of each signal peptide region":
                    st.subheader("Histograms of region length distributions")
                    MakeHistograms(tax_id+"_SP_regions_counts.tab", "TAT", group = scient_name, see_regions = True, rank = clean_rank)
                elif select_analyses == "Sequence logos":
                    st.subheader("Sequence logos aligned at region borders")
                    option1,option2 = st.columns(2)
                    option1_inf,option2_inf = st.columns(2)
                    color_scheme = option1.radio("Select a color scheme for the sequence logos", ("Default", "Hydrophobicity", "Charge", "Chemistry"))
                    with option1_inf.expander("Information about color schemes"):
                        st.image(seq_logo_color_schemes)
                        st.write("""
                            Default color scheme divides amino acids into negatively charged (red), 
                            positively charged (blue), hydrophobic (black) and others (green). This color scheme
                            is useful when examining important signal peptides properties such as
                            the positive charge of n-regions and hydrophobicity of h-regions. 
                            """)

                    keep_zero = option2.radio("Do you want to include 0 as the position before alignment site?", ("Yes", "No"))
                    with option2_inf.expander("Information about position assignments"):
                        st.write("""
                            Choosing either of these options do not influence the logos themselves, 
                            except changing whether the position before alignment site should be 
                            written as position 0 or position -1. 
                            """)
                    #st.caption("Alignment at N-H region border")
                    min_pos_n = ExtractAlignmentPosition(tax_id, "TAT", "N")[0]
                    len_seq_n = ExtractAlignmentPosition(tax_id, "TAT", "N")[1]
                    max_pos_n = int(len_seq_n) - int(min_pos_n)
                    st.write("**Select the range of positions in the logo aligned at the n-h region border**")
                    low,high = st.columns(2)
                    
                    if int(min_pos_n) >= 5:
                        lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -5, step = 1)
                    else: 
                        lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -int(min_pos_n), step = 1)
                    
                    if int(max_pos_n) >= 7:
                        upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = 7, step = 1)
                    else: 
                        upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = max_pos_n, step = 1)

                    MakeLogo(tax_id,"TAT","n", lower_n, upper_n, scient_name, color_scheme, keep_zero)
                    
                    min_pos_h = ExtractAlignmentPosition(tax_id, "TAT", "H")[0]
                    len_seq_h = ExtractAlignmentPosition(tax_id, "TAT", "H")[1]
                    max_pos_h = int(len_seq_h) - int(min_pos_h)
                    st.write("**Select the range of positions in the logo aligned at the h-c region border**")
                    low,high = st.columns(2)
                    if int(min_pos_h) >= 15:
                        lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -15, step = 1)
                    else:
                        lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -int(min_pos_h), step = 1)
                    
                    if int(max_pos_h) >= 5:
                        upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = 5, step = 1)
                    else:
                        upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = max_pos_h, step = 1)

                    MakeLogo(tax_id,"TAT","h", lower_h, upper_h, scient_name, color_scheme, keep_zero)

                    min_pos_c = ExtractAlignmentPosition(tax_id, "TAT", "C")[0]
                    len_seq_c = ExtractAlignmentPosition(tax_id, "TAT", "C")[1]
                    max_pos_c = int(len_seq_c) - int(min_pos_c)
                    st.write("**Select the range of positions in the logo aligned at cleavage site**")
                    low,high = st.columns(2)
                    
                    if int(min_pos_c) >= 7:
                        lower_c = low.number_input('Lower value', min_value = -int(min_pos_c), max_value = -1, value = -7, step = 1)
                    else: 
                        lower_c = low.number_input('Lower value', min_value = -int(min_pos_c), max_value = -1, value = -int(min_pos_c), step = 1)
                    
                    if int(max_pos_c) >= 5:
                        upper_c = high.number_input('Upper value', min_value = 1, max_value = max_pos_c, value = 5, step = 1)
                    else:
                        upper_c = high.number_input('Upper value', min_value = 1, max_value = max_pos_c, value = max_pos_c, step = 1)

                    MakeLogo(tax_id,"TAT","c", lower_c, upper_c, scient_name, color_scheme, keep_zero)

                    with st.expander("Creation of sequence logos"):
                        st.write("""
                            The sequence logos presented on this website are based on Shannon Entropy,
                            as described in "Sequence logos: a new way to display consensus sequences" 
                            by T. D. Schneider and R. M. Stephens in 1990.
                            """)
                    with st.expander("Calculation of errorbars"):
                        st.write("""
                            The magnitude of the errorbars are based on the small sample correction
                            factor. The small sample correction factor for sequence logos with 
                            amino acids is calculated as
                            """)
                        st.latex(r'''e_n = \frac{1}{\ln{2}} \cdot \frac{20-1}{2\cdot n}''')
                        st.write("""
                            The height of the errorbars correspond to twice the small sample correction
                            factor. 
                            """)
                #elif select_analyses == "Download protein sequences tagged with signal peptide":
                #    st.download_button('Download protein sequences tagged with Tat SPI', prepareDownloadSPs(group_number = str(tax_id), sp_type = select_sp))

            
            #Tat SPII analyses
            elif select_sp == "Tat SPII":
                select_analyses = st.radio("Which analyses would you like to take a look at?", ("Full signal peptide length distribution", "Length distributions of each signal peptide region", "Sequence logos", "Download protein sequences tagged with signal peptide"))
                if select_analyses == "Full signal peptide length distribution":
                    st.subheader("Histogram of signal peptide length distribution")
                    MakeHistograms(tax_id+"_SP_regions_counts.tab", "TATLIPO", group = scient_name, see_regions = False, rank = clean_rank)
                elif select_analyses == "Length distributions of each signal peptide region":
                    st.subheader("Histograms of region length distributions")
                    MakeHistograms(tax_id+"_SP_regions_counts.tab", "TATLIPO", group = scient_name, see_regions = True, rank = clean_rank)
                elif select_analyses == "Sequence logos":
                    st.subheader("Sequence logos aligned at region borders")
                    option1,option2 = st.columns(2)
                    option1_inf,option2_inf = st.columns(2)
                    color_scheme = option1.radio("Select a color scheme for the sequence logos", ("Default", "Hydrophobicity", "Charge", "Chemistry"))
                    with option1_inf.expander("Information about color schemes"):
                        st.image(seq_logo_color_schemes)
                        st.write("""
                            Default color scheme divides amino acids into negatively charged (red), 
                            positively charged (blue), hydrophobic (black) and others (green). This color scheme
                            is useful when examining important signal peptides properties such as
                            the positive charge of n-regions and hydrophobicity of h-regions. 
                            """)

                    keep_zero = option2.radio("Do you want to include 0 as the position before alignment site?", ("Yes", "No"))
                    with option2_inf.expander("Information about position assignments"):
                        st.write("""
                            Choosing either of these options do not influence the logos themselves, 
                            except changing whether the position before alignment site should be 
                            written as position 0 or position -1. 
                            """)

                    #st.caption("Alignment at N-H region border")
                    min_pos_n = ExtractAlignmentPosition(tax_id, "TATLIPO", "N")[0]
                    len_seq_n = ExtractAlignmentPosition(tax_id, "TATLIPO", "N")[1]
                    max_pos_n = int(len_seq_n) - int(min_pos_n)
                    st.write("**Select the range of positions in the logo aligned at the n-h region border**")
                    low,high = st.columns(2)
                    if int(min_pos_n) >= 5:
                        lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -5, step = 1)
                    else:
                        lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -int(min_pos_n), step = 1)
                    
                    if int(max_pos_n) >= 7:
                        upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = 7, step = 1)
                    else: 
                        upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = max_pos_n, step = 1)

                    MakeLogo(tax_id,"TATLIPO","n", lower_n, upper_n, scient_name, color_scheme, keep_zero)
                    #except ValueError:
                    #    st.warning("Your selected value for lower part of the range is too small. Please try with a greater value.")
                    min_pos_h = ExtractAlignmentPosition(tax_id, "LIPO", "H")[0]
                    len_seq_h = ExtractAlignmentPosition(tax_id, "LIPO", "H")[1]
                    max_pos_h = int(len_seq_h) - int(min_pos_h)
                    st.write("**Select the range of positions in the logo aligned at the h-c region border**")
                    low,high = st.columns(2)
                    
                    if int(min_pos_h) >= 7:
                        lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -7, step = 1)
                    else:
                        lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -int(min_pos_h), step = 1)
                    
                    if int(max_pos_h) >= 5:
                        upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = 5, step = 1)
                    else:
                        upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = max_pos_h, step = 1)

                    MakeLogo(tax_id,"TATLIPO","h", lower_h, upper_h, scient_name, color_scheme, keep_zero)

                    with st.expander("Creation of sequence logos"):
                        st.write("""
                            The sequence logos presented on this website are based on Shannon Entropy,
                            as described in "Sequence logos: a new way to display consensus sequences" 
                            by T. D. Schneider and R. M. Stephens in 1990.
                            """)
                    with st.expander("Calculation of errorbars"):
                        st.write("""
                            The magnitude of the errorbars are based on the small sample correction
                            factor. The small sample correction factor for sequence logos with 
                            amino acids is calculated as
                            """)
                        st.latex(r'''e_n = \frac{1}{\ln{2}} \cdot \frac{20-1}{2\cdot n},''')
                        st.write("""
                            with *n* being the number of amino acids at a given position.
                            The height of the errorbars correspond to twice the small sample correction
                            factor. 
                            """)
                #elif select_analyses == "Download protein sequences tagged with signal peptide":
                #    st.download_button('Download protein sequences tagged with Tat SPII', prepareDownloadSPs(group_number = str(tax_id), sp_type = select_sp))

    except NameError as err:
        st.warning('The program was not able to evoke analyses of ' + option + '. Please try with another taxonomical group or organism name.')
        st.write(str(err))


st.text("")
st.text("")
st.text("")
st.text("")
st.text("")

col1, col2 = st.columns([4,1])

col1.write("")

col2.image(dtu_logo)

col2.write("[DTU Health Tech Services](https://services.healthtech.dtu.dk/)")
col2.write("[SignalP (version 6.0)](https://services.healthtech.dtu.dk/service.php?SignalP-6.0)")