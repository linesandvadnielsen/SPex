#!/usr/bin/env python3

import streamlit as st
import matplotlib.pyplot as plt
import math
import statistics
from CountSPs import *

def MakeHistograms(filename, sp_type, group, see_regions, rank):
    """Make histograms of SP region length distributions,"""
    """given a phylogenetic group and SP type"""

    #Collect data from CountSPs() for histograms
    n_region_sp = CountSPs(filename)[0]
    h_region_sp = CountSPs(filename)[1]
    c_region_sp = CountSPs(filename)[2]
    n_region_lipo = CountSPs(filename)[3]
    h_region_lipo = CountSPs(filename)[4]
    n_region_tat = CountSPs(filename)[5]
    h_region_tat = CountSPs(filename)[6]
    c_region_tat = CountSPs(filename)[7]
    n_region_tatlipo = CountSPs(filename)[8]
    h_region_tatlipo = CountSPs(filename)[9]
    n_region_pilin = CountSPs(filename)[10]

    #Initialize
    len_sp = []
    len_lipo = []
    len_tat = []
    len_tatlipo = []
    len_pilin = []

    #Collect full SP length data
    #Sec SPI
    for i in range(len(n_region_sp)):
        sp_len = int(n_region_sp[i]) + int(h_region_sp[i]) + int(c_region_sp[i])
        len_sp.append(sp_len)

    #Sec SPII
    for i in range(len(n_region_lipo)):
        lipo_len = int(n_region_lipo[i]) + int(h_region_lipo[i])
        len_lipo.append(lipo_len)

    #Tat SPI
    for i in range(len(n_region_tat)):
        tat_len = int(n_region_tat[i]) + int(h_region_tat[i]) + int(c_region_tat[i])
        len_tat.append(tat_len)

    #Tat SPII
    for i in range(len(n_region_tatlipo)):
        tatlipo_len = int(n_region_tatlipo[i]) + int(h_region_tatlipo[i])
        len_tatlipo.append(tatlipo_len)

    #Sec SPIII
    for i in range(len(n_region_pilin)):
        pilin_len = int(n_region_pilin[i])
        len_pilin.append(pilin_len)

    
    def __Histograms(region, sp_type, group, region_lengths, w, hist_col, rank):
        """Make histograms of each reion's length, display to streamlit"""
        #Sort region lengths from shortest to longest
        region_lengths.sort()
        
        #Only produce non-empty histograms
        if len(region_lengths) != 0: 
            
            #Define number of bins (1 bin for each unique integer)
            n = math.ceil((region_lengths[-1] - region_lengths[0])/w)
            plt.style.use('ggplot')
            fig, ax = plt.subplots()
            
            #Create histogram and title, then display to streamlit
            ax.hist(region_lengths, bins = n, density = True, color = hist_col)
            ax.set_title("Distribution of "+ sp_type + " " + region + "-region lengths \n for "+ rank + " " + group)
            st.pyplot(fig)

    #Use this function if all region counts in a region 
    #for a group are the same (rare case)
    def __HistogramsOneLength(region, sp_type, group, region_lengths, hist_col, rank):
        """Make histograms of each reion's length, display to streamlit"""
        #Sort region lengths from shortest to longest
        region_lengths.sort()
        
        #Only produce non-empty histograms
        if len(region_lengths) != 0: 
        
            #Only show 1 bin, as all counts will be the same
            n = 1
            plt.style.use('ggplot')
            fig, ax = plt.subplots()
            
            #Create histogram and title, then display to streamlit
            ax.hist(region_lengths, bins = n, density = True, color = hist_col)
            ax.set_title("Distribution of "+ sp_type + " " + region + "-region lengths \n for "+ rank + " " + group)
            st.pyplot(fig)

    
    def __HistogramsCompleteSP(sp_type, group, sp_lengths, w, hist_col, rank):
        """Make histograms of the full SP lengths, display to streamlit"""
        #Sort region lengths from shortest to longest
        sp_lengths.sort()
        
        #Only produce non-empty histograms
        if len(sp_lengths) != 0: 
            
            #Define number of bins (1 bin for each integer)
            n = math.ceil((sp_lengths[-1] - sp_lengths[0])/w)
            plt.style.use('ggplot')
            fig, ax = plt.subplots()
            #Create histogram and title, then display to streamlit
            ax.hist(sp_lengths, bins = n, density = True, color = hist_col)
            ax.set_title("Distribution of "+ sp_type + " lengths \n for "+ rank + " " + group)
            st.pyplot(fig)

    def __HistogramsCompleteSPOneLength(sp_type, group, sp_lengths, hist_col, rank):
        """Make histograms of the full SP lengths, display to streamlit"""
        #Sort region lengths from shortest to longest
        sp_lengths.sort()
        
        #Only produce non-empty histograms
        if len(sp_lengths) != 0: 
        
            #Define number of bins (1 bin for each integer)
            n = 1
            plt.style.use('ggplot')
            fig, ax = plt.subplots()
        
            #Create histogram and title, then display to streamlit
            ax.hist(sp_lengths, bins = n, density = True, color = hist_col)
            ax.set_title("Distribution of "+ sp_type + " lengths \n for "+ rank + " " + group)
            st.pyplot(fig)

    #Display histograms of SP region lengths of the selected SP type
    if see_regions:
        #Displaying histograms for Sec SPI
        if sp_type == "SP":
            if len(n_region_sp) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                #n-region length histograms
                try:
                    __Histograms("n", "Sec SPI", group, n_region_sp, 1, "blue", rank)
                    st.write("**Mean n-region length**: " + str(round(statistics.mean(n_region_sp), 2)))
                    st.write("**Median n-region length**: " + str(round(statistics.median(n_region_sp), 2)))


                except ValueError as err:
                    __HistogramsOneLength("n", "Sec SPI", group, n_region_sp, "blue", rank)
                    if len(n_region_sp) > 1:
                        st.info("All Sec SPI n-regions within this group consist of " + str(n_region_sp[0]) + " residues.")
                
                #h-region length histograms
                try:
                    __Histograms("h", "Sec SPI", group, h_region_sp, 1, "green", rank)
                    st.write("**Mean h-region length**: " + str(round(statistics.mean(h_region_sp), 2)))
                    st.write("**Median h-region length**: " + str(round(statistics.median(h_region_sp), 2)))
                except ValueError as err:
                    __HistogramsOneLength("h", "Sec SPI", group, h_region_sp, "green", rank)
                    if len(h_region_sp) > 1:
                        st.info("All Sec SPI h-regions within this group consist of " + str(h_region_sp[0]) + " residues.")
                
                #c-region length histograms
                try:
                    __Histograms("c", "Sec SPI", group, c_region_sp, 1, "red", rank)
                    st.write("**Mean c-region length**: " + str(round(statistics.mean(c_region_sp), 2)))
                    st.write("**Median c-region length**: " + str(round(statistics.median(c_region_sp), 2)))
                except ValueError as err:
                    __HistogramsOneLength("c", "Sec SPI", group, c_region_sp, "red", rank)
                    if len(c_region_sp) > 1:
                        st.info("All Sec SPI c-regions within this group consist of " + str(c_region_sp[0]) + " residues.")
        
        #Displaying histograms for Sec SPII
        elif sp_type == "LIPO":
            if len(n_region_lipo) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                #n-region length histograms
                try:
                    __Histograms("n", "Sec SPII", group, n_region_lipo, 1, "blue", rank)
                    st.write("**Mean n-region length**: " + str(round(statistics.mean(n_region_lipo), 2)))
                    st.write("**Median n-region length**: " + str(round(statistics.median(n_region_lipo), 2)))
                except ValueError as err:
                    __HistogramsOneLength("n", "Sec SPII", group, n_region_lipo, "blue", rank)
                    if len(n_region_lipo) > 1:
                        st.info("All Sec SPII n-regions within this group consist of " + str(n_region_lipo[0]) + " residues.")

                #h-region length histograms
                try:
                    __Histograms("h", "Sec SPII", group, h_region_lipo, 1, "green", rank)
                    st.write("**Mean h-region length (including lipobox)**: " + str(round(statistics.mean(h_region_lipo), 2)))
                    st.write("**Median h-region length (including lipobox)**: " + str(round(statistics.median(h_region_lipo), 2)))
                except ValueError as err:
                    __HistogramsOneLength("h", "Sec SPII", group, h_region_lipo, "green", rank)
                    if len(h_region_lipo) > 1:
                        st.info("All Sec SPII h-regions within this group consist of " + str(h_region_lipo[0]) + " residues.")
        
        #Displaying histograms for Tat SPI
        elif sp_type == "TAT":
            if len(n_region_tat) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                #n-region length histograms
                try:
                    __Histograms("n", "Tat SPI", group, n_region_tat, 1, "blue", rank)
                    st.write("**Mean n-region length**: " + str(round(statistics.mean(n_region_tat), 2)))
                    st.write("**Median n-region length**: " + str(round(statistics.median(n_region_tat), 2)))
                except ValueError as err:
                    __HistogramsOneLength("n", "Tat SPI", group, n_region_tat, "blue", rank)
                    if len(n_region_tat) > 1:
                        st.info("All Tat SPI n-regions within this group consist of " + str(n_region_tat[0]) + " residues.")
                
                #h-region length histograms
                try:
                    __Histograms("h", "Tat SPI", group, h_region_tat, 1, "green", rank)
                    st.write("**Mean h-region length**: " + str(round(statistics.mean(h_region_tat), 2)))
                    st.write("**Median h-region length**: " + str(round(statistics.median(h_region_tat), 2)))
                except ValueError as err:
                    __HistogramsOneLength("h", "Tat SPI", group, h_region_tat, "green", rank)
                    if len(h_region_tat) > 1:
                        st.info("All Tat SPI h-regions within this group consist of " + str(h_region_tat[0]) + " residues.")
                
                #c-region length histograms
                try:
                    __Histograms("c", "Tat SPI", group, c_region_tat, 1, "red", rank)
                    st.write("**Mean c-region length**: " + str(round(statistics.mean(c_region_tat), 2)))
                    st.write("**Median c-region length**: " + str(round(statistics.median(c_region_tat), 2)))
                except ValueError as err:
                    __HistogramsOneLength("c", "Tat SPI", group, c_region_tat, "red", rank)
                    if len(c_region_tat) > 1:
                        st.info("All Tat SPI c-regions within this group consist of " + str(c_region_tat[0]) + " residues.")
        
        #Displaying histograms for Tat SPII
        elif sp_type == "TATLIPO":
            if len(n_region_tatlipo) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                #n-region length histograms
                try:
                    __Histograms("n", "Tat SPII", group, n_region_tatlipo, 1, "blue", rank)
                    st.write("**Mean n-region length**: " + str(round(statistics.mean(n_region_tatlipo), 2)))
                    st.write("**Median n-region length**: " + str(round(statistics.median(n_region_tatlipo), 2)))
                except ValueError as err:
                    __HistogramsOneLength("n", "Tat SPII", group, n_region_tatlipo, "blue", rank)
                    if len(n_region_tatlipo) > 1:
                        st.info("All Tat SPII n-regions within this group consist of " + str(n_region_tatlipo[0]) + " residues.")
                
                #h-region length histograms
                try:
                    __Histograms("h", "Tat SPII", group, h_region_tatlipo, 1, "green", rank)
                    st.write("**Mean h-region length (including lipobox)**: " + str(round(statistics.mean(h_region_tatlipo), 2)))
                    st.write("**Median h-region length (including lipobox)**: " + str(round(statistics.median(h_region_tatlipo), 2)))
                except ValueError as err:
                    __HistogramsOneLength("h", "Tat SPII", group, h_region_tatlipo, "green", rank)
                    if len(h_region_tatlipo) > 1:
                        st.info("All Tat SPII h-regions within this group consist of " + str(h_region_tatlipo[0]) + " residues.")
        
        #Displaying histograms for Sec SPIII
        elif sp_type == "PILIN":
            if len(n_region_pilin) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                #n-region length histograms
                try:
                    __Histograms("n", "Sec SPIII", group, n_region_pilin, 1, "blue", rank)
                    st.write("**Mean length**: " + str(round(statistics.mean(n_region_pilin), 2)))
                    st.write("**Median length**: " + str(round(statistics.median(n_region_pilin), 2)))
                except ValueError as err:
                    __HistogramsOneLength("n", "Sec SPIII", group, n_region_pilin, "blue", rank)
                    if len(n_region_pilin) > 1:
                        st.info("All Sec SPIII signal peptides within this group contains " + str(n_region_pilin[0]) + " residues.")

    #If histogram of the complete SP length should be shown
    if not see_regions:
        #Displaying histogram for Sec SPI
        if sp_type == "SP":
            if len(len_sp) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __HistogramsCompleteSP("Sec SPI", group, len_sp, 1, "purple", rank)
                    st.write("**Mean length**: " + str(round(statistics.mean(len_sp), 2)))
                    st.write("**Median length**: " + str(round(statistics.median(len_sp), 2)))
                except ValueError as err:
                    __HistogramsCompleteSPOneLength("Sec SPI", group, len_sp, "purple", rank)
                    if len(len_sp) > 1:
                        st.info("All Sec SPI signal peptides within this group consist of " + str(len_sp[0]) + " residues.")
        
        #Displaying histogram for Sec SPII
        elif sp_type == "LIPO":
            if len(len_lipo) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __HistogramsCompleteSP("Sec SPII", group, len_lipo, 1, "purple", rank)
                    st.write("**Mean length**: " + str(round(statistics.mean(len_lipo), 2)))
                    st.write("**Median length**: " + str(round(statistics.median(len_lipo), 2)))
                except ValueError as err:
                    __HistogramsCompleteSPOneLength("Sec SPII", group, len_lipo, "purple", rank)
                    if len(len_lipo) > 1:
                        st.info("All Sec SPII signal peptides within this group consist of " + str(len_lipo[0]) + " residues.")

        #Displaying histogram for Tat SPI
        elif sp_type == "TAT":
            if len(len_tat) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __HistogramsCompleteSP("Tat SPI", group, len_tat, 1, "purple", rank)
                    st.write("**Mean length**: " + str(round(statistics.mean(len_tat), 2)))
                    st.write("**Median length**: " + str(round(statistics.median(len_tat), 2)))
                except ValueError as err:
                    __HistogramsCompleteSPOneLength("Tat SPI", group, len_tat, "purple", rank)
                    if len(len_tat) > 1:
                        st.info("All Tat SPI signal peptides within this group consist of " + str(len_tat[0]) + " residues.")

        #Displaying histogram for Tat SPII
        elif sp_type == "TATLIPO":
            if len(len_tatlipo) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __HistogramsCompleteSP("Tat SPII", group, len_tatlipo, 1, "purple", rank)
                    st.write("**Mean length**: " + str(round(statistics.mean(len_tatlipo), 2)))
                    st.write("**Median length**: " + str(round(statistics.median(len_tatlipo), 2)))
                except ValueError as err:
                    __HistogramsCompleteSPOneLength("Tat SPII", group, len_tatlipo, "purple", rank)
                    if len(len_tatlipo) > 1:
                        st.info("All Tat SPII signal peptides within this group consist of " + str(len_tatlipo[0]) + " residues.")

        #Displaying histogram for Sec SPIII
        elif sp_type == "PILIN":
            if len(len_pilin) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __HistogramsCompleteSP("Sec SPIII", group, len_pilin, 1, "purple", rank)
                    st.write("**Mean length**: " + str(round(statistics.mean(len_pilin), 2)))
                    st.write("**Median length**: " + str(round(statistics.median(len_pilin), 2)))
                except ValueError as err:
                    __HistogramsCompleteSPOneLength("Sec SPIII", group, len_pilin, "purple", rank)
                    if len(len_pilin) > 1:
                        st.info("All Sec SPIII signal peptides within this group consist of " + str(len_pilin[0]) + " residues.")
