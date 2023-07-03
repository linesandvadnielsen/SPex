#!/usr/bin/env python3

import streamlit as st
import os
import sys
from CountSPs import *
from ExtractEukaryoticGroups import *
from ExtractProtCounts import *
    
def VisualizeSPFrequencies(tax_id):
    """Displays the fractions of proteins tagged with each SP type for a given group to streamlit"""
    euk_tax_groups = ExtractEukaryoticGroups()
    
    #Only show Sec SPI for eukaryotes
    if tax_id in euk_tax_groups:
        #Save region counts
        n_region_sp = CountSPs(tax_id+"_SP_regions_counts.tab")[0]
        
        #Save total protein count for phylogenetic group
        prot_count = ExtractProtCounts(tax_id+"_SP_types_counts.tab")[0]
        
        try:
            #Setup in streamlit
            sp, none_ = st.columns(2)
            
            #Display fractions of each SP type to streamlit, round to 2 significant digits
            sp.metric("Sec SPI", str(round((len(n_region_sp)/int(prot_count)*100), 2))+" %")
            none_.metric("", "")
        except ZeroDivisionError as err:
            st.info("Not Available.")

        #Write explanation to calculations to streamlit
        with st.expander("Further information"):
         st.write("""
            The fraction of proteins tagged with a signal peptide type
            is calculated as the count of signal peptides of a given type
            in this phylogenetic group,
            divided by the total number of proteins in this phylogenetic 
            group. The total number of
            proteins is """, prot_count, """ in this group.""")
    else:
        #Save region counts
        n_region_sp = CountSPs(tax_id+"_SP_regions_counts.tab")[0]
        n_region_lipo = CountSPs(tax_id+"_SP_regions_counts.tab")[3]
        n_region_tat = CountSPs(tax_id+"_SP_regions_counts.tab")[5]
        n_region_tatlipo = CountSPs(tax_id+"_SP_regions_counts.tab")[8]
        n_region_pilin = CountSPs(tax_id+"_SP_regions_counts.tab")[10]


        #Save total protein count for phylogenetic group
        prot_count = ExtractProtCounts(tax_id+"_SP_types_counts.tab")[0]
        
        try:
            #Setup in streamlit
            sp, lipo, tat, tatlipo, pilin = st.columns(5)
            
            #Display fractions of each SP type to streamlit, round to 2 significant digits
            sp.metric("Sec SPI", str(round((len(n_region_sp)/int(prot_count)*100), 2))+" %")
            lipo.metric("Sec SPII", str(round((len(n_region_lipo)/int(prot_count)*100), 2))+" %")
            tat.metric("Tat SPI", str(round((len(n_region_tat)/int(prot_count)*100), 2))+" %")
            tatlipo.metric("Tat SPII", str(round((len(n_region_tatlipo)/int(prot_count)*100), 2))+" %")
            pilin.metric("Sec SPIII", str(round((len(n_region_pilin)/int(prot_count)*100), 2))+" %")
        except ZeroDivisionError as err:
            st.write("Not Available.")

        #Write explanation to calculations to streamlit
        with st.expander("Further information"):
         st.write("""
            The fraction of proteins tagged with a signal peptide type
            is calculated as the count of signal peptides of a given type
            in this phylogenetic group,
            divided by the total number of proteins in this phylogenetic 
            group. The total number of
            proteins is """, prot_count, """ in this group.""")
