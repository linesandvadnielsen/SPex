#!/usr/bin/env python3

import streamlit as st
import os
import sys
from CountSPs import *
from ExtractEukaryoticGroups import *


def VisualizeSPCounts(tax_id):
    """Displays the counts of each SP type for a given group to streamlit"""
    euk_tax_groups = ExtractEukaryoticGroups()
    
    #Only show Sec SPI for eukaryotes
    if tax_id in euk_tax_groups:
        #Save region counts
        n_region_sp = CountSPs(tax_id+"_SP_regions_counts.tab")[0]
        n_region_sp_len = len(n_region_sp)

        if n_region_sp_len == 0:
            st.info("Not Available")

        else:
            #Setup in streamlit
            sp, none_ = st.columns(2)
            
            #Display count of each SP type to streamlit
            sp.metric("Sec SPI", n_region_sp_len)
            none_.metric("", "")
    
    else:    
        #Save region counts
        n_region_sp = CountSPs(tax_id+"_SP_regions_counts.tab")[0]
        n_region_lipo = CountSPs(tax_id+"_SP_regions_counts.tab")[3]
        n_region_tat = CountSPs(tax_id+"_SP_regions_counts.tab")[5]
        n_region_tatlipo = CountSPs(tax_id+"_SP_regions_counts.tab")[8]
        n_region_pilin = CountSPs(tax_id+"_SP_regions_counts.tab")[10]

        
        if len(n_region_sp) == 0 and len(n_region_lipo) == 0 and len(n_region_tat) == 0 and len(n_region_tatlipo) == 0 and len(n_region_pilin) == 0:
            st.info("Not Available")

        else:
            #Setup in streamlit
            sp, lipo, tat, tatlipo, pilin = st.columns(5)
        
            #Display count of each SP type to streamlit
            sp.metric("Sec SPI", len(n_region_sp))
            lipo.metric("Sec SPII", len(n_region_lipo))
            tat.metric("Tat SPI", len(n_region_tat))
            tatlipo.metric("Tat SPII", len(n_region_tatlipo))
            pilin.metric("Sec SPIII", len(n_region_pilin))

