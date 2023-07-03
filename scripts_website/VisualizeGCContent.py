#!/usr/bin/env python3

import streamlit as st
import os
import sys

def VisualizeGCContent(tax_id):
	"""Take a file with GC-content information about a tax group, """
	"""visualize it in streamlit"""

	infile_GC = open('../generated_data/data_website/GC_content_phylogenetic_groups/'+tax_id + ".tab", "r")

	#read first line
	infile_GC.readline()

	#read 2nd line with GC-content information, split and store GC-contents
	GC_contents = infile_GC.readline().split("\t")

	GC_coding = GC_contents[1]
	GC1 = GC_contents[2]
	GC2 = GC_contents[3]
	GC3 = GC_contents[4][:-1]

	#If GC-content is not calculated, display message to streamlit
	if GC_contents[1] == "Not Available":
		st.info("Not Available")
	#Display GC-contents
	else:
		gc_coding, gc1, gc2, gc3 = st.columns(4)

		gc_coding.metric("Coding regions", GC_coding)
		gc1.metric("Codon position 1", GC1)
		gc2.metric("Codon position 2", GC2)
		gc3.metric("Codon position 3", GC3)


	infile_GC.close()

