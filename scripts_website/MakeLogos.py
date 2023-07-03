#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import math
import streamlit as st
import sys
from PIL import Image
import logomaker
import gzip
import io

def ExtractAlignmentPosition(tax_id,sp_type,region):
	"""Extracts aligned sequences from gzipped file"""

	#Ensure correct case for input file
	region_lowercase = region.lower()
	region_uppercase = region.upper()
	sp_type = sp_type.upper()

	#Read gzipped content, save in object
	input_ = gzip.GzipFile('../generated_data/data_website/Aligned_sequences/'+tax_id+"/"+sp_type+"_align_"+region_lowercase+".txt.gz", 'rb')
	zipfile_content = input_.read()
	input_.close()

	#Write saved object to file to unzip it
	output = open('../generated_data/data_website/Aligned_sequences/'+tax_id+"/"+sp_type+"_align_"+region_lowercase+".txt", 'wb')
	output.write(zipfile_content)
	output.close()

	#Open alignment-file
	infile = open('../generated_data/data_website/Aligned_sequences/'+tax_id+"/"+sp_type+"_align_"+region_lowercase+".txt", "r")
	line = infile.readline()
	
	#Save the length of the aligned sequences
	len_seqs = len(line[:-1])

	infile.close()

	#Open file showing positions for alignments
	infile_positions = open('../generated_data/data_website/Aligned_sequences/'+tax_id+"/"+"align_positions.tab", "r")
	region_matched = False

	try:
		for line in infile_positions:
			splits = line.split("\t")

			#Find information about alignment positions for the chosen SP type
			if splits[0] == sp_type:
				
				#Select the border alignments should be based on (N, H, C)
				if splits[1] == region_uppercase:
					
					#Extract position 1 of logo (selected region border)
					alignment_position = splits[2] 	
					region_matched = True

		if not region_matched:
			raise ValueError
	except ValueError as err:
		print("Error! The signal peptide type {0} does not have an {1}-region.".format(sp_type,region))

	infile_positions.close()
	#Remove the unzipped file (only store gzipped version)
	os.remove('../generated_data/data_website/Aligned_sequences/'+tax_id+"/"+sp_type+"_align_"+region_lowercase+".txt")

	return alignment_position, len_seqs


def MakeLogo(tax_id,sp_type,region, low_limit, high_limit, scient_name, color_scheme, keep_zero, max_bits_plot):
	"""Create sequence logos to be display using Streamlit introducing dynamic elements"""
	
	#Define default color scheme 
	color_1990 = {"K": "blue", "R": "blue", "H": "blue",
	"D": "red", "E": "red",
	"A": "black", "V": "black", "L": "black", "I": "black",
	"P": "black", "W": "black", "F": "black", "M": "black",
	"N": "green", "G": "green", "Q": "green", "S": "green",
	"C": "green", "Y": "green", "T": "green",
	"X": "magenta"}

	#Define colorscheme-codes
	if color_scheme == "Default":
		color_logo = color_1990
	elif color_scheme == "Hydrophobicity":
		color_logo = "hydrophobicity"
	elif color_scheme == "Charge":
		color_logo = "charge"
	elif color_scheme == "Chemistry":
		color_logo = "chemistry"

	#Ensure correct case for input file
	region_lowercase = region.lower()
	region_uppercase = region.upper()
	sp_type = sp_type.upper()
	
	#Read gzipped content, save in object
	input_ = gzip.GzipFile('../generated_data/data_website/Aligned_sequences/'+tax_id+"/"+sp_type+"_align_"+region_lowercase+".txt.gz", 'rb')
	zipfile_content = input_.read()
	input_.close()

	#Write saved object to file to unzip it
	output = open('../generated_data/data_website/Aligned_sequences/'+tax_id+"/"+sp_type+"_align_"+region_lowercase+".txt", 'wb')
	output.write(zipfile_content)
	output.close()
	
	#Open alignment-file
	infile = open('../generated_data/data_website/Aligned_sequences/'+tax_id+"/"+sp_type+"_align_"+region_lowercase+".txt", "r")

	seqs = []

	#Store all SP sequences as strings in a list
	for line in infile:
		seqs.append(line[:-1])

	infile.close()

	#Open file showing positions for alignments
	infile_positions = open('../generated_data/data_website/Aligned_sequences/'+tax_id+"/"+"align_positions.tab", "r")
	region_matched = False

	try:
		for line in infile_positions:
			splits = line.split("\t")

			#Find information about alignment positions for the chosen SP type
			if splits[0] == sp_type:
				
				#Select the border alignments should be based on (N, H, C)
				if splits[1] == region_uppercase:
					
					#Extract position 1 of logo (selected region border)
					alignment_position = splits[2] 	
					region_matched = True

		if not region_matched:
			raise ValueError
	except ValueError as err:
		print("Error! The signal peptide type {0} does not have an {1}-region.".format(sp_type,region))

	infile_positions.close()

	#Remove unzipped file
	os.remove('../generated_data/data_website/Aligned_sequences/'+tax_id+"/"+sp_type+"_align_"+region_lowercase+".txt")

	if int(alignment_position) == 0:
		st.info("No signal peptides of type " + sp_type + " were detected for this group.")

	else:
		def __SetPositions(low_limit, high_limit, keep_zero):
			"""Set lower limit and upper limit of position range for logo,"""
			"""Make sure input is correctly given"""
			try:
				if int(low_limit) > 0:
					raise ValueError
			except ValueError as err:
				print("The lower limit of the selected range must be a negative number. Your input lower input number {0} will be converted to -{0}.".format(low_limit))
				low_limit = -int(low_limit)
			try:
				if int(high_limit) < 0:
					raise ValueError
			except ValueError as err:
				print("The high limit of the selected range must be a positive number. Your input upper input number {0} will be converted to {1}.".format(str(high_limit), str(high_limit)[1:]))
				high_limit = str(high_limit)[1:]
				high_limit = int(high_limit)

			#Store positions in a list
			positions = list(range(low_limit,high_limit+1))

			#Determine if position "0" should be displayed
			if keep_zero == "Yes":
				positions.remove(int(low_limit))
			elif keep_zero == "No":
				positions.remove(0)

			#Save number of positions to be in logo
			number_of_pos = int(high_limit) - int(low_limit)

			#Plot-window settings
			set_boundary = -low_limit - 0.5

			return positions, low_limit, high_limit, number_of_pos, set_boundary

		#Extract output from inner function 
		func_output = __SetPositions(low_limit, high_limit, keep_zero)
		pos_range = func_output[0]
		lower = func_output[1] + int(alignment_position)
		upper = func_output[2] + int(alignment_position) - 1
		number_of_pos = func_output[3]
		set_boundary = func_output[4]

		#Make a dataframe of amino acid counts in each position in each sequence
		count_df = logomaker.alignment_to_matrix(sequences = seqs, to_type = "counts", characters_to_ignore = "-")

		Rseq = []
		inf_df = count_df.copy()

		#Compute preparation of position heights
		for i in range(len(count_df)):
			n = sum(count_df.loc[i])
			if n != 0:
				e_n = (20-1)/(2*math.log(2)*n)
				Hl = 0
				for j in range(len(count_df.loc[i])):
					freq = (count_df.loc[i][j]/n)
					if freq > 0:
						logfreq = math.log2(freq)
					else:
						logfreq = 0
					Hl += freq*logfreq
			else:
				e_n = 0
				Hl = 0
				freq = 0
				logfreq = 0
			H = -Hl
			if math.log2(20)-(H+e_n) < 0:
				Rseq.append(0)
			else:
				Rseq.append(math.log2(20)-(H+e_n))

		#Calculate information matrix based on count matrix
		for i in range(len(count_df)):
			n = sum(count_df.loc[i])
			for j in range(len(count_df.loc[i])):
				if n != 0:
					freq = count_df.loc[i][j]/n
					inf_df.loc[i][j] = freq*Rseq[i]
				else:
					freq = 0
					inf_df.loc[i][j] = freq*Rseq[i]

		#Only keep positions to be visualized
		inf_df = inf_df.loc[lower:upper]
		inf_df_cleaned = inf_df.copy()

		inf_df_cleaned.index = list(range(0,number_of_pos))

		count_df = count_df.loc[lower:upper]
		count_df_cleaned = count_df.copy()
		count_df_cleaned.index = list(range(0,number_of_pos))

		E_n = []
		heights = []

		#Calculate small sample correction factor
		for i in range(len(count_df_cleaned)):
		    n = sum(count_df_cleaned.loc[i])
		    if n != 0:
		    	e_n = (20-1)/(2*math.log(2)*n)
		    
		    else:
		    	e_n = 0
		    
		    #Calculation of error bar magnitude
		    errorbar = 2*e_n
		    E_n.append(errorbar)

		    #Append position height information to list
		    heights.append(sum(inf_df_cleaned.loc[i]))

		##Preparation for logo styling depending on SP type and region border
		if sp_type == "SP":
			sp_type_desc = "Sec SPI"
			if region_lowercase == "n":
				region_border = "n-h region border"
			elif region_lowercase == "h":
				region_border = "h-c region border"
			elif region_lowercase == "c":
				region_border = "cleavage site"
		
		elif sp_type == "LIPO":
			sp_type_desc = "Sec SPII"
			if region_lowercase == "n":
				region_border = "n-h region border"
			elif region_lowercase == "h":
				region_border = "cleavage site"

		elif sp_type == "PILIN":
			sp_type_desc = "Sec SPIII"
			if region_lowercase == "n":
				region_border = "cleavage site"
		
		elif sp_type == "TAT":
			sp_type_desc = "Tat SPI"
			if region_lowercase == "n":
				region_border = "n-h region border"
			elif region_lowercase == "h":
				region_border = "h-c region border"
			elif region_lowercase == "c":
				region_border = "cleavage site"
		
		elif sp_type == "TATLIPO":
			sp_type_desc = "Tat SPII"
			if region_lowercase == "n":
				region_border = "n-h region border"
			elif region_lowercase == "h":
				region_border = "cleavage site"

		##Make the Logo
		fig, ax = plt.subplots(1,1,figsize=[4,2])

		seq_logo = logomaker.Logo(inf_df_cleaned,
			stack_order = "big_on_top",
			font_name='DejaVu Sans',
			color_scheme=color_logo,
			vpad=0)

		#Style using Logo methods
		seq_logo.style_xticks(anchor=0, spacing=30, rotation=0)

		seq_logo.ax.set_xticks(range(len(pos_range)))
		seq_logo.ax.set_xticklabels(x for x in pos_range)
		seq_logo.ax.axvline(set_boundary, color='k', linewidth=1, linestyle=':')
		seq_logo.ax.set_title(sp_type_desc + " " + region_border + " in " + scient_name)
		seq_logo.ax.set_ylabel('Information (bits)')
		seq_logo.ax.set_xlabel('Position')
		seq_logo.ax.set_xlim([-1, len(inf_df_cleaned)])
		seq_logo.ax.set_ylim([0, max_bits_plot])	##Customize here

		#Prepare errorbar visualization; make error bar for each position
		x = np.linspace(0,number_of_pos-1,number_of_pos)
		y = heights

		#Define error for each position
		xerr = 0
		yerr = E_n
		
		#Insert errorbars
		seq_logo.ax.errorbar(x, y,
		            xerr=xerr,
		            yerr=yerr,
		            fmt='.',
		           color='grey',
		           ecolor='grey')
		seq_logo.ax.grid(False)

		logo_bytes = io.BytesIO()
		plt.savefig(logo_bytes, format='png')
		logo_bytes.seek(0)
		im = Image.open(logo_bytes)

		#Display the logo to streamlit
		st.image(im)