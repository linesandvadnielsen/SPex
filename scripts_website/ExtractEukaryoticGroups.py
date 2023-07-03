#!/usr/bin/env python3

import os
import sys

def ExtractEukaryoticGroups():
	"""Make a file of all taxonomic group numbers belonging to the eukaryotic domain"""

	#Open file with information about eukaryotic domain taxonomies
	infile_eukarya = open("../generated_data/data_website/eukarya_tax_info.tab", "r")

	euk_tax_groups = []

	#Extract all taxonomical group IDs belonging to eukaryotic domain
	for line in infile_eukarya:
		#Extract taxonomical information from file
		tax_data = line.split("\t")[8:]
		for element in tax_data:
			#Extract all tax IDs once
			if element.isnumeric():
				if element not in euk_tax_groups:
					euk_tax_groups.append(element)

	infile_eukarya.close()

	return euk_tax_groups
