#!/usr/bin/env python3

import os
import streamlit as st
import lmdb
import pickle


def prepareDownloadSPs(group_number, sp_type):
	#sp_type: write either "Sec SPI", "Sec SPII", "Tat SPI", "Tat SPII" or "Sec SPIII"
	
	#Open information about organisms in a domain

	org_tax_info = open("../generated_data/data_website/all_domains_tax_info.tab", "r")
	
	for line in org_tax_info:
		linesplit = line.split("\t")
		if group_number in linesplit[8:]:
			Domain = linesplit[-2]
			break

	if Domain == "Eukaryota":
		domain = "eukarya"
	else:
		domain = Domain.lower()

	org_tax_info.close()

	#Extract all proteome IDs belonging to a domain
	#proteome_IDs_domain = os.listdir("../generated_data/data_prep/Processed_entries/"+domain)
	with open("../generated_data/data_website/Proteome_IDs/"+domain, "rb") as fp:
		proteome_IDs_domain = pickle.load(fp)


	org_tax_info = open("../generated_data/data_website/"+domain+"_tax_info.tab", "r")

	proteome_IDs = []

	#Extract the proteome IDs of organisms belonging to the group
	for line in org_tax_info:
		linesplit = line.split("\t")
		if group_number in linesplit[8:]:
			if linesplit[3]+".fasta" in proteome_IDs_domain:
				proteome_IDs.append(linesplit[3])
				continue

	org_tax_info.close()

	if sp_type == "Sec SPI":
		sp = "SP"
	elif sp_type == "Sec SPII":
		sp = "LIPO"
	elif sp_type == "Sec SPIII":
		sp = "PILIN"
	elif sp_type == "Tat SPI":
		sp = "TAT"
	elif sp_type == "Tat SPII":
		sp = "TATLIPO"

	env_db = lmdb.Environment("../generated_data/data_website/Download_SPs_DB/db_download_" + sp) 

	db = env_db.begin()

	sequences = ""

	for i in range(len(proteome_IDs)):
		sequences += db.get(str(proteome_IDs[i]).encode()).decode()
	
	st.download_button("Download protein sequences from group " + group_number + " tagged with " + sp_type, sequences)
	
	env_db.close()


