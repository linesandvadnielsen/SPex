#!/usr/bin/env python3

import pickle
import os
import sys

def TreeSearchDict(domain):
    """Open a: 1. dictionary with supgroups to each group"""
    """2. dictionary with parent groups to each group"""
    

    #Open subgroup dictionary
    infile_dict = open('../generated_data/data_website/'+domain+"_tax_groups_children.pk", "rb")
    tree_search_file_subgroups = pickle.load(infile_dict)
    infile_dict.close()

    #Open parent-group dictionary
    infile_dict_parents = open('../generated_data/data_website/'+domain+"_tax_groups_parents.pk", "rb")
    tree_search_file_parents = pickle.load(infile_dict_parents)
    infile_dict_parents.close()

    
    return tree_search_file_subgroups, tree_search_file_parents
