import numpy as np
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import os
import argparse

__author__="Sevvalli Thavapalan"

"""
Takes the detailed hits files from BLAST and cross checks the files for possible bifunctional
enzymes.
"""

def get_files():
	parser = argparse.ArgumentParser(description='')
	# order of input is important 
	parser.add_argument('-i', '--input', help="Files for the bifunctional enzymes. Can be found in detailed hits from BLAST", required=True, nargs='+') 
	parser.add_argument('-r', '--rest', required=True, help = "All other files to filter out mulifunctional enzymes",nargs='+')
	parser.add_argument('-o', '--out', help="Path to output files", required=True)
	parser.add_argument('-c', '--combinations', help="Combination of fusion enzymes", nargs='+', required=True ,type=str)
	
	args = parser.parse_args()
	arguments = args.__dict__

	return arguments


def filter_data(intersection, files):
	"""
	The intersection of the two enzymes are compared against all other enzymes to filter out 
	possible multifunctional enzymes
	"""
	dfs = []
	remove_list = set()
	for file in files:
		print("Filtering: ",file)
		df = pd.read_csv(file, engine='python', sep='\t', header=None) 
		dfs.append(df)
		
		accession_list = df[1].tolist()
		update = (set(intersection).intersection(accession_list))
		print(len(update))
		remove_list |= update
	final = intersection.difference(remove_list)
	return final
	

def get_domain_information(accessions, enzyme1, enzyme2):
	"""
	Information of the bifunctional enzymes are calculated by comparing query and hits. 
	order of information: "Accession", "Start Domain 1","End Domain 1", "Linker", "Start Domain 2","End Domain 2", 
	"Length Enzyme", "Back Overhang"
	"""
	dom12 = []
	dom21 = []
	for entry in accessions:
		cross_1 = enzyme1[enzyme1[1].str.contains(entry)].values.tolist()
		cross_2 = enzyme2[enzyme2[1].str.contains(entry)].values.tolist()
		if  len(cross_1) !=0 and len(cross_2) !=0:
			begin_query1 = cross_1[0][7] -1
			end_query1 = cross_1[0][5]- cross_1[0][8]
		
			begin_dom1 = cross_1[0][9] - begin_query1
			end_dom1 = cross_1[0][10] + end_query1
		
			begin_query2 = cross_2[0][7] -1
			end_query2= cross_2[0][5]- cross_2[0][8]
		
			begin_dom2 = cross_2[0][9] - begin_query2
			end_dom2 = cross_2[0][10] + end_query2 
			
			if begin_dom1 < begin_dom2:
				dom12.append([cross_1[0][1], begin_dom1,1,begin_dom1-1,end_dom1,begin_dom2-end_dom1-1, begin_dom2,end_dom2, cross_1[0][6],end_dom2-cross_1[0][6]])
			else:
				dom21.append([cross_2[0][1],begin_dom2,1,begin_dom2-1,end_dom2,begin_dom1-end_dom2-1,begin_dom1,end_dom1, cross_1[0][6],end_dom1 -cross_1[0][6]])
	
	return dom12, dom21
		

def main():
	infiles = get_files()
	out_path = infiles["out"]
	# query from one file could be present as hit in the other file if they are fused together
	enzyme1 =  pd.read_csv(infiles["input"][0], engine='python', sep='\t', header=None) 
	enzyme2 =  pd.read_csv(infiles["input"][1], engine='python', sep='\t', header=None)
	
	enzyme1_accessions = enzyme1[1].tolist()
	enzyme2_accessions = enzyme2[1].tolist()

	intersection = set(enzyme1_accessions).intersection(enzyme2_accessions) # unfiltered intersection between the two enzymes
	print("Number of intersecting enzymes (unfiltered): ")
	print(len(intersection))
	final_accessions = filter_data(intersection, infiles["rest"])
	print("Total number of bifunctional enzymes:", len(final_accessions))
	
	fusion12, fusion21 = get_domain_information(final_accessions, enzyme1, enzyme2)
	fusion12_df = pd.DataFrame(fusion12)
	# Order: ["Accession", "Start Domain 1","End Domain 1", "Linker", "Start Domain 2","End Domain 2", "Length Enzyme","Back Overhang"] 
	fusion21_df = pd.DataFrame(fusion21)
	out_file1 = fusion12_df.to_csv(str(out_path) + '/' + infiles["combinations"][0] + '.tsv', sep='\t', index = False, header=False)
	out_file2 = fusion21_df.to_csv(str(out_path) + '/' + infiles["combinations"][1] + '.tsv', sep='\t', index = False, header=False)
	
	#print(fusion12_df.head(10))
	print("Files saved in: ")
	print(out_path)
	print("Done")
	

if __name__ == "__main__":
    main()	
