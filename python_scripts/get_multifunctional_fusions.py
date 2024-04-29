import numpy as np
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import os
import argparse
import itertools as it

def get_files():
	parser = argparse.ArgumentParser(description='Multifunctional fusion enzymes can be identified in the detailed hits files.')
	parser.add_argument('-i', '--input', help="files for the multifunctional enzymes. Can be found in detailed hits from blast", required=True, nargs='+') # order of input is important 
	parser.add_argument('-n', '--number', required=True, help = "what type of multifunctional enzyme, bi. tri, tetra etc",type=int)
	parser.add_argument('-o', '--out', help="Path to output files", required=True)
	#parser.add_argument('-c', '--combinations', help="Fusion combinations", nargs='+', required=True ,type=str)
	
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


def get_domain_order(accessions,dfs):

	order_dict = {}
	for entry in accessions:
		starts = []
		order = ""
		for df in dfs:
			row = df[0][df[0][1].str.contains(entry)].values.tolist()
			starts.append((row[0][9],df[1]))
		starts.sort()
		#print(starts)
		for i in range(len(starts)):

			order += starts[i][1]
		#print(order)	

		if order not in order_dict:
			order_dict[order] = [entry]
		else:
			order_dict[order].append(entry)
	return(order_dict)
	
def main():
# query from one file could be present as hit in the other file
	infiles = get_files()
	out_path = infiles["out"]
	combinations =  list(it.combinations(infiles["input"],infiles["number"])) # 
	counter = 0
	for i in range(len(combinations)):
		enzyme_accessions = []
		dfs = []
		rest = []
		enzymes = ""
		for j in combinations[i]:
		
			df = pd.read_csv(j, engine='python', sep='\t', header=None)
			enzyme_accessions.append(df[1].tolist())
			dfs.append((df,j[-6]))
			enzymes += j[-6]
		print(enzymes)
		
		intersection = set()
		for j in range(len(enzyme_accessions)-1):
			if j == 0:
				intersection = set(enzyme_accessions[j]).intersection(enzyme_accessions[j+1])
			else: 
				intersection = (set(intersection).intersection(enzyme_accessions[j+1]))
			
		#intersection_2 = (set(intersection_1).intersection(enzyme_accessions[2]))
		#intersection_3 = (set(intersection_2).intersection(enzyme_accessions[3]))
		print("Number of intersecting enzymes (unfiltered): ")
		print(len(intersection))
		for j in infiles['input']:
			if j not in combinations[i]:
				rest.append(j)
		final_accessions = filter_data(intersection, rest)
		print("Total number of multifunctional enzymes:", len(final_accessions))

		
		if len(final_accessions) != 0:
			info = get_domain_order(final_accessions,dfs)
			print(info.keys())
			print(len(final_accessions))
			for key, value in info.items():
				with open(infiles['out'] + key +".txt", "w") as outfile:
					outfile.write("\n".join(value))
			counter += 1
			print("\n")
			
			
		else: 
			print("Skip\n")
	print("Total:", counter, "\n")

	

if __name__ == "__main__":
    main()
