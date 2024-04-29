from __future__ import print_function
import pandas as pd
import argparse
import numpy as np
import math
import matplotlib.pyplot as plt
from pathlib import Path



def get_files():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-i', '--input', help="File with distances which was created using calculate_operon_distances.py", required=True)
	parser.add_argument('-f', '--fasta', help="File created by FlaGs called *.fasta_cluster_out", required=True)
	parser.add_argument('-o', '--out', help="Path to output files", required=True)
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments


def get_accessions_distances(results):
	print("Determining Distances: ")
	queries = results.Query.unique()

	left_accessions = []
	right_accessions = []
	distances_right = []
	distances_left = []
	
	check = []
	for query in queries:
		i = results[results.Query==query]
	
		for j in range(len(i)):
			if (len(i) > 9) and (len(i)<= 17) :
				if i.iloc[j]['Start_position_relative_to_query'] == 1:
					l = i.iloc[j-1]['Prot_accession_and_order_of_output']
					split_l = l.split('#')
					left_accessions.append(split_l[0])
					distances_left.append(i.iloc[j-1]['Distances'])
					r = i.iloc[j+1]['Prot_accession_and_order_of_output'] # j+1
					split_r = r.split('#')
					right_accessions.append(split_r[0])
					distances_right.append(i.iloc[j+1]['Distances'])
					check.append(query)
					break			
			else:
				if (i.iloc[j]['Start_position_relative_to_query'] ==1) and (j == 0):
					r = i.iloc[j+1]['Prot_accession_and_order_of_output']
					split_r = r.split('#')
					right_accessions.append(split_r[0])
					distances_right.append(i.iloc[j+1]['Distances'])
					left_accessions.append(np.nan)
					distances_left.append(np.nan)
					check.append(query)
					break
				elif (i.iloc[j]['Start_position_relative_to_query'] ==1) and (j == len(i)-1):
					l = i.iloc[j-1]['Prot_accession_and_order_of_output']
					split_l = l.split('#')
					left_accessions.append(split_l[0])
					distances_left.append(i.iloc[j-1]['Distances'])
					right_accessions.append(np.nan)
					distances_right.append(np.nan)
					check.append(query)
					break
				elif (i.iloc[j]['Start_position_relative_to_query'] ==1) and (j!=0):
					l = i.iloc[j-1]['Prot_accession_and_order_of_output']
					split_l = l.split('#')
					left_accessions.append(split_l[0])
					distances_left.append(i.iloc[j-1]['Distances'])
					r = i.iloc[j+1]['Prot_accession_and_order_of_output']
					split_r = r.split('#')
					right_accessions.append(split_r[0])
					distances_right.append(i.iloc[j+1]['Distances'])
					check.append(query)
					break
	print("Check if the following numbers are equal:")				
	print(len(right_accessions))
	print(len(left_accessions))
	print(len(queries))
	print(len(set(check)))
	print("Done")
	return left_accessions, right_accessions, distances_right, distances_left



def get_protein_names(fasta, right_accessions, left_accessions):

	prot_list = []

	file1 = open(fasta, "r")
	for line in file1:
		if line.startswith('>'):
			for entry in right_accessions:
				if str(entry) in line:
					prot_list.append(line)
			for entry in left_accessions:
				if str(entry) in line:
					prot_list.append(line)
	file1.close()
	
	return prot_list
			

def map_accessions_to_proteins(right_accessions, left_accessions, prot_list):
	print("Mapping accessions to proteins:...")
	right_proteins = []
	left_proteins = []
	for entry_r in right_accessions:
		for prot in prot_list:
			x = prot[1:-1].split("|")	
			if str(entry_r) in x[0]:
				name = x[1].split(" ",1)
				right_proteins.append(name[1])
				break
		else:
			right_proteins.append(entry_r)
	

	for entry_l in left_accessions:
		#print(entry_l)
		for prot in prot_list:
			#print(prot)
			x = prot[1:-1].split("|")	
			if str(entry_l) in x[0]:
				name = x[1].split(" ",1)
				left_proteins.append(name[1])
				break
		else:
			left_proteins.append(entry_l)
	print("Done")
	return right_proteins, left_proteins



def calc_stats(result, out_path):
	print("Calculating abundance: ")

	unique_enzymes_right = result["Enzyme_right"].unique()
	unique_enzymes_left = result["Enzyme_left"].unique()
	with open(str(out_path) +"/stats.tsv", "w") as f:

		print("Average distances of right flanking genes: \n", file=f)
		for entry in unique_enzymes_right:
			sum_distances = 0
			i = result[result.Enzyme_right==entry]
			if str(entry) != "nan":
				for j in range(len(i)):
					sum_distances += i.iloc[j]["Distance_right"]
				print(len(i), "\t", entry, "\t", round(sum_distances/len(i),3), file=f)
				print(len(i), "\t", entry, "\t", round(sum_distances/len(i),3))
	
	
		print("\nAverage distances of left flanking genes: \n", file=f)
		for entry in unique_enzymes_left:
			sum_distances = 0
			i = result[result.Enzyme_left==entry]
			if str(entry) != "nan":
				for j in range(len(i)):
					sum_distances += i.iloc[j]["Distance_left"]
				print(len(i), "\t", entry, "\t", round(sum_distances/len(i), 3), file=f)
	print("Done")
					
		

def main():
	infiles = get_files()
	results = pd.read_csv(infiles['input'], engine='python', sep='\t')
	fasta_out = infiles['fasta']
	queries = results.Query.unique()
	left_accessions, right_accessions, distances_right, distances_left = get_accessions_distances(results)
	prot_list = get_protein_names(fasta_out, right_accessions, left_accessions)
	right_proteins, left_proteins = map_accessions_to_proteins(right_accessions, left_accessions, prot_list)
	
	
	result = pd.DataFrame({
		   "Query": queries ,
		   "Accession_right": right_accessions,
		   "Enzyme_right": right_proteins,
		   "Distance_right" : distances_right,
		   "Accession_left": left_accessions,
		   "Enzyme_left": left_proteins,
		   "Distance_left": distances_left
	
	})
	
	out_path = Path(infiles['out'])
	result.to_csv(str(out_path) + '/Right_and_left_distances.tsv', sep='\t', index=False)
	stats = calc_stats(result, out_path)


	

if __name__ == "__main__":
    main()



