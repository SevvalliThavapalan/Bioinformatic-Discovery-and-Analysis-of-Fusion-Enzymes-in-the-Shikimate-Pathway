from __future__ import print_function
import pandas as pd
import argparse
import numpy as np
import math
import matplotlib.pyplot as plt
from pathlib import Path
from Bio import SeqIO

__author__="Sevvalli Thavapalan"

def get_files():
	parser = argparse.ArgumentParser(description='Calculates the average distances for each flanking enzyme. Maps protein to the gene using the '
												'.fasta_cluster_out file. Results are saved as Right_and_left_distances.tsv.')
	parser.add_argument('-i', '--input', help="File with distances which was created using calculate_operon_distances.py", required=True, nargs='+')
	parser.add_argument('-f', '--fasta', help="File created by FlaGs called *.fasta_cluster_out", required=True, nargs='+')
	parser.add_argument('-o', '--out', help="Path to output files", required=True, nargs='+')
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments


def get_accessions_distances(results):
	'''
	The distances between query and the flanking genes is determined. As the unique enzymes are called individualy there are different cases 
	based on how many flanking genes on each side were found.
	'''
	print("Determining Distances...")
	queries = results.Query.unique()

	left_accessions = []
	right_accessions = []
	distances_right = []
	distances_left = []
	
	check = []
	for query in queries:
		i = results[results.Query==query]
	
		for j in range(len(i)):
			if (len(i) > 9) and (len(i)<= 17) : # max 8 flanking genes on each side
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
	print("Check if the following numbers are equal or the data frame can not be created:")				
	print(len(right_accessions))
	print(len(left_accessions))
	print(len(queries))
	print(len(set(check)))
	print("Done")
	return left_accessions, right_accessions, distances_right, distances_left



def get_protein_names(fasta, right_accessions, left_accessions):
	'''
	The given fasta file is used to create a list of possible protein names for the flanking genes.
	As we are only interested in the direct neighbours its easier to filter the file first before mapping
	each accession to the protein. 
	'''

	prot_list = []
	FastaFile = open(fasta, 'rU')
	for rec in SeqIO.parse(FastaFile, 'fasta'):
		name = rec.description
		seq = rec.seq
		seqLen = len(rec)
		for entry in right_accessions:
			if str(entry) in name:
				#print(name)
				prot_list.append((name, seqLen))
				break
		for entry in left_accessions:
			if str(entry) in name:
				prot_list.append((name, seqLen))
				break
			
	#print(prot_list)	
	return prot_list
			

def map_accessions_to_proteins(right_accessions, left_accessions, prot_list):
	'''
	Maps the accessions of the flanking proteins to the protein names.
	'''
	print("Mapping accessions to proteins:...")
	right_proteins = []
	right_len = []
	left_proteins = []
	left_len=[]
	for entry_r in right_accessions:
		#print("entry")
		#print(entry_r)
		for prot in prot_list:
			x = prot[0].split("|")
			#print(x)	# access only the first part containing the name of the protein
			if str(entry_r) in x[1]:
				name = x[1].split(" ",1)
				right_proteins.append(name[1])
				right_len.append(prot[1])
				break
		else:
			right_proteins.append(entry_r)
			right_len.append(np.nan)
				#break
	print(len(right_proteins))
	for entry_l in left_accessions:
		#print(entry_l)
		for prot in prot_list:
			#print(prot)
			x = prot[0].split("|")	
			if str(entry_l) in x[1]:
				name = x[1].split(" ",1)
				left_proteins.append(name[1])
				left_len.append(prot[1])
				break
		else:
			left_proteins.append(entry_l)
			left_len.append(np.nan)
				
	print(len(left_proteins))
	print("Done")
	return right_proteins, left_proteins, right_len, left_len



def calc_stats(result, out_path):
	'''
	Calculates the abundance of each protein and the average distance.
	'''
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
	for i, fasta_out, out in zip(infiles['input'], infiles['fasta'], infiles['out']):
		print(i, fasta_out, out)
		results = pd.read_csv(i, engine='python', sep='\t')
		queries = results.Query.unique()
		left_accessions, right_accessions, distances_right, distances_left = get_accessions_distances(results)
		prot_list = get_protein_names(fasta_out, right_accessions, left_accessions)
		right_proteins, left_proteins, right_len, left_len = map_accessions_to_proteins(right_accessions, left_accessions, prot_list)
	
	
		result = pd.DataFrame({
		   	"Query": queries ,
		   	"Accession_right": right_accessions,
		   
			"Enzyme_right": right_proteins,
		   	"Distance_right" : distances_right,
		   	"Length_right": right_len,
		   	"Accession_left": left_accessions,
		   	"Enzyme_left": left_proteins,
		   	"Distance_left": distances_left,
		   	"Length_left": left_len
	
		})
	
		out_path = Path(out)
		result.to_csv(str(out_path) + '/Right_and_left_distances.tsv', sep='\t', index=False)
		#stats = calc_stats(result, out_path) # Is not needed 
	

if __name__ == "__main__":
    main()



