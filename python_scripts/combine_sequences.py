import pandas as pd
import argparse
from pathlib import Path
import sys
import numpy as np
import time
from collections import Counter, defaultdict
from Bio import Entrez
from analyze_cluster_distances import read_accessions



def parse_arguments():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-i', '--input', help="Right and left distances", required=True, nargs='+')
	parser.add_argument('-f', '--fusions', required=True, help='Linker and length info of fusions')
	parser.add_argument('-c', '--cluster_accessions', required=True, help='')
	parser.add_argument('-a', '--accessions', required=True, help = "Accessions of the flanking protein")
	parser.add_argument('-s', '--side', required=True, type=str)
	#parser.add_argument('-o', '--out', help="Path to output files", required=True)

	args = parser.parse_args()
	arguments = args.__dict__
	return arguments

def main(): 
	final_cluster_accessions = []
	final_flags_accessions = []
	infiles = parse_arguments()
	flanking_proteins = []
	entrezDBName = 'protein'
	Entrez.email = 'sevvalli.thavapalan@student.uni-tuebingen.de'
	API_KEY = '76340e02587e2d5027b1cb302cf0b11fc808' 
	
	# Flanking proteins 

	# Accessions obtained from cluster 
	cluster_accessions = read_accessions(infiles['cluster_accessions']) 
	print("Cluster accessions read")
	#print(cluster_accessions)

	for line in cluster_accessions:
		final_cluster_accessions.append(line[0:12])
	#print(final_cluster_accessions)

	# Accessions obtained from FlagS
	flags_accessions = read_accessions(infiles['accessions'])
	print("FlagS accessions read")


	for line in flags_accessions:
		y = line.rstrip("\n")
		final_flags_accessions.append(y[:-2])
	#print(final_flags_accessions)

	final_accessions = [x for x in final_cluster_accessions  if x in final_flags_accessions]
	print(len(final_accessions))
	
	for file in infiles['input']:
		print("Determining distances")
		results = pd.read_csv(file, engine='python', sep='\t')

		for line in final_accessions:
			i = results[results["Query"].str.contains(line)]
			if infiles["side"] == "r":
				if len(i) != 0:
					accession_right = i["Accession_right"].values
					distance_right = i["Distance_right"].values
					#if (distance_right[0] == -4.0) or (distance_right[0] == -1.0) :
					flanking_proteins.append((line,accession_right[0][:-2],distance_right[0] )) #first enzyme then neighbour
			else:
				if len(i) != 0:
					accession_left = i["Accession_left"].values
					distance_left = i["Distance_left"].values
					#if (distance_left[0] == -4.0) or (distance_left[0] == -1.0):
					flanking_proteins.append((line,accession_left[0][:-2],distance_left[0]))
	print(len(flanking_proteins))
	
	fusions = pd.read_csv(infiles['fusions'], engine='python', sep='\t', header=None)
	fusion_accession = []
	linker = []
	for lines in final_cluster_accessions:
		#print(lines[0:12])
		i = fusions[fusions.iloc[:,0].str.contains(lines[0:12])].values
		if len(i) != 0:
			linker.append(i[0][5]*3) #nucleotides
			fusion_accession.append(i[0][0])

	counted_linker = Counter(linker)
	
	c = [fusion_accession, linker] 
	with open("/ebio/abt1_share/small_projects/sthavapalan/results/AroK/cluster_info/accessions_fusion.txt", "w") as file:
		for x in zip(*c):
			file.write("{0}\t{1}\n".format(*x))
	
	
	
	fasta_file = []
	
	print("Neighbours:")
	for i in flanking_proteins:
		
		time.sleep(3)
		main_fasta = Entrez.efetch(db=entrezDBName, id=i[0], rettype='fasta', api_key=API_KEY).read() # include outfiles into parser, make a function out of it
		main_str = main_fasta.rstrip()
		print("Getting fasta for: ",i[0])
		neighbour = Entrez.efetch(db=entrezDBName, id=i[1], rettype='fasta', api_key=API_KEY).read()
		print("Getting fasta for: ",i[1])
		#print(main,neighbour)
		fasta_str = neighbour.split("]")
		main_str += fasta_str[1].rstrip()
		fasta_file.append(main_str)
		
	print("\nFusions:")	
	for entry in fusion_accession:
		time.sleep(3)
		print("Getting fasta for: ", entry)
		fusion_fasta = Entrez.efetch(db=entrezDBName, id=entry, rettype='fasta', api_key=API_KEY).read()
		fasta_file.append(fusion_fasta.rstrip())
		
	with open("/ebio/abt1_share/small_projects/sthavapalan/results/AroB/cluster_info/cluster1.fasta", "w") as outfile:
		outfile.write("\n".join(fasta_file))
	
	with open("/ebio/abt1_share/small_projects/sthavapalan/results/AroE/cluster_info/accession_neighbours.txt", "w") as outfile:
		outfile.write('\n'.join('{} {}'.format(x[0],x[2]) for x in flanking_proteins))
	
		





if __name__ == "__main__":
    main()	



