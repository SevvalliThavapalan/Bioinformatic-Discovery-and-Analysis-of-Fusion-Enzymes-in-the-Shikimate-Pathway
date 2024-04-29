import pandas as pd
import argparse
from pathlib import Path
import sys
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import pprint
import numpy as np

__author__="Sevvalli Thavapalan"

def parse_arguments():
	parser = argparse.ArgumentParser(description='Takes distances from neighbors and fusions inside a cluster from clans and writes all relevant information into files. List for PyMOL is alos included.')
	parser.add_argument('-i', '--input', help="Right and left distances", required=True, nargs='+')
	parser.add_argument('-f', '--fusions', required=True, help='Linker and length info of fusions', nargs='+')
	parser.add_argument('-c', '--cluster_accessions', required=True, help='', nargs='+')
	parser.add_argument('-a', '--accessions', required=True, help = "Accessions of the flanking protein")
	parser.add_argument('-s', '--side', required=True, type=str)
	#parser.add_argument('-o', '--out', help="Path to output files", required=True)

	args = parser.parse_args()
	arguments = args.__dict__
	return arguments


def read_accessions(accession_file):
	accession_list = []
	with open(accession_file) as file:
		accession_list = file.readlines()
	file.close()

	return accession_list


def main():
	distances = []
	linker = []
	all_accessions = []
	final_cluster_accessions = []
	final_flags_accessions = []
	infiles = parse_arguments()
	cluster_accessions = []
	# Flanking proteins 

	# Accessions obtained from cluster 
	for i in infiles['cluster_accessions']:
		cluster_accessions.extend(read_accessions(i))
	print("Cluster accessions read")

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
	#print(len(final_accessions))



	for file in infiles['input']:
		print("Determining distances")
		results = pd.read_csv(file, engine='python', sep='\t')

		for line in final_accessions:
			i = results[results["Query"].str.contains(line)]
			if infiles["side"] == "r":
				if len(i) != 0:
					distances.extend(i["Distance_right"].values)
			else:
				if len(i) != 0:
					distances.extend(i["Distance_left"].values)


	counted_distances = Counter(distances)
	print(counted_distances)
	

	fusion_accession = []
	# Fusion enzymes
	for file in infiles['fusions']:
		fusions = pd.read_csv(file, engine='python', sep='\t', header=None)
	
		for lines in final_cluster_accessions:
		#print(lines[0:12])
			i = fusions[fusions.iloc[:,0].str.contains(lines[0:12])].values
			if len(i) != 0:
				linker.append(i[0][5]*3) #nucleotides
				fusion_accession.append(i[0][0])

	counted_linker = Counter(linker)
	print(fusion_accession)	
	fusion_pymol = []

	# hardcoded paths should be changed...
	
	for accession in fusion_accession:
		fusion_pymol.append("/ebio/global2/mhartmann/SK-pathway/structures/onB/"+accession + "*ed_0*.pdb") # add this to parser
	
	c = [fusion_accession, linker] 
	with open("/ebio/abt1_share/small_projects/sthavapalan/results/AroB/cluster_info/accessions_fusion.txt", "w") as file:
		for x in zip(*c):
			file.write("{0}\t{1}\n".format(*x))
	
	with open("/ebio/abt1_share/small_projects/sthavapalan/results/AroB/cluster_info/Distances.txt", "w") as f:
		f.write(pprint.pformat(counted_distances))
	
	with open("/ebio/abt1_share/small_projects/sthavapalan/results/AroB/cluster_info/Linker.txt", "w") as f:
		f.write(pprint.pformat(counted_linker))
	
	with open("/ebio/abt1_share/small_projects/sthavapalan/results/AroB/cluster_info/Accessions_fusion.txt", "w") as outfile:
		outfile.write("\n".join(fusion_accession))
	
	with open("/ebio/abt1_share/small_projects/sthavapalan/results/AroB/cluster_info/Accessions_pymol.txt", "w") as outfile:
		outfile.write("\n".join(fusion_pymol))
	
		
	
	plt.bar(list(counted_distances.keys()), counted_distances.values(), label="Neighbours")
	plt.bar(list(counted_linker.keys()), counted_linker.values(), label= "Fusions", color='C1', alpha=0.8)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.title("AroB Fusions (BA, BE, KB, BQ)",fontsize=20)
	plt.ylabel("Abundance", fontsize=20)
	plt.xlabel("Length [bp]", fontsize=20)
	plt.legend()
	plt.show()
	
	#print(Counter(distances))
					

if __name__ == "__main__":
    main()	
	
