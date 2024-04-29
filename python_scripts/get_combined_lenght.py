import pandas as pd
import argparse
from pathlib import Path
import sys
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import pprint

__author__="Sevvalli Thavapalan"

"""
Gets information about linker an flanking genes and plots them together into a barchart.
Files with information about both type of distances can be saved if needed.
Code needs to be adjusted.
"""


def get_files():
	parser = argparse.ArgumentParser(description='Gets information about linker an flanking genes and plots them together into a barchart.')
	parser.add_argument('-i', '--input', help="operon distances", required=True, nargs='+')
	parser.add_argument('-a', '--accessions', required=True, help = "Accessions of the flanking protein")
	parser.add_argument('-f', '--fusions', required=True, help="File with information of the fusion enzymes from get_linker_lengths.py")
	#parser.add_argument('-o', '--out', help="Path to output files", required=True)
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments
	


def main():
	lines= []
	length_combined = []
	linker_combined = []
	infiles = get_files()
	with open(infiles['accessions']) as file:
		lines = file.readlines()
	file.close()
	print("Accessions read")
	final_lines = []
	for line in lines:
		final_lines.append(line.rstrip("\n"))


	fusions = pd.read_csv(infiles['fusions'], engine='python', sep='\t')
	counted_length_fusions = Counter(fusions.iloc[:,8]*3)
	counted_linker_fusions = Counter(fusions.iloc[:,5]*3)
	
	for file in infiles['input']:
		results = pd.read_csv(file, engine='python', sep='\t')
		print(file)
		neighbour = []
		queries = []
		
		for lines in final_lines: #checks for the accessions in the dataframe
			i = results[results["Prot_accession_and_order_of_output"].str.contains(lines)].values.tolist()
			if len(i) != 0:
				query = results.loc[results["Query"] == i[0][0]]
				j = query.loc[query["Distances"] == 0].values.tolist() # returns a nested list of each row

				if i[0][12] < 0:
					combined_length = i[0][1]+j[0][1]-i[0][12]
					neighbour.append((j[0][1],i[0][0],i[0][1],i[0][9],i[0][12], combined_length))
					# len of query enzyme, accession of query, len of neighbouring enzyme,accession, and distance
				else:
					combined_length = i[0][1]+j[0][1] + i[0][12]
					neighbour.append((j[0][1],i[0][0],i[0][1],i[0][9],i[0][12], combined_length))
		
		for entry in neighbour:
			#print(entry)
			length_combined.append(entry[5])
			linker_combined.append(entry[4])


	counted_length_combined= Counter(length_combined)
	counted_linker_combined = Counter(linker_combined)
	
	
	# files can be created using the next lines if needed; Overview with plots was enough

	#with open("test.txt", "w") as f:
	#	f.write(pprint.pformat(counted_length, indent=4))
	#with open("test2.txt", "w") as f:
		#f.write(pprint.pformat(counted_linker, indent=4))
		    
	plt.bar(list(counted_length_combined.keys()), counted_length_combined.values(), label="Neighbours")
	plt.bar(list(counted_length_fusions.keys()), counted_length_fusions.values(), label= "Fusions")
	plt.title("AroB + AroA Length",fontsize=20)
	plt.ylabel("Abundance")
	plt.xlabel("Length")
	plt.legend()
	plt.show()
	
	print("Done")


if __name__ == "__main__":
    main()	
