import pandas as pd
import sys
import argparse
import matplotlib.pyplot as plt
import os.path
import numpy as np
from pathlib import Path
from collections import Counter

__author__="Sevvalli Thavapalan"

"""
Abundances and distances of all neighbors are calculated and printed. Plot with distances of specific enzymes can be 
created needs to adjusted in the code if needed. Else all distances of enzymes are plotted.
"""

def get_files():
	parser = argparse.ArgumentParser(description='Gets overall abundances of flanking enzymes from Right_and_left_distances.tsv ')
	parser.add_argument('-i', '--input', help="Path to Right_and_left_distances.tsv which was created using calc_stat.py", required=True, nargs='+')
	parser.add_argument('-o', '--out', help="Path to output files", required=True)
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments



def create_histogram(enzyme, title, out_path):
	print("Creating histogram...")
	fig = plt.figure()
	
	plt.hist(enzyme, bins=np.arange(min(enzyme), max(enzyme) + 5), color='#1f77b4') # 50
	plt.title(title)
	plt.xlabel("Distance [bp]")
	plt.ylabel("Frequency")
	fig.savefig(str(out_path) + "/" + title + ".png", bbox_inches='tight')
	
	plt.hist(enzyme, bins=np.arange(min(enzyme), max(enzyme) + 5), color='#1f77b4') # 50
	plt.title(title)
	plt.xlabel("Distance [bp]")
	plt.xlim([-100,200])
	plt.ylabel("Frequency")
	fig.savefig(str(out_path) + "/" + title + "_zoom1.png", bbox_inches='tight')

	#plt.hist(enzyme, bins=np.arange(min(enzyme), max(enzyme) + 5), color='#1f77b4') 
	#plt.title(title)
	#plt.xlabel("Distance [bp]")
	#plt.xlim([-50,100])
	#plt.ylabel("Frequency")
	#fig.savefig(str(out_path) + "/" + title + "_zoom2.png", bbox_inches='tight')

def small_frame(enzyme, title, out_path):
	fig = plt.figure()
	plt.hist(enzyme, bins=np.arange(min(enzyme), max(enzyme) + 1), color='#1f77b4') # 50
	plt.title(title)
	plt.xlabel("Distance [bp]")
	plt.xlim([-20,20])
	plt.ylabel("Frequency")
	fig.savefig(str(out_path) + "/" + title + ".png", bbox_inches='tight')



def get_abundances(enzyme_dict, name):
	list_merged_enzyme = []
	for i in enzyme_dict:
		if i.find(name) !=-1:
			list_merged_enzyme.extend(enzyme_dict[i])

	return list_merged_enzyme


def sort_dict(enzyme_dict):
	# sort the dictionary after getting the frequency of each enzyme
	len_dict= {key: len(value) for key, value in enzyme_dict.items()}
	sorted_dict = sorted(len_dict.items(), key=lambda kv: kv[1], reverse=True)
	
	# filter out the most frequent enzymes for each side; all above 50
	#print(len(sorted_dict))
	abundances = []
	for k in range(len(sorted_dict)):
		if sorted_dict[k][1] > 1:
			abundances.append(sorted_dict[k])
	return abundances, sorted_dict
	
	
def main():
	infiles = get_files()
	enzymes_right = {}
	enzymes_left = {}
	enzymes_right_accessions = {}
	enzymes_left_accessions = {}
	out_path = Path(infiles['out'])
	total = 0
	#print(infiles['input'])
	for file in infiles['input']:
		results = pd.read_csv(file, engine='python', sep='\t')
		total += results.shape[0]
		# get unique enzymes on left and right side
		unique_enzymes_right = results["Enzyme_right"].unique()
		unique_enzymes_left = results["Enzyme_left"].unique()
		# loop through each entry and get the calculated distance
		for entry in unique_enzymes_right:
			i = results[results.Enzyme_right==entry]
			if str(entry) != "nan":
				for j in range(len(i)):
					#check if enzyme is already in the dict; distances are added as a list to the dictionary
					if entry in enzymes_right:
						enzymes_right[entry].append(i.iloc[j]["Distance_right"])
					else:
						enzymes_right[entry] = []
						enzymes_right[entry].append(i.iloc[j]["Distance_right"])
					if entry in enzymes_right_accessions:
						enzymes_right_accessions[entry].append((i.iloc[j]["Distance_right"],i.iloc[j]["Accession_right"]))
					else:
						enzymes_right_accessions[entry] = []
						enzymes_right_accessions[entry].append((i.iloc[j]["Distance_right"],i.iloc[j]["Accession_right"]))
		
		# same for the left side 
		for entry in unique_enzymes_left:
			i = results[results.Enzyme_left==entry]
			if str(entry) != "nan":
				for j in range(len(i)):
					if entry in enzymes_left:
						enzymes_left[entry].append(i.iloc[j]["Distance_left"])
					else:
						enzymes_left[entry] = []
						enzymes_left[entry].append(i.iloc[j]["Distance_left"])
					if entry in enzymes_left_accessions:
						enzymes_left_accessions[entry].append((i.iloc[j]["Distance_left"],i.iloc[j]["Accession_left"]))
					else:
						enzymes_left_accessions[entry] = []
						enzymes_left_accessions[entry].append((i.iloc[j]["Distance_left"],i.iloc[j]["Accession_left"]))
	print("Total number of queries: ", total)

	print("Right Enzymes:\n")
	right_abundant, sorted_right = sort_dict(enzymes_right)		
	for name in right_abundant:
		print(name)
		title = name[0].replace(" ", "_")
		plot = create_histogram(enzymes_right[name[0]],title, out_path )
	
	# Some have different names but are the same proteins; manually look them up with 

	
	print("\n\nLeft Enzymes:\n")
	left_abundant, sorted_left = sort_dict(enzymes_left)		
	for name in left_abundant:
		print(name)
		title = name[0].replace(" ", "_")
		plot = create_histogram(enzymes_left[name[0]],title, out_path)


	
	#concat1 = small_frame(list_3DS, '3-dehydroquinate_synthase_all', out_path)
	#concat = small_frame(list_3DD, '3-dehydroquinate_dehydratase_all', out_path)
	#concat3 = small_frame(list_chor, 'chorismate_synthase_all', out_path)
	#concat5 = small_frame(list_aroa, 'EPSP_synthase_all', out_path)
	#concat4 = small_frame(list_sd, 'shikimate_dehydrogenase_all', out_path)




if __name__ == "__main__":
    main()
