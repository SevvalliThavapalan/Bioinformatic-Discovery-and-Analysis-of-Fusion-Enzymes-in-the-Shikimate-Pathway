import pandas as pd
import numpy as np
import argparse
from collections import Counter, defaultdict


"""
Cluster information from Clans without accessions is taken and the
individual taxonomic ranks are counted to get an overview of all 
taxa inside the given cluster. All taxa which occur more than two 
times are printed.
"""




def get_files():
	parser = argparse.ArgumentParser(description="The modified files from clans without the accessions are needed to count" \
	 											"the different taxa. Files are named _taxa.txt")
	parser.add_argument('-i', '--input', help="Taxa in one cluster", required=True)
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments


def main():
	lines = []
	infiles = get_files()
	with open(infiles['input']) as file:
		lines = file.readlines()
	file.close()
	final_lines = []
	for line in lines:
		final_lines.append(line.rstrip("\n"))
	
	all_ranks = []
	
	for line in final_lines:
		taxa = line.split(";")
		for rank in taxa:
			all_ranks.append(rank)

	counted_ranks = Counter(all_ranks)



	for key, value in counted_ranks.items():
		if value > 2:
			print(key, value)


if __name__ == "__main__":
    main()	
    
