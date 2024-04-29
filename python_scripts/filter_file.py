from Bio import Entrez
import time
from pathlib import Path
import argparse
import pandas as pd

__author__="Sevvalli Thavapalan"

"""
Cross checks accession list with Right and left distances files. Checks for the query.
Results are saved into a csv file.
"""

def get_files():
	parser = argparse.ArgumentParser(description='If results with distances are already present, you can filter by accessions.')
	parser.add_argument('-i', '--input', help="", required=True, nargs='+') #right and left distances file
	parser.add_argument('-a', '--accessions', required=True)
	parser.add_argument('-o', '--output', help='path and name to outfile, should be a .tsv', required=True)
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments
	
def main():
	lines = []
	df_list = []
	infiles = get_files()
	out_path = Path(infiles['output'])
	with open(infiles['accessions']) as file:
		lines = file.readlines()
	file.close()
	
	for files in infiles['input']:
		print("Filtering:" + files)
		results = pd.read_csv(files, engine='python', sep='\t')
		queries = results['Query']
		no_fusions = []
		for i in queries:
			for entry in lines:
				if i.find(entry.rstrip("\n")) !=-1:
					no_fusions.append(i)
					break

		print(len(no_fusions))
		df_list.append(results[results['Query'].isin(no_fusions)])
	df = pd.concat(df_list)	
	df.to_csv(str(out_path), sep='\t', index=False)

if __name__ == "__main__":
    main()	

