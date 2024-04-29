import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
from pathlib import Path

__author__="Sevvalli Thavapalan"

"""
For all possible fusion enzymes which are not part of the pathway.
The proteins are filtered by their length. Sequence lengths have to be adapted 
for each enzyme manually. 
"""

def get_files():
	parser = argparse.ArgumentParser(description='Finds possible fusions inside the FlaGs results by filtering the proteins by length.')
	parser.add_argument('-i', '--input', help="operon.tsv is needed", required=True, nargs='+')
	parser.add_argument('-o', '--output', help='path to outfile', required=True)
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments



def main():
	infiles = get_files()
	genome_accession = []
	genome_pos = []
	query = []
	len_query = []
	list_prots = []
	single_enzymes = []
	new_df = pd.DataFrame()
	out_path = Path(infiles['output'])
	for file in infiles['input']:
		Distances = pd.read_csv(file, engine='python', sep='\t', error_bad_lines=False, header=0)
		queries = Distances.loc[Distances['Start_position_relative_to_query']==1]

		pot_fusions = queries.loc[queries.Length_in_aa > 750] # are given in nucleotide length; *3 to convert
		for i in pot_fusions['Prot_accession_and_order_of_output']:
			a = i.split("#") # order of output is still included in table
			list_prots.append(a[0])
		no_fusions_1 =  queries.loc[queries.Length_in_aa < 750]
		no_fusions_2 = no_fusions_1.loc[no_fusions_1.Length_in_aa > 420 ]    #TODO: specific for kinases, adapt for other enzymes
		
		for i in no_fusions_2['Prot_accession_and_order_of_output']:
			a = i.split("#") # order of output is still included in table
			single_enzymes.append(a[0])
			
			
	
	with open(str(out_path) + '/fusions.txt', 'w') as f:
		for item in list_prots:
			f.write("%s\n" % item)
	f.close()
	
	with open(str(out_path) + '/no_fusions_kinase.txt', 'w') as w:
		for item in single_enzymes:
			w.write("%s\n" % item)
	w.close()
	
	
	
	#print(len(list_prots))		


if __name__ == "__main__":
    main()	
