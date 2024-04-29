import pandas as pd
import numpy as np
import argparse

__author__="Sevvalli Thavapalan"

"""
Checks the files created by calc_stat.py for specific neighbours on each side.
The names of the enzymes have to be specified. 
"""

def parse_arguments():
	parser = argparse.ArgumentParser(description='Checks for specific right and left enzymes in given files created by calc_stat.py.')
	parser.add_argument('-i', '--input', help="Right and left distances", required=True, nargs='+')
	parser.add_argument('-r', '--right_enzyme', required=True, type=str)
	parser.add_argument('-l', '--left_enzyme', required=True, type=str)
	parser.add_argument('-o', '--out', help="Path to output files", required=True)

	args = parser.parse_args()
	arguments = args.__dict__
	return arguments


def main():
	
	infiles = parse_arguments()
	accessions = []
	
	for file in infiles['input']:
		print("Reading:" , file)
		results = pd.read_csv(file, engine='python', sep='\t')
		for i in range(len(results)) :
  			if results.loc[i, "Enzyme_right"] == infiles['right_enzyme'] and results.loc[i, "Enzyme_left"] == infiles['left_enzyme']:
  				query = results.loc[i,"Query"]
  				accessions.append(query[0:12]) # order of output from FlaGs is still present
  	
	print(len(accessions))
	if len(accessions) != 0: # in case there a no hits
		with open(infiles['out'], "w") as outfile:
			outfile.write("\n".join(accessions))
	else:
		print("Skip")

if __name__ == "__main__":
    main()	
