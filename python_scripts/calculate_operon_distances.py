import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import os
import sys

__author__="Sevvalli Thavapalan"

def parse_commandline():
	"""
	Takes multiple files as input (*_operon.tsv). Calculates the distances between each query and its flanking genes.
	Distances between each falnking gene and query is added into the data frame as an additional column.
	"""
	
	parser = argparse.ArgumentParser(description='Calculates the distances between each query and its flanking genes. Output is '
										'saved into a new tsv file which includes header and the distances.')
	parser.add_argument('-i', '--input', help='Input file of FlaGs results named *_operon.tsv to calculate distances between genes.', required=True, nargs='+')
	parser.add_argument('-o', '--out', help='Path to output file called *_with_Distances.tsv.', required=True, nargs='+')
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments


def calculate_distances(df):
	"""
	For each input file the distance between query and flanking gene is calculated and added to the dataframe
	"""
	queries = df.Query.unique()
	distances = []
	#check= []
	for query in queries:
		i = df[df.Query == query]
		end_pos = 0
		pairs = []
		for j in range(len(i)):
			pairs.append((i.iloc[j]['Start_position_relative_to_query'], i.iloc[j]['End_position_relative_to_query']))
			if i.iloc[j]['Start_position_relative_to_query'] == 1:
				end_pos += i.iloc[j]['End_position_relative_to_query']
		for j in range(len(pairs)):
			if pairs[j][0] < 1:
				distances.append(0 - pairs[j][1])
			elif pairs[j][0] > 1:
				distances.append(pairs[j][0] - end_pos -1)
			else:
				distances.append(0)
	df['Distances'] = np.array(distances)
	return df


def write_out_file(df, arguments, file, out_path):
	"""
	DataFrame is written into a tsv file for further use
	"""
	name = os.path.basename(file)
	out_file = df.to_csv(str(out_path) + '/' + name[:-4] + '_with_Distances.tsv', sep='\t', index = False)
	print('outfile written to: ', out_path)
		

def main():
	arguments = parse_commandline()
	for i, o  in  zip(arguments['input'], arguments['out']):
		print("Calculating distances for: " + i)
		results = pd.read_csv(i, sep='\t', engine='python', header=None, error_bad_lines=False)
		results.columns = ['Query', 'Length_in_bp', 'Query_strand', 'Direction_to_query', 'Cluster', 
							'Start_position_relative_to_query', 'End_position_relative_to_query','Start_in_genome', 
							'End_in_genome', 'Prot_accession_and_order_of_output', 'Genome_accession', 
							'Genome_assembly_identifier']
		new_results = calculate_distances(results)
		final = write_out_file(new_results, arguments,i, Path(o))


if __name__ == "__main__":
    main()
