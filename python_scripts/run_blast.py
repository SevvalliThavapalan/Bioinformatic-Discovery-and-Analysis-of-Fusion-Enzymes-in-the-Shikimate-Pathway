import argparse
import sys
import pandas as pd
from pathlib import Path
import glob, os
import time

__author__="Sevvalli Thavapalan"

# use nohup ... > .log 2>&1& to redirect
# paths need to be adjusted if the location of the script is changed
# also check if the right flaGs file is used (two versions, linux and normal)

def parse_commandline() -> dict:
	"""
	Saves commandline arguments for a blastp call in a dictionary. Most values are set to recommended values
	from the blast documentation. Blastp version: 2.11
	"""
	
	
	parser = argparse.ArgumentParser(description='Blastp search for given queries is performed. Most parameter are set to default but can be changed if neccessary. The Blast output ist then used to run FlaGs using the protein accessions of the hits. Parameters for FlaGs are hardcoded and need to be changed inside the script if neccessary.')
	# Required parameters; need to be provided
	parser.add_argument('-query', help='Query file, should contain accession number or Fasta without gi', required=True)
	parser.add_argument('-db', help='Path to database', required=True) #custom  sequences as db, needs to be adjusted if needed
	parser.add_argument('-out', help='Path to output file', required=True)
	parser.add_argument('-evalue', default=0.001, help='E value threshold for hits; default: 0.001', required = True)
	parser.add_argument('-max_target_seqs', default = 10, help='Maximum number of aligned sequences')
	# Optional parameters, can be adjusted from the user
	parser.add_argument('-task', default='blastp', choices=['blastp', 'blastp-fast', 'blastp-short'], help='Task to execute; default blastp; default values are adjusted for use of blastp.') # only for blastp
	parser.add_argument('-outfmt', default='"7 qacc saccver slen pident nident sstart send evalue"', help='Formatting options; Look at Blast documentation for details. These results are created to work with FlaGs. Needs to be put into a string')
	parser.add_argument('-word_size', default=3, help='Word size for wordfinder algorithm; default: 3')
	parser.add_argument('-gapopen', default=11, help='Cost to open a gap') # gap costs are adjusted to the matrix
	parser.add_argument('-gapextend', default=1, help='Cost to extend a gap')
	parser.add_argument('-matrix', default='BLOSUM62', help='Scoring matrix; default BLOSUM62, recommended for Proteins')
	parser.add_argument('-threshold', default=11, help='Minimum score to add to BLAST lookup table; default = 11')
	parser.add_argument('-num_threads', default=16, help='Number of threads which should be used; default = 16')
	parser.add_argument('-comp_based_stats', default=2, choices=['D','d','0','F','f','1','2','T','t','3'], help='Use composition based statistics; ; D and d = default = 2; 0,F and f = no composition based statistics; 1 refers to  NAR 29:2994-3005, 2001, 2 refers to' 							 'Bioinformatics 21:902-911 and 3 to 21:902-911, 2005, unconditionally')
	parser.add_argument('-soft_masking', default=False, help='Apply filtering locations as soft masks; default = False')
	parser.add_argument('-xdrop_gap_final', default=25, help='Heuristic value for final gapped alignment; default = 25' )
	parser.add_argument('-window_size', default=40, help='Multiple hits window size, use 0 to specify 1-hit algorithm; default = 40')
	#parser.add_argument('-num_iterations', default=3) #only for psiblast
	
	
	args = parser.parse_args()
	
	# save arguments into a dictionary
	arguments = args.__dict__
	
	return arguments


def convert_dict_to_string(arguments: dict) -> str:
	"""
	Converts the dictionary containing the commandline arguments into a string.
	This string is nedded for the BLASTP call.
	"""
	dict_string = ""
	for key, value in arguments.items():
		dict_string +=  "-"+ key + " " + str(value) + " " 
	
	return dict_string


def parse_blast_results(arguments):
	"""
	The Blastp output is parsed and the refseq accessions of the hits are extracted and 
	saved into a new file. A seperate file for each query is created for the FlagS input.
	"""

	results = pd.read_csv(arguments['out'], sep='\t', engine='python', comment = '#', header=None)
	out_path = Path(arguments['out'])
	queries = results.loc[:,0].unique()
	#print("Queries: ")
	#print(queries)
	list_acc =  [] #modified
	for i in queries:
		list_acc.extend(results.loc[results[0]==i,1].values)
		#print(accessions)
	
	with open(str(out_path.parent) + '/cluster_hits.txt', "w") as outfile:
		outfile.write("\n".join(list_acc))
	
	
	#out_file = list_acc.to_csv(str(out_path.parent) + '/ cluster1_hits.txt', index = False, header = None, sep='\n')			
	#print('Protein accessions saved in ' + str(out_path.parent) + '/')
	#return out_path


def run_flagS(out_path):
	"""
	Takes the _hits.txt files as input and performs the flanking gene analysis using FlaGs.
	A seperate folder for each query (from blast) is created and the hits file is moved in there too.
	"""
	os.chdir(out_path.parent)
	for file in glob.glob("*_hits.txt"):
		prefix = file[:-9]
		os.system('mkdir ' + prefix)
		flags_dict =  {'e': 1e-10, 'n': 3, 'g': 8, 'u': 'sevvalli.thavapalan@tuebingen.mpg.de', 'c': 16}
		flags_dict['p'] = prefix + '.txt'
		flags_dict['o']=  prefix + '/' + prefix
		flags_string = convert_dict_to_string(flags_dict)
		print(flags_string)
		os.system('python3 FlaGs.py ' +  flags_string[:-1]) #  
		os.system('mv ' + file + ' ' + prefix + '/' + file)

	
	
def main():
	
	t1 = time.perf_counter()
	arguments = parse_commandline()
	arguments_string = convert_dict_to_string(arguments)
	print('blastp ' + arguments_string)
	os.system('/ebio/abt1_toolkit/bioprogs/tools/blastplus/bin/blastp ' + arguments_string[:-1])
	print('Blast results saved in ' + arguments['out'])
	accession_list = parse_blast_results(arguments)
	#flags = run_flagS(accession_list) # 
	t2 = time.perf_counter()
	print('Total time in s:')
	print(t2-t1)

if __name__ == "__main__":
    main()
	
