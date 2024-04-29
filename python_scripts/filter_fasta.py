from Bio import SeqIO
import pandas as pd
import argparse
from pathlib import Path

__author__="Sevvalli Thavapalan"

"""
A given fasta file is filtered by length. Sequence length hs to be manually added 
in the script at line 32
"""

def get_file():
	parser = argparse.ArgumentParser(description='A given fasta file is filtered by lenght.')
	parser.add_argument('-i', '--input', help="fastafile", required=True) 
	parser.add_argument('-o', '--output', help='path and name to outfile', required=True)
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments


def main():
	infiles = get_file()
	out_path = Path(infiles['output'])
	filtered_accessions = []
	FastaFile = open(infiles['input'], 'rU')
	for rec in SeqIO.parse(FastaFile, 'fasta'):
		name = rec.id
		seq = rec.seq
		seqLen = len(rec)
		if seqLen > 430: # and seqLen < 180: # try adding this to the parser; more convinient
			#print(seqLen)
			filtered_accessions.append(name)

	FastaFile.close()
	with open(str(out_path), 'w') as w:
		for item in filtered_accessions :
			w.write("%s\n" % item)
	w.close()

if __name__ == "__main__":
    main()	

