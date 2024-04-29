from Bio import SeqIO
import pandas as pd
import argparse
from pathlib import Path

__author__="Sevvalli Thavapalan"

def get_file():
	parser = argparse.ArgumentParser(description='Compares a list of protein accessions to a fasta file.' 
												'The matches are written into a new file')
	parser.add_argument('-i', '--input', help="fastafile", required=True)
	parser.add_argument('-a', '--accessions', required=True)
	parser.add_argument('-o', '--out', help='path and name to outfile', required=True)

	args = parser.parse_args()
	arguments = args.__dict__
	return arguments

def main():
	infiles = get_file()
	out = Path(infiles['out'])
	fasta_info = []


	with open(infiles['accessions']) as file:
		print('Getting accessions')
		lines = file.readlines()
	file.close()
	unique_lines = list(set(lines))
	print(unique_lines)

	FastaFile = open(infiles['input'], 'rU')
	for rec in SeqIO.parse(FastaFile, 'fasta'):
		name = rec.id
		seq = rec.seq
		seqLen = len(rec)
		for accession in unique_lines:
			if accession[:-1] in name:
				fasta_info.append((rec.description,seqLen,str(seq)))
	with open(str(out), 'w') as outfile:
		outfile.write('\n'.join('>{} {}\n{}'.format(x[0],x[1],x[2]) for x in fasta_info))


if __name__ == "__main__":
    main()	
