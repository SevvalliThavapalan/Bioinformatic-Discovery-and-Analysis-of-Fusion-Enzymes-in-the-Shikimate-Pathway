from Bio import Entrez
import time
from pathlib import Path
import argparse

__author__="Sevvalli Thavapalan"

"""
Takes a file containing protein accessions and retrieves the fasta files
"""

def get_files():
	parser = argparse.ArgumentParser(description='Takes a txt file containing protein accessions and will get the corresponding fasta file.')
	parser.add_argument('-i', '--input', help=".txt file with accessions", required=True, nargs="+")
	parser.add_argument('-o', '--output', help='path to outfile', required=True , nargs="+")
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments


def main():
	print(__author__)
	infiles = get_files()
	out_path = infiles['output']
	entrezDBName = 'protein'
	Entrez.email = 'sevvalli.thavapalan@student.uni-tuebingen.de'
	API_KEY = '76340e02587e2d5027b1cb302cf0b11fc808' 
	lines = []
	for files, out in zip(infiles['input'], out_path):
		with open(files) as file:
			print('Getting accessions')
			lines = file.readlines()
		file.close()
		for i in lines:
			print(i[:-1])
		print("Lines read")
		unique_lines = list(set(lines))
		print(len(unique_lines))
		with open(str(out), 'w') as f:
			print("Writing to file... ")
			for i in unique_lines:
				time.sleep(3)
				entryData = Entrez.efetch(db=entrezDBName, id=i[:-1], rettype='fasta', api_key=API_KEY).read()
				print("Getting fasta for: ",i)
				f.write(entryData)
		f.close()

if __name__ == "__main__":
    main()	
