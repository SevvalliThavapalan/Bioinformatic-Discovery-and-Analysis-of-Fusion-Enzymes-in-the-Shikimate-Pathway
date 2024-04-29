from Bio import Entrez
import time

__author__="Sevvalli Thavapalan"

def get_files():
	parser = argparse.ArgumentParser(description='Takes a txt file containing protein accessions and will get the corresponding fasta file.')
	parser.add_argument('-i', '--input', help=".txt file with accessions", required=True)
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments


def main():
	print(__author__)
	infiles = get_files()
	entrezDBName = 'protein'
	Entrez.email = 'sevvalli.thavapalan@tuebingen.mpg.de'
	lines = []
	with open(infiles['input']) as file:
		print('Getting accessions')
		lines = file.readlines()
	file.close()
	for i in lines:
		print(i[:-1])
	print("Lines read")
	with open('fusions.fasta', 'w') as f:
		print("Writing to file... ")
		for i in lines:
			entryData = Entrez.efetch(db=entrezDBName, id=i[:-1], rettype='fasta').read()
			time.sleep(2.0)
			print("Getting fasta for: ",i)
			f.write(entryData)
	f.close()

if __name__ == "__main__":
    main()	
