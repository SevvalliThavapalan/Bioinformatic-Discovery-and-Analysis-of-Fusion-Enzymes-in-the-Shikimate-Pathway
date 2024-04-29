import pandas as pd
import argparse
import re
from pathlib import Path
import sys
from collections import Counter


def get_files():
	parser = argparse.ArgumentParser(description='Plots histograms using the information from Right_and_left_distances.tsv ')
	parser.add_argument('-i', '--input', help="Path to Right_and_left_distances.tsv which was created using calc_stat.py", required=True, nargs='+')
	parser.add_argument('-o', '--out', help="Path to output files", required=True)
	parser.add_argument('-e', '--enzyme', required=True, type=str)
	parser.add_argument('-s', '--side', required=True, type=str)
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments


def write_table(list_enzyme, name, out_path):
	counted_list= Counter(list_enzyme)
	sorted_list = sorted(counted_list.items(), key=lambda kv: kv[1], reverse=True)
	orig_stdout = sys.stdout
	abundances_file =  open(str(out_path) + '/' + name + '_abundances.txt', 'w')
	sys.stdout = abundances_file
	print(name, '\nTotal: ', len(list_enzyme))
	print('Distance\tAbundance')
	for i in sorted_list:
		print(i[0], "\t",i[1])
	sys.stdout = orig_stdout
	abundances_file.close()

def write_accession_file(list_accessions, name, out_path):
	orig_stdout = sys.stdout
	f =  open(str(out_path) + '/' + name + '.txt', 'w')
	sys.stdout = f	
	for i in list_accessions:
		print(i)
	sys.stdout = orig_stdout
	f.close()
	

def main():
	infiles = get_files()
	out_path = Path(infiles['out'])
	neighbour = infiles['enzyme']
	accessions = []
	abundances = []
	for file in infiles['input']:
		results = pd.read_csv(file, engine='python', sep='\t')
		if infiles['side'] == 'r':	 
			unique_enzymes_right = results["Enzyme_right"].unique()
			for i in unique_enzymes_right:
				enzyme = results[results.Enzyme_right==i]
				if str(i) != "nan":
					for j in range(len(enzyme)):
						if neighbour in i:
							#if enzyme.iloc[j]["Distance_left"] ==-4:
							query = enzyme.iloc[j]["Query"]
							sp = query.split('#')
							#print(sp[0])
							accessions.append(sp[0])
							abundances.append(enzyme.iloc[j]["Distance_right"])
		if infiles['side'] == 'l':	 
			unique_enzymes_left = results["Enzyme_left"].unique()
			for i in unique_enzymes_left:
				enzyme = results[results.Enzyme_left==i]
				if str(i) != "nan":
					for j in range(len(enzyme)):
						if neighbour in i:
							#if enzyme.iloc[j]["Distance_left"] ==-4:
							query = enzyme.iloc[j]["Query"]
							sp = query.split('#')
							#print(sp[0])
							accessions.append(sp[0])
							abundances.append(enzyme.iloc[j]["Distance_left"])
	
	name = re.sub(r"\s+", "_",neighbour)
	print(neighbour,":", len(accessions))

	abundance_file = write_table(abundances, name, out_path)
	#accession_file = write_accession_file(accessions, name, out_path)


if __name__ == "__main__":
    main()






			
