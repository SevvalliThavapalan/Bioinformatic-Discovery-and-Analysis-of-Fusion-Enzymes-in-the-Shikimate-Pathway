import pandas as pd
import argparse
import re
from pathlib import Path
import sys
from collections import Counter, defaultdict

__author__="Sevvalli Thavapalan"

"""
Individual abundances are written to files. Enzyme name and side need to be specified.
"""

def get_files():
	parser = argparse.ArgumentParser(description='Determines individual abundances of enzymes and writes the accessions in to files ')
	parser.add_argument('-i', '--input', help="Path to Right_and_left_distances.tsv which was created using calc_stat.py", required=True, nargs='+')
	parser.add_argument('-o', '--out', help="Path to output files", required=True)
	parser.add_argument('-e', '--enzyme', required=True, type=str, nargs='+')
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
	neighbours = infiles['enzyme']
	accessions = []
	accessions_all = []
	abundances = []
	flanking_accessions = []
	for file in infiles['input']:
		results = pd.read_csv(file, engine='python', sep='\t') # This part could be put into a function; looks nicer :)
		for neighbour in neighbours:
		
			if infiles['side'] == 'r':	 
				unique_enzymes_right = results["Enzyme_right"].unique()
				for i in unique_enzymes_right:
					enzyme = results[results.Enzyme_right==i]
					if str(i) != "nan":
						for j in range(len(enzyme)):
							if neighbour == i or "MULTISPECIES: " + neighbour in i:  #bei abundances; bei Accessions mit Leerzeichen...
							#if neighbour in i:
								query = enzyme.iloc[j]["Query"]
								flanking_accessions.append(enzyme.iloc[j]["Accession_right"])
								sp = query.split('#')
								accessions_all.append(sp[0])
								accessions.append((enzyme.iloc[j]["Distance_right"],sp[0]))
								abundances.append(enzyme.iloc[j]["Distance_right"])
			if infiles['side'] == 'l':	 
				unique_enzymes_left = results["Enzyme_left"].unique()
				for i in unique_enzymes_left:
					enzyme = results[results.Enzyme_left==i]
					if str(i) != "nan":
						for j in range(len(enzyme)):
							if neighbour == i or "MULTISPECIES: " + neighbour in i:
							#if neighbour in i:
								query = enzyme.iloc[j]["Query"]
								flanking_accessions.append(enzyme.iloc[j]["Accession_left"])
								sp = query.split('#')
								accessions_all.append(sp[0])
								accessions.append((enzyme.iloc[j]["Distance_left"],sp[0]))
								abundances.append(enzyme.iloc[j]["Distance_left"])
	
	name1 = re.sub(r"\s+", "_",neighbours[0]) #replace spaces with underscore before saving the file
	name = name1.replace('/', ' ')
	print(neighbours[0],":", len(accessions))
	print(len(flanking_accessions))
	#abundance_file = write_table(abundances, name, out_path)
	accession_file = write_accession_file(accessions_all, name, out_path)
	accession_dict = defaultdict(list)
	for i, j in accessions:
		accession_dict[i].append(j)
	
	for key, value in accession_dict.items():
		if len(accession_dict[key]) >= 10:
			print(key, len(accession_dict[key]))
			title= name + str(key)
			distance_file = write_accession_file(accession_dict[key], title, out_path)
	
	
"""If accessions of flanking proteins are needed"""

	#unique_accessions = list(set(flanking_accessions))
	#print(len(unique_accessions))
	#with open(str(out_path), 'w') as w:
	#	for item in flanking_accessions :
	#		w.write("%s\n" % item)
	
	#w.close()



if __name__ == "__main__":
    main()






			
