import pandas as pd
import argparse
from pathlib import Path

__author__="Sevvalli Thavapalan"

def get_files():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-i', '--input', help="Path to Right_and_left_distances.tsv which was created using calc_stat.py", required=True, nargs='+')
	parser.add_argument('-o', '--out', help="Path to output file as txt", required=True)
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments


def main():
	infiles = get_files()
	out_path = Path(infiles['out'])
	accessions = []
	
	for file in infiles['input']:
		results = pd.read_csv(file, engine='python', sep='\t')
		unique_enzymes_right = results["Enzyme_right"].unique()
		for i in unique_enzymes_right:
			enzyme = results[results.Enzyme_right==i]
			if str(i) != "nan":
				for j in range(len(enzyme)):
							#if neighbour == i or "MULTISPECIES: " + neighbour == i:  #bei abundances; bei Accessions mit Leerzeichen...
					if "hypothetical protein" in i:
						#abundances.append(enzyme.iloc[j]["Distance_right"])
						print(enzyme.iloc[i]["Accessions_right"])
						
							
							
								#unique_enzymes_right = results["Enzyme_right"].unique()
		#for i in results:
		#	if "hypothetical protein" in i["Enzyme_right"]:
		#		accessions.extend(i["Accession_right"].values.tolist())
		
		#accessions.extend(i["Accession_left"].values.tolist())
		#accessions.extend(j["Accession_left"].values.tolist())
	#print(len(accessions))
	#unique_accessions = list(set(accessions))
	#print(len(unique_accessions))
	#with open(str(out_path), 'w') as w:
	#	for item in unique_accessions :
	#		w.write("%s\n" % item)
	#w.close()



if __name__ == "__main__":
    main()
	
