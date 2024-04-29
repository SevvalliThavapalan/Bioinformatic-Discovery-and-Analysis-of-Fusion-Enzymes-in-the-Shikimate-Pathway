import argparse
from pathlib import Path
import glob, os

__author__="Sevvalli Thavapalan"
#TODO use nohup ... > .log 2>&1& to redirect
#TODO paths need to be adjusted if the location of the script is changed
#TODO also check if the right flaGs file is used (two versions, linux and normal)

def get_args():
	parser = argparse.ArgumentParser(description='Splits the query list obtained using BLAST into smaller files (Default: 1000). FlaGs is used to perform flanking gene analysis.')
	parser.add_argument('-i', '--input', help="blast result file will be split", required=True)
	parser.add_argument('-o','--outpath', help='Path to output file', required=True)
	
	args = parser.parse_args()
	arguments = args.__dict__
	return arguments
	
	
def split_file(blast_result):
	'''
	Gets list of lines from the blast result and splits them into 1000 lines for each file. Number can be adjusted.
	'''	
	test = open(blast_result, 'r')
	listoflines = test.readlines()
	test.close()
	lines  = 1000
	prefix = blast_result[:-4]
	final = [listoflines[i* lines:(i+1)*lines] for i in range((len(listoflines)+lines-1) // lines)]
	for element in range(len(final)):
		textfile = open(prefix + "_" + str(element) + '_hits.txt', 'w')
		for line in final[element]:
			textfile.write(line)
		textfile.close()
		
		
def convert_dict_to_string(arguments):
	"""
	Converts the dictionary containing the commandline arguments into a string.
	This string is nedded for the BLASTP call.
	"""
	dict_string = ""
	for key, value in arguments.items():
		dict_string +=  "-"+ key + " " + str(value) + " " 
	print(dict_string)
	return dict_string

		
def run_flagS(out_path):
	"""
	Takes the _hits.txt files as input and performs the flanking gene analysis using FlaGs.
	A seperate folder for each query (from blast) is created.
	"""
	os.chdir(out_path.parent)
	for file in glob.glob(str(out_path) + "/*_hits.txt"):
		prefix = file[:-9]
		os.system('mkdir ' + prefix)
		flags_dict =  {'e': 1e-10, 'n': 3, 'g': 3, 'u': 'sevvalli.thavapalan@tuebingen.mpg.de', 'c': 4}
		flags_dict['p'] = prefix + '.txt'
		flags_dict['o']=  prefix + '/' + os.path.basename(prefix)
		flags_string = convert_dict_to_string(flags_dict)
		os.system('nohup python3 /ebio/abt1_share/small_projects/sthavapalan/FlaGs/FlaGs-FlaGs_LatestVersion/FlaGs.py ' +  flags_string[:-1]+ ' >' +prefix + '/' + os.path.basename(prefix) +'.log 2>&1&') 
		# Path needs to be adjusted if needed



def main():
	args = get_args()
	blast_result = args['input']
	out_path = Path(args['outpath'])
	splitted_files = split_file(blast_result)
	flags = run_flagS(out_path)



if __name__ == "__main__":
    main()
