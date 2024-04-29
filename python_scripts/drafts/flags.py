import pandas as pd
results = pd.read_csv("kinase_20000.tsv", sep='\t', engine='python', comment = '#', header=None)
queries = results.loc[:,0].unique()
	#print("Queries: ")
	#print(queries)
for i in queries:
	accessions = results.loc[results[0]==i,3]
		#print(accessions)
	out_file = accessions.to_csv('kinase'+ i +'_hits.txt', index = False, header = None, sep='\n')
	

nohup python3 /ebio/abt1_share/small_projects/sthavapalan/FlaGs/FlaGs-FlaGs_LatestVersion/FlaGs.py -p kinase_hits_0_hits.txt -g 8 -u sevvalli.thavapalan@tuebingen.mpg.de -c 4 -k -o kinase_1/kinase_1 > kinase_1/kinase1.log 2>&1&




