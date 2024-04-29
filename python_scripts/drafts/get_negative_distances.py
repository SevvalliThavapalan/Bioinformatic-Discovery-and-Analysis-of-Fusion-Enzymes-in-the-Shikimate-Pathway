import pandas as pd

results = pd.read_csv('../results/08_07_21/kinase_10_operon_with_Distances.tsv', sep='\t', engine='python', error_bad_lines=False)
query = []
prot_accession = []
start_genome = []
gen_accession=[]
dist = []
for i in range(len(results)):
	if results.loc[i,"Distances"] < 0:
		query.append(results.loc[i,"Query"])
		prot_accession.append(results.loc[i,"Prot_accession_and_order_of_output"])
		start_genome.append(results.loc[i,"Start_in_genome"])
		gen_accession.append(results.loc[i,"Genome_accession"])
		dist.append(results.loc[i,"Distances"])

negative_distances = pd.DataFrame({
		   "Query": query ,
		   "Accession": prot_accession,
		   "Start_in_genome": start_genome,
		   "Genome_accession" : gen_accession,
		   "Distances": dist,
	})


out_file = negative_distances.to_csv('../results/08_07_21/negative_distances.tsv', sep='\t', index = False)
