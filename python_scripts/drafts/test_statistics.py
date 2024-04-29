import pandas as pd

results = pd.read_csv('kinase1_operon_with_Distances.tsv', engine='python', sep='\t')
queries = results.Query.unique()
left_accessions = []
right_accessions = []
cluster_left = []
cluster_right = [] # als set ?
for query in queries:
	i = results[results.Query == query]
	
	sorted_columns = i.sort_values(by=['Distances'])
	if len(sorted_columns) == 9:
		if sorted_columns.iloc[1]['Start_position_relative_to_query'] < 1:
			accession = sorted_columns.iloc[1]['Prot_accession_and_order_of_output']
			splitted = accession.split('#')
			left_accessions.append((splitted[0], sorted_columns.iloc[1]['Distances']))
			cluster_left.append(sorted_columns.iloc[1]['Cluster'])
		else:
			accession = sorted_columns.iloc[1]['Prot_accession_and_order_of_output']
			splitted = accession.split('#')
			right_accessions.append((splitted[0], sorted_columns.iloc[1]['Distances']))
			cluster_right.append(sorted_columns.iloc[1]['Cluster'])

		if sorted_columns.iloc[2]['Start_position_relative_to_query'] < 1:
			accession = sorted_columns.iloc[2]['Prot_accession_and_order_of_output']
			splitted = accession.split('#')
			left_accessions.append((splitted[0], sorted_columns.iloc[2]['Distances']))
			cluster_left.append(sorted_columns.iloc[2]['Cluster'])
		else:
			accessions = sorted_columns.iloc[2]['Prot_accession_and_order_of_output']
			splitted = accession.split('#')
			right_accessions.append((splitted[0], sorted_columns.iloc[2]['Distances']))
			cluster_right.append(sorted_columns.iloc[2]['Cluster'])
	else:
		print(query)
		right_accessions.append((0,0))
		left_accessions.append((0,0))
clusters = []
clusters.append(set(cluster_right))
clusters.append(set(cluster_left))
accessions_right = [accession[0] for accession in right_accessions]
print(len(right_accessions))
print(len(left_accessions))
print(len(queries))
#right = pd.DataFrame({
		#"Query" : queries,
		#"Accessions_Right" : accessions_right

	#})
#print(right)
