
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import concurrent.futures

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Set the GENE_SELECTED and load data.')
parser.add_argument('GENE_SELECTED', type=str, help='The selected gene name')
parser.add_argument('--map_file', type=str, default='../data/Achilles_guide_map.csv', help='Path to the map file')
parser.add_argument('--data_file', type=str, default='../data/20Q4v2_Achilles_gene_effect_.csv', help='Path to the data file')

args = parser.parse_args()

GENE_SELECTED = args.GENE_SELECTED
map_file = args.map_file
data_file = args.data_file

map = pd.read_csv(map_file)
achilles = pd.read_csv(data_file)

achilles = achilles.dropna()
achilles.columns = achilles.columns.str.replace(" ", "_").str.replace("[()]", "", regex=True)
achilles_clean= pd.melt(achilles, id_vars=["DepMap_ID"], var_name="genes", value_name="value_gene_selected")


genes = pd.Series(achilles_clean.genes.unique())
achilles_clean[['gene', 'gene_id']] = achilles_clean['genes'].str.split('_', expand=True)
achilles_clean.drop(columns=['genes'], inplace = True)

map[['chromosome', 'avg_pos','_']] = map['genome_alignment'].str.split('_', expand=True)
map[['gene', 'gene_id']] = map['gene'].str.split(' ', expand=True)
map.gene_id = map.gene_id.str.replace("[()]", "", regex=True)


data = achilles_clean.copy()
selected_gene_data = data[data['gene'] == GENE_SELECTED]['value_gene_selected']

# Create an empty dictionary to store correlations
correlations = []

gene_count = 0
def calculate_correlation(selected_gene, other_gene_data):
    correlation, pvalue = pearsonr(selected_gene, other_gene_data)
    return correlation, pvalue, -np.log10(pvalue)


# Get unique genes in the dataset
unique_genes = data['gene'].unique()

correlations = []
print("Starting to compute gene dependency for " + GENE_SELECTED)
# correlations for a single gene
def calculate_gene_correlations(other_gene):
    other_gene_data = data[data['gene'] == other_gene]['value_gene_selected']
    correlation, pvalue, logpvalue = calculate_correlation(selected_gene_data, other_gene_data)
    return {'SELECTED_GENE': GENE_SELECTED, 'EVALUATED_GENE': other_gene, 'CORRELATION': correlation, 'PVALUE': pvalue, 'LOGPVALUE': logpvalue}

# concurrent.futures to parallelize the computation
with concurrent.futures.ProcessPoolExecutor() as executor:
    for gene_count, result in enumerate(executor.map(calculate_gene_correlations, unique_genes), 1):
        correlations.append(result)
        
        # Print progress 
        if gene_count % 1000 == 0:
            print(f"Computed correlations for {gene_count} / {len(unique_genes)} genes.")


results = pd.DataFrame(correlations)
FILENAME = "../output/correlations/correlations_" + GENE_SELECTED + "_by_chr.csv"

# Save the results to a CSV file
results.to_csv(FILENAME, index=False)