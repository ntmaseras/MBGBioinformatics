import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Set the GENE_SELECTED and load data.')
parser.add_argument('GENE_SELECTED', type=str, help='The selected gene name')
args = parser.parse_args()

GENE_SELECTED = args.GENE_SELECTED
#GENE_SELECTED = "DIS3L2"
FILENAME = "../output/correlations/correlations_" + GENE_SELECTED + "_by_chr.csv"
correlations = pd.read_csv(FILENAME)
threshold = 0.05

filtered_correlations = correlations[correlations['LOGPVALUE'] < threshold]

# Specify the output filename for the filtered data
filtered_filename = "../output/correlations/filtered/filtered_correlations_" + GENE_SELECTED + "_by_chr.csv"

# Write the filtered data to a CSV file
filtered_correlations.to_csv(filtered_filename, index=False)

# Display the first few rows of the filtered data
print('Filtered file: ' + filtered_filename)