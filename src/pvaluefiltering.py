import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Set the GENE_SELECTED and load data.')
parser.add_argument('-G',dest='GENE_SELECTED', type=str, help='The selected gene name')
#parser.add_argument('--LOG',type=float, default=0.05, help='The log threshold value')
#parser.add_argument('-T',dest='THRESHOLD PVALUE', type=float, default=0.05, help='The threshold value')


args = parser.parse_args()

GENE_SELECTED = args.GENE_SELECTED
#GENE_SELECTED = "DIS3L2"
FILENAME = "../output/correlations/correlations_" + GENE_SELECTED + "_by_chr.csv"
correlations = pd.read_csv(FILENAME)
threshold = 0.05#-np.log10(args.LOG)

filtered_correlations = correlations[correlations['PVALUE'] <= threshold]
print("{} significant genes".format(len(filtered_correlations)))
# Specify the output filename for the filtered data
filtered_filename = "../output/correlations/filtered/filtered_correlations_" + GENE_SELECTED + "_by_chr.csv"

# Write the filtered data to a CSV file
filtered_correlations.to_csv(filtered_filename, index=False)

# Display the first few rows of the filtered data
print('Filtered file: ' + filtered_filename)