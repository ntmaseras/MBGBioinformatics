#!/bin/bash
#SBATCH --partition normal
#SBATCH --account RNAMetabolism
#SBATCH --mem 80g
#SBATCH -t 10:00:00
#SBATCH --cpus-per-task 4 

# Define paths
GENE_LIST="data/testgenes.txt"
LOG_FILE="log/execution.log"
COMPLETED_FILE="log/completed_genes.txt"

# Check if the log file exists; if not, create it
touch "$LOG_FILE"

# Create the completed genes file if it doesn't exist
touch "$COMPLETED_FILE"

# Process genes from the list
while read -r GENE_NAME; do
    # Check if the gene has already been processed
    if grep -q "$GENE_NAME" "$COMPLETED_FILE"; then
        echo "Gene $GENE_NAME already processed. Skipping."
    else
        echo "Running script for gene: $GENE_NAME"
        python3 src/correlations.py "$GENE_NAME"
        
        # Log the processed gene
        echo "$GENE_NAME" >> "$COMPLETED_FILE"
        
        echo "Finished script for gene: $GENE_NAME"
    fi
done < "$GENE_LIST"