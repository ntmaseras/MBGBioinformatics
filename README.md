# Study of the Role of RNA Metabolism in Tumors

This project aims to establish a map of genetic interactions focused on RNA-binding proteins and RNA metabolism enzymes. The study involves a multi-step approach, starting with an extensive review of the experimental and statistical methodology used to construct gene dependency scores. It will be followed by an investigation into the use of gene-gene correlation scores to identify genetic interactions. The project will leverage data from DepMap, a genome-scale CRISPR screen, conducted on a panel of more than 800 cancer cell lines.

## Week 1

- **Correlations from file:** 
    ```bash
    python3 correlations.py DROSHA
    ```

- **Filtering per p-value (Threshold: -t):**
    ```bash
    python3 pvaluefiltering.py -G DROSHA -t 0.01 
    ```