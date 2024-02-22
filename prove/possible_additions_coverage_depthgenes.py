import pandas as pd
import numpy as np

def flatten_exons(gene, cov_map):
    contig = gene.contig.values[0]
    
    start = min(gene.start.values) - 1
    end = max(gene.end.values)
    size = end - start
    
    mask = np.full(size, False)
    
    # Create an array of indices for the entire range and then flatten it
    indices = np.array([np.arange(start, end) for start, end in zip(gene.start.values - 1, gene.end.values)]).flatten()
    
    # Use these indices to set the corresponding mask values to True
    mask[indices - start] = True
    
    sub_cov_map = cov_map[contig][start:start + size]
    
    try:
        return sub_cov_map[mask].mean()
    except Exception as e:
        return np.nan

def calculate_coverage(df, cov_map):
    return df.groupby("ID_gene").apply(flatten_exons, cov_map=cov_map)

# Load dataframes and coverage map
BUSCO_exons = pd.read_csv(snakemake.input["busco_exons"], sep="\t", names=["contig", "start", "end", "ID_gene"])
tot_genes = pd.read_csv(snakemake.input["genes_id"], sep="\t", names=["contig", "start", "end", "ID_gene"])
cov = pd.read_csv(snakemake.input["coverage"], sep="\t", names=["contig", "position", "coverage"])
cov_map = cov.groupby("contig")["coverage"].apply(np.array).to_dict()

# Calculate coverage
coverage_genes_busco = calculate_coverage(BUSCO_exons, cov_map)
coverage_tot_genes = calculate_coverage(tot_genes, cov_map)

# Rest of the code remains the same
THRSH = coverage_genes_busco.median()/8	#setting of a threshold value for the coverage to consider a gene "absent"

print(f"threshold: {THRSH}")

with open(snakemake.output["results_lista"], "w") as of:	# open the output file for writing
                                                            # the id each gene under threshold 
    for gene in coverage_tot_genes[coverage_tot_genes<THRSH].index:
        of.write("{}\n".format(gene))

coverage_tot_genes.to_csv(snakemake.output["results_tab"])