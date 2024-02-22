import pandas as pd
import numpy as np
import gc  # Garbage Collector interface

def flatten_exons(gene, cov_map):
    contig = gene.contig.values[0]

    start = min(gene.start.values) - 1
    end = max(gene.end.values)
    size = end - start

    mask = np.full(size, False, dtype=bool)

    # Vectorized operation to create the mask

    starts = gene.start.values - start
    ends = gene.end.values - start

    mask[starts[None, :] <= np.arange(size)[:, None]] = True
    mask[np.arange(size)[:, None] < ends[None, :]] = True

    sub_cov_map = cov_map[contig][start:end]

    # In-place operation to apply the mask
    
    np.putmask(sub_cov_map, ~mask, np.nan)

    return np.nanmean(sub_cov_map)

def calculate_coverage(df, cov_map):
    # This function will replace the repeated groupby and apply pattern
    return df.groupby("ID_gene").apply(flatten_exons, cov_map=cov_map)

def main():
    # Read input files
    BUSCO_exons = pd.read_csv(snakemake.input["busco_exons"], sep="\t", names=["contig", "start", "end", "ID_gene"])
    tot_genes = pd.read_csv(snakemake.input["genes_id"], sep="\t", names=["contig", "start", "end", "ID_gene"])
    cov = pd.read_csv(snakemake.input["coverage"], sep="\t", names=["contig", "position", "coverage"])

    # Create coverage map
    cov_map = cov.groupby("contig")["coverage"].apply(np.array).to_dict()

    # Calculate coverage
    coverage_genes_busco = calculate_coverage(BUSCO_exons, cov_map)
    coverage_tot_genes = calculate_coverage(tot_genes, cov_map)

    THRSH = coverage_genes_busco.median()/8	#setting of a threshold value for the coverage to consider a gene "absent"

    print(f"threshold: {THRSH}")

    with open(snakemake.output["results_lista"], "w") as of:	# open the output file for writing
                                                                # the id each gene under threshold 
        for gene in coverage_tot_genes[coverage_tot_genes<THRSH].index:
            of.write("{}\n".format(gene))

    # Once cov_map is no longer needed, explicitly delete it
    del cov_map
    gc.collect()  # Explicitly call the garbage collector

    # Save results
    coverage_tot_genes.to_csv(snakemake.output["results_tab"])

if __name__ == "__main__":
    main()