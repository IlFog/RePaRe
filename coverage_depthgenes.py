import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# returns coverage of contiguous regions of the genome 
# that are covered by at least one annotation as "exon"
# all other regions are excluded.
# The gene parameter is from a coordinate dataframe
# The cov map is global and is derived from the mapping
def flatten_exons(gene):
    global cov_map
    contig = gene.contig.values[0]
    start = min(gene.start.values)-1
    end = max(gene.end.values)
    size = end - start
    mask = np.full([size], False)
    pairs= np.array([gene.start.values, gene.end.values])
    startendarray = pairs.flatten('F')
    for i in range(0, startendarray.size, 2):
        this_exon = np.arange(startendarray[i]-start,startendarray[i+1]-start)
        mask[this_exon] = True
    sub_cov_map = cov_map[contig][start:start+size]
    try:
        return sub_cov_map[mask].mean()
    except Exception as e:
        return np.nan



BUSCO_exons = pd.read_csv(snakemake.input["busco_exons"], sep="\t" , names= ["contig", "start", "end", "ID_gene"])	#creates dataframe of busco exons from tsv file

#print(snakemake.input["busco_exons"])	#prints the input file to check the correct execution of the script
#print(BUSCO_exons.head())

tot_genes = pd.read_csv(snakemake.input["genes_id"], sep = "\t", names = ["contig", "start", "end", "ID_gene"])		# creates dataframe of genes from tsv file
                                                                                                                    # the columns are: contig id, gene start, gene end, gene id 

#print(snakemake.input["genes_id"])	#prints to check the correct execution of the script
#print(tot_genes.head())

cov = pd.read_csv(snakemake.input["coverage"], sep = "\t", names = ["contig", "position", "coverage"])		#creates dataframe from mapping depth tsv file
base_name = snakemake.input["coverage"].split("/")[-1]	# obtain the base name of the coverage file

#print(snakemake.input["coverage"])	#prints to check the correct execution of the script
#print(cov.head())

cov_map = cov.groupby("contig")["coverage"].apply(np.array).to_dict()	#groups cov dataframe by contig, turning the coverage value into a dictionary of arrays (key = contig id) where each positions indicates the coverage of each contig at that postion

print(f"length of the coverage map: {len(cov_map)}")	#print to check the correct execution of the script, should be the number of contigs in the genome

coverage_genes_busco = BUSCO_exons.groupby("ID_gene").apply(flatten_exons)	#groups busco_exons dataframe by ID and applies the flatten_exons function

#print(coverage_genes_busco.head())	#print to check the correct execution of the script

THRSH = coverage_genes_busco.median()/8	#setting of a threshold value for the coverage to consider a gene "absent"

print(f"threshold: {THRSH}")		#print to check the correct execution of the script

coverage_tot_genes = tot_genes.groupby("ID_gene").apply(flatten_exons)		#groups tot_genes by ID and applies the previously defined function

#print(coverage_tot_genes.head())	#print to check the correct execution of the script


with open(snakemake.output["results_lista"], "w") as of:	# open the output file for writing
                                                            # the id each gene under threshold 
    for gene in coverage_tot_genes[coverage_tot_genes<THRSH].index:
        of.write("{}\n".format(gene))
del cov_map

coverage_tot_genes.to_csv(snakemake.output["results_tab"])	#turns coverage_tot_genes into a csv file

