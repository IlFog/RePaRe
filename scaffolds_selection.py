from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pds
import matplotlib as mpl
import seaborn as sbn
import matplotlib.pyplot as plt
import gzip


mpl.use('Agg')

#CODE TO CALCULATE LIMITS FOR GC COUNT TO SELECT DATA

GC_dict = {} #creation of a void dictionary to store GC count

#opening of gzip file of assembly and GC count of the different scaffolds
with gzip.open(snakemake.input["gzp"], "rt") as handle:
    for seq_record in SeqIO.parse(handle, "fasta"):
        count = 100*float(seq_record.seq.count("G")+seq_record.seq.count("C")+seq_record.seq.count("g")+seq_record.seq.count("c"))/(len(seq_record.seq)-(seq_record.seq.count("N")+seq_record.seq.count("n")))
        GC_dict[seq_record.id] = count

#transformation of the dictionary to a pandas DataFrame
GC_items = GC_dict.items()

GC_list = list(GC_items)

GC_df = pds.DataFrame (GC_list)

#renaming of DataFrame columns
GC_df.columns = ["ID", "Count"]

#calculating limits for the GC count for later use
sup = GC_df.Count.quantile(0.975)
inf = GC_df.Count.quantile(0.025)

#CODE TO CALCULATE LIMITS FOR DEPTH COVERAGE SELECTION
#BUSCO GENES DEFINE HOMOZYGOUS PEAK

#opening and readin .csv file with depth coverage for genes
genes_op = open (snakemake.input["bcsv"])
genes_cover = pds.read_csv (genes_op)

#creating void list
buscor = []

#opening and reading file with busco genes names and inserting them in list
busco = open (snakemake.input["blist"])
buscor = [i.rstrip() for i in busco.readlines ()]

#selecting genes with a cover inferior to 200
genes_cover1 = genes_cover[genes_cover["0"]<200]

#selecting busco genes
busco_cover = genes_cover1[genes_cover1["ID_gene"].isin(buscor)]

#calculating median to locate the homozygous peak
busco_median = busco_cover.median()

#extracting median number from DataFrame
bmed = busco_median.iloc[0]

#calculating upper and lower limit based on homo- and heterozygous peaks
superior_busco = bmed * 1.35
inferior_busco = (bmed / 2)*0.65


#CODE TO SELECT SCAFFOLDS I AM INTERESTED IN

#creating void dictionary to store GC count values
GC_dictionary = {}

#opening fasta file containing scaffolds sequences and calculating GC count for them, storing values in dictionary
scaf_seq = SeqIO.parse (open(snakemake.input["scaf"]), 'fasta')
for seq_records in scaf_seq:
    count = 100*float(seq_records.seq.count("G")+seq_records.seq.count("C"))/(len(seq_records.seq)-seq_records.seq.count("N"))
    GC_dictionary[seq_records.id] = count

#opening coverage .csv file and storing it in a DataFrame
cov_csv = open(snakemake.input["cov"])
cov_read = pds.read_csv(cov_csv, sep="\t")

#creating a smaller DataFrame containing just the values I am interested in
sub_cov_df = cov_read[["#rname","endpos","meandepth"]]

#adding a new column to the DataFrame containing GC counts, pairing the values based on the scaffold ID
sub_cov_df["GC_count"] = sub_cov_df["#rname"].apply(lambda x : GC_dictionary[x])

#renaming the columns
sub_cov_df.columns= ["Name", "Length", "Depth", "Count"]

#selection of the scaffolds based on length
smaller_len = sub_cov_df[sub_cov_df.Length > 1000]

#selection of the scaffolds based on depth of coverage
inflim = smaller_len[smaller_len.Depth > inferior_busco]
smaller_len_cov = inflim[inflim.Depth < superior_busco]

#selection of the scaffolds based on GC count
inf_count = smaller_len_cov[smaller_len_cov.Count > inf]
final_scaffolds_list = inf_count[inf_count.Count < sup]

#creates file to store fasta sequences
fi = open(snakemake.output["selected"], "a")

#writes IDs from DataFrame to a list
ID_list = final_scaffolds_list["Name"].values.tolist()

#controls if sequences ID is in the ID list from the DataFrame and if so it writes the sequence into the new file
all_sc = open(snakemake.input["pgm"], "a")
scaf_seq = SeqIO.parse (open(snakemake.input["scaf"]), 'fasta')
for sequences in scaf_seq:
    if sequences.id in ID_list:
        r = SeqIO.write(sequences, fi, "fasta")
        p = SeqIO.write(sequences, all_sc, "fasta")

#creates plot from not-selected DataFrame
grafico = sbn.jointplot(sub_cov_df.Count, sub_cov_df.Depth, kind="kde")
#creates plot from selected DataFrame
grafico1 = sbn.jointplot(final_scaffolds_list.Count, final_scaffolds_list.Depth, kind="kde")

#saves both plots into different .svg files
grafico.savefig(snakemake.output["graf_nonsel"], format="svg")
grafico1.savefig(snakemake.output["graf_sel"], format="svg")
