#!/usr/bin/env python
# coding: utf-8


from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pds
import matplotlib as mpl
import seaborn as sbn
import matplotlib.pyplot as plt
import csv
from mdutils.mdutils import MdUtils as mdu
from mdutils import Html


#CODE TO CALCULATE LIIMITS FOR DEPTH COVERAGE SELECTION
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

#calculating upper and lower limit based on homo- and hemizygous peaks
superior_busco = bmed * 1.35
inferior_busco = (bmed / 2)*0.65
sup = str(round(superior_busco, 2))
inf = str (round(inferior_busco, 2))



# Will close handle cleanly
with open(snakemake.input["scaffolds"]) as fasta_file:  
    identifiers = []
    lengths = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        identifiers.append(seq_record.id)
        lengths.append(len(seq_record.seq))

#converting lists to pandas Series    
s1 = pds.Series(identifiers, name='ID')
s2 = pds.Series(lengths, name='Length')

#Gathering Series into a pandas DataFrame and rename index as ID column
scaffolds_df = pds.DataFrame(dict(ID=s1, Length=s2)).set_index(['ID'])

# Calculates number of total scaffolds built
num_scaffolds = len(scaffolds_df.index)
num_scaf = str(num_scaffolds)


# Will close handle cleanly
with open(snakemake.input["selec"]) as selec:  
    idf = []
    leng = []
    for seq_record in SeqIO.parse(selec, 'fasta'):  # (generator)
        idf.append(seq_record.id)
        leng.append(len(seq_record.seq))

#converting lists to pandas Series    
ser_1 = pds.Series(idf, name='ID')
ser_2 = pds.Series(leng, name='Length')

#Gathering Series into a pandas DataFrame and rename index as ID column
selected_df = pds.DataFrame(dict(ID=ser_1, Length=ser_2)).set_index(['ID'])

# Calculates number of selected scaffolds
num_selected = len(selected_df.index)
num_sel = str(num_selected)

# calculates number of Mb added
megabasi = str(sum(leng)/1000000)

# takes in reads count
with open(snakemake.input["calc_read"], "r") as file:
    read_counts = file.readline()
    reads_count = [float(i) for i in read_counts.strip().split()]
reads_counts = [str(reads_count) for reads_count in reads_count]




# Writes the report for each individual based on the data previously calculated and analyzed
# For the first individual it should be done manually
report = mdu(file_name = snakemake.output["report"], title =  snakemake.wildcards["sample"]+" Individual Report")
report.new_header(level = 1, title = " ")
report.new_header(level = 2, title = "1. Number of unmapped reads extracted from the mapping of the individual's reads on the reference genome:")
report.new_line("    ")
report.new_line(reads_counts[0])
report.write(" reads.")
report.new_paragraph()
report.new_header(level = 2, title = "2. Number of reads that mapped on the previous individuals' selected scaffolds:")
report.new_line("    ")
report.new_line(reads_counts[1])
report.write(" reads.")
report.new_paragraph()
report.new_header(level = 2, title = "3. Mapping rate for unmapped reads from genome onto individuals' selected scaffolds and number of reads that do not map on the scaffolds, respectively:")
report.new_line("    ")
report.new_line(reads_counts[2])
report.write(" mapping rate and ")
report.write(reads_counts[3])
report.write(" reads.")
report.new_paragraph()
report.new_header(level = 2, title = "4. Number of scaffolds assembled by SPAdes assembly program, with associated graph:")
report.new_line("    ")
report.new_line(num_scaf)
report.write(" scaffolds.")
report.new_line(report.new_inline_image(text = "", path = snakemake.input["svg_nonsel"]))
report.new_paragraph()
report.new_header(level = 2, title = "5. Parameters for scaffold selection:")
report.new_line("    ")
report.new_line("Scaffolds length: 1000 bp.\n")
report.new_line("    ")
report.new_line("Coverage: based on the coverage value found for Busco genes from the reference genome. Limits considered are 135% of homozygous peak and 65% of hemizygous peak.")
report.new_line("        ")
report.write("  - Upper limit = ")
report.write(sup)
report.new_line("        ")
report.write("  - Lower limit = ")
report.write(inf)
report.write("\n")

# WARNING!! 
# Modify this line with your second assembly, from which you calculate GC count. 
# Remember to also modify the upper and lower limit values.
report.new_line("GC count: based on GC count from " + report.new_inline_link(link = "https://www.ncbi.nlm.nih.gov/assembly/GCA_000297895.2", text = "ASM29789v2"))
report.write(" assembly from NCBI database. Considered 95th percentile as valid.")
report.new_line("        ")
report.write(" - Upper limit = 40,04") # CHANGE THIS
report.new_line("        ")
report.write(" - Lower limit = 28,89") # CHANGE THIS
report.new_paragraph()
report.new_header(level = 2, title = "6. Number of scaffolds after selection through parameters, with associated graph:")
report.new_line("    ")
report.new_line(num_sel)
report.write(" scaffolds. The number of megabases for the selected scaffolds are: ")
report.write(megabasi)
report.write(".")
report.new_line(report.new_inline_image(text = "", path = snakemake.input["svg_sel"]))
report.create_md_file()

