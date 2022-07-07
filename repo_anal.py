#!/usr/bin/env python
# coding: utf-8

from collections import defaultdict
import plotly.express as px
import os
import pandas as pd
import numpy as np
from mdutils.mdutils import MdUtils as mdu
from mdutils import Html


reports = {}


dirs = [snakemake.params["dir1"], snakemake.params["dir2"]]		#list with the directories used (fastqc 1 / 2)
for one_dir in dirs:		#takes info from fastqc.txt files and extracts them into a reports vairable
    in_file = f'/mnt/HADES/Fogal/Crassostrea_gigas/genomic/again/qc_reports/{one_dir}/fastqc_data.txt'
    report = defaultdict(dict)
    name = in_file.split('/')[-2]
    rows = []
    with open(in_file) as handle:
        for line in handle.readlines():
            if line[:2] == '##':
                continue
            if line[:2] == ">>":
                if line.rstrip() != '>>END_MODULE':
                    module_name, status = line.rstrip()[2:].split("\t")
                    report[module_name]['status'] = status
                else:
                    report[module_name]['data'] = pd.DataFrame(rows, columns=header)
                    report[module_name]['data']['sample'] = one_dir
            else:
                if line[0] == '#':
                    header = line.rstrip().split("\t")
                    rows = []
                    report[module_name]['data'] = pd.DataFrame(columns=line.rstrip().split())
                else:
                    rows.append(line.rstrip().split("\t"))
                    #print(name, module_name,line.rstrip())
    reports[one_dir] = report

dc_1 = reports[snakemake.params["dir1"]]		
dc_1 = dc_1["Per base sequence content"]
df_1 = pd.DataFrame.from_dict(dc_1["data"])


dc_2 = reports[snakemake.params["dir2"]]
dc_2 = dc_2["Per base sequence content"]
df_2 = pd.DataFrame.from_dict(dc_2["data"])


df_1["G"] = df_1["G"].astype(float)
df_2["G"] = df_2["G"].astype(float)


bigger_1 = df_1["G"].median()+0.5
bigger_2 = df_2["G"].median()+0.5
smaller_1 = df_1["G"].median()-0.5
smaller_2 = df_2["G"].median()-0.5


for i in range(0, len(df_1.G)):
    med_1 = df_1.G[i:i+7].median() #calculating median on selected window
    st_dev_1 = df_1.G[i:i+7].std()
    if (med_1 < bigger_1 or med_1 > smaller_1) and (st_dev_1 < 0.3): #condition: if the value of the median is bigger or smaller than a certain value
        trim_pos_1 = i
        print(trim_pos_1)
        break


for i in range(0, len(df_2.G)):
    med_2 = df_2.G[i:i+7].median() #calculating median on selected window
    st_dev_2 = df_2.G[i:i+7].std()
    if (med_2 < bigger_2 or med_2 > smaller_2) and (st_dev_2 < 0.3): #condition: if the value of the median is bigger or smaller than a certain value
        trim_pos_2 = i
        print(trim_pos_2)
        break


if trim_pos_1 > trim_pos_2 :
    trim_pos = trim_pos_1
else:
    trim_pos = trim_pos_2
 

# In[21]:

#fl = open("/mnt/HADES/Fogal/Crassostrea_gigas/snakemake_files/again/config.yaml", "w")
#fl.write("trim_pos : " + str(trim_pos))
#fl.close()

bigger_1 = str(bigger_1)
bigger_2 = str(bigger_2)
smaller_1 = str(smaller_1)
smaller_2 = str(smaller_2)
trim_pos = str(trim_pos)

report = mdu(file_name = snakemake.output["reporto"], title =  snakemake.wildcards["sample"]+" Reads Report")
report.new_header(level = 1, title = " ")
report.new_header(level = 2, title = "1. Upper limit chosen for the median composition based on appearance of base G:")
report.new_line("    ")
report.new_line("Read 1: ") 
report.write(bigger_1)
report.new_line("Read 2: ")
report.write(bigger_2)
report.new_line("    ")
report.new_header(level = 2, title = "2. Lower limit chosen for the median composition based on appearance of base G:")
report.new_line("    ")
report.new_line("Read 1: ")
report.write(smaller_1)
report.new_line("Read 2: ")
report.write(smaller_2)
report.new_line("    ")
report.new_header(level = 2, title = "3. Bases trimmed from the start of the reads: ")
report.new_line("    ")
report.new_line(trim_pos)
report.write(" bases.")
report.create_md_file()

