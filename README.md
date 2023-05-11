# RePaRe - Recursive Pangenome Reassembly and PAV analysis pipeline
<!--
da aggiungere istruzioni per:
- installazione conda environment
- modifica config file
- costruzione paths ed eventualmente directories
- utilizzo snakefile 1
- utilizzo snakefile 2
- utilizzo e modifica snake_bash
-->

RePaRe is a bioinformatic pipeline that aims to analyze resequencing data from as many individuals as possible coming from same species, in order to access the phenomenon of presence-absence variation in the species and conseuqently extract and collect the dispensable section of the genome.
The pipeline is composed of two sections: a parallel section, analyzing all the individuals at once, and an iterative section, analyzing and extracting data from one individual at a time.
These need to run separately at the moment.

**ATTENTION:** as of now the pipeline is tested on a bivalve organism, therefore it was built to analyze diployd species. While it would be possible to adapt it to analyze organisms with different rates of ploydy <!--esiste come parola? controllare-->, the code would need some changes and adaptations. Unless you yourself want to change the code and adapt it to work on different ploydy, **ONLY USE RePaRe TO WORK ON DIPLOYD ORGANISMS.**


## Installation

<!--
controllare che siano tutti i software necessari!
-->
I advise creating a conda environment with all the necessary softwares inside, in which the pipeline can run smoothly.
A list of links to the installation of the necessary softwares follows.

Snakemake installation:
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

Snakemake tutorial:
https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html

Bioawk installation:
https://anaconda.org/bioconda/bioawk

biopython installation:
https://biopython.org/wiki/Packages
https://anaconda.org/conda-forge/biopython

bwa installation:
https://anaconda.org/bioconda/bwa

fastp installation:
https://anaconda.org/bioconda/fastp

fastqc installation:
https://anaconda.org/bioconda/fastqc

mdutils installation:
https://pypi.org/project/mdutils/

numpy installation:
https://numpy.org/install/

pandas installation:
https://pandas.pydata.org/docs/getting_started/install.html

plotly installation:
https://plotly.com/python/getting-started/

samtools installation:
https://anaconda.org/bioconda/samtools

seaborn installation:
https://seaborn.pydata.org/installing.html

spades installation:
https://anaconda.org/bioconda/spades


## Usage of RePaRe

### Config file

The config file contains all the paths to the files needed by the two parts of the pipeline to work properly.
To have the pipeline working, make sure a list of the names of the samples is available. For example:

> SRRxxxxx_sample1
> 
> SRRxxxxx_sample2
> 
> SRRxxxxx_sample3

These will be used during the analysis as names for all the subsequent files and directories. If you want to change them, make sure to have the names of the reads matching the names in the list.

As can be seen in the config file, the needed paths and files are:

- list with name of the samples (first part of the pipeline)
- main working directory
- reference genome as FASTA file
- list with IDs of BUSCO exons
- list with IDs of exons of genes
- path to temporary list that will contain one sample at a time (second part of the pipeline)
- path to file in the reads directory (placeholder as output of the reads removal rule)
- file with names of BUSCO genes
- second assembly (to calculate GC content, can also be used the reference genome if a second assembly is not available)
- threads number for the mapping rules

If one or more of these files are not available, make sure to have them ready before starting up the pipeline.

### First section: parallel analysis

The first section of the pipeline runs in parallel. To start it, once everything has been installed and prepared, type in the terminal:
```bash
snakemake -s snakefile_optimized_part1
```

If you need to add other flags, follow the snakemake manual for instructions on which flags does what.
The pipeline should run smoothly on its own.

### Second section: iterative analysis

This section runs through a bash script that loops the analysis.

First, though, the first individual needs to run on its own. I'd advise to use this individual to make sure this section of the pipeline runs smoothly, rule by rule. 

**ATTENTION:** THE FOLLOWING RULES ARE NOT NEEDED FOR THE FIRST INDIVIDUAL, DO **NOT** RUN THEM WITH THE FIRST INDIVIDUAL

- rule bwa_index_1
- rule bwa_remapping_sort
- rule samtools_extraction_2
- rule samtools_conversion_2

**ATTENTION:** TAKE A LOOK AT report_script.py BEFORE RUNNING THE PIPELINE FOR THE FIRST TIME! THIS SCRIPT NEEDS A SLIGHT MODIFICATION FOR THE FIRST INDIVIDUAL (since some rules and analysis are not performed)

Be aware that consequently the input files for the `rule bwa_remapping_2` need to be changed to the unmapped reads extracted from the mapping on the reference genome.
As for the report script, a few data will not be available for the first individual, so make sure to change the script a bit before running it.
<!--sarebbe forse più semplice creare una pipeline adattata per il primo individuo? pensarci-->

Once the analysis of the first individual has been completed, the pipeline can run as it has been written.
To run it, use the command:

```bash
bash snake_bash.sh
```

If the config file was adequately compiled, everything should run smoothly wihtout the need of input from the operator.


Once the analysis is completed, an annotation analysis is advised, to discover how many genes have actually been collected with the dispensable part of the pangenome.
