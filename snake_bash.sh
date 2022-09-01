#!/bin/bash

# Short bash code to run the second part of the pipeline recursively for each individual.
# CREATED ONLY FOR THE SECOND PART OF THE PIPELINE! THE FIRST PART DOES NOT NEED TO RUN RECURSIVELY
# DO NOT USE THIS IF YOU WANT TO RUN THE ANALYSIS IN PARALLEL

# change the paths inside the code to suit your own directories
cat /path/to/list/of/reads.txt | while read lines; do
  
  # the temporary file will store the individuals' names one at a time and run the
  # pipeline only on that individual, then change the individual name inside the file
  touch /path/to/a/temporary/file.txt
  
  echo "$lines" >> /path/to/a/temporary/file.txt
   
   cat /path/to/a/temporary/file.txt | while read line; do
    
    # change/add the flag you need to the snakemake command, change
    # the core numbers to suit your machine architecture.
    time snakemake --cores 16 -F
    
    # this passage is not mandatory, but useful to check that all the
    # individuals have been analyzed and the pipeline and this script worked properly
    echo "$line" >> /path/to/list/of/analyzed/individuals.txt
  done
  
  # removing temporary file to create another, clean one
  # with the next individual (we don't want to analyze the
  # same individual two or more times)
  rm /path/to/a/temporary/file.txt
done
