#!/bin/bash

cat /mnt/HADES/Fogal/Crassostrea_gigas/snakemake_files/again/list_of_new_reads | while read lines; do
  touch /mnt/HADES/Fogal/Crassostrea_gigas/snakemake_files/again/lista_1
  echo "$lines" >> /mnt/HADES/Fogal/Crassostrea_gigas/snakemake_files/again/lista_1
   cat /mnt/HADES/Fogal/Crassostrea_gigas/snakemake_files/again/lista_1 | while read line; do
    time snakemake --cores 16 -F
    echo "$line" >> /mnt/HADES/Fogal/Crassostrea_gigas/snakemake_files/again/lista_done
  done
  rm /mnt/HADES/Fogal/Crassostrea_gigas/snakemake_files/again/lista_1
done
