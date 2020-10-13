#!/bin/bash

# Get the fastq sequences tar from Gluster
# This should be a folder with all fwd and rev reads in it (and nothing else)
# The whole thing is zipped as a .tar.gz
cp /mnt/gluster/twhitman/Braus_pH/Braus_pH.tar ./

# unzip the .fastq sequencing files
tar -xzvf Braus_pH.tar
ls
ls Braus_pH
# Should yield a folder  with all the fwd and rev reads in it
# in the form of .fastq.gz

# Get rid of the .tar.gz file so it's not transferred back automatically 
rm *.tar

# Get the R installation tar from Gluster
cp /mnt/gluster/twhitman/R/R.tar.gz ./

# untar R installation and remove .tar.gz file
tar -xzf R.tar.gz
rm R.tar.gz

# Make sure script will use R installation
export PATH=$(pwd)/R/bin:$PATH
export RHOME=$(pwd)/R

#run R script for dada2 analysis
Rscript dada2.R

# Zip up the relevant files
mkdir Dada2_Results
mv OTUtab.rds Dada2_Results/
mv OTUtab.nochim.rds Dada2_Results/
mv DADA2_seqs.fasta Dada2_Results/
mv DADA2_seqs_nochim.fasta Dada2_Results/
mv track.rds Dada2_Results/
mv pErrR.rds Dada2_Results/
mv pErrF.rds Dada2_Results/

tar -czvf Dada2_Results.tar.gz Dada2_Results/

# Rename it to something specific
mv Dada2_Results.tar.gz Dada2_Results_Braus_pH.tar.gz

# Move it out to Gluster
cp Dada2_Results_Braus_pH.tar.gz /mnt/gluster/twhitman/Braus_pH

# Remove remaining files
rm -r Braus_pH
rm -r Dada2_Results
rm R.gz
rm *.fastq
