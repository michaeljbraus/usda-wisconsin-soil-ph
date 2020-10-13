# 2020-03-12-hatch-final-data-merged

# Let's load some packages and the data set. 

library(ggplot2)
library(leaps)
library(cowplot)
dat <- read.csv(file = "metrology-compiled/2020-03-12-hatch-data-merged.csv", header = T)
str(dat)

## Biom + Metadata

# Intall:
# BiocManager::install("phyloseq")
# Load libraries
library(phyloseq)
#library(dada2)
library(ggplot2)
library(RColorBrewer)
library(vegan)

# Here we'll load in the microbial data. 

# Import sequence processing table and check it out
# How many sequences were retained at each step?
track = readRDS("microbial-compiled/track.rds")
track
# Import final otu table
OTUs = readRDS("microbial-compiled/OTUtab.nochim.rds")
# Change rows and columns so taxa are rows
OTUs = t(OTUs)
# Import as phyloseq object table
otutab = otu_table(OTUs, taxa_are_rows=TRUE)
head(otutab)
# Look at lowest-abundance samples
head(sort(sample_sums(otutab)))
name = row.names(otutab)[1]
name
hist(nchar(row.names(otutab)))
colnames(otutab)
# The sample names have suffixes appended to them - we just want sample IDs
# Get names
names = colnames(otutab)
# Remove last 17 characters
newnames = gsub('.{17}$', '', names)
# Check out names
newnames

# Assign new names to samples
colnames(otutab) = newnames

# Same thing if we wanted to simplify OTU names to IDs
# Instead of sequences (although that's useful)
otunames=row.names(otutab)
length(otunames)
newotunames = paste("OTU",rep(1:length(otunames)),sep="")
head(newotunames)
row.names(otutab) = newotunames
head(otutab)

# Cut off the sampling ID to match metadata
for (i in 1:length(colnames(otutab))){
  newname = strsplit(colnames(otutab[i]),"_")[[i]][1]
  colnames(otutab)[i]=newname
}

# Import sample data
samdat = read.csv("metrology-compiled/2020-03-12-hatch-data-merged-as-metadata.csv", header=T)
# samdat$Sample.ID = as.character(samdat$Sample.ID)
# samdat$Sample.ID[1:9] = paste("00",samdat$Sample.ID[1:9],sep="")
# samdat$Sample.ID[10:length(samdat$Sample.ID)] = paste("0",samdat$Sample.ID[10:length(samdat$Sample.ID)],sep="")

# Check we have all the same sample names now
samdat$Sample.ID[1:65] == colnames(otutab)[1:65]
samdat$Sample.ID[66:125] == colnames(otutab)[66:125]
samdat = sample_data(samdat)

# Exponentiate all pH values
samdat$DNA.Extr.Hplus.After.C1 <- 10^-samdat$DNA.Extr.pH.After.C1
samdat$DNA.Extr.Hplus.After.C2 <- 10^-samdat$DNA.Extr.pH.After.C2
samdat$Lab.CO2.Hplus.one2one <- 10^-samdat$Lab.CO2.pH.one2one
samdat$Lab.CO2.Hplus.one2two <- 10^-samdat$Lab.CO2.pH.one2two
samdat$Lab.CO2.Hplus.one2three <- 10^-samdat$Lab.CO2.pH.one2three
samdat$Lab.CO2.Hplus.one2four <- 10^-samdat$Lab.CO2.pH.one2four
samdat$High.CO2.Hplus.one2one <- 10^-samdat$High.CO2.pH.one2one
samdat$High.CO2.Hplus.one2two <- 10^-samdat$High.CO2.pH.one2two
samdat$High.CO2.Hplus.one2three <- 10^-samdat$High.CO2.pH.one2three
samdat$High.CO2.Hplus.one2four <- 10^-samdat$High.CO2.pH.one2four

row.names(samdat) = samdat$Sample.ID

# Create phyloseq object
ps = phyloseq(otu_table=otutab,sample_data=samdat)
ps = subset_samples(ps, Sample.ID != "049-W3-Compost" & Sample.ID!="004-K2-Muck" & Sample.ID!="028-S-0-30" & Sample.ID!="029-S-30-60")

# Melt
ps.melt = psmelt(ps)
ps.melt

head(ps.melt)

hist(sample_sums(ps))

# Relative abundances
ps.norm = transform_sample_counts(ps,function(x) x/sum(x))
head(otu_table(ps.norm))
sample_sums(ps.norm)
# Looks good

colnames(sample_data(ps.norm))

# Spooner Vs. Wisconsin

ps.wisc <- subset_samples(ps, Study=="Wisconsin")
ps.spoon <- subset_samples(ps, Study=="Spooner")

ps.norm.wisc <- subset_samples(ps.norm, Study=="Wisconsin")
ps.norm.spoon <- subset_samples(ps.norm, Study=="Spooner")

# [...]
