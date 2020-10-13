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

library(dplyr)

# Checking total read depth
Depth = data.frame(sample_sums(ps))
Depth$Sample.ID = row.names(Depth)
colnames(Depth)[1] = "SampleSums"
Depth = Depth %>%
    arrange(SampleSums)
head(Depth)

# It looked to me like this has already been done (no Compost or Muck)
# We don't want the Spooner samples removed anyway, I think
#ps = subset_samples(ps, Sample.ID != "049-W3-Compost" & Sample.ID!="004-K2-Muck" & Sample.ID!="028-S-0-30" & Sample.ID!="029-S-30-60")

# Melt
ps.melt = psmelt(ps)
head(ps.melt)

hist(sample_sums(ps))

# Relative abundances
ps.norm = transform_sample_counts(ps,function(x) x/sum(x))
head(otu_table(ps.norm))
sample_sums(ps.norm)
# Looks good

# Spooner Vs. Wisconsin
ps.wisc <- subset_samples(ps, Study=="Wisconsin")
ps.spoon <- subset_samples(ps, Study=="Spooner")

ps.norm.wisc <- subset_samples(ps.norm, Study=="Wisconsin")
ps.norm.spoon <- subset_samples(ps.norm, Study=="Spooner")

sample_names(ps.wisc)

# [...]

# Reminding myself what the data columns are called
colnames(sample_data(ps.norm.wisc))

ORD = ordinate(ps.norm.wisc,method="PCoA",distance="bray")
plot_ordination(ps.norm.wisc,ORD,color="Lab.CO2.pH.one2one")

# When I looked at these plots, I was thinking that the 
# 1:1 wasn't so different from the 1:2,3, or 4. So, we want
# to test that with PERMANOVA

### Running permanova using the vegan "adonis" function
# See here for a nice overview
# https://sites.google.com/site/mb3gustame/hypothesis-tests/manova/npmanova

# Set the phyloseq object, starting with Spooner
ps = ps.norm.spoon

# Collect the sample data
df = as(sample_data(ps), "data.frame")

# Calculate Bray-Curtis dissimilarity
d = distance(ps, method = "bray")

# Test whether the dissimilarity between samples (d)
# is significantly correlated with one given variable
s.adonis.Lp1 = adonis(d ~ Lab.CO2.pH.one2one, data=df)
s.adonis.Lp2 = adonis(d ~ Lab.CO2.pH.one2two, data=df)
s.adonis.Lp3 = adonis(d ~ Lab.CO2.pH.one2three, data=df)
s.adonis.Lp4 = adonis(d ~ Lab.CO2.pH.one2four, data=df)

# Looking at the results
s.adonis.Lp1
s.adonis.Lp2
s.adonis.Lp3
s.adonis.Lp4

# For the Spooner dataset at lab CO2, the 1:1 value has the lowest R2.
# They aren't huge changes, but I think we could say they are consistent
# with the possibility that the in situ measurements better represent
# microbial conditions.

# There are many possible variations on the model - in particular,
# including other variables, and then seeing what pH explains
# after controlling for those. E.g., for Wisc:

# Set the phyloseq object
ps = ps.norm.wisc
# Collect the sample data
df = as(sample_data(ps), "data.frame")
# Calculate Bray-Curtis dissimilarity
d = distance(ps, method = "bray")

# Test whether the dissimilarity between samples (d)
# is significantly correlated with one given variable
w.adonis.Lp1.C = adonis(d ~ Perc.Clay+Lab.CO2.pH.one2one, data=df)
w.adonis.Lp2.C = adonis(d ~ Perc.Clay+Lab.CO2.pH.one2two, data=df)
w.adonis.Lp3.C = adonis(d ~ Perc.Clay+Lab.CO2.pH.one2three, data=df)
w.adonis.Lp4.C = adonis(d ~ Perc.Clay+Lab.CO2.pH.one2four, data=df)

# Looking at the results,
# We can see percent clay explains ~5% of the variation in the community data
# and after controlling for that, pH explains ~8-9%, without huge differences
# between methods, although 1:1 is marginally the worst predictor

w.adonis.Lp1.C
w.adonis.Lp2.C
w.adonis.Lp3.C
w.adonis.Lp4.C



###### Estimating Richness #########
# Install breakaway from github, not CRAN:
# install.packages("devtools")
# devtools::install_github("adw96/breakaway")
library(breakaway)

# breakaway also depends on the purrr package
library(purrr)

# Set phyloseq object
# NOTE: Richness estimates require non-normalized data!
ps = ps.spoon

# Drop the OTUs that have zeros across all remaining samples
# For some reason, it needs a taxonomy table to do this

# Making a fake taxonomy table
FakeTT = data.frame(matrix(nrow = dim(otu_table(ps))[1],ncol=2))
FakeTT = tax_table(FakeTT)
row.names(FakeTT)=taxa_names(ps)
colnames(FakeTT)=c("Taxonomy1","Taxonomy2")
head(FakeTT)

ps = phyloseq(sample_data(ps),otu_table(ps),FakeTT)
ps
# Drop the OTUs that have zeros across all remaining samples
ps.nozero = subset_taxa(ps,taxa_sums(ps)>0)
ps.nozero

###### Estimating Richness ######
##### Going through current breakaway tutorial ####
# https://adw96.github.io/breakaway/articles/breakaway.html

# Set OTU table
otu_data = otu_table(ps.nozero)
# Set sample data
SamDat = data.frame(sample_data(ps.nozero))
# Checking rownames match sample data
head(colnames(otu_data) == rownames(meta_data))

# Run Breakaway's frequency table list function
frequencytablelist = build_frequency_count_tables(otu_data)

# Check out one of them (#63)
head(frequencytablelist[[12]])

# Try Breakaway on a couple of samples
breakaway(frequencytablelist[[1]])
breakaway(frequencytablelist[[4]])

# Note: the error estimates are likely too low - without singletons
# breakaway runs a Poisson model, as implemented in CatchAll
# 

# Run the richness estimator (breakaway) on all our samples (lists of frequency tables)
RichEsts = lapply(frequencytablelist,breakaway)

# Pull out the estimates, errors, and the model
Estimate = as.matrix(map(RichEsts, "estimate"))
Error = as.matrix(map(RichEsts, "error"))
Model = as.matrix(map(RichEsts, "model"))
df = data.frame(Estimate,Error,Model)

# Add sample ID column, estimate, and error
df$Sample.ID = row.names(df)
df$Estimate=as.numeric(df$Estimate)
df$Error=as.numeric(df$Error)

SamDat$Sample.ID %in% df$Sample.ID

# Merge the estimates with the sample data
RichPlot = merge(SamDat,df,by="Sample.ID")
head(RichPlot)

# Plot the Richness estimates
p = ggplot(RichPlot,aes(y=Estimate,x=Lab.CO2.pH.one2one))
p = p + geom_point() + geom_errorbar(aes(ymin=Estimate-Error,ymax=Estimate+Error))
p