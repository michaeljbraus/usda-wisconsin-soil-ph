library(phyloseq)
library(dada2)
library(ggplot2)
library(RColorBrewer)
# Load libraries

# Import sequence processing table and check it out
# How many sequences were retained at each step?
track = readRDS("track.rds")
track

# Import final otu table
OTUs = readRDS("OTUtab.nochim.rds")
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
samdat = read.csv("Hatch-Metadata.csv")
samdat$Sample.ID = as.character(samdat$Sample.ID)
samdat$Sample.ID[1:9] = paste("00",samdat$Sample.ID[1:9],sep="")
samdat$Sample.ID[10:length(samdat$Sample.ID)] = paste("0",samdat$Sample.ID[10:length(samdat$Sample.ID)],sep="")

# Check we have all the same sample names now
samdat$Sample.ID == colnames(otutab)[1:65]
samdat = sample_data(samdat)

row.names(samdat) = samdat$Sample.ID

# Create phyloseq object
ps = phyloseq(otu_table=otutab,sample_data=samdat)

ps.melt = psmelt(ps)
ps.melt

head(ps.melt)

hist(sample_sums(ps))

# Relative abundances
ps.norm = transform_sample_counts(ps,function(x) x/sum(x))
head(otu_table(ps.norm))
sample_sums(ps.norm)
# Looks good

ps.norm.noC = prune_samples(sample_names(ps.norm) != "049-W3-Compost",ps.norm)
ps.norm.noC = prune_samples(sample_names(ps.norm.nocompost) != "004-K2-Muck",ps.norm.nocompost)
ps.norm.noC

colnames(sample_data(ps.norm))

ps.ord.nmds = ordinate(ps.norm.noC, "NMDS", "bray",trymax=1000,k=2)
p = plot_ordination(ps.norm.noC, ps.ord.nmds, type="samples")
palette = brewer.pal(8, "Spectral")
p = p + aes(colour=Soil.pH) + geom_point(size=3) + scale_colour_gradientn(colors=palette)
p

row.names(samdat)
ps.ord.pcoa = ordinate(ps, "PCoA", "bray")
p = plot_ordination(ps, ps.ord.pcoa, type="samples",label="Names")
p
