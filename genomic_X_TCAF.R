load("genomic_X_TCAF.RData")
ls()

metadata <- read.csv("metadata.csv", header = T, stringsAsFactors=F)

####################################################################
# Check data from metadata and genomic data

str(metadata)

head(metadata)

names(metadata)

str(genotype)

str(sample_id)

str(POS)

mode(metadata)

mode(genotype)

mode(sample_id)

mode(POS)

####################################################################
# Take the genotype from samples of coluzzi Burkina Faso 2012

sample_id_col_2012<-metadata[metadata$country=="Burkina Faso" & metadata$year==2012 & metadata$aim_species=='coluzzii',]$sample_id

length(sample_id_col_2012)

sample_id %in% sample_id_col_2012

genotype_col_2012<-genotype[sample_id %in% sample_id_col_2012,]

dim(genotype_col_2012)

####################################################################
# Take the genotype from samples of coluzzi Burkina Faso 2014

sample_id_col_2014<-metadata[metadata$country=="Burkina Faso" & metadata$year==2014 & metadata$aim_species=='coluzzii',]$sample_id

length(sample_id_col_2014)

sample_id %in% sample_id_col_2014

genotype_col_2014<-genotype[sample_id %in% sample_id_col_2014,]

dim(genotype_col_2014)

####################################################################
# Take the genotype from samples of gambiae Burkina Faso 2012

sample_id_gamb_2012<-metadata[metadata$country=="Burkina Faso" & metadata$year==2012 & metadata$aim_species=='gambiae',]$sample_id

length(sample_id_gamb_2012)

sample_id %in% sample_id_gamb_2012

genotype_gamb_2012<-genotype[sample_id %in% sample_id_gamb_2012,]

dim(genotype_gamb_2012)

####################################################################
# Take the genotype from samples of gambiae Burkina Faso 2014

sample_id_gamb_2014<-metadata[metadata$country=="Burkina Faso" & metadata$year==2014 & metadata$aim_species=='gambiae',]$sample_id

length(sample_id_gamb_2014)

sample_id %in% sample_id_gamb_2014

genotype_gamb_2014<-genotype[sample_id %in% sample_id_gamb_2014,]

dim(genotype_gamb_2014)

####################################################################
# Calculate number of individuals of each sample, number of loci and convert
# genotype data to numeric data

fsum <- function(x) {.Primitive('sum')(x)}

s_col_2012 <- nrow(genotype_col_2012)

s_gamb_2012 <- nrow(genotype_gamb_2012)

s_col_2014 <- nrow(genotype_col_2014)

s_gamb_2014 <- nrow(genotype_gamb_2014)

K <- ncol(genotype_col_2012)

c(s_col_2012,s_col_2014, s_gamb_2012, s_gamb_2014, K)

mode(genotype_col_2012) <- 'integer'
mode(genotype_col_2014) <- 'integer'
mode(genotype_gamb_2012) <- 'integer'
mode(genotype_gamb_2014) <- 'integer'

####################################################################
# Calculate allele frequencies of each sample using the genotype,
# the number of individuals and the sum function

af_col_2012 <- apply(genotype_col_2012, 2, fsum)/(2*s_col_2012)
af_col_2014 <- apply(genotype_col_2014, 2, fsum)/(2*s_col_2014)
af_gamb_2012 <- apply(genotype_gamb_2012, 2, fsum)/(2*s_gamb_2012)
af_gamb_2014 <- apply(genotype_gamb_2014, 2, fsum)/(2*s_gamb_2014)

####################################################################
# Filter the allele frequencies to obtain only those whose MAF > 0.05

af_col_2012_TF <- af_col_2012 > 0.05 & af_col_2012 < 0.95
genotype_col_2012_MAF<-genotype_col_2012[,af_col_2012_TF]

af_col_2014_TF <- af_col_2014 > 0.05 & af_col_2014 < 0.95
genotype_col_2014_MAF<-genotype_col_2014[,af_col_2014_TF]

af_gamb_2012_TF <- af_gamb_2012 > 0.05 & af_gamb_2012 < 0.95
genotype_gamb_2012_MAF<-genotype_gamb_2012[,af_gamb_2012_TF]

af_gamb_2014_TF <- af_gamb_2014 > 0.05 & af_gamb_2014 < 0.95
genotype_gamb_2014_MAF<-genotype_gamb_2014[,af_gamb_2014_TF]

####################################################################
# Filter the genotypes which allele frequencies have a MAF > 0.05

genotype_col_2012_MAF_c<-genotype_col_2012[,af_col_2012_TF & af_col_2014_TF]
genotype_col_2014_MAF_c<-genotype_col_2014[,af_col_2012_TF & af_col_2014_TF]
genotype_gamb_2012_MAF_c<-genotype_gamb_2012[,af_gamb_2012_TF & af_gamb_2014_TF]
genotype_gamb_2014_MAF_c<-genotype_gamb_2014[,af_gamb_2012_TF & af_gamb_2014_TF]

####################################################################
# Take the positions of the genotypes which allele frequencies have a MAF > 0.05

POS_col <- POS[af_col_2012_TF & af_col_2014_TF]
POS_gamb <- POS[af_gamb_2012_TF & af_gamb_2014_TF]

####################################################################
# Remove data to store space

rm(genotype)
rm(genotype_col_2012, genotype_col_2014, genotype_gamb_2012, genotype_gamb_2014)
rm(genotype_col_2012_MAF, genotype_col_2014_MAF, genotype_gamb_2012_MAF, genotype_gamb_2014_MAF)
rm(af_col_2012, af_col_2014, af_gamb_2012, af_gamb_2014)
rm(af_col_2012_TF, af_col_2014_TF, af_gamb_2012_TF, af_gamb_2014_TF)
rm(POS)
rm(sample_id, sample_id_col_2012, sample_id_col_2014, sample_id_gamb_2012, sample_id_gamb_2014)
rm(metadata)
rm(s_col_2012, s_col_2014, s_gamb_2012, s_gamb_2014, K)

####################################################################
# Create a window sequence to divide the chromosome positions in intervals.
# Then, we will take only one position randomly for each interval to avoid linkage
# disequilibrium between SNPs positions that are really closer. This is done
# using only specie coluzzii.

# Notice that we are making windows until the position 12,5 Mb.
# Now, we are only taking the positions located in the higher recombination regions.
# In the chr X, the higher recombination regions are located between the
# position 1 and position 12,5 Mb. In the chr 3R are located between postion
# 1 and 37 Mb. Then, select the genotypes using the positions filter.
# This is done for coluzzii specie.

window_size <- 1000
window_start_col <- seq(from = 1, to = 12500000, by = window_size)
window_stop_col <- window_start_col + window_size -1

POS_col_link <- rep(NA, length(window_start_col))

for (x in 1:length(window_start_col)) {
  POS_col_TF <- POS_col < window_stop_col[x] & POS_col > window_start_col[x]
  POS_col_f <- POS_col[POS_col_TF]
  if(length(POS_col_f) == 0){
    
  }else if(length(POS_col_f) == 1){
    POS_col_link[x] = POS_col_f
  }else if(length(POS_col_f) > 1){
    POS_col_def <- sample(POS_col_f,size=1)
    POS_col_link[x] = POS_col_def
  }
}

POS_col_link <- na.omit(POS_col_link)

####################################################################
# Take only the genotypes of the coluzzii specie whose positions
# have been selected before

POS_col %in% POS_col_link
genotype_col_2012_link<-genotype_col_2012_MAF_c[,POS_col %in% POS_col_link]
genotype_col_2014_link<-genotype_col_2014_MAF_c[,POS_col %in% POS_col_link]

####################################################################
# Same process but with the specie gambiae

window_start_gamb <- seq(from = 1, to = 12500000, by = window_size)
window_stop_gamb <- window_start_gamb + window_size -1

POS_gamb_link <- rep(NA, length(window_start_gamb))

for (x in 1:length(window_start_gamb)) {
  POS_gamb_TF <- POS_gamb < window_stop_gamb[x] & POS_gamb > window_start_gamb[x]
  POS_gamb_f <- POS_gamb[POS_gamb_TF]
  if(length(POS_gamb_f) == 1){
    POS_gamb_link[x] = POS_gamb_f
  }else if(length(POS_gamb_f) > 1){
    POS_gamb_def <- sample(POS_gamb_f,size=1)
    POS_gamb_link[x] = POS_gamb_def
  }
}

POS_gamb_link <- na.omit(POS_gamb_link)

####################################################################
# Take only the genotypes of the gambiae specie whose positions
# have been selected before

POS_gamb %in% POS_gamb_link
genotype_gamb_2012_link<-genotype_gamb_2012_MAF_c[,POS_gamb %in% POS_gamb_link]
genotype_gamb_2014_link<-genotype_gamb_2014_MAF_c[,POS_gamb %in% POS_gamb_link]

####################################################################
# Give data appropiate names

genotype_col_2012 <- genotype_col_2012_link
genotype_col_2014 <- genotype_col_2014_link
genotype_gamb_2012 <- genotype_gamb_2012_link
genotype_gamb_2014 <- genotype_gamb_2014_link
rm(POS_col, POS_gamb)
POS_col <- POS_col_link
POS_gamb <- POS_gamb_link

####################################################################
# Remove data to store space

rm(genotype_col_2012_MAF_c, genotype_col_2014_MAF_c, genotype_gamb_2012_MAF_c, genotype_gamb_2014_MAF_c)
rm(POS_col_def, POS_col_f, POS_col_TF, POS_gamb_def, POS_gamb_f, POS_gamb_TF)
rm(window_size, window_start_col, window_stop_col, window_start_gamb, window_stop_gamb)
rm(genotype_col_2012_link, genotype_col_2014_link, genotype_gamb_2012_link, genotype_gamb_2014_link)
rm(POS_col_link, POS_gamb_link)

####################################################################
# Export data

save(genotype_col_2012,genotype_col_2014,genotype_gamb_2012,genotype_gamb_2014, POS_col, POS_gamb, file = "genomic_X_TCAF_filter.RData")