load("genomic_X_pi.RData")
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
# Take the genotype from samples of gambiae Burkina Faso 2012

sample_id_gamb_2012<-metadata[metadata$country=="Burkina Faso" & metadata$year==2012 & metadata$aim_species=='gambiae',]$sample_id

length(sample_id_gamb_2012)

sample_id %in% sample_id_gamb_2012

genotype_gamb_2012<-genotype[sample_id %in% sample_id_gamb_2012,]

dim(genotype_gamb_2012)

####################################################################
# Calculate number of individuals of each sample, number of loci and convert
# genotype data to numeric data

fsum <- function(x) {.Primitive('sum')(x)}

s_col_2012 <- nrow(genotype_col_2012)
s_gamb_2012 <- nrow(genotype_gamb_2012)
K <- ncol(genotype_col_2012)

c(s_col_2012, s_gamb_2012, K)

mode(genotype_col_2012) <- 'integer'
mode(genotype_gamb_2012) <- 'integer'

####################################################################
# Remove data to store space

rm(genotype)
rm(sample_id, sample_id_col_2012, sample_id_gamb_2012)
rm(metadata)

####################################################################
# Sliding window analysis in coluzzi

window_size <- 100000
window_stop <- seq(from = window_size, to = max(POS), by = window_size)
window_start <- window_stop - window_size
mid = window_start + (window_stop-window_start)/2

pi_vector <- rep(NA, length(window_start))

for (x in 1:length(window_start)){
  POS_TF <- POS < window_stop[x] & POS > window_start[x]
  gen <- genotype_col_2012[,POS_TF]
  if(ncol(gen)>0){
    af_p <- apply(gen, 2, fsum)/(2*s_col_2012)
    af_q <- 1- af_p
    pi <- fsum(2*af_p*af_q)/ncol(gen)
    pi_vector[x] = pi
  }}

plot(mid, pi_vector, xlab = "Position in bp", ylab = "Nucleotide diversity (pi)", type="l", ylim=c(0, 0.0125), main="Sliding window nucleotide diversity coluzzi Chr X")
# We can see that after position 12,5 Mb nucleotide diversity or pi
# decreases, so the higher recombination regions are located between the
# position 1 and position 12,5 Mb

####################################################################
# Sliding window analysis in gambiae

pi_vector <- rep(NA, length(window_start))

for (x in 1:length(window_start)){
  POS_TF <- POS < window_stop[x] & POS > window_start[x]
  gen <- genotype_gamb_2012[,POS_TF]
  if(ncol(gen)>0){
    af_p <- apply(gen, 2, fsum)/(2*s_gamb_2012)
    af_q <- 1- af_p
    pi <- fsum(2*af_p*af_q)/ncol(gen)
    pi_vector[x] = pi
  }}

plot(mid, pi_vector, xlab = "Position in bp", ylab = "Nucleotide diversity (pi)", type="l", ylim=c(0, 0.0125), main="Sliding window nucleotide diversity gambiae Chr X")
# Same observation as before, after position 12,5 Mb nucleotide diversity or pi
# decreases, so the higher recombination regions are located between the
# position 1 and position 12,5 Mb  

####################################################################
# Create a window sequence to divide the chromosome positions in intervals.
# Then, we will take only one position randomly for each interval to avoid linkage
# disequilibrium between SNPs positions that are really closer. This is done
# using both species because they share the same positions. Then, select the
# genotypes using the positions filter.

window_size <- 1000
window_start <- seq(from = 1, to = 12500000, by = window_size)
window_stop <- window_start + window_size -1

POS_link <- rep(NA, length(window_start))

for (x in 1:length(window_start)) {
  POS_TF <- POS < window_stop[x] & POS > window_start[x]
  POS_f <- POS[POS_TF]
  if(length(POS_f) == 0){
    
  }else if(length(POS_f) == 1){
    POS_link[x] = POS_f
  }else if(length(POS_f) > 1){
    POS_def <- sample(POS_f,size=1)
    POS_link[x] = POS_def
  }
}

POS_link <- na.omit(POS_link)

####################################################################
# Take only the genotypes of both species whose positions
# have been selected before

POS %in% POS_link
genotype_col_2012_link<-genotype_col_2012[,POS %in% POS_link]
genotype_gamb_2012_link<-genotype_gamb_2012[,POS %in% POS_link]

####################################################################
# Remove data to store space

rm(genotype_col_2012, genotype_gamb_2012, gen)
rm(POS_def, POS_f, POS_TF, P, POS)
rm(window_size, window_start, window_stop, mid)
rm(af_p, af_q, pi, pi_vector)
rm(s_col_2012, s_gamb_2012, K)

####################################################################
# Calculate nucleotide diversity or pi using all the alleles from the
# higher recombination rates

s_col_2012 <- nrow(genotype_col_2012_link)
s_gamb_2012 <- nrow(genotype_gamb_2012_link)
K <- ncol(genotype_col_2012_link)

c(s_col_2012, s_gamb_2012, K)

af_col_2012_1 <- apply(genotype_col_2012_link, 2, fsum)/(2*s_col_2012)
af_col_2012_TF <- af_col_2012_1 > 0 & af_col_2012_1 < 1
af_col_2012_p <- af_col_2012_1[af_col_2012_TF]
af_col_2012_q <- 1 - af_col_2012_p

pi_col <- fsum(2*af_col_2012_p*af_col_2012_q)/K
mu <- 2.8e-9
N_col <- pi_col/(4*mu)
N_col

af_gamb_2012_1 <- apply(genotype_gamb_2012_link, 2, fsum)/(2*s_gamb_2012)
af_gamb_2012_TF <- af_gamb_2012_1 > 0 & af_gamb_2012_1 < 1
af_gamb_2012_p <- af_gamb_2012_1[af_gamb_2012_TF]
af_gamb_2012_q <- 1 - af_gamb_2012_p

pi_gamb <- fsum(2*af_gamb_2012_p*af_gamb_2012_q)/K
mu <- 2.8e-9
N_gamb <- pi_gamb/(4*mu)
N_gamb

####################################################################
# Remove data to store space

rm(genotype_col_2012_link, genotype_gamb_2012_link)
rm(af_col_2012_1, af_col_2012_p, af_col_2012_q, af_col_2012_TF)
rm(af_gamb_2012_1, af_gamb_2012_p, af_gamb_2012_q, af_gamb_2012_TF)
rm(K, mu, pi_col, pi_gamb, s_col_2012, s_gamb_2012)
rm(POS_link)
