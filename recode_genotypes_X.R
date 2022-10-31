# LOAD THE PACKAGE AND OTHERS
#install.packages("hdf5r")
require(hdf5r)
require(compiler)

# MOUNT THE hdf5 FILE TO R, 'r' MEANS READ ONLY
dat<-H5File$new('X_BFAB_pass_biallelic_intergenic_complete.hdf5', 'r')
dat

dat[['pos']]

dat[['genotype']]

dat[['sample_id']][]
# R FUNCTIONS ARE ALLOWED AS THEY BECOME STANDARD R OBJECTS
length(dat[['sample_id']][])

# FOR LARGE DATASET YOU CAN LOAD A SMALL CHUNK OF IT
dat[['pos']][1:20]
# IT IS NOW A STANDARD R OBJECT, AND IT IS A INTEGER VECTOR, SO WE CAN CALCULATE THE MEAN, SAY
mean(dat[['pos']][1:20])

# OR EVEN PULL THE ENTIRE POS OUT
pos<-dat[['pos']][]
min(pos)
max(pos)

# THE FIRST SNP
dat[['genotype']][,,1]
dim(dat[['genotype']][,,1])

# THE FIRST INDIVIDUAL, SNPs 1:20
dat[['genotype']][,1,1:20]

names(dat[['genotype']])
dat[['genotype']]$dims

dyn.load('Rside.dll')

# CREATE A 2D GENOTYPE MATRIX
genotype<-matrix(0, nr=dat[['genotype']]$dims[2], nc=dat[['genotype']]$dims[3])
mode(genotype)<-'raw'
dim(genotype)

# START POINT AND END POINT
s1<-seq(1, ncol(genotype), 20000)
s2<-s1+20000-1
s2<-sapply(s2, min, ncol(genotype))

# FUNCTION TO CALCULATE UNIQUE ALLELE
f<-function(j)
{
  temp<-unique(c(g1[,j], g2[,j]))
  return(as.raw(sort(temp)))
}
f<-cmpfun(f)

for (i in 1:length(s1))
{
  # ONE CHUNK OF GENOTYPE
  g1<-dat[['genotype']][1,,s1[i]:s2[i]]
  g2<-dat[['genotype']][2,,s1[i]:s2[i]]
  # CREATE UNIQUE ALLELE LIST
  unique_list<-lapply(1:ncol(g1), f)
  mode(g1)<-'raw'
  mode(g2)<-'raw'
  temp<-.Call('recode', g1, g2, unique_list, 12)
  genotype[,s1[i]:s2[i]]<-temp
}
gc()

dyn.unload('Rside.dll')

# THEY ARE STILL 32BIT INTEGER HERE
dim(genotype)
object.size(genotype)

# CHECKING
sum(genotype==as.raw(0))+sum(genotype==as.raw(1))+sum(genotype==as.raw(2))
nrow(genotype)*ncol(genotype)

# SAVE EVERYTHING AS RData, HIGHLY COMPRESSED
POS<-dat[['pos']][]
sample_id<-dat[['sample_id']][]
save(genotype, POS, sample_id, file='genomic_X_TCAF.RData')

