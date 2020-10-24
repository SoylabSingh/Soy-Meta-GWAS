#This code was developed to transfer the format and imputate the missing datapoints 
#of the SNP data downloaded from the Soybase to form GD and GM files that work for 
#GAPIT. The only thing needs to do is to install the packages and modify the work 
#directory to tell R where to read the input file and where to save the created files.

#Install needed packages
install.packages("synbreed")
install.packages("gdata")

#Load needed libraries
library(synbreed)
library(gdata)

#Read in VCF from Soybase, shuffle, label, examine
setwd("/Users/jmshook/Desktop")
GD <- read.vcf2matrix("fullpanel.vcf.gz", cores=20)
dim(GD)
head(GD)[,1:5]
tail(GD)[,1:5]
GD <- GD[c(3,1,2,4:ncol(GD))]
names(GD)[1:3] <- c("SNP","Chromosome","Position")
dim(GD)

#Remove markers on scaffolds, shift columns, convert chromosome keys to numbers, examine
GD <- GD[which(nchar(as.character(GD$Chrom)) < 5),]
dim(GD)
table(GD$Chrom)
GD$chrom <- as.numeric(substr(GD$Chrom,3,4))
GD <- GD[c(1,ncol(GD),3:(ncol(GD)-1))]
GD <- GD[order(GD$chrom,GD$Position),]
head(GD)[,1:5]
tail(GD)[,1:5]

#Convert letter-based SNPs to numeric for analysis
convert.snp <- function(x) 
  {
  alleles <- sort(unique(x))
  y <- rep(NA, length(x))
  y[which(x == "A")] <- "A/A"
  y[which(x == "T")] <- "T/T"
  y[which(x == "C")] <- "C/C"
  y[which(x == "G")] <- "G/G"
  
  if ("H" %in% alleles)
{
  if (which(alleles=="H")==3) 
  {
    y[which(x == "H")] <- paste(alleles[1], alleles[2], sep = "/")
  }
  else
  {
    if ((which(alleles=="H")==2) & ("T" %in% alleles))
    {
      y[which(x == "H")] <- paste(alleles[1], alleles[3], sep = "/")
    }
    else
    {
      if ((which(alleles=="H")==2) & !("T" %in% alleles))
      {y[which(x == "H")] <- paste(alleles[1], "T", sep = "/")}
      else {y[which(x == "H")] <- "A/T"}  
     }
   }
 }
  return(y)
}
geno <- apply(GD[1:nrow(GD),-(1:3)], 1, convert.snp)

#Examine format of geno
dim(geno)
class(geno)

#Label geno rows and columns
rownames(geno) <- colnames(GD[4:ncol(GD)])
colnames(geno) <- GD[1:nrow(GD),1]

#View dataframes
GD[1:5,1:8]
geno[1:5,1:5]

#Create Mapfile of SNPs
Map <- GD[2:3]
class(Map)
colnames(Map) <- c("chr", "pos")
Map[1:4,1:2]
rownames(Map) <- colnames(geno)


##############filter and imputation Map and Geno###################################
#Perform Beagle imputation from the comfort of R
setwd("/Users/jmshook/Desktop/")
gData<- create.gpData(pheno=NULL, geno=geno, map=Map, covar= NULL, reorderMap=F,map.unit="bp")
gData_coded <- codeGeno(gData, maf= 0.01, nmiss= 0.1, impute=T, impute.type="beagle", label.heter="alleleCoding")

#Examine data and subset out the needed information to build a map
dim(gData$map)
dim(gData_coded$map)
gData_coded$map[1:5,1:2]
gData_coded$map$SNP <- colnames(gData_coded$geno)
dim(gData_coded$map)
filted_Map <- gData_coded$map[c(3,1,2)]
filted_Map[1:5,1:3]

#Create map file for use in GWAS
setwd("/Users/jmshook/Desktop/")
write.table(filted_Map, "all_GM.txt", row.names=F,quote=F, sep="\t")

#Examine data and subset out the needed information to compile all markers
dim(gData$geno)
dim(gData_coded$geno)
gData$geno[180:nrow(gData$geno),1:5]
gData_coded$geno[180:nrow(gData$geno),1:5]
str(gData_coded$geno)
class(gData_coded$geno)

#Create Marker file
geno_df <- as.data.frame(gData_coded$geno)
geno_df$PI <- rownames(geno_df)
cn <- ncol(geno_df)
geno_df[1:5,(cn-5):cn]
geno_df <- geno_df[c(cn,1:(cn-1))]
geno_df[1:5,1:5]
write.table(geno_df,row.names=F, quote=F, "all_GD.txt",sep="\t")
