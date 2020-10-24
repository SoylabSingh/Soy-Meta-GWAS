##########Processing VCF.gz file after Beagle imputation (v5.0, 28Sep18.793, available from JMS)###########
memory.limit(1000000) #Due to large size of full panel SNP file, we need to increase RAM availability.  Virtual RAM is utilized in the form of a 1 TB pagefile on a Seagate External Hard Drive

setwd('/Users/User/Desktop') #Set working directory
GD <- read.table('fullimputed.vcf.gz') #read in the marker file.  Unfortunately, headers do not read in properly with the current system.
dim(GD)
GH <- read.delim('fullimputed.vcf.gz') #reads in the complete marker file as a ridicuously long df with one row.
GH[1:15,]
GF <- as.data.frame(GH[9:20104,]) #after skipping VCF headers, the cells that should be column names are subsetted from the long df
colnames(GD) <- GF[1:20096,] #column names are set 

convert.snp <- function(x) #a function is created to convert from VCF style marker encoding to 0,1,2 encoding, with 0 homozygous major and 1 heterozygous
{
  alleles <- sort(unique(x))
  y <- rep(NA, length(x))
  y[which(x == "0|0")] <- 0
  y[which(x == "0|1")] <- 1
  y[which(x == "1|0")] <- 1
  y[which(x == "1|1")] <- 2
  return(y)
}
GD$QUAL <- NULL
GD$INFO <- NULL
GD$FILTER <- NULL
GD$FORMAT <- NULL
geno <- apply(GD[1:nrow(GD),-(1:5)], 1, convert.snp) # this function is applied to the columns containing encoded markers

GZ <- t(geno) #the resultant df is transposed. I'm not sure why it inverted in the first place.
rownames(GZ) <- rownames(GD) #set rownames for the new df
colnames(GZ) <- colnames(GD)[6:ncol(GD)]#set colnames for the new df
GP <- cbind(GD[,1:5],GZ)#bind dfs 
rownames(GP) <- GP$ID#set rownames for new df
#write.csv(GP[,1:3],"left.csv")#write csv to manually remove the pound sign from the first column and convert chromosome number to numeric
GQ <- read.csv("left.csv")#read the edited csv back in
GQ$X <- NULL #remove rowname column that didn't read in as a rowname
rownames(GQ) <- GQ$ID #set rownames
GD <- cbind(GQ,GP[,4:20092]) #combine dfs
GM <- GD[,1:3] #create map file
myGD <- GD[,6:ncol(GD)] #create preliminary marker file
MyGD <- t(myGD) #transpose marker file
QQ <- as.data.frame(rownames(MyGD))
markers <- cbind(QQ,MyGD) #create new column for genotype names
names(markers)[1:2] <- c("PI","ss715578788")
GM <- GM[c(3,1,2)]
names(GM)[1:3] <- c("SNP","Chr","Position")
rm(geno,GQ,QQ,GP,GF,GH,GZ)
gc()
gc()
###########Install packages#######
install.packages("BiocManager")
BiocManager::install(c("Biobase","multtest"))
source("https://bioconductor.org/biocLite.R")
biocLite("Biobase", "multtest")
library(Biobase)
library(multtest)
install.packages('Rcpp', dependencies = TRUE)
library(Rcpp)
install.packages('gplots')
install.packages('LDheatmap')
install.packages('genetics')
#install.packages('compiler')
install.packages('scatterplot3d')
library(gplots)
library(LDheatmap)
library(genetics)
library(compiler)
library(scatterplot3d)
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
source("http://www.zzlab.net/GAPIT/emma.txt")
#########GWAS########
results_dir <- "/Users/User/Desktop/JMS/Mat Results 2018/"

markers$PI <- NULL
setwd("/Users/User/Desktop/Maturity Paper")
Y <- read.csv('BLUPs2018pt2.csv')
Y[is.na(Y)] <- ""
  
for(i in 2:ncol(Y)){
  y <- colnames(Y)[i]
  line = Y[names(Y) %in% c('gid', y)]
  line1 <- line[(line$gid %in% colnames(GD)),]#select lines in phenotype file which have been genotyped.
  myY1 <- line1[which(line1$gid %in% colnames(GD)),]#select lines in phenotype file which have been genotyped
  myY1[,2] <-as.numeric(myY1[,2])#turn the second column (phenotype) back into numeric type, as this was lost with line2
  myGD <- as.data.frame(t(GD[,which(colnames(GD) %in% myY1$gid)]))#limit GD file to lines which have observations available
  g= nrow(myY1)#set number of groups equal to number of lines for MLM
  dir.create(paste(results_dir,y))#create individual folders for each test
  setwd(paste(results_dir,y))#set working directory to newly created folder
  myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = GM, PCA.total =3, group.from =g, group.to = g, SNP.fraction = 1, SNP.MAF = 0.01)#Regular mixed linear model with no compression by setting group number = g (the number of PI).
 #myGAPIT <- GAPIT(Y = myY1, GD = myGD, GM = GM, PCA.total =10, Model.selection = TRUE)
}

#######Compile Results#######
setwd("/Users/User/Desktop/JMS/GWAS MSel Results/")
files <- list.files(pattern = "\\.GWAS.Results.csv$", recursive = TRUE)
data2=lapply(files, read.table, header=T, sep=",")
for (i in 1:length(data2)){data2[[i]]<-cbind(data2[[i]],files[i])}
data_rbind <- do.call("rbind", data2) 
colnames(data_rbind)[c(1:11)]<-c("SNP", "Chromosome", "Position", "P.value", "maf", "nobs", "Rsquare.of.Model.without.SNP", "Rsquare.of.Model.with.SNP", "FDR_Adjusted_P-values","effect","file.name")
signif <- data_rbind[(data_rbind[,4]<0.0001),]
write.csv(signif,'SigSNPs_2122019.csv')
