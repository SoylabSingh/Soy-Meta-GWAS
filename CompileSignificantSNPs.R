setwd("D:/JMS/GWAS Results/MLM/Comp")
files <- list.files(pattern = "\\.GWAS.Results.csv$", recursive = TRUE)
data2=lapply(files, read.table, header=T, sep=",")
for (i in 1:length(data2)){
  markers = dim(data2[[i]])
  data2[[i]]<-cbind(data2[[i]],files[i],markers[1])
  }
data_rbind <- do.call("rbind", data2) 
colnames(data_rbind)[c(1:12)]<-c("SNP", "Chromosome", "Position", "P.value", "maf", "nobs", "Rsquare.of.Model.without.SNP", "Rsquare.of.Model.with.SNP", "FDR_Adjusted_P-values","effect","file.name","markers")
signif <- data_rbind[(data_rbind[,4]<0.000005),]
write.csv(signif,'SigSNPs_07212020Comp.csv')
rm(list=ls())