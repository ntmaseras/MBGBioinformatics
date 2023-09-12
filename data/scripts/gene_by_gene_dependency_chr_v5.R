gene_selected <- "DYRK1A"

setwd("~/Desktop/P3 miR-9 Follow up/TCGA_miRNA_mutants/DICER1 depmap/original_data/")
require(data.table)
achilles <- fread(file = "20Q4v2_Achilles_gene_effect.csv", header = T)
colnames(achilles)[1] <- "Depmap.Id"

colunm_names <- colnames(achilles) 
colunm_names<- gsub(" ","_",colunm_names)
colunm_names<- gsub("\\(","",colunm_names)
colunm_names<- gsub("\\)","",colunm_names)
colnames(achilles)  <- colunm_names 
rm(colunm_names)

#install.packages("reshape2")
library(reshape2)
achilles2 <- melt(achilles, id.vars = "Depmap.Id")

names <- as.data.frame(colnames(achilles))
names$gene_cor <- NA
names <- names[c(2:18120),]
names$gene_p <- NA
colnames(names)[1] <- "variable"
names$variable <- substr(names$variable,1,regexpr("_",names$variable)-1)

achilles2$variable <- substr(achilles2$variable,1,regexpr("_",achilles2$variable)-1)

gene_set <- achilles2[which(achilles2$variable==gene_selected),]
colnames(gene_set) <- c("Depmap.Id","gene_selected","value_gene_selected")

names_empty <- names
rm(achilles)

names <- names_empty

names$gene_selected <- gene_selected

i <- 1
for (i in 1:nrow(names)) {
  gene_test <- names$variable[i]
  
  test_set <- achilles2[which(achilles2$variable==gene_test),]
  test_set <- merge(test_set,gene_set,by="Depmap.Id")
  pearson_test <- cor.test(test_set$value_gene_selected, test_set$value, method = "pearson")
  names$gene_cor[i] <- as.numeric(unlist(pearson_test)[4])
  p_test <- as.numeric(unlist(pearson_test)[3])
  names$gene_p[i] <- -log10(p_test)

  print(paste0(i/nrow(names)*100, " ", gene_selected))
  rm(gene_test,test_set,p_test,pearson_test)
}


plot(names$gene_cor,names$gene_p)

ArrayAtt <- fread(file = "Achilles_guide_map.csv", header = T)
ArrayAtt$genome_alignment <- substr(ArrayAtt$genome_alignment, 1,regexpr("_",ArrayAtt$genome_alignment)-1)
ArrayAtt$gene <- substr(ArrayAtt$gene, 1,regexpr(" ",ArrayAtt$gene)-1)
ArrayAtt <- unique(ArrayAtt[,c(3,2)])
ArrayAtt$chromosome_name <- substr(ArrayAtt$genome_alignment,4,nchar(ArrayAtt$genome_alignment))
ArrayAtt$chromosome_name <- gsub("X","24",ArrayAtt$chromosome_name)
ArrayAtt$chromosome_name <- as.numeric(ArrayAtt$chromosome_name)


colnames(names)[1] <- "gene"
names <- merge(names, ArrayAtt, by="gene")


names <- names[which(names$gene_p>=-log10(0.05)),]

names <- names[order(-names$gene_p,-names$chromosome_name),]
plot(names$chromosome_name,names$gene_p)
colnames(names)[3] <- paste0(gene_selected,"_p")

setwd("~/Desktop/P3 miR-9 Follow up/TCGA_miRNA_mutants/DICER1 depmap/manhattan_gene_dependency_correlations/")
write.table(x = names, file = paste0("correlations_",gene_selected,"_by_chr.tsv"),row.names = F,col.names = T, sep = "\t")
setwd("~/Desktop/P3 miR-9 Follow up/TCGA_miRNA_mutants/DICER1 depmap/")
