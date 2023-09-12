#set directory
local_directory <- "~/Desktop/P3 miR-9 Follow up/TCGA_miRNA_mutants/DICER1 depmap/original_data/"
local_output <- "~/Desktop/P3 miR-9 Follow up/TCGA_miRNA_mutants/DICER1 depmap/original_data/"
cloud_start <- "/sbgenomics/project-files/"
cloud_end <- "/sbgenomics/output-files"

input <- local_directory
output <- local_output
setwd(input)

require(data.table)
ArrayAtt <- fread(file = "Achilles_guide_map.csv", header = T)
ArrayAtt$position <- substr(ArrayAtt$genome_alignment, regexpr("_",ArrayAtt$genome_alignment)+1, nchar(ArrayAtt$genome_alignment))
ArrayAtt$position <- as.numeric(substr(ArrayAtt$position, 1, regexpr("_",ArrayAtt$position)-1))
ArrayAtt$genome_alignment <- substr(ArrayAtt$genome_alignment, 1,regexpr("_",ArrayAtt$genome_alignment)-1)
ArrayAtt$gene <- substr(ArrayAtt$gene, 1,regexpr(" ",ArrayAtt$gene)-1)
chr <- unique(ArrayAtt[,c(3,2)])
pos <- aggregate(ArrayAtt$position, by=list(gene=ArrayAtt$gene),mean)
colnames(pos) <- c("gene","avg_pos")
ArrayAtt <- merge(chr, pos,by="gene")
ArrayAtt$chromosome_name <- substr(ArrayAtt$genome_alignment,4,nchar(ArrayAtt$genome_alignment))
ArrayAtt$chromosome_name <- gsub("X","24",ArrayAtt$chromosome_name)
ArrayAtt$chromosome_name <- as.numeric(ArrayAtt$chromosome_name)
rm(pos,chr)

achilles <- fread(file = "20Q4v2_Achilles_gene_effect.csv", header = T)
colnames(achilles)[1] <- "Depmap.Id"

colunm_names <- colnames(achilles)
colunm_names<- gsub(" ","_",colunm_names)
colunm_names<- gsub("\\(","",colunm_names)
colunm_names<- gsub("\\)","",colunm_names)
colnames(achilles)  <- colunm_names 
gene_names <- colunm_names[2:length(colunm_names)]
gene_names <- as.data.frame(gene_names)
colnames(gene_names) <- "gene"
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

#gene_names$gene <- substr(gene_names$gene,1,regexpr("_",gene_names$gene)-1)
#gene_names$tested <- F
#write.table(x = gene_names, file = "list_of_genes_completed.tsv",row.names = F,col.names = T, sep = "\t")
gene_names <- fread(file = "list_of_genes_completed.tsv", header = T)
setwd(output)
write.table(x = gene_names, file = "list_of_genes_completed.tsv",row.names = F,col.names = T, sep = "\t")

gene_names <- gene_names[gene_names$gene %in% achilles2$variable,]
ArrayAtt <- ArrayAtt[ArrayAtt$gene %in% achilles2$variable,]

sel_chr <- "chr21"
selected_chr <- ArrayAtt[which(ArrayAtt$genome_alignment==sel_chr),]

gene_names <- gene_names[which(gene_names$tested==F),]
gene_names <- gene_names[gene_names$gene %in% selected_chr$gene,]

#downsampling
set.seed(123)
index <- sample(x = 1:nrow(gene_names), size = 5)
gene_names <- gene_names[index, ]
#specific gene list
# index <- c("DICER1","AGO2","TUT4")
# gene_names <- gene_names[gene_names$gene %in% index, ]
rm(index)

rm(achilles)
names_empty <- names
gene_names$time <- NA
a <- 1
for (a in 1:nrow(gene_names)) {
  start.time <- Sys.time()
  
  gene_selected <- gene_names$gene[a]
  gene_set <- achilles2[which(achilles2$variable==gene_selected),c(3)] #vectorized gene_set
  
  names <- names_empty
  names$gene_selected <- gene_selected
  
  gene_cor <- character (nrow(names)) # initialize output vector
  gene_p <- character (nrow(names)) # initialize output vector
  
  i <- 1
  for (i in 1:nrow(names)) {
    gene_test <- names$variable[i]
    test_set <- achilles2[which(achilles2$variable==gene_test),c(3)] #vectorized test_set
    pearson_test <- cor.test(gene_set, test_set, method = "pearson")
    
    gene_cor[i] <- as.numeric(unlist(pearson_test)[4])
    gene_p[i] <- -log10(as.numeric(unlist(pearson_test)[3]))
    
    print(paste0(i/nrow(names)*100, " ", gene_selected," ",a,"/",nrow(gene_names)))
    rm(gene_test,test_set,p_test,pearson_test)
  }
  end.time <- Sys.time()
  names$gene_cor <- gene_cor  # finally assign to data frame
  names$gene_p <- gene_p  # finally assign to data frame
  
  #plot(names$gene_cor,names$gene_p)
  colnames(names)[1] <- "gene"
  names <- merge(names, ArrayAtt, by="gene")
  names$gene_p <- as.numeric(names$gene_p)
  names <- names[which(names$gene_p>=-log10(0.05)),]
  names <- names[order(-names$gene_p,-names$chromosome_name),]
  #plot(names$chromosome_name,names$gene_p)
  colnames(names)[3] <- paste0(gene_selected,"_p")
  
  setwd(output)
  gene_names$time[a] <- difftime(end.time, start.time, units = "mins")
  gene_names$tested[a] <- T
  write.table(x = names, file = paste0("correlations_",gene_selected,"_by_chr.tsv"),row.names = F,col.names = T, sep = "\t")
  write.table(x = gene_names, file = "session_times.tsv",row.names = F,col.names = T, sep = "\t",append = F)
  
  setwd(output)
  gene_done <- fread(file = "list_of_genes_completed.tsv", header = T)
  gene_done[which(gene_done$gene==gene_selected),2] <- T
  write.table(x = gene_done, file = "list_of_genes_completed.tsv",row.names = F,col.names = T, sep = "\t")
  rm(gene_done)
  rm(gene_selected)
}