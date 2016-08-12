exp_dat <- read.csv("~/Documents/HEK293expression/Data/rna_celline.csv", header=TRUE)
gene_list <- read.csv("~/Documents/HEK293expression/Data/TADA_top65ASDgenes_June2016.csv", header=TRUE)
head(exp_dat)
head(gene_list)
n <- nrow(gene_list)

HEK293_dat <- exp_dat[exp_dat$Sample == "HEK 293", c("Gene.name", "Sample", "Value")]

gene_fpkm <- merge(gene_list, HEK293_dat, by.x = "RefSeqName", by.y = "Gene.name")
head(gene_fpkm)
gene_fpkm <- gene_fpkm[, c("RefSeqName", "AltName", "Value", "qvalue.combined")]
colnames(gene_fpkm) <- c("GeneID", "AltGeneID", "FPKM", "FDR")
head(gene_fpkm)
gene_fpkm <- gene_fpkm[order(-gene_fpkm$FPKM),]
head(gene_fpkm)
write.table(gene_fpkm, "HEK293ExpTopAsdGenes.txt", quote = FALSE, col.names = TRUE, sep = "\t", row.names=FALSE)

#plot q/FDR and HEK fpkm 
plot((-log10(gene_fpkm$FDR)), gene_fpkm$FPKM, xlab = "-log10 of FDR", ylab = "Gene expression value (FPKM)", main = "Top 65 ASD Genes" , type = "p")

#for only dna reg genes
reg_genes <- read.csv("~/Documents/HEK293expression/Data/TADA_top179ASDgenes_July2016_annotatedByFunction.csv", header=TRUE)
head(reg_genes)
n <- nrow(reg_genes)
reg_fpkm <- merge(reg_genes, HEK293_dat, by.x = "ID", by.y = "Gene.name")
head(reg_fpkm)
reg_fpkm <- reg_fpkm[, c("ID", "Gene.Name", "Value", "TADA.FDR")]
colnames(reg_fpkm) <- c("GeneID", "Func", "FPKM", "FDR")
head(reg_fpkm)
reg_fpkm <- reg_fpkm[order(-reg_fpkm$FPKM),]
head(reg_fpkm)
write.table(reg_fpkm, "HEK293ExpTopRegAsdGenes.txt", quote = FALSE, col.names = TRUE, sep = "\t", row.names=F)

#where FPKM is at least 6
reg_fpkm6 <- reg_fpkm[reg_fpkm$FPKM > 6,]
tail(reg_fpkm6)
write.table(reg_fpkm6, "HEK293ExpTopRegAsdGenesGT6.txt", quote=FALSE, col.names=TRUE, sep = "\t", row.names=F)

#plot q/FDR and HEK fpkm 
plot((-log10(gene_fpkm$FDR)), gene_fpkm$FPKM, xlab = "-log10 of FDR", ylab = "Gene expression value (FPKM)", main = "Top 65 ASD Genes" , type = "p")
plot((-log10(reg_fpkm$FDR)), reg_fpkm$FPKM, xlab = "-log10 of FDR", ylab = "Gene expression value (FPKM)", main = "Top DNA Reg ASD Genes" , type = "h")

#plot with log of fpkm
plot((-log10(gene_fpkm$FDR)), (log(gene_fpkm$FPKM)), xlab = "-log10 of FDR", ylab = "log Gene expression value (FPKM)", main = "Top 65 ASD Genes" , type = "p")
plot((-log10(reg_fpkm$FDR)), (log(reg_fpkm$FPKM)), xlab = "-log10 of FDR", ylab = "log Gene expression value (FPKM)", main = "Top DNA Reg ASD Genes" , type = "p")


#check correlation
cor(gene_fpkm$FDR, gene_fpkm$FPKM)
cor(reg_fpkm$FDR, reg_fpkm$FPKM)


#create heat map
  sample.types <- unique(exp_dat$Sample)
  sample <- exp_dat[exp_dat$Sample == sample.types[1], c("Gene.name", "Value")]
  asd_hmap <- merge(gene_list, sample, by.x = "RefSeqName", by.y = "Gene.name")
  j <- 4
  colnames(asd_hmap)[j] <- sample.types[1]
  head(asd_hmap)


sample.types <- unique(exp_dat$Sample)
sample <- exp_dat[exp_dat$Sample == sample.types[1], c("Gene.name", "Value")]
asd_hmap <- merge(gene_list, sample, by = intersect("RefSeqName", "Gene.name"), all.x = FALSE, all.y = FALSE)
                  
                  
                  
                  = "RefSeqName", by.y = "Gene.name")

for (i in 6:length(sample.types)) {
  sample <- exp_dat[exp_dat$Sample == sample.types[i], c("Gene.name", "Value")]
  asd_hmap <- merge(asd_hmap, sample, by.x = "RefSeqName", by.y = "Gene.name")
  j <- 3 + i
  colnames(asd_hmap)[j] <- i
  head(asd_hmap)
}












j=3

reform <- function(sampletypes) {
  sample <- exp_dat[exp_dat$Sample == sample.types[i], c("Gene.name", "Value")]
  asd_hmap <- merge(gene_list, sample, by.x = "RefSeqName", by.y = "Gene.name")
}

reform(sample.types)


