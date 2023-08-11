library(DESeq2)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(VennDiagram)
library(pheatmap)

#Load counts
rawcounts<-read.csv("/gpfs/data/giannarellilab/Covidproject/RNAseq/Macrophages/DESeq/counts_notinfected.csv", header = F)
#for some reason row 8686 has no counts...
#Load metadata
sample_table<-read.csv("/gpfs/data/giannarellilab/Covidproject/RNAseq/Macrophages/DESeq/metadata_notinfected.csv")
rawcounts<-na.omit(rawcounts)
colnames(rawcounts) <- rawcounts[1,]
rawcounts <- rawcounts[-1,]

names <- make.unique(rawcounts$Geneid)
rownames(rawcounts) <- names
rawcounts <- rawcounts[,-1] 

sample_table$Donor<-as.factor(sample_table$Donor)
sample_table$Cell<-as.factor(sample_table$Cell)
sample_table$Condition<-as.factor(sample_table$Condition)
sample_table$Timepoint<-as.factor(sample_table$Timepoint)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"), verbose = T)
genes <- data.frame(rownames(rawcounts), ensembl = rownames(rawcounts))
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
annot <- merge(
  x = genelist,
  y = genes,
  by.x="ensembl_gene_id",
  by.y="ensembl")

################################################################################################################################
#Figure 2G
#Infected macrophages and foam cells vs their respective non-infected samples. Timepoint included as co-variate.
##macrophages
counts_mac <- rawcounts[,which(sample_table$Cell == 'primary mac' & sample_table$Timepoint != 0)]
metadata_mac <- sample_table[which(sample_table$Cell == 'primary mac' & sample_table$Timepoint != 0),]
metadata_mac$Condition_Timepoint <- paste(metadata_mac$Condition,metadata_mac$Timepoint, sep = '_')
counts_mac <- counts_mac[-8686,]
counts_mac <- mutate_all(counts_mac, function(x) as.numeric(as.character(x)))

dds_mac<- DESeqDataSetFromMatrix(countData = counts_mac, colData = metadata_mac, design = ~Condition+Timepoint+Donor)
dds_mac<- DESeq(dds_mac)
resultsNames(dds_mac)
res <- results(dds_mac, name = "Condition_SARSCov2.infected_vs_not.infected", cooksCutoff = Inf)
res <- as.data.frame(res)
genes <- rownames(res)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
res$genes<-rownames(res)
annot <- merge(x = genelist, y = res, by.x="ensembl_gene_id", by.y="genes", all.y = T)
length(which(res$padj < 0.05))
res_mac <- annot
write.csv(res_mac[which(res_mac$padj < 0.05),], file = 'inf_vs_notinf_mac_sig.csv')

vsd <- vst(dds_mac, blind=FALSE)
pcaData<-plotPCA(vsd, intgroup=c("Cell", "Condition","Timepoint", "Donor", "Batch"),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Timepoint)) +
  geom_point(size=3)
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Donor)) +
  geom_point(size=3)
ggplot(pcaData, aes(PC1, PC2, color=Batch, shape=Timepoint)) +
  geom_point(size=3)

## foam cells
counts_foam <- rawcounts[,which(sample_table$Cell == 'foam cells' & sample_table$Timepoint != 0)]
metadata_foam <- sample_table[which(sample_table$Cell == 'foam cells' & sample_table$Timepoint != 0),]
metadata_foam$Condition_Timepoint <- paste(metadata_foam$Condition,metadata_foam$Timepoint, sep = '_')
counts_foam <- counts_foam[-8686,]
counts_foam <- mutate_all(counts_foam, function(x) as.numeric(as.character(x)))

dds_foam <- DESeqDataSetFromMatrix(countData = counts_foam, colData = metadata_foam, design = ~Condition+Timepoint+Donor)
dds_foam <- DESeq(dds_foam)
resultsNames(dds_foam)
res <- results(dds_foam, name = "Condition_SARSCov2.infected_vs_not.infected", cooksCutoff = Inf)
res <- as.data.frame(res)
genes <- rownames(res)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
res$genes<-rownames(res)
annot <- merge(x = genelist, y = res, by.x="ensembl_gene_id", by.y="genes", all.y = T)
length(which(res$padj < 0.05))
res_foam <- annot
write.csv(res_foam[which(res_foam$padj < 0.05),], file = 'inf_vs_notinf_foam_sig.csv')

vsd <- vst(dds_foam, blind=FALSE)
pcaData<-plotPCA(vsd, intgroup=c("Cell", "Condition","Timepoint", "Donor", "Batch"),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Timepoint, shape=Condition)) +
  geom_point(size=3)
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Donor)) +
  geom_point(size=3)
ggplot(pcaData, aes(PC1, PC2, color=Batch, shape=Timepoint)) +
  geom_point(size=3)

## Venn
foam_genes <- res_foam$ensembl_gene_id[which(res_foam$padj < 0.05)]
mac_genes <- res_mac$ensembl_gene_id[which(res_mac$padj < 0.05)]
venn.diagram(
  x = list(na.omit(foam_genes), na.omit(mac_genes)),
  category.names = c("Foam genes" , "mac genes"),
  output=F,
  filename = 'mac_foam_venn.png'
)
#Identify and save unique and shared genes
shared <- intersect(foam_genes, mac_genes)
mac_shared <- res_mac[which(res_mac$ensembl_gene_id %in% shared),]
foam_shared <- res_foam[which(res_foam$ensembl_gene_id %in% shared),]
shared <- cbind(mac_shared, foam_shared)
col <- c(rep('macrophages', ncol(mac_shared)), rep('foam cells', ncol(foam_shared)))
colnames(shared) <- paste(colnames(shared), col, sep = '_')
write.csv(shared, file = 'inf_vs_notinf_shared.csv', quote = F)
write.csv(res_foam[which(res_foam$ensembl_gene_id %in% foam_genes[which(!foam_genes %in% mac_genes)]),], file = 'inf_vs_notinf_foam_specific.csv', quote = F)
write.csv(res_mac[which(res_mac$ensembl_gene_id %in% mac_genes[which(!mac_genes %in% foam_genes)]),], file = 'inf_vs_notinf_mac_specific.csv', quote = F)

################################################################################################################################
#Figures 2E, 2I, Extended 3F, Extended 3I
#Infected macrophages vs infected foam cells
##DESeq comparisons
counts_infected <- rawcounts[,which(sample_table$Condition == 'SARSCov2 infected' | sample_table$Timepoint == '0')]
metadata_infected <- sample_table[which(sample_table$Condition == 'SARSCov2 infected' | sample_table$Timepoint == '0'),]
counts_infected <- counts_infected[-8686,]
counts_infected <- mutate_all(counts_infected, function(x) as.numeric(as.character(x)))
metadata_infected$Cell_Timepoint <- paste(metadata_infected$Cell,metadata_infected$Timepoint, sep = '_')

metadata_infected$Cell_Timepoint <- factor(metadata_infected$Cell_Timepoint)
dds_infected <- DESeqDataSetFromMatrix(countData = counts_infected, colData = metadata_infected, design = ~Cell_Timepoint+Donor)
dds_infected <- DESeq(dds_infected)
resultsNames(dds_infected)

res_0 <- results(dds_infected, name = "Cell_Timepoint_primary.mac_0_vs_foam.cells_0")
res_0 <- lfcShrink(dds_infected, coef = "Cell_Timepoint_primary.mac_0_vs_foam.cells_0", res = res_0)
res_0 <- as.data.frame(res_0)
genes <- rownames(res_0)
res_0$genes<-rownames(res_0)
annot <- merge(x = genelist, y = res_0, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_0<- annot

metadata_infected$Cell_Timepoint <- relevel(metadata_infected$Cell_Timepoint, ref = 'foam cells_2')
dds_infected <- DESeqDataSetFromMatrix(countData = counts_infected, colData = metadata_infected, design = ~Cell_Timepoint+Donor)
dds_infected <- DESeq(dds_infected)
resultsNames(dds_infected)
res_2 <- results(dds_infected, name = "Cell_Timepoint_primary.mac_2_vs_foam.cells_2")
res_2 <- lfcShrink(dds_infected, coef = "Cell_Timepoint_primary.mac_2_vs_foam.cells_2", res = res_2)
res_2 <- as.data.frame(res_2)
genes <- rownames(res_2)
res_2$genes<-rownames(res_2)
annot <- merge(x = genelist, y = res_2, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_2<- annot

metadata_infected$Cell_Timepoint <- relevel(metadata_infected$Cell_Timepoint, ref = 'foam cells_8')
dds_infected <- DESeqDataSetFromMatrix(countData = counts_infected, colData = metadata_infected, design = ~Cell_Timepoint+Donor)
dds_infected <- DESeq(dds_infected)
resultsNames(dds_infected)
res_8 <- results(dds_infected, name = "Cell_Timepoint_primary.mac_8_vs_foam.cells_8")
res_8 <- lfcShrink(dds_infected, coef = "Cell_Timepoint_primary.mac_8_vs_foam.cells_8", res = res_8)
res_8 <- as.data.frame(res_8)
genes <- rownames(res_8)
res_8$genes<-rownames(res_8)
annot <- merge(x = genelist, y = res_8, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_8<- annot

metadata_infected$Cell_Timepoint <- relevel(metadata_infected$Cell_Timepoint, ref = 'foam cells_24')
dds_infected <- DESeqDataSetFromMatrix(countData = counts_infected, colData = metadata_infected, design = ~Cell_Timepoint+Donor)
dds_infected <- DESeq(dds_infected)
resultsNames(dds_infected)
res_24 <- results(dds_infected, name = "Cell_Timepoint_primary.mac_24_vs_foam.cells_24")
res_24 <- lfcShrink(dds_infected, coef = "Cell_Timepoint_primary.mac_24_vs_foam.cells_24", res = res_24)
res_24 <- as.data.frame(res_24)
genes <- rownames(res_24)
res_24$genes<-rownames(res_24)
annot <- merge(x = genelist, y = res_24, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_24<- annot

metadata_infected$Cell_Timepoint <- relevel(metadata_infected$Cell_Timepoint, ref = 'foam cells_48')
dds_infected <- DESeqDataSetFromMatrix(countData = counts_infected, colData = metadata_infected, design = ~Cell_Timepoint+Donor)
dds_infected <- DESeq(dds_infected)
resultsNames(dds_infected)
res_48 <- results(dds_infected, name = "Cell_Timepoint_primary.mac_48_vs_foam.cells_48")
res_48 <- lfcShrink(dds_infected, coef = "Cell_Timepoint_primary.mac_48_vs_foam.cells_48", res = res_48)
res_48 <- as.data.frame(res_48)
genes <- rownames(res_48)
res_48$genes<-rownames(res_48)
annot <- merge(x = genelist, y = res_48, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_48<- annot

###Write ouput
write.csv(file = 'infected_mac_vs_foam/Time_0_sig_0.1FDR.csv', res_0[which(res_0$padj < 0.1),])
write.csv(file = 'infected_mac_vs_foam/Time_2_sig_0.1FDR.csv', res_2[which(res_2$padj < 0.1),])
write.csv(file = 'infected_mac_vs_foam/Time_8_sig_0.1FDR.csv', res_8[which(res_8$padj < 0.1),])
write.csv(file = 'infected_mac_vs_foam/Time_24_sig_0.1FDR.csv', res_24[which(res_24$padj < 0.1),])
write.csv(file = 'infected_mac_vs_foam/Time_48_sig_0.1FDR.csv', res_48[which(res_48$padj < 0.1),])

###load gene sets
covid <- c('E', 'S', 'orf8', 'orf7a', 'orf6', 'orf3a', 'orf1ab', 'orf10', 'N', 'M')
interferon<-c("JAK1", "IRF1", "IFNAR2", "IFNA1", "IFNAR1", "IRF4", "SOCS3", "IRF5", "SOCS1", "STAT2", "MX1", "IFI35", "IRF7", "BST2", "IFITM1", "STAT1", "TYK2", "PSMB8", "PTPN6", "IFNA7","IFNA2", "IFNG", "IFNL1", "IFNB1","IFNA6", "IFNA10", "IFNA16", "IFNLR1", "ZBP1", "NFKB1", "MAVS", "IRF3")
remove <- c("IFNG", "IFNA6", "IFNA16", "IFNL1", "IFNA10", "IFNA2", "IFNA1", "IFNA7", "IFNB1")
interferon_revised <- interferon[-which(interferon %in% remove)]
comp_genes <- c('CD46','CD55','CD59','CD93','CFH','CFHR1','CFHR2','CFHR3','CFHR4','CFHR5','CFI','CLU','CR1','CR2','CSMD1','CSMD2','CSMD3','C3AR1','C4BPA','C4BPB','C5AR1','ELANE','F2','ITGAM','ITGAX','ITGB2','SERPING1','VSIG4','CFB','CFD','CFP','COLEC10','COLEC11','C1QA','C1QB','C1QC','C1R','C1S','C2','C3','C4A','C4B','C5','C6','C7','C8A','C8B','C8G','C9','FCN1','FCN2','FCN3','MASP1','MASP2','MBL2')
lipids <- c("CD36", "SCARB1", "MRC1", "CD68", "TREM2", "ABCA1", "ABCG1", "LDLR", "PPARG")

###Make heatmaps
fc_int <- data.frame(int_0 = res_0$log2FoldChange[which(res_0$hgnc_symbol %in% interferon_revised)],int_2 = res_2$log2FoldChange[which(res_2$hgnc_symbol %in% interferon_revised)],int_8 = res_8$log2FoldChange[which(res_8$hgnc_symbol %in% interferon_revised)], int_24 = res_24$log2FoldChange[which(res_24$hgnc_symbol %in% interferon_revised)], int_48 = res_48$log2FoldChange[which(res_48$hgnc_symbol %in% interferon_revised)], row.names = res_2$hgnc_symbol[which(res_2$hgnc_symbol %in% interferon_revised)])
fc_int[is.na(fc_int)] <- 0

p_int <- data.frame(int_0 = res_0$padj[which(res_0$hgnc_symbol %in% interferon_revised)],int_2 = res_2$padj[which(res_2$hgnc_symbol %in% interferon_revised)],int_8 = res_8$padj[which(res_8$hgnc_symbol %in% interferon_revised)], int_24 = res_24$padj[which(res_24$hgnc_symbol %in% interferon_revised)], int_48 = res_48$padj[which(res_48$hgnc_symbol %in% interferon_revised)], row.names = res_2$hgnc_symbol[which(res_2$hgnc_symbol %in% interferon_revised)])
p_int[p_int < 0.1] <- "*"
p_int[p_int != "*"] <- ""
p_int[is.na(p_int)] <- ""

fc_covid <- data.frame(int_0 = res_0$log2FoldChange[which(res_0$ensembl_gene_id %in% covid)],int_2 = res_2$log2FoldChange[which(res_2$ensembl_gene_id %in% covid)],int_8 = res_8$log2FoldChange[which(res_8$ensembl_gene_id %in% covid)], int_24 = res_24$log2FoldChange[which(res_24$ensembl_gene_id %in% covid)], int_48 = res_48$log2FoldChange[which(res_48$ensembl_gene_id %in% covid)], row.names = res_2$ensembl_gene_id[which(res_2$ensembl_gene_id %in% covid)])
fc_covid[is.na(fc_covid)] <- 0

p_covid <- data.frame(int_0 = res_0$padj[which(res_0$ensembl_gene_id %in% covid)],int_2 = res_2$padj[which(res_2$ensembl_gene_id %in% covid)],int_8 = res_8$padj[which(res_8$ensembl_gene_id %in% covid)], int_24 = res_24$padj[which(res_24$ensembl_gene_id %in% covid)], int_48 = res_48$padj[which(res_48$ensembl_gene_id %in% covid)], row.names = res_2$ensembl_gene_id[which(res_2$ensembl_gene_id %in% covid)])
p_covid[p_covid < 0.1] <- "*"
p_covid[p_covid != "*"] <- ""
p_covid[is.na(p_covid)] <- ""

fc_comp <- data.frame(int_0 = res_0$log2FoldChange[which(res_0$hgnc_symbol %in% comp_genes)],int_2 = res_2$log2FoldChange[which(res_2$hgnc_symbol %in% comp_genes)],int_8 = res_8$log2FoldChange[which(res_8$hgnc_symbol %in% comp_genes)], int_24 = res_24$log2FoldChange[which(res_24$hgnc_symbol %in% comp_genes)], int_48 = res_48$log2FoldChange[which(res_48$hgnc_symbol %in% comp_genes)], row.names = res_2$hgnc_symbol[which(res_2$hgnc_symbol %in% comp_genes)])
fc_comp[is.na(fc_comp)] <- 0

p_comp <- data.frame(int_0 = res_0$padj[which(res_0$hgnc_symbol %in% comp_genes)],int_2 = res_2$padj[which(res_2$hgnc_symbol %in% comp_genes)],int_8 = res_8$padj[which(res_8$hgnc_symbol %in% comp_genes)], int_24 = res_24$padj[which(res_24$hgnc_symbol %in% comp_genes)], int_48 = res_48$padj[which(res_48$hgnc_symbol %in% comp_genes)], row.names = res_2$hgnc_symbol[which(res_2$hgnc_symbol %in% comp_genes)])
p_comp <- p_comp[rowSums(is.na(p_comp)) != 5,]
p_comp[p_comp < 0.1] <- "*"
p_comp[p_comp != "*"] <- ""
p_comp[is.na(p_comp)] <- ""
fc_comp <- fc_comp[which(row.names(fc_comp) %in% row.names(p_comp)),]

fc_lipids <- data.frame(int_0 = res_0$log2FoldChange[which(res_0$hgnc_symbol %in% lipids)],int_2 = res_2$log2FoldChange[which(res_2$hgnc_symbol %in% lipids)],int_8 = res_8$log2FoldChange[which(res_8$hgnc_symbol %in% lipids)], int_24 = res_24$log2FoldChange[which(res_24$hgnc_symbol %in% lipids)], int_48 = res_48$log2FoldChange[which(res_48$hgnc_symbol %in% lipids)], row.names = res_2$hgnc_symbol[which(res_2$hgnc_symbol %in% lipids)])
fc_lipids[is.na(fc_lipids)] <- 0

p_lipids <- data.frame(int_0 = res_0$padj[which(res_0$hgnc_symbol %in% lipids)],int_2 = res_2$padj[which(res_2$hgnc_symbol %in% lipids)],int_8 = res_8$padj[which(res_8$hgnc_symbol %in% lipids)], int_24 = res_24$padj[which(res_24$hgnc_symbol %in% lipids)], int_48 = res_48$padj[which(res_48$hgnc_symbol %in% lipids)], row.names = res_2$hgnc_symbol[which(res_2$hgnc_symbol %in% lipids)])
p_lipids[p_lipids < 0.1] <- "*"
p_lipids[p_lipids != "*"] <- ""
p_lipids[is.na(p_lipids)] <- ""

#Figure 2E
myBreaks <- c(seq(-0.5, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, 0.5, length.out=floor(paletteLength/2)))
pdf(file = 'infected_mac_vs_foam/mac_vs_foam_covid_0.1FDR.pdf', height = 4, width = 3)
pheatmap(fc_covid, cluster_cols = F, breaks = myBreaks, color = myColor, display_numbers = p_covid)
dev.off()
#Figure 2I
pdf(file = 'infected_mac_vs_foam/mac_vs_foam_Interferon_0.1FDR.pdf', height = 6, width = 3)
pheatmap(fc_int, cluster_cols = F, breaks = myBreaks, color = myColor, display_numbers = p_int)
dev.off()
#Extended 3F
myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, 1, length.out=floor(paletteLength/2)))
pdf(file = 'infected_mac_vs_foam/mac_vs_foam_compliment_0.1FDR.pdf', height = 6, width = 3)
pheatmap(fc_comp, cluster_cols = F, breaks = myBreaks, color = myColor, display_numbers = p_comp)
dev.off()
#Extended 3I
pdf(file = 'infected_mac_vs_foam/mac_vs_foam_lipids_0.1FDR.pdf', height = 6, width = 3)
pheatmap(fc_lipids, cluster_cols = F, breaks = myBreaks, color = myColor, display_numbers = p_lipids)
dev.off()


################################################################################################################################
# Figure 2H
## mac
counts_mac <- rawcounts[,which(sample_table$Cell == 'primary mac' & sample_table$Timepoint != 0)]
metadata_mac <- sample_table[which(sample_table$Cell == 'primary mac' & sample_table$Timepoint != 0),]
metadata_mac$Condition_Timepoint <- paste(metadata_mac$Condition,metadata_mac$Timepoint, sep = '_')
counts_mac <- counts_mac[-8686,]
counts_mac <- mutate_all(counts_mac, function(x) as.numeric(as.character(x)))

dds_mac<- DESeqDataSetFromMatrix(countData = counts_mac, colData = metadata_mac, design = ~Condition+Timepoint+Donor+Condition:Timepoint)
dds_mac<- DESeq(dds_mac)
resultsNames(dds_mac)
res <- results(dds_mac, name = "ConditionSARSCov2.infected.Timepoint8" , cooksCutoff = Inf)
res <- as.data.frame(res)
genes <- rownames(res)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
res$genes<-rownames(res)
annot <- merge(x = genelist, y = res, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_mac_int_8 <- annot
res <- results(dds_mac, name = "ConditionSARSCov2.infected.Timepoint24" , cooksCutoff = Inf)
res <- as.data.frame(res)
genes <- rownames(res)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
res$genes<-rownames(res)
annot <- merge(x = genelist, y = res, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_mac_int_24 <- annot
res <- results(dds_mac, name = "ConditionSARSCov2.infected.Timepoint48" , cooksCutoff = Inf)
res <- as.data.frame(res)
genes <- rownames(res)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
res$genes<-rownames(res)
annot <- merge(x = genelist, y = res, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_mac_int_48 <- annot


p_interactions_int <- data.frame(int_8 = res_mac_int_8$padj[which(res_mac_int_8$hgnc_symbol %in% interferon_revised)],int_24 = res_mac_int_24$padj[which(res_mac_int_24$hgnc_symbol %in% interferon_revised)],int_48 = res_mac_int_48$padj[which(res_mac_int_48$hgnc_symbol %in% interferon_revised)], row.names = res_mac_int_8$hgnc_symbol[which(res_mac_int_8$hgnc_symbol %in% interferon_revised)])

metadata_mac$Condition_Timepoint <- factor(metadata_mac$Condition_Timepoint)
dds_mac<- DESeqDataSetFromMatrix(countData = counts_mac, colData = metadata_mac, design = ~Condition_Timepoint+Donor)
dds_mac<- DESeq(dds_mac)
resultsNames(dds_mac)

res_2 <- results(dds_mac, name = "Condition_Timepoint_SARSCov2.infected_2_vs_not.infected_2" , cooksCutoff = Inf)
res_2 <- lfcShrink(dds_mac, coef = "Condition_Timepoint_SARSCov2.infected_2_vs_not.infected_2" , res = res_2)
res_2 <- as.data.frame(res_2)
genes <- rownames(res_2)
res_2$genes<-rownames(res_2)
annot <- merge(x = genelist, y = res_2, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_2 <- annot

metadata_mac$Condition_Timepoint <- relevel(metadata_mac$Condition_Timepoint, ref = 'not infected_8')
dds_mac<- DESeqDataSetFromMatrix(countData = counts_mac, colData = metadata_mac, design = ~Condition_Timepoint+Donor)
dds_mac<- DESeq(dds_mac)
res_8 <- results(dds_mac, name = "Condition_Timepoint_SARSCov2.infected_8_vs_not.infected_8" , cooksCutoff = Inf)
res_8 <- lfcShrink(dds_mac, coef = "Condition_Timepoint_SARSCov2.infected_8_vs_not.infected_8" , res = res_8)
res_8 <- as.data.frame(res_8)
genes <- rownames(res_8)
res_8$genes<-rownames(res_8)
annot <- merge(x = genelist, y = res_8, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_8 <- annot

metadata_mac$Condition_Timepoint <- relevel(metadata_mac$Condition_Timepoint, ref = 'not infected_24')
dds_mac<- DESeqDataSetFromMatrix(countData = counts_mac, colData = metadata_mac, design = ~Condition_Timepoint+Donor)
dds_mac<- DESeq(dds_mac)
res_24 <- results(dds_mac, name = "Condition_Timepoint_SARSCov2.infected_24_vs_not.infected_24" , cooksCutoff = Inf)
res_24 <- lfcShrink(dds_mac, coef = "Condition_Timepoint_SARSCov2.infected_24_vs_not.infected_24" , res = res_24)
res_24 <- as.data.frame(res_24)
genes <- rownames(res_24)
res_24$genes<-rownames(res_24)
annot <- merge(x = genelist, y = res_24, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_24 <- annot

metadata_mac$Condition_Timepoint <- relevel(metadata_mac$Condition_Timepoint, ref = 'not infected_48')
dds_mac<- DESeqDataSetFromMatrix(countData = counts_mac, colData = metadata_mac, design = ~Condition_Timepoint+Donor)
dds_mac<- DESeq(dds_mac)
res_48 <- results(dds_mac, name = "Condition_Timepoint_SARSCov2.infected_48_vs_not.infected_48" , cooksCutoff = Inf)
res_48 <- lfcShrink(dds_mac, coef = "Condition_Timepoint_SARSCov2.infected_48_vs_not.infected_48" , res = res_48)
res_48 <- as.data.frame(res_48)
genes <- rownames(res_48)
res_48$genes<-rownames(res_48)
annot <- merge(x = genelist, y = res_48, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_48 <- annot

fc_interactions_int <- data.frame(int_2 = res_2$log2FoldChange[which(res_2$hgnc_symbol %in% interferon_revised)],int_8 = res_8$log2FoldChange[which(res_8$hgnc_symbol %in% interferon_revised)],int_24 = res_24$log2FoldChange[which(res_24$hgnc_symbol %in% interferon_revised)], int_48 = res_48$log2FoldChange[which(res_48$hgnc_symbol %in% interferon_revised)], row.names = res_2$hgnc_symbol[which(res_2$hgnc_symbol %in% interferon_revised)])
fc_interactions_int[is.na(fc_interactions_int)] <- 0

p_tp_int <- data.frame(res_2 = res_2$padj[which(res_2$hgnc_symbol %in% interferon_revised)], res_8 = res_8$padj[which(res_8$hgnc_symbol %in% interferon_revised)], res_24 = res_24$padj[which(res_24$hgnc_symbol %in% interferon_revised)], res_48 = res_48$padj[which(res_48$hgnc_symbol %in% interferon_revised)], row.names = res_8$hgnc_symbol[which(res_8$hgnc_symbol %in% interferon_revised)])
p_tp_int[p_tp_int < 0.05] <- "*"
p_tp_int[p_tp_int != "*"] <- ""
p_tp_int[is.na(p_tp_int)] <- ""

p_interactions_int[p_interactions_int < 0.05] <- "(*)"
p_interactions_int[p_interactions_int != "(*)"] <- ""
p_interactions_int[is.na(p_interactions_int)] <- ""
p_interactions_int <- cbind(rep("", nrow(p_interactions_int)), p_interactions_int)
for(i in 1:ncol(p_interactions_int)){
  p_tp_int <- cbind(p_tp_int, paste(p_tp_int[,i], p_interactions_int[,i]))
}
p_tp_int <- p_tp_int[,5:8]

paletteLength <- 50 
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(-2, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(2/paletteLength, 2, length.out=floor(paletteLength/2)))
pdf(file = 'Mac_time_Interferon.pdf', height = 6, width = 3)
pheatmap(fc_interactions_int, cluster_cols = F, breaks = myBreaks, color = myColor, display_numbers = p_tp_int)
dev.off()
hm <- pheatmap(fc_interactions_int, cluster_cols = F, breaks = myBreaks, color = myColor, display_numbers = p_tp_int)
order_int <- hm$tree_row$order


##Foam cells
counts_foam <- rawcounts[,which(sample_table$Cell == 'foam cells' & sample_table$Timepoint != 0)]
metadata_foam <- sample_table[which(sample_table$Cell == 'foam cells' & sample_table$Timepoint != 0),]
metadata_foam$Condition_Timepoint <- paste(metadata_foam$Condition,metadata_foam$Timepoint, sep = '_')
counts_foam <- counts_foam[-8686,]
counts_foam <- mutate_all(counts_foam, function(x) as.numeric(as.character(x)))

dds_foam<- DESeqDataSetFromMatrix(countData = counts_foam, colData = metadata_foam, design = ~Condition+Timepoint+Donor+Condition:Timepoint)
dds_foam<- DESeq(dds_foam)
resultsNames(dds_foam)
res <- results(dds_foam, name = "ConditionSARSCov2.infected.Timepoint8" , cooksCutoff = Inf)
res <- as.data.frame(res)
genes <- rownames(res)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
res$genes<-rownames(res)
annot <- merge(x = genelist, y = res, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_foam_int_8 <- annot
res <- results(dds_foam, name = "ConditionSARSCov2.infected.Timepoint24" , cooksCutoff = Inf)
res <- as.data.frame(res)
genes <- rownames(res)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
res$genes<-rownames(res)
annot <- merge(x = genelist, y = res, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_foam_int_24 <- annot
res <- results(dds_foam, name = "ConditionSARSCov2.infected.Timepoint48" , cooksCutoff = Inf)
res <- as.data.frame(res)
genes <- rownames(res)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
res$genes<-rownames(res)
annot <- merge(x = genelist, y = res, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_foam_int_48 <- annot


p_interactions_int <- data.frame(int_8 = res_foam_int_8$padj[which(res_foam_int_8$hgnc_symbol %in% interferon_revised)],int_24 = res_foam_int_24$padj[which(res_foam_int_24$hgnc_symbol %in% interferon_revised)],int_48 = res_foam_int_48$padj[which(res_foam_int_48$hgnc_symbol %in% interferon_revised)], row.names = res_foam_int_8$hgnc_symbol[which(res_foam_int_8$hgnc_symbol %in% interferon_revised)])

metadata_foam$Condition_Timepoint <- factor(metadata_foam$Condition_Timepoint)
dds_foam<- DESeqDataSetFromMatrix(countData = counts_foam, colData = metadata_foam, design = ~Condition_Timepoint+Donor)
dds_foam<- DESeq(dds_foam)
resultsNames(dds_foam)

res_2 <- results(dds_foam, name = "Condition_Timepoint_SARSCov2.infected_2_vs_not.infected_2" , cooksCutoff = Inf)
res_2 <- lfcShrink(dds_foam, coef = "Condition_Timepoint_SARSCov2.infected_2_vs_not.infected_2" , res = res_2)
res_2 <- as.data.frame(res_2)
genes <- rownames(res_2)
res_2$genes<-rownames(res_2)
annot <- merge(x = genelist, y = res_2, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_2 <- annot

metadata_foam$Condition_Timepoint <- relevel(metadata_foam$Condition_Timepoint, ref = 'not infected_8')
dds_foam<- DESeqDataSetFromMatrix(countData = counts_foam, colData = metadata_foam, design = ~Condition_Timepoint+Donor)
dds_foam<- DESeq(dds_foam)
res_8 <- results(dds_foam, name = "Condition_Timepoint_SARSCov2.infected_8_vs_not.infected_8" , cooksCutoff = Inf)
res_8 <- lfcShrink(dds_foam, coef = "Condition_Timepoint_SARSCov2.infected_8_vs_not.infected_8" , res = res_8)
res_8 <- as.data.frame(res_8)
genes <- rownames(res_8)
res_8$genes<-rownames(res_8)
annot <- merge(x = genelist, y = res_8, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_8 <- annot

metadata_foam$Condition_Timepoint <- relevel(metadata_foam$Condition_Timepoint, ref = 'not infected_24')
dds_foam<- DESeqDataSetFromMatrix(countData = counts_foam, colData = metadata_foam, design = ~Condition_Timepoint+Donor)
dds_foam<- DESeq(dds_foam)
res_24 <- results(dds_foam, name = "Condition_Timepoint_SARSCov2.infected_24_vs_not.infected_24" , cooksCutoff = Inf)
res_24 <- lfcShrink(dds_foam, coef = "Condition_Timepoint_SARSCov2.infected_24_vs_not.infected_24" , res = res_24)
res_24 <- as.data.frame(res_24)
genes <- rownames(res_24)
res_24$genes<-rownames(res_24)
annot <- merge(x = genelist, y = res_24, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_24 <- annot

metadata_foam$Condition_Timepoint <- relevel(metadata_foam$Condition_Timepoint, ref = 'not infected_48')
dds_foam<- DESeqDataSetFromMatrix(countData = counts_foam, colData = metadata_foam, design = ~Condition_Timepoint+Donor)
dds_foam<- DESeq(dds_foam)
res_48 <- results(dds_foam, name = "Condition_Timepoint_SARSCov2.infected_48_vs_not.infected_48" , cooksCutoff = Inf)
res_48 <- lfcShrink(dds_foam, coef = "Condition_Timepoint_SARSCov2.infected_48_vs_not.infected_48" , res = res_48)
res_48 <- as.data.frame(res_48)
genes <- rownames(res_48)
res_48$genes<-rownames(res_48)
annot <- merge(x = genelist, y = res_48, by.x="ensembl_gene_id", by.y="genes", all.y = T)
res_48 <- annot

fc_interactions_int <- data.frame(int_2 = res_2$log2FoldChange[which(res_2$hgnc_symbol %in% interferon_revised)],int_8 = res_8$log2FoldChange[which(res_8$hgnc_symbol %in% interferon_revised)],int_24 = res_24$log2FoldChange[which(res_24$hgnc_symbol %in% interferon_revised)], int_48 = res_48$log2FoldChange[which(res_48$hgnc_symbol %in% interferon_revised)], row.names = res_2$hgnc_symbol[which(res_2$hgnc_symbol %in% interferon_revised)])
fc_interactions_int[is.na(fc_interactions_int)] <- 0
fc_interactions_lysosomal <- data.frame(int_2 = res_2$log2FoldChange[which(res_2$hgnc_symbol %in% lysosomal)],int_8 = res_8$log2FoldChange[which(res_8$hgnc_symbol %in% lysosomal)],int_24 = res_24$log2FoldChange[which(res_24$hgnc_symbol %in% lysosomal)], int_48 = res_48$log2FoldChange[which(res_48$hgnc_symbol %in% lysosomal)], row.names = res_2$hgnc_symbol[which(res_2$hgnc_symbol %in% lysosomal)])
fc_interactions_lysosomal[is.na(fc_interactions_lysosomal)] <- 0


p_tp_int <- data.frame(res_2 = res_2$padj[which(res_2$hgnc_symbol %in% interferon_revised)], res_8 = res_8$padj[which(res_8$hgnc_symbol %in% interferon_revised)], res_24 = res_24$padj[which(res_24$hgnc_symbol %in% interferon_revised)], res_48 = res_48$padj[which(res_48$hgnc_symbol %in% interferon_revised)], row.names = res_8$hgnc_symbol[which(res_8$hgnc_symbol %in% interferon_revised)])
p_tp_int[p_tp_int < 0.05] <- "*"
p_tp_int[p_tp_int != "*"] <- ""
p_tp_int[is.na(p_tp_int)] <- ""

p_interactions_int[p_interactions_int < 0.05] <- "(*)"
p_interactions_int[p_interactions_int != "(*)"] <- ""
p_interactions_int[is.na(p_interactions_int)] <- ""
p_interactions_int <- cbind(rep("", nrow(p_interactions_int)), p_interactions_int)
for(i in 1:ncol(p_interactions_int)){
  p_tp_int <- cbind(p_tp_int, paste(p_tp_int[,i], p_interactions_int[,i]))
}
p_tp_int <- p_tp_int[,5:8]

paletteLength <- 50 
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(-2, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(2/paletteLength, 2, length.out=floor(paletteLength/2)))
pdf(file = 'foam_time_Interferon.pdf', height = 6, width = 3)
pheatmap(fc_interactions_int[order_int,], cluster_cols = F, breaks = myBreaks, color = myColor, display_numbers = p_tp_int[order_int,])
dev.off()


################################################################################################################################
##Figure 2J
### mac
counts_mac <- rawcounts[,which(sample_table$Cell == 'primary mac')]
metadata_mac <- sample_table[which(sample_table$Cell == 'primary mac'),]
metadata_mac$Condition_Timepoint <- paste(metadata_mac$Condition,metadata_mac$Timepoint, sep = '_')
counts_mac <- counts_mac[-8686,]
counts_mac <- mutate_all(counts_mac, function(x) as.numeric(as.character(x)))

dds_mac<- DESeqDataSetFromMatrix(countData = counts_mac, colData = metadata_mac, design = ~Condition_Timepoint+Donor)
dds_mac<- DESeq(dds_mac)
resultsNames(dds_mac)

res_2 <- results(dds_mac, name = "Condition_Timepoint_SARSCov2.infected_2_vs_not.infected_0", cooksCutoff = Inf)
res_2 <- lfcShrink(dds_mac, coef = "Condition_Timepoint_SARSCov2.infected_2_vs_not.infected_0", res = res_2)
res_2 <- as.data.frame(res_2)
res_2$Infected_2 <- res_2$log2FoldChange
res_8 <- results(dds_mac, name = "Condition_Timepoint_SARSCov2.infected_8_vs_not.infected_0", cooksCutoff = Inf)
res_8 <- lfcShrink(dds_mac, coef = "Condition_Timepoint_SARSCov2.infected_8_vs_not.infected_0", res = res_8)
res_8 <- as.data.frame(res_8)
res_2$Infected_8 <- res_8$log2FoldChange
res_24 <- results(dds_mac, name = "Condition_Timepoint_SARSCov2.infected_24_vs_not.infected_0", cooksCutoff = Inf)
res_24 <- lfcShrink(dds_mac, coef = "Condition_Timepoint_SARSCov2.infected_24_vs_not.infected_0", res = res_24)
res_24 <- as.data.frame(res_24)
res_2$Infected_24 <- res_24$log2FoldChange
res_48 <- results(dds_mac, name = "Condition_Timepoint_SARSCov2.infected_48_vs_not.infected_0", cooksCutoff = Inf)
res_48 <- lfcShrink(dds_mac, coef = "Condition_Timepoint_SARSCov2.infected_48_vs_not.infected_0", res = res_48)
res_48 <- as.data.frame(res_48)
res_2$Infected_48 <- res_48$log2FoldChange

res_2$Infected_2_adj <- res_2$padj
res_2$Infected_8_adj <- res_8$padj
res_2$Infected_24_adj <- res_24$padj
res_2$Infected_48_adj <- res_48$padj

res_2_ni <- results(dds_mac, name = "Condition_Timepoint_not.infected_2_vs_not.infected_0", cooksCutoff = Inf)
res_2_ni <- lfcShrink(dds_mac, coef = "Condition_Timepoint_not.infected_2_vs_not.infected_0", res = res_2_ni)
res_2_ni <- as.data.frame(res_2_ni)
res_2$Not_Infected_2 <- res_2_ni$log2FoldChange
res_8 <- results(dds_mac, name = "Condition_Timepoint_not.infected_8_vs_not.infected_0", cooksCutoff = Inf)
res_8 <- lfcShrink(dds_mac, coef = "Condition_Timepoint_not.infected_8_vs_not.infected_0", res = res_8)
res_8 <- as.data.frame(res_8)
res_2$Not_nfected_8 <- res_8$log2FoldChange
res_24 <- results(dds_mac, name = "Condition_Timepoint_not.infected_24_vs_not.infected_0", cooksCutoff = Inf)
res_24 <- lfcShrink(dds_mac, coef = "Condition_Timepoint_not.infected_24_vs_not.infected_0", res = res_24)
res_24 <- as.data.frame(res_24)
res_2$Not_Infected_24 <- res_24$log2FoldChange
res_48 <- results(dds_mac, name = "Condition_Timepoint_not.infected_48_vs_not.infected_0", cooksCutoff = Inf)
res_48 <- lfcShrink(dds_mac, coef = "Condition_Timepoint_not.infected_48_vs_not.infected_0", res = res_48)
res_48 <- as.data.frame(res_48)
res_2$Not_Infected_48 <- res_48$log2FoldChange

genes <- rownames(res_2)
res_2$genes<-rownames(res_2)
annot_mac <- merge(x = genelist, y = res_2, by.x="ensembl_gene_id", by.y="genes", all.y = T)

annot_sub <- annot_mac[which(annot_mac$hgnc_symbol %in% interferon_revised),]
mean_frame_mac_int <- data.frame(median = c(median(na.omit(annot_sub[,8])), median(na.omit(annot_sub[,9])),median(na.omit(annot_sub[,10])),
                                            median(na.omit(annot_sub[,11])),median(na.omit(annot_sub[,16])), median(na.omit(annot_sub[,17])),
                                            median(na.omit(annot_sub[,18])),median(na.omit(annot_sub[,19]))),
                                 Time = factor(rep(c('2', '8', '24', '48'),2), levels = c('2', '8', '24', '48')),
                                 Infection = c(rep('Infected', 4), rep('Not Infected', 4)),
                                 Cell = rep('Macrophage', 8))
mean_frame_mac_int$upper_IQR_upper <- c(boxplot.stats(annot_sub[,8])$stats[4], boxplot.stats(annot_sub[,9])$stats[4], boxplot.stats(annot_sub[,10])$stats[4], boxplot.stats(annot_sub[,11])$stats[4], boxplot.stats(annot_sub[,16])$stats[4], boxplot.stats(annot_sub[,17])$stats[4], boxplot.stats(annot_sub[,18])$stats[4], boxplot.stats(annot_sub[,19])$stats[4])
mean_frame_mac_int$upper_IQR_lower <- c(boxplot.stats(annot_sub[,8])$stats[2], boxplot.stats(annot_sub[,9])$stats[2], boxplot.stats(annot_sub[,10])$stats[2], boxplot.stats(annot_sub[,11])$stats[2],boxplot.stats(annot_sub[,16])$stats[2], boxplot.stats(annot_sub[,17])$stats[2], boxplot.stats(annot_sub[,18])$stats[2], boxplot.stats(annot_sub[,19])$stats[2])
gg_box_mac_int <- reshape2::melt(annot_sub[,8:11])
gg_box_mac_int$Cell <- rep('Macrophage', nrow(gg_box_mac_int))

annot_sub <- annot_mac[which(annot_mac$ensembl_gene_id %in% covid),]
mean_frame_mac_covid <- data.frame(median = c(median(na.omit(annot_sub[,8])), median(na.omit(annot_sub[,9])),median(na.omit(annot_sub[,10])),median(na.omit(annot_sub[,11])),median(na.omit(annot_sub[,16])), median(na.omit(annot_sub[,17])),median(na.omit(annot_sub[,18])),median(na.omit(annot_sub[,19]))),
                                   Time = factor(rep(c('2', '8', '24', '48'),2), levels = c('2', '8', '24', '48')),
                                   Infection = c(rep('Infected', 4), rep('Not Infected', 4)),
                                   Cell = rep('Macrophage', 8))
mean_frame_mac_covid$upper_IQR_upper <- c(boxplot.stats(annot_sub[,8])$stats[4], boxplot.stats(annot_sub[,9])$stats[4], boxplot.stats(annot_sub[,10])$stats[4], boxplot.stats(annot_sub[,11])$stats[4], boxplot.stats(annot_sub[,16])$stats[4], boxplot.stats(annot_sub[,17])$stats[4], boxplot.stats(annot_sub[,18])$stats[4], boxplot.stats(annot_sub[,19])$stats[4])
mean_frame_mac_covid$upper_IQR_lower <- c(boxplot.stats(annot_sub[,8])$stats[2], boxplot.stats(annot_sub[,9])$stats[2], boxplot.stats(annot_sub[,10])$stats[2], boxplot.stats(annot_sub[,11])$stats[2],boxplot.stats(annot_sub[,16])$stats[2], boxplot.stats(annot_sub[,17])$stats[2], boxplot.stats(annot_sub[,18])$stats[2], boxplot.stats(annot_sub[,19])$stats[2])


### foam
counts_foam <- rawcounts[,which(sample_table$Cell == 'foam cells')]
metadata_foam <- sample_table[which(sample_table$Cell == 'foam cells'),]
metadata_foam$Condition_Timepoint <- paste(metadata_foam$Condition,metadata_foam$Timepoint, sep = '_')
counts_foam <- counts_foam[-8686,]
counts_foam <- mutate_all(counts_foam, function(x) as.numeric(as.character(x)))

dds_foam<- DESeqDataSetFromMatrix(countData = counts_foam, colData = metadata_foam, design = ~Condition_Timepoint+Donor)
dds_foam<- DESeq(dds_foam)
resultsNames(dds_foam)

res_2 <- results(dds_foam, name = "Condition_Timepoint_SARSCov2.infected_2_vs_not.infected_0", cooksCutoff = Inf)
res_2 <- lfcShrink(dds_foam, coef = "Condition_Timepoint_SARSCov2.infected_2_vs_not.infected_0", res = res_2)
res_2 <- as.data.frame(res_2)
res_2$Infected_2 <- res_2$log2FoldChange
res_8 <- results(dds_foam, name = "Condition_Timepoint_SARSCov2.infected_8_vs_not.infected_0", cooksCutoff = Inf)
res_8 <- lfcShrink(dds_foam, coef = "Condition_Timepoint_SARSCov2.infected_8_vs_not.infected_0", res = res_8)
res_8 <- as.data.frame(res_8)
res_2$Infected_8 <- res_8$log2FoldChange
res_24 <- results(dds_foam, name = "Condition_Timepoint_SARSCov2.infected_24_vs_not.infected_0", cooksCutoff = Inf)
res_24 <- lfcShrink(dds_foam, coef = "Condition_Timepoint_SARSCov2.infected_24_vs_not.infected_0", res = res_24)
res_24 <- as.data.frame(res_24)
res_2$Infected_24 <- res_24$log2FoldChange
res_48 <- results(dds_foam, name = "Condition_Timepoint_SARSCov2.infected_48_vs_not.infected_0", cooksCutoff = Inf)
res_48 <- lfcShrink(dds_foam, coef = "Condition_Timepoint_SARSCov2.infected_48_vs_not.infected_0", res = res_48)
res_48 <- as.data.frame(res_48)
res_2$Infected_48 <- res_48$log2FoldChange

res_2$Infected_2_adj <- res_2$padj
res_2$Infected_8_adj <- res_8$padj
res_2$Infected_24_adj <- res_24$padj
res_2$Infected_48_adj <- res_48$padj

res_2_ni <- results(dds_foam, name = "Condition_Timepoint_not.infected_2_vs_not.infected_0", cooksCutoff = Inf)
res_2_ni <- lfcShrink(dds_foam, coef = "Condition_Timepoint_not.infected_2_vs_not.infected_0", res = res_2_ni)
res_2_ni <- as.data.frame(res_2_ni)
res_2$Not_Infected_2 <- res_2_ni$log2FoldChange
res_8 <- results(dds_foam, name = "Condition_Timepoint_not.infected_8_vs_not.infected_0", cooksCutoff = Inf)
res_8 <- lfcShrink(dds_foam, coef = "Condition_Timepoint_not.infected_8_vs_not.infected_0", res = res_8)
res_8 <- as.data.frame(res_8)
res_2$Not_nfected_8 <- res_8$log2FoldChange
res_24 <- results(dds_foam, name = "Condition_Timepoint_not.infected_24_vs_not.infected_0", cooksCutoff = Inf)
res_24 <- lfcShrink(dds_foam, coef = "Condition_Timepoint_not.infected_24_vs_not.infected_0", res = res_24)
res_24 <- as.data.frame(res_24)
res_2$Not_Infected_24 <- res_24$log2FoldChange
res_48 <- results(dds_foam, name = "Condition_Timepoint_not.infected_48_vs_not.infected_0", cooksCutoff = Inf)
res_48 <- lfcShrink(dds_foam, coef = "Condition_Timepoint_not.infected_48_vs_not.infected_0", res = res_48)
res_48 <- as.data.frame(res_48)
res_2$Not_Infected_48 <- res_48$log2FoldChange


genes <- rownames(res_2)
res_2$genes<-rownames(res_2)
annot_foam <- merge(x = genelist, y = res_2, by.x="ensembl_gene_id", by.y="genes", all.y = T)
annot_sub <- annot_foam[which(annot_foam$hgnc_symbol %in% interferon_revised),]
mean_frame_foam_int <- data.frame(median = c(median(na.omit(annot_sub[,8])), median(na.omit(annot_sub[,9])),median(na.omit(annot_sub[,10])),
                                             median(na.omit(annot_sub[,11])),median(na.omit(annot_sub[,16])), median(na.omit(annot_sub[,17])),
                                             median(na.omit(annot_sub[,18])),median(na.omit(annot_sub[,19]))),
                                  Time = factor(rep(c('2', '8', '24', '48'),2), levels = c('2', '8', '24', '48')),
                                  Infection = c(rep('Infected', 4), rep('Not Infected', 4)),
                                  Cell = rep('Foam', 8))
mean_frame_foam_int$upper_IQR_upper <- c(boxplot.stats(annot_sub[,8])$stats[4], boxplot.stats(annot_sub[,9])$stats[4], boxplot.stats(annot_sub[,10])$stats[4], boxplot.stats(annot_sub[,11])$stats[4], boxplot.stats(annot_sub[,16])$stats[4], boxplot.stats(annot_sub[,17])$stats[4], boxplot.stats(annot_sub[,18])$stats[4], boxplot.stats(annot_sub[,19])$stats[4])
mean_frame_foam_int$upper_IQR_lower <- c(boxplot.stats(annot_sub[,8])$stats[2], boxplot.stats(annot_sub[,9])$stats[2], boxplot.stats(annot_sub[,10])$stats[2], boxplot.stats(annot_sub[,11])$stats[2],boxplot.stats(annot_sub[,16])$stats[2], boxplot.stats(annot_sub[,17])$stats[2], boxplot.stats(annot_sub[,18])$stats[2], boxplot.stats(annot_sub[,19])$stats[2])
gg_box_foam_int <- reshape2::melt(annot_sub[,8:11])
gg_box_foam_int$Cell <- rep('Foam', nrow(gg_box_foam_int))

annot_sub <- annot_foam[which(annot_foam$ensembl_gene_id %in% covid),]
mean_frame_foam_covid <- data.frame(median = c(median(na.omit(annot_sub[,8])), median(na.omit(annot_sub[,9])),median(na.omit(annot_sub[,10])),median(na.omit(annot_sub[,11])),median(na.omit(annot_sub[,16])), median(na.omit(annot_sub[,17])),median(na.omit(annot_sub[,18])),median(na.omit(annot_sub[,19]))),
                                    Time = factor(rep(c('2', '8', '24', '48'),2), levels = c('2', '8', '24', '48')),
                                    Infection = c(rep('Infected', 4), rep('Not Infected', 4)),
                                    Cell = rep('Foam', 8))
mean_frame_foam_covid$upper_IQR_upper <- c(boxplot.stats(annot_sub[,8])$stats[4], boxplot.stats(annot_sub[,9])$stats[4], boxplot.stats(annot_sub[,10])$stats[4], boxplot.stats(annot_sub[,11])$stats[4], boxplot.stats(annot_sub[,16])$stats[4], boxplot.stats(annot_sub[,17])$stats[4], boxplot.stats(annot_sub[,18])$stats[4], boxplot.stats(annot_sub[,19])$stats[4])
mean_frame_foam_covid$upper_IQR_lower <- c(boxplot.stats(annot_sub[,8])$stats[2], boxplot.stats(annot_sub[,9])$stats[2], boxplot.stats(annot_sub[,10])$stats[2], boxplot.stats(annot_sub[,11])$stats[2],boxplot.stats(annot_sub[,16])$stats[2], boxplot.stats(annot_sub[,17])$stats[2], boxplot.stats(annot_sub[,18])$stats[2], boxplot.stats(annot_sub[,19])$stats[2])


### plots
###Type I interferon median of log2FC over baseline 

mean_frame_total_int <- rbind(mean_frame_mac_int, mean_frame_foam_int)
mean_frame_total_int$Cell_infection <- paste(mean_frame_total_int$Cell, mean_frame_total_int$Infection, sep = '_')
mean_frame_total_int$Time <- as.numeric(as.vector(mean_frame_total_int$Time))
pdf(file = 'infected_mac_vs_foam/median_log2FC_over_baseline_ylim_1.25_time_continuous_int.pdf', height = 4, width = 6)
ggplot(mean_frame_total_int[which(mean_frame_total_int$Infection == 'Infected'),], aes(y = median, x = Time, group = Cell_infection, color = Cell_infection))+
  geom_point(size = 2)+
  geom_line()+
  geom_errorbar(aes(ymin = upper_IQR_lower, ymax = upper_IQR_upper), width = 1, alpha = 0.5)+
  theme_bw()+
  ylab('median log2fc over baseline')+
  ylim(c(-.5,1.25))+
  ggtitle('Interferon I\nError bars: 25th and 75th quantiles')+
  scale_x_continuous(breaks=c(2,8, 24, 48))
dev.off()

###Covid genes median of log2FC over baseline 
mean_frame_total_covid <- rbind(mean_frame_mac_covid, mean_frame_foam_covid)
mean_frame_total_covid$Cell_infection <- paste(mean_frame_total_covid$Cell, mean_frame_total_covid$Infection, sep = '_')
mean_frame_total_covid$Time <- as.numeric(as.vector(mean_frame_total_covid$Time))

pdf(file = 'infected_mac_vs_foam/median_log2FC_over_baseline_time_continuous_covid.pdf', height = 4, width = 6)
ggplot(mean_frame_total_covid[which(mean_frame_total_covid$Infection == 'Infected'),], aes(y = median, x = Time, group = Cell_infection, color = Cell_infection))+
  geom_point(size = 2)+
  geom_line()+
  geom_errorbar(aes(ymin = upper_IQR_lower, ymax = upper_IQR_upper), width = 1, alpha = 0.5)+
  theme_bw()+
  ylab('median log2fc over baseline')+
  ggtitle('COVID\nError bars: 25th and 75th quantiles')+
  scale_x_continuous(breaks=c(2,8, 24, 48))
dev.off()


################################################################################################################
## Extend. Fig 3G&H

#Compare mac vs foam cells at infected timepoint
rawcounts_infected<-read.csv("~/DESeq/counts_infected.csv",header=T)
sample_table_infected<-read.csv("~/DESeq/metadata_infected.csv",header=T)

sample_table_infected$Donor<-as.factor(sample_table_infected$Donor)
sample_table_infected$Cell<-as.factor(sample_table_infected$Cell)
sample_table_infected$Timepoint<-as.factor(sample_table_infected$Timepoint)
sample_table_infected$cell_timepoint <- as.factor(paste0(sample_table_infected$Cell, "_", sample_table_infected$Timepoint))

dds <- DESeqDataSetFromMatrix(countData = rawcounts_infected,
                              colData = sample_table_infected,
                              design = ~Donor+cell_timepoint)

dds <-DESeq(dds)
resultsNames(dds)

time0 <- results(dds,name=c("cell_timepoint_primary.mac_0_vs_foam.cells_0"))
time0 <- as.data.frame(lfcShrink(dds, coef = 7, res = time0, type='apeglm'))
time0  <- time0[order(time0$padj),]
sum(time0$padj < 0.05, na.rm=TRUE)
time0  <- time0[order(time0$pvalue),]
sum(time0$pvalue < 0.05, na.rm=TRUE)
time0 <- na.omit(time0) 

dds$cell_timepoint<-relevel(dds$cell_timepoint,"foam cells_2")
dds <-DESeq(dds)
resultsNames(dds)
time2 <- results(dds,name=c("cell_timepoint_primary.mac_2_vs_foam.cells_2"))
time2 <- lfcShrink(dds, coef=8,res = time2, type='apeglm')
time2  <- time2[order(time2$padj),]
sum(time2$padj < 0.05, na.rm=TRUE)
time2  <- time2[order(time2$pvalue),]
sum(time2$pvalue < 0.05, na.rm=TRUE)
time2 <- na.omit(time2) 

dds$cell_timepoint<-relevel(dds$cell_timepoint,"foam cells_8")
dds <-DESeq(dds)
resultsNames(dds)
time8 <- results(dds,name=c("cell_timepoint_primary.mac_8_vs_foam.cells_8"))
time8 <- lfcShrink(dds, coef=11,res = time8, type='apeglm')
time8  <- time8[order(time8$padj),]
sum(time8$padj < 0.05, na.rm=TRUE)
time8  <- time8[order(time8$pvalue),]
#time8 
sum(time8$pvalue < 0.05, na.rm=TRUE)
time8 <- na.omit(time8) 

dds$cell_timepoint<-relevel(dds$cell_timepoint,"foam cells_24")
dds <-DESeq(dds)
resultsNames(dds)
time24 <- results(dds,name=c("cell_timepoint_primary.mac_24_vs_foam.cells_24"))
time24 <- lfcShrink(dds, coef=9,res = time24, type='apeglm')
time24  <- time24[order(time24$padj),]
#time24
sum(time24$padj < 0.05, na.rm=TRUE)
time24  <- time24[order(time24$pvalue),]
sum(time24$pvalue < 0.05, na.rm=TRUE)
time24 <- na.omit(time24) 

dds$cell_timepoint<-relevel(dds$cell_timepoint,"foam cells_48")
dds <-DESeq(dds)
resultsNames(dds)
time48 <- results(dds,name=c("cell_timepoint_primary.mac_48_vs_foam.cells_48"))
time48 <- lfcShrink(dds, coef=10,res = time48, type='apeglm')
time48  <- time48[order(time48$padj),]
sum(time48$padj < 0.05, na.rm=TRUE)
time48  <- time48[order(time48$pvalue),]
sum(time48$pvalue < 0.05, na.rm=TRUE)
time48 <- na.omit(time48) 

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
time0<-as.data.frame(time0)
genes <- rownames(time0)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
time0$genes<-rownames(time0)
annot <- merge(
  x = genelist,
  y = time0,
  by.x="ensembl_gene_id",
  by.y="genes")
time0<-annot

time2<-as.data.frame(time2)
genes <- rownames(time2)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
time2$genes<-rownames(time2)
annot <- merge(
  x = genelist,
  y = time2,
  by.x="ensembl_gene_id",
  by.y="genes")
time2<-annot

time8<-as.data.frame(time8)
genes <- rownames(time8)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
time8$genes<-rownames(time8)
annot <- merge(
  x = genelist,
  y = time8,
  by.x="ensembl_gene_id",
  by.y="genes")
time8<-annot

time24<-as.data.frame(time24)
genes <- rownames(time24)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
time24$genes<-rownames(time24)
annot <- merge(
  x = genelist,
  y = time24,
  by.x="ensembl_gene_id",
  by.y="genes")
time24<-annot

time48<-as.data.frame(time48)
genes <- rownames(time48)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
time48$genes<-rownames(time48)
annot <- merge(
  x = genelist,
  y = time48,
  by.x="ensembl_gene_id",
  by.y="genes")
time48<-annot

lysos<-c("CTSA","CTSB","CTSD","CTSF","GLA","GLB1","GNS","HEXA","ARSA","ARSB","LAMP1","LAMP2")

colnames(time0)[2]<-"geneid"
time0<-time0[,-c(1)]
names <- make.unique(time0$geneid)
rownames(time0) <- names
time0 <- time0[,-1] 

colnames(time2)[2]<-"geneid"
time2<-time2[,-c(1)]
names <- make.unique(time2$geneid)
rownames(time2) <- names
time2 <- time2[,-1] 

colnames(time8)[2]<-"geneid"
time8<-time8[,-c(1)]
names <- make.unique(time8$geneid)
rownames(time8) <- names
time8 <- time8[,-1] 

colnames(time24)[2]<-"geneid"
time24<-time24[,-c(1)]
names <- make.unique(time24$geneid)
rownames(time24) <- names
time24 <- time24[,-1] 

colnames(time48)[2]<-"geneid"
time48<-time48[,-c(1)]
names <- make.unique(time48$geneid)
rownames(time48) <- names
time48 <- time48[,-1] 

padj_vals_0 = as.data.frame(time0[lysos, "pvalue"])
colnames(padj_vals_0) <- c('Cell_primary.mac_vs_foam.cells at Timepoint 0')
sig_genes_lysos<- time0[lysos,]
rownames(padj_vals_0) <- rownames(sig_genes_lysos)

padj_vals_2 = as.data.frame(time2[lysos, "pvalue"])
colnames(padj_vals_2) <- c('Cell_primary.mac_vs_foam.cells at Timepoint 2')
sig_genes_lysos<- time2[lysos,]
rownames(padj_vals_2) <- rownames(sig_genes_lysos)

padj_vals_8 = as.data.frame(time8[lysos, "pvalue"])
colnames(padj_vals_8) <- c('Cell_primary.mac_vs_foam.cells at Timepoint 8')
sig_genes_lysos<- time8[lysos,]
rownames(padj_vals_8) <- rownames(sig_genes_lysos)

padj_vals_24 = as.data.frame(time24[lysos, "pvalue"])
colnames(padj_vals_24) <- c('Cell_primary.mac_vs_foam.cells at Timepoint 24')
sig_genes_lysos<- time24[lysos,]
rownames(padj_vals_24) <- rownames(sig_genes_lysos)

padj_vals_48 = as.data.frame(time48[lysos, "pvalue"])
colnames(padj_vals_48) <- c('Cell_primary.mac_vs_foam.cells at Timepoint 48')
sig_genes_lysos<- time48[lysos,]
rownames(padj_vals_48) <- rownames(sig_genes_lysos)

pvals_BH<-data.frame(padj_vals_0$`Cell_primary.mac_vs_foam.cells at Timepoint 0`,
                     padj_vals_2$`Cell_primary.mac_vs_foam.cells at Timepoint 2`,
                     padj_vals_8$`Cell_primary.mac_vs_foam.cells at Timepoint 8`,
                     padj_vals_24$`Cell_primary.mac_vs_foam.cells at Timepoint 24`,
                     padj_vals_48$`Cell_primary.mac_vs_foam.cells at Timepoint 48`)
rownames(pvals_BH)<-rownames(sig_genes_lysos)
pvals_BH[is.na(pvals_BH)] = 0.0
pvals_BH <- pvals_BH[grepl("^NA", rownames(pvals_BH))==F,]

log_fc_0 = as.data.frame(time0[lysos, "log2FoldChange"])
colnames(log_fc_0) <- c('Primary mac vs foam cells at Timepoint 0')
sig_genes_lysos<- time0[lysos,]
rownames(log_fc_0) <- rownames(sig_genes_lysos)

log_fc_2 = as.data.frame(time2[lysos, "log2FoldChange"])
colnames(log_fc_2) <- c('Primary mac vs foam cells at Timepoint 2')
sig_genes_lysos<- time2[lysos,]
rownames(log_fc_2) <- rownames(sig_genes_lysos)

log_fc_8 = as.data.frame(time8[lysos, "log2FoldChange"])
colnames(log_fc_8) <- c('Primary mac vs foam cells at Timepoint 8')
sig_genes_lysos<- time8[lysos,]
rownames(log_fc_8) <- rownames(sig_genes_lysos)

log_fc_24 = as.data.frame(time24[lysos, "log2FoldChange"])
colnames(log_fc_24) <- c('Primary mac vs foam cells at Timepoint 24')
sig_genes_lysos<- time24[lysos,]
rownames(log_fc_24) <- rownames(sig_genes_lysos)

log_fc_48 = as.data.frame(time48[lysos, "log2FoldChange"])
colnames(log_fc_48) <- c('Primary mac vs foam cells at Timepoint 48')
sig_genes_lysos<- time48[lysos,]
rownames(log_fc_48) <- rownames(sig_genes_lysos)

foldchange_lysos<-data.frame(log_fc_0$`Primary mac vs foam cells at Timepoint 0`,
                              log_fc_2$`Primary mac vs foam cells at Timepoint 2`,
                              log_fc_8$`Primary mac vs foam cells at Timepoint 8`,
                              log_fc_24$`Primary mac vs foam cells at Timepoint 24`,
                              log_fc_48$`Primary mac vs foam cells at Timepoint 48`)
rownames(foldchange_lysos)<-rownames(sig_genes_lysos)

colnames(foldchange_lysos)<-c("Primary mac vs foam cells at Timepoint 0","Primary mac vs foam cells at Timepoint 2","Primary mac vs foam cells at Timepoint 8","Primary mac vs foam cells at Timepoint 24",
                               "Primary mac vs foam cells at Timepoint 48")
foldchange_lysos[is.na(foldchange_lysos)] = 0.0
foldchange_lysos[foldchange_lysos > 10] = 10
foldchange_lysos[foldchange_lysos < -10] = -10
foldchange_lysos <- foldchange_lysos[grepl("^NA", rownames(foldchange_lysos))==F,]

foldchange_lysos = data.matrix(foldchange_lysos)
noted = matrix("", nrow=nrow(foldchange_lysos), ncol=ncol(foldchange_lysos))
noted[pvals_BH < 0.1] = " "
noted[pvals_BH < 0.05] = "*"
noted[pvals_BH < 0.01] = "**"
noted[pvals_BH < 0.001] = "***"

pdf("DESeq/Mac_vs_foam_heatmap_lysosomal_genes_nominal_pvalue_shrunken.pdf", width=10,height=10)
heatmap.2(foldchange_lysos,
          Colv=FALSE,
          Rowv=FALSE,
          mar=c(6, 14),
          margins=c(20,8),
          trace="none",
          breaks=seq(-1,1, length.out=101),
          col=colorRampPalette(c("blue", "white", "red"))(100),
          key.xlab=expression("Primary Mac vs foam cells notinfected"("log"[2] * " FC")), key.ylab=NA,
          key.title="", 
          cellnote=noted,
          notecol="black",
          cexCol=1.0,
          cexRow=1.0,
          notecex=1.5,
          dendrogram="none"
)
dev.off()

# Differential expression analysis of Infected vs Not infected in Primary Mac and Foam cells at each Timepoint
sample_table_infected<-sample_table[which(sample_table$Cell == "primary mac"),]
sample_table_primarymac<-sample_table_infected[which(sample_table_infected$Timepoint != "0"),]
rawcounts_primarymac<-rawcounts[,-c(9:10,15:30,34:36,40:42,46:48,52:54)]

sample_table_primarymac$Donor<-as.factor(sample_table_primarymac$Donor)
sample_table_primarymac$Cell<-as.factor(sample_table_primarymac$Cell)
sample_table_primarymac$Condition<-as.factor(sample_table_primarymac$Condition)
sample_table_primarymac$Timepoint<-as.factor(sample_table_primarymac$Timepoint)
sample_table_primarymac$condition_time <- as.factor(paste0(sample_table_primarymac$Condition, "_", sample_table_primarymac$Timepoint))

dds_primarymac <- DESeqDataSetFromMatrix(countData = rawcounts_primarymac,
                                         colData = sample_table_primarymac,
                                         design = ~Donor+condition_time)

dds_primarymac <-DESeq(dds_primarymac)
resultsNames(dds_primarymac)

time2 <- results(dds_primarymac,name=c("condition_time_SARSCov2.infected_2_vs_not.infected_2"))
time2 <- as.data.frame(lfcShrink(dds_primarymac, coef = 6, res = time2, type='apeglm'))
time2  <- time2[order(time2$padj),]
sum(time2$padj < 0.05, na.rm=TRUE)
time2  <- time2[order(time2$pvalue),]
sum(time2$pvalue < 0.05, na.rm=TRUE)
time2 <- na.omit(time2) 

dds_primarymac$condition_time<-relevel(dds_primarymac$condition_time,"not infected_8")
dds_primarymac <-DESeq(dds_primarymac)
resultsNames(dds_primarymac)
time8 <- results(dds_primarymac,name=c("condition_time_SARSCov2.infected_8_vs_not.infected_8"))
time8 <- as.data.frame(lfcShrink(dds_primarymac, coef=9, res = time8, type='apeglm'))
time8  <- time8[order(time8$padj),]
sum(time8$padj < 0.05, na.rm=TRUE)
time8  <- time8[order(time8$pvalue),]
sum(time8$pvalue < 0.05, na.rm=TRUE)
time8 <- na.omit(time8) 

dds_primarymac$condition_time<-relevel(dds_primarymac$condition_time,"not infected_24")
dds_primarymac <-DESeq(dds_primarymac)
resultsNames(dds_primarymac)
time24 <- results(dds_primarymac,name=c("condition_time_SARSCov2.infected_24_vs_not.infected_24"))
time24 <- as.data.frame(lfcShrink(dds_primarymac, coef=7, res = time24, type='apeglm'))
time24  <- time24[order(time24$padj),]
sum(time24$padj < 0.05, na.rm=TRUE)
time24  <- time24[order(time24$pvalue),]
sum(time24$pvalue < 0.05, na.rm=TRUE)
time24 <- na.omit(time24) 

dds_primarymac$condition_time<-relevel(dds_primarymac$condition_time,"not infected_48")
dds_primarymac <-DESeq(dds_primarymac)
resultsNames(dds_primarymac)
time48 <- results(dds_primarymac,name=c("condition_time_SARSCov2.infected_48_vs_not.infected_48"))
time48 <- as.data.frame(lfcShrink(dds_primarymac, coef=8, res = time48, type='apeglm'))
time48  <- time48[order(time48$padj),]
sum(time48$padj < 0.05, na.rm=TRUE)
time48  <- time48[order(time48$pvalue),]
sum(time48$pvalue < 0.05, na.rm=TRUE)
time48 <- na.omit(time48) 

time2<-as.data.frame(time2)
genes <- rownames(time2)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
time2$genes<-rownames(time2)
annot <- merge(
  x = genelist,
  y = time2,
  by.x="ensembl_gene_id",
  by.y="genes")
time2<-annot

time8<-as.data.frame(time8)
genes <- rownames(time8)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
time8$genes<-rownames(time8)
annot <- merge(
  x = genelist,
  y = time8,
  by.x="ensembl_gene_id",
  by.y="genes")
time8<-annot

time24<-as.data.frame(time24)
genes <- rownames(time24)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
time24$genes<-rownames(time24)
annot <- merge(
  x = genelist,
  y = time24,
  by.x="ensembl_gene_id",
  by.y="genes")
time24<-annot

time48<-as.data.frame(time48)
genes <- rownames(time48)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
time48$genes<-rownames(time48)
annot <- merge(
  x = genelist,
  y = time48,
  by.x="ensembl_gene_id",
  by.y="genes")
time48<-annot

colnames(time2)[2]<-"geneid"
time2<-time2[,-c(1)]
names <- make.unique(time2$geneid)
rownames(time2) <- names
time2 <- time2[,-1] 

colnames(time8)[2]<-"geneid"
time8<-time8[,-c(1)]
names <- make.unique(time8$geneid)
rownames(time8) <- names
time8 <- time8[,-1] 

colnames(time24)[2]<-"geneid"
time24<-time24[,-c(1)]
names <- make.unique(time24$geneid)
rownames(time24) <- names
time24 <- time24[,-1] 

colnames(time48)[2]<-"geneid"
time48<-time48[,-c(1)]
names <- make.unique(time48$geneid)
rownames(time48) <- names
time48 <- time48[,-1] 

padj_vals_2 = as.data.frame(time2[lysos, "padj"])
colnames(padj_vals_2) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 2')
sig_genes_lysos<- time2[lysos,]
rownames(padj_vals_2) <- rownames(sig_genes_lysos)

padj_vals_8 = as.data.frame(time8[lysos, "padj"])
colnames(padj_vals_8) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 8')
sig_genes_lysos<- time8[lysos,]
rownames(padj_vals_8) <- rownames(sig_genes_lysos)

padj_vals_24 = as.data.frame(time24[lysos, "padj"])
colnames(padj_vals_24) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 24')
sig_genes_lysos<- time24[lysos,]
rownames(padj_vals_24) <- rownames(sig_genes_lysos)

padj_vals_48 = as.data.frame(time48[lysos, "padj"])
colnames(padj_vals_48) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 48')
sig_genes_lysos<- time48[lysos,]
rownames(padj_vals_48) <- rownames(sig_genes_lysos)

pvals_BH<-data.frame(padj_vals_2$`SARSCov2 Infected_vs_notinfected at Timepoint 2`,
                     padj_vals_8$`SARSCov2 Infected_vs_notinfected at Timepoint 8`,
                     padj_vals_24$`SARSCov2 Infected_vs_notinfected at Timepoint 24`,
                     padj_vals_48$`SARSCov2 Infected_vs_notinfected at Timepoint 48`)
rownames(pvals_BH)<-rownames(sig_genes_lysos)
pvals_BH[is.na(pvals_BH)] = 0.0
pvals_BH <- pvals_BH[grepl("^NA", rownames(pvals_BH))==F,]

log_fc_2 = as.data.frame(time2[lysos, "log2FoldChange"])
colnames(log_fc_2) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 2')
sig_genes_lysos<- time2[lysos,]
rownames(log_fc_2) <- rownames(sig_genes_lysos)

log_fc_8 = as.data.frame(time8[lysos, "log2FoldChange"])
colnames(log_fc_8) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 8')
sig_genes_lysos<- time8[lysos,]
rownames(log_fc_8) <- rownames(sig_genes_lysos)

log_fc_24 = as.data.frame(time24[lysos, "log2FoldChange"])
colnames(log_fc_24) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 24')
sig_genes_lysos<- time24[lysos,]
rownames(log_fc_24) <- rownames(sig_genes_lysos)

log_fc_48 = as.data.frame(time48[lysos, "log2FoldChange"])
colnames(log_fc_48) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 48')
sig_genes_lysos<- time48[lysos,]
rownames(log_fc_48) <- rownames(sig_genes_lysos)

foldchange_lysos<-data.frame(log_fc_2$`SARSCov2 Infected_vs_notinfected at Timepoint 2`,
                              log_fc_8$`SARSCov2 Infected_vs_notinfected at Timepoint 8`,
                              log_fc_24$`SARSCov2 Infected_vs_notinfected at Timepoint 24`,
                              log_fc_48$`SARSCov2 Infected_vs_notinfected at Timepoint 48`)
rownames(foldchange_lysos)<-rownames(sig_genes_lysos)

colnames(foldchange_lysos)<-c("SARSCov2 Infected_vs_notinfected at Timepoint 2","SARSCov2 Infected_vs_notinfected at Timepoint 8",
                               "SARSCov2 Infected_vs_notinfected at Timepoint 24","SARSCov2 Infected_vs_notinfected at Timepoint 48")
foldchange_lysos[is.na(foldchange_lysos)] = 0.0
foldchange_lysos[foldchange_lysos > 10] = 10
foldchange_lysos[foldchange_lysos < -10] = -10
foldchange_lysos <- foldchange_lysos[grepl("^NA", rownames(foldchange_lysos))==F,]

foldchange_lysos = data.matrix(foldchange_lysos)
noted = matrix("", nrow=nrow(foldchange_lysos), ncol=ncol(foldchange_lysos))
noted[pvals_BH < 0.1] = " "
noted[pvals_BH < 0.05] = "*"
noted[pvals_BH < 0.01] = "**"
noted[pvals_BH < 0.001] = "***"


pdf("DESeq/Mac_foldchange_heatmap_lysosomal_genes_shrunken_v1.pdf", width=10,height=10)
heatmap.2(foldchange_lysos,
          Colv=FALSE,
          Rowv=FALSE,
          mar=c(6, 14),
          margins=c(14,34),
          trace="none",
          # col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
          breaks=seq(-0.3,0.3, length.out=101),
          col=colorRampPalette(c("blue", "white", "red"))(100),
          key.xlab=expression("Infected vs Not infected/Primary Mac"("log"[2] * " FC")), key.ylab=NA,
          key.title="", 
          cellnote=noted,
          notecol="black",
          cexCol=1.0,
          cexRow=1.0,
          notecex=1.5,
          dendrogram="none"
)
dev.off()


sample_table_notinfected<-sample_table[which(sample_table$Cell == "foam cells"),]
sample_table_foam<-sample_table_notinfected[which(sample_table_notinfected$Timepoint != "0"),]
rawcounts_foam<-rawcounts[,-c(1:15,24:25,30:33,37:39,43:45,49:51)]

sample_table_foam$Donor<-as.factor(sample_table_foam$Donor)
sample_table_foam$Cell<-as.factor(sample_table_foam$Cell)
sample_table_foam$Condition<-as.factor(sample_table_foam$Condition)
sample_table_foam$Timepoint<-as.factor(sample_table_foam$Timepoint)
sample_table_foam$condition_timepoint <- as.factor(paste0(sample_table_foam$Condition, "_", sample_table_foam$Timepoint))
sample_table_foam$condition_timepoint<-as.factor(sample_table_foam$condition_timepoint)

dds_foam <- DESeqDataSetFromMatrix(countData = rawcounts_foam,
                                   colData = sample_table_foam,
                                   design = ~Donor+condition_timepoint)

dds_foam <-DESeq(dds_foam)
resultsNames(dds_foam)

time2 <- results(dds_foam,name=c("condition_timepoint_SARSCov2.infected_2_vs_not.infected_2"))
time2 <- as.data.frame(lfcShrink(dds_foam, coef=6, res = time2, type='apeglm'))
time2  <- time2[order(time2$padj),]
sum(time2$padj < 0.05, na.rm=TRUE)
time2  <- time2[order(time2$pvalue),]
sum(time2$pvalue < 0.05, na.rm=TRUE)
time2 <- na.omit(time2) 

dds_foam$condition_timepoint<-relevel(dds_foam$condition_timepoint,"not infected_8")
dds_foam <-DESeq(dds_foam)
resultsNames(dds_foam)
time8 <- results(dds_foam,name=c("condition_timepoint_SARSCov2.infected_8_vs_not.infected_8"))
time8 <- as.data.frame(lfcShrink(dds_foam, coef=9, res = time8, type='apeglm'))
time8  <- time8[order(time8$padj),]
sum(time8$padj < 0.05, na.rm=TRUE)
time8  <- time8[order(time8$pvalue),]
sum(time8$pvalue < 0.05, na.rm=TRUE)
time8 <- na.omit(time8) 

dds_foam$condition_timepoint<-relevel(dds_foam$condition_timepoint,"not infected_24")
dds_foam <-DESeq(dds_foam)
resultsNames(dds_foam)
time24 <- results(dds_foam,name=c("condition_timepoint_SARSCov2.infected_24_vs_not.infected_24"))
time24 <- as.data.frame(lfcShrink(dds_foam, coef=7, res = time24, type='apeglm'))
time24  <- time24[order(time24$padj),]
sum(time24$padj < 0.05, na.rm=TRUE)
time24  <- time24[order(time24$pvalue),]
sum(time24$pvalue < 0.05, na.rm=TRUE)
time24 <- na.omit(time24) 

dds_foam$condition_timepoint<-relevel(dds_foam$condition_timepoint,"not infected_48")
dds_foam <-DESeq(dds_foam)
resultsNames(dds_foam)
time48 <- results(dds_foam,name=c("condition_timepoint_SARSCov2.infected_48_vs_not.infected_48"))
time48 <- as.data.frame(lfcShrink(dds_foam, coef=8, res = time48, type='apeglm'))
time48  <- time48[order(time48$padj),]
sum(time48$padj < 0.05, na.rm=TRUE)
time48  <- time48[order(time48$pvalue),]
sum(time48$pvalue < 0.05, na.rm=TRUE)
time48 <- na.omit(time48) 

time2<-as.data.frame(time2)
genes <- rownames(time2)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
time2$genes<-rownames(time2)
annot <- merge(
  x = genelist,
  y = time2,
  by.x="ensembl_gene_id",
  by.y="genes")
time2<-annot

time8<-as.data.frame(time8)
genes <- rownames(time8)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
time8$genes<-rownames(time8)
annot <- merge(
  x = genelist,
  y = time8,
  by.x="ensembl_gene_id",
  by.y="genes")
time8<-annot

time24<-as.data.frame(time24)
genes <- rownames(time24)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
time24$genes<-rownames(time24)
annot <- merge(
  x = genelist,
  y = time24,
  by.x="ensembl_gene_id",
  by.y="genes")
time24<-annot

time48<-as.data.frame(time48)
genes <- rownames(time48)
genelist <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
time48$genes<-rownames(time48)
annot <- merge(
  x = genelist,
  y = time48,
  by.x="ensembl_gene_id",
  by.y="genes")
time48<-annot

colnames(time2)[2]<-"geneid"
time2<-time2[,-c(1)]
names <- make.unique(time2$geneid)
rownames(time2) <- names
time2 <- time2[,-1] 

colnames(time8)[2]<-"geneid"
time8<-time8[,-c(1)]
names <- make.unique(time8$geneid)
rownames(time8) <- names
time8 <- time8[,-1] 

colnames(time24)[2]<-"geneid"
time24<-time24[,-c(1)]
names <- make.unique(time24$geneid)
rownames(time24) <- names
time24 <- time24[,-1] 

colnames(time48)[2]<-"geneid"
time48<-time48[,-c(1)]
names <- make.unique(time48$geneid)
rownames(time48) <- names
time48 <- time48[,-1] 

padj_vals_2 = as.data.frame(time2[lysos, "padj"])
colnames(padj_vals_2) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 2')
sig_genes_lysos<- time2[lysos,]
rownames(padj_vals_2) <- rownames(sig_genes_lysos)

padj_vals_8 = as.data.frame(time8[lysos, "padj"])
colnames(padj_vals_8) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 8')
sig_genes_lysos<- time8[lysos,]
rownames(padj_vals_8) <- rownames(sig_genes_lysos)

padj_vals_24 = as.data.frame(time24[lysos, "padj"])
colnames(padj_vals_24) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 24')
sig_genes_lysos<- time24[lysos,]
rownames(padj_vals_24) <- rownames(sig_genes_lysos)

padj_vals_48 = as.data.frame(time48[lysos, "padj"])
colnames(padj_vals_48) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 48')
sig_genes_lysos<- time48[lysos,]
rownames(padj_vals_48) <- rownames(sig_genes_lysos)

pvals_BH<-data.frame(padj_vals_2$`SARSCov2 Infected_vs_notinfected at Timepoint 2`,
                     padj_vals_8$`SARSCov2 Infected_vs_notinfected at Timepoint 8`,
                     padj_vals_24$`SARSCov2 Infected_vs_notinfected at Timepoint 24`,
                     padj_vals_48$`SARSCov2 Infected_vs_notinfected at Timepoint 48`)
rownames(pvals_BH)<-rownames(sig_genes_lysos)
pvals_BH[is.na(pvals_BH)] = 0.0
pvals_BH <- pvals_BH[grepl("^NA", rownames(pvals_BH))==F,]

log_fc_2 = as.data.frame(time2[lysos, "log2FoldChange"])
colnames(log_fc_2) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 2')
sig_genes_lysos<- time2[lysos,]
rownames(log_fc_2) <- rownames(sig_genes_lysos)

log_fc_8 = as.data.frame(time8[lysos, "log2FoldChange"])
colnames(log_fc_8) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 8')
sig_genes_lysos<- time8[lysos,]
rownames(log_fc_8) <- rownames(sig_genes_lysos)

log_fc_24 = as.data.frame(time24[lysos, "log2FoldChange"])
colnames(log_fc_24) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 24')
sig_genes_lysos<- time24[lysos,]
rownames(log_fc_24) <- rownames(sig_genes_lysos)

log_fc_48 = as.data.frame(time48[lysos, "log2FoldChange"])
colnames(log_fc_48) <- c('SARSCov2 Infected_vs_notinfected at Timepoint 48')
sig_genes_lysos<- time48[lysos,]
rownames(log_fc_48) <- rownames(sig_genes_lysos)

foldchange_lysos<-data.frame(log_fc_2$`SARSCov2 Infected_vs_notinfected at Timepoint 2`,
                              log_fc_8$`SARSCov2 Infected_vs_notinfected at Timepoint 8`,
                              log_fc_24$`SARSCov2 Infected_vs_notinfected at Timepoint 24`,
                              log_fc_48$`SARSCov2 Infected_vs_notinfected at Timepoint 48`)
rownames(foldchange_lysos)<-rownames(sig_genes_lysos)

colnames(foldchange_lysos)<-c("SARSCov2 Infected_vs_notinfected at Timepoint 2","SARSCov2 Infected_vs_notinfected at Timepoint 8",
                               "SARSCov2 Infected_vs_notinfected at Timepoint 24","SARSCov2 Infected_vs_notinfected at Timepoint 48")
foldchange_lysos[is.na(foldchange_lysos)] = 0.0
foldchange_lysos[foldchange_lysos > 10] = 10
foldchange_lysos[foldchange_lysos < -10] = -10
foldchange_lysos <- foldchange_lysos[grepl("^NA", rownames(foldchange_lysos))==F,]

foldchange_lysos = data.matrix(foldchange_lysos)
noted = matrix("", nrow=nrow(foldchange_lysos), ncol=ncol(foldchange_lysos))
noted[pvals_BH < 0.1] = " "
noted[pvals_BH < 0.05] = "*"
noted[pvals_BH < 0.01] = "**"
noted[pvals_BH < 0.001] = "***"


pdf("DESeq/Foam_foldchange_heatmap_lysosomal_genes_shrunken_v1.pdf", width=10,height=10)
heatmap.2(foldchange_lysos,
          Colv=FALSE,
          Rowv=FALSE,
          mar=c(6, 14),
          margins=c(14,34),
          trace="none",
          # col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
          breaks=seq(-0.3,0.3, length.out=101),
          col=colorRampPalette(c("blue", "white", "red"))(100),
          key.xlab=expression("Infected vs Not infected/Foam"("log"[2] * " FC")), key.ylab=NA,
          key.title="", 
          cellnote=noted,
          notecol="black",
          cexCol=1.0,
          cexRow=1.0,
          notecex=1.5,
          dendrogram="none"
)
dev.off()



