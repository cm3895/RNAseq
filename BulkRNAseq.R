library("DESeq2")
library("ggplot2")
library("ggrepel")

#### Load and Clean Data, Only First Time- Counts #######

# Read in Count matrix, make gene names unique
Counts <- read.table('est_counts_genes_kallisto.txt')
uninames <- make.names(rownames(Counts),unique=TRUE)
row.names(Counts) <- uninames
Info <- read.csv('Feb_bulkRNAseq_GroupV.csv',row.names='Sample')
dds <- DESeqDataSetFromMatrix(countData = Counts,
                              colData = Info,
                              design = ~ Treatment)
dds$Treatment <- relevel(dds$Treatment,"CFA")
dds <- DESeq(dds)

resultsNames(dds)
res <- results(dds,name="Treatment_SYN_vs_CFA")

sortres <- res[order(res$padj),]
head(sortres)

viewme <- as.data.frame(sortres)
write.csv(viewme, file="deseq.csv")

# p-value histogram:
hist(res$pvalue)
hist(res$padj)


##VolcanoPlot

de <- read.csv("deseq.csv")

my_theme <- function(){
  theme_bw() +
    theme(
      plot.title = element_text(face = "bold",size=10,colour = "black"),
      axis.title = element_text(face = "bold",size=10,colour = "black"),
      axis.text = element_text(face = "bold",size=10,colour = "black"),
      axis.text.x = element_text(face = "bold", hjust = 1),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

#theme_set(my_theme())

cleandata = subset(de, !is.na(padj))

cleandata$diffexpressed <- "NO"
cleandata$diffexpressed[cleandata$log2FoldChange > 0.5 & cleandata$padj < 0.05] <- "UP"
cleandata$diffexpressed[cleandata$log2FoldChange < -0.5 & cleandata$padj < 0.05] <- "DOWN"

mycolors <- c("blue", "red", "black", "darkgreen")
names(mycolors) <- c("DOWN", "UP", "NO", "PLOT")

cleandata$cleandatalabel <- NA
cleandata$cleandatalabel[cleandata$diffexpressed != "NO"] <- cleandata$X[cleandata$diffexpressed != "NO"]
write.csv(cleandata, 'cleandata.csv')

volcano_all_data <- read.csv("cleandata_ALL.csv")
volcano_innate <- read.csv("cleandata_innate.csv")
volcano_cytokine <- read.csv("cleandata_cytokine.csv")
volcano_IFNB <- read.csv("cleandata_IFNB.csv")
volcano_Lymphocyte <- read.csv("cleandata_Lymphocyte.csv")
volcano_humoral <- read.csv("cleandata_humoral.csv")

VolcanoPlot <- ggplot(data=volcano_humoral,aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label = cleandatalabel))+ 
  geom_point() + geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text_repel() +
  geom_hline(yintercept=-log10(0.05), col="red") + scale_colour_manual(values = mycolors) + my_theme()

VolcanoPlot
