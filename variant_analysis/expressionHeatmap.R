#####################################
# SARS-CoV-2 Gene Expression Heatmap
#
# Author: Jessica Luc, 2022
#####################################

library(pheatmap)
library(ggplot2)
library(viridis)



tb = read.delim("perGeneLevels.txt",sep='\t',header=T,row.names=1)

tb_trans = t(tb)

metadata = read.table("metadata-humanbat.txt",header=F,sep='\t',row.names=1,stringsAsFactors = FALSE)

annotation_col = data.frame(
                    CellType = factor(metadata[2,]), 
                    Species = factor(metadata[3,]),
                    Time = factor(metadata[5,]),
                    Treatment = factor(metadata[4,], stringsAsFactors=FALSE)
)
rownames(annotation_col) = colnames(tb_trans)
head(annotation_col)

#COLOURING - Modified further in Illustrator
ann_colors = list(
  Treatment = c("mock" = "#FFFFFF", "infected" = "#990000"),
  Time = c("12" = "#414487FF", "24" = "#FFCF40"),
  CellType = c("Calu3" = "#FDE725FF", "A549-hACE2" = "#35B779FF", "Efk3B-hACE2" ="#31688EFF", "pEfLu" = "#440154FF"),
  Species = c("human" = "black", "bat" = "grey")
)

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

pheatmap(log(tb_trans+1,base=2),cluster_cols=F,cluster_rows=F,scale="row",annotation_col=annotation_col,
         annotation_colors=ann_colors, color=myColor)
show(heatmap)
ggsave(file='Expression heatmap-new.pdf', plot=heatmap)

#creating boxplot dataframe

dfmock = data.frame(
  CellType = factor(metadata[2,1:24]), 
  Species = factor(metadata[3,1:24]),
  Treatment = factor(metadata[4,1:24]),
  Time = factor(metadata[5,1:24]),
  GeneLevel = tb[1:24,10]
)
dfinfect = data.frame(
  CellType = factor(metadata[2,25:48]), 
  Species = factor(metadata[3,25:48]),
  Treatment = factor(metadata[4,25:48]),
  Time = factor(metadata[5,25:48]),
  GeneLevel = tb[25:48,10]
)

#fix order of plots
dfinfect$fixedtype = factor(dfinfect$CellType, levels = c('Calu3','A549-hACE2','Efk3B-hACE2','pEfLu'))
dfinfect$fixedtime = factor(dfinfect$Time, levels = c('12 hr','24 hr'))

ggplot(dfinfect, aes(x = Time, y= GeneLevel, fill= Time),) + ggtitle("Infected Samples")+
  geom_boxplot(fill = "light blue", "dark red") + 
  facet_wrap(~CellType, scale="free")

image = ggplot(dfinfect, aes(x = Time, y= GeneLevel, fill= Time),) +
  geom_boxplot(show.legend = FALSE) + 
  scale_fill_manual(breaks = c("12", "24"), 
                    values=c("light blue", "red")) +
  facet_wrap(~fixedtype, scale="free", ncol = 4) +
  labs(y = "N gene expression level") +
  theme(strip.background= element_blank()) +
  theme_bw()
ggsave(file = 'N expression-new.pdf', plot = image)
#Andrew's code
op <- par(mfrow = c(5,4),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)
op <- par(mfrow = c(5,4))
#all genes
for (i in 1:11)
{
  boxplot(tb_trans[i,c(25:27)],tb_trans[i,c(37:39)], col = c("light blue","dark red"),
          main=rownames(tb_trans)[i])
  boxplot(tb_trans[i,c(28:30)],tb_trans[i,c(40:42)], col = c("light blue","dark red"),
          main=rownames(tb_trans)[i])
  boxplot(tb_trans[i,c(31:33)],tb_trans[i,c(43:45)], col = c("light blue","dark red"),
          main=rownames(tb_trans)[i])
  boxplot(tb_trans[i,c(34:36)],tb_trans[i,c(46:48)], col = c("light blue","dark red"),
          main=rownames(tb_trans)[i])
}

boxplot(tb_trans[10,c(25:27)],tb_trans[10,c(37:39)], col = c("light blue","dark red"),
        main="Calu3", ylab= "N gene expression level", xlab = c("12 hr","24 hr") )
boxplot(tb_trans[10,c(28:30)],tb_trans[10,c(40:42)], col = c("light blue","dark red"),
        main="A549-hACE2", ylab= "N gene expression level", xlab = c("12 hr","24 hr"))
boxplot(tb_trans[10,c(31:33)],tb_trans[10,c(43:45)], col = c("light blue","dark red"),
        main="Efk3B-hACE2", ylab= "N gene expression level", xlab = c("12 hr","24 hr"))
boxplot(tb_trans[10,c(34:36)],tb_trans[10,c(46:48)], col = c("light blue","dark red"),
        main="pEfLu", ylab= "N gene expression level", xlab = c("12 hr","24 hr"))
