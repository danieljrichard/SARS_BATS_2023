####
##Visualizing viral protein expression values
##April 16th 2024
####


data <- read.csv("proteomics_viral_proteins.csv")

quant_dat <- data[, 1:(which(colnames(data) == "C..Only.identified.by.site")-1)]

rownames(quant_dat) <- data[, "T..Protein.IDs"]
##going to make a line plot with error bars and individual points...

plot_set <- list()
for (prot in rownames(quant_dat)) {
    ##need one per protein, of course.
    plot_frame <- data.frame(type = do.call("c", list(rep("MOCK", length(which(grepl("Mock", colnames(quant_dat))))),
        rep("24h", length(which(grepl(".24h", colnames(quant_dat), fixed = T)))),
        rep("48h", length(which(grepl(".48h", colnames(quant_dat), fixed = T)))))),
        LFQ = do.call("c", lapply(list(quant_dat[prot, which(grepl("Mock", colnames(quant_dat)))],
            quant_dat[prot, which(grepl(".24h", colnames(quant_dat), fixed = T))],
            quant_dat[prot, which(grepl(".48h", colnames(quant_dat), fixed = T))]
        ), as.numeric)))
    library(ggplot2)
    pd <- position_dodge(0.1) # move them .05 to the left and right
    library(dplyr)

    plot_frame$LFQ <- log2(plot_frame$LFQ)
    simple_frame <- plot_frame %>% group_by(type) %>% summarize_at("LFQ", c("mean", "sd"))
  simple_frame$TIME <- 24
     simple_frame$TIME[ simple_frame$type == "48h"] <- 48
   
    plot_frame_cut <- plot_frame[plot_frame$type != "MOCK",]
    curr_plot <- ggplot(plot_frame_cut, aes(x = type, y = LFQ, color = type)) + geom_boxplot() + geom_jitter() + theme_classic() + xlab("") + ylab("log2 LFQ") + ggtitle(tail(unlist(strsplit(prot, "|", fixed = T)),1))
    plot_set[[prot]] <- curr_plot
    print(prot)
}

library(gridExtra)
pdf("viral_protein_expression_boxplots.pdf", width = 12, height = 4)
grid.arrange(grobs = plot_set, nrow = 2, ncol = 4)
dev.off()