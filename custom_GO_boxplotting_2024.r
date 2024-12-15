####
##Custom aligned boxplots for Kaushal
##April 16th 2024
####

target_dat <- readRDS("BAT_GO_enrichment_top10_SQUASHING_ALIGN_HUMAN_AGGRESSIVE_target_dat.rds")
grand_EXP_frame <- readRDS("BAT_GO_enrichment_top10_SQUASHING_ALIGN_HUMAN_AGGRESSIVE_grand_EXP_frame.rds")

of_interest <- unique(target_dat$GO_TERM)[c(1, 3, 5)]
#of_interest <- unique(target_dat$GO_TERM)[c(2, 4)]

outfix <- "BAT_GO_enrichment_aligned_human_PROINFLAMM"

of_interest <- unique(target_dat$GO_TERM)[c(2, 4)]
outfix <- "BAT_GO_enrichment_aligned_human_ANTIVIRAL"


target_dat <- target_dat[target_dat$GO_TERM %in% of_interest,]
grand_EXP_frame <- grand_EXP_frame[grand_EXP_frame$GO_TERM %in% of_interest,]

	library(ggplot2)
	library(ggsci)
	test <- ggplot(data = target_dat, aes(x = GO_TERM, y = log2FoldChange, fill = TYPE)) + geom_boxplot(aes(fill = TYPE), position = position_dodge(0.9))
	test <- test + theme_classic() + theme(axis.title.y = element_blank()) + xlab("log2FC Expression over Mock")
	
	test <- test +  geom_hline(yintercept = 0,
    linetype = c("dotted"),
    colour = c("black"),
    size = c(0.5))
    test <- test + coord_flip()
    test <- test + theme(legend.position = "left")

    reverse_col <- T

	if(reverse_col) {
		test <- test + scale_fill_npg()
	}
	
	test_ALIGN <- ggplot(data = grand_EXP_frame, aes(x = GO_TERM, y = log2FoldChange, fill = TYPE)) + geom_boxplot(aes(fill = TYPE), position = position_dodge(0.9))
	test_ALIGN <- test_ALIGN + theme_classic() + theme(axis.title.y = element_blank()) + xlab("log2FC Expression over Mock")
	
	test_ALIGN <- test_ALIGN + geom_hline(yintercept = 0,
    linetype = c("dotted"),
    colour = c("black"),
    size = c(0.5))
    
    test_ALIGN <- test_ALIGN + coord_flip()
    ##colours
	if(!reverse_col) {
		test_ALIGN <- test_ALIGN +  scale_fill_npg()
	}
        #test_ALIGN <- test_ALIGN + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
	##this actually squashes things...
	test_ALIGN <- test_ALIGN + scale_x_discrete(position = "top") ##counterintuitive, but we flipped!
	
	##gridextra.
	
	library(gridExtra)
	library(grid)
	
	margin = theme(plot.margin = unit(c(0.92,0.02,0.02,0.02), "cm"))

	pdf(paste0(outfix, "_logFC_boxplots.pdf"), width = 20, height = 12)
	
	#title <- textGrob("GO TERM", rot = 90, vjust = 1)

	print(grid.arrange(grobs = lapply(list(test, test_ALIGN), "+", margin), ncol = 2))
	
	dev.off()
	
	##CUTTING DOWN AXES.
	
	target_dat_SUBSET <- target_dat
	target_dat_SUBSET$log2FoldChange[target_dat_SUBSET$log2FoldChange > 2] <- 2
	target_dat_SUBSET$log2FoldChange[target_dat_SUBSET$log2FoldChange < -2] <- -2
	
	test <- ggplot(data = target_dat_SUBSET, aes(x = GO_TERM, y = log2FoldChange, fill = TYPE)) + geom_boxplot(aes(fill = TYPE), position = position_dodge(0.9))
	test <- test + theme_classic() + theme(axis.title.y = element_blank()) + xlab("log2FC Expression over Mock")
	
	test <- test +  geom_hline(yintercept = 0,
    linetype = c("dotted"),
    colour = c("black"),
    size = c(0.5))
    test <- test + coord_flip() + theme(legend.position = "left")
	if(reverse_col) {
		test <- test + scale_fill_npg()
	}
	
	grand_EXP_frame_SUBSET <- grand_EXP_frame
	grand_EXP_frame_SUBSET$log2FoldChange[abs(grand_EXP_frame_SUBSET$log2FoldChange) > 2] <- 2
	
	test_ALIGN <- ggplot(data = grand_EXP_frame_SUBSET, aes(x = GO_TERM, y = log2FoldChange, fill = TYPE)) + geom_boxplot(aes(fill = TYPE), position = position_dodge(0.9))
	test_ALIGN <- test_ALIGN + theme_classic() + theme(axis.title.y = element_blank()) + xlab("log2FC Expression over Mock")
	
	test_ALIGN <- test_ALIGN + geom_hline(yintercept = 0,
    linetype = c("dotted"),
    colour = c("black"),
    size = c(0.5))
    test_ALIGN <- test_ALIGN + coord_flip()
    ##colours
    if(!reverse_col) {
    test_ALIGN <- test_ALIGN +  scale_fill_npg()
	}
    ##trying to squash things together visually...
	test_ALIGN <- test_ALIGN + scale_x_discrete(position = "top") ##counterintuitive, but we flipped!

	##gridextra.
	
	margin = theme(plot.margin = unit(c(0.92,0.02,0.02,0.02), "cm"))

	pdf(paste0(outfix, "_logFC_boxplots_LOGFCSQUASH.pdf"), width = 20, height = 12)
	
	print(grid.arrange(grobs = lapply(list(test, test_ALIGN), "+", margin), ncol = 2))
	
	dev.off()
	
