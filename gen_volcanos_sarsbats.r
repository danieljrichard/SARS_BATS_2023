##############
##Comparing human and bat genes
##but instead of boxplots, it was requested for VOLCANOES
##so let's just...see about this.
##August 21st 2023
##
###############

meta_file <- "SARS_infect_human_bat_METADATA_OCT2020_extra.csv"

meta_data <- read.csv(meta_file, header = T)
##NAME
#meta_subset$NAME <- paste0(meta_subset$cell_type, "_", meta_subset$treatment)
meta_subset <- meta_data
meta_subset$NAME <- paste0(meta_subset$species, "_", meta_subset$cell_type, "_", meta_subset$sample_num, "_", meta_subset$treatment)

##NOW. we need to convert.

meta_subset$sample <- paste0(meta_subset$sample_num, "_S", meta_subset$sample_num)

meta_subset$simple <- paste0("S", meta_subset$sample_num, "_", meta_subset$treatment)

bat_files <- readLines("bat_humandef.txt")
#find . -name "*definitions.csv" | grep -v "PROTEOMICS" > bat_humandef.txt

bat_data <- lapply(bat_files, read.csv)
#find . -name "*nameconvert.csv" | grep -v "PROTEOMICS" > human_humandef.txt

human_files <- readLines("human_humandef.txt")
human_dat <- lapply(human_files, read.csv)

##Now I'm fairly confident I can come up with scaled expression values in a lovely heatmap format.
##but I need to first define the genes I'm going to use.
##Someone suggested all proteases, but what about significant genes coming from both species?

human_sigs <- list()
for (x in 1:length(human_dat)) {
    curr_set <- human_dat[[x]]
    curr_set <- curr_set[curr_set$padj < 0.05,]
 #   curr_set <- curr_set[curr_set$log2FoldChange > log2(1.5),]
    curr_set <- curr_set[curr_set$log2FoldChange > 0,]
    human_sigs[[x]] <- curr_set
}

out_counts <- table(unlist(lapply(human_sigs, function(x) x$SYMBOL)))
##want at least in two conditions
out_count_two <- out_counts[out_counts >= 2]

##let's now do the same for bats.

bat_sigs <- list()
bat_data_forhuman <- list()
for (x in 1:length(bat_data)) {
     bat_set <- bat_data[[x]]
     bat_set <- bat_set[!is.na(bat_set$padj),]
    bat_set <- bat_set[bat_set$padj < 0.05,]
    #bat_set <- bat_set[bat_set$log2FoldChange > log2(1.5),]
    bat_set <- bat_set[bat_set$log2FoldChange > 0,]
    for (a in which(is.na(bat_set$ORTHOFIND))) {
        bat_set$ORTHOFIND[a] <- bat_set$gene[a]
    }
    bat_sigs[[x]] <- bat_set
    bat_data_forhuman[[x]] <- bat_data[[x]][!is.na(bat_data[[x]]$ORTHOFIND),]
 }


########
##Aug 28th 2023
##now going to make some really fancy vennplots
##we'll make two, one for each time-point.

bat_names <- unlist(lapply(bat_files, function(x) unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))]))
bat_names <- gsub("_NORMALIZED_counts_matrix_ALLGENES_STATS_orthofinder_human_gene_definitions.csv", "", bat_names)
bat_names <- gsub("BAT_SARS_infection_", "", bat_names)

human_names <- unlist(lapply(human_files, function(x) unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))]))
human_names <- gsub("_NORMALIZED_counts_matrix_ALLGENES_STATS_nameconvert.csv", "", human_names)
human_names <- gsub("human_sars_infection_", "", human_names)

bat_genes <- lapply(bat_sigs, function(x) x$ORTHOFIND)
human_genes <- lapply(human_sigs, function(x) x$SYMBOL)

twelve_set <- list(bat_genes[grepl("12", bat_names)], human_genes[grepl("12", human_names)])
twelve_set_full <- do.call(list, unlist(twelve_set, recursive = FALSE))
names(twelve_set_full) <- c(bat_names[grepl("12", bat_names)], human_names[grepl("12", human_names)])
names(twelve_set_full) <- gsub("_TIME_", "-", names(twelve_set_full))

library(ggvenn)

twelve_plots <- ggvenn(
  twelve_set_full, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

twentyfour_set <- list(bat_genes[grepl("24", bat_names)], human_genes[grepl("24", human_names)])
twentyfour_set_full <- do.call(list, unlist(twentyfour_set, recursive = FALSE))
names(twentyfour_set_full) <- c(bat_names[grepl("24", bat_names)], human_names[grepl("24", human_names)])
names(twentyfour_set_full) <- gsub("_TIME_", "-", names(twentyfour_set_full))

twentyfour_plots <- ggvenn(
  twentyfour_set_full, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

pdf("SARS_BATS_infection_upreg_venndiagram.pdf", width = 12, height = 12)
print(twelve_plots)
print(twentyfour_plots)
dev.off()

##################
try(dir.create("VOLCANOES"))
setwd("VOLCANOES")

lapply(1:length(human_dat), function(x) gen_volcano(human_dat[[x]], gsub("_TIME_", "-", human_names[x]), T))
lapply(1:length(bat_data), function(x) gen_volcano(bat_data[[x]], gsub("_TIME_", "-", bat_names[x]), F))

gen_volcano <- function(dat, outfix, human = F) {
	
	library(EnhancedVolcano)
	#https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
	##each comparison we do, at least for the COVID work, should be only one comparison...

    curr_results <- dat[!is.na(dat$padj),]
	curr_results$padj[curr_results$padj == 0] <- 1e-16
    curr_results <- curr_results[curr_results$log2FoldChange != "",]
    curr_results$log2FoldChange <- as.numeric(curr_results$log2FoldChange)

    if(!human) {
    for (a in which(is.na(curr_results$ORTHOFIND))) {
        curr_results$ORTHOFIND[a] <- curr_results$gene[a]
    }
    curr_results$gene <- curr_results$ORTHOFIND
}else{
    curr_results$gene <- curr_results$SYMBOL
    }
    sig_results <- curr_results[curr_results$padj < 0.05,]
    
	##Let's define some...CUSTOM COLOURS..

	colCustom = rep("grey", dim(curr_results)[1])
	colCustom[which(curr_results$padj < 0.05)] <- "blue"
	colCustom[which(abs(curr_results$log2FoldChange) > 1)] <- "red"
	colCustom[which(curr_results$padj > 0.05)] <- "grey"

	###
	##Update April 20th 2021
	##there's some crazy outliers, so..

	curr_results$log2FoldChange[curr_results$log2FoldChange > 10] <- 10
	curr_results$log2FoldChange[curr_results$log2FoldChange < -10] <- -10

	names(colCustom)[colCustom == "grey"] <- "Insignificant"
	names(colCustom)[colCustom == "blue"] <- "Significant"
	names(colCustom)[colCustom == "red"] <- "Significant-logFC > 1"

	top_genes <- sig_results$gene[order(-abs(sig_results$log2FoldChange))][1:min(c(20, dim(sig_results)[1]))]

	#outfix <- gsub(".rds", "", rds)

	#out_title <- paste0(outfix, "_differential_expression")
	out_title <- outfix
    curr_volc <- EnhancedVolcano(curr_results, lab = curr_results$gene,
	x = "log2FoldChange", y = "padj", pCutoff = 0.05, selectLab = top_genes, FCcutoff = 1, col = c("grey", "grey", "blue", "red"), ylab = bquote(~-Log[10]~adjusted~italic(p)),
	title = NULL, subtitle = out_title, caption = "", legendPosition = "none")
	pdf(paste0(outfix, "_differential_expression_VOLCANO", ".pdf"), width = 12, height = 12)
	print(curr_volc)
	dev.off()
	##############
}



gen_volcano_frame <- function(sigs, label) {

library(ggplot2)
library(dplyr)
padj_cutoff <- 0.05
log2fc_cutoff <- 0.5
sig_frame <- as.data.frame(sigs)
sig_frame <- sig_frame[sig_frame$log2FoldChange != "",]
sig_frame$log2FoldChange <- as.numeric(sig_frame$log2FoldChange)
sig_frame <- sig_frame[!is.na(sig_frame$padj),]
res_table_thres <- sig_frame[!is.na(sig_frame$padj), ] %>% 
  mutate(threshold = padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff)
min(log10(res_table_thres$padj))

res_table_thres$log2FoldChange[res_table_thres$log2FoldChange > 5] <- 5
res_table_thres$log2FoldChange[res_table_thres$log2FoldChange < -5] <- -5

res_table_thres$threshold[res_table_thres$threshold == "TRUE"] <- "red"
res_table_thres$threshold[res_table_thres$threshold == "FALSE"] <- "black"
res_table_thres$threshold[(res_table_thres$padj < padj_cutoff) & (abs(res_table_thres$log2FoldChange) < log2fc_cutoff)] <- "blue"
man_scale <- list("red" = "red", "black" = "black", "blue" = "blue")
print(table(res_table_thres$threshold))

##can we add some labels?
res_table_thres <- res_table_thres[order(res_table_thres$padj),]
res_table_thres$label <- ""
res_table_thres$label[1:10] <- res_table_thres$X[1:10]
## Generate plot
out <- ggplot(res_table_thres) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle(paste0(label, " - old vs young comparison")) +
  xlab("log2 fold change") +
  xlim(min(res_table_thres$log2FoldChange), max(res_table_thres$log2FoldChange)) +
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0, max(-log10(res_table_thres$padj)))) +
  scale_color_manual(values = c("grey60", "red3")) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.3), hjust = 0.5),
        axis.title = element_text(size = rel(1.15)))      

out <- out + geom_text(aes(x = log2FoldChange, y = -log10(padj), label = label))#, nudge_x = 0.25, nudge_y = 0.25, 
#    check_overlap = T)

out <- out + scale_colour_manual(values = man_scale)

return(out)

}



