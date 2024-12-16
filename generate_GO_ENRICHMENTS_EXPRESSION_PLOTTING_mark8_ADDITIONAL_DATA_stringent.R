#################
##
##
##Comparing GO enrichment sets
##February 9th 2021
##################

###################################

masta <- "bat_stringent_metadata.csv"
global_map <- "/extra/SARS_BATS/CUSTOM_GO_ENRICHMENTS/BAT/BAT_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
##There we go.
##one other thing we need: gene expression data:
#find $PWD/. -name "*.csv" | grep "STATS.csv" > BAT_final_timepoints_statsfiles.txt
stats_files <- "BAT_final_timepoints_statsfiles.txt"
outfix <- "BAT_GO_enrichment_top10_SQUASHING_STRINGENT"
consolidate_species_level_GO_terms(masta, global_map, stats_files, outfix, topn=10)
consolidate_species_level_GO_terms(masta, global_map, stats_files, paste0(outfix, "_AGGRESSIVE"), topn=10, AGGRESSIVE = T)


###
##Update March 13th 2021
global_map <- "/extra/SARS_BATS/CUSTOM_GO_ENRICHMENTS/HUMAN/HUMAN_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
stats_files <- "HUMAN_final_timepoints_statsfiles.txt"
#find $PWD/. -name "*STATS.csv" > HUMAN_final_timepoints_statsfiles.txt
outfix <- "HUMAN_GO_enrichment_top10_SQUASHING_STRINGENT"
masta <- "human_meta_STRINGENT.csv"
topn <- 10
consolidate_species_level_GO_terms(masta, global_map, stats_files, outfix, topn, revcol = T)
consolidate_species_level_GO_terms(masta, global_map, stats_files, paste0(outfix, "_AGGRESSIVE"), topn, AGGRESSIVE = T, revcol = T)

consolidate_species_level_GO_terms <- function(masta, global_map, stats_files, outfix, topn = 10, AGGRESSIVE = F, revcol = F) {
	
	meta <- read.csv(masta, stringsAsFactors = FALSE)
	
	csv_collapse <- unlist(c(meta$file1, meta$file2))
	csv_collapse <- csv_collapse[csv_collapse != ""]
    csv_collapse <- csv_collapse[!grepl("DOWNREG", csv_collapse)]
    topn <- topn*2
	csv_dat <- lapply(csv_collapse, read.csv)
	
	csv_filt <- lapply(csv_dat, function(x) x[x$PADJ < 0.05,])
	all_hits_ever <- unique(unlist(sapply(csv_filt, function(x) x$GO.ID)))
    library(data.table)
	if(AGGRESSIVE) {
		csv_ULTRAMEGA <- as.data.frame(rbindlist(csv_filt))
	}
	csv_filt <- lapply(csv_filt, function(x) x[order(x$PADJ),][1:topn,])
	
	library(data.table)
	csv_UPREG <- as.data.frame(rbindlist(csv_filt[grep("UPREG", csv_collapse)]))
	csv_DOWNREG <- as.data.frame(rbindlist(csv_filt[grep("DOWNREG", csv_collapse)]))
	
	csv_MEGA <- as.data.frame(rbindlist(csv_filt)) ##save for later labelling...
	
	collapse_UPREG <- table(csv_UPREG$GO.ID)
	collapse_UPREG <- collapse_UPREG[order(-collapse_UPREG)]
	collapse_DOWNREG <- table(csv_DOWNREG$GO.ID)
	collapse_DOWNREG <- collapse_DOWNREG[order(-collapse_DOWNREG)]
	
	if(AGGRESSIVE) {
	
	library(data.table)
	glob_GO <- as.data.frame(fread(global_map, header = F))
	
	stats_dat <- readLines(stats_files)
    stats_dat <- stats_dat[!grepl("EF_LU", stats_dat)]
	stats_tables <- lapply(stats_dat, read.csv)
	
	gimme_name <- function(csv1) {
	out1 <- unlist(strsplit(csv1, "/"))[length(unlist(strsplit(csv1, "/")))]
	out1 <- unlist(strsplit(out1, "_infection_"))[2]
	out1 <- unlist(strsplit(out1, "_NORMALIZED"))[1]
	return(out1)
	}
	
	stats_names <- sapply(stats_dat, gimme_name)
	
	ALL_genesets <- lapply(all_hits_ever, function(x) glob_GO[grep(x, glob_GO$V2), "V1"])
	
	names(ALL_genesets) <- all_hits_ever
	
	for (x in 1:length(stats_tables)) {
		stats_tables[[x]]$TYPE <- stats_names[x] 
	}
	
	subset_EXP_data <- lapply(ALL_genesets, function(x) rbindlist(lapply(1:length(stats_tables), function(y) 
			subset(stats_tables[[y]], stats_tables[[y]]$gene %in% x)[, c("log2FoldChange", "TYPE")])))

	for (x in 1:length(ALL_genesets)) {
		subset_EXP_data[[x]]$GO <- names(ALL_genesets)[x]
		
	}
	
	super_grand_EXP_frame <- rbindlist(subset_EXP_data)
	##let's order by expression...
	
	per_term_EXP <- unlist(lapply(unique(super_grand_EXP_frame$GO), function(x) mean(subset(super_grand_EXP_frame, super_grand_EXP_frame$GO %in% x)$log2FoldChange)))
	names(per_term_EXP) <- unique(super_grand_EXP_frame$GO)
	per_term_EXP <- sort(per_term_EXP)
	collapse_DOWNREG <- per_term_EXP[1:(topn/2)]
	collapse_UPREG <- rev(per_term_EXP)[1:(topn/2)]
	
}
	
	UPREG_genes <- setdiff(names(collapse_UPREG), names(collapse_DOWNREG))[1:(topn/2)]
	DOWNREG_genes <- setdiff(names(collapse_DOWNREG), names(collapse_UPREG))[1:(topn/2)]
	if(length(intersect(names(collapse_UPREG), names(collapse_DOWNREG))) != 0) {
		print("that's gonna be a problemo chiefo")
		#browser()
	}
	library(data.table)
	glob_GO <- as.data.frame(fread(global_map, header = F))
	
	stats_dat <- readLines(stats_files)
    stats_dat <- stats_dat[!grepl("EF_LU", stats_dat)]

	stats_tables <- lapply(stats_dat, read.csv)
	
	gimme_name <- function(csv1) {
	out1 <- unlist(strsplit(csv1, "/"))[length(unlist(strsplit(csv1, "/")))]
	out1 <- unlist(strsplit(out1, "_infection_"))[2]
	out1 <- unlist(strsplit(out1, "_NORMALIZED"))[1]
	return(out1)
	}
	
	stats_names <- sapply(stats_dat, gimme_name)
	
	UPREG_genesets <- lapply(UPREG_genes, function(x) glob_GO[grep(x, glob_GO$V2), "V1"])
	DOWNREG_genesets <- lapply(DOWNREG_genes, function(x) glob_GO[grep(x, glob_GO$V2), "V1"])
	
	ALL_genesets <- do.call("c", list(UPREG_genesets, DOWNREG_genesets))
	names(ALL_genesets) <- c(UPREG_genes, DOWNREG_genes)
	
	for (x in 1:length(stats_tables)) {
		stats_tables[[x]]$TYPE <- stats_names[x] 
	}
	
	
	ARINJAY_DATA <- T
	if(ARINJAY_DATA) {
		
		stats_tables_list <- list()
		for (x in 1:length(stats_tables)) {
			curr_thing <- data.table(stats_tables[[x]])
			setkey(curr_thing, "gene")
			stats_tables_list[[x]] <- curr_thing
		}
		
		global_genes <- unique(unlist(lapply(stats_tables_list, function(x) x$gene)))
		
		grand_collapse_stats <- do.call("cbind", lapply(stats_tables_list, function(x) x[.(global_genes)][, 2:7]))
		grand_collapse_stats$gene <- global_genes
		
		subset_EXP_data_custom <- lapply(ALL_genesets, function(x) subset(grand_collapse_stats, grand_collapse_stats$gene %in% x))
	
	##Now, gotta just label.
	for (x in 1:length(ALL_genesets)) {
		subset_EXP_data_custom[[x]]$GO <- names(ALL_genesets)[x]
		
	}
	
	grand_EXP_frame_CUSTOM <- rbindlist(subset_EXP_data_custom)
	
	
	if(AGGRESSIVE) {
		
	map_GO <- unlist(sapply(unique(grand_EXP_frame_CUSTOM$GO), function(x) csv_ULTRAMEGA[which(csv_ULTRAMEGA$GO.ID == x)[1], "Term"]))
	}else{
	map_GO <- unlist(sapply(unique(grand_EXP_frame_CUSTOM$GO), function(x) csv_MEGA[which(csv_MEGA$GO.ID == x)[1], "Term"]))
	}
	
	#names(per_term_EXP) <- unique(grand_EXP_frame$GO)
	grand_EXP_frame_CUSTOM$GO_TERM <- unlist(lapply(grand_EXP_frame_CUSTOM$GO, function(x) map_GO[[x]]))	
	
	fwrite(grand_EXP_frame_CUSTOM, paste0(outfix, "_grand_EXP_frame_CUSTOM_APRIL2021.csv"))

		
		
	}


	subset_EXP_data <- lapply(ALL_genesets, function(x) rbindlist(lapply(1:length(stats_tables), function(y) 
			subset(stats_tables[[y]], stats_tables[[y]]$gene %in% x)[, c("log2FoldChange", "TYPE")])))
	
	##Now, gotta just label.
	for (x in 1:length(ALL_genesets)) {
		subset_EXP_data[[x]]$GO <- names(ALL_genesets)[x]
		
	}
	
	grand_EXP_frame <- rbindlist(subset_EXP_data)
	grand_EXP_frame$TYPE <- factor(grand_EXP_frame$TYPE, stats_names[order(stats_names)]) ##force alphabetical order in boxplots...
		
	##and now FINALLY...geom_boxplot
	
	##let's order by expression...
	
	per_term_EXP <- unlist(lapply(unique(grand_EXP_frame$GO), function(x) mean(subset(grand_EXP_frame, grand_EXP_frame$GO %in% x)$log2FoldChange)))
	
	if(AGGRESSIVE) {
	map_GO <- unlist(sapply(unique(grand_EXP_frame$GO), function(x) csv_ULTRAMEGA[which(csv_ULTRAMEGA$GO.ID == x)[1], "Term"]))
	}else{
	map_GO <- unlist(sapply(unique(grand_EXP_frame$GO), function(x) csv_MEGA[which(csv_MEGA$GO.ID == x)[1], "Term"]))
	}
	
	#names(per_term_EXP) <- unique(grand_EXP_frame$GO)
	names(per_term_EXP) <- map_GO
	per_term_EXP <- sort(per_term_EXP, decreasing = F)
	grand_EXP_frame$GO_TERM <- unlist(lapply(grand_EXP_frame$GO, function(x) map_GO[[x]]))	
	grand_EXP_frame$GO_TERM <- factor(grand_EXP_frame$GO_TERM, names(per_term_EXP))
	
	saveRDS(grand_EXP_frame, paste0(outfix, "_grand_EXP_frame.rds"))

	
	library(ggplot2)
	library(forcats)
	test <- ggplot(data = grand_EXP_frame, aes(x = GO_TERM, y = log2FoldChange, fill = TYPE)) + geom_boxplot(aes(fill = TYPE), position = position_dodge(0.9))
	test <- test + theme_classic() + theme(axis.title.y = element_blank()) + xlab("log2FC Expression over Mock")
	
	test <- test +  geom_hline(yintercept = 0,
    linetype = c("dotted"),
    colour = c("black"),
    size = c(0.5))
    test <- test + coord_flip()
	
	if(revcol) {
		library(ggsci)
		test <- test + scale_fill_npg()
	}
	
	pdf(paste0(outfix, "_collapsed_topn10_GO_term_expression_barplots.pdf"), width = 14, height = 12)
	print(test)
	dev.off()
	
	##Let's just make a reduced-axis plot...
	
	grand_EXP_frame_REDUCED <- grand_EXP_frame
	grand_EXP_frame_REDUCED[grand_EXP_frame_REDUCED$log2FoldChange > 2, "log2FoldChange"] <- 2
	grand_EXP_frame_REDUCED[grand_EXP_frame_REDUCED$log2FoldChange < -2, "log2FoldChange"] <- -2
	
	test_REDUCED <- ggplot(data = grand_EXP_frame_REDUCED, aes(x = GO_TERM, y = log2FoldChange, fill = TYPE)) + geom_boxplot(aes(fill = TYPE), position = position_dodge(0.9))
	test_REDUCED <- test_REDUCED + theme_classic() + theme(axis.title.y = element_blank()) + xlab("log2FC Expression over Mock")
	
	test_REDUCED <- test_REDUCED +  geom_hline(yintercept = 0,
    linetype = c("dotted"),
    colour = c("black"),
    size = c(0.5))
    test_REDUCED <- test_REDUCED + coord_flip()
	
	if(revcol) {
		library(ggsci)
		test_REDUCED <- test_REDUCED + scale_fill_npg()
	}
	##I like this a lot actually...
	pdf(paste0(outfix, "_collapsed_topn10_GO_term_expression_barplots_FCSQUASH.pdf"), width = 14, height = 12)
	print(test_REDUCED)
	dev.off()
	
}

##aligning HUMAN to BAT hits:

target_RDS <- "BAT_GO_enrichment_top10_SQUASHING_STRINGENT_grand_EXP_frame.rds"
global_map <- "/extra/SARS_BATS/CUSTOM_GO_ENRICHMENTS/HUMAN/HUMAN_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
stats_files <- "HUMAN_final_timepoints_statsfiles.txt"
#find $PWD/. -name "*STATS.csv" > HUMAN_final_timepoints_statsfiles.txt
outfix <- "BAT_GO_enrichment_STRINGENT_top10_SQUASHING_ALIGN_HUMAN"

ALIGN_CROSS_SPECIES_FANCYSAUCE(target_RDS, global_map, stats_files, outfix)

target_RDS <- "HUMAN_GO_enrichment_top10_SQUASHING_STRINGENT_grand_EXP_frame.rds"
if(FALSE) {
	test <- readRDS("HUMAN_GO_enrichment_top10_SQUASHING_STRINGENT_grand_EXP_frame.rds")
	library(dplyr)
	test2 <- test %>% group_by(GO_TERM) %>% summarise(mean = mean(log2FoldChange), sd = sd(log2FoldChange))
	##on average upregulated
	test3 <- test[test$GO_TERM %in% test2$GO_TERM[test2$mean > 0],]
	test3$GO_TERM <- droplevels(test3$GO_TERM)
	saveRDS(test3, "HUMAN_GO_enrichment_top10_SQUASHING_STRINGENT_grand_EXP_frame_upreg.rds")
}

target_RDS <- "HUMAN_GO_enrichment_top10_SQUASHING_STRINGENT_grand_EXP_frame_upreg.rds"
global_map <- "/extra/SARS_BATS/CUSTOM_GO_ENRICHMENTS/BAT/BAT_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
#grep -v  "EF_LU" BAT_final_timepoints_statsfiles.txt > BAT_final_timepoints_statsfiles_noEF_LU.txt
stats_files <- "BAT_final_timepoints_statsfiles_noEF_LU.txt"
outfix <- "HUMAN_GO_enrichment_STRINGENT_upreg_top10_SQUASHING_ALIGN_BAT"
ALIGN_CROSS_SPECIES_FANCYSAUCE(target_RDS, global_map, stats_files, outfix, reverse_col = T)

stats_files <- "HUMAN_final_timepoints_statsfiles.txt"
global_map <- "/extra/SARS_BATS/CUSTOM_GO_ENRICHMENTS/HUMAN/HUMAN_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
#find $PWD/. -name "*STATS.csv" > HUMAN_final_timepoints_statsfiles.txt
outfix <- "BAT_GO_enrichment_STRINGENT_top10_SQUASHING_ALIGN_HUMAN"


ALIGN_CROSS_SPECIES_FANCYSAUCE <- function(target_RDS, global_map, stats_files, outfix, reverse_col = F) {

	target_dat <- readRDS(target_RDS)

	library(data.table)
	glob_GO <- as.data.frame(fread(global_map, header = F))
	
	stats_dat <- readLines(stats_files)
	stats_tables <- lapply(stats_dat, read.csv)
	
	gimme_name <- function(csv1) {
	out1 <- unlist(strsplit(csv1, "/"))[length(unlist(strsplit(csv1, "/")))]
	out1 <- unlist(strsplit(out1, "_infection_"))[2]
	out1 <- unlist(strsplit(out1, "_NORMALIZED"))[1]
	return(out1)
	}
	
	stats_names <- sapply(stats_dat, gimme_name)
	
	target_GENES <- unique(target_dat$GO)
	
	ALL_genesets <- lapply(target_GENES, function(x) glob_GO[grep(x, glob_GO$V2), "V1"])
	names(ALL_genesets) <- target_GENES
	
	for (x in 1:length(stats_tables)) {
		stats_tables[[x]]$TYPE <- stats_names[x] 
	}
	
	
	ARINJAY_DATA <- T
	if(ARINJAY_DATA) {
		
		stats_tables_list <- list()
		for (x in 1:length(stats_tables)) {
			curr_thing <- data.table(stats_tables[[x]])
			setkey(curr_thing, "gene")
			stats_tables_list[[x]] <- curr_thing
		}
		
		global_genes <- unique(unlist(lapply(stats_tables_list, function(x) x$gene)))
		
		grand_collapse_stats <- do.call("cbind", lapply(stats_tables_list, function(x) x[.(global_genes)][, 2:7]))
		grand_collapse_stats$gene <- global_genes
		
		subset_EXP_data_custom <- lapply(ALL_genesets, function(x) subset(grand_collapse_stats, grand_collapse_stats$gene %in% x))
	
	##Now, gotta just label.
	for (x in 1:length(ALL_genesets)) {
		subset_EXP_data_custom[[x]]$GO <- names(ALL_genesets)[x]
		
	}
	
	grand_EXP_frame_CUSTOM <- rbindlist(subset_EXP_data_custom)
	
	
	convert_frame <- unique(target_dat[, c("GO", "GO_TERM")])
	convert_vect <- convert_frame$GO_TERM
	names(convert_vect) <- convert_frame$GO
	
	grand_EXP_frame_CUSTOM$GO_TERM <- unlist(sapply(grand_EXP_frame_CUSTOM$GO, function(x) convert_vect[[x]]))
	
	grand_EXP_frame_CUSTOM$GO_TERM <- factor(grand_EXP_frame_CUSTOM$GO_TERM, levels = levels(target_dat$GO_TERM)) ##this may be the first time in history I've actually used levels()
	
	fwrite(grand_EXP_frame_CUSTOM, paste0(outfix, "_grand_EXP_frame_CUSTOM_APRIL2021.csv"))
		
	}


	
	subset_EXP_data <- lapply(ALL_genesets, function(x) rbindlist(lapply(1:length(stats_tables), function(y) 
			subset(stats_tables[[y]], stats_tables[[y]]$gene %in% x)[, c("log2FoldChange", "TYPE")])))
	
	for (x in 1:length(ALL_genesets)) {
		subset_EXP_data[[x]]$GO <- names(ALL_genesets)[x]
		
	}
	
	grand_EXP_frame <- rbindlist(subset_EXP_data)
	grand_EXP_frame$TYPE <- factor(grand_EXP_frame$TYPE, stats_names[order(stats_names)]) ##force alphabetical order in boxplots...
		
	convert_frame <- unique(target_dat[, c("GO", "GO_TERM")])
	convert_vect <- convert_frame$GO_TERM
	names(convert_vect) <- convert_frame$GO
	
	grand_EXP_frame$GO_TERM <- unlist(sapply(grand_EXP_frame$GO, function(x) convert_vect[[x]]))
	
	grand_EXP_frame$GO_TERM <- factor(grand_EXP_frame$GO_TERM, levels = levels(target_dat$GO_TERM)) ##this may be the first time in history I've actually used levels()
	
saveRDS(target_dat, paste0(outfix, "_target_dat.rds"))

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
	if(reverse_col) {
		test <- test + scale_fill_npg()
	}
	
	saveRDS(grand_EXP_frame, paste0(outfix, "_grand_EXP_frame.rds"))

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
	

}

##cleanup human output

if(FALSE) {
	test <- readRDS("HUMAN_GO_enrichment_top10_SQUASHING_STRINGENT_grand_EXP_frame_upreg.rds")
	human_dat <- read.csv("HUMAN_GO_enrichment_top10_SQUASHING_STRINGENT_grand_EXP_frame_CUSTOM_APRIL2021.csv")
	human_dat_cut <- human_dat[human_dat$GO %in% test$GO,]
	write.csv(human_dat_cut, "HUMAN_GO_enrichment_top10_SQUASHING_STRINGENT_grand_EXP_frame_CUSTOM_APRIL2021_cut.csv", row.names = F)
}

##OK...
convert_genes("HUMAN_GO_enrichment_top10_SQUASHING_grand_EXP_frame_CUSTOM_APRIL2021.csv")

convert_genes("BAT_GO_enrichment_top10_SQUASHING_ALIGN_HUMAN_grand_EXP_frame_CUSTOM_APRIL2021.csv")

convert_genes <- function(csv) {
	library(clusterProfiler)
	dat <- read.csv(csv)
	conv_genes <- bitr(dat$gene, fromType = "ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
	dat <- subset(dat, dat$gene %in% conv_genes$ENSEMBL)
	colnames(conv_genes)[1] <- "gene"
	dat <- merge(dat, conv_genes, by = "gene")
	dat <- dat[order(dat[, "GO_TERM"]),]
	write.csv(dat, gsub(".csv", "_convert.csv", csv), row.names = F)
}

