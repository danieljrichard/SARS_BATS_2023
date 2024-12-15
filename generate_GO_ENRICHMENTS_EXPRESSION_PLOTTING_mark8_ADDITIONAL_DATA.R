#################
##
##
##Comparing GO enrichment sets b/c I'm awesome.
##February 9th 2021
##################

###################################
##UPDATE MARCH 13th 2021
##while this shit is good, and I think I'll retain it for if I ever need to directly compare GO enrichments
##between two sets of gene-expression data (which is actually a pretty cool method...ACTUALLY.
##Now I'm going to try to weave in more of the expression data itself...

masta <- "BAT_BP_ENRICHMENT_COMPARISONS_METATABLE.csv"
global_map <- "/extra/SARS_BATS/CUSTOM_GO_ENRICHMENTS/BAT/BAT_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
##There we go.
##one other thing we need: gene expression data:
#find $PWD/. -name "*.csv" | grep "STATS.csv" > BAT_final_timepoints_statsfiles.txt
stats_files <- "BAT_final_timepoints_statsfiles.txt"
outfix <- "BAT_GO_enrichment_top10_SQUASHING"
consolidate_species_level_GO_terms(masta, global_map, stats_files, outfix, topn=10)
consolidate_species_level_GO_terms(masta, global_map, stats_files, paste0(outfix, "_AGGRESSIVE"), topn=10, AGGRESSIVE = T)

##..let's just first take 'any GO term from any timepoint',
##rather than the more complicated 'terms only found at 24H timepoint' or something like that...

###
##Update March 13th 2021
##Now doing das humans
global_map <- "/extra/SARS_BATS/CUSTOM_GO_ENRICHMENTS/HUMAN/HUMAN_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
stats_files <- "HUMAN_final_timepoints_statsfiles.txt"
#find $PWD/. -name "*STATS.csv" > HUMAN_final_timepoints_statsfiles.txt
outfix <- "HUMAN_GO_enrichment_top10_SQUASHING"
masta <- "HUMAN_BP_MEGATABLE.csv"
topn <- 10
consolidate_species_level_GO_terms(masta, global_map, stats_files, outfix, topn, revcol = T)
consolidate_species_level_GO_terms(masta, global_map, stats_files, paste0(outfix, "_AGGRESSIVE"), topn, AGGRESSIVE = T, revcol = T)

consolidate_species_level_GO_terms <- function(masta, global_map, stats_files, outfix, topn = 10, AGGRESSIVE = F, revcol = F) {
	
	meta <- read.csv(masta, stringsAsFactors = FALSE)
	
	csv_collapse <- unlist(c(meta$file1, meta$file2))
	
	csv_dat <- lapply(csv_collapse, read.csv)
	
	##note to future self: can also just take the top 5 GO terms from the global set.
	##or b/c I'm fancy:
	csv_filt <- lapply(csv_dat, function(x) x[x$PADJ < 0.05,])
	all_hits_ever <- unique(unlist(sapply(csv_filt, function(x) x$GO.ID)))
	if(AGGRESSIVE) {
		csv_ULTRAMEGA <- as.data.frame(rbindlist(csv_filt))
	}
	csv_filt <- lapply(csv_filt, function(x) x[order(x$PADJ),][1:topn,])
	
	library(data.table)
	csv_UPREG <- as.data.frame(rbindlist(csv_filt[grep("UPREG", csv_collapse)]))
	csv_DOWNREG <- as.data.frame(rbindlist(csv_filt[grep("DOWNREG", csv_collapse)]))
	
	csv_MEGA <- as.data.frame(rbindlist(csv_filt)) ##save for later labelling...
	
	#collapse_UPREG <- table(unlist(lapply(csv_UPREG, function(x) x$GO.ID)))
	#collapse_DOWNREG <- table(unlist(lapply(csv_DOWNREG, function(x) x$GO.ID)))
	
	collapse_UPREG <- table(csv_UPREG$GO.ID)
	collapse_UPREG <- collapse_UPREG[order(-collapse_UPREG)]
	collapse_DOWNREG <- table(csv_DOWNREG$GO.ID)
	collapse_DOWNREG <- collapse_DOWNREG[order(-collapse_DOWNREG)]
	
	if(AGGRESSIVE) {
	
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
	##...do I force the same genes regardless of significance? My gut says YES...
	#...just grab all the GO terms first.
	##I'm just defining things.
	
	ALL_genesets <- lapply(all_hits_ever, function(x) glob_GO[grep(x, glob_GO$V2), "V1"])
	
	names(ALL_genesets) <- all_hits_ever
	
	##this is actually a lot, A LOT, more complicated than I thought. 
	##gonna have to switch to ggplot2's paired boxplots, maybe with coord_flip()
	##das OK.
	
	for (x in 1:length(stats_tables)) {
		stats_tables[[x]]$TYPE <- stats_names[x] 
	}
	
	subset_EXP_data <- lapply(ALL_genesets, function(x) rbindlist(lapply(1:length(stats_tables), function(y) 
			subset(stats_tables[[y]], stats_tables[[y]]$gene %in% x)[, c("log2FoldChange", "TYPE")])))
	
	##can we get a sense for gene drop-out in different sets?
	#sapply(subset_EXP_data, function(x) table(x$TYPE))
	##it's not awful. keep in mind that gene dropout can happen because of quality/coverage in some samples...
	
	##Now, gotta just label.
	for (x in 1:length(ALL_genesets)) {
		subset_EXP_data[[x]]$GO <- names(ALL_genesets)[x]
		
	}
	
	super_grand_EXP_frame <- rbindlist(subset_EXP_data)
	#super_grand_EXP_frame$TYPE <- factor(grand_EXP_frame$TYPE, stats_names[order(stats_names)]) ##force alphabetical order in boxplots...
		
	##and now FINALLY...geom_boxplot
	
	##let's order by expression...
	
	per_term_EXP <- unlist(lapply(unique(super_grand_EXP_frame$GO), function(x) mean(subset(super_grand_EXP_frame, super_grand_EXP_frame$GO %in% x)$log2FoldChange)))
	names(per_term_EXP) <- unique(super_grand_EXP_frame$GO)
	per_term_EXP <- sort(per_term_EXP)
	collapse_DOWNREG <- per_term_EXP[1:(topn/2)]
	collapse_UPREG <- rev(per_term_EXP)[1:(topn/2)]
	
}
	
	##...um. let's prioritize GO terms which are shared amongst the majority of datasets...shall we?
	
	#length(intersect(UPREG_genes, DOWNREG_genes))
	##..that's pretty...unexpected...
	
	#UPREG_genes <- names(collapse_UPREG)[1:(topn/2)]
	#DOWNREG_genes <- names(collapse_DOWNREG)[1:(topn/2)]
	
	UPREG_genes <- setdiff(names(collapse_UPREG), names(collapse_DOWNREG))[1:(topn/2)]
	DOWNREG_genes <- setdiff(names(collapse_DOWNREG), names(collapse_UPREG))[1:(topn/2)]
	##^^paranoid. but we did catch a me-writing-copypasta error so that's nice.
	if(length(intersect(names(collapse_UPREG), names(collapse_DOWNREG))) != 0) {
		print("that's gonna be a problemo chiefo")
		#browser()
	#I'm over it.
	}
	##majority are only shared with two...weird...whatever.
	##hmm...hmm.
	##I could...plot the p-values for all the terms on the boxplot as well...with a red/blue colour scale? Or is that too much?
	##TOO MUCH INFO. MAKE DECISIONS.
	##I wonder if we should force 5 upreg and 5 downreg? I think maybe we should...
	## I don't know what genes are assoc...
	##let's pull out all the genes associated with terms, then go into the gene expression csv files and pull them out...
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
	##...do I force the same genes regardless of significance? My gut says YES...
	#...just grab all the GO terms first.
	##I'm just defining things.
	
	UPREG_genesets <- lapply(UPREG_genes, function(x) glob_GO[grep(x, glob_GO$V2), "V1"])
	DOWNREG_genesets <- lapply(DOWNREG_genes, function(x) glob_GO[grep(x, glob_GO$V2), "V1"])
	
	##at this point I can probably collapse these down...
	ALL_genesets <- do.call("c", list(UPREG_genesets, DOWNREG_genesets))
	names(ALL_genesets) <- c(UPREG_genes, DOWNREG_genes)
	
	##this is actually a lot, A LOT, more complicated than I thought. 
	##gonna have to switch to ggplot2's paired boxplots, maybe with coord_flip()
	##das OK.
	
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
	
	##can we get a sense for gene drop-out in different sets?
	#sapply(subset_EXP_data, function(x) table(x$TYPE))
	##it's not awful. keep in mind that gene dropout can happen because of quality/coverage in some samples...
	
	##Now, gotta just label.
	for (x in 1:length(ALL_genesets)) {
		subset_EXP_data_custom[[x]]$GO <- names(ALL_genesets)[x]
		
	}
	
	grand_EXP_frame_CUSTOM <- rbindlist(subset_EXP_data_custom)
	
	
	##one last thingy...
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
	
	##can we get a sense for gene drop-out in different sets?
	#sapply(subset_EXP_data, function(x) table(x$TYPE))
	##it's not awful. keep in mind that gene dropout can happen because of quality/coverage in some samples...
	
	##Now, gotta just label.
	for (x in 1:length(ALL_genesets)) {
		subset_EXP_data[[x]]$GO <- names(ALL_genesets)[x]
		
	}
	
	grand_EXP_frame <- rbindlist(subset_EXP_data)
	grand_EXP_frame$TYPE <- factor(grand_EXP_frame$TYPE, stats_names[order(stats_names)]) ##force alphabetical order in boxplots...
		
	##and now FINALLY...geom_boxplot
	
	##let's order by expression...
	
	per_term_EXP <- unlist(lapply(unique(grand_EXP_frame$GO), function(x) mean(subset(grand_EXP_frame, grand_EXP_frame$GO %in% x)$log2FoldChange)))
	
	##one last thingy...
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
	
	##I like this a lot actually...
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
	
	##Note to self: I'll want to definitely retain these RDS data for later when I'm directly plotting the two species relative to one another...
	##though I guess not since they'll have different contents...
	##I'll retain the code used to subset the expression data based on...GO terms. YEAH.
}

##MARCH13th 2021
##OK. The next function is to align one species to the output from the above on another species...
##so I can save a bit of time.
##yes one could make the criticism that I could make a more generalized function plotting all the data from both species
##for a single set of GO terms (especially given that all the species formatting is the exact same)...but this isn't a super-critical thing
#right now.


##aligning HUMAN to BAT hits:

target_RDS <- "BAT_GO_enrichment_top10_SQUASHING_grand_EXP_frame.rds"
global_map <- "/extra/SARS_BATS/CUSTOM_GO_ENRICHMENTS/HUMAN/HUMAN_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
stats_files <- "HUMAN_final_timepoints_statsfiles.txt"
#find $PWD/. -name "*STATS.csv" > HUMAN_final_timepoints_statsfiles.txt
outfix <- "BAT_GO_enrichment_top10_SQUASHING_ALIGN_HUMAN"

ALIGN_CROSS_SPECIES_FANCYSAUCE(target_RDS, global_map, stats_files, outfix)
ALIGN_CROSS_SPECIES_FANCYSAUCE("BAT_GO_enrichment_top10_SQUASHING_AGGRESSIVE_grand_EXP_frame.rds", global_map, stats_files, outfix = paste0(outfix, "_AGGRESSIVE"))

target_RDS <- "HUMAN_GO_enrichment_top10_SQUASHING_grand_EXP_frame.rds"
global_map <- "/extra/SARS_BATS/CUSTOM_GO_ENRICHMENTS/BAT/BAT_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
stats_files <- "BAT_final_timepoints_statsfiles.txt"
outfix <- "HUMAN_GO_enrichment_top10_SQUASHING_ALIGN_BAT"
ALIGN_CROSS_SPECIES_FANCYSAUCE(target_RDS, global_map, stats_files, outfix, reverse_col = T)
ALIGN_CROSS_SPECIES_FANCYSAUCE(target_RDS= "HUMAN_GO_enrichment_top10_SQUASHING_AGGRESSIVE_grand_EXP_frame.rds", global_map, stats_files, outfix = paste0(outfix, "_AGGRESSIVE"), reverse_col = T)

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
	##...do I force the same genes regardless of significance? My gut says YES...
	#...just grab all the GO terms first.
	##I'm just defining things.
	
	target_GENES <- unique(target_dat$GO)
	
	ALL_genesets <- lapply(target_GENES, function(x) glob_GO[grep(x, glob_GO$V2), "V1"])
	names(ALL_genesets) <- target_GENES
	
	##this is actually a lot, A LOT, more complicated than I thought. 
	##gonna have to switch to ggplot2's paired boxplots, maybe with coord_flip()
	##das OK.
	
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
	
	##can we get a sense for gene drop-out in different sets?
	#sapply(subset_EXP_data, function(x) table(x$TYPE))
	##it's not awful. keep in mind that gene dropout can happen because of quality/coverage in some samples...
	
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
	
	##can we get a sense for gene drop-out in different sets?
	#sapply(subset_EXP_data, function(x) table(x$TYPE))
	##it's not awful. keep in mind that gene dropout can happen because of quality/coverage in some samples...
	
	##Now, gotta just label.
	for (x in 1:length(ALL_genesets)) {
		subset_EXP_data[[x]]$GO <- names(ALL_genesets)[x]
		
	}
	
	grand_EXP_frame <- rbindlist(subset_EXP_data)
	grand_EXP_frame$TYPE <- factor(grand_EXP_frame$TYPE, stats_names[order(stats_names)]) ##force alphabetical order in boxplots...
		
	##and now FINALLY...geom_boxplot
	
	##almost. One more term for grand_EXP_frame
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

