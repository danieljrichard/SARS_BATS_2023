##############
##Comparing human and bat genes
##but instead of boxplots, it was requested for HEATMAPS
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
    curr_set <- curr_set[curr_set$log2FoldChange > log2(1.5),]
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
    bat_set <- bat_set[bat_set$log2FoldChange > log2(1.5),]
    bat_sigs[[x]] <- bat_set
    bat_data_forhuman[[x]] <- bat_data[[x]][!is.na(bat_data[[x]]$ORTHOFIND),]
 }

out_counts_bats <- table(unlist(lapply(bat_sigs, function(x) x$gene)))
##want at least in two conditions
out_count_two_bats <- out_counts_bats[out_counts_bats >= 2]
##oh....k.

massive_convert <- do.call("rbind", lapply(bat_data, function(x) x[, c("gene", "ORTHOFIND")]))
massive_convert <- massive_convert[!duplicated(massive_convert),]

converted_bat_counts <- unlist(lapply(names(out_counts_bats), function(x) massive_convert[massive_convert$gene == x, "ORTHOFIND"]))

shared_human_bat <- intersect(names(out_counts), converted_bat_counts)
human_unique <- setdiff(names(out_counts), converted_bat_counts)
##that's a lot of human unique.
human_unique_counts <- out_counts[human_unique]
human_unique_string <- names(human_unique_counts)[human_unique_counts > 1]
##...let's do it better.
human_unique_foldchanges <- unlist(lapply(human_unique, function(x) max(unlist(lapply(human_sigs, function(y) y[y$SYMBOL == x, "log2FoldChange"])))))
human_unique_nonzeroes <- unlist(lapply(human_unique, function(x) max(unlist(lapply(human_sigs, function(y) y[y$SYMBOL == x, "log2FoldChange"])))))
names(human_unique_foldchanges) <- human_unique
human_unique_foldchanges <- sort(human_unique_foldchanges, decreasing = T)
best_twenty_humans <- names(human_unique_foldchanges)[1:20]

human_unique_nonzeroes <- unlist(lapply(human_unique, function(x) length(which(!is.na(unlist(lapply(human_sigs, function(y) y[y$SYMBOL == x, "log2FoldChange"])))))))
names(human_unique_nonzeroes) <- human_unique
human_unique_nonzeroes <- sort(human_unique_nonzeroes, decreasing = T)
best_twenty_humans_nonzero <- names(human_unique_nonzeroes)[1:20]

bat_unique <- setdiff(converted_bat_counts, names(out_counts))
additional_bats <- names(out_counts_bats)[is.na(converted_bat_counts)]

total_bat_unique <- c(bat_unique, additional_bats)
total_bat_unique <- total_bat_unique[!is.na(total_bat_unique)]
##Let's pick top 20 for bats.
bat_unique_foldchanges <- unlist(lapply(total_bat_unique, function(x) max(unlist(lapply(bat_sigs, function(y) y[y$gene == x, "log2FoldChange"])))))
bat_unique_foldchanges <- as.numeric(bat_unique_foldchanges)
names(bat_unique_foldchanges) <- total_bat_unique
bat_unique_foldchanges <- bat_unique_foldchanges[!is.na(bat_unique_foldchanges)]
bat_unique_foldchanges <- sort(bat_unique_foldchanges, decreasing = T)
best_twenty_bats <- names(bat_unique_foldchanges)[1:20]

bat_unique_nonzeroes <- unlist(lapply(total_bat_unique, function(x) length(which(!is.na(unlist(lapply(bat_sigs, function(y) y[y$gene == x, "log2FoldChange"])))))))
bat_unique_nonzeroes <- as.numeric(bat_unique_nonzeroes)
names(bat_unique_nonzeroes) <- total_bat_unique
bat_unique_nonzeroes <- bat_unique_nonzeroes[!is.na(bat_unique_nonzeroes)]
bat_unique_nonzeroes <- sort(bat_unique_nonzeroes, decreasing = T)
best_twenty_bats_nozeroes <- names(bat_unique_nonzeroes)[1:20]


##now, for each part of the heatmap...

##start with SHARED...

shared_rowset <- list()

for (gene in shared_human_bat) {
    curr_bat_dat <- lapply(1:length(bat_data_forhuman), function(x) bat_data_forhuman[[x]][bat_data_forhuman[[x]]$ORTHOFIND == gene,])

    curr_human_dat <- lapply(1:length(human_dat), function(x) human_dat[[x]][human_dat[[x]]$SYMBOL == gene,])

    ##oh gods...what about multi-mappers?
    #...throw them out?
    size_set <- unlist(lapply(curr_bat_dat, function(x) dim(x)[1]))
    if(length(which(size_set > 1)) != 0) {
        print("multi-mapper")
        next()
    }else{
        print("will do one-to-one")
    }

    human_expression <- lapply(curr_human_dat, function(x) x[, 2:(grep("baseMean", colnames(x))-1)])
    bat_expression <- lapply(curr_bat_dat, function(x) x[, 2:(grep("baseMean", colnames(x))-1)])

    for (x in 1:length(bat_expression)) {
        if(dim(bat_expression[[x]])[1] == 0) {
            fake_frame <- data.frame(t(rep(0, dim(bat_expression[[x]])[2])))
            colnames(fake_frame) <- colnames(bat_expression[[x]])
            bat_expression[[x]] <- fake_frame
        }
        curr_scale <- data.frame(t(scale(as.numeric(bat_expression[[x]][1,]))))
        colnames(curr_scale) <- colnames(bat_expression[[x]])
        bat_expression[[x]] <- curr_scale
    }

    for(x in 1:length(human_expression)) {
        hum_curr_scale <- data.frame(t(scale(as.numeric(human_expression[[x]][1,]))))
        colnames(hum_curr_scale) <- colnames(human_expression[[x]])
        human_expression[[x]] <- hum_curr_scale
    }

    bat_row <- do.call("cbind", bat_expression)
    human_row <- do.call("cbind", human_expression)

    super_row <- cbind(human_row, bat_row)
   # super_row2 <- as.matrix(super_row)
   # super_row2[is.nan(super_row2)] <- 0
   # super_row <- as.data.frame(super_row2)
    colnames(super_row) <- unlist(lapply(colnames(super_row), function(x) meta_subset[meta_subset$simple == x, "NAME"]))
shared_rowset[[gene]] <- super_row
}

shared_rowset_frame <- do.call("rbind", shared_rowset)

meta_subset_reorg <- meta_subset[unlist(lapply(colnames(shared_rowset_frame), function(x) which(meta_subset$NAME == x))),]
library(ComplexHeatmap)
#bat_data <- readLines("bat_human_orthofinder.txt")

col_set_geno <- list(TREAT = c("mock" = "blue", "infected" = "yellow"))
treat_col <- HeatmapAnnotation(TREAT = meta_subset_reorg$TREAT, col = col_set_geno, which = "column", show_annotation_name = FALSE)

col_set_SPEC <- list(SPECIES = c("human" = "lightblue", "bat" = "orange"))
spec_col <- HeatmapAnnotation(SPECIES = meta_subset_reorg$species, col = col_set_SPEC, which = "column", show_annotation_name = FALSE)

col_set_cell_type <- list(CELL = c("calu3" = "red", "a549" = "purple", "efk3B" = "brown", "EF_LU" = "gold"))
spec_cell <- HeatmapAnnotation(CELL = meta_subset_reorg$cell_type, col = col_set_cell_type, which = "column", show_annotation_name = FALSE)

col_set_time <- list(TIME = c("12" = "pink", "24" = "magenta"))
spec_time <- HeatmapAnnotation(TIME = meta_subset_reorg$TIME, col = col_set_time, which = "column", show_annotation_name = FALSE)

SHARED_thresh_map = Heatmap(shared_rowset_frame, column_title= "Z-scored",# title = "Z-scored", #paste0(outfix, "_all_samples_all_genes_STAR_heatmap"),
	cluster_columns=FALSE, row_names_gp = gpar(fontsize = 6), show_row_names=T, show_column_names = F, top_annotation = c(treat_col,
		spec_col, spec_cell, spec_time))


############
############
############
############
##HUMAN SPECIFIC YALL

human_best_rowset <- list()

for (gene in best_twenty_humans) {
#for (gene in best_twenty_humans_nonzero) {
    curr_bat_dat <- lapply(1:length(bat_data_forhuman), function(x) bat_data_forhuman[[x]][bat_data_forhuman[[x]]$ORTHOFIND == gene,])

    curr_human_dat <- lapply(1:length(human_dat), function(x) human_dat[[x]][human_dat[[x]]$SYMBOL == gene,])

    ##oh gods...what about multi-mappers?
    #...throw them out?
    size_set <- unlist(lapply(curr_bat_dat, function(x) dim(x)[1]))
    if(length(which(size_set > 1)) != 0) {
        print("multi-mapper")
        next()
    }else{
        print("will do one-to-one")
    }

    human_expression <- lapply(curr_human_dat, function(x) x[, 2:(grep("baseMean", colnames(x))-1)])
    bat_expression <- lapply(curr_bat_dat, function(x) x[, 2:(grep("baseMean", colnames(x))-1)])

    for (x in 1:length(bat_expression)) {
        if(dim(bat_expression[[x]])[1] == 0) {
            fake_frame <- data.frame(t(rep(0, dim(bat_expression[[x]])[2])))
            colnames(fake_frame) <- colnames(bat_expression[[x]])
            bat_expression[[x]] <- fake_frame
        }
        curr_scale <- data.frame(t(scale(as.numeric(bat_expression[[x]][1,]))))
        colnames(curr_scale) <- colnames(bat_expression[[x]])
        bat_expression[[x]] <- curr_scale
    }

    for(x in 1:length(human_expression)) {
        hum_curr_scale <- data.frame(t(scale(as.numeric(human_expression[[x]][1,]))))
        colnames(hum_curr_scale) <- colnames(human_expression[[x]])
        human_expression[[x]] <- hum_curr_scale
    }

    bat_row <- do.call("cbind", bat_expression)
    human_row <- do.call("cbind", human_expression)

    super_row <- cbind(human_row, bat_row)
   # super_row2 <- as.matrix(super_row)
   # super_row2[is.nan(super_row2)] <- 0
   # super_row <- as.data.frame(super_row2)
    colnames(super_row) <- unlist(lapply(colnames(super_row), function(x) meta_subset[meta_subset$simple == x, "NAME"]))
human_best_rowset[[gene]] <- super_row
}

human_best_rowset_frame <- do.call("rbind", human_best_rowset)
human_best_rowset_frame[is.na(human_best_rowset_frame)] <- 0

human_best_thresh_map = Heatmap(human_best_rowset_frame, column_title= "Z-scored",# title = "Z-scored", #paste0(outfix, "_all_samples_all_genes_STAR_heatmap"),
	cluster_columns=FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize = 6), show_row_names=T, show_column_names = F, top_annotation = c(treat_col,
		spec_col, spec_cell, spec_time))

##########################
##########################
##########################
##########################
##and last but not least, the bat-specific

bat_best_rowset <- list()

for (gene in best_twenty_bats) {
#for (gene in best_twenty_bats_nozeroes) {
#    curr_bat_dat <- lapply(1:length(bat_data_forhuman), function(x) bat_data_forhuman[[x]][bat_data_forhuman[[x]]$gene == gene,])
    curr_bat_dat <- lapply(1:length(bat_data_forhuman), function(x) bat_data[[x]][bat_data[[x]]$gene == gene,])

    curr_human_dat <- lapply(1:length(human_dat), function(x) human_dat[[x]][human_dat[[x]]$SYMBOL == gene,])

    ##oh gods...what about multi-mappers?
    #...throw them out?
    size_set <- unlist(lapply(curr_bat_dat, function(x) dim(x)[1]))
    if(length(which(size_set > 1)) != 0) {
        print("multi-mapper")
        next()
    }else{
        print("will do one-to-one")
    }

    human_expression <- lapply(curr_human_dat, function(x) x[, 2:(grep("baseMean", colnames(x))-1)])
    bat_expression <- lapply(curr_bat_dat, function(x) x[, 2:(grep("baseMean", colnames(x))-1)])

    for (x in 1:length(bat_expression)) {
        if(dim(bat_expression[[x]])[1] == 0) {
            fake_frame <- data.frame(t(rep(0, dim(bat_expression[[x]])[2])))
            colnames(fake_frame) <- colnames(bat_expression[[x]])
            bat_expression[[x]] <- fake_frame
        }
        curr_scale <- data.frame(t(scale(as.numeric(bat_expression[[x]][1,]))))
        colnames(curr_scale) <- colnames(bat_expression[[x]])
        bat_expression[[x]] <- curr_scale
    }

    for(x in 1:length(human_expression)) {
        hum_curr_scale <- data.frame(t(scale(as.numeric(human_expression[[x]][1,]))))
        colnames(hum_curr_scale) <- colnames(human_expression[[x]])
        human_expression[[x]] <- hum_curr_scale
    }

    bat_row <- do.call("cbind", bat_expression)
    human_row <- do.call("cbind", human_expression)

    super_row <- cbind(human_row, bat_row)
   # super_row2 <- as.matrix(super_row)
   # super_row2[is.nan(super_row2)] <- 0
   # super_row <- as.data.frame(super_row2)
    colnames(super_row) <- unlist(lapply(colnames(super_row), function(x) meta_subset[meta_subset$simple == x, "NAME"]))
bat_best_rowset[[gene]] <- super_row
}

bat_best_rowset_frame <- do.call("rbind", bat_best_rowset)
bat_best_rowset_frame[is.na(bat_best_rowset_frame)] <- 0

bat_best_thresh_map = Heatmap(bat_best_rowset_frame, column_title= "Z-scored",# title = "Z-scored", #paste0(outfix, "_all_samples_all_genes_STAR_heatmap"),
	cluster_columns=FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize = 6), show_row_names=T, show_column_names = F, top_annotation = c(treat_col,
		spec_col, spec_cell, spec_time))

		
pdf("COMBINED_species_heatmaps_AUG2023.pdf", width = 12, height = 12)
pdf("COMBINED_species_heatmaps_AUG2023_NONZEROES.pdf", width = 12, height = 12)
print(SHARED_thresh_map)
print(human_best_thresh_map)
print(bat_best_thresh_map)
dev.off()

##pre-defined proteases:

proteases <- readLines("PROTEASE.txt")
proteases <- unique(proteases)

protease_rowset <- list()

for (gene in proteases) {
    curr_bat_dat <- lapply(1:length(bat_data_forhuman), function(x) bat_data_forhuman[[x]][bat_data_forhuman[[x]]$ORTHOFIND == gene,])

    curr_human_dat <- lapply(1:length(human_dat), function(x) human_dat[[x]][human_dat[[x]]$SYMBOL == gene,])

    ##oh gods...what about multi-mappers?
    #...throw them out?
    size_set <- unlist(lapply(curr_bat_dat, function(x) dim(x)[1]))
    if(length(which(size_set > 1)) != 0) {
        print("multi-mapper")
        next()
    }else{
        print("will do one-to-one")
    }

    human_expression <- lapply(curr_human_dat, function(x) x[, 2:(grep("baseMean", colnames(x))-1)])
    bat_expression <- lapply(curr_bat_dat, function(x) x[, 2:(grep("baseMean", colnames(x))-1)])

    for (x in 1:length(bat_expression)) {
        if(dim(bat_expression[[x]])[1] == 0) {
            fake_frame <- data.frame(t(rep(0, dim(bat_expression[[x]])[2])))
            colnames(fake_frame) <- colnames(bat_expression[[x]])
            bat_expression[[x]] <- fake_frame
        }
        curr_scale <- data.frame(t(scale(as.numeric(bat_expression[[x]][1,]))))
        colnames(curr_scale) <- colnames(bat_expression[[x]])
        bat_expression[[x]] <- curr_scale
    }

    for(x in 1:length(human_expression)) {
        hum_curr_scale <- data.frame(t(scale(as.numeric(human_expression[[x]][1,]))))
        colnames(hum_curr_scale) <- colnames(human_expression[[x]])
        human_expression[[x]] <- hum_curr_scale
    }

    bat_row <- do.call("cbind", bat_expression)
    human_row <- do.call("cbind", human_expression)

    super_row <- cbind(human_row, bat_row)
   # super_row2 <- as.matrix(super_row)
   # super_row2[is.nan(super_row2)] <- 0
   # super_row <- as.data.frame(super_row2)
    colnames(super_row) <- unlist(lapply(colnames(super_row), function(x) meta_subset[meta_subset$simple == x, "NAME"]))
protease_rowset[[gene]] <- super_row
}

protease_rowset_frame <- do.call("rbind", protease_rowset)
protease_rowset_frame[is.na(protease_rowset_frame)] <- 0

meta_subset_reorg <- meta_subset[unlist(lapply(colnames(protease_rowset_frame), function(x) which(meta_subset$NAME == x))),]
library(ComplexHeatmap)
#bat_data <- readLines("bat_human_orthofinder.txt")

col_set_geno <- list(TREAT = c("mock" = "blue", "infected" = "yellow"))
treat_col <- HeatmapAnnotation(TREAT = meta_subset_reorg$TREAT, col = col_set_geno, which = "column", show_annotation_name = FALSE)

col_set_SPEC <- list(SPECIES = c("human" = "lightblue", "bat" = "orange"))
spec_col <- HeatmapAnnotation(SPECIES = meta_subset_reorg$species, col = col_set_SPEC, which = "column", show_annotation_name = FALSE)

col_set_cell_type <- list(CELL = c("calu3" = "red", "a549" = "purple", "efk3B" = "brown", "EF_LU" = "gold"))
spec_cell <- HeatmapAnnotation(CELL = meta_subset_reorg$cell_type, col = col_set_cell_type, which = "column", show_annotation_name = FALSE)

col_set_time <- list(TIME = c("12" = "pink", "24" = "magenta"))
spec_time <- HeatmapAnnotation(TIME = meta_subset_reorg$TIME, col = col_set_time, which = "column", show_annotation_name = FALSE)

PROTEASE_thresh_map = Heatmap(protease_rowset_frame, column_title= "Z-scored",# title = "Z-scored", #paste0(outfix, "_all_samples_all_genes_STAR_heatmap"),
	cluster_columns=FALSE, row_names_gp = gpar(fontsize = 6), show_row_names=T, show_column_names = F, top_annotation = c(treat_col,
		spec_col, spec_cell, spec_time))

protease_rowset_frame_sums <- apply(protease_rowset_frame, 1, function(x) length(which(x != 0)))

protease_rowset_frame_cut <- protease_rowset_frame[protease_rowset_frame_sums >= (dim(protease_rowset_frame)[2]*0.25),]

PROTEASE_thresh_map_cut = Heatmap(protease_rowset_frame_cut, column_title= "Z-scored",# title = "Z-scored", #paste0(outfix, "_all_samples_all_genes_STAR_heatmap"),
	cluster_columns=FALSE, row_names_gp = gpar(fontsize = 6), show_row_names=T, show_column_names = F, top_annotation = c(treat_col,
		spec_col, spec_cell, spec_time))
pdf("COMBINED_species_heatmaps_AUG2023_PROTEASES.pdf", width = 12, height = 12)
print(PROTEASE_thresh_map)
print(PROTEASE_thresh_map_cut)
dev.off()


#####################################
#####################################
#####################################
#####################################
###
##Now working with customized GO terms


proteases <- readLines("PROTEASE.txt")
proteases <- unique(proteases)

protease_rowset <- list()

for (gene in proteases) {
    curr_bat_dat <- lapply(1:length(bat_data_forhuman), function(x) bat_data_forhuman[[x]][bat_data_forhuman[[x]]$ORTHOFIND == gene,])

    curr_human_dat <- lapply(1:length(human_dat), function(x) human_dat[[x]][human_dat[[x]]$SYMBOL == gene,])

    ##oh gods...what about multi-mappers?
    #...throw them out?
    size_set <- unlist(lapply(curr_bat_dat, function(x) dim(x)[1]))
    if(length(which(size_set > 1)) != 0) {
        print("multi-mapper")
        next()
    }else{
        print("will do one-to-one")
    }

    human_expression <- lapply(curr_human_dat, function(x) x[, 2:(grep("baseMean", colnames(x))-1)])
    bat_expression <- lapply(curr_bat_dat, function(x) x[, 2:(grep("baseMean", colnames(x))-1)])

    for (x in 1:length(bat_expression)) {
        if(dim(bat_expression[[x]])[1] == 0) {
            fake_frame <- data.frame(t(rep(0, dim(bat_expression[[x]])[2])))
            colnames(fake_frame) <- colnames(bat_expression[[x]])
            bat_expression[[x]] <- fake_frame
        }
        curr_scale <- data.frame(t(scale(as.numeric(bat_expression[[x]][1,]))))
        colnames(curr_scale) <- colnames(bat_expression[[x]])
        bat_expression[[x]] <- curr_scale
    }

    for(x in 1:length(human_expression)) {
        hum_curr_scale <- data.frame(t(scale(as.numeric(human_expression[[x]][1,]))))
        colnames(hum_curr_scale) <- colnames(human_expression[[x]])
        human_expression[[x]] <- hum_curr_scale
    }

    bat_row <- do.call("cbind", bat_expression)
    human_row <- do.call("cbind", human_expression)

    super_row <- cbind(human_row, bat_row)
   # super_row2 <- as.matrix(super_row)
   # super_row2[is.nan(super_row2)] <- 0
   # super_row <- as.data.frame(super_row2)
    colnames(super_row) <- unlist(lapply(colnames(super_row), function(x) meta_subset[meta_subset$simple == x, "NAME"]))
protease_rowset[[gene]] <- super_row
}

protease_rowset_frame <- do.call("rbind", protease_rowset)
protease_rowset_frame[is.na(protease_rowset_frame)] <- 0

meta_subset_reorg <- meta_subset[unlist(lapply(colnames(protease_rowset_frame), function(x) which(meta_subset$NAME == x))),]
library(ComplexHeatmap)
#bat_data <- readLines("bat_human_orthofinder.txt")

col_set_geno <- list(TREAT = c("mock" = "blue", "infected" = "yellow"))
treat_col <- HeatmapAnnotation(TREAT = meta_subset_reorg$TREAT, col = col_set_geno, which = "column", show_annotation_name = FALSE)

col_set_SPEC <- list(SPECIES = c("human" = "lightblue", "bat" = "orange"))
spec_col <- HeatmapAnnotation(SPECIES = meta_subset_reorg$species, col = col_set_SPEC, which = "column", show_annotation_name = FALSE)

col_set_cell_type <- list(CELL = c("calu3" = "red", "a549" = "purple", "efk3B" = "brown", "EF_LU" = "gold"))
spec_cell <- HeatmapAnnotation(CELL = meta_subset_reorg$cell_type, col = col_set_cell_type, which = "column", show_annotation_name = FALSE)

col_set_time <- list(TIME = c("12" = "pink", "24" = "magenta"))
spec_time <- HeatmapAnnotation(TIME = meta_subset_reorg$TIME, col = col_set_time, which = "column", show_annotation_name = FALSE)

PROTEASE_thresh_map = Heatmap(protease_rowset_frame, column_title= "Z-scored",# title = "Z-scored", #paste0(outfix, "_all_samples_all_genes_STAR_heatmap"),
	cluster_columns=FALSE, row_names_gp = gpar(fontsize = 6), show_row_names=T, show_column_names = F, top_annotation = c(treat_col,
		spec_col, spec_cell, spec_time))

protease_rowset_frame_sums <- apply(protease_rowset_frame, 1, function(x) length(which(x != 0)))

protease_rowset_frame_cut <- protease_rowset_frame[protease_rowset_frame_sums >= (dim(protease_rowset_frame)[2]*0.25),]

PROTEASE_thresh_map_cut = Heatmap(protease_rowset_frame_cut, column_title= "Z-scored",# title = "Z-scored", #paste0(outfix, "_all_samples_all_genes_STAR_heatmap"),
	cluster_columns=FALSE, row_names_gp = gpar(fontsize = 6), show_row_names=T, show_column_names = F, top_annotation = c(treat_col,
		spec_col, spec_cell, spec_time))
pdf("COMBINED_species_heatmaps_AUG2023_PROTEASES.pdf", width = 12, height = 12)
print(PROTEASE_thresh_map)
print(PROTEASE_thresh_map_cut)
dev.off()
