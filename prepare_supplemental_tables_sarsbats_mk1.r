##############
##Generating outputs for Supplemental Table S1
##April 16th 2024
###

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

bat_names <- unlist(lapply(bat_files, function(x) unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))]))
bat_names <- gsub("_NORMALIZED_counts_matrix_ALLGENES_STATS_orthofinder_human_gene_definitions.csv", "", bat_names)
bat_names <- gsub("BAT_SARS_infection_", "", bat_names)

human_names <- unlist(lapply(human_files, function(x) unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))]))
human_names <- gsub("_NORMALIZED_counts_matrix_ALLGENES_STATS_nameconvert.csv", "", human_names)
human_names <- gsub("human_sars_infection_", "", human_names)

for (x in 1:length(human_dat)) {
    human_set <- human_dat[[x]]
    human_set <- human_set[!is.na(human_set$padj),]
    human_set$TREATMENT <- human_names[x]
    human_dat[[x]] <- human_set
    print(x)
}

for (x in 1:length(bat_data)) {
     bat_set <- bat_data[[x]]
    bat_set <- bat_set[!is.na(bat_set$padj),]
    for (a in which(is.na(bat_set$ORTHOFIND))) {
        bat_set$ORTHOFIND[a] <- bat_set$gene[a]
    }
    bat_set$TREATMENT <- bat_names[x]
    bat_data[[x]] <- bat_set
}

library(writexl)
names(human_dat) <- human_names
names(bat_data) <- bat_names

super_list <- do.call("c", list(human_dat, bat_data))

try(dir.create("SUPPLEMENTAL_TABLES"))
write_xlsx(super_list, "SUPPLEMENTAL_TABLES/Supp_table_S1_degs.xlsx")

saveRDS(super_list, "ALL_DEG_res_April2024.rds")

#####
##Update April 16th 2024
##going to also output the sets of genes associated with TOPGO enriched GO terms

all_res <- readRDS("/extra/SARS_BATS/ALL_DEG_res_April2024.rds")

setwd("/extra/SARS_BATS/BAT_ALIGN/FINAL/BAT_BP_ENRICHMENTS_SLIM")

res <- system("ls *SLICED.csv", intern = T)

bat_names <- unlist(lapply(res, function(x) unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))]))
bat_names <- gsub("_TOPGO_enrichments_FDR0.05_SLICED.csv", "", bat_names)
bat_names <- gsub("BAT_SARS_infection_", "", bat_names)

bat_data <- lapply(res, read.csv)
global_map <- "/extra/SARS_BATS/CUSTOM_GO_ENRICHMENTS/BAT/BAT_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"

##need to go through each result and input the genes associated.
##being careful to watch directionality.
library(data.table)
global_dat <- as.data.frame(fread(global_map, header = F))
###looks like every gene with every possible GO term. Cool.

bat_output <- list()

for (x in 1:length(bat_data)) {
    bat_res <- bat_data[[x]]
    if(dim(bat_res)[1] ==  0) {
        print("no significant terms")
        next()
    }
    bat_file <- bat_names[x]
    bat_file_short <- paste0(unlist(strsplit(bat_file, "_"))[1:(length(unlist(strsplit(bat_file, "_"))) - 1)], collapse = "_")
    bat_deg <- all_res[[bat_file_short]]
    bat_deg$log2FoldChange <- as.numeric(bat_deg$log2FoldChange)

    direction <- tail(unlist(strsplit(bat_file, "_")), 1)

    if(direction == "UPREG") {
        bat_deg <- bat_deg[bat_deg$log2FoldChange > 0,]
        }else{
        bat_deg <- bat_deg[bat_deg$log2FoldChange < 0,]
    }
    bat_deg <- bat_deg[!is.na(bat_deg$padj),]
    bat_deg <- bat_deg[bat_deg$padj < 0.05,]
    ##for each term, output collapsed set of genes.
    curr_sub_mat <- global_dat[global_dat[,1] %in% bat_deg$gene,]
    ##now everything in here will be relevant.
    go_collapses <- unlist(lapply(bat_res[, "GO.ID"], function(x) paste0(curr_sub_mat[which(grepl(x, curr_sub_mat[,2])), 1], collapse = ",")))
    bat_res <- bat_res[, c("GO.ID", "Term", "Annotated", "Significant", "Expected", "topgoFisher", "PADJ")]
    bat_res <- data.frame(TREATMENT = bat_file, bat_res)
    bat_res$HIT_GENES <- go_collapses
    bat_output[[bat_file]] <- bat_res
    print(bat_file)
}

saveRDS(bat_output, "SUPP_TABLE_2_BAT_GO_res.rds")

library(writexl)
write_xlsx(bat_output, "SUPP_TABLE_2_BAT_GO_res.xlsx")


all_res <- readRDS("/extra/SARS_BATS/ALL_DEG_res_April2024.rds")

setwd("/extra/SARS_BATS/HUMAN_ALIGN/FINAL_TEST/HUMAN_BP_ENRICHMENTS_SLIM")

res <- system("ls *SLICED.csv", intern = T)

human_names <- unlist(lapply(res, function(x) unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))]))
human_names <- gsub("_TOPGO_enrichments_FDR0.05_SLICED.csv", "", human_names)
human_names <- gsub("human_sars_infection_", "", human_names)

global_map <- "/extra/SARS_BATS/CUSTOM_GO_ENRICHMENTS/HUMAN/HUMAN_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"

bat_data <- lapply(res, read.csv)

##need to go through each result and input the genes associated.
##being careful to watch directionality.
library(data.table)
global_dat <- as.data.frame(fread(global_map, header = F))
###looks like every gene with every possible GO term. Cool.

bat_output <- list()

for (x in 1:length(bat_data)) {
    bat_res <- bat_data[[x]]
    if(dim(bat_res)[1] ==  0) {
        print("no significant terms")
        next()
    }
    bat_file <- human_names[x]
    bat_file_short <- paste0(unlist(strsplit(bat_file, "_"))[1:(length(unlist(strsplit(bat_file, "_"))) - 1)], collapse = "_")
    bat_deg <- all_res[[bat_file_short]]
    bat_deg$log2FoldChange <- as.numeric(bat_deg$log2FoldChange)

    direction <- tail(unlist(strsplit(bat_file, "_")), 1)

    if(direction == "UPREG") {
        bat_deg <- bat_deg[bat_deg$log2FoldChange > 0,]
        }else{
        bat_deg <- bat_deg[bat_deg$log2FoldChange < 0,]
    }
    bat_deg <- bat_deg[!is.na(bat_deg$padj),]
    bat_deg <- bat_deg[bat_deg$padj < 0.05,]
    ##for each term, output collapsed set of genes.
    curr_sub_mat <- global_dat[global_dat[,1] %in% bat_deg$gene,]
    curr_sub_mat$SYMBOL <- unlist(lapply(curr_sub_mat[,1], function(x) bat_deg$SYMBOL[which(bat_deg$gene == x)[1]]))
    ##now everything in here will be relevant.
    go_collapses <- unlist(lapply(bat_res[, "GO.ID"], function(x) paste0(curr_sub_mat[which(grepl(x, curr_sub_mat[,2])), "SYMBOL"], collapse = ",")))
    bat_res <- bat_res[, c("GO.ID", "Term", "Annotated", "Significant", "Expected", "topgoFisher", "PADJ")]
    bat_res <- data.frame(TREATMENT = bat_file, bat_res)
    bat_res$HIT_GENES <- go_collapses
    bat_output[[bat_file]] <- bat_res
    print(bat_file)
}

saveRDS(bat_output, "SUPP_TABLE_2_HUMAN_GO_res.rds")

library(writexl)
write_xlsx(bat_output, "SUPP_TABLE_2_HUMAN_GO_res.xlsx")

####

bat <- readRDS("./BAT_ALIGN/FINAL/BAT_BP_ENRICHMENTS_SLIM/SUPP_TABLE_2_BAT_GO_res.rds")
human <- readRDS("./HUMAN_ALIGN/FINAL_TEST/HUMAN_BP_ENRICHMENTS_SLIM/SUPP_TABLE_2_HUMAN_GO_res.rds")
library(writexl)
super_list <- do.call("c", list(human, bat))
write_xlsx(super_list, "SUPPLEMENTAL_TABLES/TABLE_S2_human_bat_GO_enrich.xlsx")