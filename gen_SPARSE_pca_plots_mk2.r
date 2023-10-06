####
##Generating a custom-go-term-loaded-PCA for Arinjay et al.
##
##

##I'm going to do this three ways.
##First, a 'global' transcriptional response as indicated in the tenOever paper Kushal sent over
##
##then separate PCAs with different GO terms loaded.
##GO.db should be appropriate for this?

##and I think also we're going to put in both species, nyeah?

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

##we're first doing a global response. So all genes...ever.
##ah. sparse pca, so we will have to come up with a global list of genes
##and sub in NAs where appropriate.

for (x in 1:length(bat_data)) {
    curr_na <- which(is.na(bat_data[[x]]$ORTHOFIND))
    for (ind in curr_na) {
        bat_data[[x]]$ORTHOFIND[ind] <- bat_data[[x]]$gene[ind]
    }
}

bat_genes <- unique(unlist(lapply(bat_data, function(x) x$ORTHOFIND)))
human_genes <- unique(unlist(lapply(human_dat, function(x) x$SYMBOL)))

all_genes <- union(bat_genes, human_genes)
all_genes <- all_genes[!is.na(all_genes)]

huge_matrix <- matrix(rep(rep(NA, length(all_genes)), length(c(human_names, bat_names))), nrow = length(all_genes))
huge_matrix <- as.data.frame(huge_matrix)
rownames(huge_matrix) <- all_genes
colnames(huge_matrix) <- c(human_names, bat_names)

for (x in 1:length(human_dat)) {
    curr_set <- human_dat[[x]]
    for (y in 1:dim(curr_set)[1]) {
        huge_matrix[curr_set$SYMBOL[y], human_names[x]] <- as.numeric(curr_set$log2FoldChange[y])
    }
    print(human_names[x])
}

for (x in 1:length(bat_data)) {
    curr_set <- bat_data[[x]]
    curr_set <- curr_set[!is.na(curr_set$ORTHOFIND),]
    for (y in 1:dim(curr_set)[1]) {
        huge_matrix[curr_set$ORTHOFIND[y], bat_names[x]] <- as.numeric(curr_set$log2FoldChange[y])
    }
    print(bat_names[x])
}

#neat.
saveRDS(huge_matrix, "logFC_matrix_global_for_PCA.rds")

library(PMA)

huge_matrix2 <- as.matrix(huge_matrix)

test <- apply(huge_matrix2, 1, function(x) mean(x, na.rm =T))

huge_matrix2 <- huge_matrix2[!is.na(test),]

curr_pca <- SPC(t(huge_matrix2), K = 2)

#plot(curr_pca$u[,1], curr_pca$u[,2])
plotframe <- as.data.frame(curr_pca$u)
colnames(plotframe) <- c("x", "y")
plotframe$sample <- colnames(huge_matrix2)

library(ggplot2)
thing <- ggplot(plotframe, aes(x = x, y = y, fill = sample)) + geom_point(size = 2, shape = 23) ##shape by species maybe? Or time-point.
##regardless, calu 12H really behaves fundamentally different to everything else.

try(dir.create("PCA"))
pdf("PCA/whole_transcriptome_PCA.pdf", width = 12, height = 12)
print(thing)
dev.off()

#####and next up.

##...geo terms.
#Fig 2D.: GO:0035457, GO:0035458, GO:0035455, GO:0035456, GO:0034340
##^^So just the ifn-1 response terms.

#Fig 2E: GO: 0005125, GO: 0008009
##So just the chemokine terms...

response_terms <- c("GO:0034097", "GO:0045087", "GO:0009615", "GO:0006954")

##chemokine/cytokine
chemokine <- c("GO:0005125", "GO:0008009")

##IFN-1 response
IFN <- c("GO:0035457","GO:0035458", "GO:0035455", "GO:0035456", "GO:0034340")

##Cell death
death <- c("GO:0008219")

##Leukocyte activation
#GO:0045431
##they got it wrong
activate <- c("GO:0045321")

##manually found this out-dated term
activate <- c("AZU1","CADM1","CD1D","CD2","CD24","CD276","CD28","CD3D","CD3E","CD4","CD40LG","CD47","CD7","CD79A","CEBPG","CLEC7A","CRTAM","CX3CL1","EBI3","ELF4","FOXP3","GLMN","HDAC4","HDAC5","HDAC7","HDAC9","HELLS","ICOSLG","IL12A","IL12B","IL18","IL2","IL21","IL27","IL4","IL7","CXCL8","CXCR2","INHA","INHBA","INS","JAG2","LAT","LAT2","LAX1","LCK","LST1","NCK1","NCK2","NFAM1","NHEJ1","NLRC3","PREX1","PRG3","PTPRC","SART1","SFTPD","SIRPG","SIT1","SLA2","SOCS5","SPACA3","SPINK5","TGFB1","THY1","TLR4","TNFSF13","TPD52","ZAP70")

##Now, let's do it...

library(org.Hs.eg.db)

ret <- AnnotationDbi::select(org.Hs.eg.db, keytype = "GOALL", keys = "GO:0045321", columns="SYMBOL")

library(msigdb)
library(msigdbr)
gsc = getMsigdb('hs', 'SYM')

all_gene_sets = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

response <- gen_sub_matrix(response_terms, huge_matrix2, "Immune Response")
#chemo <- gen_sub_matrix("CHEMOKINE", huge_matrix2, "Chemokine", dogrep = "CHEMOKINE")

####
##Update October 4th 2023
##turns out I didn't need to do all that grepping....just searched msigdb a bit more deeply...whoops.
chemokine_one <- c("ADIPOQ","AIMP1","ALKAL1","ALKAL2","AREG","BMP1","BMP10","BMP15","BMP2","BMP3","BMP4","BMP5","BMP6","BMP7","BMP8A","BMP8B","C17orf99","C1QTNF4","C5","CCL1","CCL11","CCL13","CCL14","CCL15","CCL16","CCL17","CCL18","CCL19","CCL2","CCL20","CCL21","CCL22","CCL23","CCL24","CCL25","CCL26","CCL27","CCL28","CCL3","CCL3L1","CCL3L3","CCL4","CCL5","CCL7","CCL8","CD40LG","CD70","CER1","CKLF","CLCF1","CMTM1","CMTM2","CMTM3","CMTM5","CMTM7","CMTM8","CNTF","CRLF1","CSF1","CSF2","CSF3","CX3CL1","CXCL1","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14","CXCL16","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8","CXCL9","EBI3","EDN1","EPO","FAM3B","FAM3C","FAM3D","FASLG","FGF2","FLT3LG","GDF1","GDF10","GDF11","GDF15","GDF2","GDF3","GDF5","GDF6","GDF7","GDF9","GPI","GPR15LG","GREM1","GREM2","GRN","HMGB1","IFNA1","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNB1","IFNE","IFNG","IFNK","IFNL1","IFNL2","IFNL3","IFNL4","IFNW1","IL10","IL11","IL12A","IL12B","IL13","IL15","IL16","IL17A","IL17B","IL17C","IL17D","IL17F","IL18","IL19","IL1A","IL1B","IL1F10","IL1RN","IL2","IL20","IL21","IL22","IL23A","IL24","IL25","IL26","IL27","IL3","IL31","IL32","IL33","IL34","IL36A","IL36B","IL36G","IL36RN","IL37","IL4","IL5","IL6","IL7","INHA","INHBA","INHBB","INHBC","INHBE","KITLG","LEFTY1","LEFTY2","LIF","LTA","LTB","MIF","MSMP","MSTN","NAMPT","NDP","NODAL","NRG1","OSM","PF4","PF4V1","PPBP","SCG2","SCGB3A1","SECTM1","SLURP1","SPP1","TAFA5","TGFB1","TGFB2","TGFB3","THNSL2","THPO","TIMP1","TNF","TNFRSF11B","TNFSF10","TNFSF11","TNFSF12","TNFSF13","TNFSF13B","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9","TSLP","VEGFA","VSTM1","WNT1","WNT10A","WNT10B","WNT11","WNT16","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B","WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT9A","WNT9B","XCL1","XCL2")
chemokine_two <- c("C5","CCL1","CCL11","CCL13","CCL14","CCL15","CCL16","CCL17","CCL18","CCL19","CCL2","CCL20","CCL21","CCL22","CCL23","CCL24","CCL25","CCL26","CCL27","CCL28","CCL3","CCL3L1","CCL3L3","CCL4","CCL5","CCL7","CCL8","CKLF","CX3CL1","CXCL1","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14","CXCL16","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8","CXCL9","GPR15LG","PF4","PF4V1","PPBP","XCL1","XCL2")

all_chemo <- unique(c(chemokine_one, chemokine_two))
chemo_res <- gen_sub_matrix("bla", huge_matrix2, "Chemokine", manual = all_chemo)
IFN_res <- gen_sub_matrix(IFN, huge_matrix2, "IFN-I Response")
death_res <- gen_sub_matrix(death, huge_matrix2, "Cell_Death", "CELL_DEATH")
activate_res <- gen_sub_matrix(activate, huge_matrix2, "Leukocyte_activation", manual = activate)

pdf("custom_PCA_plots_sars_bats_mk2_UPD.pdf", width = 12, height = 12)
print(response)
print(chemo_res)
print(IFN_res)
print(death_res)
print(activate_res)
dev.off()

gen_sub_matrix <- function(goterms, huge_matrix2, label, dogrep= NULL, manual = NULL) {
    curr_geneset <- all_gene_sets[all_gene_sets$gs_exact_source %in% goterms,]

    if(!is.null(dogrep)) {
        print('grep grep grep')
        curr_geneset <- all_gene_sets[grep(dogrep, all_gene_sets$gs_name),]
    }

    if(!is.null(manual)) {
        curr_geneset <- data.frame(gene_symbol = manual)
    }
    ##whatever. I'll just take it and call it a day.
    cut_matrix <- huge_matrix[rownames(huge_matrix) %in% curr_geneset$gene_symbol,]
    sub_pca <- SPC(t(cut_matrix), K = 2)

#plot(curr_pca$u[,1], curr_pca$u[,2])
plotframe <- as.data.frame(sub_pca$u)
colnames(plotframe) <- c("y", "x")
plotframe$sample <- colnames(huge_matrix2)

plotframe$cell <- factor(unlist(lapply(plotframe$sample, function(x) unlist(strsplit(x, "_"))[1])))
plotframe$time <- factor(unlist(lapply(plotframe$sample, function(x) unlist(strsplit(x, "_"))[length(unlist(strsplit(x, "_")))])))

library(ggplot2)
curr_thing <- ggplot(plotframe, aes(x = x, y = y, fill = cell, color = cell, shape = time, label = sample)) + geom_point(size = 6) ##shape by species maybe? Or time-point.
##regardless, calu 12H really behaves fundamentally different to everything else.

curr_thing <- curr_thing + geom_text() + theme_classic() + xlab(paste0("SPC1 (", round(sub_pca$prop.var.explained[2]*100, 2), "%)")) + ylab(paste0("SPC2 (", round(sub_pca$prop.var.explained[1]*100, 2), "%)"))

curr_thing <- curr_thing + ggtitle(label)

return(curr_thing)
}


response <- gen_sub_matrix_tiny(response_terms, huge_matrix2, "Immune Response")

chemo <- gen_sub_matrix_tiny("CHEMOKINE", huge_matrix2, "Chemokine", manual = all_chemo)

gen_sub_matrix_tiny <- function(goterms, huge_matrix2, label, dogrep= NULL, manual = NULL) {
    curr_geneset <- all_gene_sets[all_gene_sets$gs_exact_source %in% goterms,]

    if(!is.null(dogrep)) {
        print('grep grep grep')
        curr_geneset <- all_gene_sets[grep(dogrep, all_gene_sets$gs_name),]
    }

    if(!is.null(manual)) {
        curr_geneset <- data.frame(gene_symbol = manual)
    }
    ##whatever. I'll just take it and call it a day.
    cut_matrix <- huge_matrix[rownames(huge_matrix) %in% curr_geneset$gene_symbol,]
    if(!is.null(dogrep)) {
        print("used grep, return GO terms used")
        ret_list <- list()
        ret_list[[1]] <- cut_matrix
        ret_list[[2]] <- unique(curr_geneset$gs_exact_source)
        return(ret_list)
    }else{
    return(cut_matrix)
    }
}

###now let's make some heatmaps based on these....
setwd("PCA")
gen_heatmap(response, "GO_IFNI-response_terms", "GO IFN-I RESPONSE")
gen_heatmap(response, "GO_IFNI-response_terms_NOLU", "GO IFN-I RESPONSE", T)

##chemokine terms
gen_heatmap(chemo, "GO_Chemokine_terms", "GO CHEMOKINE")
gen_heatmap(chemo, "GO_Chemokine_terms_NOLU", "GO CHEMOKINE", T)

gen_heatmap <- function(log_frame, OUTFIX, plot_name, remove_LU = F) {
if(remove_LU) {
    print("removing LU!")
    log_frame <- log_frame[, !grepl("EF_LU", colnames(log_frame))]
}
sampleset <- colnames(log_frame)

na_counts <- apply(log_frame, 1, function(x) length(which(is.na(x))))
log_frame <- log_frame[na_counts < 4,]
##we'll allow half at most.

plotframe <- data.frame(sample = sampleset)

plotframe$cell <- factor(unlist(lapply(sampleset, function(x) unlist(strsplit(x, "_"))[1])))
plotframe$time <- factor(unlist(lapply(sampleset, function(x) unlist(strsplit(x, "_"))[length(unlist(strsplit(x, "_")))])))

plotframe$species <- "human"
plotframe$species[plotframe$cell %in% c("EFK3B", "EF")] <- "bat"

library(ComplexHeatmap)
#bat_data <- readLines("bat_human_orthofinder.txt")

#col_set_geno <- list(TREAT = c("mock" = "blue", "infected" = "yellow"))
#treat_col <- HeatmapAnnotation(TREAT = meta_subset_reorg$TREAT, col = col_set_geno, which = "column", show_annotation_name = FALSE)

col_set_SPEC <- list(SPECIES = c("human" = "lightblue", "bat" = "orange"))
spec_col <- HeatmapAnnotation(SPECIES = plotframe$species, col = col_set_SPEC, which = "column", show_annotation_name = FALSE)

col_set_cell_type <- list(CELL = c("CALU3" = "red", "A549" = "purple", "EFK3B" = "brown", "EF" = "gold"))
spec_cell <- HeatmapAnnotation(CELL = plotframe$cell, col = col_set_cell_type, which = "column", show_annotation_name = FALSE)

col_set_time <- list(TIME = c("12" = "pink", "24" = "magenta"))
spec_time <- HeatmapAnnotation(TIME = plotframe$time, col = col_set_time, which = "column", show_annotation_name = FALSE)

PROTEASE_thresh_map = Heatmap(log_frame, column_title= plot_name, name = "log2FC", #paste0(outfix, "_all_samples_all_genes_STAR_heatmap"),
	cluster_columns=FALSE, row_names_gp = gpar(fontsize = 6), show_row_names=F, show_column_names = F, top_annotation = c(spec_col, spec_cell, spec_time))

#log_frame_sums <- apply(log_frame, 1, function(x) length(which(x != 0)))

na_counts2 <- apply(log_frame, 1, function(x) length(which(is.na(x))))

#log_frame_sums_cut <- log_frame[log_frame_sums >= (dim(log_frame_sums)[2]*0.25),]
log_frame_sums_cut <- log_frame[na_counts2 < 3,]

PROTEASE_thresh_map_cut = Heatmap(log_frame_sums_cut, column_title= plot_name,  name = "log2FC", #paste0(outfix, "_all_samples_all_genes_STAR_heatmap"),
	cluster_columns=FALSE, row_names_gp = gpar(fontsize = 6), show_row_names=F, show_column_names = F, top_annotation = c(spec_col, spec_cell, spec_time))

pdf(paste0(OUTFIX, "_cross_species_heatmaps.pdf"), width = 12, height = 12)
print(PROTEASE_thresh_map)
print(PROTEASE_thresh_map_cut)
dev.off()

}


##^^it's that one 

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
