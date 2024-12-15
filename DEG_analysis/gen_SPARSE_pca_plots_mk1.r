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
##^^So just 

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

##Now, let's do it...

library(org.Hs.eg.db)

ret <- AnnotationDbi::select(org.Hs.eg.db, keytype = "GOALL", keys = "GO:0045321", columns="SYMBOL")

library(msigdb)
gsc = getMsigdb('hs', 'SYM')

all_gene_sets = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

response <- gen_sub_matrix(response_terms, huge_matrix2, "Immune Response")
chemo <- gen_sub_matrix("CHEMOKINE", huge_matrix2, "Chemokine", dogrep = "CHEMOKINE")
IFN_res <- gen_sub_matrix(IFN, huge_matrix2, "IFN-I Response")
death_res <- gen_sub_matrix(death, huge_matrix2, "Cell_Death", "CELL_DEATH")
activate_res <- gen_sub_matrix(activate, huge_matrix2, "Leukocyte_activation", "LEUKOCYTE_ACTIVATION")

pdf("custom_PCA_plots_sars_bats_mk1.pdf", width = 12, height = 12)
print(response)
print(chemo)
print(IFN_res)
print(death_res)
print(activate_res)
dev.off()

gen_sub_matrix <- function(goterms, huge_matrix2, label, dogrep= NULL) {
    curr_geneset <- all_gene_sets[all_gene_sets$gs_exact_source %in% goterms,]

    if(!is.null(dogrep)) {
        print('grep grep grep')
        curr_geneset <- all_gene_sets[grep(dogrep, all_gene_sets$gs_name),]
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
