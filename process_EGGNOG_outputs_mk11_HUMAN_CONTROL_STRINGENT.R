####################
##Processing EGGNOG OUTPUT.
##

eggnog <- "query_seqs.fa.emapper.annotations"
CDS_file <- "GCF_000308155.1_EptFus1.0_cds_from_genomic.fna"

transcript_file <- "bat_sars_infection_trimgalore_pass2_STAR_RAW_counts_matrix_genenames.csv"

outfix <- "BAT_EGGNOG"


eggnog <- "query_seqs.fa.emapper.annotations"
CDS_file <- "Homo_sapiens.GRCh37.cds.all.fa"
##human annotations have gene names right in the peptide file.
CDS_file <- "Homo_sapiens.GRCh37.pep.all.fa"
transcript_file <- "human_sars_infection_trimgalore_pass2_STAR_RAW_counts_matrix_genenames.csv"
outfix <- "HUMAN_EGGNOG"
human <- T

collapse_gene_GO_annotations <- function(eggnog, CDS_file, transcript_file, outfix, human = F) {
	library(data.table)
	egg <- fread(eggnog)
	relevant_egg <- as.data.frame(egg[, c("#query_name", "seed_eggNOG_ortholog", "seed_ortholog_evalue", "seed_ortholog_score", "Preferred_name", "GOs", "V21", "V22")])
	fwrite(relevant_egg, paste0(eggnog, ".cleanup"), row.names = F, sep = "\t")
	
	if(length(which(duplicated(relevant_egg[,1]))) != 0) {
		print("no duplicates - NICE!")
	}
	
	##next step: marrying backup with the initial genes...
	
	trans_dat <- read.csv(transcript_file)
	
	trans_genes <- as.character(trans_dat$gene)
	
	if(human) {
		pept_grab <- system(paste0("grep -F -f - ", CDS_file), input = paste0(">", relevant_egg[,1]), intern=T)
		
		if(length(pept_grab) != length(unique(relevant_egg[,1]))) {
			print("...how unfortunate!")
			browser()
		}
		super_genes <- sapply(pept_grab, function(x) unlist(strsplit(x, "gene:", fixed = T))[2])
		super_genes <- sapply(super_genes, function(x) unlist(strsplit(x, " trans", fixed = T))[1])
		
		##let's just...do the following
		
		super_prots <- sapply(pept_grab, function(x) unlist(strsplit(x, ">", fixed = T))[2])
		super_prots <- sapply(super_prots, function(x) unlist(strsplit(x, " pep", fixed = T))[1])
		
		lookup_table <- data.frame(gene = super_genes, protein_ID = super_prots)
		colnames(relevant_egg)[1] <- "protein_ID"
		lookup_SUPER <- merge(lookup_table, relevant_egg, by = "protein_ID")
		
		##...it's all just...formatting.
		
		fwrite(lookup_SUPER, paste0(outfix, ".GENEID_CDS_MAPPING"), row.names = F, sep = "\t")
	}else if(!is.null(CDS_file)) {
	
	pept_grab <- system(paste0("grep -F -f - ", CDS_file), input = paste0("protein_id=", relevant_egg[,1]), intern=T)
		
		if(length(pept_grab) != length(unique(relevant_egg[,1]))) {
			print("...how unfortunate!")
			browser()
		}
		##I'm just trying to make this as...seamless as possible.
		
		super_genes <- sapply(pept_grab, function(x) unlist(strsplit(x, "[gene=", fixed = T))[2])
		super_genes <- sapply(super_genes, function(x) unlist(strsplit(x, "]", fixed = T))[1])
		
		##let's just...do the following
		
		super_prots <- sapply(pept_grab, function(x) unlist(strsplit(x, "[protein_id=", fixed = T))[2])
		super_prots <- sapply(super_prots, function(x) unlist(strsplit(x, "]", fixed = T))[1])
		
		lookup_table <- data.frame(gene = super_genes, protein_ID = super_prots)
		colnames(relevant_egg)[1] <- "protein_ID"
		lookup_SUPER <- merge(lookup_table, relevant_egg, by = "protein_ID")
		
		##...it's all just...formatting.
		
		fwrite(lookup_SUPER, paste0(outfix, ".GENEID_CDS_MAPPING"), row.names = F, sep = "\t")
	}else{
		
		print("human") ##might just use the human CDS for streamlining.
		browser()
	}
	
	if(human) {
		lookup_SUPER$gene <- sapply(lookup_SUPER$gene, function(x) unlist(strsplit(x, ".", fixed=T))[1])
	}
	
	subset_SUPER <- subset(lookup_SUPER, lookup_SUPER$gene %in% trans_genes)
	
	frame <- subset(subset_SUPER, subset_SUPER$gene %in% "WNK1")
	frame <- subset(subset_SUPER, subset_SUPER$gene %in% "TRIP6")
	
	collapse_gene_cat <- function(frame) {
		go_set <- unique(unlist(lapply(frame$GOs, function(x) unlist(strsplit(x, ",", fixed=T))))) ##confirmed with table, this looks good.
		if(length(go_set) == 0) {
			go_set <- "NA"
		}
		curr_gene_row <- data.frame(gene = frame$gene[1], GO = paste0(go_set, collapse = ","))
		return(curr_gene_row)
	}
	GENES_SUPER <- lapply(unique(subset_SUPER$gene), function(x) collapse_gene_cat(subset(subset_SUPER, subset_SUPER$gene %in% x)))
	
	GENE_SUPER_FRAME <- as.data.frame(rbindlist(GENES_SUPER))
	super_clean <- subset_SUPER[!duplicated(subset_SUPER$gene), c("gene", "Preferred_name", "V22")]
	
	informative_SUPER_GENE <- merge(GENE_SUPER_FRAME, super_clean, by = "gene")
	informative_SUPER_GENE$GO <- as.character(informative_SUPER_GENE$GO)
	
	fwrite(informative_SUPER_GENE, paste0(outfix, ".TRANSCRIPT_ALIGNED_USEFUL.output"), sep = "\t", row.names = F, quote = F)
	
	gene_lines <- paste0(informative_SUPER_GENE$gene, "\t", informative_SUPER_GENE$GO)
#	gene_lines <- gsub(",", "\t", gene_lines)
	writeLines(gene_lines, paste0(outfix, "_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"))
		
	if(length(which(duplicated(informative_SUPER_GENE$gene))) != 0) {
		print("oh snappo")
		browser()
	}
	
	clean_informative_SUPER_GENE <- informative_SUPER_GENE[informative_SUPER_GENE$GO != "NA",]
	GO_freq_table <- table(unlist(lapply(clean_informative_SUPER_GENE$GO, function(x) unlist(strsplit(x, ",", fixed=T)))))
	png(paste0(outfix, "_GO_term_freq_histogram.png"), width = 800, height = 800, units = "px")
	hist(GO_freq_table, breaks = 200)
	dev.off()
	
	##let's just...write this to file?
	GO_freq_out <- data.frame(GO_TERM = names(GO_freq_table), gene_counts = as.numeric(GO_freq_table))
	fwrite(GO_freq_out, paste0(outfix, "_TRANSCRIPT_ALIGNED_GO_TERM_FREQUENCY_BACKGROUNDS.tsv"), sep = "\t", row.names = F, quote = F)
	}

##########

global_map <- "BAT_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
targets <- "24H_infection_de_set.txt"
outfix <- "EFK3B_24H_TEST"

targets <- "24H_UPREG_EFK3B_TEST.txt"
global_map <- "BAT_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
outfix <- "EFK3B_24H_UPREG"
cutoff = 0.05

targets <- "24H_DOWNREG_EFK3B_TEST.txt"
global_map <- "BAT_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
outfix <- "EFK3B_24H_DOWNREG"
cutoff = 0.05

targets <- "HUMAN_CALU3_24H_INF_UPREG_genes.txt"
global_map <- "HUMAN_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
outfix <- "HUMAN_CALU3_24H_UPREG"
cutoff = 0.05

targets <- "HUMAN_CALU3_24H_INF_DOWNREG_genes.txt"
global_map <- "HUMAN_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
outfix <- "HUMAN_CALU3_24H_DOWNREG"
cutoff = 0.05
test_GO_with_TOPGO(targets, global_map, outfix, cutoff)

human_map <- "/extra/SARS_BATS/CUSTOM_GO_ENRICHMENTS/HUMAN/HUMAN_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"

human_wrapper <- function() {
	human_files <- readLines("HUMAN_sig_hits.txt")
	try(dir.create("HUMAN_BP_ENRICHMENTS_SLIM_STRINGENT"))
	for (file in human_files) {
		end <- unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
		prefix <- unlist(strsplit(end, "_NORMALIZED"))[1]
		curr_dat <- read.csv(file)
		curr_sig <- curr_dat[curr_dat$padj < 0.05,] ##just in case.
		try(unlink(paste0("/thing/temp_UPREG_HUMAN_", prefix)))
		try(unlink(paste0("/thing/temp_DOWNREG_HUMAN_", prefix)))
		writeLines(curr_sig$gene[curr_sig$log2FoldChange > 1], paste0("/thing/temp_UPREG_HUMAN_", prefix))
		writeLines(curr_sig$gene[curr_sig$log2FoldChange < -1], paste0("/thing/temp_DOWNREG_HUMAN_", prefix))
		
		cutoff <- 0.05
		
		if(file.exists(paste0("HUMAN_BP_ENRICHMENTS_SLIM_STRINGENT/", prefix, "_UPREG_LOGFC1", "_TOPGO_enrichments_UNADJUSTED", cutoff, ".csv"))) {
			print("skip")
		}else{
            if(length(curr_sig$gene[curr_sig$log2FoldChange > 1]) < 5) {
                print("not doing")
                next()
            }
		test_GO_with_TOPGO(paste0("/thing/temp_UPREG_HUMAN_", prefix), human_map, paste0("HUMAN_BP_ENRICHMENTS_SLIM_STRINGENT/", prefix, "_UPREG_LOGFC1"), 0.05)
		}
		print(paste0(prefix, "_UPREG"))
		if(file.exists(paste0("HUMAN_BP_ENRICHMENTS_SLIM_STRINGENT/", prefix, "_DOWNREG_LOGFC1", "_TOPGO_enrichments_UNADJUSTED", cutoff, ".csv"))) {
			print("skip")
		}else{
            if(length(curr_sig$gene[curr_sig$log2FoldChange < -1]) < 5) {
                print("not doing")
                next()
            }
			test_GO_with_TOPGO(paste0("/thing/temp_DOWNREG_HUMAN_", prefix), human_map, paste0("HUMAN_BP_ENRICHMENTS_SLIM_STRINGENT/", prefix, "_DOWNREG_LOGFC1"), 0.05)
		}
		print(paste0(prefix, "_DOWNREG"))
		
		print("additional GO validation")
		
if(FALSE) {

		if(file.exists(paste0("HUMAN_BP_ENRICHMENTS_SLIM/", prefix, "_UPREG_VALIDATE", "_TOPGO_enrichments_UNADJUSTED", cutoff, ".csv"))) {
			print("skip")
		}else{
		test_GO_with_TOPGO_ANNOTATION_DATABASE("/thing/temp_UPREG_HUMAN", human_map, paste0("HUMAN_BP_ENRICHMENTS_SLIM/", prefix, "_UPREG_VALIDATE"), 0.05)
		}
		print(paste0(prefix, "_UPREG"))
		if(file.exists(paste0("HUMAN_BP_ENRICHMENTS_SLIM/", prefix, "_DOWNREG_VALIDATE", "_TOPGO_enrichments_UNADJUSTED", cutoff, ".csv"))) {
			print("skip")
		}else{
		test_GO_with_TOPGO_ANNOTATION_DATABASE("/thing/temp_DOWNREG_HUMAN", human_map, paste0("HUMAN_BP_ENRICHMENTS_SLIM/", prefix, "_DOWNREG_VALIDATE"), 0.05)
		}
		print(paste0(prefix, "_DOWNREG"))
	}
    }
	
}


bat_map <- "/extra/SARS_BATS/CUSTOM_GO_ENRICHMENTS/BAT/BAT_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"

bat_wrapper <- function() {
	bat_files <- readLines("BAT_sig_hits.txt")
    bat_files <- bat_files[!grepl("EF_LU", bat_files)]
	try(dir.create("BAT_BP_ENRICHMENTS_SLIM_STRINGENT"))
	for (file in bat_files) {
		end <- unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
		prefix <- unlist(strsplit(end, "_NORMALIZED"))[1]
		curr_dat <- read.csv(file)
		curr_sig <- curr_dat[curr_dat$padj < 0.05,] ##just in case.
		try(unlink(paste0("/thing/temp_UPREG_BAT_", prefix)))
		try(unlink(paste0("/thing/temp_DOWNREG_BAT_", prefix)))
		writeLines(curr_sig$gene[curr_sig$log2FoldChange > 1], paste0("/thing/temp_UPREG_BAT_", prefix))
		writeLines(curr_sig$gene[curr_sig$log2FoldChange < -1], paste0("/thing/temp_DOWNREG_BAT_", prefix))
   		cutoff <- 0.05
		if(file.exists(paste0("BAT_BP_ENRICHMENTS_SLIM_STRINGENT/", prefix, "_UPREG_LOGFC1", "_TOPGO_enrichments_UNADJUSTED", cutoff, ".csv"))) {
			print("skip")
		}else{
            if(length(curr_sig$gene[curr_sig$log2FoldChange > 1]) < 5) {
                print("not doing")
                next()
            }
		test_GO_with_TOPGO(paste0("/thing/temp_UPREG_BAT_", prefix), bat_map, paste0("BAT_BP_ENRICHMENTS_SLIM_STRINGENT/", prefix, "_UPREG_LOGFC1"), 0.05)
		}
		
		print(paste0(prefix, "_UPREG"))
		if(file.exists(paste0("BAT_BP_ENRICHMENTS_SLIM_STRINGENT/", prefix, "_DOWNREG_LOGFC1", "_TOPGO_enrichments_UNADJUSTED", cutoff, ".csv"))) {
			print("skip")
		}else{
            if(length(curr_sig$gene[curr_sig$log2FoldChange < -1]) < 5) {
                print("not doing")
                next()
            }
		test_GO_with_TOPGO(paste0("/thing/temp_DOWNREG_BAT_", prefix), bat_map, paste0("BAT_BP_ENRICHMENTS_SLIM_STRINGENT/", prefix, "_DOWNREG_LOGFC1"), 0.05)
		}
		print(paste0(prefix, "_DOWNREG"))
		
	}
	
}


test_GO_with_TOPGO <- function(targets, global_map, outfix, cutoff = 0.05) {
	
	library(topGO)
	geneID2GO <- readMappings(file = global_map)

	all_gene_background <- unique(names(geneID2GO))
	
	target_genes <- readLines(targets)
	
	geneList <- factor(as.integer(all_gene_background %in% target_genes))
	names(geneList) <- all_gene_background
		
	myGOdata <- new("topGOdata", description = 'thing', ontology = "BP", allGenes = geneList,
		annot = annFUN.gene2GO, gene2GO = geneID2GO)
	##optional 'nodeSize' to prune small GO terms?
	
	#resultFisher <- runTest(myGOData, algorithm = "classic", statistic = "fisher")
	
	##weighing hierarchy:
	#resultFisher_weight <- runTest(myGOData, algorithm = "weight01", statistic = "fisher")
	
	##....
	resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
mysummary <- summary(attributes(resultTopgo)$score <= cutoff)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001

allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)
write.csv(allRes, paste0(outfix, "_TOPGO_enrichments_UNADJUSTED", cutoff, ".csv"), row.names = F)

##padjustments though..
allGO = usedGO(object = myGOdata) 
BIG_allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = length(allGO))
BIG_allRes$PADJ <- p.adjust(BIG_allRes$topgoFisher, method = "BH")
write.csv(BIG_allRes, paste0(outfix, "_TOPGO_enrichments_FDR_APPLIED", ".csv"), row.names = F)
write.csv(BIG_allRes[BIG_allRes$PADJ < 0.05,], paste0(outfix, "_TOPGO_enrichments_FDR", cutoff, "_SLICED.csv"), row.names = F)

#output_file2 <- paste0(outfix, "_TOPGO_FDR", cutoff)
#printGraph(myGOdata, resultTopgo, firstSigNodes = length(which(BIG_allRes$PADJ < cutoff)), fn.prefix = output_file2, useInfo = "all", pdfSW = TRUE)

#output_file2 <- paste0(outfix, "_TOPGO_SLICED5")
#printGraph(myGOdata, resultTopgo, firstSigNodes = 10, fn.prefix = output_file2, useInfo = "all", pdfSW = TRUE)

##additional output?
term_matches <- list()

####

if(length(which(BIG_allRes$PADJ < cutoff)) == 0) {
	print("nothing doing")
	return(NA)
}else {
	print("grab enriches")

#myterms <- allRes$GO.ID
myterms <- BIG_allRes$GO.ID[BIG_allRes$PADJ < cutoff]
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms))
{
   myterm <- myterms[i]
   mygenesforterm <- mygenes[myterm][[1]]
   myfactor <- mygenesforterm %in% target_genes # find the genes that are in the list of genes of interest
   mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
   mygenesforterm2 <- paste(mygenesforterm2, collapse='@')
   #print(paste("Term",myterm,"genes:",mygenesforterm2))
	term_matches[[myterm]] <- mygenesforterm2
	print(myterm)
}
out_frame <- data.frame(term = names(term_matches), terms = unlist(term_matches), stringsAsFactors = F)
write.csv(out_frame, paste0(outfix, "_enriched_TOPGO_FDR", cutoff, ".csv"), row.names = F)


fancy_graph <- function() {
	
require(ggplot2)
library(scales)

ggdata <- BIG_allRes[BIG_allRes$PADJ < 0.05,]
ggdata <- ggdata[!duplicated(ggdata$Term),]
ggdata <- ggdata[1:min(c(dim(ggdata)[1], 20)),]
ggdata <- ggdata[!is.na(ggdata[,1]),]

ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
gg1 <- ggplot(ggdata,
  aes(x = Term, y = -log10(PADJ), size = Significant, fill = -log10(PADJ))) +

  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO Biological processes',
    subtitle = 'Significant terms by weighted-Fisher adjusted p-value.',
    caption = 'Cut-off lines drawn at equivalents of padj. p=0.05, p=0.01, p=0.001') +

  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    linetype = c("dotted", "longdash", "solid"),
    colour = c("black", "black", "black"),
    size = c(0.5, 1.5, 3)) +

  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),

    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),

    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +

  coord_flip()
  
	pdf(paste0(outfix, "_TOPGO_FDR", cutoff, ".pdf"), width = 16, height = 12)
	print(gg1)
	dev.off()
	
}

fancy_graph()

}

}

##HUMAN CONTROL
targets <- "HUMAN_CALU3_24H_INF_UPREG_genes.txt"
global_map <- "HUMAN_EGGNOG_TRANSCRIPT_ALIGNED_TOPGO_INPUT.output"
outfix <- "HUMAN_CALU3_24H_UPREG_EXTERNAL_VAL"
cutoff = 0.05

test_GO_with_TOPGO_ANNOTATION_DATABASE <- function(targets, global_map, outfix, cutoff = 0.05) {
	
	all_gene_background <- system(paste0("cut -f1 ", global_map), intern=T)
	
	library(topGO)
	library(org.Hs.eg.db)
	
	target_genes <- readLines(targets)
	
	geneList <- factor(as.integer(all_gene_background %in% target_genes))
	names(geneList) <- all_gene_background
	
	myGOdata <- new("topGOdata", description = 'thing', ontology = "BP", allGenes = geneList, annot = annFUN.org,
		mapping = "org.Hs.eg.db", ID = "ensembl")

	#resultFisher <- runTest(myGOData, algorithm = "classic", statistic = "fisher")
	
	##weighing hierarchy:
	#resultFisher_weight <- runTest(myGOData, algorithm = "weight01", statistic = "fisher")
	
	##....
	resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
mysummary <- summary(attributes(resultTopgo)$score <= cutoff)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001

allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = numsignif)
write.csv(allRes, paste0(outfix, "_TOPGO_enrichments_UNADJUSTED", cutoff, ".csv"), row.names = F)

##padjustments though..
allGO = usedGO(object = myGOdata) 
BIG_allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", topNodes = length(allGO))
BIG_allRes$PADJ <- p.adjust(BIG_allRes$topgoFisher, method = "BH")
write.csv(BIG_allRes, paste0(outfix, "_TOPGO_enrichments_FDR_APPLIED", ".csv"), row.names = F)
write.csv(BIG_allRes[BIG_allRes$PADJ < 0.05,], paste0(outfix, "_TOPGO_enrichments_FDR", cutoff, "_SLICED.csv"), row.names = F)

#output_file2 <- paste0(outfix, "_TOPGO_FDR", cutoff)
#printGraph(myGOdata, resultTopgo, firstSigNodes = length(which(BIG_allRes$PADJ < cutoff)), fn.prefix = output_file2, useInfo = "all", pdfSW = TRUE)

#output_file2 <- paste0(outfix, "_TOPGO_SLICED5")
#printGraph(myGOdata, resultTopgo, firstSigNodes = 10, fn.prefix = output_file2, useInfo = "all", pdfSW = TRUE)

##additional output?
term_matches <- list()

if(length(which(BIG_allRes$PADJ < cutoff)) == 0) {
	print("nothing doing")
	return(NA)
}else {
	
#myterms <- allRes$GO.ID
myterms <- BIG_allRes$GO.ID[BIG_allRes$PADJ < cutoff]
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms))
{
   myterm <- myterms[i]
   mygenesforterm <- mygenes[myterm][[1]]
   myfactor <- mygenesforterm %in% target_genes # find the genes that are in the list of genes of interest
   mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
   mygenesforterm2 <- paste(mygenesforterm2, collapse='@')
   #print(paste("Term",myterm,"genes:",mygenesforterm2))
	term_matches[[myterm]] <- mygenesforterm2
	print(myterm)
}
out_frame <- data.frame(term = names(term_matches), terms = unlist(term_matches))
write.csv(out_frame, paste0(outfix, "_enriched_TOPGO_FDR", cutoff, ".csv"), row.names = F)

fancy_graph <- function() {
		
	require(ggplot2)
library(scales)

ggdata <- BIG_allRes[BIG_allRes$PADJ < 0.05,]

ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
gg1 <- ggplot(ggdata,
  aes(x = Term, y = -log10(PADJ), size = Significant, fill = -log10(PADJ))) +

  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +

  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'GO Biological processes',
    subtitle = 'Significant terms by weighted-Fisher adjusted p-value.',
    caption = 'Cut-off lines drawn at equivalents of padj. p=0.05, p=0.01, p=0.001') +

  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    linetype = c("dotted", "longdash", "solid"),
    colour = c("black", "black", "black"),
    size = c(0.5, 1.5, 3)) +

  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),

    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),

    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +

  coord_flip()
  
	pdf(paste0(outfix, "_TOPGO_FDR", cutoff, ".pdf"), width = 16, height = 12)
	print(gg1)
	dev.off()
	
}

fancy_graph()

}

}
