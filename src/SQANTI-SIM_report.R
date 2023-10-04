#######################################
#                                     #
#    SQANTI-SIM report generation     #
#                                     #
#######################################

#Author: Jorge Mestre Tomas (jorge.mestre.tomas@csic.es)

#######################################
#                                     #
#      PACKAGES AND LIBRARIES         #
#                                     #
#######################################

suppressMessages(library(dplyr))
suppressMessages(library(DT))
suppressMessages(library(fmsb))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(gridExtra))
suppressMessages(library(knitr))
suppressMessages(library(rmarkdown))
suppressMessages(library(tidyr))


#######################################
#                                     #
#      FUNCTIONS  AND CLASSES         #
#                                     #
#######################################

gene_level_metrics <- function(data.query, data.index, MAX_TSS_TTS_DIFF){
  data.index <- data.index[data.index$sim_counts > 1, ]
  TP_genes <- unique(data.index$gene_id[data.index$pipeline_performance != "FN"])
  FN_genes <- unique(data.index$gene_id[!data.index$gene_id %in% TP_genes])
  FP_genes <- unique(data.query$associated_gene[!(data.query$associated_gene %in% TP_genes | data.query$associated_gene %in% FN_genes) & !(data.query$isoform %in% data.index$pipeline_performance)])
  
  TP <-  length(TP_genes)
  FP <- length(FP_genes)
  FN <- length(FN_genes)
  PTP <- 0
  sensitivity <- TP/(TP + FN)
  precision <- TP/(TP+FP)
  
  gene.metrics <- c("Total" = (TP + FN), "TP" = TP, "PTP" = PTP,
                    "FP" = FP, "FN" = FN, "Sensitivity" = sensitivity, "Precision" = precision, 
                    "F-score" = (2*(precision*sensitivity)/(precision+sensitivity)),
                    "False_Discovery_Rate" = (FP + PTP) / (FP + PTP +  TP),
                    "Positive_Detection_Rate" = (TP + PTP) / (TP + FN),
                    "False_Detection_Rate" = (FP) / (FP + PTP + TP))
  
  return(gene.metrics)
}

isoform_level_metrics <- function(data.query, data.index, MAX_TSS_TTS_DIFF, min_supp=0){
  # Simulated ref transcripts
  if (min_supp > 0){
    idx <- data.index[which(data.index$sim_counts >= min_supp), ]
  } else {
    idx <- data.index
  }
  data.novel <- idx[which(idx$sim_type == 'novel'),]
  data.known <- idx[which(idx$sim_type == 'known'),]
  sim.sc <- unique(data.novel$structural_category)

  # Matches between simulated and reconstructed transcripts:
  # First all splice-junctions must be identical
  # Second, difference between the annotated and reconstructed TSS and TTS must be smaller than MAX_TSS_TTS_DIFF
  matches <- inner_join(data.query, idx, by=c('junctions', 'chrom'), relationship = "many-to-many") %>%
    mutate(diffTSS = abs(TSS_genomic_coord.x - TSS_genomic_coord.y), diffTTS = abs(TTS_genomic_coord.x - TTS_genomic_coord.y), difftot = diffTSS+diffTTS) %>%
    arrange(difftot) %>%
    distinct(isoform, .keep_all = T)
  matches <- matches[which(matches$sim_type != "absent"),]
  perfect.matches <- matches[which(matches$diffTSS < MAX_TSS_TTS_DIFF & matches$diffTTS < MAX_TSS_TTS_DIFF),]
  cond <- (perfect.matches$exons > 1) | (perfect.matches$strand.x == '+' & perfect.matches$TSS_genomic_coord.x <= perfect.matches$TTS_genomic_coord.y & perfect.matches$TSS_genomic_coord.y <= perfect.matches$TTS_genomic_coord.x) | (perfect.matches$strand.x == '-' & perfect.matches$TTS_genomic_coord.x <= perfect.matches$TSS_genomic_coord.y & perfect.matches$TTS_genomic_coord.y <= perfect.matches$TSS_genomic_coord.x)
  perfect.matches <- perfect.matches[cond,]
  matches <- matches[which(matches$junctions != "" | matches$isoform %in% perfect.matches$isoform),]  # Delete PTP of mono-exons which always match because there is no splice junctions
  
  known.matches <- matches[which(matches$sim_type == "known"),]
  known.perfect.matches <- perfect.matches[which(perfect.matches$sim_type == "known"),]
  novel.matches <- matches[which(matches$sim_type == "novel"),]
  novel.perfect.matches <- perfect.matches[which(perfect.matches$sim_type == "novel"),]
  
  # Summary dataframe
  plot.cols <- c("transcript_id", "gene_id", "structural_category", "exons", "length", "sim_type", "sim_counts")
  if ('within_CAGE_peak' %in% colnames(data.index)){
    plot.cols <- c(plot.cols, "within_CAGE_peak", "dist_to_CAGE_peak")
  }
  if ('min_cov' %in% colnames(data.index)) {
    plot.cols <- c(plot.cols, "min_cov")
  }
  if ('ratio_TSS' %in% colnames(data.index)) {
    plot.cols <- c(plot.cols, "ratio_TSS")
  }
  
  data.summary <- rbind(
    data.novel[, plot.cols],
    data.known[, plot.cols]
  )
  data.summary$structural_category <- as.character(data.summary$structural_category)
  data.summary$structural_category[which(data.summary$transcript_id %in% data.known$transcript_id)] <- "Known"
  data.summary$structural_category <- factor(data.summary$structural_category, levels=c("Known", "FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron"))
  
  if (nrow(data.summary) > 0){
    data.summary$sim_match_type <- "FN"
    data.summary$sim_match_type[which(data.summary$transcript_id %in% data.known$transcript_id)] <-  "FN_known"
    data.summary$sim_match_type[which(data.summary$transcript_id %in% data.novel$transcript_id)] <-  "FN_novel"
    data.summary$sim_match_type[which(data.summary$transcript_id %in% known.perfect.matches$transcript_id)] <-  "TP_known"
    data.summary$sim_match_type[which(data.summary$transcript_id %in% novel.perfect.matches$transcript_id)] <-  "TP_novel"
    data.summary$sim_match_type <- factor(data.summary$sim_match_type, levels=c("TP_known", "FN_known", "TP_novel", "FN_novel"))
    
    data.summary$match_type <- "FN"
    
    PTP <- unique(matches[!(matches$isoform %in% perfect.matches$isoform), c("transcript_id","chrom", "junctions")])
    PTP$chrom_junc <- paste(PTP$chrom, PTP$junctions, sep = "_")
    PTP <- PTP[!(PTP$chrom_junc %in% paste(perfect.matches$chrom, perfect.matches$junctions, sep = "_")), ]
    
    data.summary$match_type[data.summary$transcript_id %in% PTP$transcript_id] <- "PTP"
    data.summary$match_type[which(data.summary$transcript_id %in% known.perfect.matches$transcript_id | data.summary$transcript_id %in% novel.perfect.matches$transcript_id)] <- "TP"
  }
  
  data.summary <- merge(data.summary, perfect.matches[,c("isoform", "transcript_id")], by="transcript_id", all.x=T)
  
  #  Compute metrics
  sqantisim.stats <- data.frame(init=c())
  for (sc in xaxislabelsF1){
    
    if (sc == 'FSM') {
      n.sim <- nrow(data.known)
      TP <- sum(data.summary$match_type == "TP" & data.summary$sim_type == "known")
      PTP <- sum(data.summary$match_type == "PTP" & data.summary$sim_type == "known")
      FN <- sum(data.summary$match_type == "FN" & data.summary$sim_type == "known") + PTP
      FP <- sum(!(data.query$isoform %in% data.summary$isoform) & data.query$structural_category == "FSM")
    } else {
      n.sim <- nrow(data.novel[which(data.novel$structural_category == sc),])
      TP <- sum(data.summary$match_type == "TP" & data.summary$sim_type == "novel" & data.summary$structural_category == sc)
      PTP <- sum(data.summary$match_type == "PTP" & data.summary$sim_type == "novel" & data.summary$structural_category == sc)
      FN <- sum(data.summary$match_type == "FN" & data.summary$sim_type == "novel" & data.summary$structural_category == sc) + PTP
      FP <- sum(!(data.query$isoform %in% data.summary$isoform) & data.query$structural_category == sc)
    }
    #FN <- n.sim - TP
    
    if (sum(n.sim, TP, PTP, FP, FN) > 0){
      sqantisim.stats['Total', sc] <- n.sim
      sqantisim.stats['TP', sc] <- TP
      sqantisim.stats['PTP', sc] <- PTP
      sqantisim.stats['FP', sc] <- FP
      sqantisim.stats['FN', sc] <- FN
      sqantisim.stats['Sensitivity', sc] <- TP/ (TP + FN)
      sqantisim.stats['Precision', sc] <- TP/ (TP + FP)
      sqantisim.stats['F-score', sc] <- 2*((sqantisim.stats['Sensitivity', sc]*sqantisim.stats['Precision', sc])/(sqantisim.stats['Sensitivity', sc]+sqantisim.stats['Precision', sc]))
      sqantisim.stats['False_Discovery_Rate', sc] <- (FP) / (FP +  TP)
      sqantisim.stats['Positive_Detection_Rate', sc] <- (TP + PTP) / (TP + FN)
      sqantisim.stats['False_Detection_Rate', sc] <- (FP) / (FP + PTP + TP)
    } else {
      sqantisim.stats['Total', sc] <- 0
      sqantisim.stats['TP', sc] <- 0
      sqantisim.stats['PTP', sc] <- 0
      sqantisim.stats['FP', sc] <- 0
      sqantisim.stats['FN', sc] <- 0
      sqantisim.stats['Sensitivity', sc] <- 0
      sqantisim.stats['Precision', sc] <- 0
      sqantisim.stats['F-score', sc] <- 0
      sqantisim.stats['False_Discovery_Rate', sc] <- 0
      sqantisim.stats['Positive_Detection_Rate', sc] <- 0
      sqantisim.stats['False_Detection_Rate', sc] <- 0 
    }
  }
  
  known.TP <- 0
  known.PTP <- 0
  known.FN <- 0
  known.FP <- 0
  novel.TP <- 0
  novel.PTP <- 0
  novel.FN <- 0 
  novel.FP <- 0
  total.FP <- 0
  for (sc in colnames(sqantisim.stats)){
    if (sc == "FSM"){
      known.TP <- known.TP + sqantisim.stats['TP', sc]
      known.PTP <- known.PTP + sqantisim.stats['PTP', sc]
      known.FN <- known.FN + sqantisim.stats['FN', sc]
      known.FP <- known.FP + sqantisim.stats['FP', sc]
    }else{
      novel.TP <- novel.TP + sqantisim.stats['TP', sc]
      novel.PTP <- novel.PTP + sqantisim.stats['PTP', sc]
      novel.FN <- novel.FN + sqantisim.stats['FN', sc]
      novel.FP <- novel.FP + sqantisim.stats['FP', sc]
    }
  }
  
  sqantisim.stats['Total', 'Known'] <-  nrow(data.known)
  sqantisim.stats['TP', 'Known'] <- known.TP
  sqantisim.stats['PTP', 'Known'] <- known.PTP
  sqantisim.stats['FP', 'Known'] <- known.FP
  sqantisim.stats['FN', 'Known'] <- known.FN
  sqantisim.stats['Precision', 'Known'] <- known.TP / (known.TP + known.FP)
  sqantisim.stats['Sensitivity', 'Known'] <- known.TP / (known.TP + known.FN)
  sqantisim.stats['F-score', 'Known'] <- 2*((sqantisim.stats['Sensitivity', 'Known']*sqantisim.stats['Precision', 'Known'])/(sqantisim.stats['Sensitivity', 'Known']+sqantisim.stats['Precision', 'Known']))
  sqantisim.stats['False_Discovery_Rate', 'Known'] <- (known.FP) / (known.FP +  known.TP)
  sqantisim.stats['Positive_Detection_Rate', 'Known'] <- (known.TP + known.PTP) / (known.TP + known.FN)
  sqantisim.stats['False_Detection_Rate', 'Known'] <- (known.FP) / (known.FP + known.PTP +  known.TP)
  
  sqantisim.stats['Total', 'Novel'] <-  nrow(data.novel)
  sqantisim.stats['TP', 'Novel'] <- novel.TP
  sqantisim.stats['PTP', 'Novel'] <- novel.PTP
  sqantisim.stats['FP', 'Novel'] <- novel.FP
  sqantisim.stats['FN', 'Novel'] <- novel.FN
  sqantisim.stats['Precision', 'Novel'] <- novel.TP / (novel.TP + novel.FP)
  sqantisim.stats['Sensitivity', 'Novel'] <- novel.TP / (novel.TP + novel.FN)
  sqantisim.stats['F-score', 'Novel'] <- 2*((sqantisim.stats['Sensitivity', 'Novel']*sqantisim.stats['Precision', 'Novel'])/(sqantisim.stats['Sensitivity', 'Novel']+sqantisim.stats['Precision', 'Novel']))
  sqantisim.stats['False_Discovery_Rate', 'Novel'] <- (novel.FP) / (novel.FP +  novel.TP)
  sqantisim.stats['Positive_Detection_Rate', 'Novel'] <- (novel.TP + novel.PTP) / (novel.TP + novel.FN)
  sqantisim.stats['False_Detection_Rate', 'Novel'] <- (novel.FP) / (novel.FP + novel.PTP +  novel.TP)
  
  col.order <- c("Known", "Novel", "FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
  row.order <- c('Total', 'TP', 'PTP', 'FP', 'FN', 'Sensitivity', 'Precision', 'F-score', 'False_Discovery_Rate', 'Positive_Detection_Rate', 'False_Detection_Rate')
  sqantisim.stats <- sqantisim.stats[intersect(row.order, rownames(sqantisim.stats)), intersect(col.order, colnames(sqantisim.stats))]
  
  res <- list(data.summary, known.perfect.matches, novel.perfect.matches, sqantisim.stats)
  names(res) <- c("data.summary", "known.perfect.matches", "novel.perfect.matches", "sqantisim.stats")
  return(res)
}

modify_index_file <- function(index.file, res.full, output_directory){
  modif.index <- read.table(index.file, header=T, as.is=T, sep="\t")
  n <- colnames(modif.index)
  modif.index <- merge(x = modif.index, y = res.full$data.summary[,c("transcript_id", "isoform")], by = "transcript_id", all.x = TRUE)
  colnames(modif.index) <- c(n, "pipeline_performance")
  modif.index$pipeline_performance[is.na(modif.index$pipeline_performance)] <- "FN"
  modif.index$pipeline_performance[which(modif.index$sim_counts <= 0)] <-  "absent"
  #modif.index$pipeline_performance[which(modif.index$transcript_id %in% res.full$data.summary$transcript_id[which(res.full$data.summary$match_type == "PTP")])] <- "PTP"
  #modif.index$pipeline_performance[which(modif.index$transcript_id %in% res.full$data.summary$transcript_id[which(res.full$data.summary$match_type == "TP")])] <- "TP"
  return(modif.index)
}

#######################################
#                                     #
#                MAIN                 #
#                                     #
#######################################

# -------------------- Argument parser
args <- commandArgs(trailingOnly = TRUE)
class.file <- args[1] # classification file SQANTI3
junc.file <- args[2] # junctions file SQANTI3
index.file <- args[3] # index file
min.supp <- as.numeric(args[4]) # min support reads
src.path <- args[5] # path to src utilities
quant.file <- args[6] # quantification file

output_directory <- dirname(dirname(class.file))
output_name <- basename(strsplit(class.file, "_classification.txt")[[1]][1])

# -------------------- Read files
# Read classification file
data.class <- read.table(class.file, header=T, as.is=T, sep="\t")
rownames(data.class) <- data.class$isoform
xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
data.class$structural_category = factor(data.class$structural_category,
                                        labels = xaxislabelsF1,
                                        levels = xaxislevelsF1,
                                        ordered=TRUE)

# Read juncs file
data.junction <- read.table(junc.file, header=T, as.is=T, sep="\t")
data.junction <- data.junction %>%
  group_by(isoform) %>%
  summarise(Donors=list(as.numeric(genomic_start_coord)), Acceptors=list(as.numeric(genomic_end_coord)))
data.junction$Donors <- lapply(data.junction$Donors, function(x){
  paste(sort(unlist(x)), collapse=',')
})
data.junction$Acceptors <- lapply(data.junction$Acceptors, function(x){
  paste(sort(unlist(x)), collapse=',')
})
data.junction$junctions <- paste(data.junction$Donors, data.junction$Acceptors, sep=',')

# Combine class and junc
data.query <- full_join(data.class, data.junction, by='isoform')
data.query$junctions[which(is.na(data.query$junctions))] <- ''
data.query <- data.query[,c('isoform', 'chrom', 'strand', 'structural_category', 'associated_gene', 'junctions', 'TSS_genomic_coord', 'TTS_genomic_coord', 'all_canonical', 'dist_to_CAGE_peak', 'within_CAGE_peak', 'min_cov', 'ratio_TSS')]

# Read index file
data.index <- read.table(index.file, header=T, as.is=T, sep="\t")
data.index$donors <- lapply(data.index$donors, function(x){
  paste(sort(unlist(as.numeric(unlist(strsplit(as.character(x), ",")))+1)), collapse=',')
})
data.index$acceptors <- lapply(data.index$acceptors, function(x){
  paste(sort(unlist(as.numeric(unlist(strsplit(as.character(x), ",")))-1)), collapse=',')
})
data.index$junctions <- paste(data.index$donors, data.index$acceptors, sep=',')
data.index$junctions[which(data.index$junctions == ',')] <- ''
data.index[which(data.index$strand == "-"),"TSS_genomic_coord"]  <- data.index[which(data.index$strand == "-"),"TSS_genomic_coord"] - 1
data.index[which(data.index$strand == "-"),"TTS_genomic_coord"]  <- data.index[which(data.index$strand == "-"),"TTS_genomic_coord"] + 1
data.index$structural_category = factor(data.index$structural_category,
                                      labels = xaxislabelsF1,
                                      levels = xaxislevelsF1,
                                      ordered=TRUE)
data.index$donors <- NULL
data.index$acceptors <- NULL
data.index$sim_type[which(data.index$sim_counts == 0)] <- 'absent' # Ignore not simulated

# -------------------- Performance metrics
# Matched for novel and known
MAX_TSS_TTS_DIFF = 50
res.full <- isoform_level_metrics(data.query, data.index, MAX_TSS_TTS_DIFF)
res.full[is.na(res.full)] <-  0
res.min <- isoform_level_metrics(data.query, data.index, MAX_TSS_TTS_DIFF, min.supp)
res.min$sqantisim.stats <- res.min$sqantisim.stats[c("Total", "TP", "FN", "Sensitivity"),]
res.min[is.na(res.min)] <-  0

modif.index <- modify_index_file(index.file, res.full, output_directory)
index_file_name <- basename(index.file)
modif_index_path <- file.path(output_directory, paste0(substr(index_file_name, 1, nchar(index_file_name) - 4), ".eval.tsv"))
print(paste("Modified index file written to:", modif_index_path))

write.table(modif.index, file = modif_index_path, quote = F, sep = "\t", na = "NA", row.names = F)

res.full$data.summary$match_type[which(res.full$data.summary$match_type == "PTP")] <- "FN"
res.full$data.summary$match_type <- factor(res.full$data.summary$match_type, levels = c("TP", "FN"))
res.min$data.summary$match_type[which(res.min$data.summary$match_type == "PTP")] <- "FN"
res.min$data.summary$match_type <- factor(res.min$data.summary$match_type, levels = c("TP", "FN"))

data.index <- merge(x = data.index, y = modif.index[,c("transcript_id", "pipeline_performance")], by = "transcript_id", all.x = TRUE)

res.gene <- gene_level_metrics(data.query, data.index, MAX_TSS_TTS_DIFF)

sel <- c("transcript_id", "gene_id", "structural_category", "exons", "length", "sim_type", "sim_counts", "pipeline_performance")
if ("within_CAGE_peak" %in% colnames(data.index)) {
  sel <- c(sel, "within_CAGE_peak")
}
if ("min_cov" %in% colnames(data.index)) {
  sel <- c(sel, "min_cov")
} 

data.index <- data.index[data.index[,"sim_counts"] > 0, sel]
data.index$pipeline_performance <- as.character(data.index$pipeline_performance)

data.class$pipeline_performance <- ifelse(data.class$isoform %in% data.index$pipeline_performance, "TP", "FP")

data.index$pipeline_performance[which(data.index$pipeline_performance != "FN")] <- "TP"
data.index$structural_category[data.index$sim_type == "known"] <- "FSM"

if ('min_cov' %in% colnames(data.index) & "within_CAGE_peak" %in% colnames(data.index)) {
  dataTPFNFP <- rbind(data.index[, c("structural_category", "exons", "length", "min_cov", "within_CAGE_peak", "pipeline_performance")],
                      data.class[data.class$pipeline_performance == "FP", c("structural_category", "exons", "length",  "min_cov", "within_CAGE_peak", "pipeline_performance")])
} else if ('min_cov' %in% colnames(data.index)) {
  dataTPFNFP <- rbind(data.index[, c("structural_category", "exons", "length", "min_cov", "pipeline_performance")],
                      data.class[data.class$pipeline_performance == "FP", c("structural_category", "exons", "length",  "min_cov", "pipeline_performance")])
} else if ("within_CAGE_peak" %in% colnames(data.index)) {
  dataTPFNFP <- rbind(data.index[, c("structural_category", "exons", "length", "within_CAGE_peak", "pipeline_performance")],
                      data.class[data.class$pipeline_performance == "FP", c("structural_category", "exons", "length", "within_CAGE_peak", "pipeline_performance")])
} else {
  dataTPFNFP <- rbind(data.index[, c("structural_category", "exons", "length", "pipeline_performance")],
                      data.class[data.class$pipeline_performance == "FP", c("structural_category", "exons", "length", "pipeline_performance")])
}


#######################################
#                                     #
#     TABLE AND PLOT GENERATION       #
#                                     #
#######################################

# -------------------- 
# -------------------- 
# TABLE INDEX
# t1: SQANTI-SIM metrics
# t2: SQANTI-SIM metrics above min_supp

# -------------------- 
# -------------------- 
# PLOT INDEX
# p1: Simulated and reconstructed isoform distribution
# p2: structural classification
# p3: TP vs FN - mono/multi-exon
# p4: canonical juncs
# p5: Number of exons
# p6: Transcript length
# p7: Expression levels
# p8: Transcripts per gene
# p9: within cage peak
# p10: distance to cage peak
# p11: min SJ cov by short reads
# p12: ratio TSS
# px: radar chart from perfomance metrics

print("***Generating plots for the report")

# -------------------- Global plot parameters
# SAME FORMAT AS SQANTI3 REPORT
#myPalette = c("#6BAED6","#FC8D59","#78C679","#EE6A50","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")
myPalette = c("#3A5A81", "#D31336", "#252131", "#6BAED6","#FC8D59","#78C679","#EE6A50")

cat.palette = c("Known"="#6BAED6", "FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                "NNC"="#EE6A50", "Genic\nGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                "Intergenic" = "darksalmon", "Genic\nIntron"="#41B6C4")

mytheme <- theme_classic(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", linewidth = 0.4),
        axis.line.y = element_line(color="black", linewidth = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=12) ) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size=12), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15, hjust = 0.5))

# -------------------- 
# TABLE 1: SQANTI-SIM metrics
sketch1 <- htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th(rowspan = 2, ""),
      th(colspan = 1, "Gene\nlevel"),
      th(colspan = 2, "Transcript level"),
      th(colspan = ncol(res.full$sqantisim.stats) - 2, "Structural categories")
    ),
    tr(
      lapply(c("Genes", colnames(res.full$sqantisim.stats)), th)
    )
  )
))

sketch2 <- htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th(rowspan = 2, ""),
      th(colspan = 2, "Transcript level"),
      th(colspan = sum(res.full$sqantisim.stats["Total",] > 0) - 2, "Structural categories")
    ),
    tr(
      lapply(colnames(res.full$sqantisim.stats[,res.full$sqantisim.stats["Total",] > 0]), th)
    )
  )
))


t1 <- DT::datatable(cbind("Genes"=res.gene, res.full$sqantisim.stats), 
                    container=sketch1, class = 'compact', extensions = "Buttons", 
                    options = list(pageLength = 15, dom = 't', buttons = c("copy", "csv", "pdf")),
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: center;','Table 3: ',
                      htmltools::em('SQANTI-SIM evaluation of transcriptome reconstruction.'))) %>%
  formatRound(c("Genes", colnames(res.full$sqantisim.stats)), digits = 3, rows=c(6:11), zero.print = 0)
write.table(cbind("Genes"=res.gene, res.full$sqantisim.stats), file = paste(output_directory, 'SQANTI-SIM_metrics.tsv', sep = "/"), quote = F, sep = "\t", na = "NA",row.names = F, col.names = c("Genes", gsub("\n", "_", colnames(res.full$sqantisim.stats))))

# TABLE 2: SQANTI-SIM metrics above min_supp
min.sqantisim.stats <- res.min$sqantisim.stats[, colnames(res.full$sqantisim.stats[,res.full$sqantisim.stats["Total",] > 0])]
t2 <- DT::datatable(min.sqantisim.stats, container = sketch2, class = 'compact', extensions = "Buttons", 
                    options = list(pageLength = 15, dom = 't', buttons = c("copy", "csv", "pdf")),
                    caption = htmltools::tags$caption(
                      style = 'caption-side: bottom; text-align: center;','Table 4: ',
                      htmltools::em(paste0('SQANTI-SIM sensititvity evaluation for transcripts with at least ', min.supp, " simulated reads.")))) %>%
  formatRound(colnames(min.sqantisim.stats), digits = 3, rows=c(4), zero.print = 0)
write.table(min.sqantisim.stats, file = paste(output_directory, 'SQANTI-SIM_metrics_min_supp.tsv', sep = "/"), quote = F, sep = "\t", na = "NA",row.names = F, col.names = gsub("\n", "_", colnames(min.sqantisim.stats)))

# -------------------- PLOT FULL
# PLOT 1: Simulated and reconstructed isoform distribution
p1A <- data.index[data.index$sim_counts >= 1, ] %>%
  mutate(structural_category = ifelse(sim_type == "known", "FSM", as.character(structural_category))) %>%
  mutate(structural_category = factor(structural_category, levels = xaxislabelsF1, ordered=TRUE)) %>% 
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100, fill=structural_category), color="black", linewidth=0.3, width=0.7) +
  scale_x_discrete(drop=FALSE) + 
  mytheme +
  geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size=12)) +
  scale_fill_manual(values = cat.palette, guide='none') +
  ggtitle("Simulated isoforms" ) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_blank()) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))

p1B <- data.class %>%
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100, fill=structural_category, alpha = pipeline_performance), color="black", linewidth=0.3, width=0.7) +
  scale_x_discrete(drop=FALSE) + 
  scale_alpha_manual(values = c(0.2, 1), name = "Trancript model calls") +
  mytheme +
  geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size=12)) +
  scale_fill_manual(values = cat.palette, guide='none') +
  ggtitle("Reconstructed transcript models" ) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_blank())

p1 <- annotate_figure(
  ggarrange(p1A, p1B, ncol = 2, align = "h", common.legend = T, legend = "bottom"),
  left = "Transcripts %"
)

# PLOT 2: Sn - Pr - F1
p2 <- cbind(data.frame(Gene = c(res.gene["Sensitivity"], res.gene["Precision"], res.gene["F-score"])),
      res.full$sqantisim.stats[c("Sensitivity", "Precision", "F-score"), colnames(res.full$sqantisim.stats[, res.full$sqantisim.stats["Total",] > 0])]) %>% 
  mutate(stat = rownames(.)) %>% 
  select(!FSM) %>% 
  pivot_longer(!stat, names_to = "sim_level", values_to = "value") %>% 
  mutate(sim_level = factor(sim_level,
                            levels = rev(c("Gene", colnames(res.full$sqantisim.stats))))) %>% 
  ggplot(aes(y = value, color = stat, x = sim_level)) +
  geom_point(size=3.5) +
  theme_classic() +
  mytheme +
  scale_x_discrete(labels = c("Gene" = "Gene level", "Known" = "Kwnown\ntranscript level", "Novel" = "Novel\ntranscript level")) +
  scale_color_manual(values = c("#15918A", "#F58A53", "#FDC659"), name="", labels = c("F1-score", "Precision", "Sensitivity")) +
  ylim(c(0,1)) +
  coord_flip() +
  theme(panel.grid.major.y = element_line(color = "grey", linewidth = 0.5, linetype = 2)) +
  theme(legend.position="top") + 
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())

# PLOT 3: Sim reads
data.index$TPM <- (data.index$sim_counts/sum(data.index$sim_counts))*1000000
data.index$pipeline_performance <- factor(data.index$pipeline_performance, levels = c("TP", "FN", "FP"))
p3 <- ggplot(data.index[data.index$structural_category %in% colnames(res.full$sqantisim.stats[,res.full$sqantisim.stats["Total",] > 0]), ], aes(x=structural_category, y=log(TPM+1), fill=structural_category, alpha=pipeline_performance)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), color = "#656565", width = 0.5) +
  scale_fill_manual(values = cat.palette, name = "Structural category", guide = "none") +
  scale_alpha_manual(values=c("TP" = 1, "FN" = 0.5), name="Pipeline performance", guide = "none") +
  theme_classic() +
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Log simulated reads") +
  coord_cartesian(ylim = c(quantile(log(dataTPFNFP$TPM + 1), probs = c(0.05, 0.95)))) + 
  theme(axis.title.x = element_blank(),
        axis.title.y  = element_text(size=16)) +
  theme(plot.subtitle = element_text(hjust = 0.5),
        plot.title.position = "plot") +
  theme(legend.position = "bottom") +
  guides(alpha = guide_legend(override.aes = list(linetype = 0, fill = "darkgrey"))) +
  ggtitle("Simulated reads by structural category")

# PLOT 4: Number of exons
dataTPFNFP$pipeline_performance <- factor(dataTPFNFP$pipeline_performance, levels = c("TP", "FN", "FP"))
p4 <- ggplot(dataTPFNFP[dataTPFNFP$structural_category %in% colnames(res.full$sqantisim.stats[,res.full$sqantisim.stats["Total",] > 0]), ], aes(x=structural_category, y=exons, fill=structural_category, alpha=pipeline_performance)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), color = "#656565", width = 0.5) +
  scale_fill_manual(values = cat.palette, name = "Structural category", guide = "none") +
  scale_alpha_manual(values=c("TP" = 1, "FN" = 0.5, "FP" = 0), name="Pipeline performance") +
  scale_y_continuous(limits = quantile(dataTPFNFP$exons, c(0.1, 0.9))) +
  theme_classic() +
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("") +
  ylab("Number of exons") +
  coord_cartesian(ylim = c(quantile(dataTPFNFP$exons, probs = c(0.05, 0.95)))) + 
  theme(axis.title.x = element_blank(),
        axis.title.y  = element_text(size=16),
        axis.text.x = element_blank()) +
  theme(plot.subtitle = element_text(hjust = 0.5),
        plot.title.position = "plot") +
  theme(legend.position = "bottom") +
  guides(alpha = guide_legend(override.aes = list(linetype = 0, fill = "darkgrey"))) +
  ggtitle("Number of exons by structural category")

# PLOT 5: Transcript length
p5 <- ggplot(dataTPFNFP[dataTPFNFP$structural_category %in% colnames(res.full$sqantisim.stats[,res.full$sqantisim.stats["Total",] > 0]), ], aes(x=structural_category, y=log(length), fill=structural_category, alpha=pipeline_performance)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), color = "#656565", width = 0.5) +
  scale_fill_manual(values = cat.palette, name = "Structural category", guide = "none") +
  scale_alpha_manual(values=c("TP" = 1, "FN" = 0.5, "FP" = 0), name="Pipeline performance") +
  theme_classic() +
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("") +
  ylab("Log transcript length")+
  coord_cartesian(ylim = c(quantile(log(dataTPFNFP$length), probs = c(0.00005, 0.99995)))) + 
  theme(axis.title.x = element_blank(),
        axis.title.y  = element_text(size=16),
        axis.text.x = element_blank()) +
  theme(plot.subtitle = element_text(hjust = 0.5),
        plot.title.position = "plot") +
  theme(legend.position = "bottom") +
  guides(alpha = guide_legend(override.aes = list(linetype = 0, fill = "darkgrey"))) +
  ggtitle("Transcript length by structural category")

# PLOT 6: Isoform complexity (transcripts per gene)
trans.per.gene <- res.full$data.summary %>%
  group_by(gene_id) %>%
  summarise(trans_per_gene=n())
res.full$data.summary <- merge(res.full$data.summary, trans.per.gene, by = "gene_id")

p6 <- ggplot(res.full$data.summary[res.full$data.summary$structural_category %in% colnames(res.full$sqantisim.stats[,res.full$sqantisim.stats["Total",] > 0]) ,], aes(x=structural_category, y=trans_per_gene, fill=structural_category, alpha=match_type)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), color = "#656565", width = 0.5) +
  scale_fill_manual(values = cat.palette, name = "Structural category", guide = "none") +
  scale_alpha_manual(values=c("TP" = 1, "FN" = 0.5), name="Pipeline performance") +
  theme_classic() +
  mytheme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylab("Transcriptss per gene") +
  coord_cartesian(ylim = c(quantile(res.full$data.summary$trans_per_gene, probs = c(0.01, 0.99)))) + 
  theme(axis.title.x = element_blank(),
        axis.title.y  = element_text(size=16),
        axis.text.x = element_blank()) +
  theme(plot.subtitle = element_text(hjust = 0.5),
        plot.title.position = "plot") +
  theme(legend.position = "bottom") +
  guides(alpha = guide_legend(override.aes = list(linetype = 0, fill = "darkgrey"))) +
  ggtitle("Isoform complexity by structural category")

# PLOT 7: Canonical Junctions
p7 <- data.class[data.class$pipeline_performance == "FP" & !is.na(data.class$all_canonical), ] %>% 
  ggplot(aes(1, alpha = all_canonical, fill = structural_category)) +
  geom_bar(aes(y = ..count..), position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = cat.palette, name = "Structural category") +
  scale_alpha_manual(values = c(0, 1), guide = "none") +
  coord_flip() +
  theme_classic() +
  mytheme +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom") +
  ggtitle("FP with non-canonical splice junctions")

# PLOT 8: Short-read coverage (min SJ cov)
if ('min_cov' %in% colnames(data.index)) {
  dataTPFNFP$min_cov[is.na(dataTPFNFP$min_cov)] <- 0
  dataTPFNFP <- dataTPFNFP %>%
    mutate(SJ_cov = ifelse(min_cov >= 1, "True", "False")) %>%
    mutate(structural_category2 = ifelse(structural_category %in% as.character(unique(data.index$structural_category)), as.character(structural_category), "Other\nSCs")) %>% 
    na.omit()
  p8A <- ggplot(dataTPFNFP[dataTPFNFP$structural_category2 == "FSM",], aes(x = pipeline_performance, y = (..count..), fill = structural_category, alpha=SJ_cov)) +
    geom_bar(position = "stack") +
    facet_grid(. ~ structural_category) +
    scale_alpha_manual(values = c(0.5, 1), name = "Pipeline performance", labels = c("Not supported", "SR support")) +
    scale_fill_manual(values = cat.palette, guide = "none") +
    theme_classic() +
    mytheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  p8B <- ggplot(dataTPFNFP[dataTPFNFP$structural_category2 != "FSM",], aes(x = pipeline_performance, y = (..count..), fill = structural_category, alpha=SJ_cov)) +
    geom_bar(position = "stack") +
    facet_grid(. ~ structural_category) +
    scale_alpha_manual(values = c(0.5, 1), name = "Pipeline performance", labels = c("Not supported", "SR support")) +
    scale_fill_manual(values = cat.palette, guide = "none") +
    theme_classic() +
    mytheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  p8 <- annotate_figure(
    ggarrange(p8A, p8B, ncol = 2, align = "h", widths = c(1.5, length(unique(dataTPFNFP$structural_category)) - 1), common.legend = T, legend = "bottom"),
    top = "Transcripts supported by short reads"
  )
}

# PLOT 9: CAGE support
if ('within_CAGE_peak' %in% colnames(data.index)){
  dataTPFNFP$within_CAGE_peak <- ifelse(dataTPFNFP$within_CAGE_peak == T | dataTPFNFP$within_CAGE_peak == "True", "True", "False")
  p9A <- ggplot(dataTPFNFP[dataTPFNFP$structural_category2 == "FSM",], aes(x = pipeline_performance, y = (..count..), fill = structural_category, alpha=within_CAGE_peak)) +
    geom_bar(position = "stack") +
    facet_grid(. ~ structural_category) +
    scale_alpha_manual(values = c(0.5, 1), name = "Pipeline performance", labels = c("Not supported", "SR support")) +
    scale_fill_manual(values = cat.palette, guide = "none") +
    theme_classic() +
    mytheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  p9B <- ggplot(dataTPFNFP[dataTPFNFP$structural_category2 != "FSM",], aes(x = pipeline_performance, y = (..count..), fill = structural_category, alpha=within_CAGE_peak)) +
    geom_bar(position = "stack") +
    facet_grid(. ~ structural_category) +
    scale_alpha_manual(values = c(0.5, 1), name = "Pipeline performance", labels = c("Not supported", "CAGE support")) +
    scale_fill_manual(values = cat.palette, guide = "none") +
    theme_classic() +
    mytheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  p9 <- annotate_figure(
    ggarrange(p9A, p9B, ncol = 2, align = "h", widths = c(1.5, length(unique(dataTPFNFP$structural_category)) - 1), common.legend = T,  legend = "bottom"),
    top = "Transcripts supported by CAGE peaks"
  )
}

# PLOT 10: SR and CAGE
if ('min_cov' %in% colnames(data.index) & 'within_CAGE_peak' %in% colnames(data.index)) {
  dataTPFNFP$supp <- paste(dataTPFNFP$SJ_cov, dataTPFNFP$within_CAGE_peak, sep = "_")
  dataTPFNFP$supp <- ifelse(dataTPFNFP$supp == "True_True", "True", "False")
  dataTPFNFP$supp <- factor(dataTPFNFP$supp,
                           levels = c("False", "True"),
                           labels = c("Not fully\nsupported", "SR+CAGE\nsupport"))
  p10A <- ggplot(dataTPFNFP[dataTPFNFP$structural_category2 == "FSM",], aes(x = pipeline_performance, y = (..count..), fill = structural_category, alpha=supp)) +
    geom_bar(position = "stack") +
    facet_grid(. ~ structural_category) +
    scale_alpha_manual(values = c(0.5, 1), name = "Pipeline performance", labels = c("Not supported", "SR + CAGE support")) +
    scale_fill_manual(values = cat.palette, guide = "none") +
    theme_classic() +
    mytheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  p10B <- ggplot(dataTPFNFP[dataTPFNFP$structural_category2 != "FSM",], aes(x = pipeline_performance, y = (..count..), fill = structural_category, alpha=supp)) +
    geom_bar(position = "stack") +
    facet_grid(. ~ structural_category) +
    scale_alpha_manual(values = c(0.5, 1), name = "Pipeline performance", labels = c("Not supported", "SR support")) +
    scale_fill_manual(values = cat.palette, guide = "none") +
    theme_classic() +
    mytheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  p10 <- annotate_figure(
    ggarrange(p10A, p10B, ncol = 2, align = "h", widths = c(1.5, length(unique(dataTPFNFP$structural_category)) - 1), common.legend = T,  legend = "bottom"),
    top = "Transcripts supported by short reads and CAGE peaks"
  )
}

# PLOT X: Radar chart
# Generated in RMD file

# -------------------- Output report
rmarkdown::render(
  input = paste(src.path, 'SQANTI-SIM_report.Rmd', sep = "/"),
  output_dir = output_directory,
  output_file = paste0(output_name, "_SQANTI-SIM_report.html")
)

