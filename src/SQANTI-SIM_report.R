#######################################
#                                     #
#    SQANTI-SIM report generation     #
#                                     #
#######################################

# Author: Jorge Mestre Tomas (jormart2@alumni.uv.es)

#######################################
#                                     #
#      PACKAGES AND LIBRARIES         #
#                                     #
#######################################

suppressMessages(library(dplyr))
suppressMessages(library(DT))
suppressMessages(library(fmsb))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(knitr))
suppressMessages(library(rmarkdown))
suppressMessages(library(tidyr))


#######################################
#                                     #
#      FUNCTIONS  AND CLASSES         #
#                                     #
#######################################

gene_level_metrics <- function(data.query, data.index, MAX_TSS_TTS_DIFF, min_supp=0){
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
  matches <- inner_join(data.query, idx, by=c('junctions', 'chrom')) %>%
    mutate(diffTSS = abs(TSS_genomic_coord.x - TSS_genomic_coord.y), diffTTS = abs(TTS_genomic_coord.x - TTS_genomic_coord.y), difftot = diffTSS+diffTTS) %>%
    arrange(difftot) %>%
    distinct(isoform, .keep_all = T)
  #matches <- matches[which(matches$sim_type != "absent"),]
  perfect.matches <- matches[which(matches$diffTSS < MAX_TSS_TTS_DIFF & matches$diffTTS < MAX_TSS_TTS_DIFF),]
  cond <- (perfect.matches$exons > 1) | (perfect.matches$strand.x == '+' & perfect.matches$TSS_genomic_coord.x <= perfect.matches$TTS_genomic_coord.y & perfect.matches$TSS_genomic_coord.y <= perfect.matches$TTS_genomic_coord.x) | (perfect.matches$strand.x == '-' & perfect.matches$TTS_genomic_coord.x <= perfect.matches$TSS_genomic_coord.y & perfect.matches$TTS_genomic_coord.y <= perfect.matches$TSS_genomic_coord.x)
  perfect.matches <- perfect.matches[cond,]
  
  sim_genes <- idx$gene_id[idx$sim_counts > 0]
  sim_genes <- sim_genes[!duplicated(sim_genes)]
  matched_genes <- perfect.matches$gene_id[perfect.matches$gene_id %in% sim_genes]
  matched_genes <- matched_genes[!duplicated(matched_genes)]
  
  fp_genes <- data.query[!(data.query$isoform %in% perfect.matches$isoform),]
  fp_genes <- data.query[!(data.query$associated_gene %in% matched_genes),]
  fp_genes <- fp_genes$associated_gene[!duplicated(fp_genes$associated_gene)]
  
  TP <-  length(matched_genes)
  FP <- length(fp_genes)
  FN <- length(sim_genes)-length(matched_genes)
  PTP <- 0
  sensitivity <- TP/length(sim_genes)
  precision <- TP/(TP+FP)
  
  gene.metrics <- c("Total" = length(sim_genes), "TP" = TP, "PTP" = PTP,
                    "FP" = FP, "FN" = FN, "Sensitivity" = sensitivity, "Precision" = precision, 
                    "F-score" = (2*(precision*sensitivity)/(precision+sensitivity)),
                    "False_Discovery_Rate" = (FP + PTP) / (FP + PTP +  TP),
                    "Positive_Detection_Rate" = (TP + PTP) / length(sim_genes),
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
  matches <- inner_join(data.query, idx, by=c('junctions', 'chrom')) %>%
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
  
  #  Compute metrics
  sqantisim.stats <- data.frame(init=c())
  for (sc in xaxislabelsF1){
    
    if (sc == 'FSM') {
      n.sim <- nrow(data.known)
      TP <- nrow(known.perfect.matches[which(known.perfect.matches$structural_category.x == sc),])
      PTP <- nrow(known.matches[which(known.matches$structural_category.x == sc),]) - TP
      FP <- nrow(data.query[which(data.query$structural_category == sc),]) - TP - PTP
    } else {
      n.sim <- nrow(data.novel[which(data.novel$structural_category == sc),])
      TP <- nrow(novel.perfect.matches[which(novel.perfect.matches$structural_category.y == sc),])
      PTP <- nrow(novel.matches[which(novel.matches$structural_category.y == sc),]) - TP
      FP <- nrow(data.query[which(data.query$structural_category == sc),]) - TP - (nrow(novel.matches[which(novel.matches$structural_category.x == sc),]) - TP)
    }
    FN <- n.sim - TP

    if (sum(n.sim, TP, PTP, FP, FN) > 0){
      sqantisim.stats['Total', sc] <- n.sim
      sqantisim.stats['TP', sc] <- TP
      sqantisim.stats['PTP', sc] <- PTP
      sqantisim.stats['FP', sc] <- FP
      sqantisim.stats['FN', sc] <- FN
      sqantisim.stats['Sensitivity', sc] <- TP/ (TP + FN)
      sqantisim.stats['Precision', sc] <- TP/ (TP + PTP + FP)
      sqantisim.stats['F-score', sc] <- 2*((sqantisim.stats['Sensitivity', sc]*sqantisim.stats['Precision', sc])/(sqantisim.stats['Sensitivity', sc]+sqantisim.stats['Precision', sc]))
      sqantisim.stats['False_Discovery_Rate', sc] <- (FP + PTP) / (FP + PTP +  TP)
      sqantisim.stats['Positive_Detection_Rate', sc] <- (TP + PTP) / n.sim
      sqantisim.stats['False_Detection_Rate', sc] <- (FP) / (FP + PTP + TP)
    }
  }
  
  known.TP <- 0
  known.PTP <- 0
  known.FN <- 0
  novel.TP <- 0
  novel.PTP <- 0
  novel.FN <- 0 
  total.FP <- 0
  for (sc in colnames(sqantisim.stats)){
    if (sc == "FSM"){
      known.TP <- known.TP + sqantisim.stats['TP', sc]
      known.PTP <- known.PTP + sqantisim.stats['PTP', sc]
      total.FP <- total.FP + sqantisim.stats['FP', sc]
    }else{
      novel.TP <- novel.TP + sqantisim.stats['TP', sc]
      novel.PTP <- novel.PTP + sqantisim.stats['PTP', sc]
      total.FP <- total.FP + sqantisim.stats['FP', sc]
    }
  }
  
  sqantisim.stats['Total', 'Known'] <-  nrow(data.known)
  sqantisim.stats['TP', 'Known'] <- known.TP
  sqantisim.stats['PTP', 'Known'] <- known.PTP
  sqantisim.stats['FP', 'Known'] <- total.FP
  sqantisim.stats['FN', 'Known'] <- nrow(data.known) - known.TP
  sqantisim.stats['Precision', 'Known'] <- known.TP / (known.TP + known.PTP + total.FP + novel.PTP)
  sqantisim.stats['Sensitivity', 'Known'] <- known.TP / (known.TP + sqantisim.stats['FN', 'Known'])
  sqantisim.stats['F-score', 'Known'] <- 2*((sqantisim.stats['Sensitivity', 'Known']*sqantisim.stats['Precision', 'Known'])/(sqantisim.stats['Sensitivity', 'Known']+sqantisim.stats['Precision', 'Known']))
  sqantisim.stats['False_Discovery_Rate', 'Known'] <- (total.FP + known.PTP) / (total.FP + known.PTP +  known.TP)
  sqantisim.stats['Positive_Detection_Rate', 'Known'] <- (known.TP + known.PTP) / nrow(data.known)
  sqantisim.stats['False_Detection_Rate', 'Known'] <- (total.FP) / (total.FP + known.PTP +  known.TP)
  
  sqantisim.stats['Total', 'Novel'] <-  nrow(data.novel)
  sqantisim.stats['TP', 'Novel'] <- novel.TP
  sqantisim.stats['PTP', 'Novel'] <- novel.PTP
  sqantisim.stats['FP', 'Novel'] <- total.FP
  sqantisim.stats['FN', 'Novel'] <- nrow(data.novel) - novel.TP
  sqantisim.stats['Precision', 'Novel'] <- novel.TP / (novel.TP + novel.PTP + total.FP + known.PTP)
  sqantisim.stats['Sensitivity', 'Novel'] <- novel.TP / (novel.TP + sqantisim.stats['FN', 'Novel'])
  sqantisim.stats['F-score', 'Novel'] <- 2*((sqantisim.stats['Sensitivity', 'Novel']*sqantisim.stats['Precision', 'Novel'])/(sqantisim.stats['Sensitivity', 'Novel']+sqantisim.stats['Precision', 'Novel']))
  sqantisim.stats['False_Discovery_Rate', 'Novel'] <- (total.FP + novel.PTP) / (total.FP + novel.PTP +  novel.TP)
  sqantisim.stats['Positive_Detection_Rate', 'Novel'] <- (novel.TP + novel.PTP) / nrow(data.novel)
  sqantisim.stats['False_Detection_Rate', 'Novel'] <- (total.FP) / (total.FP + novel.PTP +  novel.TP)
  
  col.order <- c("Known", "Novel", "FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
  row.order <- c('Total', 'TP', 'PTP', 'FP', 'FN', 'Sensitivity', 'Precision', 'F-score', 'False_Discovery_Rate', 'Positive_Detection_Rate', 'False_Detection_Rate')
  sqantisim.stats <- sqantisim.stats[intersect(row.order, rownames(sqantisim.stats)), intersect(col.order, colnames(sqantisim.stats))]
  
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
    data.summary$match_type[which(data.summary$transcript_id %in% known.matches$transcript_id | data.summary$transcript_id %in% novel.matches$transcript_id)] <- "PTP"
    data.summary$match_type[which(data.summary$transcript_id %in% known.perfect.matches$transcript_id | data.summary$transcript_id %in% novel.perfect.matches$transcript_id)] <- "TP"
  }
  
  data.summary <- merge(data.summary, perfect.matches[,c("isoform", "transcript_id")], by="transcript_id", all.x=T)
  
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
min.supp <- args[4] # min support reads
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
res.min <- isoform_level_metrics(data.query, data.index, MAX_TSS_TTS_DIFF, min.supp)
res.min$sqantisim.stats <- res.min$sqantisim.stats[c("Total", "TP", "FN", "Sensitivity"),]

res.gene <- gene_level_metrics(data.query, data.index, MAX_TSS_TTS_DIFF)


modif.index <- modify_index_file(index.file, res.full, output_directory)
write.table(modif.index, file = paste(output_directory, paste0(substr(index.file, 1,nchar(index.file)-4), ".eval.tsv"), sep = "/"), quote = F, sep = "\t", na = "NA",row.names = F)

res.full$data.summary$match_type[which(res.full$data.summary$match_type == "PTP")] <- "FN"
res.full$data.summary$match_type <- factor(res.full$data.summary$match_type, levels = c("TP", "FN"))
res.min$data.summary$match_type[which(res.min$data.summary$match_type == "PTP")] <- "FN"
res.min$data.summary$match_type <- factor(res.min$data.summary$match_type, levels = c("TP", "FN"))

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
# p1: simulated expression profile
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

cat.palette = c("FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                "NNC"="#EE6A50", "Genic\nGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                "Intergenic" = "darksalmon", "Genic\nIntron"="#41B6C4")

mytheme <- theme_classic(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=12) ) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size=12), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15, hjust = 0.5)) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm"))

# -------------------- 
# TABLE 1: SQANTI-SIM metrics
sketch = htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th(rowspan = 2, ""),
      th(rowspan = 2, "Genes"),
      th(colspan = ncol(res.full$sqantisim.stats), "Isoform level")
    ),
    tr(
      lapply(colnames(res.full$sqantisim.stats), th)
    )
  )
))

t1 <- DT::datatable(cbind("Genes"=res.gene, res.full$sqantisim.stats), container=sketch, class = 'compact', extensions = "Buttons", options = list(pageLength = 15, dom = 'Bfrtip', buttons = c("copy", "csv", "pdf"))) %>%
  formatRound(c("Genes", colnames(res.full$sqantisim.stats)), digits = 3, rows=c(6:11), zero.print = 0)
write.table(cbind("Genes"=res.gene, res.full$sqantisim.stats), file = paste(output_directory, 'SQANTI-SIM_metrics.tsv', sep = "/"), quote = F, sep = "\t", na = "NA",row.names = F, col.names = c("Genes", gsub("\n", "_", colnames(res.full$sqantisim.stats))))

# TABLE 2: SQANTI-SIM metrics above min_supp
t2 <- DT::datatable(res.min$sqantisim.stats, class = 'compact', extensions = "Buttons", options = list(pageLength = 15, dom = 'Bfrtip', buttons = c("copy", "csv", "pdf"))) %>%
  formatRound(colnames(res.min$sqantisim.stats), digits = 3, rows=c(4), zero.print = 0)
write.table(res.min$sqantisim.stats, file = paste(output_directory, 'SQANTI-SIM_metrics_min_supp.tsv', sep = "/"), quote = F, sep = "\t", na = "NA",row.names = F, col.names = gsub("\n", "_", colnames(res.min$sqantisim.stats)))

# -------------------- PLOT FULL
# PLOT 1: simulated expression profile
expr.dist <- data.index[which(data.index$sim_type %in% c('novel', 'known')), c('sim_type', 'sim_counts')]

p1 <- expr.dist %>%
  ggplot(aes(x=sim_counts, fill=sim_type)) +
  geom_histogram(aes(y=stat(count)/sum(count)), color="black", alpha=0.5, position = 'identity', bins = round(sqrt(nrow(expr.dist)))) +
  mytheme +
  scale_fill_manual(values = myPalette[1:2], name='Transcript type') +
  xlab('Number of simulated reads') + 
  ylab('Percentage of transcripts') +
  ggtitle('Simulated expression levels')

# PLOT 2: structural classification
p2 <- data.class %>%
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100, fill=structural_category), color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=FALSE) + 
  xlab('') + 
  ylab('Transcripts %') +
  mytheme +
  geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = cat.palette, guide='none') +
  ggtitle("Isoform Distribution Across Structural Categories\n\n" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))

# PLOT 3: TP vs FN - mono/multi-exon
p3 <- res.full$data.summary %>%
  mutate(exon_type=ifelse(exons > 1, 'multi-exon', 'mono-exon')) %>%
  group_by(structural_category, match_type, exon_type) %>%
  summarise(value=n()) %>%
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(fill=match_type, y=value, alpha=exon_type), position="fill", stat="identity") +
  scale_fill_manual(values=myPalette, name='Stats') +
  scale_alpha_manual(values=c(0.5,1), name='Exons') +
  mytheme +
  ylab('Percentage %') +
  xlab('')+
  ggtitle('Single- and Multi-exon identifications')

p3.min <- res.min$data.summary %>%
  mutate(exon_type=ifelse(exons > 1, 'multi-exon', 'mono-exon')) %>%
  group_by(structural_category, match_type, exon_type) %>%
  summarise(value=n()) %>%
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(fill=match_type, y=value, alpha=exon_type), position="fill", stat="identity") +
  scale_fill_manual(values=myPalette, name='Stats') +
  scale_alpha_manual(values=c(0.5,1), name='Exons') +
  mytheme +
  ylab('Percentage %') +
  xlab('')+
  ggtitle('Single- and Multi-exon identifications (min. supp.)')


# PLOT 4: canonical juncs
data.query$match_type <- 'FP'
data.query$match_type[which(data.query$isoform %in% res.full$novel.perfect.matches$isoform)] <- 'TP'
data.query$match_type[which(data.query$isoform %in% res.full$known.perfect.matches$isoform)] <- 'TP'
p4 <- data.query[which(!is.na(data.query$all_canonical)),] %>%
  group_by(structural_category, match_type, all_canonical) %>%
  summarise(value=n()) %>%
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(fill=match_type, y=value, alpha=all_canonical), position="fill", stat="identity") +
  scale_fill_manual(values=myPalette[1:2], name='Stats') +
  scale_alpha_manual(values=c(1, 0.5), name='Junctions') +
  mytheme +
  ylab('Percentage %') +
  xlab('') +
  ggtitle('Canonical or Non Canonical Junctions') +
  theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))

# PLOT 5: Number of exons
p5.1 <- ggplot(res.full$data.summary, aes(x=match_type, y=exons, fill=match_type)) +
  geom_boxplot(alpha=1, outlier.shape = NA) +
  mytheme +
  scale_fill_manual(values=myPalette, name='Stats') +
  scale_y_continuous(limits = quantile(res.full$data.summary$exons, c(0.1, 0.9))) +
  ylab("Number of exons") +
  xlab("") +
  ggtitle("Number of exons") +
  guides(fill="none")

p5.2 <- ggplot(res.full$data.summary, aes(x=sim_match_type, y=exons, fill=sim_match_type)) +
  geom_boxplot(alpha=1, outlier.shape = NA) +
  mytheme +
  scale_fill_manual(values=myPalette, name='Stats') +
  scale_y_continuous(limits = quantile(res.full$data.summary$exons, c(0.1, 0.9))) +
  ylab("Number of exons") +
  xlab("") +
  ggtitle("Number of exons") +
  guides(fill="none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

if (nrow(res.min$data.summary) > 0){
  p5.min.1 <- ggplot(res.min$data.summary, aes(x=match_type, y=exons, fill=match_type)) +
    geom_boxplot(alpha=1, outlier.shape = NA) +
    mytheme +
    scale_fill_manual(values=myPalette, name='Stats') +
    scale_y_continuous(limits = quantile(res.min$data.summary$exons, c(0.1, 0.9))) +
    ylab("Number of exons") +
    xlab("") +
    ggtitle("Number of exons") +
    guides(fill="none")
  
  p5.min.2 <- ggplot(res.min$data.summary, aes(x=sim_match_type, y=exons, fill=sim_match_type)) +
    geom_boxplot(alpha=1, outlier.shape = NA) +
    mytheme +
    scale_fill_manual(values=myPalette, name='Stats') +
    scale_y_continuous(limits = quantile(res.min$data.summary$exons, c(0.1, 0.9))) +
    ylab("Number of exons") +
    xlab("") +
    ggtitle("Number of exons") +
    guides(fill="none") +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
}

# PLOT 6: Length of transcripts
p6.1 <- ggplot(res.full$data.summary, aes(x=match_type, y=length, fill=match_type)) +
  geom_boxplot(alpha=1, outlier.shape = NA) +
  mytheme +
  scale_fill_manual(values=myPalette, name='Stats') +
  scale_y_continuous(limits = quantile(res.full$data.summary$length, c(0.1, 0.9))) +
  ylab("Transcript length (bp)") +
  xlab("") +
  ggtitle("Transcript length") +
  guides(fill="none")

p6.2 <- ggplot(res.full$data.summary, aes(x=sim_match_type, y=length, fill=sim_match_type)) +
  geom_boxplot(alpha=1, outlier.shape = NA) +
  mytheme +
  scale_fill_manual(values=myPalette, name='Stats') +
  scale_y_continuous(limits = quantile(res.full$data.summary$length, c(0.1, 0.9))) +
  ylab("Transcript length (bp)") +
  xlab("") +
  ggtitle("Transcript length") +
  guides(fill="none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

if (nrow(res.min$data.summary) > 0){
  p6.min.1 <- ggplot(res.min$data.summary, aes(x=match_type, y=length, fill=match_type)) +
    geom_boxplot(alpha=1, outlier.shape = NA) +
    mytheme +
    scale_fill_manual(values=myPalette, name='Stats') +
    scale_y_continuous(limits = quantile(res.min$data.summary$length, c(0.1, 0.9))) +
    ylab("Transcript length (bp)") +
    xlab("") +
    ggtitle("Transcript length") +
    guides(fill="none")
  
  p6.min.2 <- ggplot(res.min$data.summary, aes(x=sim_match_type, y=length, fill=sim_match_type)) +
    geom_boxplot(alpha=1, outlier.shape = NA) +
    mytheme +
    scale_fill_manual(values=myPalette, name='Stats') +
    scale_y_continuous(limits = quantile(res.min$data.summary$length, c(0.1, 0.9))) +
    ylab("Transcript length (bp)") +
    xlab("") +
    ggtitle("Transcript length") +
    guides(fill="none") +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
}


# PLOT 7: Expression level
p7.1 <- ggplot(res.full$data.summary, aes(x=match_type, y=sim_counts, fill=match_type)) +
  geom_boxplot(alpha=1, outlier.shape = NA) +
  mytheme +
  scale_fill_manual(values=myPalette, name='Stats') +
  scale_y_continuous(limits = quantile(res.full$data.summary$sim_counts, c(0.1, 0.9))) +
  ylab("Number of simulated reads") +
  xlab("") +
  ggtitle("Expression level") +
  guides(fill="none")

p7.2 <- ggplot(res.full$data.summary, aes(x=sim_match_type, y=sim_counts, fill=sim_match_type)) +
  geom_boxplot(alpha=1, outlier.shape = NA) +
  mytheme +
  scale_fill_manual(values=myPalette, name='Stats') +
  scale_y_continuous(limits = quantile(res.full$data.summary$sim_counts, c(0.1, 0.9))) +
  ylab("Number of simulated reads") +
  xlab("") +
  ggtitle("Expression level") +
  guides(fill="none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

if (nrow(res.min$data.summary) > 0){
  p7.min.1 <- ggplot(res.min$data.summary, aes(x=match_type, y=sim_counts, fill=match_type)) +
    geom_boxplot(alpha=1, outlier.shape = NA) +
    mytheme +
    scale_fill_manual(values=myPalette, name='Stats') +
    scale_y_continuous(limits = quantile(res.min$data.summary$sim_counts, c(0.1, 0.9))) +
    ylab("Number of simulated reads") +
    xlab("") +
    ggtitle("Expression level") +
    guides(fill="none")
  
  p7.min.2 <- ggplot(res.min$data.summary, aes(x=sim_match_type, y=sim_counts, fill=sim_match_type)) +
    geom_boxplot(alpha=1, outlier.shape = NA) +
    mytheme +
    scale_fill_manual(values=myPalette, name='Stats') +
    scale_y_continuous(limits = quantile(res.min$data.summary$sim_counts, c(0.1, 0.9))) +
    ylab("Number of simulated reads") +
    xlab("") +
    ggtitle("Expression level") +
    guides(fill="none") +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
}


# PLOT 8: Transcripts per gene
trans.per.gene <- res.full$data.summary %>%
  group_by(gene_id) %>%
  summarise(trans_per_gene=n())
res.full$data.summary <- merge(res.full$data.summary, trans.per.gene, by = "gene_id")

p8.1 <- ggplot(res.full$data.summary, aes(x=match_type, y=trans_per_gene, fill=match_type)) +
  geom_boxplot(alpha=1, outlier.shape = NA) +
  mytheme +
  scale_fill_manual(values=myPalette, name='Stats') +
  scale_y_continuous(limits = quantile(res.full$data.summary$trans_per_gene, c(0.1, 0.9))) +
  ylab("Number of transcripts per gene") +
  xlab("") +
  ggtitle("Isoform complexity") +
  guides(fill="none")

p8.2 <- ggplot(res.full$data.summary, aes(x=sim_match_type, y=trans_per_gene, fill=sim_match_type)) +
  geom_boxplot(alpha=1, outlier.shape = NA) +
  mytheme +
  scale_fill_manual(values=myPalette, name='Stats') +
  scale_y_continuous(limits = quantile(res.full$data.summary$trans_per_gene, c(0.1, 0.9))) +
  ylab("Number of transcripts per gene") +
  xlab("") +
  ggtitle("Isoform complexity") +
  guides(fill="none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

if (nrow(res.min$data.summary) > 0){
  trans.per.gene <- res.min$data.summary %>%
    group_by(gene_id) %>%
    summarise(trans_per_gene=n())
  res.min$data.summary <- merge(res.min$data.summary, trans.per.gene, by = "gene_id")
  
  p8.min.1 <- ggplot(res.min$data.summary, aes(x=match_type, y=trans_per_gene, fill=match_type)) +
    geom_boxplot(alpha=1, outlier.shape = NA) +
    mytheme +
    scale_fill_manual(values=myPalette, name='Stats') +
    scale_y_continuous(limits = quantile(res.min$data.summary$trans_per_gene, c(0.1, 0.9))) +
    ylab("Number of transcripts per gene") +
    xlab("") +
    ggtitle("Isoform complexity") +
    guides(fill="none")
  
  p8.min.2 <- ggplot(res.min$data.summary, aes(x=sim_match_type, y=trans_per_gene, fill=sim_match_type)) +
    geom_boxplot(alpha=1, outlier.shape = NA) +
    mytheme +
    scale_fill_manual(values=myPalette, name='Stats') +
    scale_y_continuous(limits = quantile(res.min$data.summary$trans_per_gene, c(0.1, 0.9))) +
    ylab("Number of transcripts per gene") +
    xlab("") +
    ggtitle("Isoform complexity") +
    guides(fill="none") +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
}


if ('within_CAGE_peak' %in% colnames(data.index)){
  # PLOT 9: within cage peak
  data.query$match_type <- 'FP'
  data.query$match_type[which(data.query$isoform %in% res.full$novel.perfect.matches$isoform)] <- 'TP_novel'
  data.query$match_type[which(data.query$isoform %in% res.full$known.perfect.matches$isoform)] <- 'TP_known'
  p9.known_FN <- data.index[which(data.index$transcript_id %in% res.full$data.summary$transcript_id[which(res.full$data.summary$sim_match_type == "FN_known")]),]
  if (nrow(p9.known_FN) > 0){
    p9.known_FN$match_type <- 'FN_known'
  } else {
    p9.known_FN$match_type <- character()
  }
  
  p9.novel_FN <- data.index[which(data.index$transcript_id %in% res.full$data.summary$transcript_id[which(res.full$data.summary$sim_match_type == "FN_novel")]),]
  if (nrow(p9.novel_FN) > 0){
    p9.novel_FN$match_type <- 'FN_novel'
  } else {
    p9.novel_FN$match_type <- character()
  }
  p9.all <- rbind(data.query[,c('structural_category', 'match_type', 'within_CAGE_peak', 'dist_to_CAGE_peak')],
                  p9.known_FN[,c('structural_category', 'match_type', 'within_CAGE_peak', 'dist_to_CAGE_peak')],
                  p9.novel_FN[,c('structural_category', 'match_type', 'within_CAGE_peak', 'dist_to_CAGE_peak')])
  
  p9 <- p9.all[which(!is.na(p9.all$within_CAGE_peak)),] %>%
    group_by(match_type, within_CAGE_peak) %>%
    summarise(value=n()) %>%
    ggplot(aes(x=match_type)) +
    geom_bar(aes(y=value, fill=within_CAGE_peak), position="fill", stat="identity") +
    scale_fill_manual(values=myPalette[1:2], name='CagePeak') +
    mytheme +
    ylab('Percentage %') +
    xlab('') +
    ggtitle('Within CAGE peak') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
  
  # PLOT 10: distance to cage peak
  p10 <- p9.all[which(!is.na(p9.all$dist_to_CAGE_peak)),] %>%
    ggplot(aes(x=dist_to_CAGE_peak, color=match_type, fill=match_type)) +
    geom_density(alpha=.3) +
    scale_color_manual(values = myPalette) +
    scale_fill_manual(values = myPalette) +
    mytheme +
    ylab('Distance to CAGE peak') +
    xlab('') +
    ggtitle('Distance To Cage Peak') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
}

if ('min_cov' %in% colnames(data.index)) {
  # PLOT 11: min SJ cov by short reads
  data.query$match_type <- 'FP'
  data.query$match_type[which(data.query$isoform %in% res.full$novel.perfect.matches$isoform)] <- 'TP_novel'
  data.query$match_type[which(data.query$isoform %in% res.full$known.perfect.matches$isoform)] <- 'TP_known'
  p11.known_FN <- data.index[which(data.index$transcript_id %in% res.full$data.summary$transcript_id[which(res.full$data.summary$sim_match_type == "FN_known")]),]
  if (nrow(p11.known_FN) > 0){
    p11.known_FN$match_type <- 'FN_known'
  } else {
    p11.known_FN$match_type <- character()
  }
  p11.novel_FN <- data.index[which(data.index$transcript_id %in% res.full$data.summary$transcript_id[which(res.full$data.summary$sim_match_type == "FN_novel")]),]
  
  if (nrow(p11.novel_FN) > 0){
    p11.novel_FN$match_type <- 'FN_novel'
  } else {
    p11.novel_FN$match_type <- character()
  }
  
  p11.all <- rbind(data.query[,c('structural_category', 'match_type', 'min_cov')],
                  p11.known_FN[,c('structural_category', 'match_type', 'min_cov')],
                  p11.novel_FN[,c('structural_category', 'match_type', 'min_cov')])
  p11.all$Coverage_SJ <- 'False'
  p11.all[which(p11.all$min_cov>0), 'Coverage_SJ'] <- 'True'
  
  p11 <- p11.all[which(!is.na(p11.all$Coverage_SJ)),] %>%
    group_by(match_type, Coverage_SJ) %>%
    summarise(value=n()) %>%
    ggplot(aes(x=match_type)) +
    geom_bar(aes(y=value, fill=Coverage_SJ), position="fill", stat="identity") +
    scale_fill_manual(values=myPalette, name='Coverage SJ') +
    mytheme +
    ylab('Percentage %') +
    xlab('') +
    ggtitle('Splice Junctions Short Reads Coverage') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
}

if ('ratio_TSS' %in% colnames(data.index)) {
  data.query$match_type <- 'FP'
  data.query$match_type[which(data.query$isoform %in% res.full$novel.perfect.matches$isoform)] <- 'TP_novel'
  data.query$match_type[which(data.query$isoform %in% res.full$known.perfect.matches$isoform)] <- 'TP_known'
  p12.known_FN <- data.index[which(data.index$transcript_id %in% res.full$data.summary$transcript_id[which(res.full$data.summary$sim_match_type == "FN_known")]),]
  if (nrow(p12.known_FN) > 0){
    p12.known_FN$match_type <- 'FN_known'
  } else {
    p12.known_FN$match_type <- character()
  }
  p12.novel_FN <- data.index[which(data.index$transcript_id %in% res.full$data.summary$transcript_id[which(res.full$data.summary$sim_match_type == "FN_novel")]),]
  
  if (nrow(p12.novel_FN) > 0){
    p12.novel_FN$match_type <- 'FN_novel'
  } else {
    p12.novel_FN$match_type <- character()
  }
  # PLOT 12: ratio TSS
  p12.all <- rbind(data.query[,c('structural_category', 'match_type', 'ratio_TSS')],
                  p12.known_FN[,c('structural_category', 'match_type', 'ratio_TSS')],
                  p12.novel_FN[,c('structural_category', 'match_type', 'ratio_TSS')])
  
  p12 <- p12.all[which(!is.na(p12.all$ratio_TSS)),] %>%
    ggplot(aes(x=log(ratio_TSS), color=match_type, fill=match_type)) +
    geom_density(alpha=.3) +
    scale_color_manual(values = myPalette) +
    scale_fill_manual(values = myPalette) +
    mytheme +
    ylab('log TSS ratio') +
    xlab('') +
    ggtitle('Ratio TSS') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
}

# Quantification descriptors
if (quant.file != "none"){
  quantification <- read.table(quant.file, header = F, sep="\t")
  colnames(quantification) <- c("transcript_id", "predicted_counts")
  #quantification <- merge(quantification, modif.index, by.x = "transcript_id", by.y = "pipeline_performance", all.x = T)
  quantification <- merge(quantification, modif.index, by.x = "transcript_id", by.y = "pipeline_performance")
  quantification$sim_counts[is.na(quantification$sim_counts)] <- 0
  quantification <- quantification[,c("transcript_id", "predicted_counts", "sim_counts")]
  
  #RMSE(RMSD)
  rmse <- sqrt(mean((quantification$sim_counts - quantification$predicted_counts)^2))
  
  #MAPE
  mape <- mean(abs((quantification$sim_counts-quantification$predicted_counts)/quantification$sim_counts)) * 100

  #Q-Q plot
  nq <- 500
  p <- (1:nq)/nq -0.5/nq
  p13 <- ggplot() +
    geom_point(aes(x=quantile(log(quantification$sim_counts), p), y=quantile(log(quantification$predicted_counts), p)), color = "#15918A") +
    geom_abline(slope=1, intercept=0, color="#F58A53") +
    annotate(geom = 'text', label = paste0("RMSE = ", round(rmse,3), "\nMAPE = ", round(mape, 3), "%"), x = -Inf, y = Inf, hjust = -0.1, vjust = 1) +
    mytheme +
    xlab("Simulated log counts") +
    ylab("Predicted log counts") +
    ggtitle("Q-Q plot of TP count levels")
}

# PLOT X: Radar chart
# Generated in RMD file

# -------------------- Output report
rmarkdown::render(
  input = paste(src.path, 'SQANTI-SIM_report.Rmd', sep = "/"),
  output_dir = output_directory,
  output_file = paste0(output_name, "_SQANTI-SIM_report.html")
)

