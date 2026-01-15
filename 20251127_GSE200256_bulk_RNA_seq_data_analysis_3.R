################################################################################
#
# 20251127 GSE200256 Bulk RNA Sequencing data analysis 
#
################################################################################
#test 1
rm(list = ls())

setwd("/storage/research/dbmr_luisierlab/projects/Wilms_tumor/public_data/bulk_RNAseq/SRP358696/Larissa/")



rsem_res <- read.table("/storage/research/dbmr_luisierlab/projects/Wilms_tumor/public_data/bulk_RNAseq/SRP358696/rsem/all.genes.expected_count.results_coding", header=TRUE, sep="\t")
counts <- rsem_res[-1:-27, ] # the first 27 columns in rsem are gene ino
counts[, 28:70] <- round(counts[,28:70])

counts$gene_name_uniq <- make.unique(counts$gene_name)


# ------------------------------------------------------------------------------
# Wrangle data for DESeq2
# ------------------------------------------------------------------------------

# DESeq expects the counts to have gene IDs as row names
head(counts)
rownames(counts) = counts$gene_name_uniq
head(counts)

counts <- counts[, 29:ncol(counts)]
# remove last column
counts <- counts[, -ncol(counts)]


metadata <-read.delim("metadata.tsv",sep = "\t", header = TRUE)
# only metadata of bulk RNA-seq
metadata <- metadata[metadata$seq_unit =="bulk",]

# Remove underscore
metadata$submitter_id <- gsub("_", "", metadata$submitter_id)
rownames(metadata) <- metadata$submitter_id

# check whether columnnames matches rownames
all(colnames(counts) == rownames(metadata)) 
# setdiff(x,y) all elements of x that are not in y
setdiff(colnames(counts), rownames(metadata))
setdiff(rownames(metadata), colnames(counts))
# SJWLM071548D1, there is no count data for this patient

#search in count_total columns
"SJWLM071548D1" %in% colnames(counts)
# not in count data

# remove metadata of patient SJWLM071548D1
metadata <- metadata[rownames(metadata) !="SJWLM071548D1",]
# check again
all(colnames(counts) == rownames(metadata))
setdiff(colnames(counts), rownames(metadata))
setdiff(rownames(metadata), colnames(counts))
# colums and rows iedntical

#correct order?
all(colnames(counts) == rownames(metadata))
# FALSE

#equal names but not same order?
setequal(colnames(counts), rownames(metadata))
#TRUE

# order of metadata 
counts_total <- counts[, match(rownames(metadata), colnames(counts))]

# correct order
all(colnames(counts_total) == rownames(metadata))
#true

#remove diaphragma sample of patient SJWLM071558 (keep kidney sample)
"SJWLM071558D2" %in% colnames(counts_total)
metadata <- metadata[rownames(metadata) !="SJWLM071558D2",]
counts_total <- counts_total[, colnames(counts_total) != "SJWLM071558D2"]

#remove liver patient
"SJWLM071564D1" %in% colnames(counts_total)
metadata <- metadata[rownames(metadata) !="SJWLM071564D1",]
counts_total <- counts_total[, colnames(counts_total) != "SJWLM071564D1"]

table(metadata$participant_id)

#remove one sample for participants SJWLM051023 SJWLM051026,as they have 2 samples

"SJWLM051023D2" %in% colnames(counts_total)
metadata <- metadata[rownames(metadata) !="SJWLM051023D2",]
counts_total <- counts_total[, colnames(counts_total) != "SJWLM051023D2"]

"SJWLM051026D2" %in% colnames(counts_total)
metadata <- metadata[rownames(metadata) !="SJWLM051026D2",]
counts_total <- counts_total[, colnames(counts_total) != "SJWLM051026D2"]

table(metadata$participant_id)
#per patient one sample of kidney

#-------------------------------------------------------------------------------
# Metadata analysis
#-------------------------------------------------------------------------------

# How many participants?
print(length(unique(metadata$participant_id)))

# Subdiagnosis
table(unique(metadata[, c("participant_id", "subdiagnosis")])$subdiagnosis)
# Anaplastic 17, Favorable 21


library(dplyr)

# Metastasis yes/no?
metadata %>%
  filter(subdiagnosis == "Favorable") %>%
  count(metastasis)
# 3 metastasis, 18 NA

metadata %>%
  filter(subdiagnosis == "Anaplastic") %>%
  count(metastasis)
# 6 metastasis, 11 NA

#Metastasen, Relpase
metadata %>%
  filter(subdiagnosis == "Favorable") %>%
  count(metastasis, relapse_status)
# 3 metastasis and relapse, 18 no metastasis, no relapse

metadata %>%
  filter(subdiagnosis == "Anaplastic") %>%
  count(metastasis, relapse_status)
# 6 metastasis and relapse, 11 no metastasis, no relapse

# treatment
table(unique(metadata[, c("participant_id", "treatment")])$treatment)
# Resection post chemotherapy 10, Upfront resection 28

# 
metadata %>%
  filter(treatment == "Resection post chemotherapy") %>%
  count(subdiagnosis,metastasis)

metadata %>%
  filter(treatment == "Upfront resection") %>%
  count(subdiagnosis,metastasis)
# participants upfront resection, 3 anaplastic + metastasis, 6 anaplastic 
#no metastasis, 1 favorable no metastasis

metadata %>%
  filter(treatment == "Resection post chemotherapy") %>%
  count(subdiagnosis,metastasis, relapse_status)
# participants Resection post chemotherapy, 3 anaplastic + metastasis + relapse,
# 6 anaplastic no metastasis no relapse, 1 favorable no metastasis no relapse


metadata %>%
  filter(treatment == "Upfront resection") %>%
  count(subdiagnosis,metastasis, relapse_status)
# participants upfront resection, 3 anaplastic + metastasis + relapse, 
# 5 anaplastic no metastasis no relapse, 3 favorable metastasis + relapse,
# 17 favorable no metastasis no relapse

library(DESeq2)

#-------------------------------------------------------------------------------
# Favorable vs. Anaplastic
#-------------------------------------------------------------------------------

metadata$subdiagnosis <- factor(metadata$subdiagnosis, levels = c("Favorable",
                                                            "Anaplastic"))
levels(metadata$subdiagnosis)


dds <-DESeqDataSetFromMatrix(countData= counts_total, colData = metadata, design = ~ subdiagnosis)

#anaplastic vs. favorable
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
head(res)
################################################################################
# PCA plot
AF <- vst(dds, blind = FALSE)
plotPCA(AF, intgroup = "subdiagnosis")
# no clear separation

#rlog
AF_rlog <- rlog(dds, blind = FALSE)
plotPCA(AF_rlog, intgroup = "subdiagnosis")
# outlier
p <- plotPCA(AF_rlog, intgroup = "subdiagnosis")
p$data[which(p$data$PC2 < -100), ]   # oder > 50 je nachdem

#check library size for outlier
colSums(counts_total)[ "SJWLM049191D2" ]
summary(colSums(counts_total))
# not an outlier because of library size

#Quality check
#distance heatmap
library(pheatmap)
r <- rlog(dds, blind = FALSE)
d <- dist(t(assay(r)))
pheatmap(as.matrix(d), annotation_col = metadata)

summary(metadata)
sapply(metadata, function(x) length(unique(x)))

r <- rlog(dds, blind = FALSE)
d <- dist(t(assay(r)))

ann_cols <- metadata[, sapply(metadata, function(x) length(unique(x)) > 1), drop = FALSE]

pheatmap(as.matrix(d), annotation_col = ann_cols)
# SJWLM049191D2 big outlier

#remove outlier
dds_no <- dds[, colnames(dds) != "SJWLM049191D2" ]
dds_no <- DESeq(dds_no)

#PCA without outlier
AF_no <- vst(dds_no, blind = FALSE)
plotPCA(AF_no, intgroup = "subdiagnosis")

#PCA rlog
AF_no_rlog <- rlog(dds_no, blind = FALSE)
plotPCA(AF_no_rlog, intgroup = "subdiagnosis")
#no clear separation

# Hochreguliert in Anaplastic (padj < 0.05 und log2FC > 1)
anaplastic_high <- subset(res, padj < 0.05 & log2FoldChange > 1)

#sort log2FoldChange
anaplastic_high_sorted <- anaplastic_high[order(anaplastic_high$log2FoldChange, decreasing = TRUE), ]
top20_anaplastic <- head(anaplastic_high_sorted, 20)
top20_anaplastic
top20_anaplastic_high_dataset<- as.data.frame(top20_anaplastic)


# downregulated in Anaplastic (padj < 0.05 und log2FC < -1)
anaplastic_down <- subset(res, padj < 0.05 & log2FoldChange < -1)

#sort log2FoldChange
anaplastic_down_sorted <- anaplastic_down[order(anaplastic_down$log2FoldChange, decreasing = FALSE),]
top20_anaplastic_down <- head(anaplastic_down_sorted, 20)
top20_anaplastic_down
top20_anaplastic_down_dataset <- as.data.frame(top20_anaplastic_down)

#-------------------------------------------------------------------------------
# Metastasis vs. no metastasis
#-------------------------------------------------------------------------------
metadata$metastasis[is.na(metadata$metastasis)] <- "unknown"


metadata$metastasis <- factor(metadata$metastasis)
levels(metadata$metastasis)
#need a new columns
metadata$metastasis_status <- ifelse(metadata$metastasis == "unknown",
                                     "Nein",
                                     "Ja")


metadata$metastasis_status <- factor(metadata$metastasis_status,
                                     levels = c("Nein", "Ja"))




dds_m <-DESeqDataSetFromMatrix(countData= counts_total, colData = metadata, design = ~ metastasis_status)


dds_m <- DESeq(dds_m)
resultsNames(dds_m)
res_m <- results(dds_m)
head(res_m)


# PCA plot
M <- vst(dds_m, blind = FALSE)
plotPCA(M, intgroup = "metastasis_status")
# no clear separation


# Hochreguliert in metastasis (padj < 0.05 und log2FC > 1)
metastasis_high <- subset(res_m, padj < 0.05 & log2FoldChange > 1)

#sort log2FoldChange
metastasis_high_sorted <- metastasis_high[order(metastasis_high$log2FoldChange, decreasing = TRUE), ]
top20_metastasis <- head(metastasis_high_sorted, 20)
top20_metastasis
top20_metastasis_high_dataset<- as.data.frame(top20_metastasis)

# downregulated in Anaplastic (padj < 0.05 und log2FC < -1)
metastasis_down <- subset(res_m, padj < 0.05 & log2FoldChange < -1)

#sort log2FoldChange
metastasis_down_sorted <- metastasis_down[order(metastasis_down$log2FoldChange, decreasing = FALSE),]
top20_metastasis_down <- head(metastasis_down_sorted, 20)
top20_metastasis_down
top20_metastasis_down_dataset <- as.data.frame(top20_metastasis_down)


#-------------------------------------------------------------------------------
# Treatment 
#-------------------------------------------------------------------------------

metadata$treatment <- factor(metadata$treatment, levels = c("Upfront resection",
                                                            "Resection post chemotherapy"))

dds_t <-DESeqDataSetFromMatrix(countData= counts_total, colData = metadata, design = ~ treatment)


dds_t <- DESeq(dds_t)
resultsNames(dds_t)
res_t <- results(dds_t)
head(res_t)

# Hochreguliert in metastasis (padj < 0.05 und log2FC > 1)
treatment_high <- subset(res_t, padj < 0.05 & log2FoldChange > 1)

#sort log2FoldChange
treatment_high_sorted <- treatment_high[order(treatment_high$log2FoldChange, decreasing = TRUE), ]
top20_treatment <- head(treatment_high_sorted, 20)
top20_treatment
top20_treatment_high_dataset<- as.data.frame(top20_treatment)


# downregulated in Anaplastic (padj < 0.05 und log2FC < -1)
treatment_down <- subset(res_t, padj < 0.05 & log2FoldChange < -1)

#sort log2FoldChange
treatment_down_sorted <- treatment_down[order(treatment_down$log2FoldChange, decreasing = FALSE),]
top20_treatment_down <- head(treatment_down_sorted, 20)
top20_treatment_down
top20_treatment_down_dataset <- as.data.frame(top20_treatment_down)

#-------------------------------------------------------------------------------
# Treatment, subdiagnosis
#-------------------------------------------------------------------------------

#check whether variables are factors
levels(metadata$subdiagnosis)
levels(metadata$treatment)

metadata$Group <- factor(paste(metadata$subdiagnosis,
                               metadata$treatment,
                               sep = "_"))

levels(metadata$Group)
dds_group <- DESeqDataSetFromMatrix(countData = counts_total,
                                    colData = metadata,
                                    design = ~ Group)
dds_group <- DESeq(dds_group)

#anaplastic Resection post chemotherapy vs. anaplastic Ufront resection
res_anaplastic <- results(dds_group,
                          contrast = c("Group",
                                       "Anaplastic_Resection post chemotherapy",
                                       "Anaplastic_Upfront resection"))


# in data.frame umwandeln
df <- as.data.frame(res_anaplastic)

# Filter setzen
df_filt <- df[df$log2FoldChange > 1 & df$pvalue < 0.05, ]

# Nach LFC sortieren (größte zuerst)
df_filt <- df_filt[order(df_filt$log2FoldChange, decreasing = TRUE), ]

# Top 20 Gene auswählen
top20 <- head(df_filt, 20)

top20

#-------------------------------------------------------------------------------
# Favorable Resection post chemotherapy vs. Fav Upfront surgery
#-------------------------------------------------------------------------------


res_f_c <- results(dds_group,
                          contrast = c("Group",
                                       "Favorable_Resection post chemotherapy",
                                       "Favorable_Upfront resection"))


# in data.frame umwandeln
df_fc <- as.data.frame(res_f_c)

# Filter setzen
df_fc_f <- df_fc[df_fc$log2FoldChange > 1 & df$pvalue < 0.05, ]

# Nach LFC sortieren (größte zuerst)
df_filt_fc <- df_fc_f[order(df_fc_f$log2FoldChange, decreasing = TRUE), ]

# Top 20 Gene auswählen
top20_fc <- head(df_filt_fc, 20)

top20_fc

#-------------------------------------------------------------------------------
# Anaplastic Upfront surgery vs. Fav Upfront surgery
#-------------------------------------------------------------------------------


res_a_u <- results(dds_group,
                   contrast = c("Group",
                                "Anaplastic_Upfront resection",
                                "Favorable_Upfront resection"))


# in data.frame umwandeln
df_au <- as.data.frame(res_a_u)

# Filter setzen
df_au_f <- df_au[df_au$log2FoldChange > 1 & df$pvalue < 0.05, ]

# Nach LFC sortieren (größte zuerst)
df_filt_au <- df_au_f[order(df_au_f$log2FoldChange, decreasing = TRUE), ]

# Top 20 Gene auswählen
top20_au <- head(df_filt_au, 20)

top20_au

#-------------------------------------------------------------------------------
# Subdiagnosis, treatment PCA
#-------------------------------------------------------------------------------

vsd <- vst(dds_group, blind = FALSE)
plotPCA(vsd, intgroup = "Group")

#-------------------------------------------------------------------------------
# Metastasis PCA
#-------------------------------------------------------------------------------

vsd_m <- vst(dds_m, blind = FALSE)
plotPCA(vsd_m, intgroup = "metastasis_status")




#-------------------------------------------------------------------------------
# Metastasis, Subdiagnosis PCA
#-------------------------------------------------------------------------------

#check whether variables are factors
levels(metadata$subdiagnosis)
levels(metadata$metastasis_status)

metadata$Group_MS <- factor(paste(metadata$subdiagnosis,
                               metadata$metastasis_status,
                               sep = "_"))

levels(metadata$Group_MS)
dds_Group_MS <- DESeqDataSetFromMatrix(countData = counts_total,
                                    colData = metadata,
                                    design = ~ Group_MS)
dds_Group_MS <- DESeq(dds_Group_MS)


vsd_MS <- vst(dds_Group_MS, blind = FALSE)
plotPCA(vsd_MS, intgroup = "Group_MS")


#-------------------------------------------------------------------------------
# Metastasis, Subdiagnosis, treatment PCA
#-------------------------------------------------------------------------------

#check whether variables are factors
levels(metadata$subdiagnosis)
levels(metadata$metastasis_status)

metadata$Group_MST <- factor(paste(metadata$Group_MS,
                                  metadata$treatment,
                                  sep = "_"))

levels(metadata$Group_MST)
dds_Group_MST <- DESeqDataSetFromMatrix(countData = counts_total,
                                       colData = metadata,
                                       design = ~ Group_MST)
dds_Group_MST <- DESeq(dds_Group_MST)


vsd_MST <- vst(dds_Group_MST, blind = FALSE)
plotPCA(vsd_MST, intgroup = "Group_MST")

#16:36

