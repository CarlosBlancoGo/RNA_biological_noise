# PCA
my_theme <- theme(
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20),
  axis.title.x = element_text(size = 25),
  axis.title.y = element_text(size = 25),
  legend.title = element_text(size = 20),
  legend.text = element_text(size = 18),
  legend.key.size = unit(0.5, 'cm'),
  panel.background = element_rect(fill = "transparent", colour = NA), 
  plot.background = element_rect(fill = "transparent", colour = NA),
  legend.background = element_rect(fill = "transparent", colour = NA),
  plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
)

library("dplyr")
library("data.table")
library("UpSetR")
library("stringr")
library(RColorConesa)
list_classification <- list.files(path="/home/carblan/Documents/R_files/TFM/filtered_classifications", pattern = "B", full.names = TRUE)

# Filenames as 'names' before replacing by actual sample
names <- basename(list_classification)

# Read all dataframes and give them corresponding name
df <- lapply(list_classification, fread)
names(df) <- names

# Add column with sample
for (x in 1:length(df) ) { df[[x]]$sample<-names(df)[x] }

# Select only interesting columns
df <- lapply(df, function(x) {
  select(x, isoform, chrom, length, exons, associated_transcript, associated_gene, structural_category, subcategory, RTS_stage, all_canonical, perc_A_downstream_TTS, filter_result, sample, FL)
})

# rbind all into a big dataframe
df <- data.table::rbindlist(df)

# Give actual sample names
df[grep("B31_", df$sample), "sample"] <- "B31"
df[grep("B32_", df$sample), "sample"] <- "B32"
df[grep("B33_", df$sample), "sample"] <- "B33"
df[grep("B34_", df$sample), "sample"] <- "B34"
df[grep("B35_", df$sample), "sample"] <- "B35"
df[grep("B151_", df$sample), "sample"] <- "B151"
df[grep("B152_", df$sample), "sample"] <- "B152"
df[grep("B153_", df$sample), "sample"] <- "B153"
df[grep("B154_", df$sample), "sample"] <- "B154"

# Order for levels
df$structural_category <- factor(df$structural_category, levels = c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic", "antisense", "fusion", "intergenic", "genic_intron"))
df$sample <- factor(df$sample, levels = c("B31", "B32", "B33", "B34", "B35", "B151", "B152", "B153", "B154"))

# Separate chromosomal transcripts and SIRV transcripts
SIRV<-df[str_detect(df$chrom,"SIRV\\d+$")]
chr_transcripts<-df[!str_detect(df$chrom,"SIRV\\d+$")]


# NIC/NNC proportions between young and old mice.

# Remove novelGenes since they don't match across samples, and also fusion associated genes ('_') because
# they can change their ID order even when running SQANTI in the same sample (they are not useful now)!
# This df shows the number of reads per gene and sample.
library("tidyverse")
library("corrr")
library("ggcorrplot")
library("FactoMineR")
library("factoextra")

test_expression <- chr_transcripts %>%
  filter(!str_detect(associated_gene, "novelGene")) %>%
  filter(!str_detect(associated_gene, "_")) %>%
  filter(!str_detect(filter_result, "Artifact")) %>% 
  group_by(associated_gene, sample) %>%
  summarise(total_FL = sum(FL)) %>%
  ungroup()

genes_to_analyze_expression<-names(which(table(test_expression$associated_gene)==9))
final_df_expression<-test_expression %>% dplyr::filter(associated_gene %in% genes_to_analyze_expression)

# Pivoteo de los datos para tener los genes como filas y las muestras como columnas
formatted_df_expression <- final_df_expression %>%
  pivot_wider(values_from = total_FL, names_from = associated_gene)

# Define un vector con la clasificación de las muestras (young u old)
grupos_muestras <- c("young", "young", "young", "young", "young", "old", "old", "old", "old")
pool_muestras <- c("1", "1", "1", "2", "2", "1", "1", "2", "2")

# Añade la columna de clasificación justo después de la columna 'sample'
formatted_df_expression <- formatted_df_expression %>%
  mutate(group = grupos_muestras, .before = 2) %>%
  mutate(pool = pool_muestras, .before = 2)
formatted_df_expression = as.data.frame(formatted_df_expression)
rownames(formatted_df_expression) = formatted_df_expression$sample

# Y ahora ya la PCA: 
exp.pca <- PCA(formatted_df_expression[, 4:ncol(formatted_df_expression)], graph = FALSE)
eig.val <- get_eigenvalue(exp.pca)
fviz_eig(exp.pca, addlabels = TRUE, ylim = c(0, 50))
var_exp <- get_pca_var(exp.pca)


# Color by age group
fviz_pca_ind(exp.pca,
             col.ind = formatted_df_expression$group,
             palette = colorConesa(2, palette = "main"),
             legend.title = "Age group",
             pointsize = 3,
             mean.point = FALSE,
             repel = TRUE,
             title = "Samples Principal Component Analysis",
             subtitle = "Read counts per gene",
             xlab = paste("PC1 -", round(eig.val[1,2], digits = 2), "%"), 
             ylab = paste("PC2 -", round(eig.val[2,2], digits = 2), "%")
)

# Color by pool
fviz_pca_ind(exp.pca, 
             col.ind = formatted_df_expression$pool,
             palette = colorConesa(2, palette = "main"),
             legend.title = "Pool",
             pointsize = 3,
             mean.point = FALSE,
             repel = TRUE,
             title = "Samples Principal Component Analysis",
             subtitle = "Read counts per gene",
             xlab = paste("PC1 -", round(eig.val[1,2], digits = 2), "%"), 
             ylab = paste("PC2 -", round(eig.val[2,2], digits = 2), "%"),
)



# This df shows the median length of the reads per gene and sample
test_length <- chr_transcripts %>%
  filter(!str_detect(associated_gene, "novelGene")) %>%
  filter(!str_detect(associated_gene, "_")) %>%
  filter(!str_detect(filter_result, "Artifact")) %>% 
  group_by(associated_gene, sample) %>%
  summarise(median_length = median(length)) %>%
  ungroup()

genes_to_analyze_length<-names(which(table(test_length$associated_gene)==9))
final_df_length<-test_length %>% dplyr::filter(associated_gene %in% genes_to_analyze_length)

# Pivoteo de los datos para tener los genes como filas y las muestras como columnas
formatted_df_length <- final_df_length %>%
  pivot_wider(values_from = median_length, names_from = associated_gene)

# Define un vector con la clasificación de las muestras (young u old)
grupos_muestras <- c("young", "young", "young", "young", "young", "old", "old", "old", "old")
pool_muestras <- c("1", "1", "1", "2", "2", "1", "1", "2", "2")

# Añade la columna de clasificación justo después de la columna 'sample'
formatted_df_length <- formatted_df_length %>%
  mutate(group = grupos_muestras, .before = 2) %>%
  mutate(pool = pool_muestras, .before = 2)
formatted_df_length = as.data.frame(formatted_df_length)
rownames(formatted_df_length) = formatted_df_length$sample

# Y ahora ya la PCA: 
exp.pca <- PCA(formatted_df_length[, 4:ncol(formatted_df_length)], graph = FALSE)
eig.val <- get_eigenvalue(exp.pca)
fviz_eig(exp.pca, addlabels = TRUE, ylim = c(0, 50))
var_exp <- get_pca_var(exp.pca)

# Color by age group
fviz_pca_ind(exp.pca,
             col.ind = formatted_df_length$group,
             palette = colorConesa(2, palette = "main"),
             legend.title = "Age group",
             pointsize = 3,
             mean.point = FALSE,
             repel = TRUE,
             title = "Samples Principal Component Analysis",
             subtitle = "Median of transcript models lengths per gene",
             xlab = paste("PC1 -", round(eig.val[1,2], digits = 2), "%"), 
             ylab = paste("PC2 -", round(eig.val[2,2], digits = 2), "%")
)

# Color by pool
fviz_pca_ind(exp.pca, 
             col.ind = formatted_df_length$pool,
             palette = colorConesa(2, palette = "main"),
             legend.title = "Pool",
             pointsize = 3,
             mean.point = FALSE,
             repel = TRUE,
             title = "Samples Principal Component Analysis",
             subtitle = "Median of transcript models lengths per gene",
             xlab = paste("PC1 -", round(eig.val[1,2], digits = 2), "%"), 
             ylab = paste("PC2 -", round(eig.val[2,2], digits = 2), "%"),
)


# This df shows the median number of exons of the reads per gene and sample
test_exons <- chr_transcripts %>%
  filter(!str_detect(associated_gene, "novelGene")) %>%
  filter(!str_detect(associated_gene, "_")) %>%
  filter(!str_detect(filter_result, "Artifact")) %>% 
  group_by(associated_gene, sample) %>%
  summarise(median_exons = median(exons)) %>%
  ungroup()

genes_to_analyze_exons<-names(which(table(test_exons$associated_gene)==9))
final_df_exons<-test_exons %>% dplyr::filter(associated_gene %in% genes_to_analyze_exons)

# Pivoteo de los datos para tener los genes como filas y las muestras como columnas
formatted_df_exons <- final_df_exons %>%
  pivot_wider(values_from = median_exons, names_from = associated_gene)

# Define un vector con la clasificación de las muestras (young u old)
grupos_muestras <- c("young", "young", "young", "young", "young", "old", "old", "old", "old")
pool_muestras <- c("1", "1", "1", "2", "2", "1", "1", "2", "2")

# Añade la columna de clasificación justo después de la columna 'sample'
formatted_df_exons <- formatted_df_exons %>%
  mutate(group = grupos_muestras, .before = 2) %>%
  mutate(pool = pool_muestras, .before = 2)
formatted_df_exons = as.data.frame(formatted_df_exons)
rownames(formatted_df_exons) = formatted_df_exons$sample

# Y ahora ya la PCA: 
exp.pca <- PCA(formatted_df_exons[, 4:ncol(formatted_df_exons)], graph = FALSE)
eig.val <- get_eigenvalue(exp.pca)
fviz_eig(exp.pca, addlabels = TRUE, ylim = c(0, 50))
var_exp <- get_pca_var(exp.pca)

# Color by age group
fviz_pca_ind(exp.pca,
             col.ind = formatted_df_exons$group,
             palette = colorConesa(2, palette = "main"),
             legend.title = "Age group",
             pointsize = 3,
             mean.point = FALSE,
             repel = TRUE,
             title = "Samples Principal Component Analysis",
             subtitle = "Median of exons per gene",
             xlab = paste("PC1 -", round(eig.val[1,2], digits = 2), "%"), 
             ylab = paste("PC2 -", round(eig.val[2,2], digits = 2), "%")
)

# Color by pool
fviz_pca_ind(exp.pca, 
             col.ind = formatted_df_exons$pool,
             palette = colorConesa(2, palette = "main"),
             legend.title = "Pool",
             pointsize = 3,
             mean.point = FALSE,
             repel = TRUE,
             title = "Samples Principal Component Analysis",
             subtitle = "Median of exons per gene",
             xlab = paste("PC1 -", round(eig.val[1,2], digits = 2), "%"), 
             ylab = paste("PC2 -", round(eig.val[2,2], digits = 2), "%"),
)

################################################################################
test_expression <- chr_transcripts %>%
  filter(!str_detect(associated_gene, "novelGene")) %>%
  filter(!str_detect(associated_gene, "_")) %>%
  filter(!str_detect(filter_result, "Artifact")) %>% 
  group_by(associated_gene, sample) %>%
  summarise(total_FL = sum(FL)) %>%
  ungroup()

genes_to_analyze_expression<-names(which(table(test_expression$associated_gene)==9))
final_df_expression<-test_expression %>% dplyr::filter(associated_gene %in% genes_to_analyze_expression)

# Calcular la profundidad de secuenciación por muestra
depth_per_sample <- final_df_expression %>%
  group_by(sample) %>%
  summarize(total_reads = sum(total_FL))

# Unirse al dataframe original para tener la profundidad de secuenciación por muestra en cada fila
final_df_expression <- final_df_expression %>%
  left_join(depth_per_sample, by = "sample")

# Agrupar los datos por gene y calcular los CPM para cada muestra
final_df_expression <- final_df_expression %>%
  group_by(associated_gene) %>%
  mutate(CPM = total_FL / (total_reads / 1e6))
rowSums(formatted_df_expression[,-1])
# Seleccionar las columnas necesarias
final_df_expression <- final_df_expression %>%
  select(associated_gene, sample, CPM)

# Pivoteo de los datos para tener los genes como filas y las muestras como columnas
formatted_df_expression <- final_df_expression %>%
  pivot_wider(values_from = CPM, names_from = associated_gene)

# Define un vector con la clasificación de las muestras (young u old)
grupos_muestras <- c("young", "young", "young", "young", "young", "old", "old", "old", "old")
pool_muestras <- c("1", "1", "1", "2", "2", "1", "1", "2", "2")

# Añade la columna de clasificación justo después de la columna 'sample'
formatted_df_expression <- formatted_df_expression %>%
  mutate(group = grupos_muestras, .before = 2) %>%
  mutate(pool = pool_muestras, .before = 2)
formatted_df_expression = as.data.frame(formatted_df_expression)
rownames(formatted_df_expression) = formatted_df_expression$sample

# Y ahora ya la PCA: 
percentage_df$group_age <- factor(percentage_df$group_age, levels = c("young", "old"))

exp.pca <- PCA(formatted_df_expression[, 4:ncol(formatted_df_expression)], graph = FALSE)
eig.val <- get_eigenvalue(exp.pca)
fviz_eig(exp.pca, addlabels = TRUE, ylim = c(0, 50))
var_exp <- get_pca_var(exp.pca)

formatted_df_expression$group <- factor(formatted_df_expression$group, levels = c("young", "old"))
# Color by age group
p <- fviz_pca_ind(exp.pca,
             col.ind = formatted_df_expression$group,
             repel = TRUE,
             palette = colorConesa(2, palette = "nature"),
             legend.title = "Age group",
             show.legend=FALSE,
             pointsize = 7,
             labelsize = 7,
             mean.point = FALSE,
             xlab = paste("PC1 -", round(eig.val[1,2], digits = 2), "%"), 
             ylab = paste("PC2 -", round(eig.val[2,2], digits = 2), "%"),
)
p$layers[[1]]$data$pool <- factor(formatted_df_expression$pool)
p$layers[[1]]$mapping <- aes(x, y, colour = Col., shape = pool)
p <- p + labs(shape = 'pool') + theme_minimal() + my_theme + labs(title = "")
png(filename = "/home/carblan/Documents/R_files/TFM/figures/PCA/PCA_all_genes.png")
p
dev.off()



###############################################################################
library(tidyr)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(RColorConesa)
# PCA with Illumina reads
illumina_read_counts = read.table("/home/carblan/Documents/R_files/TFM/counts_illumina.tsv")

# Primero transponemos el df para prepararlo para la pca (filas serán las 
# muestras, que serán los individuos, la expresión de genes serán las variables)
illumina_read_counts_transposed <- data.frame(t(illumina_read_counts))
illumina_read_counts_transposed$"__no_feature" <- NULL
illumina_read_counts_transposed$"__ambiguous" <- NULL
illumina_read_counts_transposed$"__too_low_aQual" <- NULL
illumina_read_counts_transposed$"__not_aligned" <- NULL
illumina_read_counts_transposed$"__alignment_not_unique" <- NULL
illumina_read_counts_transposed$"X__no_feature" <- NULL
illumina_read_counts_transposed$"X__ambiguous" <- NULL
illumina_read_counts_transposed$"X__too_low_aQual" <- NULL
illumina_read_counts_transposed$"X__not_aligned" <- NULL
illumina_read_counts_transposed$"X__alignment_not_unique" <- NULL

illumina_read_counts_transposed$total_counts <- rowSums(illumina_read_counts_transposed)

# Calcular CPM
illumina_cpm_df <- illumina_read_counts_transposed[, -ncol(illumina_read_counts_transposed)] / illumina_read_counts_transposed$total_counts * 10^6

# Calcular la varianza de cada gen.
gen_varianzas <- apply(illumina_cpm_df, 2, var)

# Filtrar los genes con varianza igual a 0.
genes_con_varianza_no_cero <- which(gen_varianzas != 0)

# Eliminar estos genes del dataframe.
illumina_cpm_df <- illumina_cpm_df[, genes_con_varianza_no_cero]

# Definir un vector con la clasificación de las muestras (young u old, y por pool)
grupos_muestras <- c("old", "old", "old", "old", "young", "young", "young", "young", "young")
pool_muestras <- c("1", "1", "2", "2", "1", "1", "1", "2", "2")

# Añade las columnas de grupos y pools justo después de la columna 'sample'
illumina_cpm_df <- illumina_cpm_df %>%
  mutate(group = grupos_muestras, .before = 1) %>%
  mutate(pool = pool_muestras, .before = 1)

# Y ahora ya la PCA: 
exp.pca <- PCA(illumina_cpm_df[, 3:ncol(illumina_cpm_df)], graph = FALSE)
eig.val <- get_eigenvalue(exp.pca)
fviz_eig(exp.pca, addlabels = TRUE, ylim = c(0, 50))
var_exp <- get_pca_var(exp.pca)

# Color by age group and different points by pool
p <- fviz_pca_ind(exp.pca,
                  col.ind = illumina_cpm_df$group,
                  repel = TRUE,
                  palette = colorConesa(2, palette = "main"),
                  legend.title = "Age group",
                  show.legend=FALSE,
                  pointsize = 4,
                  mean.point = FALSE,
                  title = "Samples Principal Component Analysis",
                  subtitle = "CPM reads corrected by sequencing depth",
                  xlab = paste("PC1 -", round(eig.val[1,2], digits = 2), "%"), 
                  ylab = paste("PC2 -", round(eig.val[2,2], digits = 2), "%"),
)
p$layers[[1]]$data$pool <- factor(illumina_cpm_df$pool)
p$layers[[1]]$mapping <- aes(x, y, colour = Col., shape = pool)
p <- p + labs(shape = 'pool')
p

setwd("/home/carblan/Documents/R_files/TFM/")
save(illumina_cpm_df, file = "df_illumina_genes_cpm_depth_corrected.RData")
################################################################################
library(tidyr)
library(dplyr)
library("data.table")
library("UpSetR")
library("stringr")
library("corrr")
library("ggcorrplot")
library("FactoMineR")
library("factoextra")
# Now PCAs for the marker genes (Tabula Muris Senis)

# Isoseq Data. 
genes_tsv <- read.table("/home/carblan/Documents/R_files/TFM/DE_genes_brain_15vs3.tsv")
gene_list <- genes_tsv$gene
setwd("/home/carblan/Documents/R_files/TFM/")
load("/home/carblan/Documents/R_files/TFM/df_filtered_Isoseq_genes_cpm_depth_corrected.RData")

formatted_df_marker_genes_expression <- subset(formatted_df_expression, select = names(formatted_df_expression) %in% gene_list)

# Clasificación de las muestras (young u old)
grupos_muestras <- c("young", "young", "young", "young", "young", "old", "old", "old", "old")
pool_muestras <- c("1", "1", "1", "2", "2", "1", "1", "2", "2")

# Añade las columnas de grupos (pool y grupo de edad)
formatted_df_marker_genes_expression <- formatted_df_marker_genes_expression %>%
  mutate(group = grupos_muestras, .before = 1) %>%
  mutate(pool = pool_muestras, .before = 1)

# Y ahora la PCA: 
exp.pca <- PCA(formatted_df_marker_genes_expression[, 3:ncol(formatted_df_marker_genes_expression)], graph = FALSE)
eig.val <- get_eigenvalue(exp.pca)
fviz_eig(exp.pca, addlabels = TRUE, ylim = c(0, 50))
var_exp <- get_pca_var(exp.pca)

formatted_df_marker_genes_expression$group <- factor(formatted_df_expression$group, levels = c("young", "old"))

# Color by age group and different points by pool
p <- fviz_pca_ind(exp.pca,
                  col.ind = formatted_df_marker_genes_expression$group,
                  repel = TRUE,
                  palette = colorConesa(2, palette = "nature"),
                  legend.title = "Age group",
                  show.legend=FALSE,
                  show.title=FALSE,
                  pointsize = 7,
                  labelsize = 7,
                  mean.point = FALSE,
                  xlab = paste("PC1 -", round(eig.val[1,2], digits = 2), "%"), 
                  ylab = paste("PC2 -", round(eig.val[2,2], digits = 2), "%"),
)
p$layers[[1]]$data$pool <- factor(formatted_df_marker_genes_expression$pool)
p$layers[[1]]$mapping <- aes(x, y, colour = Col., shape = pool)
p <- p + labs(shape = 'pool') + theme_minimal() + my_theme + labs(title = "")
png(filename = "/home/carblan/Documents/R_files/TFM/figures/PCA/PCA_ageing_marker_genes.png")
p
dev.off()



# Illumina Data. 
load("/home/carblan/Documents/R_files/TFM/df_illumina_genes_cpm_depth_corrected.RData")

illumina_cpm_marker_genes_df <- subset(illumina_cpm_df, select = names(illumina_cpm_df) %in% gene_list)
# Define un vector con la clasificación de las muestras (young u old)
grupos_muestras <- c("old", "old", "old", "old", "young", "young", "young", "young", "young")
pool_muestras <- c("1", "1", "2", "2", "1", "1", "1", "2", "2")

# Añade la columna de clasificación justo después de la columna 'sample'
illumina_cpm_marker_genes_df <- illumina_cpm_marker_genes_df %>%
  mutate(group = grupos_muestras, .before = 1) %>%
  mutate(pool = pool_muestras, .before = 1)

# Y ahora ya la PCA: 
exp.pca <- PCA(illumina_cpm_marker_genes_df[, 3:ncol(illumina_cpm_marker_genes_df)], graph = FALSE)
eig.val <- get_eigenvalue(exp.pca)
fviz_eig(exp.pca, addlabels = TRUE, ylim = c(0, 50))
var_exp <- get_pca_var(exp.pca)


# Color by age group and different points by pool
p <- fviz_pca_ind(exp.pca,
                  col.ind = illumina_cpm_marker_genes_df$group,
                  repel = TRUE,
                  palette = colorConesa(2, palette = "main"),
                  legend.title = "Age group",
                  show.legend=FALSE,
                  pointsize = 4,
                  mean.point = FALSE,
                  title = "Samples Principal Component Analysis",
                  subtitle = "CPM reads corrected by sequencing depth",
                  xlab = paste("PC1 -", round(eig.val[1,2], digits = 2), "%"), 
                  ylab = paste("PC2 -", round(eig.val[2,2], digits = 2), "%"),
)
p$layers[[1]]$data$pool <- factor(illumina_cpm_marker_genes_df$pool)
p$layers[[1]]$mapping <- aes(x, y, colour = Col., shape = pool)
p <- p + labs(shape = 'pool')
p
