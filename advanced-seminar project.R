library(GEOquery)
library(dplyr)
library(readr)

# load the raw count data 
raw_data <- read.table(gzfile("C:/Users/ragha/Desktop/GSE97668_raw_counts_GRCh38.p13_NCBI (2).tsv.gz"), header = FALSE, sep = "\t") 
raw_data

# load the annotation data
file_path <- "C:/Users/ragha/Desktop/Human.GRCh38.p13.annot.tsv.gz"
annotation_data <- read_tsv(file_path)

# create a mapping of geneID to symbol (get geneID and sample IFsymbol from annotation data)
geneID_to_symbol <- annotation_data %>% select(GeneID, Symbol)
geneID_to_symbol

#########


library(tibble)

# renames the columns in raw_data and sets the first row as the column names
colnames(raw_data) <- raw_data[1, ]
raw_data <- raw_data[-1, ]
colnames(raw_data)[1] <- "GeneID"
raw_data

#this is not needed
raw_data$GeneID <- as.character(raw_data$GeneID)
annotation_data$GeneID <- as.character(annotation_data$GeneID)

# Ensure GeneID columns are the same type (character) for the join
raw_data$GeneID <- as.character(raw_data$GeneID)
geneID_to_symbol$GeneID <- as.character(geneID_to_symbol$GeneID)


# Merge raw data with annotation(geneID_to_symbol) to get gene symbols
merged_data <- raw_data %>%
  inner_join(geneID_to_symbol, by = "GeneID")

# Rearrange columns so that Symbol is the first column
merged_data <- merged_data %>%
  relocate(Symbol, .after = GeneID)

# Remove the GeneID column, if desired so we have the Symbol now
merged_data <- merged_data %>% select(-GeneID)
merged_data

# Convert count columns to numeric, ignoring any potential warnings
merged_data[,-1] <- lapply(merged_data[,-1], as.numeric)

# remove duplicates by summing counts for each gene symbol
merged_data <- merged_data %>%
  group_by(Symbol) %>%
  summarize(across(everything(), sum, na.rm = TRUE))
#merged_data now has the raw data with symbol in fist column V


# Removing genes with low counts across all samples
# Determine a minimum count threshold
min_count <- 10
min_samples <- ncol(merged_data) / 2

# Filter out genes with low counts
filtered_merged_raw_counts <- merged_data[rowSums(merged_data >= min_count) >= min_samples, ]
filtered_merged_raw_counts

# rename the symbol column to geneID
colnames(filtered_merged_raw_counts)[colnames(filtered_merged_raw_counts) == "Symbol"] <- "GeneID"



##

# Load the metadata file 
lines <- readLines(gzfile("C:/Users/ragha/Desktop/GSE97668_series_matrix (1).txt.gz"))

# Extract metadata lines starting with '!Sample_'
metadata_lines <- grep("^!Sample_", lines, value = TRUE)

# Split by tabs and bind as a data frame
metadata_file <- do.call(rbind, strsplit(metadata_lines, "\t"))
metadata_file <- as.data.frame(metadata_file, stringsAsFactors = FALSE)

# View the resulting sample metadata
metadata_file


# Transpose the data frame
metadata_transposed <- t(metadata_file)
metadata_transposed <- as.data.frame(metadata_transposed, stringsAsFactors = FALSE)

# Set the first row as column names
colnames(metadata_transposed) <- metadata_transposed[1, ]
metadata_transposed <- metadata_transposed[-1, ]

# Rename columns as needed
colnames(metadata_transposed)[10:18] <- c("inventory_patient_id", "sex", "cell_type", "diagnosis", "apoe", "expired_age", "pmi", "braak_score")

# Convert relevant columns to factors variables, this is essential in R when preparing data for certain types of analyses, such as differential expression with DESeq2
metadata_transposed$diagnosis <- as.factor(metadata_transposed$diagnosis)
if ("batch" %in% colnames(metadata_transposed)) {
  metadata_transposed$batch <- as.factor(metadata_transposed$batch)
}
metadata_transposed$sex <- as.factor(metadata_transposed$sex)

# Save transformed metadata
write.csv(metadata_transposed, "transformed_metadata_file.csv", row.names = FALSE)
metadata_transposed

# to see where the file got saved in my computer
getwd()

#now:

# Column names of the count table (excluding GeneID)
colnames(filtered_merged_raw_counts)

# Sample identifiers from metadata
colnames(metadata_transposed)
metadata_transposed[["!Sample_geo_accession"]]

# Remove any unnecessary characters (quotes, spaces, etc.) from metadata sample identifiers.
metadata_transposed[["!Sample_geo_accession"]] <- gsub('"', '', metadata_transposed[["!Sample_geo_accession"]])
metadata_transposed[["!Sample_geo_accession"]] <- trimws(metadata_transposed[["!Sample_geo_accession"]])
metadata_transposed[["!Sample_geo_accession"]]

# Use double brackets to reference the column
metadata_transposed[["!Sample_geo_accession"]] <- gsub("old_pattern", "new_pattern", metadata_transposed[["!Sample_geo_accession"]])

# Rename the column in the metadata data frame
colnames(metadata_transposed)[colnames(metadata_transposed) == "!Sample_geo_accession"] <- "Sample_geo_accession"

# Now you can access it easily using $
metadata_transposed$Sample_geo_accession <- gsub("old_pattern", "new_pattern", metadata_transposed$Sample_geo_accession)
print(metadata_transposed$Sample_geo_accession)

metadata_transposed[["!Sample_title"]] <- gsub('"', '', metadata_transposed[["!Sample_title"]])
metadata_transposed[["!Sample_title"]] <- trimws(metadata_transposed[["!Sample_title"]])
metadata_transposed[["!Sample_title"]]

colnames(metadata_transposed)[colnames(metadata_transposed) == "!Sample_title"] <- "Sample_title"


metadata_transposed_new <- metadata_transposed[, c("Sample_title", "Sample_geo_accession")]


###########





# Set 'Sample_geo_accession' as row names
rownames(metadata_transposed_new) <- metadata_transposed_new$Sample_geo_accession

# List of specific 'Sample_geo_accession' values to select
selected_accessions <- c("GSM2574998", "GSM2574999", "GSM2575000", "GSM2575001", "GSM2575002", "GSM2575016", "GSM2575017", "GSM2575018", "GSM2575019", "GSM2575020")

# Filter metadata to include only the rows with matching 'Sample_geo_accession'
filtered_metadata <- metadata_transposed_new[metadata_transposed_new$Sample_geo_accession %in% selected_accessions, ]

# Update the 'Sample_title' column to the desired format
filtered_metadata$Sample_title <- gsub("Controls_d-7_GKH7", "control", filtered_metadata$Sample_title)
filtered_metadata$Sample_title <- gsub("Controls_d-7_GKH13", "control", filtered_metadata$Sample_title)
filtered_metadata$Sample_title <- gsub("Controls_d-7_GKH22", "control", filtered_metadata$Sample_title)
filtered_metadata$Sample_title <- gsub("Controls_d-7_GKH25", "control", filtered_metadata$Sample_title)
filtered_metadata$Sample_title <- gsub("Controls_d-7_GKH28", "control", filtered_metadata$Sample_title)
filtered_metadata$Sample_title <- gsub("Asthmatics_d-7_GKH1", "Asthma", filtered_metadata$Sample_title)
filtered_metadata$Sample_title <- gsub("Asthmatics_d-7_GKH4", "Asthma", filtered_metadata$Sample_title)
filtered_metadata$Sample_title <- gsub("Asthmatics_d-7_GKH10", "Asthma", filtered_metadata$Sample_title)
filtered_metadata$Sample_title <- gsub("Asthmatics_d-7_GKH16", "Asthma", filtered_metadata$Sample_title)
filtered_metadata$Sample_title <- gsub("Asthmatics_d-7_GKH19", "Asthma", filtered_metadata$Sample_title)

filtered_metadata$Sample_title <- gsub("Asthma0", "Asthma", filtered_metadata$Sample_title)

filtered_metadata$Sample_title <- gsub("Asthma6", "Asthma", filtered_metadata$Sample_title)

filtered_metadata$Sample_title <- gsub("Asthma9", "Asthma", filtered_metadata$Sample_title)

# Extract specific columns by name
raw_final <- filtered_merged_raw_counts[, c("GeneID", "GSM2574998", "GSM2574999", "GSM2575000", "GSM2575001", "GSM2575002", "GSM2575016", "GSM2575017","GSM2575018", "GSM2575019", "GSM2575020")]
raw_final

# Store GeneID values as row names
raw_final <- raw_final %>%
  column_to_rownames(var = "GeneID")


# check the dimension to be match between the raw_final and filtered_metadata 
dim(filtered_metadata)
dim(raw_final)

#this is the final preparation 

#filtered_metadata and raw_final are the files we going to use for analysis 


###################################
#MA plot
#volcano plot
#heatmap
#bar plot
#KEGG
#pca



library(DESeq2)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = raw_final,
                              colData = filtered_metadata,
                              design = ~ Sample_title)

dds$Sample_title <- relevel(dds$Sample_title, ref = "control")

dds <- DESeq(dds)

# Results for the condition of interest
res <- results(dds, contrast = c("Sample_title", "Asthma", "control"))

# Order by adjusted p-value
res <- res[order(res$padj), ]
head(res)
res

# Filter for significant genes
sig_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
sig_genes

#data visuaization
#MA Plot: Visualize log2 fold changes versus mean expression
plotMA(res, ylim = c(-5, 5), main = "MA Plot")

library(ggplot2)

# Example with ggsave (if applicable)
ggsave("MA_plot.png", plot = last_plot(), width = 8, height = 6, dpi = 300)



#Volcano Plot
library(ggplot2)

volcano_data <- as.data.frame(res)
volcano_data$significant <- volcano_data$padj < 0.05 & abs(volcano_data$log2FoldChange) > 1

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  scale_color_manual(values = c("gray", "red"))


library(pheatmap)

# Normalize counts
normalized_counts <- counts(dds, normalized = TRUE)

# Subset for significant genes
sig_gene_counts <- normalized_counts[rownames(sig_genes), ]

# Plot heatmap
pheatmap(sig_gene_counts,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         annotation_col = filtered_metadata[, "Sample_title", drop = FALSE])



#Bar Plot of Upregulated and Downregulated Genes

# Calculate the total number of significant genes
total_genes <- num_upregulated + num_downregulated

# Add a percentage column to the data frame
bar_data$Percentage <- (bar_data$Count / total_genes) * 100

# Plot the bar chart with percentages
library(ggplot2)

ggplot(bar_data, aes(x = Regulation, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            vjust = -0.5, size = 5) +  # Add percentage labels above bars
  theme_minimal() +
  labs(title = "Differentially Expressed Genes", 
       y = "Number of Genes", 
       x = "") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme(text = element_text(size = 14))  # Optional: Adjust text size




#Bar Plot for Top Significant Genes
#Plot the log2 fold change for the top 10 most significant genes:

# Get the top 10 most significant genes
top_genes <- head(res[order(res$padj), ], 10)
top_genes_df <- as.data.frame(top_genes)
top_genes_df$Gene <- rownames(top_genes_df)

# Plot the bar chart
ggplot(top_genes_df, aes(x = reorder(Gene, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Significant Genes", x = "Genes", y = "Log2 Fold Change") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "Regulation", labels = c("upregulated", "downregulated"))


library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGREST)
library(BiocManager)


# Extract significant DEGs (padj < 0.05 and log2FoldChange > 1 or < -1)
sig_genes <- rownames(res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ])

# Map gene symbols to Entrez IDs
gene_mapping <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids <- gene_mapping$ENTREZID

# Run KEGG enrichment
kegg_enrichment <- enrichKEGG(gene = entrez_ids, organism = "hsa", pvalueCutoff = 0.05)

head(kegg_enrichment)

######################################################


#Bar Plot of Top Pathways:
barplot(kegg_enrichment, showCategory = 10, title = "Top KEGG Pathways")

#Dot Plot of Pathway Enrichment:
dotplot(kegg_enrichment, showCategory = 10, title = "KEGG Pathway Enrichment")

#save the result
write.csv(as.data.frame(kegg_enrichment), "KEGG_Enrichment_Results.csv")

getwd()

#Performing PCA (Principal Component Analysis) on RNA-seq data requires normalized counts. Here's how you can add a PCA analysis step to your existing workflow:
#Step 1: Extract and Normalize Counts
#You already have normalized counts from your DESeq2 analysis:
# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

#Use the prcomp function on the transpose of normalized counts, ensuring rows correspond to genes and columns to samples:
# Perform PCA on normalized counts
pca_results <- prcomp(t(normalized_counts), scale. = TRUE)

# Create a data frame with PCA scores
pca_data <- as.data.frame(pca_results$x)

# Add metadata
pca_data$Sample_title <- filtered_metadata$Sample_title

#Visualize PCA Results
library(ggplot2)

# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = Sample_title)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA of RNA-seq Data",
    x = paste0("PC1: ", round(summary(pca_results)$importance[2, 1] * 100, 1), "% Variance"),
    y = paste0("PC2: ", round(summary(pca_results)$importance[2, 2] * 100, 1), "% Variance")
  ) +
  theme(text = element_text(size = 14))
#This approach will generate a PCA plot showing the clustering of samples based on their expression profiles

