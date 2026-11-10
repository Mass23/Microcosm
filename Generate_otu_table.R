library(ggplot2)
library(dplyr)
library(tidyverse)

if (!requireNamespace("phyloseq", quietly = TRUE)) {
  BiocManager::install("phyloseq")
}
if (!requireNamespace("microbiome", quietly = TRUE)) {
  BiocManager::install("microbiome")
}

library(phyloseq)
library(vegan)
library(microbiome)
setwd('work/epfl-altshuler/Microcosm')



# edit here to your own one drive
# Load data ----------------------------------------------------------------

# metadata
meta <- read.csv("DATA/amplicon_analysis/metadata.csv", header = T) 
meta$SWC <- gsub("%", "", meta$SWC) 



# for formatting raw data
# # sample barcodes
# barcodes <- read.csv("sample_barcode.csv", header = T)
# 
# # count table
# counts <- read.csv("coverm_abundance.tsv", sep = "\t", header = T)
# length(counts$Contig) # 26912557 otus
# # 
# rownames(counts) <- counts$Contig
# counts <- counts[-1]
# 
# 
# # format col names otu
# colnames(counts) 
# colnames(counts) <- gsub("all_derep.fasta.", "", colnames(counts))
# colnames(counts) <- gsub("_chopper.fastq.gz.Read.Count", "", colnames(counts))
# 
# 
# # column sample names
# for (i in 1:length(colnames(counts))){
#   colnames(counts)[i] <- barcodes[barcodes$Barcode == colnames(counts)[i], ]$Sample
# }
# 
# head(counts)






# Taxa subset genera only --------------------------------------------------------------------

# # load
# taxa_16s <- read.csv("taxonomy_silvaSSU.tsv", sep = "\t") # use this
# length(taxa_16s$Taxon) # 26912557
# 
# 
# # format
# rownames(taxa_16s) <- taxa_16s$Feature.ID
# taxa_16s <- taxa_16s[-1]
# head(taxa_16s)
# 
# taxa_16s_filt <- taxa_16s[grepl("g__", taxa_16s$Taxon), ]
# length(taxa_16s_filt$Taxon) # 11726726
# 
# # filter counts according to assigned taxa
# counts_taxa <- merge(counts, taxa_16s_filt, by = "row.names", all.y = T)
# head(counts_taxa)
# # write.csv(counts_taxa, "count_filt_genus.csv", row.names = T) # save filtered reads
# 
# # summarise across genus as feature
# counts_genus <- counts_taxa %>% group_by(Taxon) %>% summarise(across(where(is.numeric), sum))
# # write.csv(counts_genus, "taxa_count_filt_genus.csv", row.names = T) # save filtered reads
# # 2152 genera

# load counts
counts_genus <- read.csv("taxa_count_filt_genus.csv", check.names = F, header = T, row.names = 1)

counts_filt <- counts_genus
rownames(counts_filt) <- counts_filt$Taxon
counts_filt <- counts_filt[-c(1, 39)] # remove taxa and confidence cols




# Filter ------------------------------------------------------------------


# do not remove extraction control otus
colSums(counts_filt) # about 20,000 counts in controls

# only remove control samples (not genera)
counts_filt$ext_control_1 = NULL # remove controls
counts_filt$ext_control_2 = NULL
length(counts_filt$'0_A') # 2153 genera, 387 in controls

length(colnames(counts_filt)) # 35 samples
colSums(counts_filt) # all above 20,000




# Taxonomy ----------------------------------------------------------------

# check proportion of otus per demain
length(counts_filt[grepl("Bacteria", rownames(counts_filt)), ]$'0_A') # 1849 - 85%
length(counts_filt[grepl("Archaea", rownames(counts_filt)), ]$'0_A') # 23 - 1%
length(counts_filt[grepl("Eukaryota", rownames(counts_filt)), ]$'0_A') # 281 - 13%

# format taxa
rownames(counts_genus) <- counts_genus$Taxon
taxonomy <- counts_genus %>% dplyr::select(Taxon) %>% 
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")  %>% 
  as.matrix()
taxonomy <- gsub("d__","",
                 gsub("p__","",
                      gsub("o__","",
                           gsub("c__","",
                                gsub("g__","",
                                     gsub("s__","",
                                          gsub("f__","", taxonomy)))))))  %>% 
  as.matrix()




# format for phyloseq object
meta <- meta[match(colnames(counts_filt), meta$Sample), ] # reorder meta samples
rownames(meta) <- meta$Sample
all(rownames(meta) == colnames(counts_filt))

# create phyloseq object
otu <- otu_table(counts_filt, taxa_are_rows = TRUE)
sample_data <- sample_data(meta)
tax_table <- tax_table(taxonomy)
otu_all <- merge_phyloseq(otu, sample_data, tax_table) 
otu_all # 2153 genera


# filter chloroplast, mitochondria
otu_clean <- subset_taxa(otu_all, (Order != "Chloroplast") | is.na(Order))
otu_clean <- subset_taxa(otu_clean, (Family != "Mitochondria") | is.na(Family)) 
otu_clean # 2151 genera

# filter cnidaria (artifact of silva SSU database, these are poorly aligned sequences that get this annotation)
rownames(taxonomy)[grep("Cnidaria", taxonomy[, "Phylum"])] 
# "d__Eukaryota;p__Cnidaria;c__Myxozoa;o__Bivalvulida;f__Bivalvulida;g__Bivalvulida;s__"
otu_clean <- subset_taxa(otu_clean, (Genus != "Bivalvulida")) # 2150 genera


# Rarefaction curve after filtering
otu_clean_t <- t(otu_table(otu_clean))
class(otu_clean_t) <- "matrix"
rarecurve(otu_clean_t, step=1000, ylab="OTU", label=T, cex.axis=1.5, cex.lab=1.5)

# relative abundance
otu_rel <- transform_sample_counts(otu_clean, function(x) x / sum(x))
colSums(otu_table(otu_rel)) # check sum to 1



# Abundance plots ---------------------------------------------------------------


# Read abundance plot
filt_data <- psmelt(otu_clean) # create dataframe from phyloseq object

sum(filt_data$Abundance) # 8046739

# domain count abundance
sum(filt_data[filt_data$Kingdom == "Bacteria", ]$Abundance) / sum(filt_data$Abundance) * 100 # 95.31%
sum(filt_data[filt_data$Kingdom == "Archaea", ]$Abundance) / sum(filt_data$Abundance) * 100 # 1.00%
sum(filt_data[filt_data$Kingdom == "Eukaryota", ]$Abundance) / sum(filt_data$Abundance) * 100 # 28% with cnidaria, 3.69% filtered
sum(filt_data[filt_data$Phylum %in% c("Ascomycota", "Basidiomycota", "Chytridiomycota",
                                      "Cryptomycota", "Mucoromycota", "Blastocladiomycota",
                                      "Neocallimastigomycota", "Zoopagomycota", "Dikarya_Incertae_Sedis"), ]$Abundance) / sum(filt_data$Abundance) * 100 # 0.37% fungi


# kingdoms
# relative abundance plot
kingdom_data <- filt_data %>% group_by(Sample, Day, Treatment, Condition, Group, Kingdom) %>% summarise(relative_abundance = sum(Abundance))

kingdom_med <- kingdom_data %>% group_by(Day, Treatment, Condition, Group, Kingdom) %>% summarise(median_abundance = median(relative_abundance))
kingdom_med$Day <- factor(kingdom_med$Day, levels = c("inoculum", "42_neg", "0", "21", "42", "84", "126"))
kingdom_med$Treatment <- factor(kingdom_med$Treatment, levels = c("Initial", "Aerobic", "Flooded"))

ggplot(kingdom_med[kingdom_med$Condition == "experiment", ], aes(x = Day, y = median_abundance, fill = Kingdom)) + geom_bar(aes(fill = Kingdom), stat="identity", position="fill", colour = "black") + 
  scale_fill_manual(values = wes_palette("Zissou1", 3, type = "continuous")) + theme_bw() + ylab("Median Rel Abundance") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + facet_grid(cols = vars(Condition, Treatment), scales = "free", space = "free") 
#ggsave("kingdom_abundance.png", plot = last_plot(), height = 7, width = 10)






# phyla
# order all phyla
filt_phyla <- aggregate_taxa(otu_clean, level = "Phylum") # phylum level data
filt_data$Phylum <- factor(filt_data$Phylum, levels = top_taxa(filt_phyla)) # order by phylum abundance

# top 10 phyla only
other <- subset(filt_data, !(filt_data$Phylum %in% top_taxa(filt_phyla, 10))) #
other$Phylum <- "Other" # create other category
top_data <- rbind(other, filt_data[filt_data$Phylum %in% top_taxa(filt_phyla, 10), ]) # combine other and top 10 phyla
top_data$Phylum <- factor(top_data$Phylum, levels = c(top_taxa(filt_phyla, 10), "Other")) # order

ggplot(top_data, aes(x = Sample, y = Abundance, fill = Phylum)) + geom_bar(aes(fill = Phylum), stat="identity", position="stack") +
  scale_fill_manual(values = wes_palette("Zissou1", 11, type = "continuous")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggsave("top_phyla_abundance_raw.png", plot = last_plot(), height = 7, width = 10)


# relative abundance plot
top_data_formatted <- top_data %>% group_by(Sample, Day, Treatment, Condition, Group, Phylum) %>% summarise(relative_abundance = sum(Abundance))
ggplot(top_data_formatted, aes(x = Sample, y = relative_abundance, fill = Phylum)) + geom_bar(aes(fill = Phylum), stat="identity", position="fill", colour = "black") +
  scale_fill_manual(values = wes_palette("Zissou1", 11, type = "continuous")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip() + facet_grid(rows = vars(Condition, Day, Treatment), scales = "free", space = "free")

top_data_med <- top_data_formatted %>% group_by(Day, Treatment, Condition, Group, Phylum) %>% summarise(median_abundance = median(relative_abundance))
top_data_med$Day <- factor(top_data_med$Day, levels = c("inoculum", "42_neg", "0", "21", "42", "84", "126"))
top_data_med$Treatment <- factor(top_data_med$Treatment, levels = c("Initial", "Aerobic", "Flooded"))

ggplot(top_data_med, aes(x = Day, y = median_abundance, fill = Phylum)) + geom_bar(aes(fill = Phylum), stat="identity", position="fill", colour = "black") +
  scale_fill_manual(values = wes_palette("Zissou1", 11, type = "continuous")) + theme_bw() + ylab("Median Rel Abundance") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + facet_grid(cols = vars(Condition, Treatment), scales = "free", space = "free")
#ggsave("top_phyla_abundance.png", plot = last_plot(), height = 7, width = 10)
















# # # # # # # 





