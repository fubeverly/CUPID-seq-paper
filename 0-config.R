# Custom RProfile for DADA2 analyses
# NatProt 2026

rm(list=ls()) # clear environmental variables

# Load necessary packages
library(here) # current directory
library(dplyr) # dataframe manipulation
library(tidyverse) # transform data long/working form
library(ggplot2) # plotting
library(ggh4x) # ggplot2 extension for facet manipulation
# library(viridis) # color scale
library(phyloseq) # phylogenetic tree
library(readxl) # read excel files
# library(data.table) # manipulating matrices and tables
# library(ggbeeswarm) # beeswarm
# library(ggpubr) # add p-values to figures
# library(rstatix) # wilcoxon test
# library(factoextra) # PCA
# library(MicrobiomeStat) # Diff abundance with LinDA
# library(foreach) # Parallel and iterative processing
# library(ggrepel) # ggplot2 text labels
# library(scales) # scale labels
# library(ggbreak) # break axes
# library(Hmisc) # pearson correlation
# library(ggtext) # custom text
# library(readr)
# library(purrr)
# library(tidyr)
# library(pheatmap)
# library(stringr)
# library(see)

##########################
# Include custom functions
##########################

setup_tables <- function(directory) {
  # Default Setup
  
  # import custom family colors
  colors <- read_excel(paste0(directory, 'Family_colors.xlsx'), sheet = 1)
  colors$Color <- paste0("#", colors$Color)
  colors <- colors %>% mutate(fOTU = paste(Kingdom, Phylum, Class, Order, Family, sep=";"))
  colors <- colors %>% mutate(fOTU = factor(fOTU, levels = rev(unique(colors$fOTU)), ordered = TRUE))
  
  # Create phyloseq object that stores:
  # otu_table <- table of read counts for each ASV
  # tax_table <- taxonomic assignment of ASVs
  # sample_data <- metadata for all experiments
  # phy_tree <- phylogenetic tree from data
  
  # setup tables to populate phyloseq object
  seqtab_nochim <- readRDS(paste0(directory,'seqtab_all.RDS')) # read in ASV read counts
  taxa <- read.table(paste0(directory,'feature_data.txt'), header = TRUE, sep = "\t") # read in ASV taxon info
  taxa <- taxa %>% mutate(fOTU = paste(Kingdom, Phylum, Class, Order, Family, sep=";")) %>% # add fOTU
    mutate(Taxon = mapply(generate_taxon, Kingdom, Phylum, Class, Order, Family, Genus, Species)) %>% # add concise taxonomy name to plot
    mutate(Family2 = ifelse(is.na(Family) | Family == "f__", Taxon, Family)) %>% # add Family2 name to plot
    left_join(colors %>% select(fOTU, Color), by = c("fOTU")) # add Family colors 
  taxa_order <- taxa %>% arrange(match(fOTU, colors$fOTU)) # arrange taxa based on fOTU order in colors
  rownames(taxa) <- taxa$Sequence # row names must match species.names
  taxa <- taxa %>% select(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "fOTU", "Taxon", "Family2", "Color", "Feature"))
  taxa <- as.matrix(taxa) # taxa must be matrix
  
  phy_tree <- read_tree(paste0(directory,'dsvs_msa.tree')) # read in phylogenetic tree
  samdf <- read_excel(paste0(here(), '/1-data/Small-intestine-sample-data.xlsx'), sheet = 1) # Read in sample metadata
  samdf <- as.data.frame(samdf) # convert to dataframe
  rownames(samdf) <- sub('sample_', '', samdf$Sample) # change row names to be the sequencing name, remove "sample_"
  samdf <- samdf[2:ncol(samdf)] # remove Sample column
  
  # create phyloseq object
  ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE), sample_data(samdf), tax_table(taxa), phy_tree(phy_tree))
  
  # store DNA sequences in refseq slot of phyloseq object and rename taxa to short taxa names
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- sprintf("feature_%04d", seq_len(ntaxa(ps)))
  
  return(list(ps = ps, taxa_order = taxa_order, samdf = samdf))
}

generate_taxon <- function(Kingdom, Phylum, Class, Order, Family, Genus, Species) {
  # Check for missing or placeholder values and create Taxon
  if (is.na(Genus) | gsub("g__", "", Genus) == "") {
    if (is.na(Family) | gsub("f__", "", Family) == "") {
      if (is.na(Order) | gsub("o__", "", Order) == "") {
        if (is.na(Class) | gsub("c__", "", Class) == "") {
          if (is.na(Phylum) | gsub("p__", "", Phylum) == "") {
            return(paste("k:", gsub("k__", "", Kingdom)))  # Use Kingdom if Phylum is missing
          } else {
            return(paste("p:", gsub("p__", "", Phylum)))  # Use Phylum if Class is missing
          }
        } else {
          return(paste("c:", gsub("c__", "", Class)))  # Use Class if Order and Family are missing
        }
      } else {
        return(paste("o:", gsub("o__", "", Order)))  # Use Order if Family is missing
      }
    } else {
      return(paste("f:", gsub("f__", "", Family)))  # Use Family if Genus is missing
    }
  } else {
    if (is.na(Species) | gsub("s__", "", Species) == "") {
      return(paste(gsub("g__", "", Genus), "sp."))  # Use Genus and "sp" if Species is missing
    } else {
      return(paste(gsub("g__", "", Genus), gsub("s__", "", Species)))  # Use Genus and Species
    }
  }
}

set_family_order <- function(df2plot, taxa_order) {
  df2plot <- df2plot %>% 
    arrange(match(fOTU, taxa_order$fOTU)) %>% # rearrange row order to fOTU order
    mutate(Family2 = factor(Family2, levels = unique(taxa_order$Family2))) %>% # set new Family2 order
    return(df2plot)
}

set_ordered_factors <- function(df2plot, taxa_order) {
  # Set default ordered factors in dataframe
  df2plot$pH <- factor(df2plot$pH) # Convert pH to a factor (discrete values)
  df2plot$Mouse_no <- factor(df2plot$Mouse_no) # Convert pH to a factor (discrete values)
  df2plot$Buffer <- factor(df2plot$Buffer) # Convert pH to a factor (discrete values)
  
  df2plot <- df2plot %>% 
    mutate(Project = factor(Project, levels = c("c2m5", "c2m7", "m3g1", "m3g2", 
                                                "MouseExp1", "MouseExp2", "MouseExp3", "MouseAvatar", "Klebsiella1", "Klebsiella2", "SIBO"), ordered = TRUE)) %>%
    mutate(Platform = factor(Platform, levels = c("MiSeq", "NovaSeq", "NextSeq"), ordered = TRUE)) %>%
    mutate(Sample_type = factor(Sample_type, levels = c("InVitroCom", "Control", "Sequencing", "Gavage", "Gavage1", "Gavage2", "Isolate", 
                                                        "Saliva", "Capsule", "ProxSI", "MidSI", "DistalSI", "Cecum", "Stool", "Blood"), ordered = TRUE)) %>%
    mutate(Community = factor(Community, levels = c("NA", "Input", "GermFree", "Conventional", "Kleb", 
                                                    "DerivedSI", "DerivedStool_2_DerivedSI", "DerivedSI_DerivedStool", "DerivedStool", "DerivedSI_2_DerivedStool",
                                                    "DirectSI", "DirectSI_DirectStool", "DirectStool", 
                                                    "SynSI", "SynStool_2_SynSI", "SynStool", "SynSI_2_SynStool", "SynGut",
                                                    "SynSI2", "SynSI2_SynStool2", "SynStool2"), ordered = TRUE)) %>%
    mutate(Location = factor(Location, levels = c("NA", "Saliva", "SI", "AC", "Cecum", "Stool", "Gut"), ordered = TRUE)) %>%
    mutate(Donor = factor(Donor, levels = c("NA", "Subj1", "KD30", "Subj2", "Subj3", "Subj4",
                                            "Donor197", "Donor1062", "Donor1064", "Donor1065", "Donor1066", "Donor1067", "Donor1068", "Donor2522",
                                            "KC02", "KC21"), ordered = TRUE)) %>%
    mutate(Pathogen = factor(Pathogen, levels = c("NA", "None", "Kleb"), ordered = TRUE)) %>%
    mutate(Mouse_breed = factor(Mouse_breed, levels = c("NA", "B6", "SW", "SW1", "SW2"), ordered = TRUE)) %>%
    mutate(Mouse_sex = factor(Mouse_sex, levels = c("NA", "Male", "Female"), ordered = TRUE)) %>%
    mutate(Diet = factor(Diet, levels = c("NA", "SD", "MD"), ordered = TRUE)) %>%
    mutate(Mouse_no = factor(Mouse_no, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 
                                                  "10", "11", "12", "13", "14", "15"), ordered = TRUE))
  return(df2plot)
}

get_color_mapping <- function(df2plot){
  # Given a long form data frame, return a mapping of Color (HEX code) to Family2, sorted in Family2 order.
  family2_scale <- unique(df2plot[, c("Color", "Family2")])
  family2_scale$labels <- gsub("f__", "", family2_scale$Family2)
  
  # if any colors do not exist, print the missing family to the console.
  if (anyNA(family2_scale$Color)) { 
    indices <- which(is.na(family2_scale$Color))
    for (i in indices) {
      print(paste("Missing Family color: ", family2_scale$Family2[i]))
    }
  }
  
  names(family2_scale) <- c("values", "breaks", "labels")
  return(family2_scale)
}

# get_evals <- function(pcoa_out){
  # function to get variance percentages
  evals <- pcoa_out$values[,1]
  var_exp <- 100 * evals/sum(evals)
  return(list("evals" = evals, "variance_exp" = var_exp))
}

# split_sample <- function(sample_col) {
  # after merging mice samples, fill in key information
  split_vals <- str_split(sample_col, ';', simplify = TRUE)
  list(
    Project     = split_vals[,1],
    Platform    = split_vals[,2],
    Community   = split_vals[,3],
    Pathogen    = split_vals[,4],
    Diet        = split_vals[,5],
    Mouse_no    = split_vals[,6],
    Location    = split_vals[,7],
    Donor       = split_vals[,8],
    Mouse_breed = split_vals[,9]
  )
}

theme_custom_comp <- function() {
  # Set default plotting method for composition in ggplot2
  theme(
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    legend.position = "none", # no legend
  )
}

theme_custom_text <- function() {
  # set all fontsizes to be 6 and black
  theme(
    legend.text = element_text(size = 6, color = "black"),
    legend.title = element_text(size = 6, color = "black"), 
    axis.text = element_text(size = 6, color = "black"),    
    axis.title = element_text(size = 6, color = "black"),   
    plot.title = element_text(size = 6, color = "black"),   
    strip.text = element_text(size = 6, color = "black")
  )
}

# linda.plot2 <- function(linda.obj, variables.plot, alpha = 0.05, lfc.cut = 1) {
#   # Modification of published LiNDA code. Returns volcano plot information 
#   # so I can custom plot it.
#   bias <- linda.obj$bias
#   output <- linda.obj$output
#   otu.tab <- linda.obj$feature.dat.use
#   meta <- linda.obj$meta.dat.use
#   variables <- linda.obj$variables
#   
#   taxa <- rownames(otu.tab)
#   m <- length(taxa)
#   
#   tmp <- match(variables, variables.plot)
#   voi.ind <- order(tmp)[1 : sum(!is.na(tmp))]
#   padj.mat <- foreach(i = voi.ind, .combine = 'cbind') %do% {
#     output[[i]]$padj
#   }
#   
#   ## volcano plot
#   leg1 <- paste0('padj > ', alpha, ' & ', 'lfc <= ', lfc.cut)
#   leg2 <- paste0('padj > ', alpha, ' & ', 'lfc > ', lfc.cut)
#   leg3 <- paste0('padj <= ', alpha, ' & ', 'lfc <= ', lfc.cut)
#   leg4 <- paste0('padj <= ', alpha, ' & ', 'lfc > ', lfc.cut)
#   
#   gg_color_hue <- function(n) {
#     hues = seq(15, 375, length = n + 1)
#     hcl(h = hues, l = 65, c = 100)[1 : n]
#   }
#   color <- gg_color_hue(3)
#   
#   for(i in 1 : length(voi.ind)) {
#     output.i <- output[[voi.ind[i]]]
#     bias.i <- bias[voi.ind[i]]
#     lfc <- output.i$log2FoldChange
#     padj <- output.i$padj
#     
#     ind1 <- padj > alpha & abs(lfc) <= lfc.cut
#     ind2 <- padj > alpha & abs(lfc) > lfc.cut
#     ind3 <- padj <= alpha & abs(lfc) <= lfc.cut
#     ind4 <- padj <= alpha & abs(lfc) > lfc.cut
#     
#     leg <- rep(NA, m)
#     leg[ind1] = leg1
#     leg[ind2] = leg2
#     leg[ind3] = leg3
#     leg[ind4] = leg4
#     leg <- factor(leg, levels = c(leg1, leg2, leg3, leg4))
#     taxa.sig <- rep('', m)
#     taxa.sig[ind3 | ind4] <- taxa[ind3 | ind4]
#     
#     data.volcano <- cbind.data.frame(taxa = taxa.sig, Log2FoldChange = lfc,
#                                      Log10Padj = -log10(padj), leg = leg)
#   }
#   return(data.volcano)
# }
# 
# diff_taxa <- function(ps_all, proj, diet, commun, pathogen, filter_taxa, 
#                       num_samples = 1, alpha_val = 0.1, lfc_cutoff = 1) {
#   
#   ps2plot <- ps_all %>% subset_samples(Project == proj
#                                        & Final_sac == "TRUE"
#                                        & Location != "Cecum"
#                                        & Diet == diet
#                                        & Pathogen == pathogen
#                                        & Community == commun) %>%
#     merge_samples(group = "Merged_meta")
#   
#   if (!is.null(filter_taxa)) {
#     ps2plot <- prune_taxa(filter_taxa$OTU, ps2plot) %>% # keep ASVs in taxa2keep
#       transform_sample_counts(function(x) x/sum(x)) # renormalize
#   } else {
#     taxa2keep <- psmelt(ps2plot) %>%
#       filter(Abundance > cutoff_SI) 
#     
#     ps2plot <- prune_taxa(taxa2keep$OTU, ps2plot) %>% # keep ASVs in taxa2keep
#       transform_sample_counts(function(x) x/sum(x)) # renormalize
#   }
#   
#   taxa2keep <- psmelt(ps2plot) %>%
#     filter(Abundance > cutoff_SI) %>% 
#     mutate(across(c(Project, Platform, Community, Pathogen, Diet, Mouse_no, Location, Donor), 
#                   ~ split_sample(Sample)[[cur_column()]])) %>% # fill in metadata
#     group_by(Community, Location, Diet, OTU) %>% 
#     dplyr::summarise(num_OTU = n()) %>%
#     filter(num_OTU > num_samples) # Make it for at least three samples
#   
#   ps2plot <- prune_taxa(taxa2keep$OTU, ps2plot) %>% # keep ASVs in taxa2keep
#     transform_sample_counts(function(x) x/sum(x)) # renormalize
#   
#   # for LinDA
#   feature.dat <- t(ps2plot@otu_table)
#   meta.dat <- ps2plot@sam_data
#   meta.dat <- data.frame(sample_data(ps2plot))
#   meta.dat$Sample <- rownames(meta.dat)
#   meta.dat <- meta.dat %>%
#     mutate(across(c(Project, Platform, Community, Pathogen, Diet, Mouse_no, Location, Donor), 
#                   ~ split_sample(Sample)[[cur_column()]])) # fill in metadata
#   
#   linda.obj <- linda(feature.dat = feature.dat, meta.dat = meta.dat, formula = '~Location',
#                      feature.dat.type = 'proportion',
#                      is.winsor = TRUE, outlier.pct = 0.03,
#                      p.adj.method = "BH", alpha = alpha_val)
#   
#   plot.volcano <- linda.plot2(linda.obj, c('Location'), alpha = alpha_val, lfc.cut = lfc_cutoff)
#   return(plot.volcano)
#   
# }
# 
# diff_taxa2 <- function(ps_all, proj, location, commun, pathogen, filter_taxa, 
#                        num_samples = 1, alpha_val = 0.1, lfc_cutoff = 1) {
#   
#   ps2plot <- ps_all %>% subset_samples(Project == proj
#                                        & Final_sac == "TRUE"
#                                        & Location != "Cecum"
#                                        & Location == location
#                                        & Pathogen == pathogen
#                                        & Community == commun) %>%
#     merge_samples(group = "Merged_meta")
#   
#   if (!is.null(filter_taxa)) {
#     ps2plot <- prune_taxa(filter_taxa$OTU, ps2plot) %>% # keep ASVs in taxa2keep
#       transform_sample_counts(function(x) x/sum(x)) # renormalize
#   } else {
#     taxa2keep <- psmelt(ps2plot) %>%
#       filter(Abundance > cutoff_SI) 
#     
#     ps2plot <- prune_taxa(taxa2keep$OTU, ps2plot) %>% # keep ASVs in taxa2keep
#       transform_sample_counts(function(x) x/sum(x)) # renormalize
#   }
#   
#   taxa2keep <- psmelt(ps2plot) %>%
#     filter(Abundance > cutoff_SI) %>% 
#     mutate(across(c(Project, Platform, Community, Pathogen, Diet, Mouse_no, Location, Donor), 
#                   ~ split_sample(Sample)[[cur_column()]])) %>% # fill in metadata
#     group_by(Community, Location, Diet, OTU) %>% 
#     dplyr::summarise(num_OTU = n()) %>%
#     filter(num_OTU > num_samples) # Make it for at least three samples
#   
#   ps2plot <- prune_taxa(taxa2keep$OTU, ps2plot) %>% # keep ASVs in taxa2keep
#     transform_sample_counts(function(x) x/sum(x)) # renormalize
#   
#   # for LinDA
#   feature.dat <- t(ps2plot@otu_table)
#   meta.dat <- ps2plot@sam_data
#   meta.dat <- data.frame(sample_data(ps2plot))
#   meta.dat$Sample <- rownames(meta.dat)
#   meta.dat <- meta.dat %>%
#     mutate(across(c(Project, Platform, Community, Pathogen, Diet, Mouse_no, Location, Donor), 
#                   ~ split_sample(Sample)[[cur_column()]])) # fill in metadata
#   
#   linda.obj <- linda(feature.dat = feature.dat, meta.dat = meta.dat, formula = '~Diet',
#                      feature.dat.type = 'proportion',
#                      is.winsor = TRUE, outlier.pct = 0.03,
#                      p.adj.method = "BH", alpha = alpha_val)
#   
#   plot.volcano <- linda.plot2(linda.obj, c('Diet'), alpha = alpha_val, lfc.cut = lfc_cutoff)
#   return(plot.volcano)
#   
# }
# 
# flattenCorrMatrix <- function(cormat, pmat) {
#   ut <- upper.tri(cormat)
#   data.frame(
#     row = rownames(cormat)[row(cormat)[ut]],
#     column = colnames(cormat)[col(cormat)[ut]],
#     cor = cormat[ut],
#     p = pmat[ut]
#   )
# }