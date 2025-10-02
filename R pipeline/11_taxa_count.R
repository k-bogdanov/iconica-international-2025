# Table 1 - taxa count
# 16S
sample_data(all_16s_phyloseq_rare)$Depth_cm_Country <- with(sample_data(all_16s_phyloseq_rare), 
                                                            paste(Depth_cm, Country, sep = "_"))
Merge_by_TU.phylo <- merge_samples(all_16s_phyloseq_rare, "Depth_cm_Country")

phylum_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Phylum"), add_meta_data = F)
class_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Class"), add_meta_data = F)
order_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Order"), add_meta_data = F)
family_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Family"), add_meta_data = F)
genus_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Genus"), add_meta_data = F)
species_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Species"), add_meta_data = F)

unique_counts <- data.frame(
  Sample = "Total across country",
  Phylum = length(unique(phylum_otu_counts$Phylum)),
  Class = length(unique(class_otu_counts$Class)),
  Order = length(unique(order_otu_counts$Order)),
  Family = length(unique(family_otu_counts$Family)),
  Genus = length(unique(genus_otu_counts$Genus)),
  Species = length(unique(species_otu_counts$Species))
)

taxa_levels_summary<-phylum_otu_counts %>%
  dplyr::count(levels(Phylum),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Phylum")

taxa_levels_summary<-class_otu_counts %>%
  dplyr::count(levels(Class),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Class")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-order_otu_counts %>%
  dplyr::count(levels(Order),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Order")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-family_otu_counts %>%
  dplyr::count(levels(Family),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Family")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-genus_otu_counts %>%
  dplyr::count(levels(Genus),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Genus")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-species_otu_counts %>%
  dplyr::count(levels(Species),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Species")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary <- bind_rows(taxa_levels_summary, unique_counts)

taxa_levels_summary_16S <- taxa_levels_summary %>%
  mutate(Country = case_when(grepl("Ireland", Sample) ~ "Ireland", grepl("Netherlands", Sample) ~ "Netherlands",
                             grepl("Denmark", Sample) ~ "Denmark", grepl("Sweden", Sample) ~ "Sweden"),
         Depth_cm = case_when(grepl("0-10", Sample) ~ "0-10", grepl("10-30", Sample) ~ "10-30", grepl("30-50", Sample) ~ "30-50", TRUE ~ "Total unique taxa"),
         Site_Type = case_when(Country %in% c("Ireland", "Netherlands") ~ "Grassland", Country %in% c("Denmark", "Sweden") ~ "Arable")) %>%
  relocate(c(Site_Type, Country, Depth_cm), .before = Sample) %>%
  select(-Sample) %>%
  mutate(Depth_cm = factor(Depth_cm, levels = c("0-10", "10-30", "30-50", "Total unique taxa"))) %>%
  arrange(desc(Site_Type), Country, Depth_cm)

# ITS
sample_data(all_ITS_phyloseq_rare)$Depth_cm_Country <- with(sample_data(all_ITS_phyloseq_rare), 
                                                            paste(Depth_cm, Country, sep = "_"))
Merge_by_TU.phylo <- merge_samples(all_ITS_phyloseq_rare, "Depth_cm_Country")

phylum_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Phylum"), add_meta_data = F)
class_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Class"), add_meta_data = F)
order_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Order"), add_meta_data = F)
family_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Family"), add_meta_data = F)
genus_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Genus"), add_meta_data = F)
species_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Species"), add_meta_data = F)

unique_counts <- data.frame(
  Sample = "Total across site_type",
  Phylum = length(unique(phylum_otu_counts$Phylum)),
  Class = length(unique(class_otu_counts$Class)),
  Order = length(unique(order_otu_counts$Order)),
  Family = length(unique(family_otu_counts$Family)),
  Genus = length(unique(genus_otu_counts$Genus)),
  Species = length(unique(species_otu_counts$Species))
)

taxa_levels_summary<-phylum_otu_counts %>%
  dplyr::count(levels(Phylum),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Phylum")

taxa_levels_summary<-class_otu_counts %>%
  dplyr::count(levels(Class),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Class")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-order_otu_counts %>%
  dplyr::count(levels(Order),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Order")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-family_otu_counts %>%
  dplyr::count(levels(Family),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Family")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-genus_otu_counts %>%
  dplyr::count(levels(Genus),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Genus")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-species_otu_counts %>%
  dplyr::count(levels(Species),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Species")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary <- bind_rows(taxa_levels_summary, unique_counts)

taxa_levels_summary_ITS <- taxa_levels_summary %>%
  mutate(Country = case_when(grepl("Ireland", Sample) ~ "Ireland", grepl("Netherlands", Sample) ~ "Netherlands",
                             grepl("Denmark", Sample) ~ "Denmark", grepl("Sweden", Sample) ~ "Sweden"),
         Depth_cm = case_when(grepl("0-10", Sample) ~ "0-10", grepl("10-30", Sample) ~ "10-30", grepl("30-50", Sample) ~ "30-50", TRUE ~ "Total unique taxa"),
         Site_Type = case_when(Country %in% c("Ireland", "Netherlands") ~ "Grassland", Country %in% c("Denmark", "Sweden") ~ "Arable")) %>%
  relocate(c(Site_Type, Country, Depth_cm), .before = Sample) %>%
  select(-Sample) %>%
  mutate(Depth_cm = factor(Depth_cm, levels = c("0-10", "10-30", "30-50", "Total unique taxa"))) %>%
  arrange(desc(Site_Type), Country, Depth_cm)
