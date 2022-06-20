# ---
# title: "Generalized Decontamination Script"
# author: "Zack Gold"
# date: "01/16/2022"
# output: html_document
# ---

# Load Libraries ----------------------------------------------------------------------------------------------------------------
library (tidyverse)
library (vegan)
library (proxy)
library(reshape2)
library(microDecon)
library(stringr)
library(knitr)
library(ranacapa)
library(dplyr)
library(tibble)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(plotly)
library(optparse)
library(fitdistrplus)
library(rstan)
library(shinystan)
library(bayesplot)
library(broom)
library(analyze.stuff)

source("/Users/zackgold/Documents/UCLA_phd/Projects/anacapa/decontam/decontamination_utilities.R")
library(here)

#Test inputs ----------------------------------------------------------------------------------------------------------------
working_directory <- here("data","Anacapa_output_old","decontamination")
input_anacapa_path_1 <- here("data","Anacapa_output_old","12S_E","12S_taxonomy_tables","Summary_by_percent_confidence","60","12S_ASV_raw_taxonomy_60.txt")
input_anacapa_path_2 <- here("data","Anacapa_output_old","12S_U","12S_taxonomy_tables","Summary_by_percent_confidence","60","12S_ASV_raw_taxonomy_60.txt")

input_meta_path <- here("data","la_ports_meta_data.txt")
number_anacapa_tables <- "two"

read_type <- "merged_only"

low_dist_probability_cutoff <- 0.05
minimum_read_cutoff <- 5000

step3_runs <- 2
step3_thresh <- 0.7
step3_prop.thresh <- 5e-05
step3_regression <- 0

som_level <- "single"
som_chains<-10
som_iterations<-1000
som_filter_choice <- "som_min_threshold_1"

max_dist_probability_cutoff <- 0.9999
#---

#---

#Load Data ----------------------------------------------------------------------------------------------------------------
setwd(working_directory)
dir.create("Output_csv")
dir.create("Output_R")
dir.create("Output_plots")

metadata <- read.table(input_meta_path, header = TRUE, sep = "\t", stringsAsFactors = F)

anacapa_table_1 <- read.table(input_anacapa_path_1, header = 1, sep = "\t", stringsAsFactors = F)
anacapa_table_2 <- read.table(input_anacapa_path_2, header = 1, sep = "\t", stringsAsFactors = F)

anacapa_table_1 %>% 
  dplyr::select(-`X12S_LA11.4.R3.S36.L001`,
                -`X12S_LA3.4.R3.S12.L001`,
                -`X12S_LB1.4.R3.S60.L001`) -> anacapa_table_1

#Generate Hash.key

colnames(anacapa_table_1) %>%
  str_remove("X12S_") -> new_colnames_1

colnames(anacapa_table_1) <- new_colnames_1

anacapa_table_1 %>% 
  mutate(., seq_number = str_replace(seq_number, "12S_","12S_elas_")) -> anacapa_table_1

anacapa_table_1 %>% 
  dplyr::select(seq_number, sum.taxonomy) -> hash.key_1

colnames(anacapa_table_2) %>%
  str_remove("X12S_") -> new_colnames_2

colnames(anacapa_table_2) <- new_colnames_2

anacapa_table_2 %>% 
  mutate(., seq_number = str_replace(seq_number,  "12S_","12S_miu_")) ->anacapa_table_2

anacapa_table_2%>% 
  dplyr::select(seq_number, sum.taxonomy) -> hash.key_2


rbind(hash.key_1,hash.key_2) -> hash.key

saveRDS(hash.key, file="hash.key.RDS")

#Format Anacapa Table 1
anacapa_table_1 %>%
  filter(., str_detect(seq_number,"merged")) %>% 
  dplyr::select(seq_number) -> merged_hash_1
merged_hash_1$seq_number[[1]] %>% str_remove(.,"merged_") %>% str_sub(., end=-3) -> barcode_1
anacapa_table_1$Miseq_run <- barcode_1

###Fix Names
anacapa_table_1_names <- colnames(anacapa_table_1)
anacapa_table_1_names %>% as.data.frame() -> anacapa_table_1_names
colnames(anacapa_table_1_names) <- c("Seq_number")

left_join(anacapa_table_1_names, metadata) %>% 
  dplyr::select(Seq_number, New_name) %>% distinct() %>% 
  mutate(final_names = coalesce(New_name, Seq_number)) %>%  dplyr::select(final_names) %>% as.list() -> anacapa_table_1_names
colnames(anacapa_table_1) <- anacapa_table_1_names$final_names

#Format Anacapa Table 2
anacapa_table_2 %>%
  filter(., str_detect(seq_number,"merged")) %>% 
  dplyr::select(seq_number) -> merged_hash_2
merged_hash_2$seq_number[[1]] %>% str_remove(.,"merged_") %>% str_sub(., end=-3) -> barcode_2
anacapa_table_2$Miseq_run <- barcode_2

####Fix Names
anacapa_table_2_names <- colnames(anacapa_table_2)
anacapa_table_2_names %>% as.data.frame() -> anacapa_table_2_names
colnames(anacapa_table_2_names) <- c("Seq_number")
left_join(anacapa_table_2_names, metadata) %>%
  dplyr::select(Seq_number, New_name) %>% 
  mutate(final_names = coalesce(New_name, Seq_number)) %>%  dplyr::select(final_names) %>% as.list() -> anacapa_table_2_names
colnames(anacapa_table_2) <- anacapa_table_2_names$final_names

#Merge All Tables
ASV.table <- bind_rows(anacapa_table_1,anacapa_table_2)

ASV.table %>% 
  replace(is.na(.), 0) ->ASV.table

tail(ASV.table)
head(ASV.table)


# Convert to Long Data

###Format for Long Data
ASV.table$seq_number <- factor(ASV.table$seq_number)
ASV.table$Miseq_run <- factor(ASV.table$Miseq_run)

columns <- colnames(ASV.table)
remove <- c("seq_number","sum.taxonomy","Miseq_run")

gathercols <-  columns[! columns %in% remove] 

ASV.table %>% 
  pivot_longer(., cols=gathercols, names_to="sample", values_to="reads") -> ASV.table_long

ASV.table_long$reads <- as.numeric(ASV.table_long$reads)

saveRDS(ASV.table_long,  file="ASV.table_long")

ASV.table_long %>%
  filter(., !sample %in% c("LB1_0_1","LB1_0_2","LB1_0_3")) %>% 
  filter(., !str_detect(sample,"Blank")) %>% 
  dplyr::summarise(n_distinct(seq_number), n_distinct(sum.taxonomy), n_distinct(sample), n_distinct(Miseq_run), sum(reads))

#---

#---

# Cleaning Process 0: Remove all Forward, Reverse, and Unmerged reads & Remove Singletons ----------------------------------------------------------------------------------------------------------------

#Filter Merged only reads
ASV.table_long <- readRDS(file="ASV.table_long")

if(read_type =="merged_only"){
  #Filter out forward, reverse, and unmerged reads
  ASV.table_long %>% 
    filter(., str_detect(seq_number,"merged")) %>% 
    filter(., !str_detect(seq_number,"unmerged")) -> ASV.table_used
} else if ( read_type =="merged_and_forward"){
  #Remove Singletons
  ASV.table_long %>% 
    filter(., !str_detect(seq_number,"reverse")) %>% 
    filter(., !str_detect(seq_number,"unmerged")) -> ASV.table_used
} else {
  #Remove Singletons
  ASV.table_long %>%
    dplyr::group_by(seq_number) %>%
    mutate (TotalReadsperSample = sum(reads)) %>% 
    filter(., TotalReadsperSample > 1) %>% 
    dplyr::select(-TotalReadsperSample) -> ASV.table_used
}

#Calculate % ASVs Kept
ASV.table_long %>%  dim() -> all_dim
ASV.table_used %>%  dim() -> used_only_dim

paste0(round(used_only_dim[[1]]/all_dim[[1]]*100,2),"% ASVs retained")

ASV.table_used %>% 
  ungroup() %>% 
  filter(., sum.taxonomy !="" ) %>% 
  left_join(metadata, by =c("sample"="New_name")) %>% 
  filter(., Bio_rep!=0) %>% 
  dplyr::group_by(Miseq_run) %>% 
  dplyr::summarise(ASVs = n_distinct(seq_number), Tot_reads = sum(reads), n_distinct(sample))
#---

#---

# Cleaning Process 1: Estimation of *Tag-jumping* or sample *Cross-talk* ----------------------------------------------------------------------------------------------------------------

## Step 1: Nest the dataset by origin of ASVs

###Identify Positives, Negatives, and Samples

###Create list of control samples
metadata %>% 
  filter(Sample_control=="Control") %>% 
  dplyr::select(New_name) %>% unique() -> controls
controls <- controls$New_name

metadata %>% 
  filter(Control_Type=="Pos") %>% 
  dplyr::select(New_name) %>% unique()-> pos_controls
pos_controls <- pos_controls$New_name

metadata %>% 
  filter(Control_Type=="Blank") %>% 
  dplyr::select(New_name) %>% unique() -> neg_controls
neg_controls <- neg_controls$New_name

###New column that labels each ASV as from Positive (control) or Sample
ASV.table_used %>% 
  mutate(source = case_when(sample %in% pos_controls~"Positives",
                            sample %in% neg_controls~"Blanks",
                            TRUE ~"Samples")) -> ASV.table_used


###Convert to tibble
ASV.table_used <- as_tibble(ASV.table_used)

###Remove empty sequences
ASV.table_used %>% 
  filter(reads != 0)  -> ASV.table_used

###Rename Columns and remove seq_number
ASV.table_used %>%
  mutate(sample = as.character(sample),
         nReads = reads)  -> ASV.table_used

###ASVs in Positive Controls
ASV.table_used %>% 
  filter (source == "Positives") %>%
  dplyr::group_by(seq_number) %>% 
  dplyr::summarise(tot = sum(reads)) %>% 
  arrange(desc(tot)) %>% 
  pull(seq_number) -> all.seqs.in.positives

hash.key %>% 
  filter(seq_number %in% all.seqs.in.positives) %>% as_tibble() -> pos.contam.species

write.csv(pos.contam.species, file="pos.contam.species.csv")

###ASVs in Negative Controls

ASV.table_used %>% 
  filter (source == "Blanks") %>%
  dplyr::group_by(seq_number) %>% 
  dplyr::summarise(tot = sum(reads)) %>% 
  arrange(desc(tot)) %>% 
  pull(seq_number) -> all.seqs.in.blanks

ASV.table_used %>% 
  filter (source == "Blanks") %>%
  dplyr::group_by(sample,sum.taxonomy) %>%
  dplyr::summarise(tot = sum(reads)) %>% 
  pivot_wider(names_from = sample, values_from=tot, values_fill=0) %>% View()

hash.key %>% 
  filter(seq_number %in% all.seqs.in.blanks) %>% as_tibble() -> blank.contam.species

write.csv(blank.contam.species, file="blank.contam.species.csv")

### Visualize Read Counts Across Samples for Barcode_1
ASV.table_used %>% 
  group_by(sample) %>%
  filter(., Miseq_run== barcode_1) %>% 
  mutate (TotalReadsperSample = sum(nReads)) %>%
  arrange(desc(TotalReadsperSample)) %>%
  ggplot(., aes(x=sample, y=TotalReadsperSample, color=source)) + 
  geom_point() +ggtitle("Read Count Across Samples") + 
  theme(axis.text.x = element_text(angle = 90)) -> plot_1

ggsave(plot=plot_1,"Output_plots/Sample_Read_Depth_Barcode_1.png", device = "png", width = 12, height = 8, units = "in")

### Visualize Read Counts Across Samples for Barcode_2

ASV.table_used %>% 
  group_by(sample) %>%
  filter(., Miseq_run== barcode_2) %>% 
  mutate (TotalReadsperSample = sum(nReads)) %>%
  arrange(desc(TotalReadsperSample)) %>%
  ggplot(., aes(x=sample, y=TotalReadsperSample, color=source)) + 
  geom_point() +ggtitle("Read Count Across Samples") + 
  theme(axis.text.x = element_text(angle = 90)) -> plot_2

ggsave(plot=plot_2,"Output_plots/Sample_Read_Depth_Barcode_2.png", device = "png", width = 12, height = 8, units = "in")

###Nesting the dataset
ASV.table_used %>% 
  dplyr::group_by(Miseq_run, source) %>% 
  nest() %>% 
  pivot_wider(names_from=source, values_from=data) -> ASV.nested


####Summary.file.1
ASV.nested %>% 
  ungroup() %>% 
  dplyr::transmute(.,Miseq_run,Summary = purrr::map(Samples, ~ how.many(ASVtable = ., round = 0)))  -> ASV.summary

ASV.summary$Summary

ASV.nested %>% 
  ungroup() %>% 
  mutate(Step1.cleaned.tibble = purrr::map(Samples, function(.x){ 
    .x %>%
      dplyr::select (sample, seq_number, nReads)
  })) -> ASV.nested

###Add this step to the summary table we were creating

####Summary.file.2
ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step1.cleaned.tibble, ~ how.many(ASVtable = .,round = "1.Jump"))) %>% 
  left_join(ASV.summary) %>% #use left join when there are many miseq runs to join
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary 

ASV.summary$Summary

saveRDS(ASV.nested,file="ASV.nested_initial.RDS")
saveRDS(ASV.summary, file="ASV.summary_initial.RDS")
#---

#---

# Cleaning Process 1.5: List of Manual Sequences to Remove for Decontamination ----------------------------------------------------------------------------------------------------------------
hash.key <- readRDS(file="hash.key.RDS")
ASV.nested <- readRDS(,file="ASV.nested_initial.RDS")
ASV.summary <- readRDS(, file="ASV.summary_initial.RDS")

taxa_to_remove <- c("Eukaryota;Chordata;Mammalia;Primates;Hominidae;Homo;Homo sapiens",
                    "","Unassigned",
                    "Eukaryota;Chordata;Aves;Galliformes;Phasianidae;Meleagris;Meleagris gallopavo",
                    "NA",
                    "Eukaryota;Chordata;Actinopteri;Salmoniformes;Salmonidae;Salmo;Salmo salar",
                    "Eukaryota;Chordata;Mammalia;Artiodactyla;Bovidae;Bos;Bos taurus",
                    "Eukaryota;Chordata;Actinopteri;Salmoniformes;Salmonidae;Salmo;",
                    "Eukaryota;Chordata;Actinopteri;Cypriniformes;Cyprinidae;Ctenopharyngodon;Ctenopharyngodon idella",
                    "Eukaryota;Chordata;Aves;Galliformes;Phasianidae;Gallus;Gallus gallus",
                    "Eukaryota;Chordata;Mammalia;Artiodactyla;Suidae;Sus;Sus scrofa",
                    "Eukaryota;Chordata;Chondrichthyes;Rajiformes;Rajidae;Raja;Raja binoculata",                                    
                    "Eukaryota;Chordata;Mammalia;Artiodactyla;Bovidae;Ovis;Ovis aries",                                    
                    "Eukaryota;Chordata;Mammalia;Carnivora;Canidae;Canis;Canis lupus",                                    
                    "Eukaryota;Chordata;Mammalia;Carnivora;Felidae;Felis;Felis catus",                                  
                    "Eukaryota;Chordata;Mammalia;Primates;Hominidae;Pan;Pan paniscus",
                    "Eukaryota;Chordata;Mammalia;Primates;Aotidae;Aotus;"
)

ASV.nested %>% 
  mutate (Step1.5.cleaned.tibble = purrr::map(Step1.cleaned.tibble, 
                                              function(.x){
                                                .x %>%
                                                  left_join(hash.key) %>% 
                                                  dplyr::select(sample, seq_number, nReads, sum.taxonomy)
                                              }) ) -> ASV.nested


ASV.nested %>% 
  mutate(Step1.5.tibble = purrr::map (Step1.5.cleaned.tibble, ~ filter(.,!sum.taxonomy %in% taxa_to_remove) %>% ungroup))-> ASV.nested

ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step1.5.tibble, ~ how.many(ASVtable = .,round ="2.5.Manual_Remove"))) %>% 
  left_join(ASV.summary) %>%
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary

ASV.summary$Summary

#---

#---


# Cleaning Process 2: **Site Occupancy Modeling** ----------------------------------------------------------------------------------------------------------------

##STAN Model
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#Metadata test run issue: 2 sequencing runs used different names for the same site so metadata had repeated data. 
#Needed to make a metadata file that was distinct(), otherwise each sample gets duplicated in the spread
metadata %>%  dplyr::select(-Seq_number) %>%  distinct() -> metadata

#Determine which SOM level

####Set up data with different hierarchical levels and reduce down to unique patterns.
#Format 2 level hierarcjocal patterns
ASV.nested %>% 
  dplyr::select(Step1.5.tibble) %>% 
  unnest(Step1.5.tibble) %>% # This stage is a long table
  ungroup() %>%
  mutate(nReads = 1) %>% 
  left_join(metadata, by= c("sample"="New_name")) %>% 
  dplyr::select(seq_number, Site, Bio_rep, Tech_rep, nReads) %>% 
  group_by(Site) %>% 
  mutate(., K = n_distinct(Tech_rep)) -> data_som_int

data_som_int %>% 
  pivot_wider(values_from = nReads, names_from = Tech_rep, names_prefix = "Tech") %>% 
  ungroup() %>% 
  mutate(ndetections = rowSums(dplyr::select(., starts_with("Tech")),na.rm = TRUE)) %>% 
  dplyr::select(., !starts_with("Tech"))  %>% 
  pivot_wider(values_from = ndetections, names_from = Bio_rep, names_prefix = "Bio_rep") -> data_tech_summarized


#Now need to try to reduce the number of repeated tests. 
#Many species have identical patterns of presence (e.g. pattern of technical replicates a species was detected)
data_tech_summarized %>% 
  unite(repetition.level,starts_with("Bio_rep"), sep = ".", remove = FALSE) -> Pattern.of.presence_tech

#Select unique patterns
Pattern.of.presence_tech %>% 
  group_by(repetition.level) %>% 
  slice(1) -> unique_patterns

#Format for STAN Model
# K is number of technical replicates for each Biological replicate
# ndetections is the number of times a species was detected out of K trials

unique_patterns %>% 
  pivot_longer(., cols=starts_with("Bio_rep"), names_to="Bio_rep",values_to="ndetections") %>% 
  filter(., !is.na(ndetections)) %>% 
  mutate(Bio_rep =  str_replace_all(Bio_rep, "Bio_rep","")) -> unique_data

#Create unique identifier for species; for use in hierarchical modeling
#STAN model requires species names to be numbers.
SS_species <- unite(data = unique_data,
                    col = SS_species,
                    c("seq_number")
) %>% pull(SS_species)
unique_data$Species <- match(SS_species, unique(SS_species)) #index for unique site-species combinations

#Create unique identifier for combinations of site-biological_rep-species; for use in hierarchical modeling
SDS <- unite(data = unique_data,
             col = SDS,
             c("Site", "Bio_rep", "seq_number")
) %>% pull(SDS)
unique_data$SiteDateSpecies <- match(SDS, unique(SDS)) #index for unique site-date-species combinations

#Create unique identifier for combinations of site-species; for use in hierarchical modeling
SS <- unite(data = unique_data,
            col = SS,
            c("Site", "seq_number")
) %>% pull(SS)
unique_data$SiteSpecies <- match(SS, unique(SS)) #index for unique site-species combinations

#Fix names to match STAN model format
unique_data %>%
  mutate(., N= ndetections) %>%
  dplyr::select(-ndetections)-> unique_data

unique_data$N <- as.numeric(unique_data$N)

##Stan Model
sink("Stan_SOM_hierarchical.stan")
cat(
  "data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of samples (nrow)
    int<lower=1> Species[S];    // index of species, each of which will have a different value for p11 and p10
    int<lower=1> Nspecies;    // number of species, each of which will have a different value for p11 and p10
    int<lower=1> L[S];   // index of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> Nloc;   // number of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> K[S];   // number of replicates per site (ncol)
    int<lower=0> N[S]; // number of detections among these replicates
    int z[S];   // integer flag to help estimate psi parameter
}
parameters{/////////////////////////////////////////////////////////////////////
real<lower=0,upper=1> psi[Nloc];  //commonness parameter
real<lower=0,upper=1> p11[Nspecies]; //true positive detection rate
real<lower=0,upper=1> p10[Nspecies]; //false positive detection rate
}
transformed parameters{/////////////////////////////////////////////////////////////////////
}
model{/////////////////////////////////////////////////////////////////////
real p[S];

for (i in 1:S){
z[i] ~ bernoulli(psi[L[i]]);
p[i] = z[i]*p11[Species[i]] + (1-z[i])*p10[Species[i]];
N[i] ~ binomial(K[i], p[i]);
}; 

//priors
psi ~ beta(2,2); 
p11 ~ beta(2,2); 
p10 ~ beta(1,10);
}
generated quantities{
real<lower=0,upper=1> Occupancy_prob[S];    //after inferring parameters above, now calculate occupancy probability for each observation. Equation from Lahoz-Monfort et al. 2015

for (i in 1:S){
Occupancy_prob[i]  = (psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
/ ((psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
+ (((1-psi[L[i]])*(p10[Species[i]]^N[i]))*((1-p10[Species[i]])^(K[i]-N[i])))
);
}
}

",
fill=TRUE)
sink()

##Run Stan model
#note this will take a while the first time you run a particular model, because it needs to compile from C++

myHierarchicalModel <- stan(file = "Stan_SOM_hierarchical.stan", 
                            data = list(
                              S = nrow(unique_data),
                              Species = unique_data$Species,
                              Nspecies = length(unique(unique_data$Species)),
                              L = unique_data$SiteSpecies,
                              Nloc = length(unique(unique_data$SiteDateSpecies)),
                              K = unique_data$K,
                              N = unique_data$N,
                              z = ifelse(unique_data$N > 0, 1, 0)
                            ), 
                            chains = som_chains,   #number of chains
                            iter = som_iterations   #number of iterations per chain
)



myHierarchicalStanResults <- tidy(tibble(as.data.frame(myHierarchicalModel)))

#See histogram of occupancy probabilities for unique patterns of presence
myHierarchicalStanResults %>% 
  filter(grepl("Occupancy_prob", column)) %>% 
  ggplot() + geom_histogram( aes(mean)) +
  ggtitle(label = "Histogram of Mean Occupancy Probabilities", subtitle = "Unique Patterns of Presence") +
  ylab("Count") + xlab("Mean Occupancy Probability")-> plot_3

ggsave(plot=plot_3,"Output_plots/Histogram_of_Mean_Occupancy_Probabilities_unique_patterns_of_presence.png", device = "png", width = 12, height = 8, units = "in")

#Now Re-attach SOM to unique patterns
myHierarchicalStanResults %>% 
  filter(grepl("Occupancy_prob", column)) %>%
  separate(column, into=c("column","SiteDateSpecies"), sep="([\\[\\]])") %>% 
  pivot_wider(., names_from=column, values_from = c(mean, sd))-> myHierarchicalStanResults_wide
myHierarchicalStanResults_wide$SiteDateSpecies <- as.numeric(myHierarchicalStanResults_wide$SiteDateSpecies)

unique_data %>% 
  left_join(myHierarchicalStanResults_wide) %>% 
  group_by(seq_number,Site) %>% 
  summarise(., max_Occupancy_prob = max(mean_Occupancy_prob)) -> unique_data_SOM

unique_patterns %>% 
  left_join(unique_data_SOM, by=c("seq_number","Site")) %>% 
  dplyr::select(-seq_number,-Site)-> unique_patterns_SOM

Pattern.of.presence_tech %>% 
  left_join(unique_patterns_SOM, by=c("K","repetition.level")) %>% 
  left_join(hash.key)-> SOM_results

#See histogram of occupancy probabilities for unique patterns of presence
SOM_results %>% 
  ggplot() + geom_histogram( aes(max_Occupancy_prob)) +
  ggtitle(label = "Histogram of Mean Occupancy Probabilities", subtitle = "All ASVs") +
  ylab("Count") + xlab("Mean Occupancy Probability")-> plot_4

ggsave(plot=plot_4,"Output_plots/Histogram_of_Mean_Occupancy_Probabilities_all_ASVs.png", device = "png", width = 12, height = 8, units = "in")


unique_patterns_SOM %>% 
  ungroup() %>% 
  mutate (ndetections = rowSums(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) %>% 
  filter(., ndetections ==1) %>% 
  summarise(max(max_Occupancy_prob)) -> som_min_threshold_1

unique_patterns_SOM %>% 
  ungroup() %>% 
  mutate(maxdetections = rowMaxs(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) %>% 
  filter(., maxdetections == 1) %>% 
  summarise(max(max_Occupancy_prob)) -> som_min_threshold_2

unique_patterns_SOM %>% 
  ungroup() %>% 
  mutate(maxdetections = rowMaxs(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) %>% 
  filter(., maxdetections == 3) %>% 
  summarise(min(max_Occupancy_prob)) -> som_min_threshold_3


list(som_min_threshold_1,som_min_threshold_2,som_min_threshold_3) -> som_min_thresh
names(som_min_thresh) <- c("som_min_threshold_1","som_min_threshold_2","som_min_threshold_3")

#Summarize
SOM_results %>% 
  mutate(., to_keep = case_when(max_Occupancy_prob < som_min_threshold_1$`max(max_Occupancy_prob)` ~ "Toss",
                                max_Occupancy_prob > som_min_threshold_2$`max(max_Occupancy_prob)` ~ "Keep",
                                TRUE ~ "Borderline")) %>% 
  count(to_keep) %>% 
  mutate(., total=sum(n)) %>% 
  mutate(per=paste0(round(100*n/total,2),'%')) %>% 
  dplyr::select(-total)

#Save Data
SOM_results %>% 
  dplyr::select(-K, -repetition.level)-> SOM_results_to_keep_all

saveRDS(SOM_results_to_keep_all,"occupancy_results_all.RDS")

SOM_results_to_keep_all %>% 
  write.csv(paste0("Occupancy_all_",Sys.Date(),".csv"), row.names = F)

SOM_results_to_keep_all %>% 
  filter(., max_Occupancy_prob > som_min_thresh[[som_filter_choice]][[1]]) -> SOM_results_to_keep_filtered

SOM_results_to_keep_filtered %>% 
  filter(., sum.taxonomy =="Eukaryota;Chordata;Mammalia;Artiodactyla;Balaenopteridae;Megaptera;Megaptera novaeangliae")

SOM_results_to_keep_all %>% 
  filter(., str_detect(sum.taxonomy, "Phanerodon furcatus"))

saveRDS(SOM_results_to_keep_filtered,"occupancy_results_filtered.RDS")

SOM_results_to_keep_filtered %>% 
  write.csv(paste0("Occupancy_all_",Sys.Date(),".csv"), row.names = F)

ASV.nested %>% 
  mutate(Step2.tibble = purrr::map (Step1.5.tibble, ~ filter(.,sum.taxonomy %in% SOM_results_to_keep_filtered$sum.taxonomy) %>% ungroup))-> ASV.nested

ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step2.tibble, ~ how.many(ASVtable = .,round ="2.Occupancy"))) %>% 
  left_join(ASV.summary) %>%
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary

ASV.summary$Summary

ASV.nested$Step2.tibble[[2]] %>% 
  filter(., sum.taxonomy =="Eukaryota;Chordata;Mammalia;Artiodactyla;Ziphiidae;Ziphius;Ziphius cavirostris")

ASV.nested$Step2.tibble[[1]] %>% 
  filter(., sum.taxonomy =="Eukaryota;Chordata;Mammalia;Artiodactyla;Balaenopteridae;Megaptera;Megaptera novaeangliae")


# Cleaning Process 3: **Discarding PCR replicates with low number of reads** ----------------------------------------------------------------------------------------------------------------

###Pull out Sample Read Depth

ASV.nested %>% 
  dplyr::select(Miseq_run, Step2.tibble) %>% 
  unnest(Step2.tibble) %>% 
  group_by(Miseq_run,sample) %>%
  dplyr::summarise(tot = sum(nReads)) %>% 
  arrange(desc(tot)) %>% 
  unite(col = sample_name, c("Miseq_run","sample"))-> all.reps


all.reps %>%  
  pull(tot) -> reads.per.sample

names(reads.per.sample) <- all.reps %>% pull(sample_name)  


### Fit Normal Distribution

fit <- fitdist(reads.per.sample, "gamma", lower=c(0,0), start=list(scale=1,shape=1))

all.reps %>%  
  mutate(prob = pgamma(tot, shape = fit$estimate[[2]], scale = fit$estimate[[1]], lower.tail = TRUE,
                       log.p = FALSE)) -> all.reps

### Remove Outlier Samples
outliers <- all.reps %>% 
  filter(prob < low_dist_probability_cutoff  | tot < minimum_read_cutoff) # changed to 0.05 to save the two samples

ASV.nested %>% 
  mutate (Step3.tibble.edited = purrr::map(Step2.tibble, 
                                           function(.x){
                                             .x %>%
                                               mutate(., Miseq_run= case_when(str_detect(seq_number,barcode_1) ~ barcode_1,
                                                                              str_detect(seq_number,barcode_2) ~ barcode_2)) %>% 
                                               unite(col = sample_name, c("Miseq_run","sample"), remove = FALSE) %>% 
                                               dplyr::select(-Miseq_run)
                                           }) ) -> ASV.nested


ASV.nested %>% 
  mutate(Step.3.low.read.depth = purrr::map (Step3.tibble.edited,
                                             ~ filter(.,!sample_name %in% outliers$sample_name) %>% ungroup)) -> ASV.nested

ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step.3.low.read.depth, ~ how.many(ASVtable = .,round = "3.Low.Read.Depth"))) %>% 
  left_join(ASV.summary) %>% 
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary 

ASV.summary$Summary



#---

#---

# Cleaning Process 4: **Dissimilarity between PCR (biological) replicates** ----------------------------------------------------------------------------------------------------------------


ASV.nested %>% 
  dplyr::select(Miseq_run,Step.3.low.read.depth) %>% 
  unnest(Step.3.low.read.depth) %>% 
  as.data.frame() %>% 
  ungroup() %>% 
  left_join(metadata, by=c("sample"="New_name"))-> cleaned.tibble.post_occ


ASV.nested %>% 
  dplyr::select(Miseq_run,Step1.cleaned.tibble) %>% 
  unnest(Step1.cleaned.tibble) %>% 
  as.data.frame() %>% 
  ungroup() %>% 
  left_join(metadata, by=c("sample"="New_name"))-> cleaned.tibble.pre_occ


## How many samples, how many ASVs
cleaned.tibble.post_occ %>% 
  dplyr::summarise(n_distinct(sample),
                   n_distinct(seq_number))

cleaned.tibble.pre_occ %>% 
  dplyr::summarise(n_distinct(sample),
                   n_distinct(seq_number))


## eDNA Index
cleaned.tibble.post_occ %>%
  dplyr::group_by (sample) %>%
  mutate (Tot = sum(nReads),
          Row.sums = nReads / Tot) %>% 
  dplyr::group_by (seq_number) %>%
  mutate (Colmax = max (Row.sums),
          Normalized.reads = Row.sums / Colmax) -> cleaned.tibble.post_occ #transforms raw number of reads to eDNA index

cleaned.tibble.pre_occ %>%
  dplyr::group_by (sample) %>%
  mutate (Tot = sum(nReads),
          Row.sums = nReads / Tot) %>% 
  dplyr::group_by (seq_number) %>%
  mutate (Colmax = max (Row.sums),
          Normalized.reads = Row.sums / Colmax) -> cleaned.tibble.pre_occ #transforms raw number of reads to eDNA index

#Calculate Vegan Distances
tibble_to_vegdist_all (cleaned.tibble.post_occ) -> all.distances.full.post
tibble_to_vegdist_all (cleaned.tibble.pre_occ) -> all.distances.full.pre


as.tibble(subset(melt(as.matrix(all.distances.full.post)))) -> all.distances.melted.post
as.tibble(subset(melt(as.matrix(all.distances.full.pre)))) -> all.distances.melted.pre


# Now, create a three variables for all distances, they could be PCR replicates, BIOL replicates, or from the same site
all.distances.melted.pre %>%
  separate(Var1, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
  left_join(metadata) %>% 
  mutate(., Site1=Site,
         Bio_rep1=Bio_rep,
         Tech_rep1=Tech_rep,
         Miseq_run1=Miseq_run) %>% 
  dplyr::select(Var1,Miseq_run1,Site1,Bio_rep1,Tech_rep1,Var2,value) %>% 
  separate(Var2, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
  left_join(metadata) %>% 
  mutate(., Site2=Site,
         Bio_rep2=Bio_rep,
         Tech_rep2=Tech_rep,
         Miseq_run2=Miseq_run) %>% 
  dplyr::select(Var1,Var2,Site1,Miseq_run1,Bio_rep1,Tech_rep1,Miseq_run2,Site2,Bio_rep2,Tech_rep2,value) %>%
  unite( Site1, Bio_rep1, col= "station1", remove=F) %>% 
  unite( Site2, Bio_rep2, col= "station2", remove=F) %>% 
  unite( Miseq_run1, station1, col= "typers1", remove=F) %>% 
  unite( Miseq_run2, station2, col= "typers2", remove=F) %>% 
  unite( Miseq_run1, Site1, col= "type_site1", remove=F) %>% 
  unite( Miseq_run2, Site2, col= "type_site2", remove=F) %>% 
  ungroup() %>% 
  mutate(Distance.type = case_when(typers1 == typers2 ~ "PCR.replicates",
                                   type_site1 == type_site2 ~ "Biol.replicates",
                                   Miseq_run1 == Miseq_run2 ~ "Same Barcode Different Site", 
                                   TRUE ~ "Different Barcode")) %>%
  dplyr::select(Sample1 = Var1, Sample2 = Var2 , value , Distance.type) %>%
  filter (Sample1 != Sample2) -> all.distances.to.plot.pre



# Checking all went well
sapply(all.distances.to.plot.pre, function(x) summary(is.na(x)))

all.distances.to.plot.pre$Distance.type <- all.distances.to.plot.pre$Distance.type  %>% fct_relevel( "PCR.replicates", "Biol.replicates", "Same Barcode Different Site", "Different Barcode")

plot_5 <- ggplot (all.distances.to.plot.pre , aes (fill = Distance.type, x = value,after_stat(density))) +
  geom_histogram(stat = 'bin', alpha = 0.9, binwidth = 0.05) + #xlim(0, 1) +
  facet_wrap( ~ Distance.type) +
  labs (x = "Pairwise Dissimilarity", y = "Density" ,
        fill = "Groups", title = "eDNA Pairwise Dissimilarity Between Samples", subtitle = "Pre Occupancy") +theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave(plot=plot_5,"Output_plots/eDNA_Pairwise_Dissimilarity_Between_Samples_pre_occupancy.png", device = "png", width = 12, height = 8, units = "in")

all.distances.melted.post %>%
  separate(Var1, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
  left_join(metadata) %>% 
  mutate(., Site1=Site,
         Bio_rep1=Bio_rep,
         Tech_rep1=Tech_rep,
         Miseq_run1=Miseq_run) %>% 
  dplyr::select(Var1,Miseq_run1,Site1,Bio_rep1,Tech_rep1,Var2,value) %>% 
  separate(Var2, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
  left_join(metadata) %>% 
  mutate(., Site2=Site,
         Bio_rep2=Bio_rep,
         Tech_rep2=Tech_rep,
         Miseq_run2=Miseq_run) %>% 
  dplyr::select(Var1,Var2,Site1,Miseq_run1,Bio_rep1,Tech_rep1,Miseq_run2,Site2,Bio_rep2,Tech_rep2,value) %>%
  unite( Site1, Bio_rep1, col= "station1", remove=F) %>% 
  unite( Site2, Bio_rep2, col= "station2", remove=F) %>% 
  unite( Miseq_run1, station1, col= "typers1", remove=F) %>% 
  unite( Miseq_run2, station2, col= "typers2", remove=F) %>% 
  unite( Miseq_run1, Site1, col= "type_site1", remove=F) %>% 
  unite( Miseq_run2, Site2, col= "type_site2", remove=F) %>% 
  ungroup() %>% 
  mutate(Distance.type = case_when(typers1 == typers2 ~ "PCR.replicates",
                                   type_site1 == type_site2 ~ "Biol.replicates",
                                   Miseq_run1 == Miseq_run2 ~ "Same Barcode Different Site", 
                                   TRUE ~ "Different Barcode")) %>%
  dplyr::select(Sample1 = Var1, Sample2 = Var2 , value , Distance.type) %>%
  filter (Sample1 != Sample2) -> all.distances.to.plot.post


# Checking all went well
sapply(all.distances.to.plot.post, function(x) summary(is.na(x)))

all.distances.to.plot.post$Distance.type <- all.distances.to.plot.post$Distance.type  %>% fct_relevel( "PCR.replicates", "Biol.replicates", "Same Barcode Different Site", "Different Barcode")

plot_6 <- ggplot (all.distances.to.plot.post , aes (fill = Distance.type, x = value,after_stat(density))) +
  geom_histogram(stat = 'bin', alpha = 0.9, binwidth = 0.05) + #xlim(0, 1) +
  facet_wrap( ~ Distance.type) +
  labs (x = "Pairwise Dissimilarity", y = "Density" ,
        fill = "Groups", title = "eDNA Pairwise Dissimilarity Between Samples", subtitle = "Post Site Occupancy") +theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                      panel.background = element_blank())

ggsave(plot=plot_6,"Output_plots/eDNA_Pairwise_Dissimilarity_Between_Samples_post_occupancy.png", device = "png", width = 12, height = 8, units = "in")


# Instead of chosing based on the pairwise distances, we will use distance to centroid

#Need to Remove Single Sampels
cleaned.tibble.post_occ  %>% 
  group_by(Miseq_run,Site) %>% 
  dplyr::summarise(cases = n_distinct(Bio_rep)) %>% 
  filter(., cases == 1) %>% 
  unite(col="location_id", Miseq_run,Site, remove = FALSE, sep=":") %>% ungroup() %>% 
  dplyr::select(location_id) -> singles
#Currently this method does not filter out these single tech reps. There are only 3

#Nest Data
cleaned.tibble.post_occ  %>% 
  unite(col="location_id", Miseq_run,Site, remove = FALSE, sep=":") %>% 
  filter(., !location_id %in% singles$location_id) %>% 
  group_by(location_id) %>% nest() -> nested.cleaning

nested.cleaning %>% 
  mutate(matrix = purrr::map(data, tibble_to_vegdist_bio_rep_single)) -> nested.cleaning

# Convert
nested.cleaning %>% mutate(ncomparisons = purrr::map(matrix, length)) -> nested.cleaning

#Calculate distances to centroid
nested.cleaning <- nested.cleaning %>% mutate(distances = map2(matrix, location_id, dist_to_centroid))
nested.cleaning %>% 
  separate(location_id, c("Miseq_run","Site"), sep=":") %>% 
  unnest_longer(distances) %>% 
  dplyr::select(Miseq_run,Site, distances_id,distances) -> all_distances.groups #unnest data



#Calculate normal distribution of distances to centroid
normparams.step5 <- MASS::fitdistr(all_distances.groups$distances, "normal")$estimate

#Calculate Probability
probs.step5 <- pnorm(all_distances.groups$distances, normparams.step5[1], normparams.step5[2])

#Determine Outliers
outliers.step5 <- which(probs.step5>max_dist_probability_cutoff)

#Remove Outliers
discard.step5 <-names(all_distances.groups$distances[outliers.step5])

to_write_discarded.step5 <- tibble(distances_id = discard.step5,
                                   distance = all_distances.groups$distances[outliers.step5])

to_write_discarded.step5 <- to_write_discarded.step5 %>% left_join(tibble(distances_id = discard.step5,
                                                                          probs = probs.step5[outliers.step5]))
write_csv(to_write_discarded.step5 ,"step5.discared_samples.csv")

## Plot Final Replication Levels

cleaned.tibble.post_occ %>% 
  unite(col="distances_id", Miseq_run,Site,Bio_rep,Tech_rep, remove = FALSE, sep=":") ->cleaned.tibble.post_occ_for_plotting

cleaned.tibble.post_occ_for_plotting %>%
  filter(., !distances_id %in% discard.step5) %>% 
  dplyr::group_by(Miseq_run,Site,Bio_rep) %>% 
  dplyr::summarise(cases = n_distinct(distances_id)) %>% 
  ggplot()+
  geom_raster(aes(x= Site, y = Bio_rep, fill = cases))+
  geom_text(aes(x= Site, y = Bio_rep, label = cases),color="white") +
  ggtitle("Sample Replication Level") +facet_wrap(~Miseq_run) +
  theme(axis.text.x = element_text(angle = 90))-> plot_7

ggsave(plot=plot_7,"Output_plots/Final_sample_replication_level.png", device = "png", width = 20, height = 8, units = "in")


##Remove Samples

ASV.nested %>% 
  mutate (Step4.tibble_edited = purrr::map(Step.3.low.read.depth, 
                                           function(.x){
                                             .x %>%
                                               mutate(., Miseq_run= str_detect(seq_number,barcode_1)) %>% 
                                               mutate(distances_id =TRUE )
                                           }) ) -> ASV.nested


#Filter Sample
ASV.nested %>% 
  mutate(Step4.tibble = purrr::map (Step4.tibble_edited,  ~ filter(.,! distances_id %in% to_write_discarded.step5$distances_id))) -> ASV.nested


#Visualize Results of Clearence Process 4
ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step4.tibble, ~ how.many(ASVtable = .,round ="4.Dissimilarity"))) %>% 
  left_join(ASV.summary) %>%
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary


ASV.summary$Summary

ASV.nested$Step1.5.tibble[[2]] %>% 
  filter(., sum.taxonomy =="Eukaryota;Chordata;Mammalia;Artiodactyla;Ziphiidae;Ziphius;Ziphius cavirostris")


##Save Data
saveRDS(ASV.nested,file="ASV.nested_final.RDS")
saveRDS(ASV.summary,file="ASV.summary_final.RDS")



#---

#---

#Generate Final Outputs ----------------------------------------------------------------------------------------------------------------

ASV.nested <- readRDS(file="ASV.nested_final.RDS")
ASV.summary <- readRDS(file="ASV.summary_final.RDS")
metadata <- read.table(input_meta_path, header = TRUE, sep = "\t", stringsAsFactors = F)
hash.key <- readRDS(file="hash.key.RDS")

##Code for merging ASV tables

###Identify Unique Species for Sum.taxonomy merging

hash.key %>% 
  distinct(.,sum.taxonomy) -> hashes_unique

hashes_unique$number <- row.names(hashes_unique)
hashes_unique$number <- paste0("taxon_",hashes_unique$number)
row.names(hashes_unique)<-hashes_unique$number

hash.key %>% 
  left_join(hashes_unique, by="sum.taxonomy") -> hash.key.updated


#---

#---

###Pre Occupancy Sum by Taxonomy, All PCR Tech Reps Separate Samples

hash.key.updated$number %>% unique() -> total_taxa

ASV.nested$Step1.5.tibble[[1]] %>% 
  mutate(miseq = ASV.nested$Miseq_run[[1]]) %>% 
  unite(miseq,sample, col="Sample") %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) -> barcode_1_pre_df
barcode_1_pre_df$number %>%  unique() -> barcode_1_pre_df_taxa

ASV.nested$Step1.5.tibble[[2]] %>% 
  mutate(miseq = ASV.nested$Miseq_run[[2]]) %>% 
  unite(miseq,sample, col="Sample") %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) -> barcode_2_pre_df
barcode_2_pre_df$number %>%  unique() -> barcode_2_pre_df_taxa

total_kept_taxa <- (append(barcode_1_pre_df_taxa,barcode_2_pre_df_taxa)) %>% unique()

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_1_pre_df_taxa))) {
  barcode_1_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_1_pre_df
} else {
  barcode_1_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_1_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_1_pre_df
}

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_2_pre_df_taxa))) {
  barcode_2_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_2_pre_df
} else {
  barcode_2_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_2_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_2_pre_df
}

barcode_1_pre_df <- as.data.frame(barcode_1_pre_df)
row.names(barcode_1_pre_df) <- barcode_1_pre_df$number
barcode_1_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_1_pre_df

barcode_2_pre_df <- as.data.frame(barcode_2_pre_df)
row.names(barcode_2_pre_df) <- barcode_2_pre_df$number
barcode_2_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_2_pre_df

####First, we want to create proportions by dividing by the rowsums:
####we could do this with sweep() or mutate_all() or other ways, but using vegan:

barcode_1_prop <- decostand(barcode_1_pre_df, method = "total", MARGIN = 2)
barcode_2_prop <- decostand(barcode_2_pre_df, method = "total", MARGIN = 2)

####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
####We do this by scaling everything to the max observed in each row. 

####to do this WITHIN a dataset, we could just do (again, using vegan):
barcode_1_prop_index <- decostand(barcode_1_prop, method = "max", MARGIN = 1)
barcode_2_prop_index <- decostand(barcode_2_prop, method = "max", MARGIN = 1)

####This gives us an index between 0 and 1 for each species in each dataset.  

####But if we want to combine datasets, this second step has to happen in the combined dataset, so it all gets scaled to 0-1.  
####easy enough:

combined_index <- decostand(cbind(barcode_1_prop,barcode_2_prop), method = "max", MARGIN = 1)
####How both datasets are combined, on a common, comparable scale.

### Output Read Count Data
pre_results_reads = cbind(barcode_1_prop,barcode_2_prop)

hash.key.updated.2 <- hash.key.updated[!duplicated(hash.key.updated$number), ]

pre_results_reads$number <- rownames(pre_results_reads)

pre_results_reads %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> pre_results_reads

saveRDS(pre_results_reads,file="Output_R/pre_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.RDS")
write_csv(pre_results_reads ,"Output_csv/pre_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.csv")

### Output eDNA Index Data

combined_index$number <- rownames(combined_index)

combined_index %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_index

saveRDS(combined_index,file="Output_R/pre_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.RDS")
write_csv(combined_index ,"Output_csv/pre_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.csv")

#---

#---

###Post Occupancy Sum by Taxonomy, All PCR Tech Reps Separate Samples

hash.key.updated$number %>% unique() -> total_taxa

ASV.nested$Step4.tibble[[1]] %>% 
  mutate(miseq = ASV.nested$Miseq_run[[1]]) %>% 
  unite(miseq,sample, col="Sample") %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) -> barcode_1_post_occ_df
barcode_1_post_occ_df$number %>%  unique() -> barcode_1_post_occ_df_taxa


ASV.nested$Step4.tibble[[2]] %>% 
  mutate(miseq = ASV.nested$Miseq_run[[2]]) %>% 
  unite(miseq,sample, col="Sample") %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0) -> barcode_2_post_occ_df
barcode_2_post_occ_df$number %>%  unique() -> barcode_2_post_occ_df_taxa

total_kept_taxa <- append(barcode_1_post_occ_df_taxa,barcode_2_post_occ_df_taxa) %>% unique()

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_1_post_occ_df_taxa))) {
  barcode_1_post_occ_df %>% ungroup() %>% 
    arrange(number) -> barcode_1_post_occ_df
} else {
  barcode_1_post_occ_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_1_post_occ_df_taxa)) %>%
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_1_post_occ_df
}

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_2_post_occ_df_taxa))) {
  barcode_2_post_occ_df %>% ungroup() %>% 
    arrange(number) -> barcode_2_post_occ_df
} else {
  barcode_2_post_occ_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_2_post_occ_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_2_post_occ_df
}


barcode_1_post_occ_df <- as.data.frame(barcode_1_post_occ_df)
row.names(barcode_1_post_occ_df) <- barcode_1_post_occ_df$number
barcode_1_post_occ_df %>% ungroup() %>% dplyr::select(-number) -> barcode_1_post_occ_df

barcode_2_post_occ_df <- as.data.frame(barcode_2_post_occ_df)
row.names(barcode_2_post_occ_df) <- barcode_2_post_occ_df$number
barcode_2_post_occ_df %>% ungroup() %>% dplyr::select(-number) -> barcode_2_post_occ_df

####first, we want to create proportions by dividing by the rowsums:
####we could do this with sweep() or mutate_all() or other ways, but using vegan:

barcode_1_post_occ_df_prop <- decostand(barcode_1_post_occ_df, method = "total", MARGIN = 2)
barcode_2_post_occ_df_prop <- decostand(barcode_2_post_occ_df, method = "total", MARGIN = 2)

####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
####We do this by scaling everything to the max observed in each row. 

####to do this WITHIN a dataset, we could just do (again, using vegan):
barcode_1_post_occ_df_index <- decostand(barcode_1_post_occ_df_prop, method = "max", MARGIN = 1)
barcode_2_post_occ_df_index <- decostand(barcode_2_post_occ_df_prop, method = "max", MARGIN = 1)

####this gives us an index between 0 and 1 for each species in each dataset.  

####But if we want to combine datasets, this second step has to happen in the combined dataset, so it all gets scaled to 0-1.  
####easy enough:

combined_index <- decostand(cbind(barcode_1_post_occ_df_index,barcode_2_post_occ_df_index), method = "max", MARGIN = 1)
####now both datasets are combined, on a common, comparable scale.

### Output Post Occupancy Read Data
results_post_occ_reads = cbind(barcode_1_post_occ_df,barcode_2_post_occ_df)

hash.key.updated.2 <- hash.key.updated[!duplicated(hash.key.updated$number), ]

results_post_occ_reads$number <- rownames(results_post_occ_reads)

results_post_occ_reads %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> results_post_occ_reads

saveRDS(results_post_occ_reads,file="Output_R/post_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.RDS")
write_csv(results_post_occ_reads ,"Output_csv/post_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.csv")

### Output Post Occupancy eDNA Index Data
combined_index$number <- rownames(combined_index)

combined_index %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_index

saveRDS(combined_index,file="Output_R/post_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.RDS")
write_csv(combined_index ,"Output_csv/post_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.csv")

#Code for Merging Tech Reps

###Pre Occupancy Sum by Taxonomy, Biological Replicates Separate

ASV.nested$Step1.5.tibble[[1]] %>% 
  mutate(Sample = sample) %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> barcode_1_pre_df
barcode_1_pre_df$number %>%  unique() -> barcode_1_pre_df_taxa

metadata %>% 
  filter(., !(New_name %in% controls)) %>% 
  filter(.,!(New_name %in% colnames(barcode_1_pre_df))) %>% 
  pull(New_name) %>%  unique()-> columns2add

barcode_1_pre_df <- as.data.frame(barcode_1_pre_df)

barcode_1_pre_df %>% 
  tibble::add_column(!!!set_names(as.list(rep(NA, length(columns2add))),nm=columns2add)) %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(sort(tidyselect::peek_vars()))-> barcode_1_pre_df

row.names(barcode_1_pre_df) <- barcode_1_pre_df$number
barcode_1_pre_df %>% ungroup() -> barcode_1_pre_df

ASV.nested$Step1.5.tibble[[2]] %>% 
  mutate(Sample = sample) %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> barcode_2_pre_df
barcode_2_pre_df$number %>%  unique() -> barcode_2_pre_df_taxa

metadata %>% 
  filter(., !(New_name %in% controls)) %>% 
  filter(.,!(New_name %in% colnames(barcode_2_pre_df))) %>% 
  pull(New_name) %>%  unique()-> columns2add

barcode_2_pre_df <- as.data.frame(barcode_2_pre_df)

barcode_2_pre_df %>% 
  tibble::add_column(!!!set_names(as.list(rep(NA, length(columns2add))),nm=columns2add)) %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(sort(tidyselect::peek_vars())) -> barcode_2_pre_df

row.names(barcode_2_pre_df) <- barcode_2_pre_df$number
barcode_2_pre_df %>% ungroup() -> barcode_2_pre_df

total_kept_taxa <- (append(barcode_1_pre_df_taxa,barcode_2_pre_df_taxa)) %>% unique()

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_1_pre_df_taxa))) {
  barcode_1_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_1_pre_df
} else {
  barcode_1_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_1_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_1_pre_df
}

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_2_pre_df_taxa))) {
  barcode_2_pre_df %>% ungroup() %>% 
    arrange(number) -> barcode_2_pre_df
} else {
  barcode_2_pre_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_2_pre_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_2_pre_df
}

barcode_1_pre_df <- as.data.frame(barcode_1_pre_df)
row.names(barcode_1_pre_df) <- barcode_1_pre_df$number
barcode_1_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_1_pre_df

barcode_2_pre_df <- as.data.frame(barcode_2_pre_df)
row.names(barcode_2_pre_df) <- barcode_2_pre_df$number
barcode_2_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_2_pre_df

barcode_1_pre_df+barcode_2_pre_df -> combined_df

####First, we want to create proportions by dividing by the rowsums:
####we could do this with sweep() or mutate_all() or other ways, but using vegan:

combined_df_prop <- decostand(combined_df, method = "total", MARGIN = 2)

####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
####We do this by scaling everything to the max observed in each row. 

####To do this WITHIN a dataset, we could just do (again, using vegan):
####eDNA Index between 0-1 made by straight adding reads across each decontaminated datatable
combined_df_index <- decostand(combined_df_prop, method = "max", MARGIN = 1)

#Output read counts
hash.key.updated.2 <- hash.key.updated[!duplicated(hash.key.updated$number), ]

combined_df$number <- rownames(combined_df)

combined_df %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_df

saveRDS(combined_df,file="Output_R/pre_occupancy_results_sum.taxonomy_tech_reps_summed_read_counts.RDS")
write_csv(combined_df,"Output_csv/pre_occupancy_results_sum.taxonomy_tech_reps_summed_read_counts.csv")

### Output eDNA Index
combined_df_index$number <- rownames(combined_df_index)

combined_df_index %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_df_index

saveRDS(combined_df_index,file="Output_R/pre_occupancy_results_sum.taxonomy_tech_reps_summed_eDNA_index.RDS")
write_csv(combined_df_index,"Output_csv/pre_occupancy_results_sum.taxonomy_tech_reps_summed_eDNA_index.csv")


###Post Occupancy Merge, Bio reps sep

ASV.nested$Step4.tibble[[1]] %>% 
  mutate(Sample = sample) %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> barcode_1_post_df
barcode_1_post_df$number %>%  unique() -> barcode_1_post_df_taxa

metadata %>% 
  filter(., !(New_name %in% controls)) %>% 
  filter(.,!(New_name %in% colnames(barcode_1_post_df))) %>% 
  pull(New_name) %>%  unique()-> columns2add

barcode_1_post_df <- as.data.frame(barcode_1_post_df)

barcode_1_post_df %>% 
  tibble::add_column(!!!set_names(as.list(rep(NA, length(columns2add))),nm=columns2add)) %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(sort(tidyselect::peek_vars()))-> barcode_1_post_df

row.names(barcode_1_post_df) <- barcode_1_post_df$number
barcode_1_post_df %>% ungroup() -> barcode_1_post_df

ASV.nested$Step4.tibble[[2]] %>% 
  mutate(Sample = sample) %>% 
  left_join(hash.key.updated) %>% 
  dplyr::group_by(number,Sample) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Sample, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> barcode_2_post_df
barcode_2_post_df$number %>%  unique() -> barcode_2_post_df_taxa

metadata %>% 
  filter(., !(New_name %in% controls)) %>% 
  filter(.,!(New_name %in% colnames(barcode_2_post_df))) %>% 
  pull(New_name) %>%  unique()-> columns2add

barcode_2_post_df <- as.data.frame(barcode_2_post_df)

barcode_2_post_df %>% 
  tibble::add_column(!!!set_names(as.list(rep(NA, length(columns2add))),nm=columns2add)) %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(sort(tidyselect::peek_vars())) -> barcode_2_post_df

row.names(barcode_2_post_df) <- barcode_2_post_df$number
barcode_2_post_df %>% ungroup() -> barcode_2_post_df

total_kept_taxa <- (append(barcode_1_post_df_taxa,barcode_2_post_df_taxa)) %>% unique()

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_1_post_df_taxa))) {
  barcode_1_post_df %>% ungroup() %>% 
    arrange(number) -> barcode_1_post_df
} else {
  barcode_1_post_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_1_post_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_1_post_df
}

if (rlang::is_empty(setdiff(total_kept_taxa,barcode_2_post_df_taxa))) {
  barcode_2_post_df %>% ungroup() %>% 
    arrange(number) -> barcode_2_post_df
} else {
  barcode_2_post_df %>% ungroup() %>% 
    add_row(number=setdiff(total_kept_taxa,barcode_2_post_df_taxa)) %>% 
    arrange(number) %>% 
    replace(is.na(.), 0) -> barcode_2_post_df
}

barcode_1_post_df <- as.data.frame(barcode_1_post_df)
row.names(barcode_1_post_df) <- barcode_1_post_df$number
barcode_1_post_df %>% ungroup() %>% dplyr::select(-number) -> barcode_1_post_df

barcode_2_post_df <- as.data.frame(barcode_2_post_df)
row.names(barcode_2_post_df) <- barcode_2_post_df$number
barcode_2_post_df %>% ungroup() %>% dplyr::select(-number) -> barcode_2_post_df

barcode_1_post_df+barcode_2_post_df -> combined_df

####First, we want to create proportions by dividing by the rowsums:
####we could do this with sweep() or mutate_all() or other ways, but using vegan:

combined_df_prop <- decostand(combined_df, method = "total", MARGIN = 2)

####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
####We do this by scaling everything to the max observed in each row. 

####To do this WITHIN a dataset, we could just do (again, using vegan):
####eDNA Index between 0-1 made by straight adding reads across each decontaminated datatable
combined_df_index <- decostand(combined_df_prop, method = "max", MARGIN = 1)

#Output read counts
hash.key.updated.2 <- hash.key.updated[!duplicated(hash.key.updated$number), ]

combined_df$number <- rownames(combined_df)

combined_df %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_df

saveRDS(combined_df,file="Output_R/post_occupancy_results_sum.taxonomy_tech_reps_summed_read_counts.RDS")
write_csv(combined_df,"Output_csv/post_occupancy_results_sum.taxonomy_tech_reps_summed_read_counts.csv")

### Output eDNA Index
combined_df_index$number <- rownames(combined_df_index)

combined_df_index %>% 
  left_join(hash.key.updated.2, by="number") %>% 
  dplyr::select(-number,-seq_number) -> combined_df_index

saveRDS(combined_df_index,file="Output_R/post_occupancy_results_sum.taxonomy_tech_reps_summed_eDNA_index.RDS")
write_csv(combined_df_index,"Output_csv/post_occupancy_results_sum.taxonomy_tech_reps_summed_eDNA_index.csv")

#---

#---

#Code for Merging Sites

###Pre Occupancy Merge by Site
ASV.nested %>% 
  dplyr::select(Step1.5.tibble,Miseq_run) %>% 
  unnest(Step1.5.tibble) %>% 
  left_join(hash.key.updated) %>% 
  left_join(metadata, by=c("sample"="New_name")) %>% 
  dplyr::group_by(sum.taxonomy,Site) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Site, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> pre_wide_reads

saveRDS(pre_wide_reads,file="Output_R/pre_occ_site_sum.taxonomy_reads_summed.RDS")
write_csv(pre_wide_reads ,"Output_csv/pre_occ_site_sum.taxonomy_reads_summed.csv")

ASV.nested %>% 
  dplyr::select(Step1.5.tibble,Miseq_run) %>% 
  unnest(Step1.5.tibble) %>% 
  left_join(hash.key.updated) %>% 
  left_join(metadata, by=c("sample"="New_name")) %>% 
  ungroup() %>% 
  group_by(Site, sum.taxonomy) %>% 
  summarise(meanreads = mean(nReads)) %>% 
  dplyr::group_by(Site) %>%
  mutate (Tot = sum(meanreads),
          Row.sums = meanreads / Tot) %>% 
  dplyr::group_by (sum.taxonomy) %>%
  mutate (Colmax = max (Row.sums),
          Normalized.reads = Row.sums / Colmax) %>% 
  ungroup() %>% 
  dplyr::select(Site,Normalized.reads, sum.taxonomy) %>% 
  spread(., Site, Normalized.reads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> pre_wide

saveRDS(pre_wide,file="Output_R/pre_occ_site_averaged_sum.taxonomy_e_index.RDS")
write_csv(pre_wide,"Output_csv/pre_occ_site_averaged_sum.taxonomy_e_index.csv")

###Post Occupancy Merge by Site
ASV.nested %>% 
  dplyr::select(Step4.tibble) %>% 
  unnest(Step4.tibble) %>% 
  ungroup() %>% 
  left_join(hash.key.updated) %>% 
  left_join(metadata, by=c("sample"="New_name")) %>% 
  dplyr::group_by(sum.taxonomy,Site) %>%
  dplyr::summarise(nReads=sum(nReads)) %>% 
  spread(., Site, nReads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> post_wide_reads

saveRDS(post_wide_reads,file="Output_R/post_occ_site_sum.taxonomy_reads_summed.RDS")
write_csv(post_wide_reads ,"Output_csv/post_occ_site_sum.taxonomy_reads_summed.csv")

ASV.nested %>% 
  dplyr::select(Step4.tibble) %>% 
  unnest(Step4.tibble) %>% 
  left_join(hash.key.updated) %>% 
  left_join(metadata, by=c("sample"="New_name")) %>% 
  ungroup() %>% 
  group_by(Site, sum.taxonomy) %>% 
  summarise(meanreads = mean(nReads)) %>% 
  dplyr::group_by(Site) %>%
  mutate (Tot = sum(meanreads),
          Row.sums = meanreads / Tot) %>% 
  dplyr::group_by (sum.taxonomy) %>%
  mutate (Colmax = max(Row.sums),
          Normalized.reads = Row.sums / Colmax) %>% 
  dplyr::select(Site,Normalized.reads, sum.taxonomy) %>% 
  spread(., Site, Normalized.reads) %>% #convert to wide data format
  replace(is.na(.), 0)  -> post_wide


saveRDS(post_wide,file="Output_R/post_occ_site_averaged_sum.taxonomy_e_index.RDS")
write_csv(post_wide,"Output_csv/post_occ_site_averaged_sum.taxonomy_e_index.csv")

#Data Summary ----------------------------------------------------------------------------------------------------------------

ASV.nested %>% 
  dplyr::select(Step4.tibble) %>% 
  unnest(Step4.tibble) %>%
  filter(., !sample %in% c("LB1_0_1","LB1_0_2","LB1_0_3")) %>% 
  mutate(., barcode = case_when(str_detect(sample_name,"elas")~"Elas",
                                str_detect(sample_name,"miu")~"Miu", 
                                TRUE~"wrong")) %>% 
  group_by(sample, barcode) %>% 
  dplyr::summarise(., tot_reads = sum(nReads)) %>% 
  ungroup() %>% 
  group_by(barcode) %>% 
  dplyr::summarise(., sum(tot_reads), max(tot_reads), min(tot_reads), mean(tot_reads))

ASV.nested %>% 
  dplyr::select(Step4.tibble) %>% 
  unnest(Step4.tibble) %>%
  filter(., !sample %in% c("LB1_0_1","LB1_0_2","LB1_0_3")) %>% 
  mutate(., barcode = case_when(str_detect(sample_name,"elas")~"Elas",
                                str_detect(sample_name,"miu")~"Miu", 
                                TRUE~"wrong")) %>% 
  group_by(sample, barcode) %>% 
  dplyr::summarise(., tot_reads = sum(nReads)) %>% 
  ungroup() %>% 
  dplyr::summarise(., sum(tot_reads), max(tot_reads), min(tot_reads), mean(tot_reads))

ASV.nested %>% 
  dplyr::select(Step4.tibble) %>% 
  unnest(Step4.tibble) %>%
  filter(., !sample %in% c("LB1_0_1","LB1_0_2","LB1_0_3")) %>% 
  mutate(., barcode = case_when(str_detect(sample_name,"elas")~"Elas",
                                str_detect(sample_name,"miu")~"Miu", 
                                TRUE~"wrong")) %>% 
  group_by(barcode) %>% 
  dplyr::summarise(., n_distinct(seq_number), n_distinct(sample))


ASV.nested %>% 
  dplyr::select(Step4.tibble) %>% 
  unnest(Step4.tibble) %>%
  filter(., !sample %in% c("LB1_0_1","LB1_0_2","LB1_0_3")) %>% 
  mutate(., barcode = case_when(str_detect(sample_name,"elas")~"Elas",
                                str_detect(sample_name,"miu")~"Miu", 
                                TRUE~"wrong")) %>% 
  separate(., sum.taxonomy, into= c("D","P","C","O","F","G","S"), sep=";", remove=F) %>% 
  filter(., barcode !="Elas") -> holder


holder %>% 
  filter(., S !="") %>% 
  dplyr::summarise(., n_distinct(seq_number))

holder %>% 
  filter(., S =="") %>%
  filter(., G !="") %>% 
  dplyr::summarise(., n_distinct(seq_number))

holder %>% 
  filter(., S =="") %>%
  filter(., G =="") %>% 
  filter(., `F` !="") %>%
  dplyr::summarise(., n_distinct(seq_number))


holder %>% 
  filter(., S!="") %>% 
  summarise(n_distinct(S))

holder %>% 
  filter(., S=="") %>% 
  filter(., G!="") %>% 
  summarise(n_distinct(G))

holder %>% 
  filter(., S=="") %>% 
  filter(., G=="") %>% 
  filter(., `F`!="") %>% 
  summarise(n_distinct(`F`))


holder %>% 
  filter(., S =="") %>%
  filter(., G =="") %>% 
  filter(., `F` =="") %>%
  filter(., `O` =="") %>%
  filter(., `P` =="") %>%
  filter(., `D` =="") %>%
  dplyr::summarise(., n_distinct(seq_number))


ASV.nested %>% 
  dplyr::select(Step1.cleaned.tibble) %>% 
  unnest(Step1.cleaned.tibble) %>%
  mutate(., barcode = case_when(str_detect(sample,"elas")~"Elas",
                                str_detect(sample,"miu")~"Miu", 
                                TRUE~"wrong")) %>% 
  left_join(hash.key.updated) %>% 
  separate(., sum.taxonomy, into= c("D","P","C","O","F","G","S"), sep=";", remove=F) -> holder2

holder2 %>% 
  filter(., S!="") %>% 
  summarise(n_distinct(S))

holder %>% 
  summarise(n_distinct(S))



ASV.nested_old_inv <- readRDS(file ="/Users/zackgold/Documents/UCLA_phd/Projects/California/Port_LA/Data/r_analysis/ASV.cleaned.rds")

ASV.nested_old_inv$Step.1.low.reads %>%  View()

ASV.nested %>% 
  dplyr::select(Step2.tibble) %>% 
  unnest(Step2.tibble) %>%  
  