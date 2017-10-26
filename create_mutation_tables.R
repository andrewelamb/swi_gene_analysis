# This script was used on an EC2 instance to create gene and variant tables on
# genie data. The results were put in 

library(synapseClient)
library(tidyverse)
library(data.table)
library(magrittr)
library(stringr)

synapseLogin('andrew.lamb@sagebase.org', 'sageBlam1979')
synapseCacheDir("./tmp/")

# GENIE data ------------------------------------------------------------------

genie_bed_df <- synGet("syn7444851")@filePath %>% 
    fread(select = c("Hugo_Symbol", "SEQ_ASSAY_ID")) %>% 
    as_data_frame 

genie_mut_df <- synGet("syn5571527")@filePath %>% 
    fread(select = c("Hugo_Symbol", "Tumor_Sample_Barcode", "HGVSp_Short")) %>% 
    as_data_frame 

genie_pat_df <- synGet("syn9734573")@filePath %>% 
    fread(skip = 4, select = c("SAMPLE_ID", "SEQ_ASSAY_ID")) %>% 
    as_data_frame

# filter mutations for only hugo genes listed in a panel
all_genes    <- sort(unique(genie_bed_df$Hugo_Symbol))
genie_mut_df <- genie_mut_df %>% 
    filter(Hugo_Symbol %in% all_genes) %>% 
    left_join(genie_pat_df, by = c("Tumor_Sample_Barcode" = "SAMPLE_ID")) %>% 
    filter(!is.na(SEQ_ASSAY_ID))

gene_df <- genie_mut_df %>% 
    mutate(mutation = Hugo_Symbol) %>% 
    select(Tumor_Sample_Barcode, Hugo_Symbol, mutation, SEQ_ASSAY_ID)

variant_df <- genie_mut_df %>% 
    .[complete.cases(.),] %>% 
    mutate(HGVSp_Short = str_sub(HGVSp_Short, start = 3)) %>% 
    mutate(mutation = str_c(Hugo_Symbol, ":", HGVSp_Short)) %>% 
    select(Tumor_Sample_Barcode, Hugo_Symbol, mutation, SEQ_ASSAY_ID) 

panel_df <- distinct(genie_bed_df)

all_variants <- variant_df$mutation

create_feature_matrix <- function(mutation_df, panel_df, all_mutations){
    mutation_df %>% 
        split(.$SEQ_ASSAY_ID) %>%
        .[1:2] %>% 
        map(create_matrix_by_panel, panel_df, all_mutations) %>% 
        reduce(rbind)
}

create_matrix_by_panel <- function(mutation_df, panel_df, all_mutations){
    panel    <- mutation_df$SEQ_ASSAY_ID[[1]]
    patients <- sort(unique(mutation_df$Tumor_Sample_Barcode))
    
    all_genes  <- sort(unique(panel_df$Hugo_Symbol))
    
    panel_genes <- panel_df %>% 
        filter(SEQ_ASSAY_ID == panel) %>% 
        use_series(Hugo_Symbol) %>% 
        unique %>% 
        sort
    
    non_panel_genes <- all_genes[!all_genes %in% panel_genes] 
    panel_mutations <- mutation_df %>% 
        filter(Hugo_Symbol %in% panel_genes) %>% 
        use_series(mutation) %>% 
        unique %>% 
        sort
    
    non_panel_mutations <- all_mutations[!all_mutations %in% panel_mutations]
    
    
    matrix1 <- matrix(
        0, 
        nrow = length(patients), 
        ncol = length(panel_mutations), 
        dimnames = list(patients, panel_mutations))
    
    matrix2 <- matrix(
        NA, 
        nrow = length(patients), 
        ncol = length(non_panel_mutations), 
        dimnames = list(patients, non_panel_mutations))
    
    
    for(patient in patients){
        patient_mutations <- filter(mutation_df, Tumor_Sample_Barcode == patient)$mutations
        matrix1[patient, panel_mutations %in% patient_mutations] <- 1
    }
    
    matrix <- cbind(matrix1, matrix2)
    matrix <- matrix[,order(colnames(matrix))]
    return(matrix)
}

remove(genie_bed_df, genie_mut_df, genie_pat_df)

write.table(
    create_feature_matrix(gene_df, panel_df, all_genes),
    file  = "genie_genes.tsv",
    quote = F,
    sep   = "\t")

remove(gene_df, all_genes)

write.table(
    create_matrix_by_panel(variant_df, panel_df, all_variants),
    file  = "genie_variants.tsv",
    quote = F,
    sep   = "\t")

