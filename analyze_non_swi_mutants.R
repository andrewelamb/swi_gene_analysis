# This script creates file sused for plotting swi scores and non
# swi mutants

library(plyr)
library(doMC)
library(data.table)
library(tidyverse)
library(stringr)
library(magrittr)
library(synapseClient)


working_dir      <- "/home/aelamb/repos/swi_gene_analysis/"
boostrap_file    <- "bootstrap_results/bootstrap_scores.tsv"

setwd(working_dir)
synapseLogin()
synapseCacheDir("./tmp/")
registerDoMC(cores = detectCores())



get_barcode_sample <- function(string){
    string %>% 
        str_match("^TCGA-..-....-(..)[:print:]+$") %>% 
        extract2(2)
}


mut_df  <- synGet("syn4924181")@filePath %>% 
    str_c('zcat ', .) %>% 
    fread(
        select = c("Tumor_Sample_Barcode", "Variant_Classification", "Hugo_Symbol"),
        skip = 1) %>% 
    as_data_frame %>% 
    mutate(sample = str_sub(Tumor_Sample_Barcode, 1, 16)) %>% 
    select(-Tumor_Sample_Barcode) 

full_df <- fread(boostrap_file) %>% 
    as_data_frame %>% 
    filter(Hugo_Symbol == "none") %>% 
    select(sample, tumor_type, med_score) %>% 
    left_join(mut_df)

tumors <- sort(unique(full_df$tumor_type))


calculate_wilcox_by_gene <- function(gene, df){
    x <- filter(df, Hugo_Symbol == gene)$med_score
    y <- filter(df, !Hugo_Symbol == gene)$med_score
    pvalue <- wilcox.test(x, y, alternative = "greater")$p.value
    return(pvalue)
}

create_gene_pvalue_df <- function(df){
    df %>% 
        filter(!is.na(Hugo_Symbol)) %>%
        group_by(Hugo_Symbol) %>% 
        dplyr::summarise("n_mut" = n()) %>%
        mutate(pvalue = unlist(llply(Hugo_Symbol, calculate_wilcox_by_gene, full_df, .parallel = T)))
}


gene_df <- create_gene_pvalue_df(full_df)


create_gene_pvalue_df_by_tumor <- function(tumor, full_df){
    tumor_df <- filter(full_df, tumor_type == tumor)
    gene_df <- tumor_df %>% 
        filter(!is.na(Hugo_Symbol)) %>% 
        group_by(Hugo_Symbol, tumor_type) %>% 
        dplyr::summarise("n_mut" = n()) 
    pvalues <- unlist(llply(gene_df$Hugo_Symbol, calculate_wilcox_by_gene, tumor_df, .parallel = T))
    gene_df$pvalue <- pvalues
    return(gene_df)
}

gene_df_by_tumor <- tumors %>% 
    map(create_gene_pvalue_df_by_tumor, full_df) %>% 
    bind_rows
    


setwd(output_dir)

write_tsv(gene_df, 'bootstrap_results/non_swi_gene_pvalues.tsv')
write_tsv(gene_df_by_tumor, 'bootstrap_results/non_swi_gene_pvalues_by_tumor.tsv')
