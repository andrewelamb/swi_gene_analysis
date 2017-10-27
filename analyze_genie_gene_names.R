# This script was used on an EC2 instance to create gene and variant tables on
# genie data. The results were put in 

library(synapseClient)
library(tidyverse)
library(data.table)
library(magrittr)
library(stringr)
library(grid)

synapseLogin()
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

panels <- sort(unique(genie_pat_df$SEQ_ASSAY_ID)) 

genie_full_mut_df <-  genie_mut_df %>% 
    left_join(genie_pat_df, by = c("Tumor_Sample_Barcode" = "SAMPLE_ID"))

get_mut_genes_not_in_panel <- function(panel, genie_bed_df, genie_full_mut_df){
    panel_genes <- genie_bed_df %>% 
        filter(SEQ_ASSAY_ID == panel) %>% 
        use_series("Hugo_Symbol") %>% 
        unique %>% 
        sort
    sample_genes <- genie_full_mut_df %>% 
        filter(SEQ_ASSAY_ID == panel) %>% 
        use_series("Hugo_Symbol") %>% 
        unique %>% 
        sort
    genes <- sample_genes[!sample_genes %in% panel_genes]
}

genie_panel_df1 <- genie_bed_df %>% 
    group_by(SEQ_ASSAY_ID) %>% 
    summarise("panel_gene_count" = n())

genie_panel_df2 <- genie_full_mut_df %>% 
    select(-HGVSp_Short) %>% 
    select(-Tumor_Sample_Barcode) %>% 
    distinct %>% 
    group_by(SEQ_ASSAY_ID) %>% 
    summarise("sample_gene_count" = n())

non_panel_genes <- 
    map(panels, 
        get_mut_genes_not_in_panel, 
        genie_bed_df, 
        genie_full_mut_df) %>% 
    set_names(panels)

genie_panel_df3 <- data_frame(
    "SEQ_ASSAY_ID" = panels,
    "non_panel_genes_count" = map(non_panel_genes, length),
    "non_panel_genes" = map(non_panel_genes, str_c, collapse = ","))
               
genie_panel_df <- reduce(list(genie_panel_df1, genie_panel_df2, genie_panel_df3), full_join)

jpeg("images/genie_gene_names_tbl.jpg", quality = 100, width = 1200, height = 900)
grid.table(genie_panel_df)
dev.off()
    
