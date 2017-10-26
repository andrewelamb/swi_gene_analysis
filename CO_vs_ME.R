# This script computes Co-ocurrance vs Mututal exclusivity in swi genes in
# tcga and genie data sets

library(synapseClient)
library(tidyverse)
library(data.table)
library(magrittr)
library(stringr)
library(data.table)


working_dir     <- "/home/aelamb/repos/swi_gene_analysis"
swi_mut_file    <- "source_files/results-20170523-150114.csv"
swi_pat_file    <- "source_files/results-20170523-154238.csv"

setwd(working_dir)
synapseLogin()
synapseCacheDir("./")


# TCGA data
tcga_pat_df <- synGet("syn4983466")@filePath %>% 
    fread(select = c("bcr_patient_barcode", "acronym")) %>%
    as_data_frame
              

# SWI data --------------------------------------------------------------------
swi_mut_df <- fread(swi_mut_file)
swi_genes  <- unique(swi_mut_df$Hugo_Symbol)

# all tumors
swi_pat_df   <- swi_pat_file %>%
    fread %>% 
    as_data_frame %>% 
    rename(Tumor_Sample_Barcode = f0_) 

swi_pats <- swi_pat_df$Tumor_Sample_Barcode

# tumors with clinical data
swi_pat_df2 <- swi_pat_df %>%
    mutate(bcr_patient_barcode = str_sub(Tumor_Sample_Barcode, end = 12)) %>% 
    inner_join(tcga_pat_df) %>% 
    select(Tumor_Sample_Barcode, bcr_patient_barcode, acronym) %>% 
    arrange(acronym)

swi_pats2 <- swi_pat_df$Tumor_Sample_Barcode

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

# panels with swi genes
genie_panel_df <- genie_bed_df %>% 
    filter(Hugo_Symbol %in% swi_genes) %>% 
    select(SEQ_ASSAY_ID, Hugo_Symbol) %>% 
    distinct

genie_panel_sum_df1 <- genie_panel_df %>% 
    group_by(SEQ_ASSAY_ID) %>%
    summarise(swi_genes = n()) 

genie_panel_sum_df2 <- genie_pat_df %>% 
    group_by(SEQ_ASSAY_ID) %>% 
    summarise(patients = n()) 

# panel summary
genie_panel_sum_df <- left_join(genie_panel_sum_df1, genie_panel_sum_df2)

genie_pats <- genie_mut_df %>%
    use_series(Tumor_Sample_Barcode) %>%
    unique

genie_panels <- genie_panel_df %>%
    use_series(SEQ_ASSAY_ID) %>%
    unique

# matrix of zeros, patients, by swi genes -------------------------------------
swi_mut_matrix <- matrix(
    0, 
    nrow = length(swi_pats), 
    ncol = length(swi_genes), 
    dimnames = list(swi_pats, swi_genes))

genie_mut_matrix <- matrix(
    0, 
    nrow = length(genie_pats), 
    ncol = length(swi_genes), 
    dimnames = list(genie_pats, swi_genes))

# create a matrix of tumors associated with gene mutations
for(gene in swi_genes){
    # for each swi gene, get the tumors associated with that gene
    gene_str <- quo(gene)
    swi_patients_with_gene <- swi_mut_df %>% 
        filter(Hugo_Symbol == !!gene_str) %>% 
        use_series(Tumor_SampleBarcode)
    genie_patients_with_gene <- genie_mut_df %>% 
        filter(Hugo_Symbol == !!gene_str) %>% 
        use_series(Tumor_Sample_Barcode)
    swi_mut_matrix[swi_pats %in% swi_patients_with_gene, gene] <- 1
    genie_mut_matrix[genie_pats %in% genie_patients_with_gene, gene] <- 1
}

# replace 0s with NA's if panel doesn't include gene
for(panel in genie_panels){
    genie_tumors_with_panel <- genie_pat_df %>% 
        filter(SEQ_ASSAY_ID == panel) %>% 
        use_series(SAMPLE_ID) 
    genes_in_panel <- genie_bed_df %>% 
        filter(SEQ_ASSAY_ID == panel) %>% 
        filter(Hugo_Symbol %in% swi_genes) %>% 
        use_series(Hugo_Symbol) %>% 
        unique
    genie_mut_matrix[genie_pats %in% genie_tumors_with_panel, !swi_genes %in% genes_in_panel] <- NA
}

mut_matrix <- rbind(swi_mut_matrix, genie_mut_matrix)



# pw function -----------------------------------------------------------------


testPairwise <- function(M, threshold = 3, genes){
    ### pairwise comparison of co-occurence / mutual exclusivity
    n = length(genes)
    ME <- matrix(0, nrow=n, ncol=n, dimnames=list(genes, genes))
    for(i in 1:(n-1)){
        for(j in (i+1):n){
            if(sum(M[,i], na.rm = T) > 0 && sum(M[,j], na.rm = T) > 0){
                # co-occurance
                pval1 <- -log10(fisher.test(M[,i], M[,j],alternative="greater")$p.value)
                # exclusive
                pval2 <- -log10(fisher.test(M[,i], M[,j],alternative="less")$p.value)
                # positive values will represent occurance, negative, exclusivity
                ME[i,j] <- max(pval1, pval2) * ifelse(pval1 > pval2, 1, -1)
            }
        }
    }
    ME[abs(ME) < threshold] <- 0
    return(ME)
}

###



# all tumors 

pw_all_matrix <- testPairwise(mut_matrix, 0, swi_genes)
write.table(pw_all_matrix, "CO_ME_results/pw_all_matrix.tsv", quote = F, sep = "\t")
