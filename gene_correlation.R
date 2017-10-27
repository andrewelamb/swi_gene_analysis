library(data.table)
library(tidyverse)
library(stringr)
library(magrittr)
library(gplots)
library(grid)


working_dir      <- "/home/aelamb/repos/swi_gene_analysis/"
boostrap_file    <- "bootstrap_results/bootstrap_scores.tsv"
image_dir        <- "images/correlation_images"

setwd(working_dir)

full_df <- fread(boostrap_file) %>% 
    as_data_frame



mut_df <- filter(full_df, !Hugo_Symbol == "none")

n_mut_cor_df <- mut_df %>% 
    group_by(tumor_type, Hugo_Symbol) %>% 
    summarise('cor' = cor(n_mut, med_score), 'n' = n())



mut_types <- sort(unique(mut_df$Variant_Classification))

mut_model_df <- mut_df %>% 
    model.matrix(~Variant_Classification+0, .) %>% 
    as_data_frame %>% 
    set_colnames(mut_types)

mut_cor_df1 <- mut_df %>% 
    group_by(tumor_type, Hugo_Symbol, Variant_Classification) %>% 
    summarise(n = n())

mut_cor_df2 <- mut_df %>% 
    select(tumor_type, Hugo_Symbol, med_score) %>% 
    bind_cols(mut_model_df) %>% 
    group_by(tumor_type, Hugo_Symbol) %>% 
    summarise('Frame_Shift_Del'        = cor(Frame_Shift_Del, med_score),
              'Frame_Shift_Ins'        = cor(Frame_Shift_Ins, med_score),
              'In_Frame_Del'           = cor(In_Frame_Del, med_score),
              'In_Frame_Ins'           = cor(In_Frame_Ins, med_score),
              'Missense_Mutation'      = cor(Missense_Mutation, med_score),
              'Nonsense_Mutation'      = cor(Nonsense_Mutation, med_score),
              'Nonstop_Mutation'       = cor(Nonstop_Mutation, med_score),
              'Splice_Site'            = cor(Splice_Site, med_score),
              'Translation_Start_Site' = cor(Translation_Start_Site, med_score)
)



heatmap_df <- n_mut_cor_df %>% 
    select(-n) %>% 
    rename("Mutation_Count" = cor) %>% 
    left_join(mut_cor_df2) %>% 
    ungroup

remove_empty_rows <- function(m){
    return(m[rowSums(is.na(m))!=ncol(m), ])
}

set_na_to_zero <- function(m){
    m[is.na(m)] <- 0
    return(m)
}


setwd(image_dir)

# all tumors
jpeg("all_tumors_heatmap.jpg", quality = 100, width = 1200, height = 900)
heatmap_df %>% 
    mutate(rowname = str_c(tumor_type, Hugo_Symbol, sep = ":")) %>% 
    select(-Hugo_Symbol, -tumor_type) %>% 
    column_to_rownames("rowname") %>% 
    as.matrix %>% 
    remove_empty_rows %>% 
    set_na_to_zero %>% 
    heatmap.2(
        trace = 'none', 
        main = 'All tumors',
        srtCol = 45, 
        margins=c(8,8))
dev.off()

make_heatmap_by_tumor <- function(df){
    tumor_type <- df$tumor_type[[1]]
    jpeg(str_c(
        tumor_type, "_heatmap.jpg"), 
        quality = 100, 
        width = 1200, 
        height = 900)
    df %>% 
        select(-tumor_type) %>% 
        column_to_rownames("Hugo_Symbol") %>% 
        as.matrix %>% 
        remove_empty_rows %>% 
        set_na_to_zero %>% 
        heatmap.2(
            trace = 'none', 
            main = tumor_type, 
            srtCol = 45, 
            margins=c(8,8))
    dev.off()
}

heatmap_df %>% 
    split(.$tumor_type) %>% 
    walk(make_heatmap_by_tumor)
    

    

