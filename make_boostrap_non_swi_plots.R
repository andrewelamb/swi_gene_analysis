library(tidyverse)
library(magrittr)
library(gplots)
library(stringr)
library(data.table)
library(qvalue)
library(synapseClient)

working_dir    <- "/home/aelamb/repos/swi_gene_analysis/"
non_swi_file   <- "bootstrap_results/non_swi_gene_pvalues_by_tumor.tsv"
bootstrap_file <- "bootstrap_results/bootstrap_scores.tsv"
image_dir      <- "bootstrap_images_non_swi/"


setwd(working_dir)
synapseLogin()
synapseCacheDir("./tmp/")


non_swi_df <- fread(non_swi_file)
qvals <- non_swi_df %>% 
    use_series(pvalue) %>% 
    qvalue %>% 
    use_series(qvalues)

non_swi_df <- non_swi_df %>% 
    as_data_frame %>% 
    inset('qvalue', value = qvals) %>% 
    arrange(qvalue, pvalue)

mut_df  <- synGet("syn4924181")@filePath %>% 
    str_c('zcat ', .) %>% 
    fread(
        select = c("Tumor_Sample_Barcode", "Variant_Classification", "Hugo_Symbol"),
        skip = 1) %>% 
    as_data_frame %>% 
    mutate(sample = str_sub(Tumor_Sample_Barcode, 1, 16)) %>% 
    select(-Tumor_Sample_Barcode) 

full_df <- fread(bootstrap_file) %>% 
    as_data_frame %>% 
    filter(Hugo_Symbol == "none") %>% 
    select(sample, tumor_type, med_score) %>% 
    left_join(mut_df)

n_muts <- full_df %>% 
    group_by(sample) %>% 
    summarise(n_muts = n())

full_df <- left_join(full_df, n_muts)

tumors <- sort(unique(full_df$tumor_type))

give.n <- function(x){
    return(c(y = 1.0, label = length(x))) 
}


plot_scores_by_tumor <- function(tumor, full_df, non_swi_df){
    tumor_df <- filter(full_df, tumor_type == tumor)
    top_genes_df <- non_swi_df %>% 
        filter(tumor_type == tumor) %>% 
        arrange(pvalue) %>% 
        slice(1:25) %>% 
        select(Hugo_Symbol, pvalue, qvalue) %>% 
        bind_rows(data_frame(Hugo_Symbol = "all")) %>% 
        arrange(Hugo_Symbol) 
    top_tumor_df <- tumor_df %>% 
        inner_join(top_genes_df)
    all_df <- tumor_df %>% 
        select(sample, med_score, n_muts, Variant_Classification) %>% 
        mutate(Variant_Classification = 'mult') %>% 
        distinct %>% 
        mutate(Hugo_Symbol = 'all')
    plot_df <- bind_rows(top_tumor_df, all_df)
    p <- ggplot(plot_df, aes(x = Hugo_Symbol, y = med_score)) +
        scale_colour_gradientn(colours = rainbow(10)) +
        scale_shape_manual(values = c(15:18, 7:15 )) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(aes(color = n_muts, shape = Variant_Classification)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) +
        labs(title = tumor) +
        scale_x_discrete(labels = paste(
            top_genes_df$Hugo_Symbol,
            "(", 
            sprintf("%1.1e", top_genes_df$pvalue),
            ")", 
            sep="")) +
        scale_y_continuous(limits = c(0.0, 1.0)) +
        stat_summary(fun.data = give.n, geom = "text")
    print(p)
}

setwd(image_dir)

for (tumor in tumors){
    jpeg(str_c(tumor, "_non_swi_score_boxplot.jpg"), quality = 100, width = 1000, height = 900)
    plot_scores_by_tumor(tumor, full_df, non_swi_df)
    dev.off()
}




