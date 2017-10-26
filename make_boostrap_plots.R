# This script makes plots form the swi bootstrap models
library(tidyverse)
library(magrittr)
library(gplots)
library(stringr)


working_dir  <- "/home/aelamb/repos/swi_gene_analysis/"
input_file   <- "bootstrap_results/bootstrap_results.RDS"
tumor_file   <- "source_files/results-20170523-150114.csv"
image_dir    <- "bootstrap_images/"

setwd(working_dir)
synapseLogin()
synapseCacheDir("./tmp/")


tumor_df      <- tumor_file %>% 
    read_csv %>% 
    select(-Study) %>% 
    mutate(bcr_patient_barcode = str_sub(Tumor_SampleBarcode, end = 12))

patient_df <- synGet("syn4983466")@filePath %>% 
    fread(select = c("bcr_patient_barcode", "acronym")) %>% 
    as_data_frame %>% 
    left_join(tumor_df) %>% 
    select(-Tumor_SampleBarcode)

results_list  <- readRDS(input_file)
tumors        <- sort(names(results_list))



# table of auc scores
get_aucs <- function(tumor){
    aucs <- c()
    for(bs in 1:30){
        aucs[bs] <- results_list[[tumor]][[bs]][[2]]
    }
    data_frame('aucs' = aucs, 'tumor_type' = tumor)
}

auc_df <- map(tumors, get_aucs) %>% 
    bind_rows

# table of mutations, their score, tumor type and gene
get_scores <- function(tumor){
    scores <- list()
    for(bs in 1:30){
        scores[bs] <- results_list[[tumor]][[bs]][1]
    }
    bind_rows(scores)
}

score_df <- map(tumors, get_scores) %>% 
    bind_rows %>% 
    mutate(bcr_patient_barcode = str_sub(sample, end = 12)) %>% 
    left_join(patient_df) %>%  
    rename("score" = X1) %>% 
    rename("tumor_type" = acronym) %>% 
    select(-bcr_patient_barcode)
    


# seprate samples into ones with mutations, and wt
mut_score_df <- score_df %>% 
    filter(!is.na(Hugo_Symbol))

wt_score_df <- score_df %>% 
    filter(is.na(Hugo_Symbol))

# table of samples with at least one mutation
mutation_df <- mut_score_df %>% 
    select(-score, -tumor_type) %>% 
    distinct 

nmut_df <- mut_score_df %>% 
    select(-score) %>% 
    distinct %>% 
    group_by(sample) %>%
    summarise(n_mut = n())

med_score_df <- mut_score_df %>% 
    group_by(sample) %>% 
    summarise(med_score = median(score))

mut_sample_df <- mut_score_df %>% 
    select(sample, tumor_type) %>% 
    distinct %>% 
    left_join(nmut_df) %>%
    left_join(med_score_df)

# wt samples
wt_sample_df <- wt_score_df %>% 
    select(-Variant_Classification, -Hugo_Symbol) %>% 
    group_by(sample, tumor_type) %>% 
    summarise(med_score = median(score)) %>% 
    inset("n_mut", value = 0) %>% 
    inset("Variant_Classification", value = "none") %>% 
    inset("Hugo_Symbol", value = "none")

# all samples
df <- left_join(mut_sample_df, mutation_df) %>% 
    bind_rows(wt_sample_df) %>% 
    mutate(Hugo_Symbol = as.factor(Hugo_Symbol)) %>% 
    as_data_frame

write_tsv(df, "boostrap_results/bootstrap_scores.tsv")
write_tsv(auc_df, "boostrap_results/auc_scores.tsv")

df$Hugo_Symbol <- relevel(df$Hugo_Symbol, 'none')

levs <- levels(df$Hugo_Symbol)


# median score by gene and tumor type
genes <- sort(unique(df$Hugo_Symbol))

median_score_matrix <- matrix(
    0, 
    nrow = length(tumors),
    ncol = length(genes), 
    dimnames = list(tumors, genes))

for(tumor in tumors){
    for (gene in genes){
        score <- df %>% 
            filter(tumor_type == tumor) %>% 
            filter(Hugo_Symbol == gene) %>% 
            use_series(med_score) %>% 
            median()
        median_score_matrix[tumor, gene] <- score
    }
}

# ratio of median scores to wt scores
score_ratio_matrix <- apply(median_score_matrix, 2, function(col) col / median_score_matrix[,'none']) %>% 
    .[, -1]

# percent of tumors with mutation in gene per tumor type
genes <- sort(unique(mutation_df$Hugo_Symbol))

percent_mutated_matrix <- matrix(
    0, 
    nrow = length(tumors),
    ncol = length(genes), 
    dimnames = list(tumors, genes))

for(tumor in tumors){
    tbl <- filter(df, tumor_type == tumor)
    n_samples <- tbl %>% 
        select(sample) %>% 
        distinct %>% 
        nrow
    for (gene in genes){
        n_with_gene <- tbl %>% 
            filter(Hugo_Symbol == gene) %>% 
            nrow
        if(n_with_gene == 0){
            percent_mutated_matrix[tumor, gene] <- NA
        }
        else{
            percent_mutated_matrix[tumor, gene] <- n_with_gene/n_samples
        }
    }
}



set_wd(image_dir)

# ratio score heatmap
jpeg("score_ratio_heatmap.jpg", quality = 100, width = 700, height = 900)
heatmap.2(score_ratio_matrix, trace = 'none', main = 'Score Ratio', na.color = 'blue')
dev.off()


# percent gene by tumor heatmap
jpeg("percent_gene_heatmap.jpg", quality = 100, width = 1200, height = 900)
heatmap.2(percent_mutated_matrix, trace = 'none', main = '% Tumors with gene mutation', na.color = 'blue')
dev.off()



jpeg("aucs_boxplot.jpg", quality = 100, width = 500, height = 900)
ggplot(auc_df, aes(x = tumor_type, y = aucs)) +
    scale_colour_gradientn(colours = rainbow(10)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) +
    labs(title = "AUC") +
    scale_y_continuous(limits = c(0.0, 1.0)) %>% 
    print
dev.off()

# plot by gene

empty_df <- df %>% 
    select(Hugo_Symbol) %>% 
    distinct

# function for number of observations 
give.n <- function(x){
    return(c(y = 1.0, label = length(x))) 
}

give.p <- function(x){
    return(c(y = 0.0, label = round(length(x)/ nrow(scores_df), digits = 2)))
}

# aucs plot
plot_aucs_by_tumor <- function(df){
    p <- ggplot(df, aes(x = tumor_type, y = aucs)) +
        scale_colour_gradientn(colours = rainbow(10)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) +
        labs(title = "AUC") +
        scale_y_continuous(limits = c(0.0, 1.0))
    print(p)
}

plot_scores_by_tumor <- function(scores_df){
    p <- ggplot(scores_df, aes(x = Hugo_Symbol, y = med_score)) +
        scale_colour_gradientn(colours = rainbow(10)) +
        scale_shape_manual(values = c(15:18, 7:10 )) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(aes(color = n_mut, shape = Variant_Classification), size = 4) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20)) +
        labs(title = scores_df$tumor_type[[1]]) +
        scale_y_continuous(limits = c(0.0, 1.0)) +
        scale_x_discrete(labels = levs) + 
        stat_summary(fun.data = give.n, geom = "text") +
        stat_summary(fun.data = give.p, geom = "text")
    print(p)
}


for (tumor in tumors){
    aucs   <- filter(auc_df, tumor_type == tumor)
    scores_df <- filter(df, tumor_type == tumor) %>% 
        bind_rows(empty_df)
    jpeg(str_c(tumor, "_score_boxplot.jpg"), quality = 100, width = 1000, height = 900)
    plot_scores_by_tumor(scores_df)
    dev.off()
    jpeg(str_c(tumor, "_auc_boxplot.jpg"), quality = 100, width = 200, height = 900)
    plot_aucs_by_tumor(aucs)
    dev.off()
}






