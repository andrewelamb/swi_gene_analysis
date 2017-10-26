# This script did bootrapping on select swi models, and was run an an EC2

library(synapseClient)
library(glmnet)
library(limma)
library(ROCR)
library(tidyverse)
library(doMC)

registerDoMC(cores=10)
setwd('/home/ubuntu/')
synapseLogin()
synapseCacheDir("./")

# list of mutations
swiMuts <- read.delim2("results-20170523-150114.csv",sep=",",as.is=TRUE)
# list of tumor samples
swiPats <- read.delim2("results-20170523-154238.csv",sep=",",as.is=TRUE)[,1]

# clinical data, patient by row
print("get synapse TCGA data")
clin <- read.delim2(synGet("syn4983466")@filePath,sep="\t",header=TRUE,as.is=TRUE)
# expression data, genes by row, tumor sample by col
expr <- read.delim2(synGet("syn4976369")@filePath, as.is=TRUE, sep="\t",header=TRUE,check.names=FALSE)


print("filtering data")
# remove gene names
tmp <- data.matrix(expr[,-1])
# logical vector for any rows with atleast 1 NA
mask <- apply(tmp, 1, function(x) any(is.na(x)))
# getting gene names from names column, filtering out any with NA's in row
genes <- gsub("(.*)\\|.*", "\\1", expr[,1])[!mask]
# getting entrez numbers from names column, filtering out any with NA's in row
entrez <- gsub(".*\\|(.*)", "\\1", expr[,1])[!mask]
# filtering out na rows
em <- tmp[!mask,]

rm(tmp)
rm(expr)

types <- as.numeric(gsub("TCGA-..-....-(..).*","\\1",colnames(em)))
# filtering expresison data for cols with type under 10
em <- em[, types < 10]
# shortening expr colnames
colnames(em) <-  gsub("(TCGA-..-....-...).*", "\\1", colnames(em))


# unique genes from mutation list 
swiGenes <- unique(swiMuts$Hugo_Symbol)
n <- length(swiGenes)
# matrix of zeros, patients, by complex genes
M <- matrix(0,nrow=length(swiPats),ncol=n, dimnames=list(swiPats, swiGenes))
# create a matrix of tumors associated with gene mutations
for(gene in swiGenes){
    # for each complex gene, get the tumors associated with that gene
    pats <- swiMuts$Tumor_SampleBarcode[swiMuts$Hugo_Symbol == gene]
    M[swiPats %in% pats,gene] <- 1
}

# patients with tumors, clinical data, and expression data
commonPats <- intersect(intersect(clin$bcr_patient_barcode, 
                                  gsub("(TCGA-..-....).*","\\1", rownames(M))),
                        gsub("(TCGA-..-....).*","\\1", colnames(em)))

# index of clinical data rows to use
clin.idxs <- match(commonPats, clin$bcr_patient_barcode)
# index of mutation data rows to use
mut.idxs <- match(commonPats,gsub("(TCGA-..-....).*","\\1", rownames(M)))
# index of expression data cols to use
expr.idxs <- match(commonPats, gsub("(TCGA-..-....).*","\\1", colnames(em)))

M.m <- M[mut.idxs,]
clin.m <- clin[clin.idxs,]
expr.m <- em[, expr.idxs]

# for laptop
rm(em)
rm(clin)

#######################
####
# unique tumor types
tt <- unique(clin.m$acronym)

R <- lapply(tt, function(x){
    # filter data for tumor type
    mask <- clin.m$acronym == x
    # rows of M.m that are of this tumor and have atleast one complex gene mutation
    mut <- rowSums(M.m[mask,]) > 0
    # expression data for this tumor
    X <- expr.m[,mask]
    # Is gene differentially expressed in tumors grouped on 0/ >0 mutations in complex genes
    pvals <- eBayes(lmFit(X, model.matrix(~mut)))$p.value[,2]  
    # list n tumors, gene pvalues, n tumors with >0 complex gene mutations
    list(N=sum(mask),pvals=pvals,Nmut=sum(mut))
})
pdf("./ExprSignal.pdf",width=10,height=30)
par(mfrow=c(8,4))
for(i in 1:32){
    tmp <- R[[i]]
    hist(tmp$pvals,main=paste(tt[i],tmp$Nmut,"of",tmp$N))
}
dev.off()


#####################
# ridge modeling

boostrap_samples <- function(X, mut){
    print("bootstrapping")
    samples <- colnames(X)[colnames(X) %in% names(mut)]
    training_samples <- sample(x = samples, size = as.integer(.67 * ncol(X)), replace = T)
    test_samples <- samples[!samples %in% training_samples]
    training_exp <- t(X[,training_samples])
    test_exp <- t(X[, test_samples])
    training_mut <- mut[training_samples]
    test_mut <- mut[test_samples]
    cv.fit <- cv.glmnet(training_exp, training_mut, family="binomial",alpha = 0,type.measure="auc",parallel=TRUE)
    yhat <- predict(cv.fit, test_exp, type="response", s = "lambda.min") 
    yhat_df <- yhat %>% 
        data.frame %>% 
        rownames_to_column("sample")
    pred   <- prediction(yhat, test_mut)
    auc    <- performance(pred,"auc")@y.values[[1]]
    return(list(yhat_df, auc))
}


# pick tumor types with atleast 15 tumors and atleast 15% of total with one or more
# mutations in complex  genes
tt_to_model <- tt[which(sapply(R, function(x) x$Nmut > 15 && x$Nmut / x$N > .15))]


boostrap_tumor <- function(tumorType, M.m, clin.m, expr.m){
    print(paste("bootstrapping", tumorType))
    mask <- clin.m$acronym == tumorType
    # whether or not tumors have atleast one complex tumor
    mut <- rowSums(M.m[mask,]) > 0
    X <- expr.m[,mask]
    
    result_list <- rerun(30, boostrap_samples(X, mut)) 
    return(result_list)
}

results <- map(tt_to_model, boostrap_tumor, M.m, clin.m, expr.m)
names(results) <- tt_to_model
saveRDS(results, file = 'bootstrap_results.RDS')

