Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

display.progress = function ( index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}

# PPI prepare -------------------------------------------------------------

library(igraph)
library(rio)
library(tidyverse)

# ==== load network data ====
# 1. HumanNet-PI ====
# HumanNet_PI <- import("./Interactome2022.csv")
#
# HumanNet_PI2 <- HumanNet_PI
# HumanNet_PI2$LLS <- 1
#
# nodes <- as.data.frame(union(unique(HumanNet_PI2$EntrezA), unique(HumanNet_PI2$EntrezB)))
# names(nodes) <- "name"
# HumanNet_PI3 <- graph_from_data_frame(HumanNet_PI2, directed = FALSE, vertices = nodes)
#
# print(c(length(V(HumanNet_PI3)), gsize(HumanNet_PI3)))
#
#
# # largest connected component (LCC)
# HumanNet_PI3_LCC <- decompose.graph(HumanNet_PI3, mode = "weak")[[1]]
# print(c(length(V(HumanNet_PI3_LCC)), gsize(HumanNet_PI3_LCC)))
#
# save(HumanNet_PI3_LCC,file = 'background_net.RData')

library(dnet)
library(xlsx)

# library(clusterProfiler)
# load TNBC genes

library(clusterProfiler)

library(org.Hs.eg.db)


# load network data ====
load("./background_net.RData")

# PPI networks
net <- HumanNet_PI3_LCC
set <- V(net)$name
# RWR for TNBC_driver genes ====

tnbc_gene <- import('./tnbc_driver_final.csv')

driver_tnbc_gene <- tnbc_gene$ENTREZID

# define set of seeds
seeds_tnbc <- rep(1, length(driver_tnbc_gene))
seeds_tnbc <- data.frame(seeds_tnbc)
rownames(seeds_tnbc) <- driver_tnbc_gene

n1 <- length(intersect(rownames(seeds_tnbc), set))

# Natural Product ---------------------------------------------------------

npro <- import('tcm.csv')

length(unique(npro$Pubchem_CID))
npro_gene_entrez <- bitr(geneID = as.character(npro$Targets),
                         fromType = "SYMBOL",
                         toType = "ENTREZID",
                         OrgDb = "org.Hs.eg.db")

print(paste0("===实际匹配到entrez编码基因个数：",
             length(na.omit(npro_gene_entrez$ENTREZID)),
             "==="))

npro_final <- inner_join(npro,npro_gene_entrez,by=c('Targets'='SYMBOL'))


# npro overlap with PPI ---------------------------------------------------

npro_final2 <- npro_final[npro_final$ENTREZID%in%set,]

drug_num <- length(unique(npro_final2$Pubchem_CID))

#write.csv(npro_final,file = 'npro_final.csv',quote = F,row.names = F)

# choose α ----------------------------------------------------------------


# 改写 ----------------------------------------------------------------------


PTmatrix_disease <- suppressWarnings(suppressMessages(dRWR(g = net, normalise = "row",
                                                           setSeeds = seeds_tnbc, restart = 0.75,
                                                           parallel = FALSE)))

#
# disease_random <- lapply(1:1000, function(x){
#   seeds_disease_random_net <- rep(1, n1)
#   seeds_disease_random_net <- data.frame(seeds_disease_random_net)
#   row.names(seeds_disease_random_net) <- set[sample(1:length(set), n1)]
#
#   PTmatrix_disease_random_net <- suppressWarnings(suppressMessages(dRWR(g = net, normalise = "row",
#                                       setSeeds = seeds_disease_random_net, restart = 0.75,
#                                       parallel = FALSE)))
# }
# )
#
# save(disease_random,file = 'disease_random_rwr.Rdata')

load(file = 'disease_random_rwr.Rdata')

z_score <- NULL

# library(future.apply)
# future::plan(multiprocess(workers = 48))
#
#
# t1<-Sys.time()
#
# ll <- suppressWarnings(future_lapply(1:drug_num, function(i){
#
#   display.progress(i,drug_num)
#
#   CID = unique(npro_final$Pubchem_CID)[i]
#
#   CID_target= npro_final[npro_final$Pubchem_CID%in%CID,"ENTREZID"]
#
#   seeds_target <- rep(1, length(CID_target))
#
#   seeds_target <- data.frame(seeds_target)
#
#   rownames(seeds_target) <-  CID_target
#
#   n2 <- length(intersect(rownames(seeds_target), set))
#
#   PTmatrix_target <- suppressWarnings(suppressMessages(dRWR(g = net, normalise = "row",
#                                                              setSeeds = seeds_target, restart = 0.75,
#                                                              parallel = FALSE)))
#
#   pcc=cor(as.numeric(PTmatrix_disease), as.numeric(PTmatrix_target), method = "pearson")
#
#   drug_random=future_sapply(1:1000, function(j){
#
#     seeds_target_random_net <- rep(1, n2)
#     seeds_target_random_net <- data.frame(seeds_target_random_net)
#     row.names(seeds_target_random_net) <- set[sample(1:length(set), n2)]
#
#     PTmatrix_target_random_net <- suppressWarnings(suppressMessages(dRWR(g = net, normalise = "row",
#                                        setSeeds = seeds_target_random_net, restart = 0.75,
#                                        parallel = FALSE)))
#
#     pcc_seed_random <- cor(as.numeric(disease_random[j][[1]]), as.numeric(PTmatrix_target_random_net), method = "pearson")
#
#   })
#
#   p = sum(pcc<drug_random)/1000
#
#   mean_pcc_seed_random <- mean(drug_random)
#
#   std_pcc_seed_random <- sd(drug_random)
#
#   z_score <- (pcc - mean_pcc_seed_random)/std_pcc_seed_random
#
#   res=list(CID =CID,cor=pcc,pvalue=p,Z=z_score)
#
# })
# )
#
# t2<-Sys.time()
#
# time <- t2-t1
#
# drug_rwr <- data.frame(do.call(rbind,ll),check.rows = F,check.names = F)


library(furrr)
plan(multisession, workers = availableCores()-30)

options(future.globals.maxSize= 3221222500)

system.time(
  result <- furrr::future_map(1:drug_num, function(i) {

    CID = unique(npro_final2$Pubchem_CID)[i]

    CID_target= unique(npro_final2[npro_final2$Pubchem_CID%in%CID,"ENTREZID"])

    seeds_target <- rep(1, length(CID_target))

    seeds_target <- data.frame(seeds_target)

    rownames(seeds_target) <-  CID_target

    n2 <- length(intersect(rownames(seeds_target), set))

    PTmatrix_target <- suppressWarnings(suppressMessages(dRWR(g = net, normalise = "row",
                                                              setSeeds = seeds_target, restart = 0.75,
                                                              parallel = TRUE)))

    pcc=cor(as.numeric(PTmatrix_disease), as.numeric(PTmatrix_target), method = "pearson")

      rv <- purrr::map(1:1000, function(j) {
      seeds_target_random_net <- rep(1, n2)
      seeds_target_random_net <- data.frame(seeds_target_random_net)
      row.names(seeds_target_random_net) <- set[sample(1:length(set), n2)]

      PTmatrix_target_random_net <- suppressWarnings(suppressMessages(dRWR(g = net, normalise = "row",
                                                                           setSeeds = seeds_target_random_net, restart = 0.75,
                                                                           parallel = TRUE)))

      pcc_seed_random <- cor(as.numeric(disease_random[j][[1]]), as.numeric(PTmatrix_target_random_net), method = "pearson")
      pcc_seed_random
    })


    drug_random <- do.call(cbind,rv)

    p = sum(pcc<drug_random)/1000

    mean_pcc_seed_random <- mean(drug_random)

    std_pcc_seed_random <- sd(drug_random)

    z_score <- (pcc - mean_pcc_seed_random)/std_pcc_seed_random

    res=list(CID =CID,cor=pcc,pvalue=p,Z=z_score)
    res
  }, .progress = TRUE)
)

save(result,file = 'new_drug_rwr.Rdata')
drug_rwr <- data.frame(do.call(rbind,result),check.rows = F,check.names = F)

save(drug_rwr,file = 'update_tnbc_drug_rwr.Rdata')

drug_rwr$CID = sapply(drug_rwr$CID, function(x) x[1])
drug_rwr$cor = sapply(drug_rwr$cor, function(x) x[1])
drug_rwr$pvalue = sapply(drug_rwr$pvalue, function(x) x[1])
drug_rwr$Z = sapply(drug_rwr$Z, function(x) x[1])




load(file = 'new_drug_rwr.Rdata')

drug_rwr3 <- drug_rwr[drug_rwr$pvalue<0.05&drug_rwr$Z>3,]

drug_rwr4 <- drug_rwr[drug_rwr2$pvalue<0.05&drug_rwr2$Z>3,]

length(intersect(drug_rwr3$CID,drug_rwr4$CID))

save(drug_rwr4,file = 'tnbc_drug_rwr.Rdata')

#load(file = 'drug_rwr.Rdata')





#
#
# # calculate affinity matrix
#
# PTmatrix_tnbc_seed <- dRWR(g = net, normalise = "row",
#                            setSeeds = seeds_tnbc, restart = 0.75,
#                            parallel = FALSE)
#
#
# set <- V(net)$name
#
# n1 <- length(intersect(rownames(seeds_tnbc), set))
#
# z_score <- NULL
#
# p_value <- NULL
#
# j <- 0
#
# pcc_seed_random <- NULL
#
# drug_num <- length(unique(npro_final$Pubchem_CID))
#
# for (i in 1:2) {
#   CID = unique(npro_final$Pubchem_CID)[i]
#
#   CID_target= npro_final[npro_final$Pubchem_CID%in%CID,"ENTREZID"]
#
#   seeds_target <- rep(1, length(CID_target))
#
#   seeds_target <- data.frame(seeds_target)
#
#   rownames(seeds_target) <-  CID_target
#
#   n2 <- length(intersect(rownames(seeds_target), set))
#
#
#   PTmatrix_target_net <- dRWR(g = net, normalise = "row",
#                               setSeeds = seeds_target, restart = 0.75,
#                               parallel = FALSE)
#
#   pcc=cor(as.numeric(PTmatrix_tnbc_seed), as.numeric(PTmatrix_target_net), method = "pearson")
#
#   # repeated times
#   for (j in 1:1000) {
#
#     seeds_disease_random_net <- rep(1, n1)
#     seeds_disease_random_net <- data.frame(seeds_disease_random_net)
#     row.names(seeds_disease_random_net) <- set[sample(1:length(set), n1)]
#
#     seeds_target_random_net <- rep(1, n2)
#     seeds_target_random_net <- data.frame(seeds_target_random_net)
#     row.names(seeds_target_random_net) <- set[sample(1:length(set), n2)]
#
#     PTmatrix_disease_random_net <- dRWR(g = net, normalise = "row",
#                                         setSeeds = seeds_disease_random_net, restart = 0.75,
#                                         parallel = FALSE)
#
#     PTmatrix_target_random_net <- dRWR(g = net, normalise = "row",
#                                        setSeeds = seeds_target_random_net, restart = 0.75,
#                                        parallel = FALSE)
#
#     pcc_seed_random <- c(pcc_seed_random, cor(as.numeric(PTmatrix_disease_random_net), as.numeric(PTmatrix_target_random_net), method = "pearson"))
#
#     p_value_seed <- sum(pcc<pcc_seed_random)/1000
#
#   }
#   mean_pcc_seed_random <- mean(pcc_seed_random)
#   std_pcc_seed_random <- sd(pcc_seed_random)
#
#   z_score_seed <- (pcc - mean_pcc_seed_random)/std_pcc_seed_random
#
#   z_score <- c(z_score,z_score_seed)
#
#   p_value <- c(p_value,z_score_seed)
# }
#
#
#
#
