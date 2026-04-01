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
library(dnet)
library(xlsx)

library(clusterProfiler)

library(org.Hs.eg.db)


# load network data ====
load("./background_net.RData")

# PPI networks
net <- HumanNet_PI3_LCC
set <- V(net)$name
# RWR for TNBC_driver genes ====

tnbc_gene <- import('./TNBC_driver.csv')

driver_tnbc_gene <- tnbc_gene$ENTREZID

# define set of seeds
seeds_tnbc <- rep(1, length(driver_tnbc_gene))
seeds_tnbc <- data.frame(seeds_tnbc)
rownames(seeds_tnbc) <- driver_tnbc_gene

n1 <- length(intersect(rownames(seeds_tnbc), set))

# TCM-derived candidate compounds  ---------------------------------------------------------

npro <- import('TCM ingredient2targets.csv')

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



load(file = 'disease_random_rwr.Rdata')

z_score <- NULL

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


