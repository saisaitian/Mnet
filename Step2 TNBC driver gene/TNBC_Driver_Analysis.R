
# TNBC network ------------------------------------------------------------

library(igraph)
# load network data ====
load("./background_net.RData")

# PPI networks
net <- HumanNet_PI3_LCC
set <- V(net)$name
# RWR for TNBC_driver genes ====
library(rio)
tnbc_gene <- import('./tnbc_driver_final.csv')

driver_tnbc_gene <- tnbc_gene$ENTREZID

tnbc_PPI = driver_tnbc_gene[driver_tnbc_gene %in% V(net)$name]

require(NetSci)

gTnbc = net %>%
  induced.subgraph(., as.character(tnbc_PPI))%>%igraph::simplify()

components(gTnbc)
components(gTnbc)$csize %>% max

# 提取最大联通片网络 ---------------------------------------------------------------

V(gTnbc)$size = degree(gTnbc) %>%
  CoDiNA::normalize()
V(gTnbc)$size = (V(gTnbc)$size + 0.1)*5
V(gTnbc)$color = '#83c5be'
V(gTnbc)$frame.color = '#006d77'
V(gTnbc)$label = ifelse(V(gTnbc)$size  > 4, V(gTnbc)$name, NA )
V(gTnbc)$label.color = '#e29578'

E(gTnbc)$width = edge.betweenness(gTnbc, directed = F) %>% CoDiNA::normalize()
E(gTnbc)$width = E(gTnbc)$width + 0.01
E(gTnbc)$weight = E(gTnbc)$width
par(mar = c(0,0,0,0))
plot(gTnbc)

require(magrittr)
library(clusterProfiler)
library(dplyr)
# LCC
gTnbc_LCC <- decompose.graph(gTnbc, mode = "weak")[[2]]

components(gTnbc_LCC)$csize %>% max

V(gTnbc_LCC)$name


LCC_tnbc = LCC_Significance(N = 10000, Targets = V(gTnbc_LCC)$name,
                           G = net,bins = 100)

Histogram_LCC(LCC_tnbc)

LCC_tnbc_net <- as.data.frame(get.edgelist(gTnbc_LCC) )


LCC_tnbc_net_node <- unique(c(LCC_tnbc_net$V1,LCC_tnbc_net$V2))


eg <- bitr(LCC_tnbc_net$V1,fromType = 'ENTREZID',
           toType = c('SYMBOL'),
           OrgDb='org.Hs.eg.db',
)


eg2 <- bitr(LCC_tnbc_net$V2,fromType = 'ENTREZID',
            toType = c('SYMBOL'),
            OrgDb='org.Hs.eg.db',
)

eg3 <- bitr(LCC_tnbc_net_node,fromType = 'ENTREZID',
            toType = c('SYMBOL'),
            OrgDb='org.Hs.eg.db',
)

LCC_tnbc_net_data <- inner_join(eg,LCC_tnbc_net,by=c('ENTREZID'= 'V1'))

LCC_tnbc_net_data2 <- inner_join(eg2,LCC_tnbc_net_data,by=c('ENTREZID'= 'V2'))

LCC_tnbc_net_data3 <- LCC_tnbc_net_data2[,c(2,4)]

names(LCC_tnbc_net_data3) <- c('Source','Target')

write.csv(LCC_tnbc_net_data3,file = 'LCC_tnbc_net.csv',quote = F,row.names = F)


