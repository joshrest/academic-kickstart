#source("/Users/jrest/Documents/GitHub/My_Website/content/post/The-Unraveling-of-a-Protist-Genome/clusters.R")

 library(data.table)
 library(stringr)
 library(clusterProfiler)
 library("org.EcK12.eg.db")

trtab <- fread("/Users/jrest/Documents/personal/financial/idmapping_2023_11_02.tsv",header=TRUE)

clus1 <- fread("/Users/jrest/Downloads/targets_cluster.csv",header=TRUE)
clus2 <- lapply(clus1$targets,function(dx)strsplit(dx,split=",")[[1]])
clus3 <- lapply(clus2,function(dx) tolower(str_replace_all(dx,"ECO:","")))
clus4 <- lapply(clus2,function(dx) trtab$RefSeq[match(unlist(dx), trtab$From)])
clus5 <- lapply(clus4,function(dx){
   sapply(dx, function(ex){
    strsplit(ex,split=";")[[1]][1]
      })
    })

clus6 <- lapply(clus5,function(dx){
   sapply(dx, function(ex){
    strsplit(ex,split=".",fixed=TRUE)[[1]][1]
      })
    })
names(clus6) <- c(paste0("cluster",1:16),"all")
#x <- enrichKEGG(clus4[[17]],organism="eco",keyType="uniprot")

for(i in 1:16){
print(i)
enrichGO(clus6[[i]],keyType="REFSEQ",OrgDb="org.EcK12.eg.db")
}

xx <- compareCluster(clus6[1:16],fun="enrichGO",keyType="REFSEQ",OrgDb="org.EcK12.eg.db")
pdf("compareTargets.pdf")
dotplot(xx,includeAll=TRUE, showCategory=50,font.size = 4)
dev.off()
