library(Seurat)
library(SLICE)
#R version 3.5.1 (2018-07-02)
# Seurat_3.1.1 SLICE_0.99.0

load('data/obj.Rda')
data <- object@assays$RNA@data
name <- object@misc$fdata$name

used  <- rowSums(as.matrix(data)) > 0
data2 <- data[used, ]
name2 <- name[used]

used2 <- (! duplicated(name2) )
data3 <- data2[used2, ]
rownames(data3) <- name2[used2]

saveRDS(data3, file = 'data/data.Rds')
data <- data3 #readRDS('data/data.Rds')

metadata <- read.table("data/UMAP.data.xls", sep = "\t", header = T, row = 1)

rb.genes <- grep("^Rpl|^Rps|^Mrpl|^Mrps|_rRNA", rownames(data), value = TRUE)
es <- data[!rownames(data) %in% rb.genes, ]
str(es)

sc <- construct(exprmatrix=as.data.frame(as.matrix(es)), cellidentity=factor(metadata[colnames(es), 'Cluster']))

data(mm_kappasim)
rownames(km) <- toupper(rownames(km))
colnames(km) <- toupper(colnames(km))
str(km)

sc <- getEntropy(sc, km=km,                             # use the pre-computed kappa similarity matrix of mouse genes
                 calculation="bootstrap",               # choose the bootstrap calculation
                 B.num=100,                             # 100 iterations
                 exp.cutoff=1,                          # the threshold for expressed genes
                 B.size=1000,                           # the size of bootstrap sample
                 clustering.k=floor(sqrt(1000/2)),      # the number of functional clusters  
                 random.seed=201602)                    # set the random seed to reproduce the results in the paper

saveRDS(sc, 'sc.Rds')
entropy <- data.frame(cells = names(sc@entropy), entropy = sc@entropy)
write.table(entropy, file = 'entropy.xls', sep = "\t", row = F, quote = F)

