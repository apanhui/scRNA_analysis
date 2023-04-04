library(Seurat); library(dplyr); library(ggplot2); library(patchwork);library(DoubletFinder)

#################
### rawdata ###
#################
#"barcodes.tsv"   "genes.tsv"      "matrix.mtx"

## Create Seurat object
data_name = sample_name
object.list <- list()
for (i in seq(data_name)) {
	mat <- Read10X(data.dir = data_name[i], gene.column = 1)
	object.list[[i]] <- CreateSeuratObject(counts = mat, project = data_name[i], assay = assay)
}

## merge Seurat object
object <- merge(x = object.list[[1]], y = unlist(object.list[-1]),add.cell.ids = data_name)

## Normalize
object <- NormalizeData(object, normalization.method = "LogNormalize",scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst",nfeatures = 2000)
object <- ScaleData(object = object, vars.to.regress = NULL,features = rownames(object))

object <- RunPCA(object,npcs = 50,features = VariableFeatures(object), verbose = FALSE)
object <- RunUMAP(object, dims = 1:50, umap.method = "uwot")
object <- RunTSNE(object, dims = seq(50))
dims <- seq(object@reductions$pca)

### get pN
pN <- 0.25
### get pK
sweep.res.list <- paramSweep_v3(object, PCs = dims, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list)
bcmvn <- find.pK(sweep.stats)
pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)[1]]))

pdf(paste0("pK.", sample, ".pdf"))
plot(x = as.vector(bcmvn$pK), y = bcmvn$BCmetric, type = "b", xlab = "pK", ylab = "BCmetric", col = "blue", pch = 19)
abline(v = pK, lty = 2, col = "red")
dev.off()
WriteTable(bcmvn, file = paste0("pK.", sample, ".xls"))

### get nExp
if ( is.na(rate) )
    rate <- 7.6 * 10^-6 * ncol(object) + 5.27 * 10^-4
nExp_poi <- round(as.numeric(rate) * ncol(object))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
object <- doubletFinder_v3(object, PCs = dims, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
colnames(object@meta.data)[grep('pANN', colnames(object@meta.data))] <- "pANN"
colnames(object@meta.data)[grep('DF.classifications', colnames(object@meta.data))] <- "classifications"

embeddings <- cbind(object@reductions[["tsne"]]@cell.embeddings, object@reductions[["umap"]]@cell.embeddings)
data <- object@meta.data[,c("pANN", "classifications")]
WriteTable(cbind(Cells = rownames(data), data), file = paste0("DF.classify.", sample, ".xls"))
WriteTable(cbind(Cells = rownames(data), data, embeddings), file = paste0("DF.classify.", sample, ".tmp"))

sessionInfo()
#R version 3.5.1 (2018-07-02)
#Platform: x86_64-conda_cos6-linux-gnu (64-bit)
#Running under: CentOS release 6.9 (Final)

#Matrix products: default

# 
#locale:
# [1] LC_CTYPE=en_US.iso885915       LC_NUMERIC=C                  
# [3] LC_TIME=en_US.iso885915        LC_COLLATE=en_US.iso885915    
# [5] LC_MONETARY=en_US.iso885915    LC_MESSAGES=en_US.iso885915   
# [7] LC_PAPER=en_US.iso885915       LC_NAME=C                     
# [9] LC_ADDRESS=C                   LC_TELEPHONE=C                
#[11] LC_MEASUREMENT=en_US.iso885915 LC_IDENTIFICATION=C  
# 
#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     
# 
#other attached packages:
#[1] patchwork_1.0.0.9000 ggplot2_3.2.1        dplyr_0.8.5         
#[4] Seurat_3.1.1 DoubletFinder_2.0.3      
# 
#loaded via a namespace (and not attached):
# [1] nlme_3.1-147        tsne_0.1-3          bitops_1.0-6       
# [4] RcppAnnoy_0.0.14    RColorBrewer_1.1-2  httr_1.4.0         
# [7] sctransform_0.2.0   tools_3.5.1         R6_2.4.1           
#[10] irlba_2.3.3         KernSmooth_2.23-15  uwot_0.1.5         
#[13] lazyeval_0.2.2      colorspace_1.4-1    npsurv_0.4-0       
#[16] withr_2.1.2         tidyselect_1.0.0    gridExtra_2.3      
#[19] compiler_3.5.1      cli_2.0.2           plotly_4.8.0       
#[22] caTools_1.17.1.2    scales_1.0.0        lmtest_0.9-36      
#[25] ggridges_0.5.1      pbapply_1.4-0       stringr_1.4.0      
#[28] digest_0.6.25       R.utils_2.8.0       pkgconfig_2.0.3    
#[31] htmltools_0.5.1.1   bibtex_0.4.2        htmlwidgets_1.3    
#[34] rlang_1.0.6         zoo_1.8-5           jsonlite_1.6       
#[37] ica_1.0-2           gtools_3.8.1        R.oo_1.22.0        
#[40] magrittr_1.5        Matrix_1.2-17       Rcpp_1.0.4.6       
#[43] munsell_0.5.0       fansi_0.4.1         ape_5.3            
#[46] reticulate_1.11.1   R.methodsS3_1.7.1   stringi_1.4.3      
#[49] gbRd_0.4-11         MASS_7.3-51.3       gplots_3.0.1.1     
#[52] Rtsne_0.15          plyr_1.8.4          grid_3.5.1         
#[55] parallel_3.5.1      gdata_2.18.0        listenv_0.7.0      
#[58] ggrepel_0.8.1       crayon_1.3.4        lattice_0.20-41    
#[61] cowplot_0.9.4       splines_3.5.1       SDMTools_1.1-221   
#[64] pillar_1.4.3        igraph_1.2.4        future.apply_1.3.0 
#[67] reshape2_1.4.3      codetools_0.2-16    leiden_0.3.1       
#[70] glue_1.4.0          lsei_1.2-0          metap_1.1          
#[73] data.table_1.12.2   RcppParallel_4.4.4  png_0.1-7          
#[76] Rdpack_0.10-1       gtable_0.3.0        RANN_2.6.1         
#[79] purrr_0.3.4         tidyr_0.8.3         future_1.15.1      
#[82] assertthat_0.2.1    rsvd_1.0.2          survival_2.44-1.1  
#[85] viridisLite_0.3.0   tibble_2.1.3        cluster_2.0.7-1    
#[88] globals_0.12.4      fitdistrplus_1.0-14 ROCR_1.0-7         
