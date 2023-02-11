###heatmap for Nodal targets
normalized_rld_nodal_Sel <- read.csv("path/normalized_rld_nodal_Sel.csv",row.names = 1)
data_scale <- as.data.frame(t(apply(normalized_rld_nodal_Sel,1,scale))) ##Z-score
names(data_scale) <- names(normalized_rld_nodal_Sel)
##corrplot(cor(t(data_scale)),order = "hclust")
library(ComplexHeatmap)
library(circlize)

samples <- c("wt_cyc_4h_0pg","wt_cyc_4h_0pg","wt_cyc_4h_0pg","wt_cyc_5h_0pg","wt_cyc_5h_0pg","wt_cyc_5h_0pg","wt_cyc_6h_0pg","wt_cyc_6h_0pg","wt_cyc_6h_0pg",
             "wt_cyc_4h_2pg","wt_cyc_4h_2pg","wt_cyc_4h_2pg","wt_cyc_5h_2pg","wt_cyc_5h_2pg","wt_cyc_5h_2pg","wt_cyc_6h_2pg","wt_cyc_6h_2pg","wt_cyc_6h_2pg",
             "wt_cyc_4h_6pg","wt_cyc_4h_6pg","wt_cyc_4h_6pg","wt_cyc_5h_6pg","wt_cyc_5h_6pg","wt_cyc_5h_6pg","wt_cyc_6h_6pg","wt_cyc_6h_6pg","wt_cyc_6h_6pg",
             "wt_cyc_4h_10pg","wt_cyc_4h_10pg","wt_cyc_4h_10pg","wt_cyc_5h_10pg","wt_cyc_5h_10pg","wt_cyc_5h_10pg","wt_cyc_6h_10pg","wt_cyc_6h_10pg","wt_cyc_6h_10pg")

my_color <- c("paleturquoise1","skyblue1","slateblue1","khaki1","darkgoldenrod1","goldenrod3",
              "lightpink","tomato","orangered4","maroon1","mediumorchid1","purple1")

ha = HeatmapAnnotation(samples = samples,col = list(samples = c("wt_cyc_4h_0pg"="paleturquoise1",
                                                                "wt_cyc_5h_0pg"="skyblue1",
                                                                "wt_cyc_6h_0pg"="slateblue1",
                                                                "wt_cyc_4h_2pg"="khaki1",
                                                                "wt_cyc_5h_2pg"="darkgoldenrod1",
                                                                "wt_cyc_6h_2pg"="goldenrod3",
                                                                "wt_cyc_4h_6pg"="lavenderblush",
                                                                "wt_cyc_5h_6pg"="lavender",
                                                                "wt_cyc_6h_6pg"="lavenderblush4",
                                                                "wt_cyc_4h_10pg"="maroon1",
                                                                "wt_cyc_5h_10pg"="mediumorchid1",
                                                                "wt_cyc_6h_10pg"="purple1")))

set.seed(18)
km.res_scale <- kmeans(data_scale[1:36],4)
clusters <- as.character(km.res_scale$cluster)
genes <- names(km.res_scale$cluster)
cluste1_genes <- genes[which(clusters == "1")]
cluste2_genes <- genes[which(clusters == "2")]
cluste3_genes <- genes[which(clusters == "3")]
cluste4_genes <- genes[which(clusters == "4")]

concentration <- c("cyc_0pg","cyc_0pg","cyc_0pg","cyc_0pg","cyc_0pg","cyc_0pg","cyc_0pg","cyc_0pg","cyc_0pg",
                   "cyc_2pg","cyc_2pg","cyc_2pg","cyc_2pg","cyc_2pg","cyc_2pg","cyc_2pg","cyc_2pg","cyc_2pg",
                   "cyc_6pg","cyc_6pg","cyc_6pg","cyc_6pg","cyc_6pg","cyc_6pg","cyc_6pg","cyc_6pg","cyc_6pg",
                   "cyc_10pg","cyc_10pg","cyc_10pg","cyc_10pg","cyc_10pg","cyc_10pg","cyc_10pg","cyc_10pg","cyc_10pg")

my_color <- c("darkslateblue","darkgoldenrod","lavenderblush4","darkorchid")

top_ha = HeatmapAnnotation(concentration = concentration,col = list(concentration = c("cyc_0pg"="darkslateblue",
                                                                                      "cyc_2pg"="darkgoldenrod",
                                                                                      "cyc_6pg"="lavenderblush4",
                                                                                      "cyc_10pg"="darkorchid")))

set.seed(18)
p <- Heatmap(data_scale[1:36],name = "Z-score", 
             km = 4,
             column_names_side = "bottom",
             col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
             cluster_columns = FALSE,
             row_dend_side = "left",
             show_row_names = TRUE,
             show_column_names = FALSE,
             # km_title = "%i"
             top_annotation  = top_ha,
             bottom_annotation = ha,
             cluster_row_slices = FALSE
) ##+ Heatmap(data_scale$cluster, name = "cluster", width = unit(5, "mm"),col = c("palegreen","mediumpurple1","pink"),row_title = FALSE)


library(Seurat)
library(dplyr)
library(Matrix)

###single cell RNA-seq data analysis
## analysis all cells using seurat package
list.files("path/filtered_gene_bc_matrices/Danio_rerio.GRCz11")
nodal_explant <- Read10X(data.dir = "path//filtered_gene_bc_matrices/Danio_rerio.GRCz11")
dense.size <- object.size(x = as.matrix(x = nodal_explant))
nodal_explant_object <- CreateSeuratObject(raw.data = nodal_explant, min.cells = 3, min.genes = 100, project = "nodal_explant")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = nodal_explant_object@data), value = TRUE)
percent.mito <- Matrix::colSums(nodal_explant_object@raw.data[mito.genes, ]) / Matrix::colSums(nodal_explant_object@raw.data)
nodal_explant_object <- AddMetaData(object = nodal_explant_object, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = nodal_explant_object, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = nodal_explant_object, gene1 = "nUMI", gene2 = "nGene")
GenePlot(object = nodal_explant_object, gene1 = "nUMI", gene2 = "percent.mito")
nodal_explant_object <- FilterCells(object = nodal_explant_object, subset.names = c("nGene", "percent.mito"),     low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.05))
summary(nodal_explant_object@raw.data[,1])
nodal_explant_object <- NormalizeData(object = nodal_explant_object, normalization.method = "LogNormalize", scale.factor = 10000)
summary(nodal_explant_object@data[,1])
nodal_explant_object <- FindVariableGenes(object = nodal_explant_object, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = nodal_explant_object@var.genes)
nodal_explant_object <- ScaleData(object = nodal_explant_object, vars.to.regress = c("nUMI", "percent.mito"))
summary(nodal_explant_object@scale.data[,1])
nodal_explant_object <- RunPCA(object = nodal_explant_object, pc.genes = nodal_explant_object@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = nodal_explant_object, pcs.use = 1:2)
PCAPlot(object = nodal_explant_object, dim.1 = 1, dim.2 = 5)
nodal_explant_object <- ProjectPCA(object = nodal_explant_object, do.print = FALSE)
PCHeatmap(object = nodal_explant_object, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = nodal_explant_object, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
nodal_explant_object <- JackStraw(object = nodal_explant_object, num.replicate = 100,do.par = TRUE,num.cores = 12)
JackStrawPlot(object = nodal_explant_object, PCs = 1:12)
PCElbowPlot(object = nodal_explant_object)
nodal_explant_object_0.6 <- FindClusters(object = nodal_explant_object, reduction.type = "pca", dims.use = 1:12, resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = nodal_explant_object_0.6)
nodal_explant_object_0.6 <- RunTSNE(object = nodal_explant_object_0.6, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = nodal_explant_object_0.6)
nodal_explant_0.6.markers <- FindAllMarkers(object = nodal_explant_object_0.6, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
##nodal_explant_0.6.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top2 <- nodal_explant_0.6.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
VlnPlot(object = nodal_explant_object_0.6, features.plot = top2$gene, use.raw = TRUE, y.log = TRUE)
FeaturePlot(object = nodal_explant_object_0.6, features.plot = top2$gene, cols.use = c("grey", "red"), reduction.use = "tsne")
##nodal_explant_3d <- RunTSNE(object = nodal_explant_object, dims.use = 1:10, do.fast = TRUE,dim.embed = 3)
DoHeatmap(object = nodal_explant_object_0.6, genes.use = top20$gene, slim.col.label = TRUE, remove.key = TRUE,col.low = "darkgreen",col.mid = "lightblue",col.high = "red")
saveRDS(nodal_explant_object_0.6, file = "nodal_explant_dim_10_resolution_0.6.rds")

###construction of URD tree
library(URD)
library("rgl")
wt_cyc_count <- read.table("wt_cyc_combined_Raw_UMI_Matrix.tsv",sep="\t",header = TRUE)
wt_cyc_meta <- read.table("wt_cyc_combined_Metadata.tsv",sep="\t",header = TRUE)
wt_cyc_matrix <- as.matrix(wt_cyc_count)
wt_cyc_URD <- createURD(count.data = wt_cyc_matrix, meta = wt_cyc_meta, min.cells=3, min.counts=3)

wt_cyc_URD@group.ids$stage <- as.character(wt_cyc_URD@meta[rownames(wt_cyc_URD@group.ids),"timepoint"])

set_stage_order <- function(URD_object){
  URD_object@meta$stage.use <- NA
  time_seq <- c("A_4h","B_5h","C_6h","D_8h","E_10h")
  for (time in time_seq) {
    URD_object@meta[which(URD_object@meta$timepoint == strsplit(time,split = "_")[[1]][2]),]$stage.use <- rep(time,length(which(URD_object@meta$timepoint == strsplit(time,split = "_")[[1]][2])))
  }
  return(URD_object)
}

wt_cyc_URD <- set_stage_order(wt_cyc_URD)
wt_cyc_URD@group.ids$stage <- as.character(wt_cyc_URD@meta[rownames(wt_cyc_URD@group.ids),"stage.use"])

stages <- sort(unique(wt_cyc_URD@group.ids$stage))
var.by.stage <- lapply(seq(2,5,1), function(n) {
  findVariableGenes(wt_cyc_URD, cells.fit=cellsInCluster(wt_cyc_URD, "stage", stages[(n-1):n]),
                    set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100,
                    main.use=paste0("Stages ", stages[n-1], " to ", stages[n]), do.plot=T)
})

var.genes <- sort(unique(unlist(var.by.stage)))
wt_cyc_URD@var.genes <- var.genes
##Calculate PCA and tSNE
wt_cyc_URD <- calcPCA(wt_cyc_URD, mp.factor = 2,verbose = T)
pcSDPlot(wt_cyc_URD)
set.seed(19)
wt_cyc_URD <- calcTsne(object = wt_cyc_URD,verbose = T)
plotDim(wt_cyc_URD, "stage.use", plot.title = "tSNE: Stage")

wt_cyc_URD <- graphClustering(wt_cyc_URD, dim.use="pca", num.nn=c(15,20,30), do.jaccard=T, method="Louvain")
# plotDim(wt_cyc_URD, "stage",  legend=T, plot.title="Developmental Stage", alpha=0.5)
plotDim(wt_cyc_URD, "Louvain-15", legend=T, plot.title="Louvain-Jaccard Graph-based Clustering (15 NNs)", alpha=1,label.clusters = TRUE)
wt_cyc_URD <- calcKNN(wt_cyc_URD, nn=100)
outliers <- knnOutliers(wt_cyc_URD, nn.1=1, nn.2=20, x.max=40, slope.r=1.1, int.r=2.9, slope.b=0.85, int.b=10, title = "Identifying Outliers by k-NN Distance.")

###identify apoptotic like cells
gridExtra::grid.arrange(grobs=list(
  # Plot some apoptotic-like markers
  plotDim(wt_cyc_URD, "isg15", alpha=0.4, point.size=0.5),
  plotDim(wt_cyc_URD, "foxo3b", alpha=0.4, point.size=0.5),
  plotDim(wt_cyc_URD, "mdm2", alpha=0.4, point.size=0.5),
  # Figure out which cluster corresponds to these cells
  plotDimHighlight(wt_cyc_URD, clustering="Louvain-15", cluster="15", legend=F)
))

apoptotic.like.cells <- cellsInCluster(wt_cyc_URD, "Louvain-15", "15")
cells.keep <- setdiff(colnames(wt_cyc_URD@logupx.data), c(outliers))
wt_cyc_URD_keep <- urdSubset(wt_cyc_URD, cells.keep=cells.keep)
saveRDS(wt_cyc_URD, file="path/wt_cyc_URD_new.rds")
saveRDS(wt_cyc_URD_keep, file="path/wt_cyc_URD_keep.rds")

##need run different parameters for this section
wt_cyc_URD_keep <- calcDM(wt_cyc_URD_keep, knn=100, sigma.use=16,verbose = T)
saveRDS(wt_cyc_URD_keep, file="wt_cyc_URD_keep_dm16.rds")

plotDimArray(object = wt_cyc_URD_keep, reduction.use = "dm", dims.to.plot = 1:18, label="stage", plot.title="", outer.title="STAGE - Diffusion Map Sigma 16", legend=F, alpha=0.3)

root.cells <- rownames(wt_cyc_URD_keep@meta)[wt_cyc_URD_keep@meta$stage.use=="A_4h"]
flood.result <- floodPseudotime(wt_cyc_URD_keep, root.cells=root.cells, n=50, minimum.cells.flooded=2, verbose=T)
saveRDS(flood.result,"wt_cyc_URD_keep_flood_dm16.rds")
wt_cyc_URD_keep <- floodPseudotimeProcess(wt_cyc_URD_keep, flood.result, floods.name="pseudotime", max.frac.NA=0.4, pseudotime.fun=mean, stability.div=20)
pseudotimePlotStabilityOverall(wt_cyc_URD_keep)
plotDim(wt_cyc_URD_keep, "pseudotime", plot.title = "Pseudotime")
plotDists(wt_cyc_URD_keep, "pseudotime", "stage", plot.title="Pseudotime by stage")
saveRDS(wt_cyc_URD_keep,"wt_cyc_URD_keep_dim16_pseudotime.rds")

###determing tips 

wt_cyc_URD_keep_10hpf <- urdSubset(wt_cyc_URD_keep, cells.keep=cellsInCluster(wt_cyc_URD_keep, "stage", "E_10h"))
wt_cyc_URD_keep_10hpf@var.genes <- var.by.stage[[4]]
wt_cyc_URD_keep_10hpf <- calcPCA(wt_cyc_URD_keep_10hpf, mp.factor = 1.5)
pcSDPlot(wt_cyc_URD_keep_10hpf)
set.seed(18)
wt_cyc_URD_keep_10hpf <- calcTsne(wt_cyc_URD_keep_10hpf,dim.use="pca")
wt_cyc_URD_keep_10hpf <- graphClustering(wt_cyc_URD_keep_10hpf, num.nn =c(5,8,10,15), do.jaccard=T, method="Louvain")
wt_cyc_URD_keep_10hpf <- graphClustering(wt_cyc_URD_keep_10hpf, num.nn = c(10,15,20,30,40), method="Infomap", do.jaccard = T)

clusterings <- c(paste0("Louvain-", c(5,8,10,15)), paste0("Infomap-", c(10,15,20,30,40)))
for (c in clusterings) {
  plot(plotDim(wt_cyc_URD_keep_10hpf, c, legend=T))
}
plotDim(wt_cyc_URD_keep_10hpf, "sox17", point.size=3)
plotDim(wt_cyc_URD_keep_10hpf, "Louvain-15", plot.title = "Louvain (15 NN) graph clustering", point.size=1,label.clusters = TRUE)
clusters <- unique(wt_cyc_URD_keep_10hpf@group.ids$`Louvain-15`)
pr.markers <- lapply(clusters, function(c) markersAUCPR(wt_cyc_URD_keep_10hpf, clust.1 = c, clustering="Louvain-15", genes.use=wt_cyc_URD_keep_10hpf@var.genes))
names(pr.markers) <- clusters

###A totally independent clustering could be attempted here, but since there is so much prior knowledge in zebrafish, we use it to evaluate and annotate our clustering.
##First, we created data.frames to keep track of our cluster assignments.
Louv15.n <- length(unique(wt_cyc_URD_keep_10hpf@group.ids$`Louvain-15`))
Louv15.cluster.assignments <- data.frame(
  cluster=1:Louv15.n,
  name=rep(NA, Louv15.n),
  tip=rep(NA, Louv15.n),
  row.names=1:Louv15.n
)

##assign clusters
##endoderm 
plotDot(wt_cyc_URD_keep_10hpf,genes = c("prdx5","foxa1","foxa2","nkx2.7","cdx4","sox17","irx7"),clustering = "Louvain-15")
Louv15.cluster.assignments["24","name"] <- "Endoderm_Pharyngeal"
### axial mesoderm 
plotDot(wt_cyc_URD_keep_10hpf,genes = c("ctslb","p4ha1b","col2a1a","col8a1a","tbxta","noto","ntd5","gsc","frzb","he1.1","he1.2"),clustering = "Louvain-15")

Louv15.cluster.assignments["24","name"] <- "Endoderm_Pharyngeal"
Louv15.cluster.assignments["17","name"] <- "Prechordal_Plate_A"
Louv15.cluster.assignments["15","name"] <- "Paraxial_Mesoderm_A" ##cephalic paraxial mesoderm
Louv15.cluster.assignments["4","name"] <- "Notochord_A"
Louv15.cluster.assignments["8","name"] <- "Paraxial_Mesoderm_P"

### intermediate/lateral mesoderm absent lateral mesoderm
plotDot(wt_cyc_URD_keep_10hpf,genes = c("fsta","foxf2a","tal1","lmo2","morc3b","spi1b","pax2a","foxj1a","hand2"),clustering = "Louvain-15")
plotDim(wt_cyc_URD_keep_10hpf, "lmo2", point.size=3)

###paraxial mesoderm
plotDot(wt_cyc_URD_keep_10hpf,genes = c("pcdh8","myl10","myod1","myf5","meox1","ripply1","aldh1a2","mespba","tbx6","tbx6l","wnt8a","fgf8a","noto","tbxta","sox2"),clustering = "Louvain-15")
Louv15.cluster.assignments["18","name"] <- "Tailbud" 
Louv15.cluster.assignments["11","name"] <- "Tailbud" 
Louv15.cluster.assignments["20","name"] <- "Axis_A" ## axis and adxial mesoderm anteriror
Louv15.cluster.assignments["16","name"] <- "Axis_P" ## axis and adxial mesoderm posterior

####Neural 
plotDot(wt_cyc_URD_keep_10hpf,genes = c("sox10","sox9b","foxd3","foxj1a","shha","shhb","aldh1a2","foxa2"),clustering = "Louvain-15")
wt_cyc_URD_keep_10hpf@meta$Louvain_15 <- wt_cyc_URD_keep_10hpf@group.ids$`Louvain-15`
write.csv(wt_cyc_URD_keep_10hpf@meta,"wt_cyc_URD_keep_10hpf_meta.csv")


plotDim(wt_cyc_URD_keep_10hpf, "h1m", alpha=0.4, point.size=0.5)
plotDim(wt_cyc_URD_keep_10hpf, "nanog", alpha=0.4, point.size=0.5)
plotDim(wt_cyc_URD_keep_10hpf, "sox17", alpha=0.4, point.size=0.5)
plotDim(wt_cyc_URD_keep_10hpf, "ctslb", alpha=0.4, point.size=0.5)
plotDim(wt_cyc_URD_keep_10hpf, "krt18", alpha=0.4, point.size=0.5)

pp_cell <- rownames(wt_cyc_URD_keep_10hpf@meta[which(wt_cyc_URD_keep_10hpf@meta$Louvain_15=="17"),])
plotDimHighlight(wt_cyc_URD_keep_10hpf,clustering="Louvain-15",cluster="17", legend=F)

###annotating Louv15 clusters
Louv15.cluster.assignments["1","name"] <- "Forebrain" ##Telencephalon
Louv15.cluster.assignments["2","name"] <- "Forebrain"
Louv15.cluster.assignments["3","name"] <- "Midbrain" 
Louv15.cluster.assignments["4","name"] <- "Notochord"
Louv15.cluster.assignments["5","name"] <- "Spinal Cord"
Louv15.cluster.assignments["6","name"] <- "Epidermal" ## Anterior Neural Ridge dlx3b
Louv15.cluster.assignments["7","name"] <- "Forebrain" ##Telencephalon 
Louv15.cluster.assignments["8","name"] <- "Notochord" 
Louv15.cluster.assignments["9","name"] <- "Forebrain" ##rx3
Louv15.cluster.assignments["10","name"] <- "Hindbrain" ## dbx1a
Louv15.cluster.assignments["11","name"] <- "Tailbud" 
Louv15.cluster.assignments["12","name"] <- "Neural" ##not specific 
Louv15.cluster.assignments["13","name"] <- "Apoptotic like" ##
Louv15.cluster.assignments["14","name"] <- "Apoptotic like"
Louv15.cluster.assignments["15","name"] <- "Adaxial mesoderm"
Louv15.cluster.assignments["16","name"] <- "Notochord and Prechordal Plate hybrid"
Louv15.cluster.assignments["17","name"] <- "Prechordal Plate"
Louv15.cluster.assignments["18","name"] <- "Tailbud" 
Louv15.cluster.assignments["19","name"] <- "EVL"
Louv15.cluster.assignments["20","name"] <- "Axis Mesoderm"
Louv15.cluster.assignments["21","name"] <- "Forebrain" ## her6,her4
Louv15.cluster.assignments["22","name"] <- "EVL" 
Louv15.cluster.assignments["23","name"] <- "PGC" ##h1m, nanog
Louv15.cluster.assignments["24","name"] <- "Endoderm"
Louv15.cluster.assignments$clustering <- "Louvain-15"
wt_cyc_URD_keep_10hpf@group.ids$clusters.URD_10hpf.name <- NA
wt_cyc_URD_keep_10hpf@group.ids$clusters.URD_10hpf.num <- NA
# Copy cell identities over for each cluster
for (i in 1:nrow(Louv15.cluster.assignments)) {
  cells <- cellsInCluster(wt_cyc_URD_keep_10hpf, clustering = Louv15.cluster.assignments[i,"clustering"], cluster = Louv15.cluster.assignments[i,"cluster"])
  wt_cyc_URD_keep_10hpf@group.ids[cells,"clusters.URD_10hpf.name"] <- Louv15.cluster.assignments[i,"name"]
  wt_cyc_URD_keep_10hpf@group.ids[cells,"clusters.URD_10hpf.num"] <- as.character(Louv15.cluster.assignments[i,"cluster"])
}
plotDim(wt_cyc_URD_keep_10hpf, "clusters.URD_10hpf.name", point.size=1)

###save markers 

save_marker_list <- function(marker_list){
  for (i in names(marker_list)) {
    write.csv(marker_list[[i]],paste("cluster",i,".csv",sep = ""))
  }
}

save_marker_list(pr.markers)

# Transfer clusterings to main object
## Need to transfer cluster identities from the URD_10hpf only object to the full object.
wt_cyc_URD_keep@group.ids$`URD_10hpf-Louvain-15` <- NA
wt_cyc_URD_keep@group.ids[rownames(wt_cyc_URD_keep_10hpf@group.ids), "URD_10hpf-Louvain-15"] <- wt_cyc_URD_keep_10hpf@group.ids$`Louvain-15`
wt_cyc_URD_keep@group.ids$`URD_10hpf-Cluster` <- NA
wt_cyc_URD_keep@group.ids[rownames(wt_cyc_URD_keep_10hpf@group.ids), "URD_10hpf-Cluster"] <- wt_cyc_URD_keep_10hpf@group.ids$`clusters.URD_10hpf.name`
wt_cyc_URD_keep@group.ids$`URD_10hpf-Cluster-Num` <- NA
wt_cyc_URD_keep@group.ids[rownames(wt_cyc_URD_keep_10hpf@group.ids), "URD_10hpf-Cluster-Num"] <- wt_cyc_URD_keep_10hpf@group.ids$`clusters.URD_10hpf.num`
saveRDS(wt_cyc_URD_keep, file="wt_cyc_URD_keep_annotation.rds")
saveRDS(wt_cyc_URD_keep_10hpf, file="wt_cyc_URD_keep_10hpf_annotation.rds")
Louv15.cluster.assignments$tip <- TRUE
Louv15.cluster.assignments[Louv15.cluster.assignments$cluster == "12",]$tip <- FALSE
Louv15.cluster.assignments[Louv15.cluster.assignments$cluster == "13",]$tip <- FALSE
Louv15.cluster.assignments[Louv15.cluster.assignments$cluster == "14",]$tip <- FALSE
Louv15.cluster.assignments[Louv15.cluster.assignments$cluster == "16",]$tip <- FALSE
Louv15.cluster.assignments[Louv15.cluster.assignments$cluster == "20",]$tip <- FALSE
write.csv(Louv15.cluster.assignments, file="Louv15_tips_use.csv")

# Plot tips in diffusion map
# Not all clusters from the final stage should really comprise tips of the developmental tree
wt_cyc_URD_keep@group.ids$pop <- NA
wt_cyc_URD_keep@group.ids[cellsInCluster(wt_cyc_URD_keep, "URD_10hpf-Cluster-Num", "12"), "pop"] <- "1"
wt_cyc_URD_keep@group.ids[cellsInCluster(wt_cyc_URD_keep, "URD_10hpf-Cluster-Num", "13"), "pop"] <- "1"
wt_cyc_URD_keep@group.ids[cellsInCluster(wt_cyc_URD_keep, "URD_10hpf-Cluster-Num", "14"), "pop"] <- "1"
wt_cyc_URD_keep@group.ids[cellsInCluster(wt_cyc_URD_keep, "URD_10hpf-Cluster-Num", "16"), "pop"] <- "1"
wt_cyc_URD_keep@group.ids[cellsInCluster(wt_cyc_URD_keep, "URD_10hpf-Cluster-Num", "20"), "pop"] <- "1"
wt_cyc_URD_keep@group.ids$clusters.URD_10hpf.name <- wt_cyc_URD_keep@group.ids$`URD_10hpf-Cluster`
wt_cyc_URD_keep@group.ids$URD_10hpf_Cluster.Num <- wt_cyc_URD_keep@group.ids$`URD_10hpf-Cluster-Num`

pop_cells <- rownames(wt_cyc_URD_keep@group.ids[which(wt_cyc_URD_keep@group.ids$pop == "1"),])
wt_cyc_URD_keep@group.ids[pop_cells,]$clusters.URD_10hpf.name <- NA
wt_cyc_URD_keep@group.ids[pop_cells,]$URD_10hpf_Cluster.Num <- NA

wt_cyc_URD_keep <- loadTipCells(wt_cyc_URD_keep, tips="URD_10hpf_Cluster.Num")

library("knitr")
diffusion.logistic <- pseudotimeDetermineLogistic(wt_cyc_URD_keep, "pseudotime", optimal.cells.forward=20, max.cells.back=40, pseudotime.direction="<", do.plot=T, print.values=T)
biased.tm <- pseudotimeWeightTransitionMatrix(wt_cyc_URD_keep, pseudotime = "pseudotime", logistic.params = diffusion.logistic, pseudotime.direction = "<")
root.cells <- rownames(wt_cyc_URD_keep@meta)[wt_cyc_URD_keep@meta$stage.use=="A_4h"]
wt_cyc_URD_keep@group.ids[rownames(wt_cyc_URD_keep@group.ids), "tip.clusters"] <- wt_cyc_URD_keep@group.ids$URD_10hpf_Cluster.Num
saveRDS(wt_cyc_URD_keep,"wt_cyc_URD_keep_before_walks.rds")
wt_cyc_URD_keep_tip10000.walks <- simulateRandomWalksFromTips(wt_cyc_URD_keep, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = biased.tm, n.per.tip = 10000, root.visits = 1, max.steps = 5000, verbose = T)
saveRDS(wt_cyc_URD_keep_tip10000.walks,"wt_cyc_URD_keep_tip10000.walks.rds")
wt_cyc_URD_keep_tip10000 <- processRandomWalksFromTips(wt_cyc_URD_keep, wt_cyc_URD_keep_tip10000.walks, verbose = T)
plotDim(wt_cyc_URD_keep_tip10000, "tip.clusters", plot.title="Cells in each tip")

#ectoderm_neural_anterior forebrain
wt_cyc_URD_keep_tip10000 <- combineTipVisitation(wt_cyc_URD_keep_tip10000, "1", "2", "1")
wt_cyc_URD_keep_tip10000 <- combineTipVisitation(wt_cyc_URD_keep_tip10000, "1", "7", "1")
wt_cyc_URD_keep_tip10000 <- combineTipVisitation(wt_cyc_URD_keep_tip10000, "1", "9", "1")
wt_cyc_URD_keep_tip10000 <- combineTipVisitation(wt_cyc_URD_keep_tip10000, "1", "21", "1")
#notochord
wt_cyc_URD_keep_tip10000 <- combineTipVisitation(wt_cyc_URD_keep_tip10000, "4", "8", "4")
#tailbud
wt_cyc_URD_keep_tip10000 <- combineTipVisitation(wt_cyc_URD_keep_tip10000, "11", "18", "18")
#EVL
wt_cyc_URD_keep_tip10000 <- combineTipVisitation(wt_cyc_URD_keep_tip10000, "19", "22", "22")

##finally 19 clusters are saved as tips
### 1,2,7,9,21 anterior neural ectoderm/forebrain
### 3 middle neural ectoderm/midbrain
### 10 posterior neural ectoderm/hindbrain
### 5 spinal cord
### 11,18 tailbud
### 17 prechordal plate
### 4,8 notochord
### 15 adxial mesoderm
### 24 endoderm
### 19,22 EVL
### 23 PGC
### 6 epidermal
wt_cyc_URD_keep_tip10000.tree <- buildTree(wt_cyc_URD_keep_tip10000, pseudotime = "pseudotime", verbose = T,tips.use = c("1","3","4","5","6","10","15","17","22","23","24"),
                                           save.all.breakpoint.info = T,divergence.method = "preference", cells.per.pseudotime.bin = 80, bins.per.pseudotime.window = 5,p.thresh=0.001)

tip.names <- unique(wt_cyc_URD_keep_tip10000.tree@group.ids[,c("clusters.URD_10hpf.name", "URD_10hpf_Cluster.Num")])
tip.names <- tip.names[complete.cases(tip.names),]
wt_cyc_URD_keep_tip10000.tree <- nameSegments(wt_cyc_URD_keep_tip10000.tree, segments=tip.names$URD_10hpf_Cluster.Num, segment.names=tip.names$clusters.URD_10hpf.name)
plotTree(wt_cyc_URD_keep_tip10000.tree, "stage.use", title="Developmental Stage")
saveRDS(wt_cyc_URD_keep_tip10000.tree,"wt_cyc_URD_keep_tip10000_1.tree")


wt_cyc_URD_keep_tip10000.tree_2 <- buildTree(wt_cyc_URD_keep_tip10000, pseudotime = "pseudotime", verbose = T,tips.use = c("1","3","4","5","6","10","15","17","18","22","23","24"),
                                             save.all.breakpoint.info = T,divergence.method = "preference", cells.per.pseudotime.bin = 50, bins.per.pseudotime.window = 5,p.thresh=0.001)
tip.names <- unique(wt_cyc_URD_keep_tip10000.tree_2@group.ids[,c("clusters.URD_10hpf.name", "URD_10hpf_Cluster.Num")])
tip.names <- tip.names[complete.cases(tip.names),]
wt_cyc_URD_keep_tip10000.tree_2 <- nameSegments(wt_cyc_URD_keep_tip10000.tree_2, segments=tip.names$URD_10hpf_Cluster.Num, segment.names=tip.names$clusters.URD_10hpf.name)
plotTree(wt_cyc_URD_keep_tip10000.tree_2, "stage.use", title="Developmental Stage")
saveRDS(wt_cyc_URD_keep_tip10000.tree_2,"wt_cyc_URD_keep_tip10000_2.tree")








