load('./data/ST_heart_08032018.RData')

#Set all genes not detected oin other sections (they are NA) to 0 UMI count.
tmp<-heart$genes
tmp[which(is.na(tmp), arr.ind = TRUE)]<-0

#needs to be transposed for Seurat
tmp<-t(tmp)
#use indices instead of spot ID (e.g. 13x24).
colnames(tmp)<-seq_along(colnames(tmp))

library(Seurat)
#make Seurat object with min 1000 genes
seurat.object<-CreateSeuratObject(tmp,project="STheart", min.cells = 3, names.field = 2, min.genes = 1000)
#could add meta data here
seurat.object <- AddMetaData(object = seurat.object) #, metadata = images
#normalize
seurat.object <- NormalizeData(object = seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
#do cut off for diff expressed genes
seurat.object <- FindVariableGenes(object = seurat.object, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.2, x.high.cutoff = 5, y.cutoff = 0.5)
#scale here (can also regress out things)
seurat.object <- ScaleData(object = seurat.object) #, vars.to.regress=names(images)
#compute 100 PCAs
seurat.object <- RunPCA(object = seurat.object, pc.genes = seurat.object@var.genes, do.print = TRUE, pcs.compute = 100)
seurat.object <- JackStraw(object = seurat.object, num.replicate = 100, do.print = FALSE)
PCElbowPlot(object = seurat.object, num.pc=100)

#the first 12 CA contains most info
seurat.object <- FindClusters(object = seurat.object, reduction.type = "pca", dims.use = 1:12, resolution = 1, print.output = TRUE, save.SNN = TRUE)
#run the TSNE
seurat.object <- RunTSNE(object = seurat.object, dims.use = 1:12, do.fast = TRUE, check_duplicates = FALSE, max_iter=1000)

#plot the tSNE

par(mfrow=c(1,2))
plot(seurat.object@dr$tsne@cell.embeddings, pch=16, cex=0.9, col=paste0(heart$atlas$color2[as.integer(colnames(seurat.object@scale.data))], '80'), asp=1, xlab='', ylab='', main='', axes=FALSE)
legend.data<-data.frame(col = unique(heart$atlas$color2), acronym = na.omit(unique(heart$atlas$acronym)), name = na.omit(unique(heart$atlas$name)) )
legend.data<-legend.data[order(legend.data[,2]),]

par(mar=c(4,0,4,0))
plot(0,0, xlim=c(0,5), ylim=c(-5, nrow(legend.data)+5), col=0, axes=FALSE)
points(rep(0, nrow(legend.data)), 1:nrow(legend.data), pch=21, bg=as.character(legend.data$col))
text(rep(0, nrow(legend.data)), 1:nrow(legend.data), legend.data$acronym, pos=4)
text(rep(1, nrow(legend.data)), 1:nrow(legend.data), legend.data$name, pos=4)



########################
#
########################
#load 3D heart
load('./data/heart.RData')

library(rgl)
#open 3D plot window
open3d(windowRect = c(0, 0, 1280, 720))
#draw low-resolution heart with color coding
drawScene.rgl(organ.dwnsmp[which(names(organ.dwnsmp)%in%c('WH',  'RA', 'RV', 'LA', 'LV', 'P', 'A', 'OT'))])

#draw heart without color coding
drawScene.rgl(organ.dwnsmp[which(names(organ.dwnsmp)%in%c('WH'))])

#draw high resolution of heart
drawScene.rgl(organ[which(names(organ.dwnsmp)%in%c('WH'))])
#bounding box
box3d()

#draw all the spots with region color
drawScene.rgl(organ.dwnsmp[which(names(organ.dwnsmp)%in%c('WH'))])
spheres3d(598-heart$atlas$rostral.caudal, heart$atlas$right.left, heart$atlas$anterior.posterior, col=heart$atlas$color2, radius=5)

#save 3D image
rgl.snapshot(file='test3dheart.png')

#plot gene of interest OGN
gene.of.interest<-heart$genes[,which(colnames(heart$genes)=='ENSG00000106809')]

#color ramp palette base don gene expression
PaletteFunction <- colorRampPalette(c("blue", "white", "red"), space = "Lab")
gene.expression<-PaletteFunction(100)[as.numeric(cut(scale(log2(gene.of.interest+1)), breaks = 100))]

drawScene.rgl(organ.dwnsmp[which(names(organ.dwnsmp)%in%c('WH'))])
spheres3d(598-master.dataset$rostral.caudal, master.dataset$right.left, master.dataset$anterior.posterior, col=gene.expression, radius=5)
rgl.snapshot(filename='OGN_3d_heart.png')

