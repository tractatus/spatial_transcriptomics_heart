load('./data/ST_heart_13032018.RData')
#
mdata<-read.table('data/w6NoCluster_NrOfClusters_3_nrOfVariableGenes.tsv', sep='\t', header=TRUE)
#
meta.data<-read.table('data/Sample-info.tab.txt', sep='\t', header=TRUE)

################################################
#   Create index for clustering data set
################################################

#spot unique week by z section index
meta.dataindex1<-paste(meta.data$week, meta.data$z, sep='_')

#spot position index in meta data
mindex1<-paste(mdata$week, mdata$z, sep='_')
#get unique image name from meta data into cluster dataset
mdata$image<-substr(as.character(meta.data$Sample[match(mindex1, meta.dataindex1)]), nchar('FH6_1000L2_1'), nchar('FH5_1000L3_CN20_C1')  )
#spot position index in cluster data set
mindex2<-paste(mdata$x, mdata$y, sep='x')

#make a composite index (position and image) for clustering dataset
mdata$index<-paste(mindex2, mdata$image, sep = '_')

################################################
#   Create index Atlas data set
################################################

#make a composite index (position and image) for atlas dataset
heart$atlas$image2<-substr(heart$atlas$image, 3, 9)
heart$atlas$index<-paste(heart$atlas$spot.pos, heart$atlas$image2, sep = '_')

################
#   Final step
################
#get the color to your dataset
mdata$color<-as.character(heart$atlas$color2[match(mdata$index, heart$atlas$index)])

############
# PLOT IT
#############

#plot the tSNE
par(mfrow=c(1,2))
plot(mdata$V1, mdata$V2, pch=16, cex=0.9, col=paste0(mdata$color), asp=1, xlab='', ylab='', main='', axes=FALSE)


#plot legend is seperate plot
legend.data<-data.frame(col = unique(heart$atlas$color2), acronym = na.omit(unique(heart$atlas$acronym)), name = na.omit(unique(heart$atlas$name)) )
legend.data<-legend.data[order(legend.data[,2]),]
par(mar=c(4,0,4,0))
plot(0,0, xlim=c(0,5), ylim=c(-5, nrow(legend.data)+5), col=0, axes=FALSE, xlab='', yab='')
points(rep(0, nrow(legend.data)), 1:nrow(legend.data), pch=21, bg=as.character(legend.data$col))
text(rep(0, nrow(legend.data)), 1:nrow(legend.data), legend.data$acronym, pos=4)
text(rep(1, nrow(legend.data)), 1:nrow(legend.data), legend.data$name, pos=4)
