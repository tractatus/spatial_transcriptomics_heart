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

