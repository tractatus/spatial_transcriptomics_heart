load('/Users/danielfurth/Documents/GitHub/CarnegieAtlas/r/carnegie/CN73_C2.RData')

library(wholebrain)

load('./data/ST_heart_03042018.RData')
load('./data/heart.RData')
quartz(width=7, height=7)
plot.registration(regi)
scale.factor<-regi$transformationgrid$height/dim(regi$transformationgrid$mx)
for(i in 2:regi$atlas$numRegions){
  lines(regi$atlas$outlines[[i]]$xT*scale.factor, regi$atlas$outlines[[i]]$yT*scale.factor, col='orange')
}

quartz.save(file = 'Supplementary_FigureH.pdf', type='pdf')


quartz(width=7, height=7)
plot.registration(regi, draw.trans.grid = TRUE)

points(seg.spots$soma$x, seg.spots$soma$y, pch=21, col=rgb(1,0,0,0.5), cex=0.6)

points(heart$atlas$x[which(heart$atlas$image == '4_CN73_C2_HE')], heart$atlas$y[which(heart$atlas$image == '4_CN73_C2_HE')], bg= as.character(heart$atlas$color2[which(heart$atlas$image == '4_CN73_C2_HE')]), col='black', pch=21, cex=0.6 )


scale.factor<-regi$transformationgrid$height/dim(regi$transformationgrid$mx)
scale.factor<-scale.factor[1]
for(i in 2:regi$atlas$numRegions){
  lines(regi$atlas$outlines[[i]]$xT*scale.factor, regi$atlas$outlines[[i]]$yT*scale.factor, col='orange')
}

quartz.save(file = 'Supplementary_FigureI.pdf', type='pdf')



regi <- get.forward.warpRCPP(regi)
index <- round( cbind(seg.spots$soma$y, seg.spots$soma$x)/scale.factor)
somaX <- regi$transformationgrid$mxF[index]/scale.factor
somaY <- regi$transformationgrid$myF[index]/scale.factor

index2 <- round( cbind(heart$atlas$y[which(heart$atlas$image == '4_CN73_C2_HE')], heart$atlas$x[which(heart$atlas$image == '4_CN73_C2_HE')])/scale.factor)
somaX2 <- regi$transformationgrid$mxF[index2]/scale.factor
somaY2 <- regi$transformationgrid$myF[index2]/scale.factor

quartz()
plot(somaX, somaY, asp = 1, axes = F, xlab = "", ylab = "",
     col = 0, ylim= rev(c(0,598/2)), xlim=c(0,532/2))

numPaths <- regi$atlas$numRegions
outlines <- regi$atlas$outlines

lapply(2:numPaths, function(x) {
  polygon(outlines[[x]]$x/scale.factor, outlines[[x]]$y/scale.factor,
          border = gray(0.2), col = 'lightgray')
})


hei <- dim(regi$transformationgrid$mx)[1]
wid <- dim(regi$transformationgrid$mx)[2]
lapply(seq(1, hei, by = 50), function(x) {
  lines(regi$transformationgrid$mxF[x, ]/scale.factor,
        regi$transformationgrid$myF[x, ]/scale.factor,
        col = "lightblue")
})
lines(regi$transformationgrid$mxF[hei, ]/scale.factor,
      regi$transformationgrid$myF[hei, ]/scale.factor,
      col = "lightblue")
lapply(seq(1, wid, by = 50), function(x) {
  lines(regi$transformationgrid$mxF[, x]/scale.factor,
        regi$transformationgrid$myF[, x]/scale.factor,
        col = "lightblue")
})
lines(regi$transformationgrid$mxF[, wid]/scale.factor,
      regi$transformationgrid$myF[, wid]/scale.factor,
      col = "lightblue")


points(somaX, somaY, pch=16, col= 'black', cex=0.7)
points(somaX2, somaY2, pch=21, bg=as.character(heart$atlas$color2[which(heart$atlas$image == '4_CN73_C2_HE')]), cex=0.7)

quartz.save(file = 'Supplementary_FigureJ.pdf', type='pdf')

pp<-structure(list(FOV = 30, ignoreExtent = FALSE, listeners = 8L, 
                   mouseMode = structure(c("trackball", "zoom", "fov", "pull"
                   ), .Names = c("left", "right", "middle", "wheel")), skipRedraw = FALSE, 
                   userMatrix = structure(c(0.0179558768868446, 0.813902854919434, 
                                            0.580723106861115, 0, 0.878040909767151, -0.29065603017807, 
                                            0.380215406417847, 0, 0.478248953819275, 0.503071665763855, 
                                            -0.719858825206757, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), scale = c(1, 
                                                                                                             1, 1), viewport = structure(c(0L, 0L, 1280L, 720L), .Names = c("x", 
                                                                                                                                                                            "y", "width", "height")), zoom = 1.10249984264374, windowRect = c(2023L, 
                                                                                                                                                                                                                                              120L, 3303L, 840L), family = "sans", font = 1L, cex = 1, 
                   useFreeType = TRUE), .Names = c("FOV", "ignoreExtent", "listeners", 
                                                   "mouseMode", "skipRedraw", "userMatrix", "scale", "viewport", 
                                                   "zoom", "windowRect", "family", "font", "cex", "useFreeType"))

open3d(windowRect = c(0, 0, 1280, 720))
par3d(pp)
drawScene.rgl(organ[which(names(organ.dwnsmp)%in%c('WH'))])
radius.of.spots.in.atlas.pixels<- (100/(2383.36/532))/3
spheres3d(598-heart$atlas$rostral.caudal[which(heart$atlas$image == '4_CN73_C2_HE')], 532-heart$atlas$right.left[which(heart$atlas$image == '4_CN73_C2_HE')], heart$atlas$anterior.posterior[which(heart$atlas$image == '4_CN73_C2_HE')], col=heart$atlas$color2[which(heart$atlas$image == '4_CN73_C2_HE')], radius=radius.of.spots.in.atlas.pixels)
segments3d(c(0,0, 0,598,0,598,598,598), c(0,532,0,0,532,532,0,532), rep(c(250, 250), 4), lwd=3)

#bounding box
segments3d(x = c(0,0, 0, 0, 598, 598, 598, 598, 0, 0, 0, 0, 598, 598, 598, 598, 0, 598, 0, 598, 0, 598, 0, 598, 0, 598) , y = c(0,532, 0, 532, 0,532, 0, 532, 0, 0, 532, 532,  0, 0, 532, 532, 0, 0, 0, 0, 532, 532, 532, 532, 532, 532), z = c(0,0, 480, 480, 0,0, 480, 480, 0, 480, 0, 480,  0, 480, 0, 480, 0, 0, 480, 480, 0, 0, 480, 480))

rgl.snapshot('Supplementary_FigureK.png')


drawScene.rgl(organ[which(names(organ.dwnsmp)%in%c('WH'))])
radius.of.spots.in.atlas.pixels<- (100/(2383.36/532))/3
spheres3d(598-heart$atlas$rostral.caudal, 532-heart$atlas$right.left, heart$atlas$anterior.posterior, col=heart$atlas$color2, radius=radius.of.spots.in.atlas.pixels)
for(i in unique(heart$atlas$anterior.posterior))
segments3d(c(0, 0, 0,598,0,598,598,598), c(0,532,0,0,532,532,0,532), rep(i, 4*2), lwd=2)

#bounding box
segments3d(x = c(0,0, 0, 0, 598, 598, 598, 598, 0, 0, 0, 0, 598, 598, 598, 598, 0, 598, 0, 598, 0, 598, 0, 598, 0, 598) , y = c(0,532, 0, 532, 0,532, 0, 532, 0, 0, 532, 532,  0, 0, 532, 532, 0, 0, 0, 0, 532, 532, 532, 532, 532, 532), z = c(0,0, 480, 480, 0,0, 480, 480, 0, 480, 0, 480,  0, 480, 0, 480, 0, 0, 480, 480, 0, 0, 480, 480))

rgl.snapshot('Supplementary_FigureL.png')
