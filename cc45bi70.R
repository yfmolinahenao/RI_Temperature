setwd("~/Desktop/YFMH/BIOMOD2/")
library(biomod2)
library(raster)
library(ggplot2) 

cc45bi70_bio1 = raster("cc45bi70/cc45bi70/bio_1.tif")
bio1 = crop(cc45bi70_bio1, extent(-12,41.9,34.6,71.4))

cc45bi70_bio2 = raster("cc45bi70/cc45bi70/bio_2.tif")
bio2 = crop(cc45bi70_bio2, extent(-12,41.9,34.6,71.4))

cc45bi70_bio3 = raster("cc45bi70/cc45bi70/bio_3.tif")
bio3 = crop(cc45bi70_bio3, extent(-12,41.9,34.6,71.4))

cc45bi70_bio4 = raster("cc45bi70/cc45bi70/bio_4.tif")
bio4 = crop(cc45bi70_bio4, extent(-12,41.9,34.6,71.4))

cc45bi70_bio5 = raster("cc45bi70/cc45bi70/bio_5.tif")
bio5 = crop(cc45bi70_bio5, extent(-12,41.9,34.6,71.4))

cc45bi70_bio6 = raster("cc45bi70/cc45bi70/bio_6.tif")
bio6 = crop(cc45bi70_bio6, extent(-12,41.9,34.6,71.4))

cc45bi70_bio7 = raster("cc45bi70/cc45bi70/bio_7.tif")
bio7 = crop(cc45bi70_bio7, extent(-12,41.9,34.6,71.4))

cc45bi70_bio8 = raster("cc45bi70/cc45bi70/bio_8.tif")
bio8 = crop(cc45bi70_bio8, extent(-12,41.9,34.6,71.4))

cc45bi70_bio9 = raster("cc45bi70/cc45bi70/bio_9.tif")
bio9 = crop(cc45bi70_bio9, extent(-12,41.9,34.6,71.4))

cc45bi70_bio10 = raster("cc45bi70/cc45bi70/bio_10.tif")
bio10 = crop(cc45bi70_bio10, extent(-12,41.9,34.6,71.4))

cc45bi70_bio11 = raster("cc45bi70/cc45bi70/bio_11.tif")
bio11 = crop(cc45bi70_bio11, extent(-12,41.9,34.6,71.4))

cc45bi70_bio12 = raster("cc45bi70/cc45bi70/bio_12.tif")
bio12 = crop(cc45bi70_bio12, extent(-12,41.9,34.6,71.4))

cc45bi70_bio13 = raster("cc45bi70/cc45bi70/bio_13.tif")
bio13 = crop(cc45bi70_bio13, extent(-12,41.9,34.6,71.4))

cc45bi70_bio14 = raster("cc45bi70/cc45bi70/bio_14.tif")
bio14 = crop(cc45bi70_bio14, extent(-12,41.9,34.6,71.4))

cc45bi70_bio15 = raster("cc45bi70/cc45bi70/bio_15.tif")
bio15 = crop(cc45bi70_bio15, extent(-12,41.9,34.6,71.4))

cc45bi70_bio16 = raster("cc45bi70/cc45bi70/bio_16.tif")
bio16 = crop(cc45bi70_bio16, extent(-12,41.9,34.6,71.4))

cc45bi70_bio17 = raster("cc45bi70/cc45bi70/bio_17.tif")
bio17 = crop(cc45bi70_bio17, extent(-12,41.9,34.6,71.4))

cc45bi70_bio18 = raster("cc45bi70/cc45bi70/bio_18.tif")
bio18 = crop(cc45bi70_bio18, extent(-12,41.9,34.6,71.4))

cc45bi70_bio19 = raster("cc45bi70/cc45bi70/bio_19.tif")
bio19 = crop(cc45bi70_bio19, extent(-12,41.9,34.6,71.4))

myExplFuture_cc45bi70 = stack(bio1 , 
                              bio2 ,
                              bio3 , 
                              bio4 ,
                              bio5 ,
                              bio6 ,
                              bio7 ,
                              bio8 ,
                              bio9 ,
                              bio10, 
                              bio11,
                              bio12,
                              bio13,
                              bio14,
                              bio15,
                              bio16,
                              bio17,
                              bio18,
                              bio19)


myBiomodProjFuture_cc45bi70_1 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_1,
  new.env = myExplFuture_cc45bi70,
  proj.name = 'future_cc45bi70_1',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myFutureProj_cc45bi70_1  <- get_predictions(myBiomodProjFuture_cc45bi70_1)
Future_cc45bi70_1Result  <- calc(myFutureProj_cc45bi70_1  , fun = mean)
writeRaster(Future_cc45bi70_1Result  , filename = "Future_cc45bi70_1Result ", format = "ascii", overwrite = T)

myBiomodProjFuture_cc45bi70_10 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_10,
  new.env = myExplFuture_cc45bi70,
  proj.name = 'future_cc45bi70_10',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myFutureProj_cc45bi70_10 <- get_predictions(myBiomodProjFuture_cc45bi70_10)
Future_cc45bi70_10Result <- calc(myFutureProj_cc45bi70_10 , fun = mean)
writeRaster(Future_cc45bi70_10Result , filename = "Future_cc45bi70_10Result", format = "ascii", overwrite = T)  

myBiomodProjFuture_cc45bi70_2 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_2,
  new.env = myExplFuture_cc45bi70,
  proj.name = 'future_cc45bi70_2',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myFutureProj_cc45bi70_2  <- get_predictions(myBiomodProjFuture_cc45bi70_2)
Future_cc45bi70_2Result  <- calc(myFutureProj_cc45bi70_2  , fun = mean)
writeRaster(Future_cc45bi70_2Result  , filename = "Future_cc45bi70_2Result ", format = "ascii", overwrite = T)

myBiomodProjFuture_cc45bi70_20 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_20,
  new.env = myExplFuture_cc45bi70,
  proj.name = 'future_cc45bi70_20',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myFutureProj_cc45bi70_20 <- get_predictions(myBiomodProjFuture_cc45bi70_20)
Future_cc45bi70_20Result <- calc(myFutureProj_cc45bi70_20 , fun = mean)
writeRaster(Future_cc45bi70_20Result , filename = "Future_cc45bi70_20Result", format = "ascii", overwrite = T)