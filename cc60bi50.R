library(biomod2)
library(raster)
library(ggplot2) 

cc60bi50_bio1 = raster("cc60bi50/cc60bi50/bio_1.tif")
bio1 = crop(cc60bi50_bio1, extent(-12,41.9,34.6,71.4))

cc60bi50_bio2 = raster("cc60bi50/cc60bi50/bio_2.tif")
bio2 = crop(cc60bi50_bio2, extent(-12,41.9,34.6,71.4))

cc60bi50_bio3 = raster("cc60bi50/cc60bi50/bio_3.tif")
bio3 = crop(cc60bi50_bio3, extent(-12,41.9,34.6,71.4))

cc60bi50_bio4 = raster("cc60bi50/cc60bi50/bio_4.tif")
bio4 = crop(cc60bi50_bio4, extent(-12,41.9,34.6,71.4))

cc60bi50_bio5 = raster("cc60bi50/cc60bi50/bio_5.tif")
bio5 = crop(cc60bi50_bio5, extent(-12,41.9,34.6,71.4))

cc60bi50_bio6 = raster("cc60bi50/cc60bi50/bio_6.tif")
bio6 = crop(cc60bi50_bio6, extent(-12,41.9,34.6,71.4))

cc60bi50_bio7 = raster("cc60bi50/cc60bi50/bio_7.tif")
bio7 = crop(cc60bi50_bio7, extent(-12,41.9,34.6,71.4))

cc60bi50_bio8 = raster("cc60bi50/cc60bi50/bio_8.tif")
bio8 = crop(cc60bi50_bio8, extent(-12,41.9,34.6,71.4))

cc60bi50_bio9 = raster("cc60bi50/cc60bi50/bio_9.tif")
bio9 = crop(cc60bi50_bio9, extent(-12,41.9,34.6,71.4))

cc60bi50_bio10 = raster("cc60bi50/cc60bi50/bio_10.tif")
bio10 = crop(cc60bi50_bio10, extent(-12,41.9,34.6,71.4))

cc60bi50_bio11 = raster("cc60bi50/cc60bi50/bio_11.tif")
bio11 = crop(cc60bi50_bio11, extent(-12,41.9,34.6,71.4))

cc60bi50_bio12 = raster("cc60bi50/cc60bi50/bio_12.tif")
bio12 = crop(cc60bi50_bio12, extent(-12,41.9,34.6,71.4))

cc60bi50_bio13 = raster("cc60bi50/cc60bi50/bio_13.tif")
bio13 = crop(cc60bi50_bio13, extent(-12,41.9,34.6,71.4))

cc60bi50_bio14 = raster("cc60bi50/cc60bi50/bio_14.tif")
bio14 = crop(cc60bi50_bio14, extent(-12,41.9,34.6,71.4))

cc60bi50_bio15 = raster("cc60bi50/cc60bi50/bio_15.tif")
bio15 = crop(cc60bi50_bio15, extent(-12,41.9,34.6,71.4))

cc60bi50_bio16 = raster("cc60bi50/cc60bi50/bio_16.tif")
bio16 = crop(cc60bi50_bio16, extent(-12,41.9,34.6,71.4))

cc60bi50_bio17 = raster("cc60bi50/cc60bi50/bio_17.tif")
bio17 = crop(cc60bi50_bio17, extent(-12,41.9,34.6,71.4))

cc60bi50_bio18 = raster("cc60bi50/cc60bi50/bio_18.tif")
bio18 = crop(cc60bi50_bio18, extent(-12,41.9,34.6,71.4))

cc60bi50_bio19 = raster("cc60bi50/cc60bi50/bio_19.tif")
bio19 = crop(cc60bi50_bio19, extent(-12,41.9,34.6,71.4))

myExplFuture_cc60bi50 = stack(bio1 , 
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


myBiomodProjFuture_cc60bi50_1 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_1,
  new.env = myExplFuture_cc60bi50,
  proj.name = 'future_cc60bi50_1',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myFutureProj_cc60bi50_1  <- get_predictions(myBiomodProjFuture_cc60bi50_1)
Future_cc60bi50_1Result  <- calc(myFutureProj_cc60bi50_1  , fun = mean)
writeRaster(Future_cc60bi50_1Result  , filename = "Future_cc60bi50_1Result ", format = "ascii", overwrite = T)

myBiomodProjFuture_cc60bi50_10 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_10,
  new.env = myExplFuture_cc60bi50,
  proj.name = 'future_cc60bi50_10',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myFutureProj_cc60bi50_10 <- get_predictions(myBiomodProjFuture_cc60bi50_10)
Future_cc60bi50_10Result <- calc(myFutureProj_cc60bi50_10 , fun = mean)
writeRaster(Future_cc60bi50_10Result , filename = "Future_cc60bi50_10Result", format = "ascii", overwrite = T)  

myBiomodProjFuture_cc60bi50_2 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_2,
  new.env = myExplFuture_cc60bi50,
  proj.name = 'future_cc60bi50_2',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myFutureProj_cc60bi50_2  <- get_predictions(myBiomodProjFuture_cc60bi50_2)
Future_cc60bi50_2Result  <- calc(myFutureProj_cc60bi50_2  , fun = mean)
writeRaster(Future_cc60bi50_2Result  , filename = "Future_cc60bi50_2Result ", format = "ascii", overwrite = T)

myBiomodProjFuture_cc60bi50_20 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_20,
  new.env = myExplFuture_cc60bi50,
  proj.name = 'future_cc60bi50_20',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myFutureProj_cc60bi50_20 <- get_predictions(myBiomodProjFuture_cc60bi50_20)
Future_cc60bi50_20Result <- calc(myFutureProj_cc60bi50_20 , fun = mean)
writeRaster(Future_cc60bi50_20Result , filename = "Future_cc60bi50_20Result", format = "ascii", overwrite = T)
