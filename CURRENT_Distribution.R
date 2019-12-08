library(biomod2)
library(raster)

DataSpecies <- read.csv("Climatic_Niches.csv")
head(DataSpecies)

# the name of studied species
myRespName1 <- 'Pannonian'
myRespName2 <- 'Carpathian'

# the presence/absences data for our species
myResp1 <- as.numeric(DataSpecies[,myRespName1])
myResp2 <- as.numeric(DataSpecies[,myRespName2])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]

# load the environmental raster layers (could be .img, ArcGIS
# rasters or any supported format by the raster package)

bio1  = raster("bio_30s/bio_30s/bio_1.bil")
bio2  = raster("bio_30s/bio_30s/bio_2.bil")
bio3  = raster("bio_30s/bio_30s/bio_3.bil")
bio4  = raster("bio_30s/bio_30s/bio_4.bil")
bio5  = raster("bio_30s/bio_30s/bio_5.bil")
bio6  = raster("bio_30s/bio_30s/bio_6.bil")
bio7  = raster("bio_30s/bio_30s/bio_7.bil")
bio8  = raster("bio_30s/bio_30s/bio_8.bil")
bio9  = raster("bio_30s/bio_30s/bio_9.bil")
bio10 = raster("bio_30s/bio_30s/bio_10.bil")
bio11 = raster("bio_30s/bio_30s/bio_11.bil")
bio12 = raster("bio_30s/bio_30s/bio_12.bil")
bio13 = raster("bio_30s/bio_30s/bio_13.bil")
bio14 = raster("bio_30s/bio_30s/bio_14.bil")
bio15 = raster("bio_30s/bio_30s/bio_15.bil")
bio16 = raster("bio_30s/bio_30s/bio_16.bil")
bio17 = raster("bio_30s/bio_30s/bio_17.bil")
bio18 = raster("bio_30s/bio_30s/bio_18.bil")
bio19 = raster("bio_30s/bio_30s/bio_19.bil")

bio1  = crop(bio1  , extent(-12,41.9,34.6,71.4))
bio2  = crop(bio2  , extent(-12,41.9,34.6,71.4))
bio3  = crop(bio3  , extent(-12,41.9,34.6,71.4))
bio4  = crop(bio4  , extent(-12,41.9,34.6,71.4))
bio5  = crop(bio5  , extent(-12,41.9,34.6,71.4))
bio6  = crop(bio6  , extent(-12,41.9,34.6,71.4))
bio7  = crop(bio7  , extent(-12,41.9,34.6,71.4))
bio8  = crop(bio8  , extent(-12,41.9,34.6,71.4))
bio9  = crop(bio9  , extent(-12,41.9,34.6,71.4))
bio10 = crop(bio10 , extent(-12,41.9,34.6,71.4))
bio11 = crop(bio11 , extent(-12,41.9,34.6,71.4))
bio12 = crop(bio12 , extent(-12,41.9,34.6,71.4))
bio13 = crop(bio13 , extent(-12,41.9,34.6,71.4))
bio14 = crop(bio14 , extent(-12,41.9,34.6,71.4))
bio15 = crop(bio15 , extent(-12,41.9,34.6,71.4))
bio16 = crop(bio16 , extent(-12,41.9,34.6,71.4))
bio17 = crop(bio17 , extent(-12,41.9,34.6,71.4))
bio18 = crop(bio18 , extent(-12,41.9,34.6,71.4))
bio19 = crop(bio19 , extent(-12,41.9,34.6,71.4))

myExpl = stack(bio1 ,
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

########################################################################################

myBiomodData_1 <- BIOMOD_FormatingData(resp.var = myResp1,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName1,
                                       PA.nb.rep = 1)

myBiomodData_1

myBiomodOption <- BIOMOD_ModelingOptions()

myBiomodModelOut_1 <- BIOMOD_Modeling( myBiomodData_1, 
                                       models = c('GLM','MAXENT.Phillips'), 
                                       models.options = myBiomodOption,
                                       NbRunEval=10,
                                       DataSplit=70,
                                       Prevalence=0.5,
                                       VarImport=3,
                                       models.eval.meth = c('TSS','ROC'),
                                       SaveObj = TRUE,
                                       rescal.all.models = TRUE,
                                       do.full.models = FALSE,
                                       modeling.id = paste(myRespName1,"FirstModeling_1",sep=""))

myBiomodModelOut_1

myBiomodModelEval_1 <- get_evaluations(myBiomodModelOut_1)

dimnames(myBiomodModelEval_1)

get_variables_importance(myBiomodModelOut_1)

myBiomodEM_1 <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut_1,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.7),
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )

myBiomodEM_1

get_evaluations(myBiomodEM_1)

myBiomodProj_1 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_1,
  new.env = myExpl,
  proj.name = 'current_1',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

myBiomodProj_1

myCurrentProj_1 <- get_predictions(myBiomodProj_1)

presentResult_1 <- calc(myCurrentProj_1, fun = mean)

plot(presentResult_1)

writeRaster(presentResult_1, filename = "Pannonian_Present_1", format = "ascii", overwrite = T);

myBiomodEF_1 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM_1,
  projection.output = myBiomodProj_1)

###
myBiomodData_10 <- BIOMOD_FormatingData(resp.var = myResp1,
                                        expl.var = myExpl,
                                        resp.xy = myRespXY,
                                        resp.name = myRespName1,
                                        PA.nb.rep = 10)

myBiomodData_10

myBiomodOption <- BIOMOD_ModelingOptions()

myBiomodModelOut_10 <- BIOMOD_Modeling( myBiomodData_10, 
                                        models = c('GBM', 'CTA', 'FDA', 'MARS', 'RF'), 
                                        models.options = myBiomodOption,
                                        NbRunEval=1,
                                        DataSplit=70,
                                        Prevalence=0.5,
                                        VarImport=3,
                                        models.eval.meth = c('TSS','ROC'),
                                        SaveObj = TRUE,
                                        rescal.all.models = TRUE,
                                        do.full.models = FALSE,
                                        modeling.id = paste(myRespName1,"FirstModeling_10",sep=""))

myBiomodModelOut_10

myBiomodModelEval_10 <- get_evaluations(myBiomodModelOut_10)

dimnames(myBiomodModelEval_10)

get_variables_importance(myBiomodModelOut_10)

myBiomodEM_10 <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut_10,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.7),
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )

myBiomodEM_10

get_evaluations(myBiomodEM_10)

myBiomodProj_10 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_10,
  new.env = myExpl,
  proj.name = 'current_10',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

myBiomodProj_10

myCurrentProj_10 <- get_predictions(myBiomodProj_10)

presentResult_10 <- calc(myCurrentProj_10, fun = mean)

plot(presentResult_10)

writeRaster(presentResult_10, filename = "Pannonian_Present_10", format = "ascii", overwrite = T);

myBiomodEF_10 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM_10,
  projection.output = myBiomodProj_10)

###
myBiomodData_2 <- BIOMOD_FormatingData(resp.var = myResp2,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName2,
                                       PA.nb.rep = 1)

myBiomodData_2

myBiomodOption <- BIOMOD_ModelingOptions()

myBiomodModelOut_2 <- BIOMOD_Modeling( myBiomodData_2, 
                                       models = c('GLM','MAXENT.Phillips'), 
                                       models.options = myBiomodOption,
                                       NbRunEval=10,
                                       DataSplit=70,
                                       Prevalence=0.5,
                                       VarImport=3,
                                       models.eval.meth = c('TSS','ROC'),
                                       SaveObj = TRUE,
                                       rescal.all.models = TRUE,
                                       do.full.models = FALSE,
                                       modeling.id = paste(myRespName2,"FirstModeling_2",sep=""))

myBiomodModelOut_2

myBiomodModelEval_2 <- get_evaluations(myBiomodModelOut_2)

dimnames(myBiomodModelEval_2)

get_variables_importance(myBiomodModelOut_2)

myBiomodEM_2 <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut_2,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.7),
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )

myBiomodEM_2

get_evaluations(myBiomodEM_2)

myBiomodProj_2 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_2,
  new.env = myExpl,
  proj.name = 'current_2',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

myBiomodProj_2

myCurrentProj_2 <- get_predictions(myBiomodProj_2)

presentResult_2 <- calc(myCurrentProj_2, fun = mean)

plot(presentResult_2)

writeRaster(presentResult_2, filename = "Carpathian_Present_2", format = "ascii", overwrite = T);

myBiomodEF_2 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM_2,
  projection.output = myBiomodProj_2)

###
myBiomodData_20 <- BIOMOD_FormatingData(resp.var = myResp2,
                                        expl.var = myExpl,
                                        resp.xy = myRespXY,
                                        resp.name = myRespName2,
                                        PA.nb.rep = 10)

myBiomodData_20
myBiomodOption <- BIOMOD_ModelingOptions()

myBiomodModelOut_20 <- BIOMOD_Modeling( myBiomodData_20, 
                                        models = c('GBM', 'CTA', 'FDA', 'MARS', 'RF'), 
                                        models.options = myBiomodOption,
                                        NbRunEval=1,
                                        DataSplit=70,
                                        Prevalence=0.5,
                                        VarImport=3,
                                        models.eval.meth = c('TSS','ROC'),
                                        SaveObj = TRUE,
                                        rescal.all.models = TRUE,
                                        do.full.models = FALSE,
                                        modeling.id = paste(myRespName2,"FirstModeling_20",sep=""))




myBiomodModelOut_20

myBiomodModelEval_20 <- get_evaluations(myBiomodModelOut_20)

dimnames(myBiomodModelEval_20)

get_variables_importance(myBiomodModelOut_20)

myBiomodEM_20 <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut_20,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.7),
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )

myBiomodEM_20

get_evaluations(myBiomodEM_20)

myBiomodProj_20 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut_20,
  new.env = myExpl,
  proj.name = 'current_20',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

myBiomodProj_20

myCurrentProj_20 <- get_predictions(myBiomodProj_20)

presentResult_20 <- calc(myCurrentProj_20, fun = mean)

plot(presentResult_20)

writeRaster(presentResult_20, filename = "Carpathian_Present_20", format = "ascii", overwrite = T);

myBiomodEF_20 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM_20,
  projection.output = myBiomodProj_20)

########################################################################################
