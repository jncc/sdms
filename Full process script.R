## Currently this is set up as a script that is quite JNCC-specific and not easily generalisable. In version 0.2, the loop and output aggregation should be added to the SDMs function so that it is self contained. This should also help control the memory issues that I suspect arise from the lists that are appended to in the SDMs function, and then converted to global variables with list2env. Also need to add the package dependencies to the documentation for each function.


## Biodiversity Module
## Applying to multi taxa
## v9 21/11/2017- Main functions seperated and called from this script to allow
##                centralised changes

##Set up working directory and libraries
#.libPaths("D:/Rpackages")

###------------To do!! ---------

### Check pwd effectiveness
### Add NI to grids and vars

###---------
#library(stringi)
#library(rgdal)
#library(raster)
#library(rgeos)
#library(gbm)
#library(dismo)
#library(kernlab)
#library(randomForest)
#library(glmulti)
#library(beepr)
#library(mgcv)
#library(readr)
#library(maptools)
#library(mapdata)
#library(ff)
#library(ffbase)
#library(rnrfa)
#library(pryr)
# If you want to use parallel cores for multiple species, also
#library(foreach)
#library(doParallel)
packlist <- c("stringi", "rgdal", "raster", "rgeos", "dismo", "kernlab", "dismo", "kernlab", "randomForest", "glmulti",
              "beepr", "mgcv", "readr", "maptools", "mapdata", "ff", "ffbase", "rnrfa", "pryr")

###Set model variables
#covarResm <- 300
#lowestResm <- 10000
max_tries <- 20

ensembleSDM <- function(speciesData, envData, models, covarRes = resolution(envData)) {}

##Upload NBN data from txt file (etc.), or sp_list can be added manually
#sp_list <- read_delim("C4b_2015_specieslist_1sthalf.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
#sp_list <- unique(sp_list$X1)
#Remove any which have been done
#alreadyDone <- NULL
#for (i in 1:length(sp_list)){
#If output files already exist, add to alreadyDone list
#  sp <- sp_list[i]
#  lab <- gsub('([[:punct:]])|\\s+','_',sp)
#  if (file.exists(paste("ModuleOutputs/EAU/", sp, max_tries, ".csv", sep = ""))
#      & file.exists(paste("ModuleOutputs/EAU/", sp, max_tries, ".grd", sep = ""))
#      | file.exists(paste("ModuleOutputs/EAU/", sp, max_tries, ".tif", sep = ""))
#  ) {
#    alreadyDone <- list(alreadyDone, sp)
#  }
#  if (file.exists(paste("ModuleOutputs/EAU/", lab, max_tries, ".csv", sep = ""))
#      & file.exists(paste("ModuleOutputs/EAU/", lab, max_tries, ".grd", sep = ""))
#      | file.exists(paste("ModuleOutputs/EAU/", lab, max_tries, ".tif", sep = ""))
#  ) {
#    alreadyDone <- list(alreadyDone, sp)
#  }
#}
#alreadyDone <- unlist(alreadyDone)
#sp_list <- sp_list[!sp_list %in% alreadyDone] #remove alreadyDone species from sp_list

lab <- "zostera"

#Import environmental predictor rasters

# #setup parallel backend to use multiple processors if required
# cl <- makeCluster(detectCores()-1) #not to overload your computer
# registerDoParallel(cl)
# clusterEvalQ(cl, .libPaths("D:/Rpackages")) #guide each to libpath, as above (line 13)

ptm <- proc.time()

# foreach(i = 1:length(sp_list), .packages = packlist) %dopar% {    #use if using parallel cores
#for (i in 1:length(sp_list) ) {    #use if running one at a time
#sp <- sp_list[i]
#lab <- gsub('([[:punct:]])|\\s+','_',sp)

#Import species data
#spdat <- read_delim(file = paste("", lab, ".txt", sep = ""),
#                    "\t", escape_double = FALSE, trim_ws = TRUE)
# If BNG data and need to convert to eastings/northings
#spdat <- bngprep(speciesdf = spdat, precisionCol = "precision", bngCol = "gridReference", datafrom = "NBNgateway", mindata = 5000, covarRes = covarResm)

#Create background
#if("taxonGroup" %in% colnames(spdat)) {
#  taxon <- spdat$taxonGroup[1]
#} else {
#  #manually insert taxon group
#  taxon <- "bird"
#}
#tryCatch(load(file=paste("BGmasks", taxon, sep="")), error=function(err) NA)
#if(exists("mask1km")){
#  background <- mask1km
#  rm(mask1km)
#} else {
#  r <- vars[[1]]
#  r[!is.na(r[])] <- 1
#  background <- r
#  rm(r)
#}
background <- raster::raster("Projected rasters/bg_mask.tif")


#Repeat runs to account for random values
all_predicts <- NULL
all_models <-NULL
all_evals <- NULL
tries <- 0
occurrence <- rgdal::readOGR("DataInputs", "OSPAR2015_GB_ZosteraPoints")
occurrence <- rgdal::spTransform(occurrence, crs(dem))
pres.pts <- as.data.frame(coordinates(occurrence))
rownames(pres.pts) <- 1:nrow(pres.pts)
colnames(pres.pts) <- c("x", "y")

ptm <- proc.time()

while(tries < max_tries) {
  #rasterOptions(tmpdir= './Rtmpdir') #Change raster temp directory so can be emptied after each loop, avoid filling up disk

  # If neccessary, add random amount to eastings/northings to place coordinates
  #occurrence <- randomOcc(presframe = spdat, precisionCol = "precision")

  #pres.pts <- data.frame(occurrence[,c("easting","northing")])

  #coordinates(occurrence)<- ~ easting + northing



  #Run models
  ###
  ## !!! NB You MUST have maxent.jar in your dismo/java folder if you want to use MaxEnt !!! ##
  ###
  mods <- c("MaxEnt", "BioClim", "SVM", "RF", "GLM", "GAM", "BRT") #which models to include
  best_out <- SDMs()
  list2env(best_out ,.GlobalEnv)
  tries <- tries+1
  print(c(tries, all_evals[[2]][[8]]))
}

#Export model evaluation results
Mean_predict <- Reduce(`+`, all_predicts) / length(all_predicts)
Mean_predict <- matrix(unlist(Mean_predict[,]), nrow=nrow(vars[[1]]), ncol=ncol(vars[[1]]), byrow=FALSE)
Mean_predict <- setValues(vars[[1]], Mean_predict)
Models <- c(sort(mods), "Best")
all_evals <- cbind(data.frame(Models), data.frame(matrix(unlist(all_evals), nrow=length(mods)+1, byrow=F), stringsAsFactors=FALSE))
plot(Mean_predict, main = sp)
#beep(0)
proc.time() - ptm
mem_used()

save(all_models, file = paste("Outputs/", lab, tries, "models", sep=""))
writeRaster(Mean_predict, file = paste("Outputs/", lab, tries, sep=""), format="GTiff")
write.csv(all_evals, file = paste("Outputs/", lab, tries, ".csv", sep="") )

#}

beep()
#stop cluster if parrallel processing
stopCluster(cl)
proc.time() - ptm

