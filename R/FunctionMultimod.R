#' Set up modelling multiple species from NBN gateway or atlas for use in species distribution modelling.
#'
#' This function is to carry out modelling of multiple species using the JNCCSdms package.
#'
#'@param sp_list List of unique species names which you wish to model.
#'@param out_flder The location of the output folder for your models.
#'@param dat_flder The location of the folder containing your species occurrence data.
#'@param max_tries The number of times the number is run.
#'@param datafrom Whether it is data from the "NBNgateway" or "NBNatlas".
#'@param mult_prssr Set up a parallel backend to use multiple processors. As a default this is turned off. Need to ensure the suggested packages have been loaded in order to run this.
#'@return A copy
#'@examples

#'@export


Multi_mod <- function (sp_list = sp_list, #unique list of species names
                  out_flder = "Output/", #output folder for models
                  dat_flder = "Input/",
                  bkgd_flder = "BGmasks/", #location of the folder for the background masks
                  max_tries = 2, #number of model runs
                  datafrom = "NBNgateway",
                  bngCol = "OSGR 2km",
                  mult_prssr = FALSE, #set up multiple processors
                  vars, # predictor variables
                  rndm_occ = TRUE
)
{



alreadyDone <- NULL

## ------------------ setting up done list -----------------------##
#If output files already exist, add to alreadyDone list
for (i in 1:length(sp_list)){
  sp <- sp_list[i]
  lab <- gsub('([[:punct:]])|\\s+','_',sp)
  if (file.exists(paste(out_flder, sp, max_tries, ".csv", sep = ""))
      & file.exists(paste(out_flder, sp, max_tries, ".grd", sep = ""))
      | file.exists(paste(out_flder, sp, max_tries, ".tif", sep = ""))
  ) {
    alreadyDone <- list(alreadyDone, sp)
  }
  if (file.exists(paste(out_flder, lab, max_tries, ".csv", sep = ""))
      & file.exists(paste(out_flder, lab, max_tries, ".grd", sep = ""))
      | file.exists(paste(out_flder, lab, max_tries, ".tif", sep = ""))
  ) {
    alreadyDone <- list(alreadyDone, sp)
  }
}

alreadyDone <- unlist(alreadyDone)
sp_list <- sp_list[!sp_list %in% alreadyDone] #remove alreadyDone species from sp_list


## ------------------ setting up multiple processors -----------------------##

##setup parallel backend to use multiple processors if required
if (mult_prssr == TRUE){
  cl <- parallel::makeCluster(detectCores()-1) #not to overload your computer
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, .libPaths("D:/Rpackages")) #guide each to libpath, as above (line 13)
  packlist <- c("stringi", "rgdal", "raster", "rgeos", "dismo", "kernlab", "dismo", "kernlab", "randomForest", "glmulti", "beepr", "mgcv", "readr", "maptools", "mapdata", "ff", "ffbase", "rnrfa", "pryr")
}
ptm <- proc.time()


## ------------------ run models (non-parallel processing) -----------------------##
if (mult_prssr == FALSE){
for (i in 1:length(sp_list) ) {    #use if running one at a time
sp <- sp_list[i]
lab <- gsub('([[:punct:]])|\\s+','_',sp)

#read in and prepare species occurrence data
if (datafrom == "NBNgatway") {
  gateway_dat <- readr::read_delim(file = paste(dat_flder, sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
  spdat <- bngprep(speciesdf = gateway_dat, precisionCol = "precision", bngCol  = "gridReference", mindata = 5000, minyear = 2007, covarRes = 300)
} else if (datafrom == "NBNatlas"){
  atlas_dat<- read.csv(file=dat_flder, header=TRUE, sep=",", check.names = FALSE, strip.white = TRUE)
  spdat <- bngprep(speciesdf = atlas_dat, datafrom = "NBNatlas", mindata = 5000, minyear = 2007, covarRes = 300)
}
#convert to spatial data frame
sp::coordinates(spdat)<- ~ easting + northing

#Create background
if("taxonGroup" %in% colnames(spdat)) {
  taxon <- spdat$taxonGroup[1]
} else {
  taxon <-readline(paste("what taxon group is", sp, " ?"))#manual prompt to insert taxon group
}
tryCatch(load(file=paste(bkgd_flder, taxon, sep="")), error=function(err) NA)
if(exists("mask1km")){
  background <- mask1km
  rm(mask1km)
} else {
  r <- vars[[1]]
  r[!is.na(r[])] <- 1
  background <- r
  rm(r)
}


##run SDM function
final_out <- SDMs(occ = spdat, varstack = vars, lab = sp, rndm_occ = rndm_occ)
ptm <- proc.time()
beepr::beep()

  }
}

## ------------------ run models (parallel processing) -----------------------##
if (mult_prssr == TRUE){
foreach(i = 1:length(sp_list), .packages = packlist) %dopar% {    #use if using parallel cores

   sp <- sp_list[i]
  lab <- gsub('([[:punct:]])|\\s+','_',sp)

  #read in and prepare species occurrence data
  if (datafrom == "NBNgatway") {
    gateway_dat <- readr::read_delim(file = paste(dat_flder, sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
    spdat <- bngprep(speciesdf = gateway_dat, precisionCol = "precision", bngCol = "gridReference", mindata = 5000, minyear = 2007, covarRes = 300)
  } else if (datafrom == "NBNatlas"){
    atlas_dat<- read.csv(file=dat_flder, header=TRUE, sep=",", check.names = FALSE, strip.white = TRUE)
    spdat <- bngprep(speciesdf = atlas_dat, datafrom = "NBNatlas", mindata = 5000, minyear = 2007, covarRes = 300)
  }
  #convert to spatial data frame
  sp::coordinates(spdat)<- ~ easting + northing

  #Create background
  if("taxonGroup" %in% colnames(spdat)) {
    taxon <- spdat$taxonGroup[1]
  } else {
    #manual prompt to insert taxon group
    taxon <-readline(paste("what taxon group is", sp, " ?"))
  }
  tryCatch(load(file=paste(bkgd_flder, taxon, sep="")), error=function(err) NA)
  if(exists("mask1km")){
    background <- mask1km
    rm(mask1km)
  } else {
    r <- vars[[1]]
    r[!is.na(r[])] <- 1
    background <- r
    rm(r)
  }


  ##run SDM function
  final_out <- SDMs(occ = spdat, varstack = vars, lab = sp, rndm_occ = rndm_occ)
  ptm <- proc.time()

#stop cluster if parrallel processing
parallel::stopCluster(cl)
proc.time() - ptm
beepr::beep()

}
}
}



