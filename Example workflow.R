#### Example workflow for the JNCCsdms package ####
##
## Becky Trippier 15/05/2018
##
## -------------------------------------------------
getwd()
setwd("D:/Becky docs/SDMs")

#Import species data and convert with bngprep function - gateway
lc_gateway <- readr::read_delim(file = paste("D:/Becky docs/SDMs/Limenitis_camilla.txt", sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE)
lc_spdat <- bngprep(speciesdf = lc_gateway, precisionCol = "precision", bngCol = "gridReference", datafrom = "NBNgateway", mindata = 5000, minyear = 2007, covarRes = 300)
##Single species modelling
#Import species data and convert with bngprep function - atlas
puff_atlasdata<- read.csv(file="D:/Becky docs/SDMs/records-2018-05-17.csv", header=TRUE, sep=",", check.names = FALSE, strip.white = TRUE)

puff_spdat <- bngprep(speciesdf = puff_atlasdata,  bngCol = "OSGR 2km", datafrom = "NBNatlas", mindata = 5000, minyear = 2007, covarRes = 300)


#make file geospatial data frame
sp::coordinates(puff_spdat)<- ~ easting + northing

# get background mask
load("Z:/Prog100-TerrestrialEvidence/Innov- Natural capital models/Data/BGmasks/amphibian")

if(exists("mask1km")){
  background <- mask1km
  rm(mask1km)
} else {
  r <- vars[[1]]
  r[!is.na(r[])] <- 1
  background <- r
  rm(r)
}
## if tif is available
#background <- raster::raster("Projected rasters/bg_mask.tif")

#raster stack of predictor variables - vars
load(file = "D:/Becky docs/SDMs/varsLCM2015")

##run SDM function

final_out <- SDMs(occ = puff_spdat, max_tries = 2, lab = "puff2km", rndm_occ = TRUE)
ptm <- proc.time()

#-----------------------------------------------------------------------------------------------------

# Multiple species modelling

#get species list
sp_list <- readr::read_delim("C4b_2015_specieslist.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

#unique records
sp_list <- unique(sp_list$X1)

