## Tests to check the functionality of the Package function - SDMs
#  Becky Trippier 22/05/2018
#### ------------------------------------------------------------####

#### get test data

## test occurrence data
ng_data<- read.csv(file="Data/Inputs/Notonecta_glauca.csv", header=TRUE, sep=",", check.names = FALSE, strip.white = TRUE)
ngspdat <- bngprep(speciesdf = ng_data,  bngCol = "OSGR", datafrom = "NBNatlas", mindata = 5000, minyear = 2007, maxyear = 2012, covarRes = 300)
sp::coordinates(ngspdat)<- ~ easting + northing

## raster stack of predictor variables - vars
#get UK extent
UK <- ggplot2::map_data(map = "world", region = "UK")
max.lat <- ceiling(max(UK$lat))
min.lat <- floor(min(UK$lat))
max.lon <- ceiling(max(UK$long))
min.lon <- floor(min(UK$long))
extent <- raster::extent(x = c(min.lon, max.lon, min.lat, max.lat))

#get variables data
bio<-raster::getData('worldclim',var='bio',res=5,lon=-2,lat=40)
bio <- bio[[c("bio1","bio12")]]
names(bio) <- c("Temp","Prec")

#crop to uk
bio<-raster::crop(bio,extent)

##change to easting northing
vars <- projectRaster(bio, crs="+init=epsg:27700")

#load background mask

load("Data/BGmasks/insect - true bug (Hemiptera)")
background <- mask1km
rm(mask1km)

###### ---------------------------------------------------------------------------------------------------- ######
#### Function tests

test_that("model input error when these are incorrectly defined", {
  expect_error(SDMs(occ = ngspdat, max_tries = 1, models = "Bioclim", lab = "test", bckg = background, rndm_occ = TRUE, varstack = vars), "Model specification contains an unexpected value, please check model names input. Operation terminated.")

})

test_that("function places random occurrences when defined as TRUE", {
  expect_condition(SDMs(occ = ngspdat, max_tries = 1, models = "BioClim", lab = "test", bckg = background, rndm_occ = TRUE), "Random occurrences placed")

})

test_that("test output files produced and all models run on start up", {

library(glmulti)
SDMs(occ = ngspdat, max_tries = 1, lab = "test", bckg = background, rndm_occ = TRUE)
expect_true(exists("all_evals") == TRUE)
expect_true(exists("all_models") == TRUE)
expect_true(exists("all_predicts") == TRUE)

evals <- as.data.frame(unlist(all_evals[2]))
colnames(evals)[1] <- "score"
expect_true(all(evals$score != "NA") == TRUE)
})






