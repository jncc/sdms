## Tests to check the functionality of the Package function - SDMs
#  Becky Trippier 22/05/2018
#### ----------------------

load("D:/Github/sdms/R/sysdata.rda")

context("SDM tests - general")


#occ <- bngprep(speciesdf = sp_gatewaydata, bngCol = "gridReference", datafrom =
#                 "NBNgateway", minyear = 2007, maxyear = 2014, covarRes = 300)
#sp::coordinates(ev_spdat)<- ~ easting + northing

# get background mask
#load("Z:/Prog100-TerrestrialEvidence/Innov- Natural capital models/Data/BGmasks/amphibian")

#if(exists("mask1km")){
#  background <- mask1km
#  rm(mask1km)
#} else {
#  r <- vars[[1]]
#  r[!is.na(r[])] <- 1
#  background <- r
#  rm(r)
#}

#raster stack of predictor variables - vars
#load(file = "D:/Becky docs/SDMs/varsLCM2015")

##run SDM function
#final_out <- SDMs(occ = ev_spdat, max_tries = 2, lab = "ev", rndm_occ = TRUE)



test_that("function places random occurrences when defined as TRUE", {

})

test_that("test all models run on start up when selected", {

})

test_that("test output files produced", {

})

test_that("max tries can't exceed 100", {

})
