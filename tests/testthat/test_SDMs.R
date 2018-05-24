## Tests to check the functionality of the Package function - SDMs
#  Becky Trippier 22/05/2018
#### ----------------------

#get occurrence data for tests
load("D:/Github/sdms/R/sysdata.rda")

gatspdat <- bngprep(speciesdf = sp_gatewaydata, bngCol = "gridReference", datafrom =
                      "NBNgateway", minyear = 2007, maxyear = 2014, covarRes = 300)
sp::coordinates(gatspdat)<- ~ easting + northing

#raster stack of predictor variables - vars
load(file = "D:/Becky docs/SDMs/varsLCM2015")

######
SDMs(occ = gatspdat, max_tries = 1, models = "BioClim", lab = "test", bckg = gatbackground, rndm_occ = TRUE, varstack = bio)

test_that("model input error when these are incorrectly defined", {
  expect_error(SDMs(occ = gatspdat, max_tries = 1, models = "Bioclim", lab = "test", bckg = gatbackground, rndm_occ = TRUE), "Model specification contains an unexpected value, please check model names input. Operation terminated.")

})

test_that("function places random occurrences when defined as TRUE", {
  expect_condition(SDMs(occ = gatspdat, max_tries = 1, models = "BioClim", lab = "test", bckg = gatbackground, rndm_occ = TRUE), "Random occurrences placed")

})

library(glmulti)
SDMs(occ = gatspdat, max_tries = 1, lab = "test", bckg = gatbackground, rndm_occ = TRUE)

test_that("test output files produced", {
expect_true(exists("all_evals") == TRUE)
expect_true(exists("all_models") == TRUE)
expect_true(exists("all_predicts") == TRUE)
})

test_that("test all models run on start up when selected", {
  evals <- as.data.frame(unlist(all_evals[2]))
  colnames(evals)[1] <- "score"
  expect_true(all(evals$score != "NA") == TRUE)

})





