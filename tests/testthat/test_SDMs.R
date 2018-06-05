## Tests to check the functionality of the Package function - SDMs
#  Becky Trippier 05/06/2018
#### ------------------------------------------------------------####
context("SDMs tests")

### set working directory for tests
start <- getwd()
setwd(tempdir())
dir.create("Outputs")


#### get test data

## test occurrence data
data(ng_data)
ngspdat <- bngprep(speciesdf = ng_data,  bngCol = "OSGR", datafrom = "NBNatlas", mindata = 5000, minyear = 2007, maxyear = 2012, covarRes = 300)
sp::coordinates(ngspdat)<- ~ easting + northing

## raster stack of predictor variables - vars
data(vars)

#load background mask
data(background)

###### ---------------------------------------------------------------------------------------------------- ######
#### Function tests

test_that("model input error when these are incorrectly defined", {
  expect_error(SDMs(occ = ngspdat, max_tries = 1, models = "Bioclim", lab = "test", bckg = background, rndm_occ = TRUE, varstack = vars), "Model specification contains an unexpected value, please check model names input. Operation terminated.")

})

test_that("function places random occurrences when defined as TRUE", {
  expect_condition(SDMs(occ = ngspdat, max_tries = 1, models = "BioClim", lab = "test", bckg = background, rndm_occ = TRUE), "Random occurrences placed")

})

test_that("test output files produced", {

oldw <- getOption("warn")
options(warn=-1)

SDMs(occ = ngspdat, max_tries = 1, lab = "test", bckg = background, rndm_occ = TRUE, varstack = vars)
options(warn = oldw)

#expected output files
expect_true(exists("all_evals") == TRUE)
expect_true(exists("all_models") == TRUE)
expect_true(exists("all_predicts") == TRUE)

})

# remove temporary files and return to working directory
unlink("Outputs", recursive=TRUE)
setwd(start)


