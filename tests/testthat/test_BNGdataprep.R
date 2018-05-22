## Tests to check the functionality of the Package functions BNGdataprep and Occurrence
#  Becky Trippier 21/05/2018
#### ----------------------

load("D:/Github/sdms/R/sysdata.rda")

context("BNGdataprep - general")

gatspdat <- bngprep(speciesdf = sp_gatewaydata, bngCol = "gridReference", datafrom =
                        "NBNgateway", minyear = 2007, maxyear = 2014, covarRes = 300)

test_that("function returns dataframe of presence only records", {

  expect_is(gatspdat, "data.frame")
  expect_false(unique(gatspdat$zeroAbundance))

})

test_that("maxyear smaller than minyear warning", {

  expect_error(bngprep(speciesdf = sp_gatewaydata, bngCol = "gridReference", datafrom = "NBNgateway", minyear = 2007, maxyear = 2000, covarRes = 300),"minimum year limit is greater than maximum year")
})

test_that("where no datafrom specified, function will still run", {
  sp_atlasdata$pres <- 1000
  expect_warning(spdat <- bngprep(speciesdf = sp_atlasdata,  bngCol = "OSGR 1km", mindata = 5000, minyear = 2007, maxyear = 2014, covarRes = 300), "datafrom not specified as NBNatlas or NBNgateway." )
  expect_is(spdat, "data.frame")
})
#-----------------------------------------#

context("BNGdataprep - Gateway Example")

test_that("function returns dataframe of presence only records", {

  expect_is(gatspdat, "data.frame")
  expect_false(unique(gatspdat$zeroAbundance))

})

test_that("Subsetting dates returns within the correct range of years", {

  expect_true(min(gatspdat$year) == 2007)
  expect_true(max(gatspdat$year) == 2014)
  expect_false(unique(is.na(gatspdat$year)))

  })

test_that("easting and northing fields added and complete", {

  x<- gatspdat$easting
  y<- gatspdat$northing
  expect_true(is.numeric(x))
  expect_true(is.numeric(y))
  expect_false(unique(is.na(x)))
  expect_false(unique(is.na(y)))

})

#-----------------------------------------#
context("BNGdataprep - Atlas Example")
atspdat <- bngprep(speciesdf = sp_atlasdata,  bngCol = "OSGR 1km", datafrom = "NBNatlas", mindata = 5000, minyear = 2007, maxyear = 2014, covarRes = 300)

test_that("function returns dataframe of presence only records", {
  expect_is(atspdat, "data.frame")
  expect_true(unique(atspdat$`Occurrence status`)== "present")

})

test_that("Subsetting dates returns within the correct range of years", {
  expect_true(min(atspdat$year) == 2007)
  expect_true(max(atspdat$year) == 2014)
  expect_false(unique(is.na(atspdat$year)))

})

test_that("easting and northing fields added and complete", {

  x<- atspdat$easting
  y<- atspdat$northing
  expect_true(is.numeric(x))
  expect_true(is.numeric(y))
  expect_false(unique(is.na(x)))
  expect_false(unique(is.na(y)))

})

test_that("bngCol is supplied for NBNatlas data and is a valid column", {

  expect_error(bngprep(speciesdf = sp_atlasdata, datafrom= "NBNatlas"),'argument \"bngCol\" is missing, with no default')
 expect_error(bngprep(speciesdf = sp_atlasdata, datafrom= "NBNatlas", bngCol = "OSGR"),"Unaccepted bngCol specified.")

})

test_that("PrecionCol error generated for atlas data, but function still runs", {

expect_message(out <- bngprep(speciesdf = sp_atlasdata, datafrom= "NBNatlas", bngCol = "OSGR 1km", precisionCol = "testpres"),"unnecessary argument - do not specify precisionCol for data from NBNatlas")
expect_is(out, "data.frame")

})

