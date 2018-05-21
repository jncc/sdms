## Tests to check the functionality of the Package functions BNGdataprep and Occurrence
#  Becky Trippier 21/05/2018
#### ----------------------



context("BNGdataprep - general")

test_that("function returns dataframe of presence only records", {
  gatspdat <- bngprep(speciesdf = sp_gatewaydata, bngCol = "gridReference", datafrom =
                        "NBNgateway", minyear = 2007, maxyear = 2014, covarRes = 300)
  expect_is(gatspdat, "data.frame")
  expect_false(unique(gatspdat$zeroAbundance))

})

#test_that("maxyear smaller than minyear warning", {

#  expect_error(bngprep(speciesdf = sp_gatewaydata, bngCol = "gridReference", datafrom = "NBNgateway", #minyear = 2007, maxyear = 2000, covarRes = 300),"minimum year limit is greater than maximum year")
#})


#-----------------------------------------#

context("BNGdataprep - Gateway Example")

test_that("function returns dataframe of presence only records", {
  gatspdat <- bngprep(speciesdf = sp_gatewaydata, bngCol = "gridReference", datafrom =
                        "NBNgateway", minyear = 2007, maxyear = 2014, covarRes = 300)
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

test_that("function returns dataframe of presence only records", {
  atspdat <- bngprep(speciesdf = sp_atlasdata,  bngCol = "OSGR 1km", datafrom = "NBNatlas", mindata = 5000, minyear = 2007, maxyear = 2014, covarRes = 300)
  expect_is(atspdat, "data.frame")
  expect_true(unique(atspdat$`Occurrence status`)== "present")

})

test_that("Subsetting dates returns within the correct range of years", {
  expect_true(min(atspdat$Year) == 2007)
  expect_true(max(atspdat$Year) == 2014)
  expect_false(unique(is.na(atspdat$Year)))

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



})

test_that("Prescol", {



})

