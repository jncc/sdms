## Tests to check the functionality of the Package functions BNGdataprep
#  Becky Trippier 21/05/2018
#### ----------------------

context("BNGdataprep - Atlas Example")
data("ng_data")
ngspdat <- bngprep(speciesdf = ng_data,  bngCol = "OSGR", datafrom = "NBNatlas", mindata = 5000, minyear = 2007, maxyear = 2012, covarRes = 300)

test_that("function returns dataframe of presence only records", {
  expect_is(ngspdat, "data.frame")
  expect_true(unique(ngspdat$`Occurrence status`)== "present")

})

test_that("Subsetting dates returns within the correct range of years", {
  expect_true(min(ngspdat$year) == 2007)
  expect_true(max(ngspdat$year) == 2012)
  expect_false(unique(is.na(ngspdat$year)))

})

test_that("easting and northing fields added and complete", {

  x<- ngspdat$easting
  y<- ngspdat$northing
  expect_true(is.numeric(x))
  expect_true(is.numeric(y))
  expect_false(unique(is.na(x)))
  expect_false(unique(is.na(y)))

})


test_that("PrecionCol error generated for atlas data, but function still runs", {

expect_message(out <- bngprep(speciesdf = ng_data, datafrom= "NBNatlas", bngCol = "OSGR", precisionCol = "testpres"),"unnecessary argument - do not specify precisionCol for data from NBNatlas")
expect_is(out, "data.frame")

})

test_that("Duplicated field names are handled", {
  expect_true(unique(names(ng_data)[44:45] != names(ngspdat)[44:45]))
})


test_that("maxyear smaller than minyear warning", {
  expect_error(bngprep(speciesdf = ng_data,  bngCol = "OSGR", datafrom = "NBNatlas", minyear = 2007, maxyear = 2000, covarRes = 300),"minimum year limit is greater than maximum year")
})


test_that("where no datafrom specified, function will still run", {
  ng_data$precision <- 1000
  names(ng_data)[names(ng_data) == "Year"] <- "year"
  expect_warning(spdat <- bngprep(speciesdf = ng_data,  bngCol = "OSGR", mindata = 5000, minyear = 2007, maxyear = 2014, covarRes = 300), "datafrom not specified as NBNatlas or NBNgateway." )
  expect_is(spdat, "data.frame")
})
