## Tests to check the functionality of the Package function - Occurrence
#  Becky Trippier 22/05/2018
#### ----------------------

context("Occurrence tests - general")
ng_data<- read.csv(file="Inputs/Notonecta_glauca.csv", header=TRUE, sep=",", check.names = FALSE, strip.white = TRUE)

occ <- bngprep(speciesdf = ng_data, bngCol = "OSGR", datafrom =
                      "NBNatlas", minyear = 2007, maxyear = 2014, covarRes = 300)
occurrence <- raster::as.data.frame(occ)
occurrence <- randomOcc(presframe = occurrence, precisionCol = "precision", covarResm = 300)

test_that("function returns dataframe with easting and northing", {

  expect_is(occurrence, "data.frame")

})

test_that("function returns dataframe with different easting and northing", {

  before <- dplyr::sample_n(occ, 5)
  after <- subset(occurrence, NBNObservationID %in% before$NBNObservationID)
  before <- before[order(before$NBNObservationID),]
  after <- after[order(after$NBNObservationID),]
  expect_false(unique(before$easting == after$easting))
  expect_false(unique(before$northing == after$northing))
})
