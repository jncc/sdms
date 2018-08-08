## Tests to check the functionality of the Package function - Occurrence
#  Becky Trippier 08/08/2018
#### ----------------------

context("Occurrence tests - BNG data")
data(sd_data)
names(sd_data)[names(sd_data) == "record _ number"] <- "NBNObservationID"

occ <- bngprep(speciesdf = sd_data, bngCol = "OSGR", datafrom =
                      "NBNatlas", minyear = 1990, maxyear = 2014, covarRes = 300)
occurrence <- raster::as.data.frame(occ)
occurrence <- randomOcc(presframe = occurrence, precisionCol = "precision", covarResm = 10)

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

test_that("precision column is present warning", {
  expect_error(randomOcc(presframe = occurrence, precisionCol = "test", covarResm = 10), "precision column not found")
})

context("Occurrence tests - Lat lon data")
names(sd_data)[26] <- "Decimal longitude (WGS84)"
speciesdf <- latlonprep(speciesdf = sd_data, xCol = "Decimal latitude (WGS84)", yCol = "Decimal longitude (WGS84)", precisionCol = "Coordinate uncertainty in metres", yearCol = "Year", minyear = 2000, maxyear = 2007, GBonly = FALSE)
yCol = "longitude"
xCol = "latitude"


test_that("function returns dataframe with different easting and northing", {

  before <- dplyr::sample_n(speciesdf, 5)
  sp::coordinates(speciesdf)<- c(yCol, xCol)
  occurrence<-randomOcc(presframe = speciesdf, precisionCol = "precision", coordsys = "latlon")
  occurrence <- raster::as.data.frame(occurrence)
  after <- subset(occurrence, NBNObservationID %in% before$NBNObservationID)
  before <- before[order(before$NBNObservationID),]
  after <- after[order(after$NBNObservationID),]
  expect_false(unique(before$longitude == after$longitude))
  expect_false(unique(before$latitude == after$latitude))
})
