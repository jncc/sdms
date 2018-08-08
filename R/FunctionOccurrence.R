#' Place occurrence points at random within low-resolution grid cell
#'
#' This function is used to randomise the location of a presence point where the species occurrence data is low resolution. For example, if your species occurrence data is at 1km resolution, but your environmental predictor variables are at 200m resolution, this enables you to randomise where within the 1km cell your species record is placed. The function can handle multiple resolutions within a single dataframe, as long as these are labelled with a precision column.
#'
#' @param presframe Spatial points data frame of presence points. This should include columns titled 'easting' and 'northing', which should be x and y coordinates in metres (of the centre of the grid cell for gridded data), as well as a column giving the resolution (see \code{prescisionCol} below).
#' @param coordsys The coordinate system used to denote location, either "latlon" for decimal latitude and longitude or "bng" for british national grid easting and northings.
#' @param precisionCol The column indicating the resolution of the presence point data. This should be given as a grid cell size in metres e.g. for a 1km grid this should be 1000.
#' @param lowestResm The column indicating the lowest resolution limit for the presence point data in meters. If your presence record is in a higher resolution than this limit, then it will remain at its current state.
#' @param covarResm The resolution of the environmental covariate data layers, in metres. Data will not be discarded if it is of higher resolution than the environmental covariate layers.
#' @return A copy of the input dataframe \code{presframe} as a spatial points data frame, with new x and y values that place each point randomly within its grid cell.
#' @examples
#'# Example using coordinates in British National Grid :
#'#load in the data
#'data(ng_data)
#'
#'#prepare the data using the bngprep function in this package
#'speciesdf <- bngprep(speciesdf = ng_data, bngCol = 'OSGR', datafrom = 'NBNatlas', mindata = 5000, minyear = 2007, covarRes = 300)
#'
#'#convert to a SpatialPointsDataFrame
#'sp::coordinates(speciesdf)<- ~ easting + northing
#'
#'#run the random occurrences function
#'randomOcc(presframe = speciesdf, precisionCol = "precision", coordsys = "bng")
#'
#'Example using coordinates in latitude and longitude:
#'#load in the data
#'data(ng_data)
#'names(ng_data)[26] <- "Decimal longitude (WGS84)"
#'
#'#prepare the data using the latlonprep function in this package
#'speciesdf <- latlonprep(speciesdf = ng_data, xCol = "Decimal latitude (WGS84)", yCol = "Decimal longitude (WGS84)", precisionCol = "Coordinate uncertainty in metres", yearCol = "Year", minyear = 2000, maxyear = 2007, GBonly = TRUE)
#'
#'#convert to spatial points data frame
#'names(speciesdf)
#'yCol = "longitude"
#'xCol = "latitude"
#'sp::coordinates(speciesdf)<- c(yCol, xCol)
#'#run the random occurrences function
#'randomOcc(presframe = speciesdf, precisionCol = "precision", coordsys = "latlon")
#'
#' @export



randomOcc <- function(presframe, coordsys = "bng",
                      precisionCol = "precision", lowestResm = 10000,
                      covarResm = 300) {


  #check coordinate system used
  if (coordsys != "latlon" & coordsys != "bng")
    stop("invalid coordinate system")

  #deal with lat lon
   if (coordsys == "latlon") {
      ukgrid = "+init=epsg:27700"
      latlong = "+init=epsg:4326"
      proj4string(presframe) <- sp::CRS(latlong)
      presframe <- spTransform(presframe, CRS(ukgrid))
      presframe <- raster::as.data.frame(presframe)
      names(presframe)[which(names(presframe) == "longitude")] <-"easting"
      names(presframe)[which(names(presframe) == "latitude")] <-"northing"
   }


  #convert from spatial to data frame
  occurrence <- data.frame()
  presframe <- raster::as.data.frame(presframe)

   # Check precision column is numerical, and convert to m
  if (!(precisionCol %in% colnames(presframe))){
    stop("precision column not found")
  }
  if (!is.numeric(presframe[[precisionCol]])) {
    presframe$res <- sub("[^[:alpha:]]+", "", presframe[[precisionCol]])
    presframe[[precisionCol]] <- as.numeric(gsub("([0-9]+).*$", "\\1",
                                                 presframe[[precisionCol]]))
    presframe[[precisionCol]] <- ifelse(presframe$res == "m", presframe[[precisionCol]],
                                        presframe[[precisionCol]] * 1000)
    presframe$res <- NULL
  }

    for (j in unique(sort(presframe[[precisionCol]], decreasing = TRUE))) {
      # Run per capture resolution grid size If capture resolution is within
      # lowest resolution to use limit
      if (j %in% presframe[[precisionCol]] & lowestResm >= j) {
        df <- presframe[which(presframe[[precisionCol]] == j), ]  # Subset all data at same capture resolution
        if (j > covarResm) {
          movemax <- 0.5 * j  #Half of capture resolution grid cell size
          df$easting <- df$easting + round(stats::runif(length(df$easting),
                                                        -movemax, movemax))  #Move east/west by up to half capture grid cell size
          df$northing <- df$northing + round(stats::runif(length(df$northing),
                                                        -movemax, movemax))  #Move north/south by up to half capture grid cell size
          occurrence <- rbind(occurrence, df)
          message("jittering applied.")
        } else {
          occurrence <- rbind(occurrence, df)
        }
      } else {
        occurrence <- presframe
        message("no jittering required.")
      }
    }

  if (coordsys == "latlon") {
    names(occurrence)[which(names(occurrence) == "easting")] <-"longitude"
    names(occurrence)[which(names(occurrence) == "northing")] <-"latitude"
  }

  if (coordsys == "bng"){
      sp::coordinates(occurrence) <- ~easting + northing
  } else {
        sp::coordinates(occurrence)<- ~longitude + latitude
        proj4string(occurrence) <- sp::CRS(ukgrid)
        occurrence <- spTransform(occurrence, CRS(latlong))

  }


  return(occurrence)
}
