#' Place occurence points at random within low-resolution grid cell
#'
#' This function is used to randomise the location of a presence point where the species occurrence data is low resolution. For example, if your species occurrence data is at 1km resolution, but your environmental predictor variables are at 200m resolution, this enables you to randomise where within the 1km cell your species record is placed. The function can handle multiple resolutions within a single dataframe, as long as these are labelled with a precision column.
#'
#' @param presframe Data frame of presence points. This should include columns titled "easting" and "northing", which should be x and y coordinates in metres (of the centre of the grid cell for gridded data), as well as a column giving the resolution (see \code{prescisionCol} below).
#' @param precisionCol The column indicating the resolution of the presence point data. This should be given as a grid cell size in metres e.g. for a 1km grid this should be 1000.
#' @return A copy of the input dataframe \code{presframe}, but with new "easting" and "northing" values that place each point randomly within it's grid cell.
#' @export


randomOcc <- function(presframe, #presence data frame
                      precisionCol = "precision") { #which column defines capture resolution
  occurrence <- data.frame()
  for (j in unique(sort(presframe[[precisionCol]], decreasing = TRUE))) { # Run per capture resolution grid size
    if( j %in% presframe[[precisionCol]] & lowestResm >= j) { #If capture resolution is within lowest resolution to use limit
      df <- presframe[ which(presframe[[precisionCol]] == j), ] # Subset all data at same capture resolution
      if(j > covarResm){ #If capturee resolution is greater than co=variate resolution (random processing not needed if at same or higher resolution)
        # if(unit == "km") {movemax <- 500 * j}
        # if(unit == "m") {movemax <- 0.5 * j}
        movemax <- 0.5 * j #Half of capture resolution grid cell size
        df$easting <- df$easting + round(runif(length(df$easting), -movemax, movemax)) #Move east/west by up to half capture grid cell size
        df$northing <- df$northing + round(runif(length(df$northing), -movemax, movemax)) #Move north/south by up to half capture grid cell size
        occurrence <- rbind(occurrence, df)
      } else {
        occurrence <- rbind(occurrence, df)
      }
    }
  }
  return(occurrence)
}

