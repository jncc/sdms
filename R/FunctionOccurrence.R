## Biodiversity Module
## Function to prepare create occurrence data from presence points by adding random jitter
## Add random amount to eastings/nrothings to place coordinates
## randomly within resolution grid cell
## If capture resolution is higher than covariable resolution, 
## just take grid cell coords, no gain from using random point locality
## Originally based on script from BiodiversityMod_v8.R

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

