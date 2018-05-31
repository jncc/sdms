#' Place occurrence points at random within low-resolution grid cell
#'
#' This function is used to randomise the location of a presence point where the species occurrence data is low resolution. For example, if your species occurrence data is at 1km resolution, but your environmental predictor variables are at 200m resolution, this enables you to randomise where within the 1km cell your species record is placed. The function can handle multiple resolutions within a single dataframe, as long as these are labelled with a precision column.
#'
#' @param presframe Data frame of presence points. This should include columns titled 'easting' and 'northing', which should be x and y coordinates in metres (of the centre of the grid cell for gridded data), as well as a column giving the resolution (see \code{prescisionCol} below).
#' @param precisionCol The column indicating the resolution of the presence point data. This should be given as a grid cell size in metres e.g. for a 1km grid this should be 1000.
#' @param lowestResm The column indicating the lowest resolution limit for the presence point data in meters. If your presence record is in a higher resolution than this limit, then it will remain at its current state.
#' @param covarResm The resolution of the environmental covariate data layers, in metres. Data will not be discarded if it is of higher resolution than the environmental covariate layers.
#' @return A copy of the input dataframe \code{presframe}, with new 'easting' and 'northing' values that place each point randomly within its grid cell.
#' @examples
#'#load in the data
#'data(ng_data)
#'
#'#prepare the data using the bngprep function in this package
#'occurrence <- bngprep(speciesdf = ng_data, bngCol = 'OSGR', datafrom = 'NBNatlas', mindata = 5000, minyear = 2007, covarRes = 300)
#'
#'#run the random occurrences function
#'randomOcc(presframe = occurrence, precisionCol = 'precision')
#' @export


randomOcc <- function(presframe, precisionCol = "precision", lowestResm = 10000,
    covarResm = 1000) {

    occurrence <- data.frame()


    # Check precision column is numerical, and convert to m
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
            } else {
                occurrence <- rbind(occurrence, df)
            }
        }
    }
    return(occurrence)
}
