#' Prepare data from NBN Atlas for use in species distribution modelling.
#'
#' This function removes absence data (as the models generate their own pseudo-absences), subset the data based on year and allows modelling of GB and Northern Ireland using latitude and longitude coordinates.
#'
#'@param speciesdf Data frame exported from the NBN atlas with data for a species.
#'@param xCol The column name of the column in \code{speciesdf} giving the species record location as a decimal latitude.
#'@param yCol The column name of the column in \code{speciesdf} giving the species record location as a decimal longitude.
#'@param precisionCol The column name of the column in \code{speciesdf} denoting the precision of the species record locations. For NBNAtlas this is denoted as the "Coordinate uncertainty (m)". where denoted with km or m, this will be converted into meters.
#'@param yearCol The column name in \code{speciesdf} giving the year of the record.
#'@param minyear Numeric, the earliest year from which data should be selected. Year inclusive, data older than this will be discarded.
#'@param maxyear Numeric, the latest year from which data should be used. Year inclusive, data newer than this will be discarded.
#'@param mindata The target minimum number of data points to return. If this is specified, the lowest resolution data will be discarded if there are enough higher resolution data points available to reach this target.
#'@param covarRes The resolution of the environmental covariate data layers, in metres. Data will not be discarded if it is of higher resolution than the environmental covariate layers.
#'@param GBonly logical, TRUE if you wish to remove Northern Ireland from the records, FALSE if you wish to retain all records.

#'@return A copy of \code{speciesdf} as a data frame with absence data removed, subset to your year range and area of interest.
#'@examples
#'data(ng_data)
#'names(ng_data)[26] <- "Decimal longitude (WGS84)"
#'
#'latlonprep(speciesdf = ng_data, xCol = "Decimal latitude (WGS84)", yCol = "Decimal longitude (WGS84)", precisionCol = "Coordinate uncertainty in metres", yearCol = "Year", minyear = 2000, maxyear = 2007, GBonly = TRUE)
#'
#'@export

latlonprep <- function(speciesdf, xCol = "Latitude (WGS84)", yCol = "Longitude (WGS84)",
                       precisionCol = "Coordinate uncertainty (m)",
                       yearCol = "Year", minyear = 0, maxyear = 0, mindata = 500,
                       covarRes = 300, GBonly = TRUE) {



  #---------------------Subset data------------------------------------------------------------#
  ### incorrect input errors
  if (anyDuplicated(names(speciesdf))>0){
    names(speciesdf) <- make.unique(names(speciesdf))
    message("Duplicated field names made unique.")
  }
  if (maxyear != 0 & maxyear < minyear)
    stop("minimum year limit is greater than maximum year")

  if (!(xCol %in% colnames(speciesdf)) | !(yCol %in% colnames(speciesdf))){
    stop("check x and y column names")
    }
  if (!(yearCol %in% colnames(speciesdf))){
    stop("year column  not found")
  }
  if (!(precisionCol %in% colnames(speciesdf))){
    stop("precision column  not found")
  }

  ##--------------------------##

  ## remove any levels so it is just flat data to work with
  speciesdf <- droplevels(speciesdf)


  ## Data from the NBNatlas Inital clean to give only presence records,
  ## also subset records by bngCol as varying precisions for Grid
  ## References

  speciesdf <- speciesdf[speciesdf$`Occurrence status` == "present",
                         ]  # remove absence data

  speciesdf <- speciesdf[grepl("[[:alnum:]]", speciesdf[[xCol]]),
                         ]  #select rows with gridref records
  speciesdf <- speciesdf[grepl("[[:alnum:]]", speciesdf[[yCol]]),
                         ]  #select rows with gridref records

   if (nrow(speciesdf) == 0){
    stop("no presence records found")
  }


#rename yearCol to year
  names(speciesdf)[names(speciesdf) == yearCol] <- "year"

  # Extract by date - if this has been defined
  if (minyear > 0 & maxyear > 0) {
    speciesdf <- speciesdf[which(speciesdf$year >= minyear & speciesdf$year <=
                                   maxyear), ]  # if minyear & maxyear defined
  } else if (minyear > 0 & maxyear == 0) {
    speciesdf <- speciesdf[which(speciesdf$year >= minyear), ]  # if only minyear defined
  } else if (minyear == 0 & maxyear > 0) {
    speciesdf <- speciesdf[which(speciesdf$year <= maxyear), ]  # if only maxyear defined
  } else {
    speciesdf <- speciesdf
  }

  if (nrow(speciesdf) == 0){
    stop("no presence records found within year range")
  }

  # Check precision column is numerical, and convert to m
  if (!is.numeric(speciesdf[[precisionCol]])) {
    speciesdf$res <- sub("[^[:alpha:]]+", "", speciesdf[[precisionCol]])
    speciesdf[[precisionCol]] <- as.numeric(gsub("([0-9]+).*$", "\\1",
                                                 speciesdf[[precisionCol]]))
    speciesdf[[precisionCol]] <- ifelse(speciesdf$res == "m", speciesdf[[precisionCol]],
                                        speciesdf[[precisionCol]] * 1000)
    speciesdf$res <- NULL
    message("precisionCol converted to m")
  } else {
    message("precisionCol numerical")
  }

# find GB and Ireland
  speciesdf$GB  <- NA


  raster::rasterOptions(tmpdir = "./Rtmpdir")

  GB <- raster::getData('GADM', country="gbr", level=2)
  GB_sub <- raster::subset(GB, NAME_1 != "Northern Ireland")
  IR_sub <- raster::subset(GB, NAME_1 == "Northern Ireland")

  GBextent <- raster::extent(GB_sub)
  IRextent <- raster::extent(IR_sub)


  #highlight GB or Ireland

  IRxmin <- IRextent@ymin
  IRxmax <- IRextent@ymax
  IRymin <- IRextent@xmin
  IRymax <- IRextent@xmax

  for (i in 1:nrow(speciesdf)) {
    if (speciesdf[[xCol]][i] >=IRxmin & speciesdf[[xCol]][i] <= IRxmax & speciesdf[[yCol]][i] >=IRymin & speciesdf[[yCol]][i] <= IRymax) {
              speciesdf$GB[i] <- "IR"
    } else {
      speciesdf$GB[i] <- "GB"
    }
  }
  unlink("./Rtmpdir/*")



  #subset to GB
  if (GBonly == TRUE){
    speciesdf <- speciesdf[speciesdf$GB  == "GB",]  # subset to GB only
    message(paste(nrow(speciesdf), "occurrences after subsetting"))

  } else {
    message(paste(nrow(speciesdf), "occurrences after subsetting"))
  }

  # Remove low resolution data if sufficient data at higher resolution
  # If a minimum number of data points has been supplied
  if (!is.na(mindata)) {
    rescount <- table(speciesdf[[precisionCol]])
    startpoint <- nrow(speciesdf)
    for (j in sort(as.numeric(names(rescount)), decreasing = TRUE)) {
      # sort data by precision, with lowest resolution first
      if (j == covarRes)
      {
        break
      }  #No point refining data to a higher resolution than analysis resolution
      if (sum(rescount[!names(rescount) == as.character(j)], na.rm = TRUE) >
          mindata) {
        # If the sum of data excluding lowest resolution is greater than the
        # minimum data set,
        speciesdf <- speciesdf[which(!speciesdf[[precisionCol]] ==
                                       j), ]  #cut lowest resolution data
        rescount[[as.character(j)]] <- NA
      } else {
        break
      }
    }
    message(paste(startpoint - nrow(speciesdf)), " records with low resolution data points removed"
)
  }

  names(speciesdf)[which(names(speciesdf) == precisionCol)] <-"precision"
  names(speciesdf)[which(names(speciesdf) == xCol)] <-"latitude"
  names(speciesdf)[which(names(speciesdf) == yCol)] <-"longitude"

return(speciesdf)

}
