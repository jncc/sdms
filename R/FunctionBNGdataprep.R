#' Prepare data from NBN gateway or atlas for use in species distribution modelling.
#'
#' This function removes absence data (as the models generate their own pseudo-absences), converts British National Grid grid references into easting and northing coordinates and removes data from Northern Ireland to allow modelling using GB data layers.
#'
#'@param speciesdf Data frame exported from NBN gateway or NBN atlas with data for a species.
#'@param bngCol The column name of the column in \code{speciesdf} giving the species record location as a BNG grid reference. For NBNgateway this is the 'gridReference' field. The default is set to the 'OSGR' field which is used for data from NBNatlas and is also used to generate the location precision.
#'@param precisionCol The column name of the column in \code{speciesdf} denoting the precision of the species record locations. This only needs to be specified for NBNgateway data as this is generated from the selected \code{bngCol} for data from the NBNatlas.
#'@param datafrom Character, one of 'NBNgateway' or 'NBNatlas', indicating the data source.'NBNatlas' is the default.
#'@param minyear Numeric, the earliest year from which data should be selected. Year inclusive, data older than this will be discarded.
#'@param maxyear Numeric, the latest year from which data should be used. Year inclusive, data newer than this will be discarded.
#'@param mindata The target minimum number of data points to return. If this is specified, the lowest resolution data will be discarded if there are enough higher resolution data points available to reach this target.
#'@param covarRes The resolution of the environmental covariate data layers, in metres. Data will not be discarded if it is of higher resolution than the environmental covariate layers.
#'@return A copy of \code{speciesdf} with absence data removed, easting and northing columns generated from BNG grid references and data from Northern Ireland removed.
#'@examples
#'#Examples using data from NBN Atlas:
#'
#'data(ng_data)
#'
#'bngprep(speciesdf = ng_data, bngCol = 'OSGR', datafrom = 'NBNatlas', mindata = 5000, minyear = 2007, covarRes = 300)

#'@export

bngprep <- function(speciesdf, bngCol = "OSGR", precisionCol = "precision", datafrom = "NBNatlas",
    minyear = 0, maxyear = 0, mindata = 5000, covarRes = 300) {


    #---------------------Subset data------------------------------------------------------------#
    ### incorrect input errors
   if (anyDuplicated(names(speciesdf))>0){
     names(speciesdf) <- make.unique(names(speciesdf))
     message("Duplicated field names made unique.")
   }
    if (maxyear != 0 & maxyear < minyear)
        stop("minimum year limit is greater than maximum year")
    if (datafrom != "NBNatlas" & datafrom != "NBNgateway")
        warning("datafrom not specified as NBNatlas or NBNgateway.")
    if (datafrom == "NBNatlas" & precisionCol != "precision") {
        precisionCol <- "precision"
        message("unnecessary argument - do not specify precisionCol for data from NBNatlas")
    }

  if (!(bngCol %in% colnames(speciesdf))){
    stop("bngCol not found.")
  }

  ##----##

  ## remove any levels so it is just flat data to work with
  speciesdf <- droplevels(speciesdf)


    ## Data from the NBNgateway Inital clean to give presence records in GB
    if (datafrom == "NBNgateway") {
        speciesdf <- speciesdf[which(speciesdf$zeroAbundance == "FALSE"),
            ]  # remove absence data
        speciesdf <- speciesdf[which(speciesdf$Projection == "OSGB36"),
            ]  # subset to GB only
        speciesdf$year <- as.numeric(format(as.Date(speciesdf$startDate,
            format = "%d/%m/%Y"), "%Y"))  # convert date format and add year column
    }

    ## Data from the NBNatlas Inital clean to give presence records in GB,
    ## also subset records by bngCol as varying precisions for Grid
    ## References
    if (datafrom == "NBNatlas") {
        speciesdf <- speciesdf[speciesdf$`Occurrence status` == "present",
            ]  # remove absence data
        speciesdf <- speciesdf[!speciesdf$`State/Province` == "Northern Ireland",
            ]  # subset to GB only
        speciesdf <- speciesdf[grepl("[[:alnum:]]", speciesdf[[bngCol]]),
            ]  #select rows with gridref records
        names(speciesdf)[names(speciesdf) == "Year"] <- "year"
        speciesdf$precision <- NA  #add precision

        # ensure grid references are consistant
    speciesdf$precision <- ifelse(grepl("^[a-zA-Z]{2}[0-9]{4}$",speciesdf[[bngCol]]),'1km','NA')
    speciesdf$precision[grepl("^[a-zA-Z]{2}[0-9]{2}[a-zA-Z]{1}$",speciesdf[[bngCol]])] <- '2km'
    speciesdf$precision[grepl("^[a-zA-Z]{2}[0-9]{2}$",speciesdf[[bngCol]])] <- '10km'
    speciesdf$precision[grepl("^[a-zA-Z]{2}$",speciesdf[[bngCol]])] <- '100km'
    speciesdf$precision[grepl("^[a-zA-Z]{2}[0-9]{6}$",speciesdf[[bngCol]])] <- '100m'
    }

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

    message(paste(nrow(speciesdf), "occurrences after subsetting"))

    #---------------------Run module------------------------------------------------------------#

    # Check precision column is numerical, and convert to m
    if (!is.numeric(speciesdf[[precisionCol]])) {
        speciesdf$res <- sub("[^[:alpha:]]+", "", speciesdf[[precisionCol]])
        speciesdf[[precisionCol]] <- as.numeric(gsub("([0-9]+).*$", "\\1",
            speciesdf[[precisionCol]]))
        speciesdf[[precisionCol]] <- ifelse(speciesdf$res == "m", speciesdf[[precisionCol]],
            speciesdf[[precisionCol]] * 1000)
        speciesdf$res <- NULL
    }
    speciesdf$easting<-NA
    speciesdf$northing<-NA
    speciesdf$OSGR <- toupper(speciesdf$OSGR)

    # Split between normal grid refs and 'tetrad' (2km) grids
    nontetrad <- speciesdf[which(!speciesdf[[precisionCol]] == 2000), ]
    tetrad <- speciesdf[which(speciesdf[[precisionCol]] == 2000), ]



    if (nrow(nontetrad) > 0) {
        for (i in 1:nrow(nontetrad)) {
            ne <- rnrfa::osg_parse(nontetrad[i, bngCol])
            nontetrad$easting[i] <- ne[[1]]
            nontetrad$northing[i] <- ne[[2]]
        }
    }
    speciesdf <- nontetrad


    # Calculate tetrad grids seperately
    if (nrow(tetrad) > 0) {
        tetrad$At10kcorner <- paste(stringi::stri_sub(tetrad[[bngCol]],
            1, 4))
        tetrad$Letter <- paste(stringi::stri_sub(tetrad[[bngCol]], 5,
            -1))
        for (i in 1:nrow(tetrad)) {
            ne <- rnrfa::osg_parse(tetrad$At10kcorner[i])
            tetrad$easting[i] <- ne[[1]]
            tetrad$northing[i] <- ne[[2]]
        }
        Tetrad <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
            "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W",
            "X", "Y", "Z")
        AddEast <- c(1000, 1000, 1000, 1000, 1000, 3000, 3000, 3000, 3000,
            3000, 5000, 5000, 5000, 5000, 5000, 7000, 7000, 7000, 7000,
            7000, 9000, 9000, 9000, 9000, 9000)
        AddNorth <- c(1000, 3000, 5000, 7000, 9000, 1000, 3000, 5000,
            7000, 9000, 1000, 3000, 5000, 7000, 9000, 1000, 3000, 5000,
            7000, 9000, 1000, 3000, 5000, 7000, 9000)
        AddTetrad <- data.frame(Tetrad, AddEast, AddNorth)
        tetrad <- merge(tetrad, AddTetrad, by.x = "Letter", by.y = "Tetrad",
            all.x = T)
        tetrad$easting <- tetrad$easting + tetrad$AddEast
        tetrad$northing <- tetrad$northing + tetrad$AddNorth
        tetrad$At10kcorner <- tetrad$AddEast <- tetrad$AddNorth <- tetrad$Letter <- NULL
        colnames(tetrad)[1] <- colnames(nontetrad)[1]
        speciesdf <- rbind(nontetrad, tetrad)

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
    }
    message(paste(startpoint - nrow(speciesdf)), " records with low resolution data points removed")
    return(speciesdf)
}
