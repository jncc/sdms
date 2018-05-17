#' Prepare data from NBN gateway or atlas for use in species distribution modelling.
#'
#' This function removes absence data (as the models generate their own pseudo-absences), converts British National Grid grid references into easting and northing coordinates and removes data from Northern Ireland to allow modelling using GB data layers.
#'
#'@param speciesdf Data frame exported from NBN gateway or NBN atlas with data for a species.
#'@param bngCol The column name of the column in \code{speciesdf} giving the species record location as a BNG grid reference. For NBNatlas, this will vary with recorder precision, so you should select the most appropriate for your data, for example "OSGR 1km", "OSGR 2km" or "OSGR 10km".
#'@param precisionCol The column name of the column in \code{speciesdf} denoting the precision of the species record locations. This only needs to be specified for NBNgateway data as this is generated from the selected \code{bngCol} for data from the NBNatlas.
#'@param datafrom Character, one of "NBNgateway" or "NBNatlas", indicating the data source.
#'@param minyear Numeric, the earliest year from which data should be selected. Year inclusive, data older than this will be discarded.
#'@param maxyear Numeric, the latest year from which data should be used. Year inclusive, data newer than this will be discarded.
#'@param mindata The target minimum number of data points to return. If this is specified, the lowest resolution data will be discarded if there are enough higher resolution data points available to reach this target.
#'@param covarRes The resolution of the environmental covariate data layers, in metres. Data will not be discarded if it is of higher resolution than the environmental covariate layers.
#'@return A copy of \code{speciesdf} with absence data removed, easting and northing columns generated from BNG grid references and data from Northern Ireland removed.
#'@examples
#'#Example using data from NBN Gateway:
#'bngprep(speciesdf = sp_gatewaydata, precisionCol = "precision", bngCol = "gridReference", datafrom = "NBNgateway", mindata = 5000, minyear = 2006, maxyear = 2016, covarRes = 300)
#'
#'#Examples using data from NBN Atlas:
#'bngprep(speciesdf = sp_atlasdata, bngCol = "OSGR 1km", datafrom = "NBNatlas", mindata = 5000, minyear = 2007, covarRes = 300)
#'
#'bngprep(speciesdf = sp_atlasdata, bngCol = "OSGR 2km", datafrom = "NBNatlas", mindata = 5000, minyear = 1990, maxyear = 2010, covarRes = 500)
#'
#'bngprep(speciesdf = sp_atlasdata, bngCol = "OSGR 10km", datafrom = "NBNatlas", mindata = 5000, maxyear = 2002, covarRes = 300)
#'@export

bngprep <- function(speciesdf,
                    bngCol,
                    precisionCol = "precision",
                    datafrom = NULL,
                    minyear = 0,
                    maxyear = 0,
                    mindata = 5000,
                    covarRes) {


  #---------------------Subset data------------------------------------------------------------#

  ## Data from the NBNgateway
  #  Inital clean to give presence records in GB
  if (datafrom == "NBNgateway") {
    speciesdf <- speciesdf[ which(speciesdf$zeroAbundance == "FALSE"),] # remove absence data
    speciesdf <-speciesdf[ which(speciesdf$Projection == "OSGB36"),] # subset to GB only
    speciesdf$year <- as.numeric( format( as.Date(speciesdf$startDate, format="%d/%m/%Y"), '%Y')) # convert date format and add year column
  }

  # Extract by date - if this has been defined
  if (datafrom == "NBNgateway") {
    if(minyear > 0 & maxyear > 0) {
      speciesdf <- speciesdf[ which(speciesdf$year >= minyear & speciesdf$year <= maxyear), ] # if minyear & maxyear defined
    } else if (minyear > 0 & maxyear == 0) {
      speciesdf <- speciesdf[ which(speciesdf$year >= minyear), ] # if only minyear defined
    } else if (minyear == 0 & maxyear > 0) {
      speciesdf <- speciesdf[ which(speciesdf$year <= maxyear), ] # if only maxyear defined
    } else { speciesdf <- speciesdf}
  }



  ## Data from the NBNatlas
  #  Inital clean to give presence records in GB, also subset records by bngCol as varying precisions for Grid References
  if (datafrom == "NBNatlas") {
    speciesdf <- speciesdf[speciesdf$`Occurrence status` == "present", ] # remove absence data
    speciesdf <- speciesdf[!speciesdf$`State/Province` == "Northern Ireland",] # subset to GB only
    speciesdf <- speciesdf[ grepl("[[:alnum:]]",speciesdf[[bngCol]]),] #select rows with gridref records
    speciesdf$precision <- unlist(regmatches(bngCol, gregexpr("[[:digit:]]{1,3}km", bngCol))) #add precision


    # ensure grid references are consistant
  if (speciesdf$precision[1] == "1km") {
    speciesdf <- speciesdf[ grepl("^[a-zA-Z]{2}[0-9]{4}$",speciesdf[[bngCol]]),]
  } else if (speciesdf$precision[1] == "2km") {
    speciesdf <- speciesdf[ grepl("^[a-zA-Z]{2}[0-9]{2}[a-zA-Z]{1}$",speciesdf[[bngCol]]),]
  } else if (speciesdf$precision[1] == "10km") {
    speciesdf <- speciesdf[ grepl("^[a-zA-Z]{2}[0-9]{2}$",speciesdf[[bngCol]]),]
  } else {
    speciesdf <- speciesdf[ grepl("^[a-zA-Z]{2}$",speciesdf[[bngCol]]),]
    }
    }




  # Extract by date - if this has been defined
  if (datafrom == "NBNatlas") {
    if(minyear > 0 & maxyear > 0) {
      speciesdf <- speciesdf[ which(speciesdf$Year >= minyear & speciesdf$Year <= maxyear), ] # if minyear & maxyear defined
    } else if (minyear > 0 & maxyear == 0) {
      speciesdf <- speciesdf[ which(speciesdf$Year >= minyear), ] # if only minyear defined
    } else if (minyear == 0 & maxyear > 0) {
      speciesdf <- speciesdf[ which(speciesdf$Year <= maxyear), ] # if only maxyear defined
    } else {speciesdf <- speciesdf}
  }


  #---------------------Run module------------------------------------------------------------#

  #Check precision column is numerical, and convert to m
  if(!is.numeric(speciesdf[[precisionCol]]) ) {
    speciesdf$res <- sub("[^[:alpha:]]+", "", speciesdf[[precisionCol]])
    speciesdf[[precisionCol]] <- as.numeric(gsub("([0-9]+).*$", "\\1", speciesdf[[precisionCol]]))
    speciesdf[[precisionCol]] <- ifelse(speciesdf$res == "m", speciesdf[[precisionCol]], speciesdf[[precisionCol]]*1000)
    speciesdf$res <- NULL
  }

  # Split between normal grid refs and 'tetrad' (2km) grids
  nontetrad <- speciesdf[ which(!speciesdf[[precisionCol]] == 2000),]
  tetrad <- speciesdf[ which(speciesdf[[precisionCol]] == 2000),]

  if(nrow(nontetrad) > 0) {
    for (i in 1:nrow(nontetrad) ) {
      ne <- rnrfa::osg_parse(nontetrad[i,bngCol])
      nontetrad$easting[i] <- ne[[1]]
      nontetrad$northing[i] <- ne[[2]]
    }
  }
  speciesdf <- nontetrad

  #Calculate tetrad grids seperately
  if(nrow(tetrad) > 0){
    tetrad$At10kcorner <- paste(stringi::stri_sub(tetrad[[bngCol]], 1, 4))
    tetrad$Letter <- paste(stringi::stri_sub(tetrad[[bngCol]], 5, -1))
    for (i in 1:nrow(tetrad) ) {
      ne <- rnrfa::osg_parse(tetrad$At10kcorner[i])
      tetrad$easting[i] <- ne[[1]]
      tetrad$northing[i] <- ne[[2]]
    }
    Tetrad <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","P","Q","R","S","T","U","V","W","X","Y","Z")
    AddEast <- c(1000,1000,1000,1000,1000,3000,3000,3000,3000,3000,5000,5000,5000,5000,5000,7000,7000,7000,7000,7000,9000,9000,9000,9000,9000)
    AddNorth <- c(1000,3000,5000,7000,9000,1000,3000,5000,7000,9000,1000,3000,5000,7000,9000,1000,3000,5000,7000,9000,1000,3000,5000,7000,9000)
    AddTetrad <-data.frame(Tetrad, AddEast, AddNorth)
    tetrad <- merge(tetrad, AddTetrad, by.x="Letter", by.y="Tetrad", all.x=T)
    tetrad$easting <- tetrad$easting + tetrad$AddEast
    tetrad$northing <- tetrad$northing + tetrad$AddNorth
    tetrad$At10kcorner <- tetrad$AddEast <- tetrad$AddNorth <- tetrad$Letter <- NULL
    speciesdf <- rbind(nontetrad,tetrad)

  }

  #Remove low resolution data if sufficient data at higher resolution
  if (!is.na(mindata)) { #If a minimum number of data points has been supplied
    rescount <-table(speciesdf[[precisionCol]])
    for (j in sort(as.numeric(names(rescount)), decreasing = TRUE)) { #sort data by precision, with lowest resolution first
      if(j == covarRes) { break } #No point refining data to a higher resolution than analysis resolution
      if(sum(rescount[!names(rescount)== as.character(j)], na.rm = TRUE) > mindata) { #If the sum of data excluding lowest resolution is greater than the minimum data set,
        speciesdf <- speciesdf[ which(!speciesdf[[precisionCol]] == j),] #cut lowest resolution data
        rescount[[as.character(j)]]<-NA
      } else {
        break
      }
    }
  }

  return(speciesdf)
}
