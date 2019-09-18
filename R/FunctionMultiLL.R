#' Set up modelling multiple species from NBN atlas using locations given in latitude and longitude for use in species distribution modelling.
#'
#' This function is to carry out modelling of multiple species using the JNCCSdms package.
#'
#'@param sp_list List of unique species names which you wish to model.
#'@param out_flder The location of the output folder for your models.
#'@param dat_flder The location of the folder containing your species occurrence data, as txt or csv files exported from NBN gateway or NBN atlas. Each file should contain data for a single species and the naming convention should correspond to your species list in order to be recognised. e.g. 'Triturus cristatus' in the sp_list should have a corresponding data file named 'Triturus cristatus.csv' in the dat_folder.
#'@param bkgd_flder The location of the folder containing your background masks. These should be raster files showing the background area in which pseudo-absence points will be placed. Cells from which background points should be taken should have a value of 1 and excluded cells should be NA. This should be named after the Taxon Group e.g. 'amphibian' and if this is not found in the data by a 'taxonGroup' variable, then pseudo absences with be generated from the variables layer.
#'@param vars A RasterStack of the environmental parameters to be used as predictor variables for the species range.
#'@param max_tries The number of times the number is run.
#'@param datafrom Whether it is data from the 'NBNgateway' or 'NBNatlas'.
#'@param minyear Numeric, the earliest year from which data should be selected. Year inclusive, data older than this will be discarded.
#'@param maxyear Numeric, the latest year from which data should be used. Year inclusive, data newer than this will be discarded.
#'@param mindata The target minimum number of data points to return. If this is specified, the lowest resolution data will be discarded if there are enough higher resolution data points available to reach this target.
#'@param covarRes The resolution of the environmental covariate data layers, in metres. Data will not be discarded if it is of higher resolution than the environmental covariate layers.
#'@param models A character vector of the models to run and evaluate. This should be at least one of \code{'MaxEnt'}, \code{'BioClim'}, \code{'SVM'}, \code{'RF'}, \code{'GLM'}, \code{'GAM'}, \code{'BRT'}. Default is to run all models.
#'@param prop_test_data Numeric, the proportion of data to keep back as testing data for evaluating the models. Default is 25\%.
#'@param mult_prssr Set up a parallel backend to use multiple processors. As a default this is turned off. Need to ensure the suggested packages have been loaded in order to run this.
#'@param rndm_occ Logical, Default is TRUE and will randomise the locations of presence points where the species occurrence data is low resolution, through calling the randomOcc function.
#'@param GBonly logical, TRUE if you wish to remove Northern Ireland from the records, FALSE if you wish to retain all records.
#'@param xCol The column name of the column in \code{speciesdf} giving the species record location as a decimal latitude.
#'@param yCol The column name of the column in \code{speciesdf} giving the species record location as a decimal longitude.
#'@param precisionCol The column name of the column in \code{speciesdf} denoting the precision of the species record locations. For NBNAtlas this is denoted as the "Coordinate uncertainty (m)". where denoted with km or m, this will be converted into meters.
#'@param yearCol The column name in \code{speciesdf} giving the year of the record.
#'@return Lists containing predictions from the best models (as a raster layer showing probability of species occurrence), the best model evaluations and the best model itself for each species in a given species list.
#'@examples
#'#Provide a list of species you wish to model
#'sp_list <- c("Notonecta_glauca", "Sigara_dorsalis")
#'
#'#Organise an Input folder containing your input species files as .csv
#'dir.create("Inputs")
#'data("ng_data")
#'data("sd_data")
#'names(ng_data)[26] <- "Decimal longitude (WGS84)"
#'names(sd_data)[26] <- "Decimal longitude (WGS84)"
#'utils::write.csv(ng_data, file = "./Inputs/Notonecta_glauca.csv")
#'utils::write.csv(ng_data, file = "./Inputs/Sigara_dorsalis.csv")
#'
#'#Organise a folder containing your background masks where your pseudo absences will be generated from.
#'dir.create("BGmasks")
#'data("background")
#'latlong = "+init=epsg:4326"
#'ukgrid = "+init=epsg:27700"
#'proj4string(background) <- sp::CRS(ukgrid)
#'background = projectRaster(background, crs = latlong)
#'save(background, file = "./BGmasks/Hemiptera")
#'
#'#Create outputs folder
#'dir.create("Outputs")
#'
#'#Get variables data
#'data(vars)
#'proj4string(vars) <- sp::CRS(ukgrid)
#'vars = projectRaster(vars, crs = latlong)
#'
#'#run the function
#'output <- MultiLL(sp_list = sp_list, vars, out_flder = "Outputs/",dat_flder = "Inputs/", bkgd_flder = "BGmasks/", max_tries = 1, covarRes = 100, models = c("BioClim", "RF"), prop_test_data = 0.25, mult_prssr = FALSE, rndm_occ = TRUE, minyear =2000, maxyear = 2007, GBonly = TRUE, xCol = "Decimal latitude (WGS84)", yCol = "Decimal longitude (WGS84)", precisionCol = "Coordinate uncertainty in metres", yearCol = "Year")
#'@export


MultiLL <- function(sp_list = sp_list, out_flder = "Outputs/", dat_flder, bkgd_flder, vars, max_tries = 1,minyear = 0, maxyear = 0, mindata = 5000, covarRes = 300, models = c("MaxEnt", "BioClim", "SVM", "RF", "GLM", "GAM", "BRT"), prop_test_data = 0.25, mult_prssr = FALSE, rndm_occ = TRUE, GBonly = TRUE, xCol = "Latitude (WGS84)", yCol = "Longitude (WGS84)", precisionCol = "Coordinate uncertainty (m)", yearCol = "Year") {

  alreadyDone <- NULL
  sp_found <- 0
  sp_missing <- NULL

  ptm <- proc.time()

  ### input checks
  for (i in 1:length(sp_list)) {
    sp <- sp_list[i]
    if  (file.exists(paste(dat_flder, sp,".csv", sep = ""))) {
      sp_found <- sp_found + 1
    } else if (file.exists(paste(dat_flder, sp, ".txt", sep = ""))) {
      sp_found <- sp_found + 1
    } else {
      sp_missing <- append(sp_missing, sp)
    }
  }
  print(paste(sp_found, "out of", length(sp_list), " species files found in dat_flder."))
  if (length(sp_missing) > 0) {
    print("species files missing:")
    for (i in 1:length(sp_missing)) {
      sp <- sp_missing[i]
      print(sp)
    }
  }

  if (sp_found == 0) {
    stop("No species found. Check input data folder and file formats.")
  }

  #### ------------------ setting up done list -----------------------####
  if (any(duplicated(sp_list)) == TRUE) {
    sp_list <- unique(sp_list)
    message("Duplicate species removed.")
  }

  ## output files already exist, add to alreadyDone list
  for (i in 1:length(sp_list)) {
    sp <- sp_list[i]
    lab <- gsub("([[:punct:]])|\\s+", "_", sp)
    if (file.exists(paste(out_flder, sp, max_tries, ".csv", sep = "")) &
        file.exists(paste(out_flder, sp, max_tries, ".grd", sep = "")) |
        file.exists(paste(out_flder, sp, max_tries, ".tif", sep = ""))) {
      alreadyDone <- append(alreadyDone, sp)
    } else if (file.exists(paste(out_flder, lab, max_tries, ".csv",sep = "")) & file.exists(paste(out_flder, lab, max_tries,".grd", sep = "")) | file.exists(paste(out_flder, lab, max_tries, ".tif", sep = ""))) {
      alreadyDone <- append(alreadyDone, sp)
    }
  }

  sp_list <- sp_list[!sp_list %in% alreadyDone]  #remove alreadyDone species from sp_list

  if (length(unlist(sp_list)) < 1) {
    stop("No more species to process. Modelling terminated.")
  }

  #### ------------------ setting up multiple processors -----------------------####

  ## setup parallel backend to use multiple processors if required
  if (mult_prssr == TRUE) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)  #not to overload your computer
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl, .libPaths("D:/Rpackages"))  #guide each to libpath, as above (line 13)
    packlist <- c("stringi", "rgdal", "raster", "rgeos", "dismo",
                  "kernlab", "dismo", "kernlab", "randomForest", "glmulti",
                  "beepr", "mgcv", "readr", "maptools", "mapdata", "ff", "ffbase",
                  "rnrfa", "pryr")
    message("Parallel clustering established.")
  }
  ptm <- proc.time()


  #### ------------------ run models (non-parallel processing) -----------------------####
  if (mult_prssr == FALSE) {
    for (i in 1:length(sp_list)) {
      # use if running one at a time
      sp <- sp_list[i]
      lab <- gsub("([[:punct:]])|\\s+", "_", sp)

      # read in and prepare species occurrence data
      atlas_dat <- utils::read.csv(file = paste(dat_flder, sp,".csv", sep = ""), header = TRUE, sep = ",", check.names = FALSE,strip.white = TRUE)
      spdat<- latlonprep(atlas_dat, xCol = xCol, yCol = yCol,precisionCol = precisionCol, yearCol = yearCol, minyear = minyear, maxyear = maxyear, mindata = mindata,covarRes = covarRes, GBonly = GBonly)

      # convert to spatial data frame
      sp::coordinates(spdat)<- ~longitude + latitude

      # Create background
      if ("taxonGroup" %in% names(spdat)) {
        taxon <- as.character(spdat$taxonGroup[1])
      } else {
        print("Unable to find taxonGroup field in data.")
      }

      #try to load the background
      load_obj <- function(f)
      {
        env <- new.env()
        nm <- load(f, env)[1]
        env[[nm]]
      }

      bkgd_name <-  paste(bkgd_flder, taxon, sep = "")
      background <- tryCatch(load_obj(bkgd_name), error = function(err) NA)
      if (exists("background")) {
        print("Background mask found.")
      } else {
        r <- vars[[1]]
        r[!is.na(r[])] <- 1
        background <- r
        rm(r)
        print("Unable to obtain background mask from folder. This has been generated from the vars layer.")
      }

      ## run SDM function
      final_out <- SDMs(occ = spdat, varstack = vars, models = models,prop_test_data = prop_test_data, covarReskm = covarRes,max_tries = max_tries, lab = sp, rndm_occ = rndm_occ, out_flder = out_flder, bckg = background,  coordsys = "latlon",precisionCol="precision")
      message(paste(sp, " modelling completed."))
      print("Overall runtime:")
      print(ptm <- proc.time())
      beepr::beep()
    }
  }

  #### ------------------ run models (parallel processing) -----------------------####
  if (mult_prssr == TRUE) {
    foreach(i = 1:length(sp_list), .packages = packlist) %dopar% {
      # use if using parallel cores
      sp <- sp_list[i]
      lab <- gsub("([[:punct:]])|\\s+", "_", sp)

      # read in and prepare species occurrence data
      atlas_dat <- utils::read.csv(file = paste(dat_flder, sp,".csv", sep = ""), header = TRUE, sep = ",", check.names = FALSE,strip.white = TRUE)
      spdat<- latlonprep(atlas_dat, xCol = xCol, yCol = yCol, precisionCol = precisionCol, yearCol = yearCol, minyear = minyear, maxyear = maxyear, mindata = mindata,covarRes = covarRes, GBonly = GBonly)

      # convert to spatial data frame
      sp::coordinates(spdat)<- ~longitude + latitude

      # Create background
      if ("taxonGroup" %in% names(spdat)) {
        taxon <- as.character(spdat$taxonGroup[1])
      } else {
        print("Unable to find taxonGroup field in data.")
      }

      #try to load the background

      load_obj <- function(f)
      {
        env <- new.env()
        nm <- load(f, env)[1]
        env[[nm]]
      }

      bkgd_name <-  paste(bkgd_flder, taxon, sep = "")
      background <- tryCatch(load_obj(bkgd_name), error = function(err) NA)

      if (exists("background")) {
        print("Background mask found.")
      } else {
        r <- vars[[1]]
        r[!is.na(r[])] <- 1
        background <- r
        rm(r)
        print("Unable to obtain background mask from folder. This has been generated from the vars layer.")
      }

      ## run SDM function
      final_out <- SDMs(occ = spdat, varstack = vars, models = models,
                        prop_test_data = prop_test_data, covarReskm = covarRes,
                        max_tries = max_tries, lab = sp, rndm_occ = rndm_occ, bckg = background, out_flder = out_flder,  coordsys = "latlon", precisionCol="precision")
      ptm <- proc.time()

      # stop cluster if parrallel processing
      parallel::stopCluster(cl)
      message("Parallel clustering off.")
      print("Overall runtime:")
      print(proc.time() - ptm)
      beepr::beep()

    }
  }
}

