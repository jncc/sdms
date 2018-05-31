#' Set up modelling multiple species from NBN gateway or atlas for use in species distribution modelling.
#'
#' This function is to carry out modelling of multiple species using the JNCCSdms package.
#'
#'@param sp_list List of unique species names which you wish to model.
#'@param out_flder The location of the output folder for your models.
#'@param dat_flder The location of the folder containing your species occurrence data, as data frames exported from NBN gateway or NBN atlas with files containing data for a single species.
#'@param bkgd_flder The location of the folder containing your background masks. These should be raster files showing the background area in which pseudo-absence points will be placed. Cells from which background points should be taken should have a value of 1 and excluded cells should be NA.
#'@param vars A RasterStack of the environmental parameters to be used as predictor variables for the species range.
#'@param max_tries The number of times the number is run.
#'@param datafrom Whether it is data from the 'NBNgateway' or 'NBNatlas'.
#'@param minyear Numeric, the earliest year from which data should be selected. Year inclusive, data older than this will be discarded.
#'@param maxyear Numeric, the latest year from which data should be used. Year inclusive, data newer than this will be discarded.
#'@param mindata The target minimum number of data points to return. If this is specified, the lowest resolution data will be discarded if there are enough higher resolution data points available to reach this target.
#'@param covarRes The resolution of the environmental covariate data layers, in metres. Data will not be discarded if it is of higher resolution than the environmental covariate layers.
#'@param models A character vector of the models to run and evaluate. This should be at least one of \code{'MaxEnt'}, \code{'BioClim'}, \code{'SVM'}, \code{'RF'}, \code{'GLM'}, \code{'GAM'}, \code{'BRT'}. Default is to run all models.
#'@param prop_test_data Numeric, the proportion of data to keep back as testing data for evaluating the models. Default is 25\%.
#'@param bngCol The column name of the column in \code{speciesdf} giving the species record location as a BNG grid reference. For NBNatlas, this will vary with recorder precision, so you should select the most appropriate for your data, for example 'OSGR 1km', 'OSGR 2km' or 'OSGR 10km'.
#'@param mult_prssr Set up a parallel backend to use multiple processors. As a default this is turned off. Need to ensure the suggested packages have been loaded in order to run this.
#'@param rndm_occ Logical, Default is TRUE and will randomise the locations of presence points where the species occurrence data is low resolution, through calling the randomOcc function.
#'@return A copy
#'@examples
sp_list <- c("Notonecta_glauca", "Sigara_dorsalis")

#'@export


Multi_mod <- function(sp_list = sp_list, out_flder = "Outputs/", dat_flder = "Inputs/",
    bkgd_flder = "BGmasks/", vars, max_tries = 1, datafrom = "NBNgatweay",
    minyear = 0, maxyear = 0, mindata = 5000, covarRes = 300, models = c("MaxEnt",
        "BioClim", "SVM", "RF", "GLM", "GAM", "BRT"), prop_test_data = 0.25,
    bngCol = "gridReference", mult_prssr = FALSE, rndm_occ = TRUE) {



    alreadyDone <- NULL

    ## ------------------ setting up done list -----------------------## If
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
        } else if (file.exists(paste(out_flder, lab, max_tries, ".csv",
            sep = "")) & file.exists(paste(out_flder, lab, max_tries,
            ".grd", sep = "")) | file.exists(paste(out_flder, lab, max_tries,
            ".tif", sep = ""))) {
            alreadyDone <- append(alreadyDone, sp)
        }
    }

    sp_list <- sp_list[!sp_list %in% alreadyDone]  #remove alreadyDone species from sp_list


    ## ------------------ setting up multiple processors
    ## -----------------------##

    ## setup parallel backend to use multiple processors if required
    if (mult_prssr == TRUE) {
        cl <- parallel::makeCluster(detectCores() - 1)  #not to overload your computer
        doParallel::registerDoParallel(cl)
        parallel::clusterEvalQ(cl, .libPaths("D:/Rpackages"))  #guide each to libpath, as above (line 13)
        packlist <- c("stringi", "rgdal", "raster", "rgeos", "dismo",
            "kernlab", "dismo", "kernlab", "randomForest", "glmulti",
            "beepr", "mgcv", "readr", "maptools", "mapdata", "ff", "ffbase",
            "rnrfa", "pryr")
        message("Parallel clustering established.")
    }
    ptm <- proc.time()


    ## ------------------ run models (non-parallel processing)
    ## -----------------------##
    if (mult_prssr == FALSE) {
        for (i in 1:length(sp_list)) {
            # use if running one at a time
            sp <- sp_list[i]
            lab <- gsub("([[:punct:]])|\\s+", "_", sp)

            # read in and prepare species occurrence data
            if (datafrom == "NBNgatway") {
                gateway_dat <- readr::read_delim(file = paste(dat_flder,
                  sp, ".txt", sep = ""), "\t", escape_double = FALSE,
                  trim_ws = TRUE)
                spdat <- bngprep(speciesdf = gateway_dat, precisionCol = "precision",
                  bngCol = bngCol, mindata = mindata, minyear = minyear,
                  maxyear = maxyear, covarRes = coverRes)
            } else if (datafrom == "NBNatlas") {
                atlas_dat <- utils::read.csv(file = paste(dat_flder, sp,
                  ".csv", sep = ""), header = TRUE, sep = ",", check.names = FALSE,
                  strip.white = TRUE)
                spdat <- bngprep(speciesdf = atlas_dat, datafrom = datafrom,
                  mindata = mindata, minyear = minyear, maxyear = maxyear,
                  covarRes = covarRes, bngCol = bngCol)
            }
            # convert to spatial data frame
            sp::coordinates(spdat) <- ~easting + northing

            # Create background
            if ("taxonGroup" %in% names(spdat)) {
                taxon <- spdat$taxonGroup[1]
            } else {
                taxon <- readline(paste("unable to find background mask, please type in file name and extension.  "))  #manual prompt to insert taxon group
            }
            tryCatch(load(file = paste(bkgd_flder, taxon, sep = "")),
                error = function(err) NA)
            if (exists("mask1km")) {
                background <- mask1km
                rm(mask1km)
            } else {
                r <- vars[[1]]
                r[!is.na(r[])] <- 1
                background <- r
                rm(r)
            }


            ## run SDM function
            final_out <- SDMs(occ = spdat, varstack = vars, models = models,
                prop_test_data = prop_test_data, covarReskm = covarRes,
                max_tries = max_tries, lab = sp, rndm_occ = rndm_occ, out_flder = out_flder)
            message(paste(sp, " modelling completed."))
            ptm <- proc.time()
            beepr::beep()

        }
    }

    ## ------------------ run models (parallel processing)
    ## -----------------------##
    if (mult_prssr == TRUE) {
        foreach(i = 1:length(sp_list), .packages = packlist) %dopar% {
            # use if using parallel cores

            sp <- sp_list[i]
            lab <- gsub("([[:punct:]])|\\s+", "_", sp)

            # read in and prepare species occurrence data
            if (datafrom == "NBNgatway") {
                gateway_dat <- readr::read_delim(file = paste(dat_flder,
                  sp, ".txt", sep = ""), "\t", escape_double = FALSE,
                  trim_ws = TRUE)
                spdat <- bngprep(speciesdf = gateway_dat, precisionCol = "precision",
                  bngCol = bngCol, mindata = mindata, minyear = minyear,
                  maxyear = maxyear, covarRes = coverRes)
            } else if (datafrom == "NBNatlas") {
                atlas_dat <- utils::read.csv(file = paste(dat_flder, sp,
                  ".csv", sep = ""), header = TRUE, sep = ",", check.names = FALSE,
                  strip.white = TRUE)
                spdat <- bngprep(speciesdf = atlas_dat, datafrom = datafrom,
                  mindata = mindata, minyear = minyear, maxyear = maxyear,
                  covarRes = covarRes, bngCol = bngCol)
            }
            # convert to spatial data frame
            sp::coordinates(spdat) <- ~easting + northing

            # Create background
            if ("taxonGroup" %in% colnames(spdat)) {
                taxon <- spdat$taxonGroup[1]
            } else {
                taxon <- readline(paste("what taxon group is", sp, " ?"))  #manual prompt to insert taxon group
            }
            tryCatch(load(file = paste(bkgd_flder, taxon, sep = "")),
                error = function(err) NA)
            if (exists("mask1km")) {
                background <- mask1km
                rm(mask1km)
            } else {
                r <- vars[[1]]
                r[!is.na(r[])] <- 1
                background <- r
                rm(r)
            }


            ## run SDM function
            final_out <- SDMs(occ = spdat, varstack = vars, models = models,
                prop_test_data = prop_test_data, covarReskm = covarRes,
                max_tries = max_tries, lab = sp, rndm_occ = rndm_occ)
            ptm <- proc.time()

            # stop cluster if parrallel processing
            parallel::stopCluster(cl)
            message("Parallel clustering off.")
            proc.time() - ptm
            beepr::beep()

        }
    }
}



