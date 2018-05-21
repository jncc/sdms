#' Prepare presence and pseudo-absence data, run selected models and evaluate best model.
#'
#' This function uses presence points for a species with a background mask to generate pseudo-absences, and then uses these with environmental data layers to generate species distribution models using a range of different model algorithms. The function then selects the best perfoming model and outputs the distribution that that model predicts, an evaluation of the model performance and the model itself.
#'
#' @param occ A SpatialPointsDataFrame of presence points.
#' @param bckg A raster showing the background area in which pseudo-absence points will be placed. Cells from which background points should be taken should have a value of 1 and excluded cells should be NA.
#' @param varstack A RasterStack of the environmental parameters to be used as predictor variables for the species range.
#' @param models A character vector of the models to run and evaluate. This should be at least one of \code{'MaxEnt'}, \code{'BioClim'}, \code{'SVM'}, \code{'RF'}, \code{'GLM'}, \code{'GAM'}, \code{'BRT'}. Default is to run all models.
#' @param n_bg_points The number of pseudo-absence point to attempt to generate. Note that if a very restrictive mask is used the number actually generated may be fewer than that specified. Default is to attempt to generate the same number of pseudo-absences as presences for which there is data on the environmental parameters (this may be fewer than the number of points in \code{occ} if some of these fall in cells that are NA in one or more layers in \code{varstack}.
#' @param prop_test_data Numeric, the proportion of data to keep back as testing data for evaluating the models. Default is 25\%.
#' @param covarReskm Numeric, the resolution of the environmental covariate data layers, in metres. Data will not be discarded if it is of higher resolution than the environmental covariate layers.
#' @param max_tries Numeric, the numbers of times the model is run.
#' @param lab The name of the output files.
#' @param rndm_occ Logical, Default is TRUE and will randomise the locations of presence points where the species occurrence data is low resolution, through calling the randomOcc function.
#' @return A list containing the prediction from the best model (as a raster layer showing probability of species occurrence), the best model evaluation and the best model itself.
#' @examples
#' SDMs(occ = occurrence, bckg = background, varstack = variables, max_tries = 5, lab = 'species', rndm_occ = TRUE)
#' @export

SDMs <- function(occ = occurrence, bckg = background, varstack = vars, models = c("MaxEnt",
    "BioClim", "SVM", "RF", "GLM", "GAM", "BRT"), n_bg_points = nrow(pres_vars), prop_test_data = 0.25,
    covarReskm = 300, max_tries = 2, lab = "species", rndm_occ = TRUE) {

    all_predicts <- NULL
    all_models <- NULL
    all_evals <- NULL
    tries <- 0

    ptm <- proc.time()
    best_out <- c()

    while (tries < max_tries) {
        raster::rasterOptions(tmpdir = "./Rtmpdir")


        #---------------------generate random occurrences -------------------------------------------------#

        if (rndm_occ == TRUE) {
            occurrence <- raster::as.data.frame(occ)
            occurrence <- randomOcc(presframe = occurrence, precisionCol = "precision")
            pres.pts <- data.frame(occurrence[, c("easting", "northing")])
            sp::coordinates(occurrence) <- ~easting + northing
            message("Random occurrences placed.")
        } else {
            pres.pts <- as.data.frame(sp::coordinates(occ))
        }

        #---------------------generate pseudo-absences from mask-------------------------------------------------#

        ## Get background points (from within background mask if it exists, but excluding
        ## grids with presence) bgNonPres <- raster::mask(bckg, occ, inverse=TRUE) bg.pts <-
        ## raster::sampleRandom(bgNonPres, max_bg_points, xy = TRUE, sp=TRUE, na.rm = TRUE)
        ## ## Due to the small depth range of the mask, this function returned dramatically
        ## fewer points than max_bg_points (~10%). I've therefore substituted for
        ## randomPoints from the Dismo package below, and given it up to 50 attempts to
        ## generate the number of points requested. Also moved this to after the filtering
        ## of pres_vars for complete.cases, so that the number of background points can
        ## match the number of presence points. bg.coords <- sp::coordinates(bg.pts)

        ## Extract variable values for presence and pseudo-absence (background) points
        ppts <- as.data.frame(sp::coordinates(pres.pts))
        rownames(ppts) <- 1:nrow(ppts)
        colnames(ppts) <- c("x", "y")
        pres_vars <- data.frame(cbind(unique(ppts), raster::extract(varstack, unique(ppts))))  #length (which (is.na (pres_vars$bio1)))
        pres_vars <- pres_vars[stats::complete.cases(pres_vars), ]  #Remove any NA lines
        pres_vars$Presence <- 1
        bg.pts <- dismo::randomPoints(mask = bckg, n = n_bg_points, p = ppts, tryf = 50)
        bg_vars <- data.frame(cbind(bg.pts, raster::extract(vars, bg.pts)))
        bg_vars <- bg_vars[stats::complete.cases(bg_vars), ]  #Remove any NA lines
        bg_vars$Presence <- 0
        # Combine presence and pseudo absence training data
        combo_vars <- rbind(pres_vars, bg_vars)


        #---------------------training and test data -------------------------------------------------#


        # Split presence data into training and test datasets
        presence_random <- sample(nrow(pres_vars), ceiling(prop_test_data * (nrow(pres_vars))))
        presence_train <- pres_vars[-presence_random, ]
        presence_test <- pres_vars[presence_random, ]
        # Split pseudo-absence data into training and test datasets
        bg_random <- sample(nrow(bg_vars), ceiling(prop_test_data * (nrow(bg_vars))))
        bg_train <- bg_vars[-bg_random, ]
        bg_test <- bg_vars[bg_random, ]
        # Combine presence and pseudo absence training data
        combo_train <- rbind(presence_train, bg_train)
        # Drop variables with fewer than 12 unique values (minimum required for GAM)
        lowunique <- NULL
        for (i in 1:(raster::ncol(combo_train) - 1)) {
            p <- length(unique(combo_train[, i]))
            if (p < 12) {
                lowunique <- append(lowunique, names(combo_train[i]))
            }
        }
        if (length(lowunique) >= 1) {
            combo_vars <- combo_vars[, !(colnames(combo_vars) %in% lowunique)]  #'[[<-'(combo_vars, lowunique, value = NULL)
            combo_train <- combo_train[, !(colnames(combo_train) %in% lowunique)]  #'[[<-'(combo_train, lowunique, value = NULL)
            vars_drop <- raster::dropLayer(varstack, lowunique)
        } else {
            vars_drop <- varstack
        }
        varfactors <- names(combo_vars[, 3:(ncol(combo_vars) - 1)])

        #---------------------Define models and evaluate using test data---------------------------------------------------#

        # Define models and evaluate using test data Evaluate does not work on some models
        # (notably some GAMs) as output is array, not matrix. eval function is a reduced
        # version of dismo::evaluate which calculates auc for array outputs , x, tr
        eval <- function(p, a, model) {
            p <- kernlab::predict(model, data.frame(p))
            a <- kernlab::predict(model, data.frame(a))
            np <- length(p)
            na <- length(a)
            N <- na + np
            xc <- rJava::new("ModelEvaluation")
            xc@presence = p
            xc@absence = a
            R <- sum(rank(c(p, a))[1:np]) - (np * (np + 1)/2)
            xc@auc <- R/(as.numeric(na) * as.numeric(np))
            return(xc)
        }

        models <- sort(models)
        aucs <- NULL
        # MaxEnt - machine learning presence-only
        if ("MaxEnt" %in% models) {
            Maxent <- tryCatch({
                dismo::maxent(vars_drop, presence_train[, c("x", "y")], args = c("nolinear",
                  "noquadratic", "noproduct", "nothreshold"), path = "/out")
            }, error = function(err) NA)
            eval_test <- tryCatch(dismo::evaluate(presence_test[, c("x", "y")], bg_test[,
                c("x", "y")], Maxent, vars_drop), error = function(err) {
                tryCatch(eval(presence_test[, c("x", "y")], bg_test[, c("x", "y")],
                  Maxent, vars_drop), error = function(err) NA)
            })
            aucs <- list(aucs, MaxEnt = tryCatch(eval_test@auc, error = function(err) NA))
            message("MaxEnt model run ", tries + 1, " completed.")
        }

        # Bioclim - profile method presence-only
        if ("BioClim" %in% models) {
            BioClim <- tryCatch(dismo::bioclim(vars_drop, presence_train[, c("x", "y")]),
                error = function(err) NA)
            eval_test <- tryCatch(dismo::evaluate(presence_test[, c("x", "y")], bg_test[,
                c("x", "y")], BioClim, vars_drop), error = function(err) {
                tryCatch(eval(presence_test[, c("x", "y")], bg_test[, c("x", "y")],
                  BioClim, vars_drop), error = function(err) NA)
            })
            aucs <- list(aucs, BioClim = tryCatch(eval_test@auc, error = function(err) NA))
            message("Bioclim model run ", tries + 1, " completed.")
        }


        # Kernal SVM - machine learning presence-absence
        if ("SVM" %in% models) {
            SVM <- tryCatch(kernlab::ksvm(stats::as.formula(paste("Presence ~", paste(varfactors,
                collapse = "+"))), data = combo_train[, 3:ncol(combo_train)], na.action = stats::na.omit,
                kpar = "automatic"), error = function(err) NA)
            eval_test <- tryCatch(dismo::evaluate(presence_test[, 3:ncol(presence_test)],
                bg_test[, 3:ncol(bg_test)], SVM), error = function(err) {
                tryCatch(eval(presence_test[, 3:ncol(presence_test)], bg_test[, 3:ncol(bg_test)],
                  SVM), error = function(err) NA)
            })
            aucs <- list(aucs, SVM = tryCatch(eval_test@auc, error = function(err) NA))
            message("Kernal SVM model run ", tries + 1, " completed.")
        }

        # Random Forest - machine learning presence-absence (could modify to export
        # importance and proximity summary)
        if ("RF" %in% models) {
            mtry <- tryCatch(randomForest::tuneRF(x = combo_train[, 3:(ncol(combo_train) -
                1)], y = combo_train[, ncol(combo_train)], trace = FALSE, plot = FALSE),
                error = function(err) NA)
            if (!is.na(mtry)) {
                best.m <- tryCatch(mtry[mtry[, 2] == min(mtry[, 2]), 1], error = function(err) NA)
                RF <- tryCatch(randomForest::randomForest(stats::as.formula(paste("Presence ~",
                  paste(varfactors, collapse = "+"))), data = combo_train[, 3:ncol(combo_train)],
                  na.action = stats::na.omit, mtry = best.m), error = function(err) NA)
            } else {
                RF <- tryCatch(randomForest::randomForest(stats::as.formula(paste("Presence ~",
                  paste(varfactors, collapse = "+"))), data = combo_train[, 3:ncol(combo_train)],
                  na.action = stats::na.omit), error = function(err) NA)
            }
            eval_test <- tryCatch(dismo::evaluate(presence_test[, 3:ncol(presence_test)],
                bg_test[, 3:ncol(bg_test)], RF), error = function(err) {
                tryCatch(eval(presence_test[, 3:ncol(presence_test)], bg_test[, 3:ncol(bg_test)],
                  RF), error = function(err) NA)
            })
            aucs <- list(aucs, RF = tryCatch(eval_test@auc, error = function(err) NA))
            message("Random forest model run ", tries + 1, " completed.")
        }

        # GLM - regression presence-absence (attempts glm selectionfirst, but if too many
        # predictors does full model glm)
        if ("GLM" %in% models) {
            glm_multi <- tryCatch(do.call(glmulti::glmulti, list(y = "Presence", xr = varfactors,
                family = "binomial", data = combo_train[, 3:ncol(combo_train)], na.action = stats::na.omit,
                level = 1, method = "g", confsetsize = 1, deltaB = 0, conseq = 1, maxit = 200,
                report = FALSE, plotty = FALSE)), error = function(err) NA)
            if (class(glm_multi) == "numeric") {
                GLM <- tryCatch(stats::glm(stats::as.formula(paste("Presence ~", paste(varfactors,
                  collapse = "+"))), family = "binomial", data = combo_train[, 3:ncol(combo_train)],
                  na.action = stats::na.omit), error = function(err) NA)
            } else {
                gl_vector_string <- tryCatch(summary(glm_multi)$bestmodel, error = function(err) NA)
                gl_string <- tryCatch(paste(gl_vector_string, sep = "", collapse = ""),
                  error = function(err) NA)
                GLM <- tryCatch(stats::glm(gl_string, family = "binomial", data = combo_train[,
                  3:ncol(combo_train)], na.action = stats::na.omit), error = function(err) NA)
            }
            eval_test <- tryCatch(dismo::evaluate(presence_test[, 3:ncol(presence_test)],
                bg_test[, 3:ncol(bg_test)], GLM), error = function(err) {
                tryCatch(eval(presence_test[, 3:ncol(presence_test)], bg_test[, 3:ncol(bg_test)],
                  GLM), error = function(err) NA)
            })
            aucs <- list(aucs, GLM = tryCatch(eval_test@auc, error = function(err) NA))
            message("GLM model run ", tries + 1, " completed.")
        }

        # GAM
        if ("GAM" %in% models) {
            GAM <- tryCatch(mgcv::gam(stats::as.formula(paste("Presence ~ s(", paste(varfactors,
                collapse = ")+ s("), ")")), data = combo_train[, 3:ncol(combo_train)],
                family = stats::binomial(link = "logit"), gamma = 1.4, na.rm = TRUE),
                error = function(err) NA)
            eval_test <- tryCatch(dismo::evaluate(presence_test[, 3:ncol(presence_test)],
                bg_test[, 3:ncol(bg_test)], GAM), error = function(err) {
                tryCatch(eval(presence_test[, 3:ncol(presence_test)], bg_test[, 3:ncol(bg_test)],
                  GAM), error = function(err) NA)
            })
            aucs <- list(aucs, GAM = tryCatch(eval_test@auc, error = function(err) NA))
            message("GAM model run ", tries + 1, " completed.")
        }

        # Boosted Regression Trees - machine learning presence-absence (could modify to
        # export importance and proximity summary)
        if ("BRT" %in% models) {
            # Use tryCatch as does not work for all datasets
            BRT <- tryCatch(dismo::gbm.step(data = combo_train, gbm.x = 3:(ncol(combo_train) -
                1), gbm.y = as.numeric(ncol(combo_train)), tree.complexity = 3, learning.rate = 0.01,
                bag.fraction = 0.75, verbose = FALSE, silent = TRUE, plot.main = FALSE),
                error = function(err) NA)
            eval_test <- tryCatch(dismo::evaluate(presence_test[, 3:ncol(presence_test)],
                bg_test[, 3:ncol(bg_test)], BRT, n.trees = BRT$n.trees), error = function(err) {
                tryCatch(eval(presence_test[, 3:ncol(presence_test)], bg_test[, 3:ncol(bg_test)],
                  BRT), error = function(err) NA)
            })
            aucs <- list(aucs, BRT = tryCatch(eval_test@auc, error = function(err) NA))
            message("Boosted Regression model run ", tries + 1, " completed.")
        }

        aucs <- unlist(aucs)
        aucs <- aucs[order(names(aucs))]
        best <- ifelse(max(aucs, na.rm = TRUE) < 0.7, "Presence", names(aucs)[which.max(aucs)])
        # aucs2 <- list(auc_mx_test, auc_bc_test, auc_svm_test, auc_rf_test, auc_glm_test,
        # auc_gam_test, auc_brt_test, best)

        # Add function to return only presence observations as raster if model evaluation
        # is low
        Obs2ras <- function(points, raskm = covarReskm) {
            # nr <- 1250 / raskm nc <- 700 / raskm r <- raster::raster(nrows = nr, ncols = nc,
            # xmn = 0, xmx= 700000, ymn = 0, ymx = 1250000)
            r <- vars[[1]]
            r[!is.na(r[])] <- 0
            i <- raster::extract(r, points, cellnumbers = TRUE)[, "cells"]
            r[i] <- 1
            return(r)
        }

        # Project to entire region (GB)
        bestmod <- tryCatch(get(best), error = function(err) NA)
        prediction_best <- tryCatch(kernlab::predict(vars_drop, bestmod, type = "response"),
            error = function(err) Obs2ras(ppts[c("x", "y")]))

        # Convert raster output to ff matrix to reduce object size
        predict_ff <- ff::ff(raster::as.matrix(prediction_best))
        # Add to list of previous runs
        all_predicts <- append(all_predicts, list(predict_ff))
        # Add to list of best full models
        all_models <- append(all_models, list(bestmod))
        # Collate model evaluations for info
        all_evals <- list(all_evals, as.list(c(aucs, best)))
        close(predict_ff)
        unlink("D:/Rdata/Rtmpdir/*")
        gc()

        out <- list(all_predicts = all_predicts, all_evals = all_evals, all_models = all_models)

        best_out <- append(best_out, out)
        list2env(best_out, .GlobalEnv)
        tries <- tries + 1
        message("Run ", tries, " completed.")
        print(c(tries, all_evals[[2]][[8]]))

    }
    #---------------------Generate model evaluation results-----------------------------#

    # Export model evaluation results
    Mean_predict <- Reduce(`+`, all_predicts)/length(all_predicts)
    Mean_predict <- matrix(unlist(Mean_predict[, ]), nrow = nrow(vars[[1]]), ncol = ncol(vars[[1]]),
        byrow = FALSE)
    Mean_predict <- raster::setValues(vars[[1]], Mean_predict)
    Models <- c(sort(models), "Best")
    all_evals <- cbind(data.frame(Models), data.frame(matrix(unlist(all_evals), nrow = length(models) +
        1, byrow = F), stringsAsFactors = FALSE))
    raster::plot(Mean_predict)
    # beep(0)
    proc.time() - ptm
    pryr::mem_used()

    save(all_models, file = paste("Outputs/", lab, tries, "models", sep = ""))
    raster::writeRaster(Mean_predict, file = paste("Outputs/", lab, tries, sep = ""),
        format = "GTiff")
    utils::write.csv(all_evals, file = paste("Outputs/", lab, tries, ".csv", sep = ""))

    beepr::beep()

}
