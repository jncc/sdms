## Tests to check the functionality of the Package function - Multimod
#  Becky Trippier 05/06/2018
#### ----------------------


context("Multimod tests - general")

### set working directory for tests
start <- getwd()
setwd(tempdir())

## create test folders with data
dir.create("Outputs")
dir.create("Inputs")
dir.create("BGmasks")
data("ng_data")
data("sd_data")
data("background")
utils::write.csv(ng_data, file = "./Inputs/Notonecta_glauca.csv")
utils::write.csv(ng_data, file = "./Inputs/Sigara_dorsalis.csv")
save(background, file = "./BGmasks/Hemiptera")

#### get test data
data(vars)

#### testing

test_that("test if lists containing duplicates are handled", {
  sp_list <- c("Notonecta_glauca", "Sigara_dorsalis", "Notonecta_glauca")
  expect_message(if (any(duplicated(sp_list)) == TRUE) {
    sp_list <- unique(sp_list)
    message("Duplicate species removed.")
  },"Duplicate species removed.")

})

test_that("test if sp_lists contain names not found in input folder", {
  sp_list <- c("fakeosaurus")
  expect_error(Multi_mod(sp_list = sp_list, vars, out_flder = "Outputs/",dat_flder = "Inputs/", bkgd_flder = "BGmasks/", max_tries = 1, datafrom = "NBNatlas", covarRes = 100, models = "BioClim", prop_test_data = 0.25, bngCol = "OSGR", mult_prssr = FALSE, rndm_occ = TRUE),"No species found. Check input data folder and file formats.")

})

test_that("test if iterate through list and species complete messages generated", {
  oldw <- getOption("warn")
  options(warn=-1)

  sp_list <- c("Notonecta_glauca", "Sigara_dorsalis")

   out <- Multi_mod(sp_list = sp_list, vars, out_flder = "Outputs/",dat_flder = "Inputs/", bkgd_flder = "BGmasks/", max_tries = 1, datafrom = "NBNatlas", covarRes = 100, models = "BioClim", prop_test_data = 0.25, bngCol = "OSGR", mult_prssr = FALSE, rndm_occ = TRUE)

   expect_true(file.exists("./Outputs/Notonecta_glauca1.csv") == TRUE)
   expect_true(file.exists("./Outputs/Sigara_dorsalis1.tif") == TRUE)
   options(warn = oldw)
})

test_that("test if already done list is updated", {
  sp_list <- c("Notonecta_glauca", "Sigara_dorsalis")

   expect_error(Multi_mod(sp_list = sp_list, vars, out_flder = "Outputs/",dat_flder = "Inputs/", bkgd_flder = "BGmasks/", max_tries = 1, datafrom = "NBNatlas", covarRes = 100, models = "BioClim", prop_test_data = 0.25, bngCol = "OSGR", mult_prssr = FALSE, rndm_occ = TRUE),"No more species to process. Modelling terminated.")


})



# remove temporary files and return to working directory
unlink("Outputs", recursive=TRUE)
setwd(start)
unlink(tempdir(), recursive=TRUE)
