## Tests to check the functionality of the Package function - Multimod
#  Becky Trippier 22/05/2018
#### ----------------------


load("D:/Github/sdms/R/sysdata.rda")
load(file = "D:/Becky docs/SDMs/varsLCM2015")


#get species list
sp_list <- c("Notonecta_glauca", "Sigara_dorsalis")
 out <- Multi_mod(sp_list = sp_list, out_flder = "Outputs/", dat_flder = "Inputs/", bkgd_flder = "BGmasks/", vars, max_tries = 1, datafrom = "NBNatlas", minyear = 0, maxyear = 0, mindata = 5000, covarRes = 300, models = "BioClim", prop_test_data = 0.25, bngCol = "OSGR", mult_prssr = FALSE, rndm_occ = TRUE)

context("Multimod tests - general")

test_that("test if lists containing duplicates are handled", {
  sp_list <- c("sp_gatewaydata", "sp_atlasdata", "sp_atlasdata")

  expect_message(if (any(duplicated(sp_list)) == TRUE) {
    sp_list <- unique(sp_list)
    message("Duplicate species removed.")
  },"Duplicate species removed.")
})

