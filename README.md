# JNCCsdms
Functions to Generate Species Distribution Models using the JNCC SDM Framework

This package includes functions for ensemble species distribution modelling. The package can be used to run a suite of models on your data and select the best model based on AUC. The package also facilitates randomising the placement of pseudo-absences and the location of presence points if these are drawn from a lower-resolution grid than the environmental predictor data sets. This enables an ensemble modelling approach by repeatedly running the suite of models with the random elements varying each time, and then taking an average of the best model from each run to generate a final output.
