#' Produce a correlation matrix of all potential input variables
#'
#' @param vars A raster stack of potential input variables
#' @param na.rm `logical` Should cells with at least one NA be ignored?
#' @param order Ordering method for the correlation matrix plot. See `?corrplot::corrplot` for more details
#' @param title Title for the correlation matrix plot
#'
#' @return A matrix of correlation coefficients between layers, and a correlation matrix plot of these coefficients
#' @export
#'
#' @examples
#'
#' #get UK extent
#' UK <- ggplot2::map_data(map = "world", region = "UK")
#' max.lat <- ceiling(max(UK$lat))
#' min.lat <- floor(min(UK$lat))
#' max.lon <- ceiling(max(UK$long))
#' min.lon <- floor(min(UK$long))
#' extent <- raster::extent(x = c(min.lon, max.lon, min.lat, max.lat))
#'
#' #get variables data
#' bio<-raster::getData('worldclim',var='bio',res=5,lon=-2,lat=40)
#'
#' #crop to uk
#' bio<-raster::crop(bio,extent)
#'
#' #produce correlation matrix
#' corr <- corrVars(bio)
#'
#' #END
#'

corrVars <- function(vars, na.rm = T, method = 'color', order = 'hclust', title = 'Layer Correlation Matrix'){

  requireNamespace('raster', quietly = T)
  requireNamespace('corrplot', quietly = T)

  if(raster::nlayers(vars)<2){
    stop("Argument 'vars' must be a raster stack of 2 or more variables for a correlation matrix to be produced")
  }

  if(!is.logical(na.rm)){
    stop("Argument 'na.rm' must be logical")
  }

  if(!(method %in% c("circle", "square", "ellipse", "number", "shade", "color", "pie"))){
    stop("Argument 'method' must be one of 'circle', 'square', 'ellipse', 'number', 'shade', 'color', 'pie'. See ?corrplot::corrplot for more details")
  }

  if(!(order %in% c("original", "AOE", "FPC", "hclust", "alphabet"))){
    stop("Argument 'order' must be one of 'original', 'AOE', 'FPC', 'hclust', 'alphabet'. See ?corrplot::corrplot for more details")
  }

  temp <- raster::layerStats(vars, stat = 'pearson', na.rm = na.rm)

  plt <- corrplot::corrplot(corr = temp$`pearson correlation coefficient`, type = 'lower', method = 'color',
                            diag = F, outline = T, order = order, tl.col = 'black',
                            title = title, mar=c(0,0,3,0), is.corr = F)

  return(plt)
}
