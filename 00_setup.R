########################
# SETUP : librairies   #
########################

library(parallel)
library(taxize)
library(ncdf4)
library(dplyr)
library(RColorBrewer)
library(rvertnet)
library(spocc)
library(robis)
library(sp)
library(rgdal)
library(maptools)
library(raster)
library(gridExtra)
library(biomod2) # If you are not using a Docker container please install the Github version of the package 
library(PresenceAbsence)
library(rgeos)
library(plyr)
library(ecospat)
library(ade4)
library(viridis)
library(stringr)
library(nctools)
library(mapdata)
library(lubridate)
library(readr)
library(kali)
library(fields)
library(units)
library(R.utils)
library(scam)
library(foreach)
library(doSNOW)
library (snow)
library (doParallel)

#############
# FUNCTIONS #
#############

##### function to find bottom values
find.bottom <- function(x,depths){
  z.max <- max(depths[!(is.na(x))])
  if ( z.max==-Inf  )  {
    return(NA)
  }else{ 
    return(return(x[which(depths==z.max)]))}
}


#####
getNcTime <- function(nc) { 
  
  options(warn=1) #show warnings by default
  if (is.character(nc)) nc <- nc_open(nc)
  ncdims <- names(nc$dim) #get netcdf dimensions
  timevar <- ncdims[which(ncdims %in% c("time", "Time", "datetime", "Datetime", "date", "Date"))] #find (first) time variable
  if (length(timevar) > 1) {
    warning(paste("Found more than one time var. Using the first:", timevar[1]))
    timevar <- timevar[1]
  }
  if (length(timevar)!=1) stop("ERROR! Could not identify the correct time variable")
  times <- ncvar_get(nc, timevar) #get time data
  timeatt <- ncatt_get(nc, timevar) #get attributes
  timeunit <- timeatt$units
  units(times) <- as_units(timeunit)
  as.POSIXct(times)
}


#####
findPrimeMeridian = function(x) {
  if(any(x<0)) return("center")
  if(any(x>180)) return("left")
  warning("Indeterminate Prime Meridian from longitude values.")
  return(NULL)
}



#####
.checkVarid = function(varid, nc) {
 
  if(missing(varid)) varid = NA
  
  if(is.na(varid)) {
    if(length(nc$var)==1) varid = nc$var[[1]]$name
    msg = sprintf("Several variables found in %s, must specify 'varid'.", nc$filename)
    if(length(nc$var)>1) stop(msg)
  }
  
  if(inherits(varid, "ncvar4")) varid = varid$name
  
  if(!is.character(varid))
    stop("varid must be a string or an object of class 'ncvar4'.")
  
  varExists = varid %in% names(nc$var)
  
  msg = sprintf("Variable '%s' not found in '%s'.", varid, nc$filename)
  if(!varExists) stop(msg)
  return(varid)  
}


