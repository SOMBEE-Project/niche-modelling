############################################################################################
# I - LOOP to create background raster for each month				 	   #					    
#	1 Load Global 3D Clim data path 						   #
#	2 Mean/surface/bottom raster 							   #
#	3 Monthly 2D Backgound 							   	   #
#								   			   #
# II - MERGE monthly Background files							   #
# 											   #
############################################################################################										

###############################
#      Load GLOBAL LIMITS     #
###############################

# min & max temperature & salinity throughout the time series and all over the world  (limits of env. back)
load(file.path(Maindir,"RESULTS/Backgrounds/S_globlimits.RData"))
S_min_glob  <- S_globlimits[,"S_min_glob"]
S_max_glob <- S_globlimits[,"S_max_glob"]

load(file.path(Maindir,"RESULTS/Backgrounds/T_globlimits.RData") )
T_min_glob <- T_globlimits[,"T_min_glob"]
T_max_glob <- T_globlimits[,"T_max_glob"]


########################
#  Outputs directories #
########################

# Backgrounds direcotory : 
BackDir <- file.path(Maindir,"RESULTS/Backgrounds/") ; if(!dir.exists(BackDir)){dir.create(BackDir,recursive = TRUE)}

# Backgrounds sub-direcotories : 
BackFile <- file.path(BackDir,"bfiles",sp_layers) ; if(!dir.exists(BackFile)){dir.create(BackFile,recursive = TRUE)}
BackRast <- file.path(BackDir,"rasters",sp_layers) ; if(!dir.exists(BackRast)){dir.create(BackRast,recursive = TRUE)}


### Load background for the species min/max depth layers if it already exists (count & binary) 

BCoutput = file.path(BackDir, paste0( "Back_count_",sp_layers,".RData"   ) )
BBoutput = file.path(BackDir, paste0( "Back_binary_",sp_layers,".RData"   ) ) 


if ( file.exists(BCoutput) && file.exists(BBoutput) ) {

  Sp_back <- get(load(BBoutput))  
  Sp_back.r <- raster (resolution=c(0.1,0.1), xmn = T_min_glob, xmx = T_max_glob,ymn = S_min_glob, ymx = S_max_glob)  # empty raster res 0.1
  values(Sp_back.r) <- Sp_back 
  
  Sp_back_count <- get(load(BCoutput))
  Sp_back_count.r <- raster (resolution=c(0.1,0.1), xmn = T_min_glob, xmx = T_max_glob,ymn = S_min_glob, ymx = S_max_glob) # empty raster res 0.1
  values(Sp_back_count.r) <- Sp_back_count ; rm("count_back_all")

  cat("-----> ","\t",c("BACKGROUND & Global rasters already DONE" ,"\n") )

}else{


##########################################################
#  LOOP to create background raster for each month       #
##########################################################


# Monthly files for Mercator data : between 93 & 2018

dates = paste0 (  rep(seq(1993,2018,1), each = 12) ,   "-" ,rep(str_pad(seq(1,12,1),2,side = "left" ,"0"), length(rep(seq(1993,2018,1)))  ))


for (date in dates) {    
  
  cat(date ,";" , "\t")
  
  ### 1 Load Global 3D Clim data path : T.Month.file  && S.Month.file
  Tclim <-  Tfiles[str_detect(Tfiles,date)]
  rTemp <- brick(Tclim,lvar=4)
  Toutput = file.path(BackRast,paste0( "meanTemp_",date,".RData" ) )    

  Sclim <-  Sfiles[str_detect(Sfiles,date)]
  rSal <-   brick(Sclim,lvar=4)
  Soutput = file.path(BackRast,paste0( "meanSal_",date,".RData" ) )  


  ### 2 Mean/surface/bottom raster 

  if ( !file.exists(Toutput) | !file.exists(Soutput) ) {

	# if bottom => bottom values with function find.bottom
		if ( sp_layers == "bottom"  ) {

			rTempA <- as.array(rTemp)
			Tbottomvalues <- apply(rTempA, c(1,2), find.bottom,depths=Obs.depths)
			Tmean <- subset(rTemp, 1) ; values(Tmean) <- Tbottomvalues
			Tmean <- values(Tmean)                        
    
			rSalA <- as.array(rSal)
			Sbottomvalues <- apply(rSalA, c(1,2), find.bottom,depths=Obs.depths)
			Smean <- subset(rSal, 1) ; values(Smean) <- Sbottomvalues
			Smean <- values(Smean)

			rm ("rTemp","rTempA","rSal","rSalA","Sbottomvalues","Tbottomvalues")
		}# eo bottom
	

	# if surface => first layer
		if ( sp_layers == "surface"  ) {

		    layers = 1	   

		    Tlayers <- subset(rTemp,layers)
		    Tmean <- values(Tlayers)            
		    
		    Slayers <- subset(rSal,layers)
		    Smean <- values(Slayers) 

		    rm ("layers","Tlayers","Slayers")

	# else mean values between min/max sp depth

		}else{      

		    layers <- seq(as.numeric (unlist (str_split(sp_layers,"_") )[1]),as.numeric (unlist (str_split(sp_layers,"_") )[2]),1)
			   
		    Tlayers <- subset(rTemp,layers)
		    Tmean <- mean(Tlayers,na.rm=TRUE)
		    Tmean <- values(Tmean)                        
			    
		    Slayers <- subset(rSal,layers)
		    Smean <- mean(Slayers,na.rm=TRUE)
		    Smean <- values(Smean) 

		    rm ("layers","Tlayers","Slayers")

		} # eo ifelse

 
   ### Save rasters ###
   save( Tmean,file =Toutput , version=2 ) 
   save( Smean,file =Soutput , version=2) 
   rm ("Tmean","Smean") 

  }# if raster files exist


  ### 3 Monthly Backgound
 
  Boutput  = file.path(BackFile,paste0( "back_",date,".RData" ) )  

  if ( !file.exists(Boutput) ) {

	Tmean <- get(load(Toutput))
	Smean <- get(load(Soutput))
	  
	combTS <- data.frame( cbind("Temperature" =Tmean,"Salinity" = Smean ) )
	combTS <- na.omit(combTS)
	combTS <- round(combTS,2)
	  
	Back.r <- raster(resolution=c(0.1,0.1), xmn = T_min_glob,xmx = T_max_glob,
		         ymn = S_min_glob,ymx = S_max_glob)
	  
	combTS.r <- rasterize(combTS, Back.r, fun="count")
	  
	### 4 Save rasters ###
	back_month <- values(combTS.r)
	save(back_month,file = Boutput )

	rm ("Tmean","Smean","combTS.r","back_month")

	removeTmpFiles(h=0.2)
	gc()
	  
  }#eo Boutput


} # eo for dates


#######################################
#                                     #
# MERGE monthly Background files      #
#                                     #
#######################################

files <- dir(BackFile)  # All monthly env background rasters

back_all <- NULL  # empty object : data frame with values for all monthly env background rasters

for (file in files){
  BG <- get(load(file.path(BackFile,file)))  
  back_all <- cbind (back_all,BG)
}#eo for files


count_back_all <- apply(back_all,1,sum,na.rm=TRUE)
count_back_all[ which(count_back_all == 0 )] <- NA
Sp_back_count.r <- raster (resolution=c(0.1,0.1), xmn = T_min_glob, xmx = T_max_glob,ymn = S_min_glob, ymx = S_max_glob)
values(Sp_back_count.r) <- count_back_all

save(count_back_all,file= BCoutput)
rm("count_back_all")

count_back_bin <- ifelse(count_back_all > 0 , 1, NA )
Sp_back.r <- raster (resolution=c(0.1,0.1), xmn = T_min_glob, xmx = T_max_glob,ymn = S_min_glob, ymx = S_max_glob)  
values(Sp_back.r) <- count_back_bin

save(count_back_bin,file= BBoutput)
rm("count_back_bin")


cat("-----> ","\t",c("BACKGROUND & Global rasters DONE" ,"\n") )

}#eo else
