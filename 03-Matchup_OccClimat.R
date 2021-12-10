############################################################################################
# I - LOOP to match-up occurrences/climat for each month/year reported in the occurrences  #					    
#	1 Subset of occurrences for date j						   #
#	2 Climat data for date j : T.Month.file  && S.Month.file			   #
#	3 Extract Temp & salinity values for the occ subset : matchup climat/occ   	   #
#											   #
############################################################################################										

########################
#  Outputs directory   #
########################

### We create Output directory for Species occurrences/Climat match-up

Matchup_Dir <- file.path(Maindir,"RESULTS/ClimateOccurrences/"); if(!dir.exists(Matchup_Dir)){dir.create(Matchup_Dir,recursive = T)}


### Load match-up file if it already exists

if ( file.exists( file.path(Matchup_Dir,paste0("Occ_Clim_",SP_name,".RData") ) ) ) {
   
  sp_occ_Clim <- get( load( file.path(Matchup_Dir,paste0("Occ_Clim_",SP_name,".RData") )  ) ) 
  cat("-----> ","\t",c("MATCH-UP already DONE" ,"\n") )
  
}else{ 


###############################################################
#  LOOP to match-up occurrences/climat for each month/year    #
###############################################################

	# Occurrences with the same dates of observation

	SpeciesOCC$date <- paste0(SpeciesOCC$year,"-",str_pad(SpeciesOCC$month, 2, side ="left",pad="0") )
	samedate <- unique(SpeciesOCC$date)

	# Empty object

	sp_occ_Clim <- NULL
  
	compteur <- 0
 
	for (j in samedate){
    
	    compteur <- compteur +1 ; cat (compteur, "/", length(samedate),"; \t")
    
	    # 1 Subset of occurrences for date j
	    sub.occ <- subset( SpeciesOCC,date==j)


	    # 2 Climat data for date j : T.Month.file  && S.Month.file

	    Tclim <- file.path (BackRast, grep(dir(BackRast),pattern=paste0("Temp_",j),value=T) )
	    rTemp <-  RasterClim ; values(rTemp) <- get(load(Tclim))
        
	    Sclim <-   file.path (BackRast, grep(dir(BackRast),pattern=paste0("Sal_",j),value=T) )
	    rSal <-    RasterClim ; values(rSal) <- get(load(Sclim))
    
	    # 3 Extract Temp & salinity values for the occ subset : matchup climat/occ
	    sub.occ$longitude <- as.numeric(sub.occ$longitude ) ; sub.occ$latitude <- as.numeric(sub.occ$latitude )
	    coordinates(sub.occ) <- ~longitude + latitude
	    sub.occ.Temp <- raster::extract(rTemp,sub.occ)
	    sub.occ.Sal <- raster::extract(rSal,sub.occ)
    
	    # 4 Add matchup occ to object sp_occ_Clim
	    sp_occ_Clim <- data.frame(rbind(     sp_occ_Clim     , cbind(coordinates(sub.occ),sub.occ@data,"Temp"=sub.occ.Temp,"Sal"=sub.occ.Sal)))
    
	    rm(list=c("sub.occ.Temp","sub.occ.Sal","sub.occ","rSal","rTemp","Tclim","Sclim"))
        
	}#eo j in samedate
  
	rm("samedate","compteur")
	sp_occ_Clim <- na.omit(sp_occ_Clim)

  	### 4 Save math-up occ/climat file
 	if ( !is.null(sp_occ_Clim) && nrow(sp_occ_Clim) > 0 ){

		Name <- paste0("Occ_Clim_",SP_name,".RData")
		save(sp_occ_Clim,file=file.path(Matchup_Dir,Name),version=2)
		rm("Name")

		cat("-----> ","\t",c("MATCH-UP DONE" ,"\n") )

	}else{
		cat("-----> ","\t",c("NO OCCURRENCES after match-up" ,"\n") )
	}


	
  

}#eo else

