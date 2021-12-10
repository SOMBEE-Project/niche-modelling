####################################################################################################################
# I - LOOP to create Projection rasters / Month & Year / 3D.layer			  			   # #														   #	
#	1 Create sub directories to save Rasters						   		   #
#	2 MONTHLY RASTERS						   					   #
#		2.1 Open Global 3D Climat data path : T.Month.file  && S.Month.file				   #
#		2.2 Crop Global 3D climat files with local domain extent  (path : "RESULTS/ProjRasters/domain")    #
#		2.3 Average across layers 3D local -> 2D local						   	   #
#			2.3.1 subset layers between min and max of pre-fixed "layers3D"				   #
#			2.3.2 mean raster 						   			   #
#			2.3.3 save data.frame (coords, R values)						   #
#	3 ANNUAL RASTERS Mean of Monthly Temperature files for year y						   #
####################################################################################################################
				
########################
#  Outputs directory   #
########################

# Directory for global scale prediction rasters
ProjRastDir <- file.path(Maindir, "RESULTS/ProjectionRasters") ; if(!dir.exists(ProjRastDir)){dir.create(ProjRastDir,recursive = TRUE)}


# sub direcotry to save rasters cropped for study area
DomProjRast <-   file.path(ProjRastDir,sp.dom) ; if (!dir.exists (DomProjRast)) dir.create(DomProjRast,recursive = TRUE)


#######################################################################
#   LOOP to create Projection rasters / Month & Year / 3D.layer	      #
# from Yf = 1993 to Yt = 2018 (see LAUNCH_xxx.R)        	      #
# 2 * 12 * 25 = 624 file / 3D.layer			      	      #
#######################################################################

# Monthly & Annual Climate raster files : Mercator data : between 93 & 2018
dates = paste0 (  rep(seq(1993,2018,1), each = 12) ,   "-" ,rep(str_pad(seq(1,12,1),2,side = "left" ,"0"), length(rep(seq(1993,2018,1)))  ))
years <- seq(Yf,Yt,1)

# Define bounds of study region (minLong-1°,maxLong+1°) ; (minLat-1°,maxLat+1°) # Add 1° from each extent side
lon = (DOM.extent)[c(1,2)]
lat = (DOM.extent)[c(3,4)]
dx = 1 ; xlon = lon + c(-1, +1)*dx   
dy = 1 ; xlat = lat + c(-1, +1)*dy  


for ( layer in layers3D ){

	cat(layer ," ; ")
	regexp <- "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?"
	min_deplayer = as.numeric(str_extract_all(layer, regexp,simplify = T)[1,1])
	max_deplayer = as.numeric(str_extract_all(layer, regexp,simplify = T)[1,2])

	### 1 Create sub directories to save Rasters
		# sub direcotry to save historical climate rasters for depth layers in layers3D
		DomProjRastLayer <-  file.path(DomProjRast,paste0(min_deplayer,"-",max_deplayer) ) ; if (!dir.exists (DomProjRastLayer)) dir.create(DomProjRastLayer,recursive = TRUE)
		# sub direcotry to save historical climate rasters cropped for study area with all 55 depth layers
		DomProjRastZall = file.path(DomProjRast,"Zall") ; if (!dir.exists (DomProjRastZall)) dir.create(DomProjRastZall,recursive = TRUE)

	### 2 MONTHLY RASTERS
	# IF all MONTHLY projection rasters done
	Toutputs = file.path(DomProjRastLayer,paste0("Mercator_Temperature_" ,sp.dom,"_",dates,"_",paste0(min_deplayer,"-",max_deplayer),".RData" ) )    
	Soutputs = file.path(DomProjRastLayer,paste0("Mercator_Salinity_",sp.dom,"_",dates,"_",paste0(min_deplayer,"-",max_deplayer),".RData" ) ) 

	if (  sum(file.exists(Toutputs)) != length(file.exists(Toutputs)) | sum(file.exists(Soutputs)) != length(file.exists(Soutputs)) )  { 

		for (date in dates){

			### 2.1 Open Global 3D Climat data path : T.Month.file  && S.Month.file
			Tclim <-  Tfiles[str_detect(Tfiles,date)]
			Sclim <-  Sfiles[str_detect(Sfiles,date)]

			### 2.2 Crop Global 3D climat files with local domain extent  (path : "RESULTS/ProjRasters/domain")
	
			Toutput = file.path(DomProjRastZall,paste0("Mercator_Temperature_" ,sp.dom,"_",date,".nc" ) )    
			Soutput = file.path(DomProjRastZall,paste0("Mercator_Salinity_",sp.dom,"_",date,".nc" ) )  

		    	if ( !file.exists(Toutput) | !file.exists(Soutput) ){

				### 2.2.1 Crop 3D Temperature raster with local domain extent
				nc_subset(filename=Tclim, compression=9,output= Toutput ,latitude=xlat, longitude=xlon,ignore.case = TRUE)
				### 2.2.2 Crop 3D Salinity raster with local domain extent
				nc_subset(filename=Sclim, compression=9,output= Soutput ,latitude=xlat, longitude=xlon,ignore.case = TRUE)

		    	 }# cropped files exists

	
			### Open Local 3D climat files .nc
			rTemp <- brick(Toutput,lvar=4) ; rm("Toutput")
			rSal <- brick(Soutput,lvar=4) ; rm("Soutput")


			### 2.3 Average across layers 3D local -> 2D local
			Toutput = file.path(DomProjRastLayer,paste0("Mercator_Temperature_" ,sp.dom,"_",date,"_",paste0(min_deplayer,"-",max_deplayer),".RData" ) )    
			Soutput = file.path(DomProjRastLayer,paste0("Mercator_Salinity_",sp.dom,"_",date,"_",paste0(min_deplayer,"-",max_deplayer),".RData" ) ) 


			if ( !file.exists(Toutput) | !file.exists(Soutput) ){

				cat("\t" ,date ," : ")

				### 2.3.1 subset layers between min and max of pre-fixed "layers3D"
			
				Mercator_layers <- which (between(Obs.depths, min_deplayer,  max_deplayer ) )
				
				Tlayers <- raster::subset(rTemp,Mercator_layers)
				Slayers <- raster::subset(rSal,Mercator_layers)
	
				### 2.3.2 mean raster 
			
				Tmean <- mean(Tlayers,na.rm=TRUE) ; coords = coordinates(Tmean)
				Tmean <- cbind(coords,values(Tmean))
				
				Smean <- mean(Slayers,na.rm=TRUE)
				Smean <- cbind(coords,values(Smean))
				
				### 2.3.3 save data.frame (coords, R values)
				save(Tmean,file=Toutput,version=2)
				save(Smean,file=Soutput,version=2)
	
				rm("Tmean","Smean","Tlayers","Slayers","Mercator_layers","coords")
				gc()

			}# if mean rasters exists

		rm("Toutput","Soutput")

		}#end of for dates



	}# eo IF all MONTHLY projection rasters done

	rm("Toutputs","Soutputs")

	
	### 3 ANNUAL RASTERS
	# IF all ANNUAL projection rasters done
	Toutputs = file.path(DomProjRastLayer,paste0("Mercator_Temperature_" ,sp.dom,"_",years,"_",paste0(min_deplayer,"-",max_deplayer),".RData" ) )    
	Soutputs = file.path(DomProjRastLayer,paste0("Mercator_Salinity_",sp.dom,"_",years,"_",paste0(min_deplayer,"-",max_deplayer),".RData" ) ) 


	if (  sum(file.exists(Toutputs)) != length(file.exists(Toutputs)) | sum(file.exists(Soutputs)) != length(file.exists(Soutputs)) )  { 

		for (year in years){

		## Annual 2D local rasters
		Toutput = file.path(DomProjRastLayer,paste0("Mercator_Temperature_" ,sp.dom,"_",year,"_",paste0(min_deplayer,"-",max_deplayer),".RData" ) )    
		Soutput = file.path(DomProjRastLayer,paste0("Mercator_Salinity_",sp.dom,"_",year,"_",paste0(min_deplayer,"-",max_deplayer),".RData" ) ) 


		if ( !file.exists(Toutput) | !file.exists(Soutput) ){

			cat("\t" ,year ," : ")

		  	year.files <- dir(DomProjRastLayer,pattern = paste0(year,"-") )
	  
			# Mean of Monthly Temperature files for year y
			Tyear.files <- grep(year.files,pattern="Temperature",value=T)
			rdate <- get( load(file.path(DomProjRastLayer,Tyear.files[1])) ) ; coords = rdate[,c(1,2)]  ; rm("rdate")
			TyearVal <- sapply(Tyear.files,function(TF){ rdate <- get( load(file.path(DomProjRastLayer,TF)) ) ; return( rdate[,3] ) } )
			Tyearmean <- apply(TyearVal, 1, mean,na.rm=TRUE) ; Tyearmean <- cbind(coords,Tyearmean)
	  
			# Mean of Monthly Salinity files for year y
			Syear.files <- grep(year.files,pattern="Salinity",value=T) 
			SyearVal <- sapply(Syear.files,function(TF){ rdate <- get( load(file.path(DomProjRastLayer,TF)) ) ; return( rdate[,3] ) } )
			Syearmean <- apply(SyearVal, 1, mean,na.rm=TRUE) ; Syearmean <- cbind(coords,Syearmean)
		
			# save R values
			save(Tyearmean,file=Toutput,version=2 )
			save(Syearmean,file=Soutput,version=2 )
	  
			rm("Tyearmean","Syearmean","TyearVal","SyearVal","year.files","Tyear.files","Syear.files","coords")
			gc()


		}#eo files exists

		rm("Toutput","Soutput")


		}#eo for year

		rm ("max_deplayer","min_deplayer")

	}#eo IF all ANNUAL projection rasters done

	rm("Toutputs","Soutputs")
	rm ("max_deplayer","min_deplayer")

}#eo for layers3D


### ---- END ----- ###
cat("-----> ","\t",c("RASTERS FOR PROJECTIONS done" ,"\n") )

