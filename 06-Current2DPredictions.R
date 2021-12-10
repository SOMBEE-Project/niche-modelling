####################################################################################################################
# LOOP to predict species distributions / Month & Year / 3D.layer			  			   # 
#														   #
#	1 Create sub directories to save predictions / 3D Layer							   #
#	2 MONTHLY PREDICTIONS											   #
#		2.1 Open Global 3D Climat data path : T.Month.file  && S.Month.file				   #
#		2.2 Create a raster stack with env. predictors (Temp & Salinity)				   #
#		2.3 Biomod Prediction										   #
#		2.4 Find the selected model(s)    								   #
#		2.5 Maps of mean and standard.deviation of probability and binary result (current)		   #
#														   #
#	3 ANNUAL PREDICTIONS											   #
#		3.1 Open Global 3D Climat data path : T.Month.file  && S.Month.file				   #
#		3.2 Create a raster stack with env. predictors (Temp & Salinity)				   #
#		3.3 Biomod Prediction										   #
#		3.4 Find the selected model(s)    								   #
#		3.5 Maps of mean and standard.deviation of probability and binary result (current)		   #
####################################################################################################################

########################
#  Outputs directories #
########################


#--------# Sub-directory to save predicted distribution for each study area 
	Niche_SpDomDir <- file.path(Niche_SpDir,sp.dom) ; if(!dir.exists(Niche_SpDomDir)){dir.create(Niche_SpDomDir,recursive = T)}

#----------------# Sub-directory to save predicted distribution for historical period using Mercator Climate data / Depth layer (Layers3D)
		Niche_SpDomDirCurrent  <- file.path(Niche_SpDomDir,"CurrentPredictions")  ; if(!dir.exists(Niche_SpDomDirCurrent)){dir.create(Niche_SpDomDirCurrent,recursive=T)} 


###  LOAD nichefiles if Niche model done for species for other study area #

spnichfiles  <- file.path(nichfiles, str_replace(SP.name," ","."))   # Niche modeling Biomod files 

BiomodClibationFile = file.path(Niche_SpDir,"Biomod_cal.RData") # Biomod Calibration file


#######################################################################
#   LOOP to predict species distributions / Month & Year / 3D.layer   #
# from Yf = 1993 to Yt = 2018 (see LAUNCH_xxx.R)        	      #
# 2 * 12 * 25 = 624 file / 3D.layer			      	      #
#######################################################################

# DO PREDICTIONS ONLY FOR layers3D between species minimum and maximum depth
SPlayers3D = as.character(unique(cut(min_depth:max_depth,breaks= c(0,25,50,100,150,200,300,500,1000,max(Obs.depths)),include.lowest = T)))
 
###### monthly from 01-Yf to 12-Yt
dates = paste0 (  rep(seq(Yf,Yt,1), each = 12) ,   "-" ,rep(str_pad(seq(1,12,1),2,side = "left" ,"0"), length(rep(seq(Yf,Yt,1)))  ))

###### annual from Yf to Yt
years <- seq(Yf,Yt,1)


   ###### for each of layers3D  

	for ( layer in SPlayers3D ){

		cat("\n",layer ," ; ")
		regexp <- "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?"
		min_deplayer = as.numeric(str_extract_all(layer, regexp,simplify = T)[1,1])
		max_deplayer = as.numeric(str_extract_all(layer, regexp,simplify = T)[1,2])

		### 1 Create sub directories to save Rasters
		# sub direcotry to save historical prediction raster for depth layers in layers3D
		Niche_SpDomLayerDirCurrent <-  file.path(Niche_SpDomDirCurrent,paste0(min_deplayer,"-",max_deplayer) ) 
		if (!dir.exists (Niche_SpDomLayerDirCurrent)) 	dir.create(Niche_SpDomLayerDirCurrent,recursive = TRUE)
		# sub direcotry where historical climate rasters for depth layers in layers3D are saved (04.xxxxxxx.R)
		DomProjRastLayer <-  file.path(DomProjRast,paste0(min_deplayer,"-",max_deplayer) ) 

		# 2 MONTHLY PREDICTIONS
		### IF Files to save MONTHLY predictions for date/layer
		averaging_current_dates =  file.path(Niche_SpDomLayerDirCurrent,paste0("Averaging_",dates,".RData"))
	  
		if (  sum(file.exists(averaging_current_dates)) != length(file.exists(averaging_current_dates)) )  { 

			for (date in dates){

				cat (date , "\t")

				### 2.1 Open Global 3D Climat data path : T.Month.file  && S.Month.file
				Tfile <-  grep(dir(DomProjRastLayer),pattern=date,value=T) ; Tfile <-  grep(Tfile,pattern="Temperature",value=T)
				rTemp <- get(load( file.path(DomProjRastLayer,Tfile) ))
				TemperatureL <-  rasterFromXYZ(rTemp)

				Sfile <-  grep(dir(DomProjRastLayer),pattern=date,value=T); Sfile <-  grep(Sfile,pattern="Salinity",value=T)
				rSal <- get(load(file.path(DomProjRastLayer,Sfile))) 
				SalinityL <-  rasterFromXYZ(rSal) 


				averaging_current_date =  file.path(Niche_SpDomLayerDirCurrent,paste0("Averaging_",date,".RData"))
		  
				if( !file.exists(averaging_current_date) ) {

					# 2.2 Create a raster stack with env. predictors (Temp & Salinity)
					CurrentEnv <- stack(TemperatureL,SalinityL)
					names(CurrentEnv) <- c("Temperature", "Salinity") ; gc()  

					# 2.3 Biomod Prediction
					Biomod_current <- BIOMOD_Projection(Biomod_cal, 
				                      new.env         = CurrentEnv, 
				                      proj.name       = paste(SP_name,sp.dom,date,sep="_"),
				                      selected.models = Biomod_cal@models.computed[grep(Biomod_cal@models.computed, pattern = "4")],
				                      silent          = TRUE) ; gc()
					current_distribution <- as(Biomod_current@proj@val, "SpatialPixelsDataFrame")
		  
					# 2.4 Find the selected model(s)    
					sel.models.grep = NULL
					for (z in 1 : NROW(sel.models)){
						grep1 <- grep(sel.models[z], names(current_distribution@data), value = TRUE)
						sel.models.grep <- rbind(sel.models.grep, data.frame(grep1))  }
		  
		  
					current_distribution@data <- as.data.frame(current_distribution@data[, as.character(unlist(sel.models.grep))])
					names(current_distribution) <- sel.models

					# 2.5 Maps of mean and standard.deviation of probability and binary result (current)
					averaging.current <- current_distribution[, 1]
					averaging.current@data[, 1] <- apply(current_distribution@data, 1, FUN = function(x) weighted.mean(x, w = tss))/1000 
					averaging.current@data[, 2] <- apply(current_distribution@data, 1, FUN = sd)/1000 # standard deviation map
					if(all(is.na(averaging.current@data[, 2]))){averaging.current@data[, 2] <- 0}
					averaging.current@data[, 3] <- ifelse(averaging.current@data[, 1] > cutoff, 1, 0) # binary map
					names(averaging.current) <- c("Mean", "Standard.deviation", "Binary")
		 
		  
					# Save results 
					current_distribution_DF <-   cbind( current_distribution@coords, current_distribution@data) 
					save(current_distribution_DF,file = file.path(Niche_SpDomLayerDirCurrent,paste0("BestModels_",date,".RData"   ) ),version=2) 
					averaging_current_DF <- cbind( averaging.current@coords,averaging.current@data)
					save(averaging_current_DF,file = averaging_current_date,version=2) 

					rm("averaging_current_DF","averaging.current","current_distribution_DF","current_distribution","Biomod_current","CurrentEnv")

				}#eo file exists  

				rm("Tfile","rTemp","TemperatureL","Sfile","rSal","SalinityL")
		 
			} #eo date in dates 



	}# eo IF all MONTHLY projection rasters done



		### 3 ANNUAL PREDICTIONS

		for (year in years){

			cat(year,"\t")

			### 3.1 Open Temperature/Salinity climate rasters for year / study area / depth layer
			Tfile <-  grep(dir(DomProjRastLayer),pattern=paste0(year,"_"),value=T) ; Tfile <-  grep(Tfile,pattern="Temperature",value=T)
			rTemp <- get(load( file.path(DomProjRastLayer,Tfile) ))
			TemperatureL <-  rasterFromXYZ(rTemp) # TemperatureL <- resample(TemperatureL, osmose.grid)

			Sfile <-  grep(dir(DomProjRastLayer),pattern=paste0(year,"_"),value=T) ; Sfile <-  grep(Sfile,pattern="Salinity",value=T)
			rSal <- get(load(file.path(DomProjRastLayer,Sfile))) 
			SalinityL <-  rasterFromXYZ(rSal) 

		
			### Name the file to save predictions for year/layer (save model "Averaging"  & Best models)
			averaging_current_year =  file.path(Niche_SpDomLayerDirCurrent,paste0("Averaging_",year,".RData"))
	  

			if (!file.exists(averaging_current_year) ) {

				# 3.2 Create a raster stack with env. predictors (Temp & Salinity)
				CurrentEnv <- stack(TemperatureL,SalinityL)
				names(CurrentEnv) <- c("Temperature", "Salinity") ; gc()  

				# 3.3 Biomod Prediction
				Biomod_current <- BIOMOD_Projection(Biomod_cal, 
		                              new.env         = CurrentEnv, 
		                              proj.name       = paste(SP_name,sp.dom,year,sep="_"),
		                              selected.models = Biomod_cal@models.computed[grep(Biomod_cal@models.computed, pattern = "4")],
		                              silent          = TRUE) ; gc()
				current_distribution <- as(Biomod_current@proj@val, "SpatialPixelsDataFrame")
	  
	  
				# 3.4 Find the selected model(s)    
				sel.models.grep = NULL
				for (z in 1 : NROW(sel.models)){
				    grep1 <- grep(sel.models[z], names(current_distribution@data), value = TRUE)
				    sel.models.grep <- rbind(sel.models.grep, data.frame(grep1)) }
	  
				current_distribution@data <- as.data.frame(current_distribution@data[, as.character(unlist(sel.models.grep))])
				names(current_distribution) <- sel.models

				# 3.5 Maps of mean and standard.deviation of probability and binary result (current)
				averaging.current <- current_distribution[, 1]
				averaging.current@data[, 1] <- apply(current_distribution@data, 1, FUN = function(x) weighted.mean(x, w = tss))/1000 
				averaging.current@data[, 2] <- apply(current_distribution@data, 1, FUN = sd)/1000 # standard deviation map
				if(all(is.na(averaging.current@data[, 2]))){averaging.current@data[, 2] <- 0}
				averaging.current@data[, 3] <- ifelse(averaging.current@data[, 1] > cutoff, 1, 0) # binary map
				names(averaging.current) <- c("Mean", "Standard.deviation", "Binary")

				# Save results 
				current_distribution_DF <-   cbind( current_distribution@coords, current_distribution@data) 
				save(current_distribution_DF,file = file.path(Niche_SpDomLayerDirCurrent,paste0("BestModels_",year,".RData"   ) ),version=2) 
				averaging_current_DF <- cbind( averaging.current@coords,averaging.current@data)
				save(averaging_current_DF,file =  averaging_current_year,version=2) 
				gc()
	  
				rm("averaging_current_DF","averaging.current","current_distribution_DF","current_distribution","Biomod_current","CurrentEnv")


			}#eo file exists

		rm("Tfile","rTemp","TemperatureL","Sfile","rSal","SalinityL")
	  
		} #eo years

	}#eo SPlayers3D

		cat("-----> ","\t",c("Monthly 2D-Distribution predicted" ,"\n") )

