#########################################################################################################################
# I - Models Calibration	  											#					   
#	1 Create a date frame of Presence/ Generated Pseudo-Absences used to calibrate niche models     		#
#		1.2 remove duplicated occurrences (with the same climate conditions, to reduce the sampling bias)	#
#		1.3 Calculate Temp & Sal quantiles for species occurrences (Presences)					#
#		1.4 create a convex-hull surrounding Presences projected in the env background				#
#		1.5 Fixe Number of pseudo-absence									#
#		1.6 Pseudo-absence simulation : sampling of Abs outside of the convex hull surronding prensences	#
#		1.7 Merge Presence/ Generated Pseudo-Absences : PA file							#
#	2 Calibrate models												#
#		2.1 Formating data for biomod										#
#		2.2 Numbers of k in the  k-fold validation								#
#		2.3 Creating the SplitTable, to split the data set in k part						#
#		2.3 Biomod Calibration											#
#	3 Predict the presence of the species over the environemental background (Global Env)				#
#	4 Validate models												#
#		4.1 Create empty object to save validation results							#
#		4.2 Extract the eval values (TSS and ROC) & and add sp name column					#
#		4.3 k-fold cross validation using Boyce index								#
#															#
# II - Model averaging													#
#															#
# III - RESPONSE CURVE													#
#	1 Response Averaging												#
#	2 Mean and standard.deviation of probability and binary results							#
# IV - OPT/CRITICAL MIN/CRITICAL MAX TEMP    										#
#	1 meadian saliity when max proba de presence / Presence predicted						#
#	2 Response curve when fixed sailinity										#
#	3 Crtitical minimum /critical Maximum										#
#	4 Optimal Temp													#
#########################################################################################################################										

########################
#  Outputs directories #
########################


# Create direcotry to save nichfile : for future projections
nichfiles <- file.path(Maindir,"RESULTS/NicheModels","nichfiles") ; if(!dir.exists(nichfiles)){dir.create(nichfiles,recursive = T)}
setwd(nichfiles)

# Niche direcotry for Species to save modeling outputs
Niche_SpDir <- file.path(Maindir,"RESULTS/NicheModels",SP_name); if(!dir.exists(Niche_SpDir)){dir.create(Niche_SpDir,recursive = T)}


###  LOAD nichefiles if Niche model done for species for other study area #

spnichfiles  <- file.path(nichfiles, str_replace(SP.name," ","."))   # Niche modeling Biomod files 

BiomodClibationFile = file.path(Niche_SpDir,"Biomod_cal.RData") # Biomod Calibration file

if ( file.exists(BiomodClibationFile) && dir.exists(spnichfiles) ){

	# load(file.path(SpDir,"PA.RData"))
	load(file.path(Niche_SpDir,"ParamAverag.RData"))
	Biomod_cal = get(load(BiomodClibationFile))

	cat("-----> ","\t",c("NICHE models already calibrated" ,"\n") )

}else{

	#######################
	#  Models Calibration #
	#######################

	
	### 1 Create a date frame of Presence/ Generated Pseudo-Absences used to calibrate niche models

	# 1.1 Project all the combinaison of temperature and salinity encounter by the species occurrences (matchup) in the ENV Background (2D space)
	XY <- sp_occ_Clim[, c("Temp", "Sal")]  ; XY <- na.omit(XY)
	coordinates(XY) <- ~ Temp + Sal

	# 1.2 remove duplicated occurrences (with the same climate conditions, to reduce the sampling bias)
	XY.r <- rasterize(XY, Sp_back.r, fun = "count") ; XY.r <- XY.r > 0
	Nb.presence <- length(which(values(XY.r) > 0))
	Presence <- coordinates(XY.r)[which(values(XY.r) > 0), ]
	colnames(Presence) <- c("Temperature", "Salinity")

	# 1.3 Calculate Temp & Sal quantiles for species occurrences (Presences)
	XY.quantile <- sp_occ_Clim[, c("Temp", "Sal")]

	XY.quantile <- XY.quantile[which(XY.quantile[, 1] > quantile(XY.quantile[, 1], Q1) &
                                   XY.quantile[, 1] < quantile(XY.quantile[, 1], Q2) & 
                                   XY.quantile[, 2] > quantile(XY.quantile[, 2], Q1) & 
                                   XY.quantile[, 2] < quantile(XY.quantile[, 2], Q2)), ]

	coordinates(XY.quantile) <- ~ Temp + Sal

	# 1.4 create a convex-hull surrounding Presences projected in the env background
	hull <- rasterize(gConvexHull(XY.quantile), Sp_back.r)	
	values(hull) <- ifelse(is.na(values(hull)), 1, NA)
	hull <- Sp_back.r * hull
	XY.hull <- coordinates(hull)[which(is.na(values(hull)) == FALSE), ]

	# 1.5 Fixe Number of pseudo-absence
	Nb.absence  <- nba * nrow(Presence)
	Nb.absence <- ifelse(Nb.absence < 1000,1000,Nb.absence) # minimum of Abs = 1000

	# 1.6 Pseudo-absence simulation : sampling of Abs outside of the 
	sampledCells <- sample(1 : nrow(XY.hull), Nb.absence,replace = TRUE)
	pseudoAbs <- XY.hull[sampledCells, ]
	colnames(pseudoAbs) <- c("Temperature", "Salinity")

	# 1.7 Merge Presence/ Generated Pseudo-Absences : PA file
	PA <- as.data.frame( rbind(cbind(Presence, resp = rep(1, Nb.presence)) , cbind(pseudoAbs, resp = rep(0, Nb.absence)  )  ))
	PA$x <- 1:nrow(PA)
	PA$y <- 1:nrow(PA)


	### 2 Calibrate models

	# 2.1 Formating data for biomod
	CalibData <- BIOMOD_FormatingData(resp.var    = PA$resp, 
                                  expl.var  = PA[, c("Temperature", "Salinity")],
                                  resp.xy   = PA[, c("x","y")],
                                  resp.name = SP_name)

	# 2.2 Numbers of k in the  k-fold validation
	k <- 3 

	# 2.3 Creating the SplitTable, to split the data set in k part
	DataSplitTable <- BIOMOD_cv(CalibData, k = k, rep = 1, balance = "presences")

	# 2.3 Biomod Calibration
	Biomod_cal = BIOMOD_Modeling(data             = CalibData,
                             models           = models, 
                             models.eval.meth = c("TSS", "ROC"),
                             DataSplitTable   = DataSplitTable,
                             do.full.models   = FALSE,
                             SaveObj          = TRUE) 	; gc()


	### 3 Predict the presence of the species over the environemental background (Global Env)
	xyz<-cbind(idx=1:length(Sp_back_count),freq=Sp_back_count)
	xyz[is.na(xyz)]<-0
	M <- ifelse( sum(Sp_back_count,na.rm = T)<= 100000,sum(Sp_back_count,na.rm = T),100000) 
	idx<-sample(x=xyz[,"idx"], size=100000, replace=T,prob= xyz[,"freq"]/sum(xyz[,"freq"]))

	GlobalEnv<-coordinates(Sp_back_count.r)[idx,]
	colnames(GlobalEnv) <- c("Temperature", "Salinity")

	Biomod_current_global <- BIOMOD_Projection(modeling.output = Biomod_cal,
                                           new.env         = GlobalEnv, 
                                           selected.models = Biomod_cal@models.computed[1:(k * length(models))], 
                                           proj.name       = paste0(SP_name, "_GlobalEnv"),
                                           silent          = TRUE); gc()

	Biomod_current_global <- Biomod_current_global@proj@val
	
	### 4 Validate models

	# 4.1 Create empty object to save validation results
	Validation <- NULL

	# 4.2 Extract the eval values (TSS and ROC) & and add sp name column
	eval <- get_evaluations(Biomod_cal, as.data.frame = T)
	eval <- cbind(eval[, -4], Species = rep(SP_name, nrow(eval))) 
	Validation <- rbind(Validation, eval)

	# 4.3 k-fold cross validation using Boyce index
	names.Biomod_current_global <- sapply(Biomod_cal@models.computed,function(o){  paste ( unlist( str_split(o,"_") )[3] ,unlist( str_split(o,"_") )[4],sep="_")  })
	
	for (j in models){
		  for (kfold in 1 : k){
		    idModel <- grep(names.Biomod_current_global, pattern = paste("RUN", kfold,"_", j, sep=""))
		 if(length(idModel)>0){ 
		    val.data <- PA[DataSplitTable[, kfold] == FALSE, c("Temperature", "Salinity", "resp")]
		    val.data <- val.data[val.data$resp == 1, -3]

		    val.data <- BIOMOD_Projection(modeling.output = Biomod_cal, # Biomod projection on the climatic map
                                  new.env         = val.data, 
                                  selected.models = Biomod_cal@models.computed[idModel], 
                                  proj.name       = paste(Biomod_cal@models.computed[idModel], "_Current_Global", sep = ""),
                                  silent          = TRUE) ;gc()
		    
		    CBI <- ecospat.boyce (fit      = na.omit(Biomod_current_global[,j,paste0("RUN",kfold),]/1000), 
                          obs      = val.data@proj@val[, 1, 1, 1]/1000, 
                          nclass   = 0, 
                          window.w ="default", 
                          res      = 100, 
                          PEplot   = FALSE)
		    Validation <- rbind( Validation, 
                         data.frame(Model.name   = paste("RUN", kfold, "_", j, sep = ""),  
                                    Eval.metric  = "CBI",  
                                    Testing.data = CBI$Spearman.cor, # MAEL GERNEZ rather suggests using :  CBI$cor
                                    Cutoff       = NA, 
                                    Sensitivity  = NA,
                                    Specificity  = NA,  
                                    Species      =  SP_name)) ; gc()
			}
	  }#eo kfold
	}#eo j in models

	

	#################################
	#     Model averaging           #
	#################################

	# 1 Extract RUN 4 (final run, with all data)
	for (i in 1:4){Validation[grep(Validation$Model.name,pattern = as.character(i)), "Rep"] <- i}
	for (i in models){Validation[grep(Validation$Model.name,pattern = as.character(i)), "Algorithm"] <- i} 
	Validation <- Validation[which(Validation$Rep != 4), c(2:9)]
	Validation <- Validation[, c("Eval.metric", "Testing.data", "Cutoff", "Rep", "Species", "Algorithm")]
	Validation$Eval.metric <- as.factor(Validation$Eval.metric)
	Validation$Species <- as.factor(Validation$Species)
	Validation$Algorithm <- as.factor(Validation$Algorithm)
	bem.tmp <- Validation[which(Validation$Species ==SP_name), ] 
	tmp.cbi <-  bem.tmp[which(bem.tmp$Eval.metric == "CBI"), ]
	Mean.tmpcbi <- aggregate(tmp.cbi[, 2], list(tmp.cbi$Algorithm), mean)
	

	# if at least one good model else stop
	if (  length( which(Mean.tmpcbi$x > 0.5) ) > 0 ) { 

		# 2 Asking for the mean of the CBI of all the 3 cross validation to be higher than 0.5
		sel.models <- names(which(table(Mean.tmpcbi[which(Mean.tmpcbi$x > 0.5), "Group.1"]) == 1)) 
		sel.cbi <-  tmp.cbi[which(tmp.cbi$Algorithm %in% sel.models), c("Algorithm", "Testing.data")]

		cbi <- aggregate(sel.cbi$Testing.data, by = list(sel.cbi$Algorithm), FUN = mean)    
		cbi <- cbi$x; names(cbi) <- sel.models

		# 3 Cutoff setting from TSS 
		sel.cutoff <- na.omit(Validation[which(Validation$Species ==  SP_name & 
                                         Validation$Algorithm %in% sel.models &  
                                         Validation$Eval.metric %in% c("TSS", "CBI")), ])
		tss <- aggregate(sel.cutoff$Testing.data, by = list(sel.cutoff$Algorithm), FUN = mean)
		tss <- tss$x; names(tss) <- sel.models
		cutoff <- aggregate(sel.cutoff$Cutoff, by = list(sel.cutoff$Algorithm), FUN = mean)
		cutoff <- cutoff$x; names(cutoff) <- sel.models
		cutoff <- weighted.mean(cutoff, w = tss)/1000 # Cutoff for thr binary map from the TSS

		# SAVE MODELING OUTPUTS  in  SpDir
		save(PA,Nb.presence,file= file.path(Niche_SpDir,"PA.RData"),version=2)
		save(Biomod_current_global,file= file.path(Niche_SpDir,"Biomod_current_global.RData"),version=2)
		save(Validation,file= file.path(Niche_SpDir,"Validation.RData"),version=2)
		save(GlobalEnv,file= file.path(Niche_SpDir,"GlobalEnv.RData"),version=2)
		save(Biomod_cal,file= file.path(Niche_SpDir,"Biomod_cal.RData"),version=2)
		save(Mean.tmpcbi,tss,cbi,sel.cutoff,sel.models,sel.cbi,cutoff,file=file.path(Niche_SpDir,"ParamAverag.RData"),version=2)


		###########################
		#   RESPONSE CURVE        # 
		###########################
 
		new.env <- coordinates(Sp_back_count.r)[which(values(Sp_back_count.r) > 0),]
		colnames(new.env) <- c("Temperature", "Salinity") ; gc()

		CurveFull <- BIOMOD_Projection(Biomod_cal, 
                               new.env = new.env, 
                               selected.models = Biomod_cal@models.computed[grep(Biomod_cal@models.computed,pattern = "4")], 
                               proj.name = "_NewEnv",
                               silent=TRUE);gc()

		CurveFull <- as.data.frame(cbind(new.env,CurveFull@proj@val[, , 1, 1]/1000))

		# save Response 2D (all models)
		save(CurveFull,file=file.path(Niche_SpDir,"CurveFull.RData"),version=2)
		coordinates(CurveFull) <-~ Temperature+Salinity
		CurveFull <- as(CurveFull,"SpatialPixelsDataFrame") # Use to plot the Full responce curve
		gc()  


		# 1 Response Averaging
		response <- CurveFull
		response@data <- as.data.frame(CurveFull@data[, sel.models])
    
		# 2 Mean and standard.deviation of probability and binary results
		averaging.response <- response[, 1]
		averaging.response@data[, 1] <- apply(response@data,1, FUN = function(x) weighted.mean(x, w = tss)) 
		averaging.response@data[, 2] <- apply(response@data,1, FUN = sd) # standard deviation map
		if(all(is.na(averaging.response@data[, 2]))){averaging.response@data[, 2] <- 0}
		averaging.response@data[, 3] <- ifelse(averaging.response@data[, 1] > cutoff, 1, 0) # binary map
		names(averaging.response) <- c("Mean", "Standard.deviation", "Binary")
  

		# Save results 
		averaging.responseDF <- cbind( averaging.response@coords,averaging.response@data)
		save(averaging.responseDF,file=file.path(Niche_SpDir, paste0("Response_Averaging.RData") ),version=2 )
		gc()
    
		models.responseDF <- cbind( response@coords,response@data)
		save(models.responseDF,file=file.path(Niche_SpDir, paste0("Response_BestMods.RData") ),version=2 )

   

		#########################################
		#   OPT/CRITICAL MIN/CRITICAL MAX TEMP  #
		#########################################
  
		DF <- data.frame(cbind(round(new.env,2),averaging.responseDF[,c("Mean","Binary")]))
    
		# 1 meadian saliity when max proba de presence / Presence predicted
		test <- subset (DF,between(Mean,max(DF$Mean)-0.1,max(DF$Mean)+0.1))
		meansal = median( test[test$Binary == 1,"Salinity"] ) ; rm("test")

		# 2 Response curve when fixed sailinity
		sub.df <- subset (DF,between(Salinity,meansal-0.2,meansal+0.2))
		sub.df<- sub.df[order(sub.df$Temperature),]

		# 3 Crtitical minimum /critical Maximum
		Tcritmin <- min ( unique(DF[DF$Binary==1,"Temperature"]) )
		Tcritmax <- max ( unique(DF[DF$Binary==1,"Temperature"]) )

		# 4 Optimal Temp
		test <- subset (sub.df,between(Mean,max(sub.df$Mean)-0.1,max(DF$Mean)+0.1))
		test <- subset (test,Binary==1)
		Topt <-  round(mean(test$Temperature) ,2)   ; rm("test")

		TemperatureLimits <- list("Toptimal"=Topt , "Tmincritical" = Tcritmin , "Tmaxcritical" = Tcritmax) ; rm("Topt","Tcritmin","Tcritmax")

		# Save results 
		save(TemperatureLimits,file= file.path(Niche_SpDir,paste0(  SP_name ,"_TemperatureLimits.RData") ) ,version=2)
		save(sub.df,file= file.path(Niche_SpDir,paste0(  SP_name ,"_sub.df.RData") ) ,version=2)
		save(DF,file= file.path(Niche_SpDir,paste0(  SP_name ,"_DF.RData") ) ,version=2)

		rm("models.responseDF" , "response","CurveFull","sub.df","DF","TemperatureLimits")
		rm("averaging.responseDF" , "averaging.response") ; gc()

		rm("bem.tmp", "Biomod_current_global","CalibData","cbi","CBI","DataSplitTable","names.Biomod_current_global","Nb.absence","Nb.presence","hull",
		"idModel","idx","M","Mean.tmpcbi","sel.cbi","sel.cutoff","tmp.cbi","val.data","Validation","XY","XY.hull","XY.quantile","XY.r","xyz","sampledCells",
		"Presence","pseudoAbs","PA","GlobalEnv","eval")

		cat("-----> ","\t",c("NICHE models calibrated" ,"\n") )


	# if at least one good model 
	}else{
		cat("-----> ","\t",c("REDO NICHE calibration" ,"\n") ) } # models do not converge 

} #eo else calibration done / not done

