############################################################################################
# I - RAW DATA : Download species occurrenes at global scale from global databases	   #					    
#	1 OBIS occurrences								   #
#	2 GBIF and other database							   #
#	3 SAVE RAW DATA SpeciesOCC						   	   #
# 	    columns c("longitude","latitude","name","prov","year","month")		   #
#								   			   #
# II - PROCESSING DATA									   #
#	1 Remove record with zero in Lat and Long					   #
#	2 Remove occurrences on the land using a raster of climate data			   #
#	3 Genius.species' records (sometimes other species records in the table)	   #
#	4 Clean the month and year of the observations (remove NAs, etc)		   #
#	5 Remove duplicated records (cells with same long/lat/year/month)		   #
#	6 SAVE PROCESSED DATA								   #
############################################################################################										

########################
#  Outputs directories #
########################

### We create Outputs directories for Species raw occurrences
Sp_occ_Dir_Raw <- file.path(Maindir,"RESULTS/SpeciesOccurrences/raw/") ; if(!dir.exists(Sp_occ_Dir_Raw)){dir.create(Sp_occ_Dir_Raw,recursive = TRUE)}
# & processed occurrences 
Sp_occ_Dir <- file.path(Maindir,"RESULTS/SpeciesOccurrences")


### Load the processed occurrences file if it already exists

if ( file.exists( file.path(Sp_occ_Dir,paste0(SP_name,".RData") ) ) ) {
    
	SpeciesOCC <- get( load(  file.path(Sp_occ_Dir,  dir(Sp_occ_Dir,pattern = paste0(SP_name,".RData") ) ) ) ) 
	cat("-----> ","\t",c("OCCURRENCES already processed" ,"\n") )

}else{ 

#####################################################################
#			 	RAW DATA			    #
# "obis","gbif", "bison", "ecoengine", "vertnet"		    #
# keep columns c("longitude","latitude","name","prov","year","month"#
#####################################################################
  

### Load raw occurrences file if it already exists

if( file.exists(  paste0(Sp_occ_Dir_Raw,"occ_",SP_name,".RData") ) ) {

	SpeciesOCC <- get(load(paste0(Sp_occ_Dir_Raw,"occ_",SP_name,".RData")  ))
	cat("-----> ","\t",c("OCCURRENCES already downloded" ,"\n") )

}else{

	SpeciesOCC <- NULL

	###### 1 OBIS occurrences 
	
	# 1.1 Download species occurrenes at global scale from gloabl databases  : obis/gbif etc 

	occOBIS <- occurrence(SP.name)
	
	if ( !is.null( dim(occOBIS)[1] ) && dim(occOBIS)[1]>0  ) {
	
	# 1.2  keep columns c("longitude","latitude","name","prov","year","month")	
		year = try( pull(occOBIS,"date_year") ,silent=T)
		month = try( pull(occOBIS,"month") ,silent=T)
		if( inherits(year,"try-error") ) year = substr(pull(occOBIS,"eventDate") ,1,4)
		if(inherits(month,"try-error") ) month = substr(pull(occOBIS,"eventDate") ,6,7)
		occOBIS <- data.frame(occOBIS[,c("decimalLongitude","decimalLatitude", "species")],"prov"=rep("obis",nrow(occOBIS)),year,month)
		names(occOBIS) <- c("longitude","latitude","name","prov","year","month")
		SpeciesOCC <- rbind(SpeciesOCC,occOBIS)
	} #if there are obis records 
  

	###### 2 GBIF and other database

	# 2.1 Download occurrences from GBIF and other database 

	occOther <- occ(SP.name,limit =50000,from=c("gbif", "bison", "ecoengine", "vertnet"))
	occOther <- occ2df(occOther,"all") 
	occOther <- occOther$data

	# 2.2 keep columns c("longitude","latitude","name","prov","year","month")

	if ( !is.null( dim(occOther)[1]) && dim(occOther)[1]>0 ){
		if ( sum( str_detect(colnames(occOther),"date") ) > 0 ){
		occOther <- data.frame( occOther[,c("longitude","latitude","name","prov")], "year"= as.numeric(substr(as.character(occOther[,"date"]),1,4)),
			 		"month"= as.numeric(substr(as.character(occOther[,"date"]),6,7)) )
		SpeciesOCC <- rbind(SpeciesOCC,occOther)
		}#if there is no date for each occurrence
	}#if there are gbif records 
  

	###### 3 SAVE RAW DATA SpeciesOCC

        if ( !is.null(SpeciesOCC) && nrow(SpeciesOCC) > 0 ){

		Name <- paste0 ("occ_",SP_name,".RData")
		save(SpeciesOCC,file = file.path(Sp_occ_Dir_Raw,Name),version=2) 
		rm(list=c("occOBIS","occOther","Name")) 
		cat("-----> ","\t",c("OCCURRENCES downloded" ,"\n") )

	}else{
		cat("-----> ","\t",c("NO OCCURRENCES in global databases / CHECK SPECIES NAME" ,"\n") )
	}

} # raw data
  


###################################### 
#      	PROCESSING DATA 	     #
######################################


	###### 1 Remove record with zero in Lat and Long

	SpeciesOCC$longitude <- as.numeric(SpeciesOCC$longitude)
	SpeciesOCC$latitude <- as.numeric(SpeciesOCC$latitude)
	no.coord <- which( is.na(SpeciesOCC$longitude) | is.na(SpeciesOCC$latitude) )
	if(length(no.coord)>0) { SpeciesOCC<-SpeciesOCC[-no.coord,]  }

	
	###### 2 Remove occurrences on the land using a raster of climate data

	# 2.1 Create spatial object

	
	coordinates(SpeciesOCC) <- ~ longitude+latitude
 
	# 2.2

	#proj4string(SpeciesOCC) <- proj4string(RasterClim)
	sea <- raster::extract(RasterClim,SpeciesOCC)
	if ( length(which(is.na(sea)))>0 ) {  SpeciesOCC <- SpeciesOCC[-which(is.na(sea)),] }
  

	###### 3 sometimes we have other species records in the table -> keep only our species' records

	SpeciesOCC <- cbind(coordinates(SpeciesOCC),SpeciesOCC@data)   # re-convert to data.frame
	SpeciesOCC <- SpeciesOCC[grep(SpeciesOCC$name,pattern=SP.name),]


	###### 4 Correct the month and year of the observations

	# 4.2 Clean years
          
	SpeciesOCC$year <- as.numeric(SpeciesOCC$year)

	SpeciesOCC <- SpeciesOCC[ !is.na(SpeciesOCC$year),]		# remove NAs

	SpeciesOCC <- SpeciesOCC[ between( SpeciesOCC$year ,  1993,2018), ]     # keep only records between 1993 & 2018 (Mercator periods)

    
	# 4.3 Clean months    
          
	SpeciesOCC$month <- as.numeric(SpeciesOCC$month)

	SpeciesOCC <- SpeciesOCC[!is.na(SpeciesOCC$month),]	# remove NAs

	SpeciesOCC <- SpeciesOCC[which(SpeciesOCC$month %in% 1:12),]		# keep only months from 01 to 12


	###### 5 Remove duplicated records (cells with same long/lat/year/month)

	SpeciesOCC2 <- SpeciesOCC[,c("longitude","latitude","year","month")]
	SpeciesOCC <- SpeciesOCC[!duplicated(SpeciesOCC2),]

	rm(list=c("no.coord","sea","SpeciesOCC2"))

	###### 6 SAVE PROCESSED DATA

        if ( !is.null(SpeciesOCC) && nrow(SpeciesOCC) > 0 ){

		Name <- paste0 (SP_name,".RData")
		save(SpeciesOCC,file=file.path(Sp_occ_Dir,Name),version=2)
		rm("Name")
		cat("-----> ","\t",c("OCCURRENCES processed" ,"\n") )

	}else{
		cat("-----> ","\t",c("NO OCCURRENCES after data processing" ,"\n") )
	}

}#eo else
