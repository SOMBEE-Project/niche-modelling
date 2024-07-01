############################### Launch SDM #####################################
#
# Script from Sabrine Drira, from Sombee project with Maël Gernez and Marine Beneat modifications. 
#
# Used to build a thermal niche model in order to estimate temperature T min, max and opt. These temperatures will be used to estimate Bioen-Osmose parameters.
# 
# After running all scripts for all species on the study period, only one result is kept in RESULTs/NicheModels/"Sp_name"/"Sp_nam"_TemperatureLimits.RData
#
# Scripts modifications are highlighted with comments. Where you should add your own modification is highlighted with : ##* NEEDS TO BE ADAPTED TO YOUR DATASET.
#
##################
#									 Sabrine comments :
##															    #	
# 1 - COLLECT SPECIES OCCURRENCES DATA   : 01-Species_occurrences.R							    #
#	Download species occurrenes at global scale from global databases						    #
#															    #	
# 2 - CREATE BACKGROUND for each species : 02-Backgrounds.R								    #
# 	The environnemental background is a grid encompassing all the combinations of temperature and salinity occurring at #
#	the global scale. To generate this 2D space, we calculated mean climate raster using only layers between min max    #
#	depth (or bottom, surface) for each monthly climat 3D file (Temp/salinity)					    #
#	We used the env background for filtering of occurrences data (sampling bias) and for pseudo-absences generation	    #
# 	Mean Rasters used to generate environnemental background are used for the climat/occ match-up			    #
#															    #	
# 3 - MATCH-UP between climate (Sal & Temp) and species occurrences : 03-Matchup_OccClimat.R				    #
#	A spatiotemporal match-up between climatic climatologies and species occurrences is realised, taking into account   #
#	the geographic coordinates of occurrences as well as their corresponding decade and the vertical layer that	    #
#	corresponds to the vertical habitat of the considered species (mean between min max, bottom, surface)		    #
#															    #
# 4 - GENERATE LOCAL RASTERS for historical predictions									    #
#	We aggregated the ~55 mercator's depth layers into 9 new depth layers as following (layers3D)			    #
#	"[0,25]","(25,50]","(50,100]","(100,150]","(150,200]","(200,300]","(300,500]","(500,1000]","(1000,bottom]"	    #
#															    #
# 5 - CALIBRATE NICHE MODELS												    #
#															    #
# 6 - PREDICT historical species distributions										    #
#	We predict 3D species distribution  : monthly predictions for each of 9 deth layers				    #
#	"[0,25]","(25,50]","(50,100]","(100,150]","(150,200]","(200,300]","(300,500]","(500,1000]","(1000,bottom]"	    #
#	We agreagated 3D predictions into one final map									    #
# 															    #
# After Downscaling....													    #	
# 7 - PROJECT Future species distributions										    #
#															    #	
##############################################################################################################################


rm(list=ls())

# Main Directory for Species Distribution Modeling
Maindir <- file.path(getwd())  
setwd(Maindir)

####################################
#   Load setup file : librairies   #	
####################################

source(file.path(Maindir,"SCRIPTS/00_setup.R"))


##################################################################
#            Climate data directory & files                      #
##################################################################

climDir <- "/home/datawork-marbec-pmod/Sabrine/Mercator/new-values" # Datarmor folder

# If Glorys files with 1 variable at a time 
#Tfiles <- grep( dir(climDir) , pattern= "thetao" , value=T)
#Tfiles <- file.path ( climDir , Tfiles )
#Sfiles <- grep( dir(climDir) , pattern= "so" , value=T)
#Sfiles <- file.path ( climDir , Sfiles )

# If Glorys files have various climatic variables per file 
Climfiles <- list.files(climDir, "*.nc")
Climfiles <- file.path( climDir , Climfiles )
Tfiles <- Climfiles
Sfiles <- Climfiles
RasterClim <- raster::brick(Climfiles[1], varname="thetao")
Obs.depths <- as.numeric(as.vector(unlist(RasterClim@z)))
RasterClim <- RasterClim[[1]]

# PREDICTION/PROJECTIONS 9 depth layers : 
layers3D <- levels(cut(0:max(Obs.depths),breaks= c(0,25,50,100,150,200,300,500,1000,max(Obs.depths)),include.lowest = T))  

##################################################################
#                	STUDY AREA		                 #
##################################################################
 
# LIST OF SPATIAL EXTENT FOR EACH ECOSYSTEM (FOR SP DISTRIBUTIONS PREDICTIONS)
DOMextents <- list()
# DOMextents$BoB <- extent(-8,0,43,48) # from ICES shp...
# DOMextents$NS <- extent(-5,10.5,51,61)
# DOMextents$BS <- extent(27.56,42.04,40.775,47.975)
# DOMextents$GL <- extent(2,9,40,45)
# DOMextents$WCC <- extent(-135.4232,-122.1903,47.577,56.0678)
# DOMextents$YS <- extent(119.2,126.5,31.4,39.5)
# DOMextents$HUM <- extent(-93,-70,-20,-6)
DOMextents$MED <- extent(-6,36,29,46)


##############################################
#          Species list                      #
##############################################

# LIST OF SPECIES FOR ALL ECOSYSTEMS

ALLSp <- read.csv( file.path(Maindir, "RESULTS/ALLSpecies_List.csv" ),sep=";")

# ALLSp$Layers <- as.character(ALLSp$Layers)
# ALLSp$Species <- as.character(ALLSp$Species)


################
#  SIMULATION  #
################

# create a file .txt to follow the modeling process at different stages, for each species 

WorkProgress <- data.frame( cbind("Species"= character(nrow(ALLSp)),"Domain"=character(nrow(ALLSp)),"sp_layers"=character(nrow(ALLSp)),"Occurences"=numeric(nrow(ALLSp)),
		"Match-up"=numeric(nrow(ALLSp)),"Nichemodeling"=character(nrow(ALLSp)),"CurrentPredictions"=character(nrow(ALLSp))),stringsAsFactors=FALSE)

WorkProgress$Species <- ALLSp$Species
WorkProgress$Domain <- ALLSp$Domain

# YEARS From xxxx & xxxx for monthly predictions
Yf = 1993
Yt = 2020

# MODELS USED FOR NM
models = c( "GLM","GBM","GAM","CTA","ANN","FDA","RF")#,"MARS")

# QUANTILES TO CREATE PSEUDO-ABS ( See SCRIPT : xxxxxxxxxxxxxxxx.R ) 

Q1= 0.05 ; Q2 = 0.95

# Numbr of pseudo-abs = nba * Presence
nba = 1

# RESOLUTION OF BACKGROUND (env space used to create pseudo-abs ; see SCRIPT : xxxxxxxxxxxxxxxxx.R ; created with resolution of 0.1)
fact=3


##################################################################
#       	 LOOP for species niche modeling                 #
##################################################################


for (sp in  1:nrow(ALLSp) ){
  
  ## MAEL GERNEZ, change due to change in setwd in 05 so:
  setwd(Maindir)
    
  SP <-   ((ALLSp$Species))[sp]
  SP_name <- str_replace_all(SP," ","_")
  SP.name <- str_replace_all(SP,"_"," ")
  
  cat( paste0("\n +  ",SP_name,"\t : \n"))


  # the species’ vertical habitat
   min_depth <- ALLSp[sp,"Minimum_depth"] 
   max_depth <-  ALLSp[sp,"Maximum_depth"] 
   sp_layers <- which (between(Obs.depths, min_depth ,  max_depth ) )
   sp_layers <- paste(  sp_layers[1],  sp_layers[length(sp_layers)],sep="_")
  
	  # For Benthic species : bottom layer
   ### MAEL GERNEZ, not needed
	  #  if( ALLSp[sp,"Depth_zone"] =="benthic" )      sp_layers <- "bottom"
	  # 
	  # # For Benthic species : bottom layer
	  #  if( ALLSp[sp,"Depth_zone"] =="surface" )      sp_layers <- "surface"

  # WorkProgress table : add species layers 
  WorkProgress[sp,"sp_layers"] <- sp_layers 

  # Study area to predict species distribution
  sp.dom <- as.character(((ALLSp$Domain))[sp])
  DOM.extent  <- DOMextents[[sp.dom]]

  ##########################################
  # 1 - COLLECT SPECIES OCCURRENCES DATA   #
  ##########################################
  

  ##### 1.1 Download species occurrenes at global scale from gloabl databases  : obis/gbif etc 
  
  OccStatus <- source(file = file.path(Maindir,"SCRIPTS/01-Species_occurrences.R") )

  # if there less than 50 occ OR no records in global databases : next species
  if( nrow(SpeciesOCC)[1] < 50 ){ rm("SpeciesOCC","OccStatus" ) ; next  } ; rm("OccStatus")

  # WorkProgress table : Add Occ nbr
  WorkProgress[sp,"Occurences"] <- nrow(SpeciesOCC)


  ############################################################################
  # 2 - CREATE BACKGROUND for each species : env space to create pseudo-abs  #
  ############################################################################


  source(file = file.path(Maindir,"SCRIPTS/02-Backgrounds.R"))

  # Change RESOLUTION OF BACKGROUND from 0.1 to 0.3 (fact = 3)

  if ( exists("fact") && fact>=0 ) {
	Sp_back_count.r <- aggregate (Sp_back_count.r, fact=fact, fun=sum,na.rm=T)
	Sp_back.r <- aggregate (Sp_back.r , fact=fact, fun=max,na.rm=T)
	Sp_back = values(Sp_back.r)
	Sp_back_count = values(Sp_back_count.r) }


  #################################################################
  # 3 - MATCH-UP between climate (S & T) and species occurrences  #
  #################################################################

  source(file = file.path(Maindir,"SCRIPTS/03-Matchup_OccClimat.R"))

  if( nrow(sp_occ_Clim) < 50  ){ rm("SpeciesOCC","sp_occ_Clim") ; next  }


 # WorkProgress table : Add nbr of Occ with Matchup
  WorkProgress[sp,"Match.up"] <- nrow(sp_occ_Clim)



  ######################################################################################
  # 4- GENERATE local rasters for historical predictions using Mercator climate data   #
  ######################################################################################

  # We predicted species ditribution for 9 depth layers : layers3
  source(file = file.path(Maindir,"SCRIPTS/04-ProjectionRasters.R"))

 
   #########################
   # 5 - NICHE MODELING    #
   #########################
 
   source(file = file.path(Maindir,"SCRIPTS/05-Nichemodels.R"))
 
   if( !file.exists(BiomodClibationFile) ){  unlink(tmpdir, recursive = TRUE) ; unlink(spnichfiles, recursive = TRUE)  ; rm("SpeciesOCC","sp_occ_Clim") ; next  }

 
   # WorkProgress table
   WorkProgress[sp,"Nichemodeling"] <- ifelse ( file.exists(BiomodClibationFile) ,"DONE","REDO")
 
 
  #################################
   # 6 - Historical predictions   #
   ################################


 source(file = file.path(Maindir,"SCRIPTS/06-Current2DPredictions.R"))

 # WorkProgress table
 WorkProgress[sp,"CurrentPredictions"] <- ifelse (length(dir(Niche_SpDomDirCurrent,recursive=T))>0 ,"DONE","REDO")


  rm("sp_occ_Clim", "Sp_back_count.r" , "Sp_back.r" , "Sp_back", "Sp_back_count","SpeciesOCC" )
  
} #eo species



