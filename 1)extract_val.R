# 1)extract_val.R
	# extracts values from raster grids for location across north america

t1=Sys.time()
library(doParallel)
library(foreach)

#function to extrac values from raster grids
extract.val <- function(si,ei){



library(raster)
library(rgdal)

#list of grids
ls <- list.files("W:/EDM/Fire/DeLancey/North_America/ANALYSIS_redo_May26/input_var_tifs", pattern=".tif$")
in.dir <- "W:/EDM/Fire/DeLancey/North_America/ANALYSIS_redo_May26/input_var_tifs/"

#number of point in hexels
num_pts = c(12365,6372,12659,26496,26963,19305,20285,22788,25720,23333,18448,14914,8204,22739,21308,19081)
num_pts_fire = c(6754,4074,10560,6067,9822,1888,27287,26447,13756,20580,3705,24293,2989,45383,27252,14354)
num_pts <- num_pts + num_pts_fire
names <- c("AHM", "aspect", "bFFP", "CMD", "DD_0", "DD_18")
out.dir <- "W:/EDM/Fire/DeLancey/North_America/ANALYSIS_redo_May26/data_frame_RESULTS/"
hex.dir <- "W:/EDM/Fire/DeLancey/North_America/ANALYSIS_redo_May26/hexes"
pts.dir <- "W:/EDM/Fire/DeLancey/North_America/ANALYSIS_redo_May26/pts"

#loop through all hexels
for (i in si:ei){
	dat <- data.frame(row=1:num_pts[i])
	area <- raster(paste0(in.dir,i,".tif"))
	pts <- readOGR(pts.dir, paste0("pts",i))
	ext.fire <- extract(area,pts)
	dat <- cbind(dat,ext.fire)
	#loop through all grids
	for (j in 17:62){
		r <- raster(paste0(in.dir,ls[j]))
		ext <- extract(r,pts)
		dat <- cbind(dat,ext)
		print(j)
	}
	dat <- dat[,-1]
	colnames(dat) <- c("fire",ls[17:62])
	write.csv(dat,paste0(out.dir,"hex",i,".csv"),row.names=F)
}
}

#set up multi-core environment
cl <- makeCluster(16)
registerDoParallel(cl)


#loop through 16 iterations of the extract.val function
foreach(i=1:16) %dopar% {
	extract.val(i,i)
}


t2=Sys.time()
t2-t1