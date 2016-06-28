# Bossted_regression_trees.R
	# creates function to run boosted regression tree from extracted dataframe
	# uses multi-core processing

t1=Sys.time()
library(doParallel)
library(foreach)

#function to run boosted regression tree
boost.rt <- function(si,ei,model){



library(dismo)
library(gbm)
library(dplyr)

#number of sampleing points for every hexagon
num_pts <- c(12365,6372,12659,26496,26963,19305,20285,22788,25720,23333,18448,14914,8204,22739,21308,19081)
num_pts_fire = c(6754,4074,10560,6067,9822,1888,27287,26447,13756,20580,3705,24293,2989,45383,27252,14354)
per <- c(.03,.03,.02,.02,.02,.06,.02,.02,.02,.02,.03,.02,.04,.02,.02,.02)
in.dir <- "W:/EDM/Fire/DeLancey/North_America/ANALYSIS_redo_May26/data_frame_RESULTS/"
importance.out <- "W:/EDM/Fire/DeLancey/North_America/ANALYSIS_redo_May26/importance_plots/"
response.out <- "W:/EDM/Fire/DeLancey/North_America/ANALYSIS_redo_May26/response_curves/"

#different variables for every model
if (model == "multiple.anthro"){
	vars <- c(1,5,7,16,18,19,20,28,32,39,40,42,47)
	x <- 2:13
	y <- 1
}
if ( model == "HF"){
	vars <- c(1,5,7,16,18,19,20,28,39,42,47)
	x <- 2:11
	y <- 1
	va <- 4
}
if ( model == "pop"){
	vars <- c(1,5,7,16,19,20,28,32,39,42,47)
	x <- 2:11
	y <- 1
	va <- 7
}
if ( model == "rlv"){
	vars <- c(1,5,7,16,19,20,28,39,40,42,47)
	x <- 2:11
	y <- 1
	va <- 8
	
}



for (i in si:ei){

	# read in data for each hexel
	d <- read.csv(paste0(in.dir,"hex",i,".csv"))
	d <- na.omit(d)
	d <- d[,vars]
	df.importance <- data.frame(1:(length(vars)-1))
	df.plt <- data.frame(1:100)
	df.AUC <- data.frame(c(1))
	df.dev <- data.frame(c(1))
	AUC <- vector()
	dev <- vector()
	
	# subsampling loop done for 50 replicates
	for (j in 1:50){
		d1 <- filter(d, fire==1)
		samp1 <- sample(1:length(d1[,1]),round(num_pts_fire[i]*per[i]))
		d1 <- d1[samp1,]
		d0 <- filter(d, fire==0)
		samp <- sample(1:length(d0[,1]),length(d1[,1])*2)
		d0 <- d0[samp,]
		d.df <- rbind(d1,d0)
		
		#run boosted regression tree
		brt <- gbm.step(data=d.df, gbm.x=x ,gbm.y=y, family = "bernoulli", tree.complexity = 3, learning.rate = 0.005, bag.fraction = 0.5)
		# stats
		AUC[j] <- brt$cv.statistics$discrimination.mean
		dev[j] <- (brt$self.statistics$mean.null - brt$self.statistics$mean.resid) / brt$self.statistics$mean.null
		if (model == "multiple.anthro"){
			imp <- data.frame(summary(brt))
			imp <- arrange(imp,var)
			imp <- imp[,2]
			df.importance <- cbind(df.importance,imp)
		}else{
			plt <- plot.gbm(brt,return.grid=T,i.var=va)
			df.plt <- cbind(df.plt,plt[,2])
		}
	}
	AUC <- mean(AUC)
	dev <- mean(dev)
	stats <- cbind(AUC,dev)
	
	#output statistics
	if (model == "multiple.anthro"){
		write.csv(stats,paste0(importance.out,"hex",i,"stats.csv"),row.names=F)
		df.importance <- df.importance[,-1]
		im <- rowMeans(df.importance)
		comb <- c((im[1]+im[2]+im[6]+im[9]),(im[3]+im[7]+im[11]+im[12]),im[5],(im[4]+im[8]+im[10]))
	
		png(filename=paste0(importance.out,"hex",i,".png"),width=900,height=800,bg="transparent")
		barplot(comb,names.arg=c("Climate","Enduring","Human","Lightning"),cex.axis=4,cex.names=4,col="black",ylab="Variable group importance",cex.lab=4,ylim=c(0,90),space=c(0.05,0.05,0.05),las=2)
		par(mar=c(8,0,0,0))
		barplot(comb,names.arg=c("","","",""),cex.axis=8,cex.names=10,col="black",ylab="",cex.lab=10,ylim=c(0,90),space=c(0.05,0.05,0.05),las=2,axes=FALSE)
		axis(side=1,at=c(0.53,1.6,2.65,3.72),labels=c("C","E","L","A"),col="transparent",cex.axis=11,line=6,lwd=4,font.axis=2)
		dev.off()
	} else {
		write.csv(stats,paste0(response.out,model,"/","hex",i,"stats.csv"),row.names=F)
		df.plt <- df.plt[,-1]
		chrt <- rowMeans(df.plt)
		std <- apply(df.plt,1,sd)
		chrt.min <- chrt - std
		chrt.pls <- chrt + std
		chrt <- cbind(plt[,1],chrt,chrt.min,chrt.pls)
		write.csv(chrt,paste0(response.out,model,"/","hex",i,".csv"),row.names=F)
		
	}
}
}

#set up multi-core processing
cl <- makeCluster(16)
registerDoParallel(cl)


#loop through 16 iterations fo the boost.rt function with one iteration on each core
foreach(i=1:16) %dopar% {
	boost.rt(i,i,"rlv")
}


t2=Sys.time()
t2-t1