#	This program is to rescale the scores of assessors to make comparable scores.
#	This is for Matt White's project.

#	This program is based on the image assessment data.

#	This program is similar to score_scaling-5.r.
#	Their difference is that here we use linear mixed model (for incomplete block design)
#	rather than general linear model (for factorial design).

#	This is modified on 25 Aug 2021.

#	On 29 Aug 2021
#	This program is similar to score_scaling-6.r.
#	The difference is that here we select a new data point from a triangle distribution with the range the assessor provided and the most possible score as the mode. 

#	On 27 Apr 2022
#	This program is very similar to score_scaling-6-2.r.
#	The difference is that here we thought removing the outlier assessor is not necessary and the related parts were removed.

#	This program was written by Canran Liu.

#	source("C:/Data/U_copy/accudx/score_scaling-6-3.r")



set.seed(5510010)

###########################################################################################################################################

library(lme4)
library(triangle)

m <- 1		# dataset enlarger: select m values from the interval the assessor gives. Previously m = 5. But it is better to use m = 1.
nsm <- 1000	# times of simulations for calculating multipliers

dir0 <- "C:/Data/U_copy/stat_consult/white-5/"
dir <- "C:/Data/U_copy/stat_consult/white-5/rst1/"

fln_in <- paste(dir0,"ImageAssessment.csv",sep="")
fln_out <- paste(dir,"rst_additives_and_predictions_through_partial_response.csv",sep="")

dt0 <- read.csv(file=fln_in,header=T)

dt0$Assessor_name <- as.factor(dt0$Assessor_name)
dt0$trans_Overall_condition_score <- asin(sqrt(dt0$Overall_condition_score))


assessor <- unique(dt0$Assessor_name)

nass <- length(assessor)

md <- lmer(trans_Overall_condition_score ~ Assessor_name + (1 | ImageName), data=dt0)
summary(md)

cf <- fixef(md)

ad0 <- c(0,-cf[2:nass])				# difference from assessor 1's score in (arcsine squared root) transformed scale, here the difference for assessor 1 is 0.


#	generate many datasets by sampling using triangle distribution from the range of the original data for calculating the additives


n <- nrow(dt0)

ad <- matrix(0,nr=nsm,nc=nass)

for (j in 1:nsm) {

cat("j =",j,"\n")

dt1 <- dt0

k <- 0

for (i in 1:n) {

	rs <- rtriangle(m,dt0[i,"Condition.Rating.Lower"],dt0[i,"Condition.Rating.Upper"],dt0[i,"Overall_condition_score"])

	for (t in 1:m) {
		dt1[k+t,c(1:13,15:22)] <- dt0[i,c(1:13,15:22)]
	}

	dt1[k + 1:m,14] <- rs

	k <- k+m

}


#	Calculate the additives


dt1$trans_Overall_condition_score <- asin(sqrt(dt1$Overall_condition_score))


md <- lmer(trans_Overall_condition_score ~ Assessor_name + (1 | ImageName), data=dt1)


cf <- fixef(md)

ad[j,] <- c(0,-cf[2:nass])				# difference from assessor 1's score in transformed scale

}	# end of "for (j in 1:nsm)"

ad_mean <- apply(ad,2,mean)
ad_sd <- apply(ad,2,sd)
ad_median <- apply(ad,2,median)

qt025 <- function(x) {
	return(quantile(x,0.025))
}

qt975 <- function(x) {
	return(quantile(x,0.975))
}

ad_025 <- apply(ad,2,qt025)
ad_975 <- apply(ad,2,qt975)

adtv <- data.frame(Assessor_name=sort(assessor),ad_mean=ad_mean,ad_sd=ad_sd,ad_median=ad_median,ad_025=ad_025,ad_975=ad_975)

write.csv(adtv,fln_out,row.names=FALSE)


#----------------------------------------------------------------------------------------------------------------------------

#	generate another dataset from the original data for evaluating the additives 

library(lme4)
library(triangle)

m <- 1		# dataset enlarger: select m values from the interval the assessor gives. 
nsm <- 1000	# times of simulations for calculating additives

dir0 <- "C:/Data/U_copy/stat_consult/white-5/"
dir <- "C:/Data/U_copy/stat_consult/white-5/rst1/"

fln_in <- paste(dir0,"ImageAssessment.csv",sep="")
fln <- paste(dir,"rst_additives_and_predictions_through_partial_response.csv",sep="")
fln_fig <- paste(dir,"fig_comparison_between_rescaled_and_unrescaled.tiff",sep="")	

rst0 <- read.csv(file=fln,header=T)

dt0 <- read.csv(file=fln_in,header=T)
dt0$Assessor_name <- as.factor(dt0$Assessor_name)
dt0$trans_Overall_condition_score <- asin(sqrt(dt0$Overall_condition_score))

n <- nrow(dt0)

assessor <- unique(dt0$Assessor_name)

nass <- length(assessor)

image <- unique(dt0$ImageName)

nimg <- length(image)

prd_unscl <- matrix(0,nr=nsm,nc=nass)
mad_unscl <- numeric(nsm)
sd_unscl <- numeric(nsm)

prd_scl <- matrix(0,nr=nsm,nc=nass)
mad_scl <- numeric(nsm)
sd_scl <- numeric(nsm)

for (j in 1:nsm) {

cat("j =",j,"\n")

dt2 <- dt0

k <- 0

for (i in 1:n) {

	rs <- rtriangle(m,dt0[i,"Condition.Rating.Lower"],dt0[i,"Condition.Rating.Upper"],dt0[i,"Overall_condition_score"])

	for (t in 1:m) {
		dt2[k+t,c(1:13,15:22)] <- dt0[i,c(1:13,15:22)]
	}

	dt2[k + 1:m,14] <- rs

	k <- k+m

}


#	Calculate the partial response with the unrescaled scores


dt2$trans_Overall_condition_score <- asin(sqrt(dt2$Overall_condition_score))


md <- lmer(trans_Overall_condition_score ~ Assessor_name + (1 | ImageName), data=dt2)

cf <- fixef(md)
prd_unscl[j,] <- cf + c(0,rep(cf[1],nass-1))

mad_unscl[j] <- mad(prd_unscl[j,])
sd_unscl[j] <- sd(prd_unscl[j,])


#	Calculate the partial response with the rescaled scores

dt3 <- merge(dt2,rst0,by="Assessor_name")

dt3$calibrated_trans_Overall_condition_score <-  asin(sqrt(dt3$Overall_condition_score)) + dt3$ad_mean

md <- lmer(calibrated_trans_Overall_condition_score ~ Assessor_name + (1 | ImageName), data=dt3)

cf <- fixef(md)
prd_scl[j,] <- cf + c(0,rep(cf[1],nass-1))

mad_scl[j] <- mad(prd_scl[j,])
sd_scl[j] <- sd(prd_scl[j,])


}	# end of "for (j in 1:nsm)"



qt025 <- function(x) {
	return(quantile(x,0.025))
}

qt975 <- function(x) {
	return(quantile(x,0.975))
}


prd_unscl <- (sin(prd_unscl))^2
prd_scl <- (sin(prd_scl))^2

prd_unscl_050 <- apply(prd_unscl,2,median)
prd_unscl_025 <- apply(prd_unscl,2,qt025)
prd_unscl_975 <- apply(prd_unscl,2,qt975)


prd_scl_050 <- apply(prd_scl,2,median)
prd_scl_025 <- apply(prd_scl,2,qt025)
prd_scl_975 <- apply(prd_scl,2,qt975)


rst0$prd_unscl_050 <- prd_unscl_050
rst0$prd_unscl_025 <- prd_unscl_025
rst0$prd_unscl_975 <- prd_unscl_975


rst0$prd_scl_050 <- prd_scl_050
rst0$prd_scl_025 <- prd_scl_025
rst0$prd_scl_975 <- prd_scl_975


write.csv(rst0,fln,row.names=FALSE)

#	make a single plot for comparison between rescaled and unrescaled

tiff(file=fln_fig,width=8,height=20,units="cm",res=500,pointsize=10,compression="lzw",antialias = "cleartype")

mn <- min(min(prd_unscl),min(prd_scl))
mx <- max(max(prd_unscl),max(prd_scl))

par(mfrow=c(3,1),omi=c(0.1,0.2,0.1,0.1))

cx <- 1.2
cx1 <- 1

		par(mfg=c(1,1,3,1))
		par(cex.axis=cx1,cex.lab=cx,cex.main=cx,cex.sub=cx1)


plot(1:nass,prd_unscl_050,xlab="Assessor",ylab="Predicted score",main="(a) Unrescaled",ylim=c(mn,mx),type="p",pch=19)
points(1:nass,prd_unscl_025,pch=1)
points(1:nass,prd_unscl_975,pch=1)


		par(mfg=c(2,1,3,1))
		par(cex.axis=cx1,cex.lab=cx,cex.main=cx,cex.sub=cx1)


plot(1:nass,prd_scl_050,xlab="Assessor",ylab="Predicted score",main="(b) Rescaled",ylim=c(mn,mx),type="p",pch=19)
points(1:nass,prd_scl_025,pch=1)
points(1:nass,prd_scl_975,pch=1)


		par(mfg=c(3,1,3,1))
		par(cex.axis=cx1,cex.lab=cx,cex.main=cx,cex.sub=cx1)


dt5 <- data.frame(scaling=c(rep("Rescaled",nsm),rep("Unrescaled",nsm)),MAD=c(mad_scl,mad_unscl))
dt5$scaling <- as.factor(dt5$scaling)
dt5$scaling <- relevel(dt5$scaling, ref="Unrescaled")

boxplot(dt5$MAD ~ dt5$scaling,xlab="",ylab="MAD",main="(c) MAD comparison")

dev.off()


#----------------------------------------------------------------------------------------------------------------------------

#	display the effectiveness of the rescaling method using some images as examples


dir0 <- "C:/Data/U_copy/stat_consult/white-5/"
dir <- "C:/Data/U_copy/stat_consult/white-5/rst1/"

fln_in <- paste(dir0,"ImageAssessment.csv",sep="")
fln_in1 <- paste(dir,"rst_additives_and_predictions_through_partial_response.csv",sep="")
fln_fig <- paste(dir,"fig_comparison_between_rescaled_and_unrescaled_with_two_examples.tiff",sep="")	

dt0 <- read.csv(file=fln_in,header=T)
dt0$Assessor_name <- as.factor(dt0$Assessor_name)
add <- read.csv(file=fln_in1,header=T)
add$Assessor_name <- as.factor(add$Assessor_name)

dt0$trans_Overall_condition_score <- asin(sqrt(dt0$Overall_condition_score))

images <- unique(dt0$ImageName)
nimg <- length(images)

for (i in 1:nimg) { 
	dt1 <- dt0[dt0$ImageName == images[i],]
	cat(i,nrow(dt1),"\n")
}

#	select images 70 and 74 for examples

#	mainly use Assessor_name and Overall_condition_score

tiff(file=fln_fig,width=20,height=10,units="cm",res=500,pointsize=10,compression="lzw",antialias = "cleartype")

par(mfcol=c(1,2),omi=c(0.1,0.2,0.1,0.1))

cx <- 1.2
cx1 <- 1

k <- 0

for (image_id in c(70,74)) {

k <- k+1

par(mfg=c(1,k,1,2))
par(cex.axis=cx1,cex.lab=cx,cex.main=cx,cex.sub=cx1)

dt1 <- dt0[dt0$ImageName == images[image_id],c("Assessor_name","Overall_condition_score")]

dt1$re_score <- dt1$Overall_condition_score	# initialise a new column

assessors <- unique(dt1$Assessor_name)

for (i in 1:nrow(dt1)) {
	dt1$re_score[i] <- (sin(asin(sqrt(dt1$Overall_condition_score[i])) + add[c(1:nrow(add))[add$Assessor_name %in% dt1$Assessor_name[i]],"ad_mean"]))^2		# if do not use those from simulation, just remove "_sm"
}

id = order(as.numeric(assessors))

plot(1:nrow(dt1),dt1[id,2],xaxt="n",xlab="Assessor",ylab="Condition score",main="",ylim=c(0,1),type="p",pch=1,cex.lab=1)
axis(side=1,at=1:nrow(dt1),labels=assessors[id])
points(1:nrow(dt1),dt1[id,3],pch=19)

}

dev.off()


###################################################################################################################################

#	rescale the assessment data


dir0 <- "C:/Data/U_copy/stat_consult/white-5/"
dir1 <- "C:/Data/U_copy/stat_consult/white-5/rst1/"

fln_in0 <- paste(dir0,"SiteConditionAssessment.csv",sep="")
fln_in1 <- paste(dir1,"rst_additives_and_predictions_through_partial_response.csv",sep="")
fln_out <- paste(dir1,"rescaled_SiteConditionAssessment.csv",sep="")

dt0 <- read.csv(file=fln_in0,header=T)
dt0$Assessor_name <- as.character(dt0$Assessor_name)
add <- read.csv(file=fln_in1,header=T)
add$Assessor_name <- as.character(add$Assessor_name)

dt0$additive.recaled.Overall.Condition.Rating <- NULL

for (i in 1:nrow(dt0)) {
	additive <- add[add$Assessor_name == dt0$Assessor_name[i],"ad_mean"]
	dt0$additive.recaled.Overall.Condition.Rating[i] <- (sin(asin(sqrt(dt0$Overall.Condition.Rating[i])) + additive))^2
	
}

write.csv(dt0,fln_out,row.names=FALSE)



#---------------------------------------------------------------------------------------------------------------------------------

#	make histogram for both raw condition scores and rescaled condition scores

cx1 <- 1
cx <- 1

wd <- 15
ht <- 8

dir1 <- "C:/Data/U_copy/stat_consult/white-5/rst1/"

fln_in <- paste(dir1,"rescaled_SiteConditionAssessment.csv",sep="")
fln_out <- paste(dir1,"histgraph_for_raw_and_rescaled_condition_scores.tiff",sep="")

dt0 <- read.csv(file=fln_in,header=T)

dt <- dt0[,c("Overall.Condition.Rating","additive.recaled.Overall.Condition.Rating")]

fln_out <- paste(dir1,"histgraph_for_raw_and_rescaled_condition_scores.tiff",sep="")

tiff(file=fln_out,width=wd,height=ht,units="cm",res=500,pointsize=10,compression="lzw",antialias = "cleartype")

par(mfcol=c(1,2),omi=c(0.5,0.5,0.1,0.1))

	par(mfg=c(1,1,1,2))	# put the current figure to position (k1,k2)
	par(mai=c(0.2,0.3,0.2,0.1))
	par(adj=0.5,cex=cx1,tck=-0.02,mgp=c(2,0.5,0),lab=c(5,5,5))

hist(dt[,1],breaks=c(0,1:10/10),main="(a) Original",xlab="Condition score",xaxp=c(0,1,10),yaxp=c(0,70,7),ylim=c(0,78),labels=TRUE)

	par(mfg=c(1,2,1,2))	# put the current figure to position (k1,k2)
	par(mai=c(0.2,0.3,0.2,0.1))
	par(adj=0.5,cex=cx1,tck=-0.02,mgp=c(2,0.5,0),lab=c(5,5,5))

hist(dt[,2],breaks=c(0,1:10/10),main="(b) Rescaled",xlab="Condition score",xaxp=c(0,1,10),yaxp=c(0,70,7),ylim=c(0,78),labels=TRUE)

mtext("Condition score",side=1,outer=T,at=0.5,line=1,cex=cx)
mtext("Frequency",side=2,outer=T,at=0.51,line=1,cex=cx)
#mtext("Original",side=3,outer=T,at=0.3,line=0.5,cex=cx)
#mtext("Rescaled",side=3,outer=T,at=0.75,line=0.5,cex=cx)

dev.off()

#---------------------------------------------------------------------------------------------------------------------------------

#	make a point plot for both raw conditionscores and rescaled condition scores

cx1 <- 1
cx <- 1

wd <- 15
ht <- 8

dir1 <- "C:/Data/U_copy/stat_consult/white-5/rst1/"

fln_in <- paste(dir1,"rescaled_SiteConditionAssessment.csv",sep="")
fln_out <- paste(dir1,"point_plot_for_raw_and_rescaled_condition_scores.tiff",sep="")

dt0 <- read.csv(file=fln_in,header=T)
dt0$Assessor_name <- as.factor(dt0$Assessor_name)

dt <- dt0[,c("Overall.Condition.Rating","additive.recaled.Overall.Condition.Rating")]

id <- order(dt$Overall.Condition.Rating)


tiff(file=fln_out,width=wd,height=ht,units="cm",res=500,pointsize=10,compression="lzw",antialias = "cleartype")

plot(1:nrow(dt),dt[id,1],xlab="Data point identifier",ylab="Condition score",pch=46)
points(1:nrow(dt),dt[id,2],pch=46,col="blue")

dev.off()

#---------------------------------------------------------------------------------------------------------------------------------

#	model diagnostics

library(lme4)
library(car)

dir0 <- "C:/Data/U_copy/stat_consult/white-5/"
dir1 <- "C:/Data/U_copy/stat_consult/white-5/rst1/"

fln_in <- paste(dir0,"ImageAssessment.csv",sep="")

fln_out <- paste(dir1,"model_diagnostics.tiff",sep="")


dt0 <- read.csv(file=fln_in,header=T)

dt0$Assessor_name <- as.factor(dt0$Assessor_name)
dt0$trans_Overall_condition_score <- asin(sqrt(dt0$Overall_condition_score))


assessor <- unique(dt0$Assessor_name)

nass <- length(assessor)

md <- lmer(trans_Overall_condition_score ~ Assessor_name + (1 | ImageName), data=dt0)
summary(md)

rsd = residuals(md)
fit = fitted(md)


tiff(file=fln_out,width=10,height=20,units="cm",res=500,pointsize=10,compression="lzw",antialias = "cleartype")

par(mfrow=c(3,1),omi=c(0.1,0.2,0.1,0.1))

cx <- 1.2
cx1 <- 1

		par(mfg=c(1,1,3,1))
		par(cex.axis=cx1,cex.lab=cx,cex.main=cx,cex.sub=cx1)

boxplot(rsd ~ dt0$Assessor_name,xlab="Assessor",ylab="Residual",main="(a)")

		par(mfg=c(2,1,3,1))
		par(cex.axis=cx1,cex.lab=cx,cex.main=cx,cex.sub=cx1)

plot(fit,rsd,xlab="Predicted value",ylab="Residual",main="(b)")

		par(mfg=c(3,1,3,1))
		par(cex.axis=cx1,cex.lab=cx,cex.main=cx,cex.sub=cx1)
qqPlot(rsd,envelope=list(style="lines"),col.lines="black",lwd=1,grid=FALSE,id=FALSE,xlab="Normal quantiles",ylab="Residual",main="(c)")

dev.off()

	