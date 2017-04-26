#### Figure 1B: Diversity over time
##### Anna M. Seekatz
##### 4.12.17

##### Files needed:
	- directory: /Users/annaseekatz/Box Sync/Projects/UM.FMT/UMFMT_paper/16S_UMFMT
	- umfmt.allmeasures.txt
		- created previously in prep code
	
```

d<-read.table(file="umfmt.allmeasures.txt", header=TRUE)

data<-d
data$group<-as.numeric(gsub("D|R", "", data$patientID))
data$timepoint[data$timepoint==0 & data$unhomogenized.y=="yes"]<--4
data$timepoint[data$timepoint==0 & data$unhomogenized.y=="no"]<--3
data <- data[order(data$patientID, data$timepoint),]
df<-data[!is.na(data$nseqs), c(3, 4, 6, 13:34, 343)]			# subset for a slightly smaller df to work with
data<-droplevels(df)
data <- data[order(data$patientID, data$timepoint),]
donors<-data[which(data$human_group=="donor"), ]
recipients<-data[which(data$human_group %in% c("pre", "post")), ]
donors<-droplevels(donors)
recipients<-droplevels(recipients)
r.axis <- c(-2,-1,1,2,3,4)	# can't figure out why this is so weird, so just redefined it to include all time points for recipients

# color codes:
colors<-colorRampPalette(c("grey20","grey80"))(n = 6)
dcol<-colorRampPalette(c("chartreuse","darkgreen"))(n = 5)
rcol<-colorRampPalette(c("cadetblue1", "deepskyblue", "blue4"))(n = 6)

# graphing it (let's use a couple different diversity measures):
longiPlot<- function(n) {
	xdonors <- tapply(donors$timepoint,donors$group,function(x) return(x))
	ydonors <- tapply(donors[,n],donors$group,function(x) return(x))
	xrecipients <- tapply(recipients$timepoint,recipients$group,function(x) return(x))
	yrecipients <- tapply(recipients[,n],recipients$group,function(x) return(x))
	dmean <- tapply(donors[,n],donors$timepoint,function(x) mean(x))
	rmean <- tapply(recipients[,n],recipients$timepoint,function(x) mean(x))
	dsd <- tapply(donors[,n],donors$timepoint,function(x) sd(x)/sqrt(length(x)))
	rsd <- tapply(recipients[,n],recipients$timepoint,function(x) sd(x)/sqrt(length(x)))
	dhigh<-dmean+dsd
	rhigh<-rmean+rsd
	dlow<-dmean-dsd
	rlow<-rmean-rsd
	dcord.x<-c(unique(donors$timepoint), rev(unique(donors$timepoint)))
	dcord.y<-c(dhigh, rev(dlow))
	rcord.x<-c(r.axis, rev(r.axis))	#weirdness 1
	rcord.y<-c(rhigh, rev(rlow))
		plot(data$timepoint,data[,n],type="n", xlim=c(-4,4), ylim=c(0,max(data[,n], na.rm=TRUE)), ylab=names(data[n]), xaxt='n', xlab="sampling")
		axis(1, at=recipients$timepoint, labels=recipients$timepoint, cex.axis=0.8, tck=-0.03)
		axis(1, at=unique(donors$timepoint), labels=c("unprocessed", "processed"), cex.axis=0.8, tck=-0.03, las=2)
		abline(v=0, col = "grey80", lty=2)
		polygon(dcord.x,dcord.y,col='grey90', border=NA)
		polygon(rcord.x,rcord.y,col='grey90', border=NA)
		mapply(lines,xdonors,ydonors,col=dcol,pch=17,type="o", cex=0.8)
		mapply(lines,xrecipients,yrecipients,col=rcol,pch=19,type="o", cex=0.8)
		lines(unique(donors$timepoint), as.numeric(dmean), col="black", lwd=3)
		lines(r.axis, as.numeric(rmean), col="black", lwd=3)	#weirdness 2
		text(-3.5, max(data[,n]), "donors")
		text(-1.5, max(data[,n]), "pre-FMT")
		text(2.5, max(data[,n]), "post-FMT")
		}
par(mfrow=c(1,3))
lapply(c("invsimpson_03", "shannon_03", "sobs_03"), FUN=longiPlot)
title('Diversity measures', outer=TRUE, line=-2)

### For Fig. 1B (inverse Simpson)
lapply(c("invsimpson_03"), FUN=longiPlot)

# note: can change dot sizes above
# adjust graph size:
dev.size()
#[1] 4.669725 4.073394
dev.new(width=4.669725, height=4.073394)
lapply(c("invsimpson_03"), FUN=longiPlot)


### stats:

# comparison of unprocessed v. processed:
donors <- data[data$human_group=="donor", ]
wilcox.test(invsimpson_03~timepoint, data=donors)
#W = 14, p-value = 0.8413

# comparison of pre vs post (combined):




```