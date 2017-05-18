#### Figure 3, 4, 5 and S2: Bile acid analysis
##### Anna M. Seekatz
##### 4.27.17


##### Files needed:
	- directory: /Users/annaseekatz/Box Sync/Projects/UM.FMT/UMFMT_paper/16S_UMFMT
	- umfmt.allmeasures.txt (data in this file was merged in the prep code)
	
##### Figures:
	- Fig. 3, 4, and 5:
		- this code specifies levels of primary bile acids (Fig. 3), secondary bile acids (Fig. 4), or SCFAs (Fig. 5) over time by individual
		- the same code can be used to graph levels of any variable
	- Fig. S2
		- NMDS of bile acids (both primary and secondary) and SCFAs

#####  Fig. 3 (A-C): levels of bile acids over time

```
### Step 1: Set up the data

d<-read.table(file="umfmt.allmeasures.txt", header=TRUE)

data<-d
data$group<-as.numeric(gsub("D|R", "", data$patientID))
data$timepoint[data$timepoint==0 & data$unhomogenized.y=="yes"]<--4
data$timepoint[data$timepoint==0 & data$unhomogenized.y=="no"]<--3
data <- data[order(data$patientID, data$timepoint),]
df<-data[!is.na(data$CA), c(3, 4, 6, 35:69, 343)]			# subset for a slightly smaller df to work with--also, include only the samples that we have measurements for
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

### Step 2: create function:

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
		axis(1, at=unique(donors$timepoint), labels=c("unhomogenized", "homogenized"), cex.axis=0.8, tck=-0.03, las=2)
		abline(v=0, col = "grey80", lty=2)
		polygon(dcord.x,dcord.y,col='grey90', border=NA)
		polygon(rcord.x,rcord.y,col='grey90', border=NA)
		mapply(lines,xdonors,ydonors,col=dcol,pch=17,type="o", cex=1.2)
		mapply(lines,xrecipients,yrecipients,col=rcol,pch=19,type="o", cex=1.2)
		lines(unique(donors$timepoint), as.numeric(dmean), col="black", lwd=3)
		lines(r.axis, as.numeric(rmean), col="black", lwd=3)	#weirdness 2
		text(-3.5, max(data[,n]), "donors")
		text(-1.5, max(data[,n]), "pre-FMT")
		text(2.5, max(data[,n]), "post-FMT")
		}

### Step 3: Call out variables:

### For Fig. 3 (primary bile acids)
dev.new(width=3*4.669725, height=4.073394)
par(mfrow=c(1,3))
lapply(c("CA", "TCDCA", "TCA"), FUN=longiPlot)
title('primary bile acids', outer=TRUE, line=-2)

### For Fig. 4 (secondary bile acids)
dev.new(width=3*4.669725, height=4.073394)
par(mfrow=c(1,3))
lapply(c("DCA", "LCA", "UDCA"), FUN=longiPlot)
title('secondary bile acids', outer=TRUE, line=-2)

### For Fig. 5 (secondary bile acids)
dev.new(width=3*4.669725, height=4.073394)
par(mfrow=c(1,3))
lapply(c("Acetate_nmol", "Propionate_nmol", "Butyrate_nmol"), FUN=longiPlot)
title('short chain fatty acids', outer=TRUE, line=-2)

### For Fig. S3 (lactate and succinate)
dev.new(width=2*4.669725, height=4.073394)
par(mfrow=c(1,2))
lapply(c("lactate_uM.mg", "succinate_uM.mg"), FUN=longiPlot)

### Stats:
# donor samples:
df<-d[d$human_group %in% c("donor"), c(10,35:69)]

modelList<-list()
for(i in 2:36){
    fmla <- formula(paste(names(df)[i], " ~ unhomogenized.x"))
    modelList[[i]]<-wilcox.test(fmla, data = df, paired = FALSE)
}
modelList

# only these weresignificantly impacted by homogenization
	TCDCA W = 2, p-value = 0.03175
	TCA W = 2.5, p-value = 0.02537

# recipient samples (to donor and post/pre FMT)
data<-d[d$unhomogenized.x %in% c("yes"), c(4,35:69)]
df<-data[data$timepoint %in% c(0,-2),]

modelList<-list()
for(i in 2:36){
    fmla <- formula(paste(names(df)[i], " ~ timepoint"))
    modelList[[i]]<-wilcox.test(fmla, data = df, paired = FALSE)
}
modelList

### comparisons between donor and:
## -2 timepoint:
	succinate p = 0.01587
#valerate p-value = 0.007937
#LCA W = 3, p-value = 0.05556
#HDCA W = 3, p-value = 0.04491
	CA W = 24, p-value = 0.01587

## -1 timepoint:
	lactate W = 19, p-value = 0.03175
#valerate W = 4, p-value = 0.05195
#acetate W = 5, p-value = 0.08225
#LCA W = 5, p-value = 0.08214
#HDCA W = 4, p-value = 0.04437
	CA W = 27, p-value = 0.0303
#GDCA W = 3, p-value = 0.03534

## 2 timepoint:
#GHDCA W = 25.5, p-value = 0.04648

### between -1 timepoints in recipients and:
# 4:
#valerate W = 5, p-value = 0.04113
	acetate W = 4, p-value = 0.02597
#LCA W = 3, p-value = 0.06911
# 3:
#lactate W = 21, p-value = 0.06667
	butyrate W = 2, p-value = 0.008658
	acetate W = 2, p-value = 0.008658
	LCA W = 3, p-value = 0.02002
	CA W = 36, p-value = 0.002165
	TCDCA W = 33, p-value = 0.01515
# 2:
#acetate W = 6, p-value = 0.06494
	LCA W = 3, p-value = 0.02002
#CA W = 30, p-value = 0.06494
	TCDCA W = 31, p-value = 0.04113
# 1:
#valerate W = 1, p-value = 0.008658
#butyrate W = 4, p-value = 0.05195
	propionate W = 3, p-value = 0.0303
#acetate W = 4, p-value = 0.05195
#LCA W = 3, p-value = 0.06911

### between -2 timepoints in recipients and:
# 1:
#heptanoate W = 0, p-value = 0.007937
#octanoate W = 2, p-value = 0.03615
#hexanoate W = 0, p-value = 0.007937
#valerate W = 0, p-value = 0.007937
#butyrate W = 3, p-value = 0.05556
#TCA W = 18, p-value = 0.06169
# 2:
#LCA W = 4, p-value = 0.05195
#CA W = 26, p-value = 0.05195
	TCA W = 28, p-value = 0.02217
# 3:
#succinate W = 4, p-value = 0.05195
	butyrate W = 3, p-value = 0.0303
	CA W = 30, p-value = 0.004329
	TCDCA W = 27, p-value = 0.0303
	TCA W = 28, p-value = 0.02127
#butyrate W = 5, p-value = 0.08225
#acetate W = 5, p-value = 0.08225

```

#####  Fig. S2: NMDS/PCoA of metabolite data 

```

library(vegan)

# read in data:
d<-read.table(file="umfmt.allmeasures.txt", header=TRUE)
d<-d[order(d$human_group), ]					# data has to be ordered by the group for colors later
rownames(d)<-d$sampleID
dmatrix<-as.matrix(d[!is.na(d$CA), 35:59])		# eliminates samples with NAs, and only chooses BA axes

# use metaMDS, which uses NMDS, with standardized scaling (also adds the species to the axis)
# can specify several different types of distances here
d_NMDS=metaMDS(dmatrix, k=2, distance="horn")

# look at data:
stressplot(d_NMDS)
plot(d_NMDS)

# preliminary plot:	
ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",cex=1.25,air=0.01)

# Fig. S2: add group information:
summary(d$human_group)		# look at what levels we have
treat=c(rep("donor",10),rep("post",23),rep("pre",11))
fig <- ordiplot(d_NMDS, type = "none")
points(fig, "sites", pch=21, col="black", bg=c(rep("chartreuse4",10),rep("blue",20),rep("gold",11)))
orditorp(d_NMDS,display="sites",col=c(rep("chartreuse4",10),rep("blue",20),rep("gold",11)), air=0.01,cex=0.8, pos=1)
orditorp(d_NMDS,display="species",col="red",air=0.01, cex=0.5)

# run adonis to see if clustering is significant between groups
d<-read.table(file="umfmt.allmeasures.txt", header=TRUE)
d<-d[order(d$human_group), ]
rownames(d)<-d$sampleID
ddf<-as.data.frame(d[!is.na(d$CA), 35:59])
metadf<-as.data.frame(d[!is.na(d$CA), c(3,6)])
# use adonis:
adonis2(ddf ~ human_group, data = metadf, by = NULL)
#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999
#adonis2(formula = ddf ~ human_group, data = metadf, by = NULL)
#         Df SumOfSqs      F Pr(>F)   
#Model     2    1.711 2.9377  0.005 **
#Residual 38   11.066                 
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```