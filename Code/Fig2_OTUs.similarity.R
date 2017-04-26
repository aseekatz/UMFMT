#### Figure 2: Donor/recipient similarity
##### Anna M. Seekatz
##### 1.26.17


##### Files needed to calculate number of shared OTUs and classifications:
	- Fig. 2B: community simmilarity intra/inter and over time
		- allhumos.betasummary.txt
	- Fig. 2C and Table S2: shared OTUs across donor recipient pairs
		- umfmt.allmeasures.txt
		- allhumos4.taxonomy.names.txt


##### Fig. 2A: intra and inter-individual community similarity
	
```
# read in files:
m<-read.table(file="allhumos.betasummary.txt", header=TRUE)

## need to separate groups into 4 different categories of pairwise comparisons:
	- intra-individual, post-FMT
	- intra-individual, pre/post-FMT
	- inter-individual, post-FMT
	- inter-individual, post-FMT to donors

## set up comparisons:
	# intra-individual, post-FMT
intra.post<-m[m$s1_patientID==m$s2_patientID & m$s1_human_group %in% c("post") & m$s2_human_group %in% c("post"), ]
intra.post$COMP<-"intra-post"
	# intra-individual, pre/post-FMT
intra.pp<-m[m$s1_patientID==m$s2_patientID & m$s1_human_group %in% c("post") & m$s2_human_group %in% c("pre") | m$s1_patientID==m$s2_patientID & m$s2_human_group %in% c("post") & m$s1_human_group %in% c("pre"), ]
intra.pp$COMP<-"intra-pre-post"
	# inter-individual, pre/post-FMT
inter.post<-m[!m$s1_patientID==m$s2_patientID & m$s1_human_group %in% c("post") & m$s2_human_group %in% c("post"), ]
inter.post$COMP<-"inter-post"
	# inter-individual, post-FMT to donors
	# note: need to remove homogenized samples first
m2<-m[m$s1_unhomogenized %in% c("yes") & m$s2_unhomogenized %in% c("yes"), ]
inter.pd<-m2[!m2$s1_patientID==m2$s2_patientID & m2$s1_human_group %in% c("donor") & m2$s2_human_group %in% c("post") | !m2$s1_patientID==m2$s2_patientID & m2$s2_human_group %in% c("donor") & m2$s1_human_group %in% c("post"), ]
inter.pd$COMP<-"inter-post-donor"

## merge all together:
mcomp<-rbind(intra.post, intra.pp, inter.post, inter.pd)

## graph results
	#color scheme:
outside.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
		for (i in 1:length(n)) {
		colorvec[i] = "light grey"
		if ( n[i] %in% c("intra-post", "intra-pre-post") ) {
		colorvec[i] = "grey45"
		}
		if ( n[i] %in% c("inter-post", "inter-post-donor") ) {
		colorvec[i] = "lightcyan2"
		}
	}
	c(colorvec)
}
inside.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
		for (i in 1:length(n)) {
		colorvec[i] = "light grey"
		if ( n[i] %in% c("intra-post", "intra-pre-post") ) {
		colorvec[i] = "grey80"
		}
		if ( n[i] %in% c("inter-post", "inter-post-donor") ) {
		colorvec[i] = "lightcyan"
		}
	}
	c(colorvec)
}
	# plot:
mcomp$COMP <- factor(mcomp$COMP,levels=c("intra-post", "intra-pre-post", "inter-post", "inter-post-donor"))
df<-mcomp
	#singular:
plot<-plot(as.numeric(thetayc) ~ as.factor(COMP), data = df, ylab=expression(paste("", theta, "YC (dissimilarity)")), xlab="comparison", ylim=c(0,1), xlab="", xaxt="n", outline=F)
points(thetayc~ jitter(as.numeric(COMP, factor=0)), data = df, bg=inside.col(df$COMP), col=outside.col(df$COMP), pch=21)
names<-levels(mcomp$COMP)
text(x =  seq(1,4,by=1), y = par("usr")[3]-0.05, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)

	# compare all distances:
simplot<- function (n) {
	plot(df$COMP, df[,n], ylab=names(df[n]), xlab="comparison", xaxt="n", outline=F)
	points(jitter(as.numeric(df$COMP, factor=0)), df[,n], data = df, bg=inside.col(df$COMP), col=outside.col(df$COMP), pch=21)
	names<-levels(mcomp$COMP)
	text(x =  seq(1,4,by=1), y = par("usr")[3]-0.05, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
	}
distances<-as.character(colnames(df)[19:26])
par(mfrow=c(2, 4))
lapply(distances, FUN=simplot)

```



##### Fig. 2B: comparing shared distance between respective and random donors to recipients
	
```
# read in files:
m<-read.table(file="allhumos.betasummary.txt", header=TRUE)

## get comparison between respective donors/recipients:
# we want to compare the sharedsobs between recipient samples over time and their respective donor
# AND, the recipient samples to random donor over time
# so, only want comparisons with a donor in one of the slots

	# let's get only pairwise comparisons between respective donor and unprocessed stool
dr1<-m[m$s1_patientID %in% c("D1") & m$s2_patientID %in% c("R1") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes")| m$s2_patientID %in% c("D1") & m$s1_patientID %in% c("R1") & m$s1_unhomogenized %in% c("yes") & m$s2_unhomogenized %in% c("yes"), ]
dr1$COMP<-"DR1"
dr1$TIME<-dr1$s2_timepoint

dr2<-m[m$s1_patientID %in% c("D2") & m$s2_patientID %in% c("R2") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes")| m$s2_patientID %in% c("D2") & m$s1_patientID %in% c("R2") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes"), ]
dr2$COMP<-"DR2"
dr2$TIME<-dr2$s2_timepoint

dr3<-m[m$s1_patientID %in% c("D3") & m$s2_patientID %in% c("R3") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes") | m$s2_patientID %in% c("D3") & m$s1_patientID %in% c("R3") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes"), ]
dr3$COMP<-"DR3"
dr3$TIME<-dr3$s2_timepoint

dr4<-m[m$s1_patientID %in% c("D4") & m$s2_patientID %in% c("R4") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes") | m$s2_patientID %in% c("D4") & m$s1_patientID %in% c("R4") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes"), ]
dr4$COMP<-"DR4"
dr4$TIME<-dr4$s2_timepoint

dr5<-m[m$s1_patientID %in% c("D5") & m$s2_patientID %in% c("R5") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes") | m$s2_patientID %in% c("D5") & m$s1_patientID %in% c("R5") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes"), ]
dr5$COMP<-"DR5"
dr5$TIME<-dr5$s2_timepoint

	# donors to their own (hom/unhom):
donors<-m[m$s1_patientID==m$s2_patientID & m$s1_human_group %in% c("donor"), ]
donors$COMP<-"DONORS"
donors$TIME<--4

alldr<-rbind(dr1, dr2, dr3, dr4, dr5, donors)
	# this file has the comparison (COMP) at each time point
	
# get comparison between random donors/recipients:
	# note: this only works because all the D comparisons happen to be on one side...
#random<-m[!m$s1_patientID==m$s2_patientID & m$s1_human_group %in% c("donor") & m$s2_human_group %in% c("pre", "post"), ]
#random$COMP[random$s1_patientID==random$s2_patientID & random$s1_human_group %in% c("donor") ]<-"DONORS"
	# at the time, could not figure out how to do this nicer...
random1<-m[m$s1_patientID %in% c("D1") & m$s2_patientID %in% c("R2", "R3", "R4", "R5", "R6") & m$s1_human_group %in% c("donor") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes"), ]	
random2<-m[m$s1_patientID %in% c("D2") & m$s2_patientID %in% c("R1", "R3", "R4", "R5", "R6") & m$s1_human_group %in% c("donor") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes"), ]	
random3<-m[m$s1_patientID %in% c("D3") & m$s2_patientID %in% c("R2", "R1", "R4", "R5", "R6") & m$s1_human_group %in% c("donor") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes"), ]	
random4<-m[m$s1_patientID %in% c("D4") & m$s2_patientID %in% c("R2", "R3", "R1", "R5", "R6") & m$s1_human_group %in% c("donor") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes"), ]	
random5<-m[m$s1_patientID %in% c("D5") & m$s2_patientID %in% c("R2", "R3", "R4", "R1", "R6") & m$s1_human_group %in% c("donor") & m$s2_unhomogenized %in% c("yes") & m$s1_unhomogenized %in% c("yes"), ]	

# recipients to random donors:
allrandom<-rbind(random1, random2, random3, random4, random5)
allrandom$COMP<-"RANDOM"
allrandom$TIME<-allrandom$s2_timepoint

# donors to other random donors
donorrandom<-m[!m$s1_patientID==m$s2_patientID & m$s1_human_group %in% c("donor") & m$s2_human_group %in% c("donor"), ]
donorrandom$COMP<-"RDONORS"
donorrandom$TIME<--4

## now, combine all together:
comp<-rbind(alldr, allrandom, donorrandom)
# add a grouped name for each category:
comp$COMP2[comp$COMP %in% c("DR1", "DR2", "DR3", "DR4", "DR5")]<-"DR"
comp$COMP2[comp$COMP %in% c("RANDOM")]<-"RANDOM"
comp$COMP2[comp$COMP %in% c("DONORS")]<-"DONORS"
comp$COMP2[comp$COMP %in% c("RDONORS")]<-"RDONORS"


# I know that I want to offset the randoms by a little in my graph
#we will make 2 lists out of this:
	# set plot for variable 'n':
	# note: made two functions: one for plotting respective, and one for adding randoms
# first, separate into the respective and donor groups:
resp<-comp[comp$COMP2 %in% c("DR", "DONORS"), ]
randoms<-comp[comp$COMP2 %in% c("RANDOM", "RDONORS"), ]
# then, plot the first group:
df<-resp
distplot_overtime<- function(n) {
	# calculate some summary stats:
	xtime <- tapply(df$TIME,df$COMP2,function(x) return(x))
	yvar <- tapply(df[,n],df$COMP2,function(x) return(x))
	means <- tapply(df[,n],df$TIME,function(x) mean(x))
	sds <- tapply(df[,n],df$TIME,function(x) sd(x)/sqrt(length(x)))
	sdhigh<-means+sds
	sdlow<-means-sds
	cord.x<-c(unique(sort(df$TIME)), rev(unique(sort(df$TIME))))
	cord.y<-c(sdhigh, rev(sdlow))
	# donor labels:
	dlabels <- c("donors")
	dat <- unique(comp[comp$COMP2 %in% c("DONORS"), c("TIME")])
	# recipient labels:
	rlabels <- c("pre (~10d)", "pre(~2d)", "post (~2d)", "post (~10d)", "post (~30d)", "post (~200d)")
	rat <- unique(comp[comp$COMP2 %in% c("DR"), c("TIME")])
	xrec <- unique(sort(df$TIME))[2:7]
	# plot the data
	plot(df$TIME,df[,n],type="n", xlim=c(min(df$TIME)-1,max(df$TIME)+1), ylim=c(0,max(df[,n])), ylab=names(df[n]), xaxt='n', xlab="sampling")
	axis(1, at=dat, labels=dlabels, cex.axis=0.8, tck=-0.03, las=1)
	text(x =  rat, y = par("usr")[3]-max(df[,n])/50, srt = 45, adj = 1, labels = rlabels, xpd = TRUE, cex=0.8)
	abline(v=0, col = "grey80", lty=2)
		# if you want to add individual lines:
	colors<-colorRampPalette(c("grey20","grey80"))(n = 5)
	xreclines <- tapply(df$TIME, df$COMP, function(x) return(x))
	xreclines <- tapply(df[,n], df$COMP, function(x) return(x))
	mapply(lines,xreclines[2:6],yreclines[2:6],col=colors,pch=".",type="o")
		# add individual dots:
	points(jitter(df$TIME), df[,n], col=adjustcolor("grey45", alpha=0.5), pch=21, bg="grey80", cex=0.7)
		# add mean for recipients ONLY (and individual lines, if you want)
	lines(xrec, means[2:7], col="black", lwd=2)
	arrows(xrec, sdlow[2:7], xrec, sdhigh[2:7], length=0.05, angle=90, code=3, col="grey17", lwd=2)
		# add arrows for donor:
	arrows(dat, sdlow[1], dat, sdhigh[1], length=0.05, angle=90, code=3, col="grey17", lwd=2)
		# legend:
	}
# now, add points (offset) for random donors:
df2<-randoms
addgroup_distplot_overtime<- function(n) {
	# calculate some summary stats:
	xtime <- tapply(df2$TIME,df2$COMP2,function(x) return(x))
	yvar <- tapply(df2[,n],df2$COMP2,function(x) return(x))
	means <- tapply(df2[,n],df2$TIME,function(x) mean(x))
	sds <- tapply(df2[,n],df2$TIME,function(x) sd(x)/sqrt(length(x)))
	sdhigh<-means+sds
	sdlow<-means-sds
	cord.x<-c(unique(sort(df2$TIME)), rev(unique(sort(df2$TIME))))
	cord.y<-c(sdhigh, rev(sdlow))
	# donor labels:
	#dlabels <- c("donors")
	#dat <- unique(comp[comp$COMP2 %in% c("DONORS"), c("TIME")])
	# recipient labels--not usually necessary:
	#rlabels <- c("pre (~10d)", "pre(~2d)", "post (~2d)", "post (~10d)", "post (~30d)", "post (~200d)")
	#rat <- unique(comp[comp$COMP2 %in% c("DR"), c("TIME")])
	xrec <- unique(sort(df2$TIME))[2:7]
	# plot the data
	#plot(df2$TIME,df2[,n],type="n", xlim=c(min(df2$TIME)-1,max(df2$TIME)+1), ylim=c(0,max(df2[,n])), ylab=names(df2[n]), xaxt='n', xlab="sampling")
	#axis(1, at=dat, labels=dlabels, cex.axis=0.8, tck=-0.03, las=1)
	#text(x =  rat, y = par("usr")[3]-max(df2[,n])/50, srt = 45, adj = 1, labels = rlabels, xpd = TRUE, cex=0.8)
	#abline(v=0, col = "grey80", lty=2)
		# if you want to add individual lines--not usually necessary in added group
	#colors<-colorRampPalette(c("grey20","grey80"))(n = 5)
	#xreclines <- tapply(df2$TIME, df2$COMP, function(x) return(x))
	#yreclines <- tapply(df2[,n], df2$COMP, function(x) return(x))
	#mapply(lines,xrec[2:6],yrec[2:6],col=colors,pch=".",type="o")
		# add individual dots:
	points(jitter(df2$TIME)+0.33, df2[,n], col=adjustcolor("lightcyan2", alpha=0.5), pch=21, bg="lightcyan", cex=0.7)
		# add mean for recipients ONLY (and individual lines, if you want)
	lines(xrec+0.33, means[2:7], col="darkcyan", lwd=2)
	arrows(xrec+0.33, sdlow[2:7], xrec+0.33, sdhigh[2:7], length=0.05, angle=90, code=3, col="darkcyan", lwd=2)
		# add arrows for donor:
	arrows(dat+0.33, sdlow[1], dat+0.33, sdhigh[1], length=0.05, angle=90, code=3, col="darkcyan", lwd=2)
	}
# now, apply to distance of choice:
dev.new(width=4.669725, height=4.073394)
lapply(c("sharedsobs"), FUN=distplot_overtime)
legend("bottomleft", c("donor-recipient pairs", "random donor"), col=c("grey45", "lightcyan2"), pch=21, cex=0.7, pt.bg=c("grey80", "lightcyan"))
lapply(c("sharedsobs"), FUN=addgroup_distplot_overtime)
		# note: for some reason, mapply for individual lines does not work...
		# also, can't topple these functions together


```


##### Fig. 2C and Table S2: unique OTUs in recipient-donor samples over time


```
# let's see how many OTUs shared overall overlap between the pairs:
# lists: shared1, shared2, shared3, shared4, shared5
# (repeat from above)
d<-read.table(file="umfmt.allmeasures.txt", header=TRUE)
rownames(d)<-d$sampleID

# get only OTUs and sampleID, patientID:
df<-d[, c(1, 3, 6, 8, 70:length(d))]

# convert to 1/0 (presence/absence)
df[df > 0]<-1	

### Let's calculate how many OTUs are:
	- shared between recipients/donors (shared by all conditions-indistinguishable)
	- unique to recipient pre-FMT (only detected in recipient, pre-FMT)
	- shared OTUs unique to respective donor and found in recipient post-FMT
	- novel: never detected in recipient or donor prior to FMT

## make lists of OTU dataframes by patient:
# in donors:
listFunction <- function(n) {
	# get total OTU list for donor samples only
	otus <- df[df$patientID==(n) & df$human_group=="donor",  5:length(df)]
	# remove samples that we do not have sequences for
	otus <- otus[!is.na(otus$Otu000001), ]
	# get ONLY OTUs that are present in any of the samples associated with that patient
	otus_list<-names(which(colSums(otus == 1) > 0))
	}
	# note: still need to redefine the human_group everytime
# for these donors, apply this function to get a list of the possible OTUs ever found in each patient:
donors.list<-lapply(c("D1", "D2", "D3", "D4", "D5"), FUN=listFunction)
# how many OTUs were ever found in each donor?
lapply(donors.list, FUN=length)
# how many of these OTUs are ONLY EVER found in that donor?
unique.donors<-lapply(1:length(donors.list), function(n) setdiff(donors.list[[n]], unlist(donors.list[-n])))

# in recipients, post (OTUs ever found in patients post-FMT, and unique to each patient across post-samples):
listFunction <- function(n) {
	otus <- df[df$patientID==(n) & df$human_group=="post",  5:length(df)]
	otus <- otus[!is.na(otus$Otu000001), ]
	otus_list<-names(which(colSums(otus == 1) > 0))
	}
post.list<-lapply(c("R1", "R2", "R3", "R4", "R5"), FUN=listFunction)
lapply(post.list, FUN=length)
unique.post<-lapply(1:length(post.list), function(n) setdiff(post.list[[n]], unlist(post.list[-n])))

# pre-fmt recipients (OTUs ever found in patients pre-FMT, and unique to each patient across pre-FMT samples):
listFunction <- function(n) {
	otus <- df[df$patientID==(n) & df$human_group=="pre",  5:length(df)]
	otus <- otus[!is.na(otus$Otu000001), ]
	otus_list<-names(which(colSums(otus == 1) > 0))
	}
	# note: still need to redefine the human_group everytime
pre.list<-lapply(c("R1", "R2", "R3", "R4", "R5"), FUN=listFunction)
unique.pre<-lapply(1:length(pre.list), function(n) setdiff(pre.list[[n]], unlist(pre.list[-n])))

## Also need to calculate which OTUs are unique across the samples:
# the previous lists specify the potential OTUs found in each patient, as well as the ones unique across each category
# need to also see which OTUs are found across the conditions (donor/pre/post)

# which OTUs are shared between donors and recipients post-FMT?
shared.dr<-lapply(1:length(post.list), function(n) intersect(donors.list[[n]], post.list[[n]]))
lapply(shared.dr, FUN=length)
# unique OTUs that appear in dr pairs, but not anywhere else in dataset (so ONLY ones found in across the patient data sets)
unique.dr.otus<-lapply(1:length(shared.dr), function(n) setdiff(shared.dr[[n]], unlist(shared.dr[-n])))
lapply(unique.dr.otus, FUN=length)

## let's calculate some percentages:
	# what percent of OTUS in post samples are shared?
lapply(1:length(donors.list), function(n) length(shared.dr[[n]])/length(post.list[[n]])*100)
	# what percentage shared OTUs are unique to dr pairs?
lapply(1:length(donors.list), function(n) length(unique.dr.otus[[n]])/length(shared.dr[[n]])*100)

## need to also account for presence of pre-FMT OTUs...
# let's make these, and see how many are present before/after, and how many compare to donors

# ok, what OTUs overlap with post/pre?
pre.post.comp<-lapply(1:length(post.list), function(n) intersect(pre.list[[n]], post.list[[n]]))
lapply(pre.post.comp, FUN=length)
	# before calculating percent that is still detected in post, let's see how many overlap with our shared...
overlap.shared<-lapply(1:length(shared.dr), function(n) intersect(shared.dr[[n]], pre.list[[n]]))
	# what percent of sequences in post samples are detected prior to FMT?
lapply(1:length(donors.list), function(n) length(overlap.shared[[n]])/length(shared.dr[[n]])*100)

# let's see if any of our unique OTUs overlap with the previously detected OTUs...
overlap.shared.pre<-lapply(1:length(overlap.shared), function(n) intersect(overlap.shared[[n]], unique.dr.otus[[n]]))
# very few!!! Let's make a final 'rare' list:
unique.dr.filtered<-lapply(1:length(overlap.shared), function(n) setdiff(unique.dr.otus[[n]], pre.list[[n]]))
lapply(unique.dr.filtered, FUN=length)

#### redo of above (systematically):

# 1) get list of post OTUs that were not detected before, calculate percentages of pre in post
post.notinpre<-lapply(1:length(donors.list), function(n) setdiff(post.list[[n]], pre.list[[n]]))
lapply(post.notinpre, FUN=length)
	# percent of OTUs in post not detected in pre:
lapply(1:length(donors.list), function(n) length(post.notinpre[[n]])/length(post.list[[n]])*100)
	# percent of pre in post:
pre.in.post<-lapply(1:length(donors.list), function(n) intersect(pre.list[[n]], post.list[[n]]))
percent.pre<-lapply(1:length(donors.list), function(n) length(pre.in.post[[n]])/length(post.list[[n]])*100)
# this represents the % of the sequences that are detected beforehand (note: could include shared--will need to make a group that is detected in all three)

# more specifics:
# 1) get list of shared OTUs between dr pairs, calculate % shared by dr pairs that are NOT detected in pre
	# OTUs shared between donors and post:
shared.dr<-lapply(1:length(donors.list), function(n) intersect(post.list[[n]], donors.list[[n]]))
lapply(shared.dr, FUN=length)
# percent total shared (includes the pre)
percent.allshared<-lapply(1:length(donors.list), function(n) length(shared.dr[[n]])/length(post.list[[n]])*100)
	# remove OTUs that appear in pre:
shared.dr.nopre<-lapply(1:length(donors.list), function(n) setdiff(shared.dr[[n]], pre.list[[n]]))
lapply(shared.dr.nopre, FUN=length)
percent.shared.notinpre<-lapply(1:length(donors.list), function(n) length(shared.dr.nopre[[n]])/length(post.list[[n]])*100)
# this % is slightly lowered after accounting for the pre
mean(unlist(percent.shared.notinpre))
#[1] 64.49372

# 2) get list of OTUs detected in pre samples:
	# OTUs shared between donors and post:
shared.pp<-lapply(1:length(donors.list), function(n) intersect(post.list[[n]], pre.list[[n]]))
lapply(shared.pp, FUN=length)
	# remove OTUs that appear in donor:
shared.pp.notind<-lapply(1:length(donors.list), function(n) setdiff(shared.pp[[n]], donors.list[[n]]))
lapply(shared.pp.notind, FUN=length)
percent.pp.notind<-lapply(1:length(donors.list), function(n) length(shared.pp.notind[[n]])/length(post.list[[n]])*100)
mean(unlist(percent.pp.notind))
#[1] 6.11206

# 3) get list of OTUs detected in BOTH (within post, so detected in all 3):
detected.inall<-lapply(1:length(donors.list), function(n) intersect(shared.dr[[n]], shared.pp[[n]]))
lapply(detected.inall, FUN=length)
percent.det.inall<-lapply(1:length(donors.list), function(n) length(detected.inall[[n]])/length(post.list[[n]])*100)
mean(unlist(percent.det.inall))
#[1] 19.19625

# 4) what OTUs are not detected in either?
# the remaining sequences are 'new'--so coming from environment
percent.newly.detected<-lapply(1:length(donors.list), function(n) 100-(percent.det.inall[[n]]+percent.pp.notind[[n]]+percent.shared.notinpre[[n]]))
#[1] 10.19797

# 5) finally, what OTUs are 'rare'--only identified in individual donor-post samples, not detected in ANY pres across the board?
	# identify what shared OTUs are rare intra-individually
unique.dr.otus<-lapply(1:length(shared.dr), function(n) setdiff(shared.dr[[n]], unlist(shared.dr[-n])))
lapply(unique.dr.otus, FUN=length)
	# remove all pre OTUs from this list:
unique.dr.nopre<-lapply(1:length(shared.dr), function(n) setdiff(unique.dr.otus[[n]], unlist(pre.list)))
lapply(unique.dr.nopre, FUN=length)
percent.unique<-lapply(1:length(donors.list), function(n) length(unique.dr.nopre[[n]])/length(post.list[[n]])*100)
mean(unlist(percent.unique))
#[1] 5.896061

####
# now, let's make a bargraph of the percent matching overall in each:
dev.new(width=4.669725, height=4.073394)
shared<-do.call("rbind", percent.det.inall)		# if you just wanted one list

alldf<-do.call(rbind, Map(data.frame, donor=percent.shared.notinpre, shared=percent.det.inall, pre=percent.pp.notind, novel=percent.newly.detected))
rownames(alldf)<-c("R1", "R2", "R3", "R4", "R5")

# plot it --> Fig. 2B:
barplot(as.matrix(t(alldf)), col=c("dodgerblue2", "turquoise3", "goldenrod2", "plum1"), ylab="% OTUs detected", xlim=c(0,8))
legend("topright", col=c("dodgerblue2", "turquoise3", "goldenrod2", "plum1"), colnames(alldf), pch=19, cex=0.8)
# add text of percentages (unique dr, inter-individually)
text(x=c(0.7, 1.9, 3.1, 4.3, 5.5), y=rep(3, each=5), as.character(round(unlist(percent.unique), digits=1)))



#### finally, we can get a list of which types of OTUs are unique to donor-recipient pairs, and never found in other patients OR in the recipient's pre-sample:
# this list has the OTUs for the shared pairs
unique.dr.nopre

# let's merge these using our taxonomic classification:
taxa<-read.table(file="allhumos4.taxonomy.names.txt", header=TRUE)
# add family names:
tax <- taxa$Taxonomy
tax <- gsub("\\(\\d*\\)", "", tax)
tax <- gsub(";unclassified", "", tax)
tax <- gsub("_1", "", tax)
tax <- gsub(";$", "", tax)
tax <- gsub("/.*", "", tax)
tax <- gsub(".*les;", "", tax)
tax <- gsub(";.*", "", tax)
taxa$family<-tax
ftaxa<-taxa[, c("OTU", "taxname", "phylum", "family")]

# if you wanted to merge this individually per patient (not done):
#u1<-unlist(unique.dr.nopre[1])
#test<-data.frame("Otu" = u1, "Pair" = c("dr1"))

# you can merge the names within EACH of these to get a full list of OTUs
	# get a dataframe of these OTUs as a list:
unique.dr.nopre.classified<-lapply(1:length(unique.dr.nopre), function(n,d) merge(as.data.frame(unique.dr.nopre[[n]]), ftaxa, by.x=names(as.data.frame(unique.dr.nopre[[n]])), by.y="OTU"))
names(unique.dr.nopre.classified)<-c("dr1", "dr2", "dr3", "dr4", "dr5")
# collapse them together, and add the pair identity:
all.unique.otus <- do.call("rbind", unique.dr.nopre.classified)
all.unique.otus$pair <- rep(names(unique.dr.nopre.classified), sapply(unique.dr.nopre.classified, nrow))
names(all.unique.otus)[1]<-"OTU"
all.unique.otus[order(all.unique.otus$pair, all.unique.otus$phylum, all.unique.otus$family), ]
#write.table(all.unique.otus, file="figures/paper_figsTableS2_unique.dr.otus.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# this is Table S2


```