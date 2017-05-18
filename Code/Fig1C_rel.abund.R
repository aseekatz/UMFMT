#### Figure 1C and S1: Relative abundance of pre, post, and donor samples
##### Anna M. Seekatz
##### 4.12.17

##### Files needed:
	- directory: /Users/annaseekatz/Box Sync/Projects/UM.FMT/UMFMT_paper/16S_UMFMT
	- Fig. 1C: relative abundance plots by group
		- (files created in prep code):
		- relative abundance of top 98% genera: umfmt_genfrac2p.all_w.meta.txt
		- color assignment: umfmt_genfrac2p.txt
	- Fig. S1: OTU heatmap
		- umfmt_metadata.txt
		- mothurfiles/allhumos4.final.0.03.pick.0.03.pick.0.03.filter.shared
		- mothurfiles/allhumos4.final.0.03.cons.taxonomy

#####  Fig. 1C: comparison of top genera in pre, post, and donor samples

```

# read in phylotype information:
data<-read.table(file="umfmt_genfrac2p.all_w.meta.txt", header=TRUE)
data2<-read.table(file="umfmt_genfrac2p.txt", header=TRUE)
colors<-data2$color

# figure with top 20 (?) genera:
# if you want to organize genera by top genera present:
#df<-data[, 10:56]					# want other at end
#df<-df[ , order(colSums(df)/42, decreasing=TRUE)]
#df<-cbind(data[, c(1:9)], df, other=data[, c("other")])
df<-cbind(data[, c(1:9)], data[, c(10:54, )], unclassified_Bacteria=data[, c("Bacteria_unclassified")], other=data[, c("other")])

# group within an array, then calculate means for numeric columns
df<-data
df$human_group <- ordered(df$human_group, levels = c("pre", "post", "donor"))

# to get a matrix of a summary statistic by a certain variable across several columns:
apply(df[,10:ncol(df)], 2, function(x) tapply(x, df$human_group, mean))

# can also do this:
# get medians/min/max:
median<- sapply(levels(df$human_group), function(x) sapply(df[df$human_group == x, c(10:57)], median))
min<- sapply(levels(df$human_group), function(x) sapply(df[df$human_group == x, c(10:57)], function(x) quantile(x)[2]))
max<- sapply(levels(df$human_group), function(x) sapply(df[df$human_group == x, c(10:57)], function(x) quantile(x)[4]))
means<- sapply(levels(df$human_group), function(x) sapply(df[df$human_group == x, c(10:57)], mean))
sds<- sapply(levels(df$human_group), function(x) sapply(df[df$human_group == x, c(10:57)], sd))
ses<- sapply(levels(df$human_group), function(x) sapply(df[df$human_group == x, c(10:57)], function(x) sd(x)/sqrt(length(x))))


# then plot:
bar.col<-c(as.character(colors), "grey47")
par(mfrow=c(3,1))
par(mar=c(4,4,1,1))
par(oma=c(0,0,0,0))
par(xpd=TRUE)
	# for plotting median:
mp1<-barplot(median[,c("pre")], beside=TRUE,col=bar.col,las=2, ylim=c(0,35), ylab="Relative Abundance")
segments(mp1, max[,c("pre")], mp1, min[,c("pre")])
mp2<-barplot(median[,c("post")], beside=TRUE,col=bar.col,las=2, ylim=c(0,35), ylab="Relative Abundance")
segments(mp2, max[,c("post")], mp2, min[,c("post")])
mp3<-barplot(median[,c("donor")], beside=TRUE,col=bar.col,las=2, ylim=c(0,35), ylab="Relative Abundance")
segments(mp3, max[,c("donor")], mp3, min[,c("donor")])
#legend(36, 70, c(), col=bar.col, pch=15, cex=0.8)

	# Fig. 1C: for plotting mean + sd:
mp1<-barplot(means[,c("pre")], beside=TRUE,col=bar.col,las=2, ylim=c(0,35), ylab="Relative Abundance")
segments(mp1, means[,c("pre")]-sds[,c("pre")], mp1, means[,c("pre")]+sds[,c("pre")])
mp2<-barplot(means[,c("post")], beside=TRUE,col=bar.col,las=2, ylim=c(0,35), ylab="Relative Abundance")
segments(mp2, means[,c("post")]-sds[,c("post")], mp2, means[,c("post")]+sds[,c("post")])
mp3<-barplot(means[,c("donor")], beside=TRUE,col=bar.col,las=2, ylim=c(0,35), ylab="Relative Abundance")
segments(mp3, means[,c("donor")]-sds[,c("donor")], mp3, means[,c("donor")]+sds[,c("donor")])

	# Fig. 1C: for plotting mean + se:
bar.col<-c(as.character(colors), "grey47")
par(mfrow=c(3,1))
par(mar=c(0,4,0.5,1))
par(oma=c(4,0,2,0))
par(xpd=TRUE)
mp1<-barplot(means[,c("pre")], beside=TRUE,col=bar.col,las=2, ylim=c(0,30), ylab="Relative Abundance", xaxt="n")
segments(mp1, means[,c("pre")]-ses[,c("pre")], mp1, means[,c("pre")]+ses[,c("pre")])
mp2<-barplot(means[,c("post")], beside=TRUE,col=bar.col,las=2, ylim=c(0,30), ylab="Relative Abundance", xaxt="n")
segments(mp2, means[,c("post")]-ses[,c("post")], mp2, means[,c("post")]+ses[,c("post")])
mp3<-barplot(means[,c("donor")], beside=TRUE,col=bar.col,las=2, ylim=c(0,30), ylab="Relative Abundance")
segments(mp3, means[,c("donor")]-ses[,c("donor")], mp3, means[,c("donor")]+ses[,c("donor")])
dev.size()
#[1] 4.825688 5.825688

### stats:
data<-df
df<-data[!data$unhomogenized %in% c("no"), c(6, 10:57)]

modelList<-list()
for(i in 2:49){
    fmla <- formula(paste(names(df)[i], " ~human_group"))
    modelList[[i]]<-kruskal.test(fmla, data = df)
}
modelList

# significant:
data:  Fusobacterium by human_group
Kruskal-Wallis chi-squared = 6.4286, df = 2, p-value = 0.04018

data:  Collinsella by human_group
Kruskal-Wallis chi-squared = 8.4363, df = 2, p-value = 0.01473

data:  Bifidobacterium by human_group
Kruskal-Wallis chi-squared = 11.982, df = 2, p-value = 0.002501

data:  Campylobacter by human_group
Kruskal-Wallis chi-squared = 8.8524, df = 2, p-value = 0.01196

data:  Enterobacteriaceae_unclassified by human_group
Kruskal-Wallis chi-squared = 8.6289, df = 2, p-value = 0.01337

data:  Escherichia.Shigella by human_group
Kruskal-Wallis chi-squared = 6.5105, df = 2, p-value = 0.03857

data:  Anaeroglobus by human_group
Kruskal-Wallis chi-squared = 9.2128, df = 2, p-value = 0.009988

data:  Megasphaera by human_group
Kruskal-Wallis chi-squared = 12.434, df = 2, p-value = 0.001996

data:  Pediococcus by human_group
Kruskal-Wallis chi-squared = 10.382, df = 2, p-value = 0.005567

data:  Veillonellaceae_unclassified by human_group
Kruskal-Wallis chi-squared = 13.304, df = 2, p-value = 0.001291

data:  Coprococcus by human_group
Kruskal-Wallis chi-squared = 21.636, df = 2, p-value = 2.004e-05

data:  Dorea by human_group
Kruskal-Wallis chi-squared = 12.822, df = 2, p-value = 0.001644

data:  Veillonella by human_group
Kruskal-Wallis chi-squared = 8.3022, df = 2, p-value = 0.01575

data:  Clostridium_IV by human_group
Kruskal-Wallis chi-squared = 15.3, df = 2, p-value = 0.0004761

data:  Ruminococcus2 by human_group
Kruskal-Wallis chi-squared = 18.13, df = 2, p-value = 0.0001156

data:  Anaerostipes by human_group
Kruskal-Wallis chi-squared = 21.158, df = 2, p-value = 2.544e-05

data:  Oscillibacter by human_group
Kruskal-Wallis chi-squared = 19.725, df = 2, p-value = 5.21e-05

data:  Phascolarctobacterium by human_group

data:  Ruminococcus by human_group
Kruskal-Wallis chi-squared = 14.207, df = 2, p-value = 0.0008221

data:  Faecalibacterium by human_group
Kruskal-Wallis chi-squared = 14.982, df = 2, p-value = 0.000558

data:  Clostridium_XlVa by human_group
Kruskal-Wallis chi-squared = 9.5027, df = 2, p-value = 0.00864

data:  Clostridium_XI by human_group
Kruskal-Wallis chi-squared = 17.769, df = 2, p-value = 0.0001385

data:  Ruminococcaceae_unclassified by human_group
Kruskal-Wallis chi-squared = 22.003, df = 2, p-value = 1.668e-05

data:  Blautia by human_group
Kruskal-Wallis chi-squared = 14.541, df = 2, p-value = 0.0006958

data:  Lactobacillus by human_group
Kruskal-Wallis chi-squared = 15.901, df = 2, p-value = 0.0003525

data:  Lachnospiraceae_unclassified by human_group

data:  Clostridiales_unclassified by human_group
Kruskal-Wallis chi-squared = 17.109, df = 2, p-value = 0.0001926

data:  Firmicutes_unclassified by human_group
Kruskal-Wallis chi-squared = 15.201, df = 2, p-value = 0.0005003

data:  Bacteroides by human_group
Kruskal-Wallis chi-squared = 21.337, df = 2, p-value = 2.327e-05

data:  Alistipes by human_group
Kruskal-Wallis chi-squared = 9.5656, df = 2, p-value = 0.008373

data:  Parabacteroides by human_group
Kruskal-Wallis chi-squared = 15.985, df = 2, p-value = 0.0003379


# wilcoxon:
data<-df
df<-data[data$unhomogenized %in% c("yes") & data$human_group %in% c("pre", "post"), c(6, 10:57)]

modelList<-list()
for(i in 2:49){
    fmla <- formula(paste(names(df)[i], " ~human_group"))
    modelList[[i]]<-wilcox.test(fmla, data = df, paired=FALSE)
}
modelList

# p < 0.05:
Fusobacterium
Campylobacter
Enterobacteriaceae
e.shigella
Anaeroglobus
Veillonella

# p < 0.005:
Collinsella
	Bifidobacterium
	Megasphaera
	Pediococcus
	Veillonellaceae
	Dorea
	Ruminococcus
	Clostridium_XlVa
	Lachnospiraceae_unclassified
	Lactobacillus
	Alistipes

# p < 0.0005
	Coprococcus
	Roseburia
	Clostridium_IV
	Ruminococcus2
	Anaerostipes
	Faecalibacterium
	Oscillibacter
	Blautia
	Ruminococcaceae_unclassified
	Firmicutes_unclassified
	Clostridium_XI
	Clostridiales_unclassified
	Parabacteroides
	Bacteroides

```


#####  Fig. S1: heatmap of top OTUs

```
library(gplots)
library(vegan)
library(plyr)

# read in file:
meta<-read.table(file="umfmt_metadata.txt", header=TRUE)
shared<-read.table(file="mothurfiles/allhumos4.final.0.03.pick.0.03.pick.0.03.filter.shared", header=TRUE, row.names=2)
dim(shared)
shared$seqID<-rownames(shared)
	# merge with meta:
sum.shared<-merge(meta, shared, by.x="seqID", by.y="seqID")
sum.shared<-subset(sum.shared, select =-c(label, numOtus) )
sum.shared<-droplevels(sum.shared)
#write.table(sum.shared, file="umfmt_otus.w.meta.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# create a taxonomy file:
taxonomy_file<-read.table(file="mothurfiles/allhumos4.final.0.03.cons.taxonomy", header=TRUE)
tax <- taxonomy_file$Taxonomy
tax <- gsub("\\(\\d*\\)", "", tax)
tax <- gsub(";unclassified", "", tax)
tax <- gsub("_1", "", tax)
tax <- gsub(";$", "", tax)
tax <- gsub("/.*", "", tax)
tax <- gsub(".*;", "", tax)
tax.names <-paste(taxonomy_file$OTU, tax, sep="_")
tax.names <-gsub("000", "", tax.names)
taxonomy_file$taxname<-tax.names
phylum <- taxonomy_file$Taxonomy
phylum <- gsub("\\(\\d*\\)", "", phylum)
phylum <- gsub("Bacteria;", "", phylum)
phylum <- gsub(";$", "", phylum)
phylum <- gsub(";.*", "", phylum)
taxonomy_file$phylum<-phylum
tax<-taxonomy_file
tax<-tax[order(tax$OTU),]
#write.table(tax, file="allhumos4.taxonomy.names.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#tax<-read.table(file="allhumos4.taxonomy.names.txt", header=TRUE)

# get relative abundance and split meta (needs to remain in same order, though)
otu.matrix<-sum.shared[, 10:ncol(sum.shared)]
otu.rel<-otu.matrix/rowSums(otu.matrix)
rownames(otu.rel)<-sum.shared$sampleID
meta<-sum.shared[, 1:9]

	# add some colors based on the metadata:
meta$group.col<-mapvalues(meta$human_group, from = c("donor", "post", "pre"), to = c("chartreuse4", "blue", "gold"))
	# check that the names match:
cbind(row.names(otu.rel), as.character(meta$sampleID))

# clustering by all OTUs:
otu.rel.horn<-vegdist(otu.rel, method="horn")
row.cluster.horn<-hclust(otu.rel.horn, "average")
otu.rel.bray<-vegdist(otu.rel, method="bray")
row.cluster.bray<-hclust(otu.rel.bray, "average")

# let's get the top 50 OTUs, but cluster by the whole data set
topotus<- otu.rel[, order(-colSums(otu.rel))]	
top50<-topotus[, 1:50]
top40<-topotus[, 1:40]
otus<-top50

# filter this list to only include those in the heatmap:
#tax<-read.table(file="allhumos3.taxonomy.names.txt", header=TRUE)
topotu<-colnames(otus)									#define the OTUs you want to keep (those in the top 98%)
filtered.tax<-tax[tax$OTU %in% topotu , ]							#filter out the taxonomy metadata file
filtered.tax<-droplevels(filtered.tax)											#sometimes R keeps levels that were discarded--this gets rid of them completely

tax.col<-filtered.tax[,c("OTU", "taxname", "Size", "phylum")] #take the column with the group
tax.col$color<-mapvalues(tax.col$phylum, from = c("Actinobacteria","Bacteroidetes", "Firmicutes", "Fusobacteria", "Proteobacteria", "Synergistetes", "Verrucomicrobia"), to = c("darkmagenta", "darkgreen","skyblue", "darkred", "yellow", "cyan", "hotpink"))
tax.col<- tax.col[order(tax.col$phylum, -tax.col$Size) , ]
tax.col<-t(tax.col)
phyla.col<-tax.col[5,]

# order heatmap according to topotus, then bind to column variable:
col.order<-as.character(tax.col[1, ])								#convert the order of the OTUs into a list
otus<-otus[,col.order]											#order your matrix by the ordered list
rbind(colnames(otus), tax.col)									#check order of the matrix and list to ensure that colors will be correct (they should match!)

# make heatmaps:
my.col <- colorRampPalette(c("aliceblue","lightskyblue", "deepskyblue", "dodgerblue2", "blue4"))(n = 249)			#you can put in as many colors as you want here
my.breaks = c(seq(0,0.001,length=50),  								#just remember: your breaks must always = 1 more than the number of colors you are defining!
               seq(0.001,0.01,length=51),
               seq(0.01,0.1,length=51),
               seq(0.10,0.5,length=51),
               seq(0.50,1,length=51))
# does R 3.3.0 make these differently?!
edited<-unique(my.breaks)
	# with clustering by morisita-horn:
heatmap.2(as.matrix(otus),
                notecol="black",     
                density.info="none",
                key.xlab="",
                key.title="",  
                trace="none",         
                margins =c(4,4),    
                col=my.col,        
                breaks=edited,    
                dendrogram="none",    
                Colv=F,
                symkey=FALSE,
                Rowv=as.dendrogram(row.cluster.horn),
                srtCol=45,
                cexRow= 0.7,
                cexCol = 0.8,
                key=FALSE,
                lwid = c(3,5),
                lhei = c(2,5),
                ylab='sample',
                RowSideColors = as.character(meta$group.col),
  		        ColSideColors = as.character(phyla.col),
  		        labCol=tax.col[2,]
                    )
#legend("bottomleft",legend=c("healthy_mouse", "mFMT", "noFMT", "hFMT", "hFMT-pre", "human: donor", "human: pre", "human: post"), col=c("darkolivegreen", "chartreuse4", "magenta", "darkcyan", "yellow", "deepskyblue", "gold", "deepskyblue4"), cex=0.6, pch=19)
legend("topleft",legend=c("donor", "post", "pre"), pt.bg=c("chartreuse4", "blue", "gold"), cex=1, pch=22, col="black")
legend("top",legend=c("0-0.001", "0.001-0.01", "0.01-0.1", "0.1-0.5", "0.5-1"), col=c("aliceblue","lightskyblue", "deepskyblue", "dodgerblue2", "blue4"), cex=1, pch=15)
legend("topright",legend=c("Actinobacteria","Bacteroidetes", "Firmicutes", "Fusobacteria", "Proteobacteria", "Synergistetes", "Verrucomicrobia"), pt.bg=c("darkmagenta", "darkgreen","skyblue", "darkred", "yellow", "cyan", "hotpink"), cex=1, pch=22, col="black")


```