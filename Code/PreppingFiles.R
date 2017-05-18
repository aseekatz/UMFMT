#### Prepping mothur output files
##### Anna M. Seekatz
##### 2.3.17
######## directory: /Users/annaseekatz/Box Sync/Projects/UM.FMT/UMFMT_paper/16S_UMFMT

##### Mothur output files and resulting modified files:
	- alpha diversity measures --> umfmt_summary.txt
		- input file is a summary file of alpha diversity measures created in mothur
		- output file combines these with metadata
	- allhumos4.final.0.03.cons.taxonomy --> allhumos4.taxonomy.names.txt
		- input file lists the taxonomic classifications for each of the OTUs at 0.03 cutoff
		- output file was modified to reflect new OTU names that includes the taxonomy, as well as cleaner classifications
	- allhumos4.final.0.03.pick.0.03.filter.0.03.pick.summary --> umfmt.betasummary.txt
		- input file represents pairwise distances (beta diversity) generated in mothur, filtered to include only relevant samples
		- output file adds sample meta data to sampleIDs
	- allhumos4.final.0.03.pick.0.03.pick.0.03.filter.shared, umfmt_metadata.txt --> umfmt_otus.w.meta.txt
		- input .shared file was previously filtered using the specified measures
		- output file combines metadata with OTU counts
	- combining all data (metadata, alpha metrics, relative abundance of OTUs, and metabolites) --> umfmt_allmeasures.txt
		- all data input files
		- output file has all variables, combined in one file
	- --> umfmt_genfrac2p.all_w.meta.txt
		- input file is a file produced in mothur classifying sequences directly to the RDP database
		- output file is a file with phylotype information (relative abundance of genus-level sequence assignments) 
	
##### Other files relevant to project:
	- compiling GEE results from Krishna's code --> fmtrecip_allresults.txt
		- this represents Table S3

##### Data files compiled from experiments or from cores:
	- umfmt_metadata.txt
	- umfmt_SCFA.txt
	- umfmt_suclac.txt
	- 	

###### alpha summaries --> umfmt_summary.txt
	- allhumos4.final.0.03.pick.0.03.pick.groups.summary
	- allhumos4.final.0.03.pick.0.03.pick.thetayc.0.03.lt.pcoa.axes
	- allhumos4.final.0.03.pick.0.03.pick.thetayc.0.03.lt.nmds.axes
	- umfmt_metadata.txt
	
```
# read in files and merge together:
meta<-read.table(file="umfmt_metadata.txt", header=TRUE)
# add a '_hom' to those samples that were unhomogenized:
meta$group2<-meta$human_group
meta$group2<-as.character(meta$group2)
meta$group2[meta$unhomogenized==c("no")]<-"donor_hom"
meta$group2<-as.factor(meta$group2)

pcoa<-read.table(file="mothurfiles/allhumos4.final.0.03.pick.0.03.pick.thetayc.0.03.lt.pcoa.axes", header=TRUE)
	pcoa<-pcoa[,1:4]
	colnames(pcoa)[2:4] <- paste("pcoa03", colnames(pcoa)[2:4], sep = "_")
	colnames(pcoa)[1]<-"sampleID"
nmds<-read.table(file="mothurfiles/allhumos4.final.0.03.pick.0.03.pick.thetayc.0.03.lt.nmds.axes", header=TRUE)
	nmds<-nmds[1:4]
	colnames(nmds)[2:4] <- paste("nmds03", colnames(nmds)[2:4], sep = "_")
	colnames(nmds)[1]<-"sampleID"
sum<-read.table(file="mothurfiles/allhumos4.final.0.03.pick.0.03.pick.groups.summary", header=TRUE)
	sum<-subset(sum, select=-c(label))
	colnames(sum)[2:16] <- paste(colnames(sum)[2:16], "03", sep = "_")
	colnames(sum)[1]<-"sampleID"

combined.pcoa<-merge(meta, pcoa, by.x=c("seqID"), by.y=c("sampleID"))
combined.nmds<-merge(combined.pcoa, nmds, by.x=c("seqID"), by.y=c("sampleID"))
combined.sum<-merge(combined.nmds, sum, by.x=c("seqID"), by.y=c("sampleID"))
#write.table(combined.sum, 'umfmt_summary.txt',quote=FALSE,sep="\t", col.names=TRUE, row.names=FALSE)

```



###### allhumos4.final.0.03.cons.taxonomy --> allhumos4.taxonomy.names.txt

```
taxonomy_file<-read.table(file="mothurfiles/allhumos4.final.0.03.cons.taxonomy", header=TRUE)
# taxname to OTU:
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
# phylum variable:
phylum <- taxonomy_file$Taxonomy
phylum <- gsub("\\(\\d*\\)", "", phylum)
phylum <- gsub("Bacteria;", "", phylum)
phylum <- gsub(";$", "", phylum)
phylum <- gsub(";.*", "", phylum)
taxonomy_file$phylum<-phylum
# family variable:
fam <- taxonomy_file$Taxonomy
fam <- gsub("\\(\\d*\\)", "", tax)
fam <- gsub(";unclassified", "", tax)
fam <- gsub("_1", "", tax)
fam <- gsub(";$", "", tax)
fam <- gsub("/.*", "", tax)
fam <- gsub(".*les;", "", tax)
fam <- gsub(";.*", "", tax)
taxonomy_file$family<-fam
#write.table(tax, file="allhumos4.taxonomy.names.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# note: some family level marks may still be messed up due to naming scheme
# these can be addressed when necessary

```


###### Subsetting pairwise distances specific to human files
	-  allhumos4.final.0.03.pick.0.03.filter.0.03.pick.summary,  umfmt_metadata.txt --> umfmt.betasummary.txt
	
```
# read in file:
mdist<-read.table(file="/Users/aseekatz/Box Sync/Projects/UM.FMT/UMFMT_paper/16S_UMFMT/mothurfiles/previously_filtered/allhumos4.final.0.03.pick.0.03.filter.0.03.pick.summary", header=TRUE)
	# note: the headers in this file had to be edited a bit (weird thing in mothur)
#write.table(mdist, file="umfmt.betasummary.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# now add metadata for each sample comparison:
var<-read.table(file='umfmt_metadata.txt', header=TRUE)
	# merge files:
	#merge by sample1
m1<-merge(var, mdist, by.x=c("seqID"), by.y=c("s1")) #this merges the data based on the sampleID/group1 match
colnames(m1)[1:9] <- paste("s1", colnames(m1)[1:9], sep = "_")	
	#add sample information for 2nd sample: ('group2'):
m2<-merge(var, m1, by.x=c("seqID"), by.y=c("s2"))
colnames(m2)[1:9] <- paste("s2", colnames(m2)[1:9], sep = "_")
#write.table(m2, file="umfmt.betasummary.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

```


###### Creating an OTU count file:
	- allhumos4.final.0.03.pick.0.03.pick.0.03.filter.shared, umfmt_metadata.txt --> umfmt_otus.w.meta.txt
	
```
## heatmap with metadata and taxonomy, and clustering:
	# read in files:
meta<-read.table(file="umfmt_metadata.txt", header=TRUE)
shared<-read.table(file="mothurfiles/allhumos4.final.0.03.pick.0.03.pick.0.03.filter.shared", header=TRUE, row.names=2)
dim(shared)
shared$seqID<-rownames(shared)
	# merge with meta:
sum.shared<-merge(meta, shared, by.x="seqID", by.y="seqID")
sum.shared<-subset(sum.shared, select =-c(label, numOtus) )
sum.shared<-droplevels(sum.shared)
#write.table(sum.shared, file="umfmt_otus.w.meta.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

```

###### Combining bile acid data with metadata:
	- umfmt_summary.txt
	- umfmt_metadata.txt
	- umfmt_bile_acids.txt

```
# combine with meta:
combined<-read.table(file="umfmt_summary.txt", header=TRUE)
meta<-read.table(file="umfmt_metadata.txt", header=TRUE)
bas<-read.table(file="supervised_bile_acids.txt", header=TRUE)

div<-combined[, c(2, 11:32)]
meta2<-merge(meta, div, by="sampleID", all.x=TRUE)
data<-merge(meta2, bas, by.x="sampleID", by.y="Samples")
#write.table(data, file="umfmt_summary.w.bile.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

```

###### combining all data (metadata, alpha metrics, relative abundance of OTUs, and metabolites) --> umfmt_allmeasures.txt
	- umfmt_summary.txt
	- umfmt_metadata.txt
	- umfmt_summary.w.bile.txt
	- umfmt_SCFA.txt
	- umfmt_suclac.txt
	- mothurfiles/allhumos4.final.0.03.pick.0.03.pick.0.03.filter.shared

```
combined.raw<-read.table(file="umfmt_summary.txt", header=TRUE)
	combined<-combined.raw[, c(2, 8:ncol(combined.raw))]
meta.raw<-read.table(file="umfmt_metadata.txt", header=TRUE)
	meta<-meta.raw[!is.na(meta.raw$patientID), ]
bas.raw<-read.table(file="umfmt_summary.w.bile.txt", header=TRUE)
	bas<-bas.raw[, c(1, 32:56)]
scfa.raw<-read.table(file="../SCFA/umfmt_SCFA.txt", header=T)
	scfa<-scfa.raw[, c(2, 12:ncol(scfa.raw))]
scfa2.raw<-read.table(file="../SCFA/umfmt_suclac.txt", header=T)
	scfa2<-scfa2.raw[, c(3, 5, 8)]
shared<-read.table(file="mothurfiles/allhumos4.final.0.03.pick.0.03.pick.0.03.filter.shared", header=TRUE, row.names=2)
	otu<-subset(shared, select =-c(label, numOtus) )
	otu.rel<-otu/rowSums(otu)
	otu.rel$sampleID<-rownames(otu.rel)
	
# combine all together (note: some info is missing for some...)
combo<-merge(meta, combined, by="sampleID", all.x=TRUE)
combo2<-merge(combo, bas, by="sampleID", all.x=TRUE)
combo3<-merge(combo2, scfa, by="sampleID", all.x=TRUE)
combo4<-merge(combo3, scfa2, by="sampleID", all.x=TRUE)
combo5<-merge(combo4, otu.rel, by="sampleID", all.x=TRUE)
#write.table(combo5, file="umfmt.allmeasures.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

```

###### umfmt_genfrac2p.all_w.meta.txt
	- allhumos4.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.tax.summary
	- allhumos4_all.genera.txt	#created previously with mouse data)
	- umfmt_summary.txt (as metadata)
	
```
# step 1: create a 'phylotype' file with phylum levels

# read in mothur file; get genus-level assignments and assign phyla
tax<-read.table(file="mothurfiles/allhumos4.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.tax.summary", header=TRUE)

# get phylum designations for level 6 (genera) rows, and curate levels for graphing (later):
tax2<-tax[which(tax$taxlevel==2), ]
tax2[, c("rankID", "taxon")]
tax6<-tax[which(tax$taxlevel==6), ]
tax6$rankID<-gsub("^0.1.1.*", "20_Archaea_unclassified", tax6$rankID)
tax6$rankID<-gsub("^0.1.3.*", "20_Euryarchaeota", tax6$rankID)
tax6$rankID<-gsub("^0.1.2.*", "20_Crenarchaeota", tax6$rankID)
tax6$rankID<-gsub("^0.1.4.*", "20_Thaumarchaeota", tax6$rankID)
tax6$rankID<-gsub("^0.2.1\\..*", "10_Acidobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.2\\..*", "04_Actinobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.3\\..*", "20_Atribacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.4\\..*", "20_BRC1", tax6$rankID)
tax6$rankID<-gsub("^0.2.5\\..*", "11_Unclassified", tax6$rankID)
tax6$rankID<-gsub("^0.2.6\\..*", "01_Bacteroidetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.7\\..*", "20_Candidatus_Saccharibacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.8\\..*", "20_Chlamydiae", tax6$rankID)
tax6$rankID<-gsub("^0.2.9\\..*", "20_Chloroflexi", tax6$rankID)
tax6$rankID<-gsub("^0.2.10..*", "20_Cloacimonetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.11..*", "20_Cyanobacteria/Chloroplast", tax6$rankID)
tax6$rankID<-gsub("^0.2.12..*", "20_Deferribacteres", tax6$rankID)
tax6$rankID<-gsub("^0.2.13..*", "20_Deinococcus-Thermus", tax6$rankID)
tax6$rankID<-gsub("^0.2.14..*", "20_Fibrobacteres", tax6$rankID)
tax6$rankID<-gsub("^0.2.15..*", "02_Firmicutes", tax6$rankID)
tax6$rankID<-gsub("^0.2.16..*", "06_Fusobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.17..*", "20_Gemmatimonadetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.18..*", "20_Ignavibacteriae", tax6$rankID)
tax6$rankID<-gsub("^0.2.19..*", "20_Latescibacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.20..*", "20_Lentisphaerae", tax6$rankID)
tax6$rankID<-gsub("^0.2.21..*", "20_Microgenomates", tax6$rankID)
tax6$rankID<-gsub("^0.2.22..*", "20_Nitrospinae", tax6$rankID)
tax6$rankID<-gsub("^0.2.23..*", "20_Nitrospinae", tax6$rankID)
tax6$rankID<-gsub("^0.2.24..*", "20_Parcubacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.25..*", "20_Planctomycetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.26..*", "03_Proteobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.27..*", "09_Spirochaetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.28..*", "08_Synergistetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.29..*", "07_Tenericutes", tax6$rankID)
tax6$rankID<-gsub("^0.2.30..*", "20_Thermotogae", tax6$rankID)
tax6$rankID<-gsub("^0.2.31..*", "05_Verrucomicrobia", tax6$rankID)
tax6$rankID<-gsub("^0.2.32..*", "20_candidate_division_WPS-2", tax6$rankID)
tax6$rankID<-gsub("^0.3.1..*", "11_unknown_unclassified", tax6$rankID)
colnames(tax6)[2]<-"phylum"
	# remove samples w/ <5000:
subtax6<-subset(tax6, select=-c(taxlevel, daughterlevels))
subtax6<-subtax6[order(subtax6$phylum, -subtax6$total), ]
taxmatrix<-subtax6[, c(4:ncol(subtax6))]
duplicated(subtax6$taxon)			#identify any duplicated taxon names
subtax6$taxon<-as.character(subtax6$taxon)
subtax6$taxon[402]<-"Actinobacteria_unclassified2"
subtax6$taxon<-as.factor(subtax6$taxon)
rownames(taxmatrix)<-subtax6$taxon
genera<- taxmatrix[, colSums(taxmatrix)>5000,]
	# get rel. abund fraction:
genmatrix<-as.data.frame(t(genera))
genera.fr<-genmatrix/rowSums(genmatrix)*100
genus.fr<-t(genera.fr)
all.genera<-cbind(subtax6[1:3], genus.fr)
#write.table(all.genera, file="allhumos4_all.genera.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# note: this file included a larger dataset, and will be filtered to show the samples relevant to this project

### step 2: combine with metadata and filter out data relevant to human stuff:

# read in files
combined<-read.table(file="umfmt_summary.txt", header=TRUE)
genbar<-read.table(file="allhumos4_all.genera.txt", header=TRUE, row.names=2)
meta<-combined[, 1:9]

# filter out human samples only:
rownames(genbar)<-genbar$taxon
humans<-genbar[ , grep("^D.*|^R.*", colnames(genbar))]	# pretty close--just remove the one mouse sample
humans<-humans[, !colnames(humans) %in% c('R_WT_1', 'Rwt_2')]

# now filter to 1 or 2%:
phyla<-subtax6[1:3]
genus1<- humans[rowSums(humans>=1)>=1,]
namelist<-as.character(rownames(genus1))
phyla1p<-phyla[phyla$taxon %in% namelist, ]
genera1<-cbind(phyla1p, genus1)
	# get top 2%
genus2<- humans[rowSums(humans>=2)>=2,]
namelist<-as.character(rownames(genus2))
phyla2p<-phyla[phyla$taxon %in% namelist, ]
genera2<-cbind(phyla2p, genus2)

# add some colors to the levels, by phyla:
summary(as.factor(genera2$phylum))
color<-c("darkgreen","green3","lightgreen","seagreen",
		"midnightblue","mediumblue","blue3","blue","dodgerblue4","dodgerblue1","deepskyblue4","deepskyblue1","skyblue3","skyblue","steelblue4","steelblue1","royalblue4","royalblue1","slateblue4","purple3","orchid3","plum4","plum1","pink3","pink","lightpink1","lightpink3","palevioletred4","palevioletred1","magenta4","deeppink4","mediumvioletred","magenta3","magenta1","thistle",
		"yellow2","darkgoldenrod3","goldenrod2","orange2","yellow4", 
		"maroon", "red4", 
		"hotpink", "red", "cyan", "black", "grey67")
genera2<-cbind(phyla2p, color, genus2)
#write.table(genera2, file="umfmt_genfrac2p.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# read in file and combine with meta:
genbar<-read.table(file="umfmt_genfrac2p.txt", header=TRUE, row.names=2)
	rm_g<-subset(genbar, select =-c(phylum, color, total) )
	barg<-as.data.frame(t(rm_g))
	barg$other<-100-rowSums(barg)
	others<-100-colSums(barg)
	barg$sampleID<-rownames(barg)
	col.gen<-c(as.character(genbar$color), "grey47")
	#barg$sampleID<-gsub("X19_", "19_", barg$sampleID)
bar<-merge(meta, barg, by.x=c("seqID"), by.y=c("sampleID"))
#write.table(bar, 'umfmt_genfrac2p.all_w.meta.txt',quote=FALSE,sep="\t", col.names=NA)

# if you want all genera (including the rarer guys), do this:
meta<-read.table(file="../allhumos3/allhumos_metadata_corrected.txt", header=TRUE)
genbar<-read.table(file="allhumos4_all.genera.txt", header=TRUE, row.names=2)
	#genbar5<- genbar[rowSums(genbar[ ,3:ncol(genbar)]>=5)>=5,]
	rm_g<-subset(genbar, select =-c(phylum, total) )
	barg<-as.data.frame(t(rm_g))
	barg$sampleID<-rownames(barg)
bar<-merge(meta, barg, by.x=c("seqID"), by.y=c("sampleID"))
#write.table(bar, 'allhumos4_allgenera_w.meta.txt',quote=FALSE,sep="\t", col.names=NA)

```

###### compiling GEE results from Krishna's code --> fmtrecip_allresults.txt

```
##### Notes on Krishna's results:

# copied following code to produce a 'results' list, listing top OTUs correlating with Butyrate:
library(TeachingDemos)
library(dplyr)
library(tidyr)
library(geepack)
library(pROC)
library(QICpack) #this should be run everytime, like the other packages

# (krishna's code) create the data file with respective donor/recipient pairs was copied to replicate the results per metabolite

#read in the main file that includes everyone, donors/recipients
fmt<-read.table(file = "umfmt.allmeasures.txt", header=T)

#this data frame limits to just recipients
fmtrecip<-fmt[!fmt$human_group=="donor", ]
#need to arrange it
fmtrecip<-arrange(fmtrecip, sampleID)

#now let's try to duplicate the donor columns in each recipient with merges
fmtdonor<-fmt[fmt$human_group=="donor", ]
#need to arrange it
fmtdonor<-arrange(fmtdonor, sampleID)

#first will need to assign a study number that's the same between donors/recipients
fmtrecip$ID<-gsub("R","",fmtrecip$patientID)
fmtdonor$ID<-gsub("D","",fmtdonor$patientID)

#now merge initial donor samples (unhomogenized) with the recipient
fmtdonorunhomogenized<-fmtdonor[fmtdonor$unhomogenized.x=="yes",]
fmtwide<-merge(fmtrecip, fmtdonorunhomogenized, by="ID", all.x=T, all.y=F)
#de-duplicate this just in case
fmtwide<-fmtwide[!duplicated(fmtwide$sampleID.x),]
#rearrange so in order by sample
fmtwide<-arrange(fmtwide, sampleID.x)

### then, for butyrate, must create filtered dataframe excluding the OTUs that do not converge:
drops<-c("Otu000512", "Otu000529", "Otu000728")
fmtrecip<-fmtrecip[,!(names(fmtrecip) %in% drops)]

start<-match("Otu000001", names(fmtrecip))
end<-(length(names(fmtrecip)))

#now let's create the absence / presence of OTUs
catotus<-dplyr::select(fmtrecip, starts_with("Otu"))
catotus<-catotus>0
temp<-dplyr::select(fmtrecip, -contains("Otu"))
temp<-cbind(temp, catotus)
catfmtrecip<-temp
remove(temp)

#will need to rearrange since this scrambled the order
catfmtrecip<-arrange(catfmtrecip, sampleID)

## can use this file going forward

####### Butyrate:


# Butyrate first:
results<-matrix(data=NA, nrow = end-start+1, ncol = 3)
colnames(results)<-c("variable","P_value","coefficient")

for(j in start:end)
{
  variable<-names(catfmtrecip)[j]
  temp<-catfmtrecip[!(is.na(catfmtrecip[,"Butyrate_nmol"]) | is.na(catfmtrecip[,j])),]
  
  results[j-start+1, 1]<-variable
  results[j-start+1, 2]<-summary(geeglm(as.formula(paste("Butyrate_nmol", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Pr(>|W|)"]
  results[j-start+1, 3]<-summary(geeglm(as.formula(paste("Butyrate_nmol", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Estimate"] 
}
#butyrate predicted by individual OTUs
results

sigresults<-as.data.frame(results)
sigresults$P_value<-as.numeric(as.character(sigresults$P_value))
sigresults<-sigresults[sigresults$P_value<.0001, ]

#let's take the top significant results
topresults<-sigresults
topresults$abscoef<-abs(as.numeric(as.character(topresults$coefficient)))
topresults<-arrange(topresults, desc(abscoef))

# krishna goes on to categorize the model--will chat with him later
# for now, create a file with these topresults to merge with  taxonomy file later

topresults$metabolite<-"butyrate"

# read in taxonomic file to merge with this:
taxa<-read.table(file="../16S_UMFMT/allhumos4.taxonomy.names.txt", header=TRUE)
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
taxanames<-taxa[, c("OTU", "taxname", "phylum", "family")]
butresults_classified<-merge(topresults, taxanames, by.x="variable", by.y="OTU")
#write.table(butresults_classified, file="Anna_test/fmtrecip_but.results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

### let's model these for the best relationship:
#now create a frame just with those
subcatfmtrecip<-dplyr::select(catfmtrecip, patientID, sampleID, Butyrate_nmol, Otu000086,	Otu000476,	Otu000374,	Otu000554,	Otu000045,	Otu000312,	Otu000267,	Otu000083,	Otu000332,	Otu000227)
subcatfmtrecip<-subcatfmtrecip[complete.cases(subcatfmtrecip),]
#as always, will need to rearrange by sampleID
subcatfmtrecip<-arrange(subcatfmtrecip, sampleID)

#tried stepwise addition from there to arrive at the following "best model"
catmodel<-geeglm(Butyrate_nmol~Otu000476+	Otu000374 + Otu000312, family=gaussian(), data=subcatfmtrecip, id=patientID, corstr="ar1")
summary(catmodel)

#so this model includes three OTUs each increasing butyrate levels by an
#additive amount when present, in each individual on average, accounting
#for repeated measures
## get modedeled OTUs in a file:

#### not sure how to get outcome from this; will ask Krishna


####### Acetate:

# acetate:
results_acetate<-matrix(data=NA, nrow = end-start+1, ncol = 3)
colnames(results_acetate)<-c("variable","P_value","coefficient")

for(j in start:end)
{
  variable<-names(catfmtrecip)[j]
  temp<-catfmtrecip[!(is.na(catfmtrecip[,"Acetate_nmol"]) | is.na(catfmtrecip[,j])),]
  
  results_acetate[j-start+1, 1]<-variable
  results_acetate[j-start+1, 2]<-summary(geeglm(as.formula(paste("Acetate_nmol", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Pr(>|W|)"]
  results_acetate[j-start+1, 3]<-summary(geeglm(as.formula(paste("Acetate_nmol", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Estimate"] 
}
# predicted by individual OTUs
results_acetate

sigresults<-as.data.frame(results_acetate)
sigresults$P_value<-as.numeric(as.character(sigresults$P_value))
sigresults<-sigresults[sigresults$P_value<.0001, ]

#let's take the top significant results
topresults<-sigresults
topresults$abscoef<-abs(as.numeric(as.character(topresults$coefficient)))
topresults<-arrange(topresults, desc(abscoef))

# add to file:
topresults$metabolite<-"acetate"
aceresults_classified<-merge(topresults, taxanames, by.x="variable", by.y="OTU")
#write.table(aceresults_classified, file="Anna_test/fmtrecip_ace.results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


####### Propionate:

# propionate

results_propionate<-matrix(data=NA, nrow = end-start+1, ncol = 3)
colnames(results_propionate)<-c("variable","P_value","coefficient")

for(j in start:end)
{
  variable<-names(catfmtrecip)[j]
  temp<-catfmtrecip[!(is.na(catfmtrecip[,"Propionate_nmol"]) | is.na(catfmtrecip[,j])),]
  
  results_propionate[j-start+1, 1]<-variable
  results_propionate[j-start+1, 2]<-summary(geeglm(as.formula(paste("Propionate_nmol", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Pr(>|W|)"]
  results_propionate[j-start+1, 3]<-summary(geeglm(as.formula(paste("Propionate_nmol", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Estimate"] 
}
# predicted by individual OTUs
results_propionate

sigresults<-as.data.frame(results_propionate)
sigresults$P_value<-as.numeric(as.character(sigresults$P_value))
sigresults<-sigresults[sigresults$P_value<.0001, ]

#let's take the top significant results
topresults<-sigresults
topresults$abscoef<-abs(as.numeric(as.character(topresults$coefficient)))
topresults<-arrange(topresults, desc(abscoef))

# add to file:
topresults$metabolite<-"propionate"
proresults_classified<-merge(topresults, taxanames, by.x="variable", by.y="OTU")
#write.table(proresults_classified, file="Anna_test/fmtrecip_pro.results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


####### CA:

# CA:

results_CA<-matrix(data=NA, nrow = end-start+1, ncol = 3)
colnames(results_CA)<-c("variable","P_value","coefficient")

for(j in start:end)
{
  #Otu000086 halts, so let's skip it.
    variable<-names(catfmtrecip)[j]
  temp<-catfmtrecip[!(is.na(catfmtrecip[,"CA"]) | is.na(catfmtrecip[,j])),]
  
  if(!variable=="Otu000086")
  {results_CA[j-start+1, 1]<-variable
  results_CA[j-start+1, 2]<-summary(geeglm(as.formula(paste("CA", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Pr(>|W|)"]
  results_CA[j-start+1, 3]<-summary(geeglm(as.formula(paste("CA", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Estimate"] }
}
# predicted by individual OTUs
results_CA

sigresults<-as.data.frame(results_CA)
sigresults$P_value<-as.numeric(as.character(sigresults$P_value))
sigresults<-sigresults[sigresults$P_value<.0001, ]

#let's take the top significant results
topresults<-sigresults
topresults$abscoef<-abs(as.numeric(as.character(topresults$coefficient)))
topresults<-arrange(topresults, desc(abscoef))

# add to file:
topresults$metabolite<-"CA"
CAresults_classified<-merge(topresults, taxanames, by.x="variable", by.y="OTU")
#write.table(CAresults_classified, file="Anna_test/fmtrecip_CA.results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


####### TCA:

# TCA:

results_TCA<-matrix(data=NA, nrow = end-start+1, ncol = 3)
colnames(results_TCA)<-c("variable","P_value","coefficient")

for(j in start:end)
{
  #Otu000086 halts, so let's skip it.
  variable<-names(catfmtrecip)[j]
  temp<-catfmtrecip[!(is.na(catfmtrecip[,"TCA"]) | is.na(catfmtrecip[,j])),]
  
  if(!variable=="Otu000086")
  {results_TCA[j-start+1, 1]<-variable
  results_TCA[j-start+1, 2]<-summary(geeglm(as.formula(paste("TCA", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Pr(>|W|)"]
  results_TCA[j-start+1, 3]<-summary(geeglm(as.formula(paste("TCA", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Estimate"] }
}
# predicted by individual OTUs
results_TCA

sigresults<-as.data.frame(results_TCA)
sigresults$P_value<-as.numeric(as.character(sigresults$P_value))
sigresults<-sigresults[sigresults$P_value<.01, ]

#let's take the top significant results
topresults<-sigresults
topresults$abscoef<-abs(as.numeric(as.character(topresults$coefficient)))
topresults<-arrange(topresults, desc(abscoef))

# add to file:
topresults$metabolite<-"TCA"
TCAresults_classified<-merge(topresults, taxanames, by.x="variable", by.y="OTU")
#write.table(TCAresults_classified, file="Anna_test/fmtrecip_TCA.results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


####### TCDCA

# TCDCA

results_TCDCA<-matrix(data=NA, nrow = end-start+1, ncol = 3)
colnames(results_TCDCA)<-c("variable","P_value","coefficient")

for(j in start:end)
{
  #Otu000086 halts, so let's skip it.
  variable<-names(catfmtrecip)[j]
  temp<-catfmtrecip[!(is.na(catfmtrecip[,"TCDCA"]) | is.na(catfmtrecip[,j])),]
  
  if(!variable=="Otu000086")
  {results_TCDCA[j-start+1, 1]<-variable
  results_TCDCA[j-start+1, 2]<-summary(geeglm(as.formula(paste("TCDCA", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Pr(>|W|)"]
  results_TCDCA[j-start+1, 3]<-summary(geeglm(as.formula(paste("TCDCA", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Estimate"] }
}
# predicted by individual OTUs
results_TCDCA

sigresults<-as.data.frame(results_TCDCA)
sigresults$P_value<-as.numeric(as.character(sigresults$P_value))
sigresults<-sigresults[sigresults$P_value<.01, ]

#let's take the top significant results
topresults<-sigresults
topresults$abscoef<-abs(as.numeric(as.character(topresults$coefficient)))
topresults<-arrange(topresults, desc(abscoef))

# add to file:
topresults$metabolite<-"TCDCA"
TCDCAresults_classified<-merge(topresults, taxanames, by.x="variable", by.y="OTU")
#write.table(TCDCAresults_classified, file="Anna_test/fmtrecip_TCDCA.results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


####### DCA

# DCA:

results_DCA<-matrix(data=NA, nrow = end-start+1, ncol = 3)
colnames(results_DCA)<-c("variable","P_value","coefficient")

for(j in start:end)
{
  #Otu000086 halts, so let's skip it.
  variable<-names(catfmtrecip)[j]
  temp<-catfmtrecip[!(is.na(catfmtrecip[,"DCA"]) | is.na(catfmtrecip[,j])),]
  
  if(!variable=="Otu000086")
  {results_DCA[j-start+1, 1]<-variable
  results_DCA[j-start+1, 2]<-summary(geeglm(as.formula(paste("DCA", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Pr(>|W|)"]
  results_DCA[j-start+1, 3]<-summary(geeglm(as.formula(paste("DCA", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Estimate"] }
}
# predicted by individual OTUs
results_DCA

sigresults<-as.data.frame(results_DCA)
sigresults$P_value<-as.numeric(as.character(sigresults$P_value))
sigresults<-sigresults[sigresults$P_value<.01, ]

#let's take the top significant results
topresults<-sigresults
topresults$abscoef<-abs(as.numeric(as.character(topresults$coefficient)))
topresults<-arrange(topresults, desc(abscoef))

# add to file:
topresults$metabolite<-"DCA"
DCAresults_classified<-merge(topresults, taxanames, by.x="variable", by.y="OTU")
#write.table(DCAresults_classified, file="Anna_test/fmtrecip_DCA.results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


####### LCA

# LCA:

results_LCA<-matrix(data=NA, nrow = end-start+1, ncol = 3)
colnames(results_LCA)<-c("variable","P_value","coefficient")

for(j in start:end)
{
  #Otu000086 halts, so let's skip it.
  variable<-names(catfmtrecip)[j]
  temp<-catfmtrecip[!(is.na(catfmtrecip[,"LCA"]) | is.na(catfmtrecip[,j])),]
  
  if(!variable=="Otu000086")
  {results_LCA[j-start+1, 1]<-variable
  results_LCA[j-start+1, 2]<-summary(geeglm(as.formula(paste("LCA", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Pr(>|W|)"]
  results_LCA[j-start+1, 3]<-summary(geeglm(as.formula(paste("LCA", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Estimate"] }
}
# predicted by individual OTUs
results_LCA

sigresults<-as.data.frame(results_LCA)
sigresults$P_value<-as.numeric(as.character(sigresults$P_value))
sigresults<-sigresults[sigresults$P_value<.01, ]

#let's take the top significant results
topresults<-sigresults
topresults$abscoef<-abs(as.numeric(as.character(topresults$coefficient)))
topresults<-arrange(topresults, desc(abscoef))

# add to file:
topresults$metabolite<-"LCA"
LCAresults_classified<-merge(topresults, taxanames, by.x="variable", by.y="OTU")
#write.table(LCAresults_classified, file="Anna_test/fmtrecip_LCA.results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


####### UDCA

# UDCA:

results_UDCA<-matrix(data=NA, nrow = end-start+1, ncol = 3)
colnames(results_UDCA)<-c("variable","P_value","coefficient")

for(j in start:end)
{
  #Otu000086 halts, so let's skip it.
  variable<-names(catfmtrecip)[j]
  temp<-catfmtrecip[!(is.na(catfmtrecip[,"UDCA"]) | is.na(catfmtrecip[,j])),]
  
  if(!variable=="Otu000086")
  {results_UDCA[j-start+1, 1]<-variable
  results_UDCA[j-start+1, 2]<-summary(geeglm(as.formula(paste("UDCA", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Pr(>|W|)"]
  results_UDCA[j-start+1, 3]<-summary(geeglm(as.formula(paste("UDCA", "~", variable)), family=gaussian(), data=temp, id=patientID, corstr="ar1"))$coefficients[2,"Estimate"] }
}
# predicted by individual OTUs
results_UDCA

sigresults<-as.data.frame(results_UDCA)
sigresults$P_value<-as.numeric(as.character(sigresults$P_value))
sigresults<-sigresults[sigresults$P_value<.05, ]
#only get 8!

#let's take the top significant results
topresults<-sigresults
topresults$abscoef<-abs(as.numeric(as.character(topresults$coefficient)))
topresults<-arrange(topresults, desc(abscoef))

# add to file:
topresults$metabolite<-"UDCA"
UDCAresults_classified<-merge(topresults, taxanames, by.x="variable", by.y="OTU")
#write.table(UDCAresults_classified, file="Anna_test/fmtrecip_UDCA.results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


####### Combine all the results together for a complete file:

allresults<-rbind(butresults_classified, aceresults_classified, proresults_classified, CAresults_classified, TCAresults_classified, TCDCAresults_classified, DCAresults_classified, LCAresults_classified, UDCAresults_classified)
#write.table(allresults, file="Anna_test/fmtrecip_allresults.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#### This final file was used for Table S3

```

