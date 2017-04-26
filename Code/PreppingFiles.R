#### Prepping mothur output files
##### Anna M. Seekatz
##### 2.3.17
######## directory: /Users/annaseekatz/Box Sync/Projects/UM.FMT/UMFMT_paper/16S_UMFMT

##### Mothur output files and resulting modified files:
	- allhumos4.final.0.03.cons.taxonomy --> allhumos4.taxonomy.names.txt
		- input file lists the taxonomic classifications for each of the OTUs at 0.03 cutoff
		- output file was modified to reflect new OTU names that includes the taxonomy, as well as cleaner classifications
	- allhumos4.final.0.03.pick.0.03.filter.0.03.pick.summary --> umfmt.betasummary.txt
		- input file represents pairwise distances generated in mothur, filtered to include only relevant samples
		- output file adds sample meta data to sampleIDs
	- allhumos4.final.0.03.pick.0.03.pick.0.03.filter.shared, umfmt_metadata.txt --> umfmt_otus.w.meta.txt
		- input .shared file was previously filtered using the specified measures
		- output file combines metadata with OTU counts
	- compiling GEE results from Krishna's code --> fmtrecip_allresults.txt
		
		
		
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

```

####### Butyrate:


```
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

```

####### Acetate:

```
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


```

####### Propionate:

```
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

```

####### CA:

```
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

```

####### TCA:

```
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

```

####### TCDCA

```
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

```

####### DCA

```
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

```

####### LCA

```
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

```

####### UDCA

```
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

```

####### Combine all the results together for a complete file:

```
allresults<-rbind(butresults_classified, aceresults_classified, proresults_classified, CAresults_classified, TCAresults_classified, TCDCAresults_classified, DCAresults_classified, LCAresults_classified, UDCAresults_classified)
#write.table(allresults, file="Anna_test/fmtrecip_allresults.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

```

