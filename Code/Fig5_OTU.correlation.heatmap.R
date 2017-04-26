#### Figure 6: Correlating OTUs w/ metabolites
##### Anna M. Seekatz
##### 2.3.17


##### Files needed to calculate number of shared OTUs and classifications:
	- directory: /Users/annaseekatz/Box Sync/Projects/UM.FMT/UMFMT_paper/16S_UMFMT
	- fmtrecip_allresults.txt  #this compiled all the OTUs that were found to correlate with the selected metabolites
	- allhumos4.taxonomy.names.txt	#taxonomic information for all OTUs in dataset
	
##### Basic steps:
	- get list of OTUs that correlate (R)
	- select OTU representative sequences from .fasta file of all sequences (mothur)
	- parse out the desired OTUs based on your list (samtools)
	- align OTU sequences to silva (mothur)
	- create a tree out of the aligned sequences (raxml)
	- convert tree to ultrametric tree as a dendrogram (R)
	- make a heatmap of metabolite coefficient values (R)
	- cluster heatmap by phylogenetic dendrogram and add taxonomic classification information (R) 
	

###### Step 1: get the list of OTUs you want to include in your tree

```
# read in file that has list of OTUs
df<-read.table(file="../Krishna Results/Anna_test/fmtrecip_allresults.txt", header=TRUE)

# make a list of all the UNIQUE OTUs that we want to pull into the tree:
otu.list<-unique(df$variable)
#write.table(otu.list, file="../Krishna Results/Anna_test/corrOTUs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
	# this file includes all the OTUs correlated with each metabolite

# if you wanted to select only certain groups of metabolites, you could define each group (SCFA, BA1, or BA2)
scfa.list<-df[df$metabolite %in% c("acetate", "propionate", "butyrate"), c("variable")]
scfa.list<-sort(unique(scfa.list))
#write.table(scfa.list, file="../Krishna Results/Anna_test/scfa.corrOTUs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
ba1.list<-df[df$metabolite %in% c("CA", "TCA", "TCDCA"), unique("variable")]
ba1.list<-sort(unique(ba1.list))
#write.table(ba1.list, file="../Krishna Results/Anna_test/ba1.corrOTUs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
ba2.list<-df[df$metabolite %in% c("DCA", "LCA", "UDCA"), unique("variable")]
ba2.list<-sort(unique(ba2.list))
#write.table(ba2.list, file="../Krishna Results/Anna_test/ba2.corrOTUs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

```

###### Step 2: select representative OTU sequences and combine into a file

######## Note: some of the terminal commands here are specific to HPC used at U of Michigan

```
## ensure that you have all the files you will need to select OTUs and align
## in terminal:
module load mothur
module load samtools
module load raxml
cd /home/aseekatz/humoFMT/allhums4/UMFMT/corrOTU_tree
ln -s /home/aseekatz/humoFMT/allhums4/allhumos4.trim.contigs.good.unique.fasta ./			# this file is the fasta file created immediately before alignment in the mothur pipeline
ln -s /home/aseekatz/humoFMT/allhums4/allhumos4.trim.contigs.good.names ./					# needed for mothur to select OTUs
ln -s /home/aseekatz/humoFMT/allhums4/allhumos4.final.list ./								# the file that lists the sequences within your defined OTU clusters

# if you want to taxonomically classify the OTUs (again):
cp /scratch/youngvi_fluxm/aseekatz/Joyce_Wang/silva.v4.fasta.gz ./ 	
ln -s /home/aseekatz/mothur_files/trainset10_082014.rdp/trainset10* ./

## in mothur:
#mothur v.1.39.0
get.oturep(list=allhumos4.final.list, method=abundance, label=0.03, name=allhumos4.trim.contigs.good.names, fasta=allhumos4.trim.contigs.good.unique.fasta)
	# output: allhumos4.final.0.03.rep.fasta
	# this can be used for ANY of the subsets ot OTU lists below (samtools picks from this list)

## back in terminal:
# the command to choose the seqs that we want cannot have white spaces in the name, so replace them
cp allhumos4.final.0.03.rep.fasta corrOTU.fasta
perl -pi -e "s/M0.*Otu/Otu/g" corrOTU.fasta
perl -pi -e "s/\|.*$//g" corrOTU.fasta

# now, use samtools to select out the desired OTU sequences:
samtools faidx corrOTU.fasta Otu000002 Otu000026 Otu000029 Otu000036 Otu000045 Otu000052 Otu000053 Otu000060 Otu000062 Otu000083 Otu000086 Otu000096 Otu000102 Otu000107 Otu000124 Otu000142 Otu000157 Otu000162 Otu000175 Otu000204 Otu000227 Otu000232 Otu000237 Otu000239 Otu000243 Otu000248 Otu000251 Otu000260 Otu000264 Otu000267 Otu000275 Otu000284 Otu000312 Otu000323 Otu000332 Otu000334 Otu000339 Otu000355 Otu000374 Otu000390 Otu000473 Otu000476 Otu000493 Otu000514 Otu000521 Otu000526 Otu000554 Otu000565 Otu000588 Otu000613 Otu000625 Otu000631 Otu000688 Otu000722 Otu000004 Otu000057 Otu000087 Otu000093 Otu000283 Otu000303 Otu000313 Otu000464 Otu000494 Otu000496 Otu000566 Otu000567 Otu000579 Otu000595 Otu000697 Otu000784 Otu000187 Otu000190 Otu000268 Otu000269 Otu000343 Otu000350 Otu000400 Otu000412 Otu000416 Otu000436 Otu000447 Otu000449 Otu000504 Otu000509 Otu000553 Otu000561 Otu000573 Otu000580 Otu000582 Otu000589 Otu000592 Otu000610 Otu000614 Otu000617 Otu000622 Otu000628 Otu000632 Otu000654 Otu000663 Otu000671 Otu000699 Otu000701 Otu000711 Otu000717 Otu000731 Otu000749 Otu000769 Otu000778 Otu000782 Otu000792 Otu000049 Otu000064 Otu000137 Otu000139 Otu000152 Otu000194 Otu000195 Otu000294 Otu000333 Otu000387 Otu000389 Otu000523 Otu000751 Otu000007 Otu000066 Otu000070 Otu000131 Otu000218 Otu000245 Otu000265 Otu000308 Otu000356 Otu000375 Otu000379 Otu000382 Otu000401 Otu000411 Otu000445 Otu000497 Otu000508 Otu000532 Otu000538 Otu000557 Otu000563 Otu000630 Otu000647 Otu000659 Otu000668 Otu000684 Otu000704 Otu000027 Otu000100 Otu000318 Otu000367 Otu000576 Otu000010 Otu000043 Otu000249 Otu000373 Otu000393 Otu000450 Otu000472 Otu000073 Otu000079 Otu000119 Otu000127 Otu000141 Otu000176 Otu000193 Otu000200 Otu000278 Otu000302 Otu000320 Otu000422 Otu000429 Otu000483 Otu000513 Otu000594 Otu000642 Otu000738 > all.corrOTU.fasta
# for the grouped OTU lists (based on your selected list above):
#samtools faidx corrOTU.fasta Otu000002 Otu000004 Otu000026 Otu000029 Otu000036 Otu000045 Otu000052 Otu000053 Otu000057 Otu000060 Otu000062 Otu000083 Otu000086 Otu000087 Otu000093 Otu000096 Otu000102 Otu000107 Otu000124 Otu000142 Otu000157 Otu000162 Otu000175 Otu000187 Otu000190 Otu000204 Otu000227 Otu000232 Otu000237 Otu000239 Otu000243 Otu000248 Otu000251 Otu000260 Otu000264 Otu000267 Otu000268 Otu000269 Otu000275 Otu000283 Otu000284 Otu000303 Otu000312 Otu000313 Otu000323 Otu000332 Otu000334 Otu000339 Otu000343 Otu000350 Otu000355 Otu000374 Otu000390 Otu000400 Otu000412 Otu000416 Otu000436 Otu000447 Otu000449 Otu000464 Otu000473 Otu000476 Otu000493 Otu000494 Otu000496 Otu000504 Otu000509 Otu000514 Otu000521 Otu000526 Otu000553 Otu000554 Otu000561 Otu000565 Otu000566 Otu000567 Otu000573 Otu000579 Otu000580 Otu000582 Otu000588 Otu000589 Otu000592 Otu000595 Otu000610 Otu000613 Otu000614 Otu000617 Otu000622 Otu000625 Otu000628 Otu000631 Otu000632 Otu000654 Otu000663 Otu000671 Otu000688 Otu000697 Otu000699 Otu000701 Otu000711 Otu000717 Otu000722 Otu000731 Otu000749 Otu000769 Otu000778 Otu000782 Otu000784 Otu000792 > scfa.corrOTU.fasta
#samtools faidx corrOTU.fasta Otu000004 Otu000007 Otu000027 Otu000049 Otu000053 Otu000064 Otu000066 Otu000070 Otu000093 Otu000100 Otu000131 Otu000137 Otu000139 Otu000152 Otu000175 Otu000187 Otu000190 Otu000194 Otu000195 Otu000204 Otu000218 Otu000232 Otu000245 Otu000248 Otu000264 Otu000265 Otu000267 Otu000268 Otu000294 Otu000308 Otu000312 Otu000313 Otu000318 Otu000332 Otu000333 Otu000334 Otu000339 Otu000343 Otu000350 Otu000355 Otu000356 Otu000367 Otu000375 Otu000379 Otu000382 Otu000387 Otu000389 Otu000390 Otu000401 Otu000411 Otu000416 Otu000445 Otu000447 Otu000449 Otu000464 Otu000476 Otu000496 Otu000497 Otu000504 Otu000508 Otu000509 Otu000521 Otu000523 Otu000526 Otu000532 Otu000538 Otu000553 Otu000557 Otu000561 Otu000563 Otu000566 Otu000567 Otu000573 Otu000576 Otu000579 Otu000582 Otu000589 Otu000595 Otu000610 Otu000614 Otu000617 Otu000622 Otu000625 Otu000628 Otu000630 Otu000632 Otu000647 Otu000654 Otu000659 Otu000663 Otu000668 Otu000671 Otu000684 Otu000688 Otu000697 Otu000699 Otu000701 Otu000704 Otu000717 Otu000731 Otu000749 Otu000751 Otu000769 Otu000782 Otu000784 Otu000792 > ba1.corrOTU.fasta
#samtools faidx corrOTU.fasta Otu000010 Otu000026 Otu000029 Otu000036 Otu000043 Otu000053 Otu000062 Otu000064 Otu000070 Otu000073 Otu000079 Otu000102 Otu000119 Otu000124 Otu000127 Otu000139 Otu000141 Otu000142 Otu000152 Otu000157 Otu000162 Otu000175 Otu000176 Otu000187 Otu000193 Otu000194 Otu000200 Otu000204 Otu000227 Otu000232 Otu000237 Otu000249 Otu000251 Otu000260 Otu000264 Otu000265 Otu000275 Otu000278 Otu000284 Otu000294 Otu000302 Otu000312 Otu000320 Otu000332 Otu000333 Otu000343 Otu000350 Otu000373 Otu000374 Otu000379 Otu000387 Otu000389 Otu000390 Otu000393 Otu000400 Otu000401 Otu000411 Otu000412 Otu000416 Otu000422 Otu000429 Otu000445 Otu000447 Otu000449 Otu000450 Otu000472 Otu000476 Otu000483 Otu000493 Otu000496 Otu000508 Otu000509 Otu000513 Otu000514 Otu000521 Otu000523 Otu000553 Otu000554 Otu000561 Otu000566 Otu000567 Otu000573 Otu000579 Otu000582 Otu000589 Otu000592 Otu000594 Otu000595 Otu000610 Otu000613 Otu000614 Otu000617 Otu000622 Otu000631 Otu000632 Otu000642 Otu000654 Otu000663 Otu000671 Otu000684 Otu000697 Otu000699 Otu000701 Otu000704 Otu000717 Otu000731 Otu000738 Otu000749 Otu000751 Otu000769 Otu000782 Otu000784 Otu000792 > ba2.corrOTU.fasta

# you now have a file that has unaligned sequences of your selected OTUs!

## Note: at this point, you could use mothur to classify your OTUs:
classify.seqs(fasta=all.corrOTU.fasta, template=trainset10_082014.rdp.fasta, taxonomy=trainset10_082014.rdp.tax)
	# output: all.corrOTU.filter.rdp.wang.taxonomy
	# this will have the taxonomic classifications of each OTU

```

##### Step 3: Align sequences
	- you have multiple options for aligning your sequences
	- we will use both mothur and a program called mafft

```
module load mafft
module load mothur

## if using mothur to align (in mothur):
align.seqs(fasta=all.corrOTU.fasta, reference=silva.v4.fasta, processors=2)
filter.seqs(fasta=all.corrOTU.align, vertical=T, trump=., processors=2)
	# output: all.corrOTU.filter.fasta
	
# if you want to do it for specific types of metabolites (smaller lists):
align.seqs(fasta=scfa.corrOTU.fasta, reference=silva.v4.fasta, processors=2)
filter.seqs(fasta=scfa.corrOTU.align, vertical=T, trump=., processors=2)
align.seqs(fasta=ba1.corrOTU.fasta, reference=silva.v4.fasta, processors=2)
filter.seqs(fasta=ba1.corrOTU.align, vertical=T, trump=., processors=2)
align.seqs(fasta=ba2.corrOTU.fasta, reference=silva.v4.fasta, processors=2)
filter.seqs(fasta=ba2.corrOTU.align, vertical=T, trump=., processors=2)

## OR, align with mafft (in terminal):
mafft all.corrOTU.fasta > all.corrOTU.aln.fasta
# or
mafft --op 2 all.corrOTU.fasta > all.corrOTU.aln_op2fasta

# obviously, you can play around with the parameters as much as you want

```

##### Step 4: Calculate phylogenetic tree:
	- we will use the program RAxML to create our tree from the aligned sequences

```
module load raxml
module load openmpi/1.10.2/gcc/4.8.5

# on FLUX, need to call it out as this:
raxmlHPC -s all.corrOTU.filter.fasta -n all.corrOTU_mothuraln.out -m GTRCAT -f a -x 123 -N autoMRE -p 456 -T 2


```

##### Step 4: Create heatmap in R

```
# in R
library(phangorn)
library(phytools)
library(plyr)
library(gplots)

## You should now have a best tree that raxml created, so read that in:
raxml.tree<-read.tree("trees/RAxML_bestTree.all.corrOTU_mothuraln.out")
plot(raxml.tree, cex=0.5)

### in order to include this in a heatmap, you must convert it into an ultrametric tree, as a dendrogram
row_dend<-midpoint(raxml.tree) 	#root tree
row_phylo<-row_dend 			# for plotting tree
row_dend<-chronopl(row_dend, lambda = 0.01, tol = 0)
row_dend<- as.dendrogram(as.hclust.phylo(row_dend))
	# this dendrogram clustering can now be used to cluster your heatmap values (the OTUs)
    
### you also want to make a file with the coefficient values by each of the metabolites of interest from the OTU correlation results that you originally made a list with:
# scfa
but<-df[df$metabolite %in% c("butyrate"), c("variable", "coefficient")]
names(but)[2]<-"butyrate"
ace<-df[df$metabolite %in% c("acetate"), c("variable", "coefficient")]
names(ace)[2]<-"acetate"
pro<-df[df$metabolite %in% c("propionate"), c("variable", "coefficient")]
names(pro)[2]<-"propionate"
# ba1
CA<-df[df$metabolite %in% c("CA"), c("variable", "coefficient")]
names(CA)[2]<-"CA"
TCA<-df[df$metabolite %in% c("TCA"), c("variable", "coefficient")]
names(TCA)[2]<-"TCA"
TCDCA<-df[df$metabolite %in% c("TCDCA"), c("variable", "coefficient")]
names(TCDCA)[2]<-"TCDCA"
# ba2
LCA<-df[df$metabolite %in% c("LCA"), c("variable", "coefficient")]
names(LCA)[2]<-"LCA"
DCA<-df[df$metabolite %in% c("DCA"), c("variable", "coefficient")]
names(DCA)[2]<-"DCA"
UDCA<-df[df$metabolite %in% c("UDCA"), c("variable", "coefficient")]
names(UDCA)[2]<-"UDCA"

# now merge, keeping all OTUs:
scfa1<-merge(but, ace, by="variable", all.x=TRUE, all.y=TRUE)
scfa<-merge(scfa1, pro, by="variable", all.x=TRUE, all.y=TRUE)
ba1.1<-merge(CA, TCA, by="variable", all.x=TRUE, all.y=TRUE)
ba1<-merge(ba1.1, TCDCA, by="variable", all.x=TRUE, all.y=TRUE)
ba2.1<-merge(DCA, LCA, by="variable", all.x=TRUE, all.y=TRUE)
#ba2<-merge(ba2.1, UDCA, by="variable", all.x=TRUE, all.y=TRUE)
# merge all
combo1<-merge(scfa, ba1, by="variable", all.x=TRUE, all.y=TRUE)
combo2<-merge(combo1, ba2.1, by="variable", all.x=TRUE, all.y=TRUE)
#write.table(combo2, file="gee.allcorr.matrix.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

# order by scfa.tree OTU tip labels:
map.order<-raxml.tree$tip.label
# ensure that it matches
#cbind(sort(as.character(combo2$variable)), sort(map.order))
m<-combo2[order(match(as.character(combo2$variable), map.order)), ]
#cbind(as.character(m$variable), map.order)

# convert heatmap values to matrix:
rownames(m)<-m$variable
mmatrix<-as.matrix(m[2:9])

# add taxonomic info:
taxa<-read.table(file="allhumos4.taxonomy.names.txt", header=TRUE)
taxanames<-taxa[, c("OTU", "taxname", "phylum", "family")]

## add taxonomic info to row variables and make a color scheme:
otunames<-as.character(rownames(mmatrix))
taxvar<- taxanames[taxanames$OTU %in% otunames,]
	# let's add colors based on taxonomic phyla
	# sort by phyla, family:
filtered.taxa<-taxvar[order(taxvar$phylum, taxvar$family),]
filtered.taxa<-droplevels(filtered.taxa)
	# check that names are correct
filtered.taxa[, c("taxname", "phylum", "family")]
	# change, if necessary:
filtered.taxa$family<-as.character(filtered.taxa$family)
filtered.taxa$family[filtered.taxa$taxname=="Otu374_Actinobacteria_unclassified"]<-"Actinobacteria"
filtered.taxa$family[filtered.taxa$taxname=="Otu476_Proteobacteria_unclassified"]<-"Proteobacteria"
filtered.taxa$family[filtered.taxa$taxname %in% c("Otu237_Firmicutes_unclassified", "Otu275_Firmicutes_unclassified", "Otu449_Firmicutes_unclassified", "Otu521_Clostridia_unclassified", "Otu573_Firmicutes_unclassified", "Otu610_Firmicutes_unclassified", "Otu614_Firmicutes_unclassified", "Otu738_Firmicutes_unclassified")]<-"Firmicutes"
filtered.taxa$family[filtered.taxa$family=="Bacteria"]<-"unclassified"
filtered.taxa$family<-as.factor(filtered.taxa$family)
filtered.taxa<-droplevels(filtered.taxa)
	# then, let's add colors based on phyla, family:
# phyla:
tax.phylum<-unique(filtered.taxa$phylum)
#[1] Actinobacteria        Bacteria_unclassified Bacteroidetes         Firmicutes            Fusobacteria          Proteobacteria        Synergistetes         unclassified         
phy.cols<-c("pink", "black", "chartreuse4", "blue", "red", "yellow", "cyan", "black")
# family:
tax.family<-unique(filtered.taxa$family)
actinos<-colorRampPalette(c("pink", "hotpink1", "hotpink3", "hotpink4"))(4)
uncl<-colorRampPalette(c("black"))(1)
bact<-colorRampPalette(c("green4", "green3", "palegreen3", "palegreen1"))(4)
firm<-colorRampPalette(c("lightskyblue", "skyblue", "royalblue", "mediumpurple"))(14)
fus<-colorRampPalette(c("red"))(1)
prot<-colorRampPalette(c("yellow", "darkgoldenrod1", "orange1"))(7)
other<-colorRampPalette(c("cyan"))(1)
fam.cols<-c(actinos, uncl, bact, firm, fus, prot, other)
# then change the value of the phyla/fam to the appropriate color:
filtered.taxa$famcol<-mapvalues(filtered.taxa$family, from=tax.family, to=fam.cols)
filtered.taxa$phycol<-mapvalues(filtered.taxa$phylum, from=tax.phylum, to=phy.cols)

# note: it is extremely important that your taxonomic variables are in the same order as your heatmap variables!
# match order
taxvar<-filtered.taxa[order(match(as.character(filtered.taxa$OTU), map.order)), ]

# get max/min values for colorscheme:
# summary(mmatrix, na.rm=T)
mmax<-max(mmatrix, na.rm=TRUE)
mmin<-min(mmatrix, na.rm=TRUE)
	# -12 to 25
negbreaks = seq(mmin,0,length.out=25)
posbreaks = seq(0,mmax,length.out=25)

my.breaks = c(seq(mmin,-100,length=3),  								#just remember: your breaks must always = 1 more than the number of colors you are defining!
               seq(-50,-10,length=3),
               seq(-9, -3,length=3),
               seq(-2, 2,length=9),
               seq(3,9,length=3),
               seq(10,50,length=3),
               seq(100,mmax,length=3))
gradient1 = colorpanel( sum( my.breaks[-1]<=0 ), "firebrick", "orange1", "goldenrod2")
gradient2 = colorpanel( sum( my.breaks[-1]> 0 ), "skyblue", "steelblue3", "dodgerblue4")
hm.colors = c(gradient1,gradient2)

# plot
result <- heatmap.2(mmatrix,
                    notecol="black",     
                    density.info="none",
                    key.xlab="",
                    key=FALSE,  
                    trace="none",         
                    margins =c(5,5),
                    col=hm.colors,
                    scale="none",
                    na.color="grey87",   
                    breaks=my.breaks,
                    dendrogram="row",    
                    Colv=F,
                    Rowv=row_dend,
                    keysize=0.5,
                    symkey=FALSE,
                    symm=FALSE,
                    symbreaks=FALSE,
                    srtCol=45,
                    cexRow= 0.5,
                    cexCol = 1,
                    lwid = c(2,2),
                    lhei = c(1,5),
                    labRow=rownames(mmatrix),
                    RowSideColors = as.character(taxvar$famcol)
)
legend("top",legend=tax.family, pt.bg=fam.cols, cex=0.7, pch=22, col="black")
legend("bottom",legend=my.breaks[-1], col=hm.colors, cex=0.7, pch=19)

```