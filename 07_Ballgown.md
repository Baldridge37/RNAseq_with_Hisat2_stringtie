# Ballgown on R

To finish the RNAseq analysis we need to switch from the command line to R.
This can sometimes be a bit tricky if you haven't got much R experience or it's your first time running it so I'm going to go into a bit more detail than previously.

## Packages

Firstly we need to make sure we have all the packages we need to run the analysis

ballgown is installed through BioConductor https://bioconductor.org/packages/release/bioc/html/ballgown.html 

```
install.packages("BiocManager")
BiocManager::install("ballgown", version = "3.8")
library(ballgown)
```

We will also need some other packages to perform the analysis and handle the data. 
If not shown otherwise packages can be installed using install.packages(c("package_1", "package_2"))

```
library(tidyverse)
library(devtools)

devtools::install_github('alyssafrazee/RSkittleBrewer')
library(RSkittleBrewer)

BiocManager::install("genefilter", version = "3.8")
library(genefilter)
```

We will now be ready to input the data

## Ballgown object

```
#set working directory
setwd("directory/")

#Load sample data
read_tsv("./ballgown/pheno.txt", col_names =  T) -> pheno_data %>% as.data.frame()

#Build ballgown object
bg<-ballgown(dataDir = "ballgown", samplePattern = "sample_", pData = pheno_data)

#filter by expression
bg_filt = subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)
```

Now, in this section you can see that we import data on the samples in the pheno_data.
Basically this can be whatever you want it to be. It should detail your sample names, tissue type, and any other differences between them.
This is important for performing statistical tests later on where we need to say how to group the samples.

As an example here is one from an experiment I performed. Here all samples start with Tap_ and I simply have the sample name and tissue type.
```
Sample	Tissue
Tap_Meiocyte_HGXF003C	Meiocyte
Tap_Meiocyte_HGXF003E	Meiocyte
Tap_Meiocyte_HGXF006A	Meiocyte
Tap_Rosette_HGXF002A	Rosette
Tap_Rosette_HGXF002E	Rosette
Tap_Rosette_HGXF003A	Rosette
Tap_Tapetum_JZXF003A	Tapetum
Tap_Tapetum_JZXF003B	Tapetum
Tap_Tapetum_JZXF003C	Tapetum
```

Now that we've created the ballgown object we can find both differentially expressed transcripts and genes. 
This is done with the inbuilt function stattest. The covariate corresponds to a column of the pheno_data.
If we had extra detail, such as a time series, we could find genes that vary between different time points, or that vary due to library method/size.
Collecting extra information on your samples could be very useful.

```
results_transcripts = stattest(bg_filt, feature="transcript", covariate = "Tissue", getFC=F, meas="FPKM")
#getFC (fold change) must be set to F when you are dealing with more than two samples as fold changes can't be calculated.

results_genes =stattest(bg_filt, feature="gene", covariate = "Tissue", getFC = F, meas="FPKM")

#Add IDs to the transcript results
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filt),
                                 geneIDs=ballgown::geneIDs(bg_filt), results_transcripts)
#Arrange by p value
results_transcripts = arrange(results_transcripts, pval)
results_genes = arrange(results_genes, pval)

#Get gene names for gene results. I don't know why this is more complicated than for transcripts but it is.
indices <- match(results_genes$id, texpr(bg_filt, 'all')$gene_id)
gene_names_for_result <- texpr(bg_filt, 'all')$gene_name[indices]
results_genes <- data.frame(geneNames=gene_names_for_result, results_genes)

#The full table will contain all the genes contained within bg_filt
full_table <- texpr(bg_filt , 'all')

#To get gene-level data use
full_table_gene <- gexpr(bg_filt)
```

We should now have a list of genes/transcripts and their associated p and q values.
We can now very easily extract those that are significant using subset.

```
subset(results_transcripts, results_transcripts$qval<0.05) -> results_transcripts_filt
subset(results_genes, results_genes$qval<0.05) -> results_genes_filt

write_tsv(results_genes_filt, "./ballgown/genes_results_qval.txt")
write_tsv(results_transcripts_filt, "./ballgown/transcripts_results_qval.txt")
```

A list of significantly different genes is a great place to start for RNAseq analysis, but what if we want a bit more information than that?
Currently the results_genes_filt doesn't tell us what tissues the gene is expressed in.
To make the file I'll use the left_join fucntion from tidyverse. Any genes in the significant list will have the data from the full_table joined to it.

```
results_genes_filt %>% left_join(full_table,by=c('id'='gene_id')) -> sig_gene_results_all
```

So we will now have the significant genes with all their data attached. We can look for a gene and find its expression in each sample.
But what if we're only really interested in a gene's average expression in a tissue? Well we can use mutate from the tidyverse to add and calculate new columns.

```
sig_gene_results_all %>% mutate(
  Meiocyte_mean_FPKM = ((FPKM.Tap_Meiocyte_HGXF003C+ FPKM.Tap_Meiocyte_HGXF003E+ FPKM.Tap_Meiocyte_HGXF006A)/3),
  Rosette_mean_FPKM = ((FPKM.Tap_Rosette_HGXF002A + FPKM.Tap_Rosette_HGXF002E + FPKM.Tap_Rosette_HGXF003A)/3),
  Tapetum_mean_FPKM = ((FPKM.Tap_Tapetum_JZXF003A + FPKM.Tap_Tapetum_JZXF003B + FPKM.Tap_Tapetum_JZXF003C)/3)
  ) %>% select(1:14,33:35) %>% mutate(Expressed_in = case_when(
    Meiocyte_mean_FPKM > Rosette_mean_FPKM | Meiocyte_mean_FPKM > Tapetum_mean_FPKM ~ "Meiocyte",
    Rosette_mean_FPKM >  Meiocyte_mean_FPKM | Rosette_mean_FPKM > Tapetum_mean_FPKM ~ "Rosette",
    Tapetum_mean_FPKM >  Meiocyte_mean_FPKM | Tapetum_mean_FPKM > Rosette_mean_FPKM ~ "Tapetum"
  )) -> sig_gene_results_means
  
  write_tsv(sig_gene_results_means, "./ballgown/Tsig_gene_expression_qval.txt")
```

This may look complicated but it's actually performing a series of simple functions before being piped (%>%) to the next.
If we break it down we can see what the file we're building looks like.
The combined data of the significant genes is passed to mutate where new columns are made that have the mean expression.
Then I'm selecting just the columns I want, that is the columns with gene information (1:14) and the new means columns (33:35).
I am then using mutate again to create another column called Expressed_in, to make searching for tissue-enriched genes easier.
By using case_when I fill in the tissue name of the highest expression. So if meiocyte expression is > leaf and (|) tapetum then the row is filled with "Meiocyte".

Now we have much more user friendly lists of significant genes. We can look through these in Excel or whatever program you want and then explore transcripts further back in ballgown.

## Other functions for transcript visualisation

```
# palette
tropical<- c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)

#get FPKMs

fpkm <- texpr(bg_TML, meas="FPKM")
fpkm <- log2(fpkm+1)

#plot distribution of FPKMs
boxplot(fpkm,col=as.numeric(pheno_data$Tissue),las=2,ylab='log2(FPKM+1)')

##get transcript level information
ballgown::transcriptNames(bg)[12]

ballgown::geneNames(bg)[12]

##Get transcript information for a gene

plot(fpkm[12,] ~pheno_data$Tissue, border=c(1,2),
     main=paste(ballgown::geneNames(bg)[12],' : ',
                ballgown::transcriptNames(bg)[12]), pch=19, xlab="Tissue",
     ylab='log2(FPKM+1)')

points(fpkm[12,] ~ jitter(as.numeric(pheno_data$Tissue)),
       col=as.numeric(pheno_data$Tissue))

###Transcript level plots between samples
plotTranscripts(ballgown::geneIDs(bg)[gene_number], bg, main=c('Gene X in samples'), 
                sample=c(''))

#Mean plots #must be MSTRG number
plotMeans('MSTRG.1875', bg_filt, groupvar = "Tissue", legend=FALSE)
```

