require(synapseClient)
require(devtools)

source_gist("https://gist.github.com/brian-bot/6117476")
source_gist("https://gist.github.com/brian-bot/6197439")

## SET THE ROOT WIKI FOR THE TCGA FOLDER
fTcga <- synGet("syn2023932")
wTcga <- WikiPage(title="TCGA Overview", owner=fTcga, markdown="## Overview of the TCGA data\nExpression data from the Illumina Hi-seq machines (RNAseqV2 data frozen by the TCGA Pan Cancer Analysis Group) and their corresponding clinical data were consolidated across the Colon Adenocarcinoma (COAD) and Rectum Adenocarcinoma (READ) cohorts from TCGA.")
tmp <- try(synStore(wTcga), silent=TRUE)
if( class(tmp) == "try-error" ){
  tmp <- synGetWiki(fTcga)
  tmp <- synDelete(tmp)
  wTcga <- synStore(wTcga)
} else{
  wTcga <- tmp
}

## POPULATE SUBPAGES
wTcgaExpr <- knit2synapse("dataQc/tcgaCrcRNAseq.Rmd", owner="syn2023932", parentWikiId=wTcga@properties$id, wikiName="TCGA CRC RNA-Seq Expression Data")
wTcgaClin <- knit2synapse("dataQc/tcgaCrcClinical.Rmd", owner="syn2023932", parentWikiId=wTcga@properties$id, wikiName="TCGA CRC Clinical Data")

