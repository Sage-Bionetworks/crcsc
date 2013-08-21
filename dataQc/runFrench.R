require(synapseClient)
require(devtools)

source_gist("https://gist.github.com/brian-bot/6117476")
source_gist("https://gist.github.com/brian-bot/6197439")

## SET THE ROOT WIKI FOLDER
fFrench <- synGet("syn2019116")
wFrench <- WikiPage(title="FRENCH Overview", owner=fFrench, markdown="## Summary from GEO\nFrom a clinical and molecular perspective, colon cancer (CC) is a heterogeneous disease but to date no classification based on high-density transcriptome data has been established. The aim of this study was to build up a robust molecular classification of mRNA expression profiles (Affymetrix U133Plus2) of a large series of 443 CC and to validate it on an independent serie of 123 CC and 906 public dataset. We identified and validated six molecular subtypes in this large cohort as a combination of multiple molecular processes that complement current disease stratification based on clinicopathological variables and molecular markers. The biological relevance of these subtypes was consolidated by significant differences in survival. These insights open new perspectives for improving prognostic models and targeted therapies.\n\n## Design\nAmong a large series of colon cancers collected for the Cartes d'Identite des Tumeurs (CIT) program from the French Ligue Nationale Contre le Cancer (http://cit.ligue-cancer.net), 566 were analyzed for mRNA expression profiles using Affymetrix U133plus2 chip and, among theses, 463 could also be analyzed for DNA alteration profiles using the CGH Array ( CIT-CGHarray V6). The 566 tumors was divided into a discovery dataset of 443 CC and a validation dataset of 123 CC.\n\n## Citation\n###### Marisa L, de Reynies A, Duval A, Selves J et al. Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value. PLoS Med 2013;10(5):e1001453. PMID: 23700391\n\n.")
tmp <- try(synStore(wFrench), silent=TRUE)
if( class(tmp) == "try-error" ){
  tmp <- synGetWiki(fFrench)
  tmp <- synDelete(tmp)
  wFrench <- synStore(wFrench)
} else{
  wFrench <- tmp
}

## POPULATE SUBPAGES
wFrenchExpr <- knit2synapse("dataQc/frenchExpression.Rmd", owner=fFrench, parentWikiId=wFrench@properties$id, wikiName="FRENCH Expression Data")
wFrenchClin <- knit2synapse("dataQc/frenchClinical.Rmd", owner=fFrench, parentWikiId=wFrench@properties$id, wikiName="FRENCH Clinical Data")

