require(synapseClient)
require(devtools)

source_gist("https://gist.github.com/brian-bot/6117476")
source_gist("https://gist.github.com/brian-bot/6197439")

## SET THE ROOT WIKI FOR THE KFSYSCC FOLDER
fKf <- synGet("syn2019114")
wKf <- WikiPage(title="KFSYSCC Overview", owner=fKf, markdown="## Overview of the KFSYSCC data\nKFSYSCC generously provided 322 CRC samples with expression profiling on the Affymetrix U133 Plus 2 platform.")
tmp <- try(synStore(wKf), silent=TRUE)
if( class(tmp) == "try-error" ){
  tmp <- synGetWiki(fKf)
  tmp <- synDelete(tmp)
  wKf <- synStore(wKf)
} else{
  wKf <- tmp
}

## POPULATE SUBPAGES
wKfExpr <- knit2synapse("dataQc/kfsysccExpression.Rmd", owner=fKf, parentWikiId=wKf@properties$id, wikiName="KFSYSCC Expression Data")
wKfClin <- knit2synapse("dataQc/kfsysccClinical.Rmd", owner=fKf, parentWikiId=wKf@properties$id, wikiName="KFSYSCC Clinical Data")

