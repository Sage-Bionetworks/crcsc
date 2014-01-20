require(synapseClient)
require(rGithubClient)

## SOURCE IN BACKGROUND FUNCTIONS FROM JG
if( !exists("crcRepo") ){
  crcRepo <- getRepo("Sage-Bionetworks/crcsc")
}
sourceRepoFile(crcRepo, "groups/G/pipeline/JGLibrary.R")
sourceRepoFile(crcRepo, "groups/G/pipeline/subtypePipelineFuncs.R")


#####
## SET UP THE DATA THAT WE MIGHT WANT TO PULL
#####
dataset <- setRefClass("crcDataset", fields=list(exprSynId="character",
                                                 phenoSynId="character"))
phenoObj <- setRefClass("phenoObj", fields=list(data="data.frame",
                                                discreteFields="list",
                                                continuousFields="list",
                                                censoredFields="list"))

## TCGA UPDATED TO USE MERGED DATA (GA AND HISEQ)
coreDatasets <- list(amc_ajccii=dataset(exprSynId="syn2159423",phenoSynId="syn2159427"),
                     tcga_rnaseq=dataset(exprSynId="syn2325328",phenoSynId="syn2325330"),
                     kfsyscc=dataset(exprSynId="syn2169565",phenoSynId="syn2171240"),
                     french=dataset(exprSynId="syn2171434",phenoSynId="syn2171548"),
                     petacc=dataset(exprSynId="syn2175581",phenoSynId="syn2280515"),
                     nki_az=dataset(exprSynId="syn2176657",phenoSynId="syn2176653"),
                     agendia_gse42284=dataset(exprSynId="syn2192792",phenoSynId="syn2192794"),
                     agendia_ico208=dataset(exprSynId="syn2192796",phenoSynId="syn2289240"),
                     agendia_vhb70=dataset(exprSynId="syn2192799",phenoSynId="syn2289239"),
                     mdanderson=dataset(exprSynId="syn2233387",phenoSynId="syn2290781"))

publicDatasets <- list(gse10961=dataset(exprSynId="syn2177194",phenoSynId="syn2177195"),
                       gse13067=dataset(exprSynId="syn2177888",phenoSynId="syn2177889"),
                       gse13294=dataset(exprSynId="syn2177894",phenoSynId="syn2177895"), 
                       gse14333=dataset(exprSynId="syn2181079",phenoSynId="syn2181006"),
                       gse15960=dataset(exprSynId="syn2177199",phenoSynId="syn2177200"),
                       gse17536=dataset(exprSynId="syn2178137",phenoSynId="syn2178136"),
                       gse17537=dataset(exprSynId="syn2178128",phenoSynId="syn2178129"),
                       gse20916=dataset(exprSynId="syn2177899",phenoSynId="syn2177898"), 
                       gse2109=dataset(exprSynId="syn2177169",phenoSynId="syn2177168"),
                       gse23878=dataset(exprSynId="syn2177902",phenoSynId="syn2178063"), 
                       gse37892=dataset(exprSynId="syn2178082",phenoSynId="syn2178089"),
                       gse4107=dataset(exprSynId="syn2177179",phenoSynId="syn2177180"), 
                       gse4183=dataset(exprSynId="syn2177187",phenoSynId="syn2177188"),
                       gse8671=dataset(exprSynId="syn2181088",phenoSynId="syn2181090"))

cellLineDatasets <- list(ccle=dataset(exprSynId="syn2292137"),
                      sanger=dataset(exprSynId="syn2181097"),
                      gsk=dataset(exprSynId="syn2181084"));

patientDatasets <- c(coreDatasets, publicDatasets)
allDatasets <- c(patientDatasets, cellLineDatasets)

getDatanameForExprSynId <- function(synId){
  idx <- which(sapply(allDatasets, function(ds) { ds$exprSynId==synId } ))
  if(length(idx) == 0){
    warning(paste("Unable to find dataset with expr syn id:", synId))
    return("unknown")
  }else{
    return(names(allDatasets)[idx])
  }
}

groupFolders <- list(GroupA="syn2274064",
                     GroupB="syn2274065",
#                      GroupC="syn2274066",
                     GroupD="syn2319015",
                     GroupE="syn2274069",
                     GroupF="syn2274068",
                     GroupG="syn2274063")

#####
## FUNCTION TO GET PHENOTYPE DATA BY PROVIDING DATASET NAME
#####
getPhenoObjs <- function(ds){
  switch(ds,
         agendia_gse42284=(function(){
           tmp <- read.table(synGet(allDatasets$agendia_gse42284$phenoSynId)@filePath, sep=",",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           rownames(tmp) <- tmp$"Sample name"
           new("phenoObj",
               data=tmp,
               discreteFields=list(location="ch1: characteristics: Location",
                                   kras="ch1: characteristics: KRAS",
                                   braf="ch1: characteristics: BRAF",
                                   pik3ca="ch1: characteristics: PIK3CA",
                                   gender="ch1: characteristics: Gender",
                                   stage="ch1: characteristics: Stage",
                                   msi="ch1: characteristics: Microsattelite"),
               continuousFields=list(age="ch1: characteristics: Age at diagnosis"))
         })(),
         
         agendia_ico208=(function(){
           tmp <- read.table(synGet(allDatasets$agendia_ico208$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           rownames(tmp) <- tmp$"AgendiaID"
           new("phenoObj",
               data=tmp,
               discreteFields=list(location="Location",
                                   kras="KRAS",
                                   braf="BRAF",
                                   gender="Gender",
                                   stage="Stage",
                                   msi="MSI-profile"),
               continuousFields=list(age="Age at diagnosis"))
         })(),
         
         agendia_vhb70=(function(){
           tmp <- read.table(synGet(allDatasets$agendia_vhb70$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           rownames(tmp) <- tmp$"AgendiaID"
           new("phenoObj",
               data=tmp,
               discreteFields=list(location="Location",
                                   kras="KRAS",
                                   braf="BRAF",
                                   gender="Gender",
                                   stage="Stage",
                                   msi="MSI-profile"),
               continuousFields=list(age="Age at diagnosis"))
         })(),
         
         french=(function(){
           tmp <- read.table(synGet(allDatasets$french$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           rownames(tmp) <- tmp$id
           new("phenoObj", 
               data=tmp,
               discreteFields=list(location="tumorLocation",
                                   kras="krasStatus",
                                   braf="brafStatus",
                                   tp53="tp53Status",
                                   gender="gender",
                                   stage="stage",
                                   msi="microsatelite"),
               continuousFields=list(age="age"),
               censoredFields=list(rfs=c("rfsMo","rfsStat")))
         })(),
         
         amc_ajccii=(function(){
           tmp <- read.table(synGet(allDatasets$amc_ajccii$phenoSynId)@filePath, sep=",",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           rownames(tmp) <- tmp$ID
           tmp$dfs.event <- tmp$"DFS (event)" == "yes"
           new("phenoObj",
               data=tmp,
               discreteFields=list(location="Location (distal/prox)" ,
                                   site="Location.CRC",
                                   kras="Ras",
                                   braf="BRAF",
                                   tp53="p53",
                                   gender="Sex",
                                   msi="MSI"),
               continuousFields=list(age="age_at_operation (years)"),
               censoredFields=list(rfs=c("DFS (days)","dfs.event")))
         })(),
         
         kfsyscc=(function(){
           tmp <- read.table(synGet(allDatasets$kfsyscc$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           rownames(tmp) <- tmp$id
           tmp$stage <- gsub("(\\d+).*","\\1",tmp$stage)
           tmp$kras <- tmp$krasCodon12 != "None" | tmp$krasCodon13 != "None"
           new("phenoObj",
               data=tmp,
               discreteFields=list(location="tumorLocation" ,
                                   kras="kras",
                                   stage="stage"),
               continuousFields=list(age="age"),
               censoredFields=list(os=c("osMo","osStat")))
         })(),
         
         petacc=(function(){
           tmp <- read.table(synGet(allDatasets$petacc$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           rownames(tmp) <- tmp$ID
           tmp$os.event <- tmp$os.status == "dead"
           tmp$rfs.event <- tmp$rfs.status == "recurred"
           tmp$kras = c("mut","wt")[factor(tmp$KRAS.mut == 'wt')]
           tmp$msi = c("MSI-H","MSS")[factor(tmp$Microsatellite != "MSI")]
           new("phenoObj",
               data=tmp,
               discreteFields=list(gender="Gender",
                                   location="Tumor.location2",
                                   msi="msi",
                                   kras="kras",
                                   braf="BRAF.mut"),
               continuousFields=list(age="Age"),
               censoredFields=list(os=c("os.time","os.event"),
                                   rfs=c("rfs.time","rfs.event")))
         })(),
         
         nki_az=(function(){
           tmp <- read.table(synGet(allDatasets$nki_az$phenoSynId)@filePath, row.names=1,sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE,)
           new("phenoObj", data=tmp,
               discreteFields=list(gender="characteristics_ch1.8",
                                   kras="characteristics_ch1.1",
                                   braf="characteristics_ch1.2",
                                   pik3ca="characteristics_ch1.5",
                                   tp53="characteristics_ch1.4",
                                   pten="characteristics_ch1.6",
                                   msi="characteristics_ch1.7")) 
         })(),
         
         tcga_rnaseq=(function(){
           # TODO add mutation, MSI information
           tmp <- read.table(synGet(allDatasets$tcga$phenoSynId)@filePath, row.names=1,sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE,)
           new("phenoObj", data=tmp,
               discreteFields=list(gender="gender",
                                   stage="stage",
                                   location="tumorLocation",
                                   msi="microsatelite"),
               continuousFields=list(age="age"),
               censoredFields=list(os=c("osMo","osStat")))
         })(),
         
         mdanderson=(function(){
           tmp <- read.table(synGet(allDatasets$mdanderson$phenoSynId)@filePath,sep="\t",header=TRUE,comment="",as.is=TRUE,
                             na.strings=c("NA","NaN","N/A","n/a"))
           tmp <- tmp[!duplicated(tmp$Agendia.ID),]
           rownames(tmp) <- tmp$Agendia.ID
           tmp$os.event <- tmp$"Death.censor..1.dead..2.censored." == 1
           tmp$rfs.event <- tmp$"Relapse.censor..1.yes.2.censored." == 1
           new("phenoObj", data=tmp,
               discreteFields=list(gender="Sex..Male.1..Female.2.",
                                   stage="Stage",
                                   site="Tumor.site",
                                   location="Localization",
                                   pik3ca="PIK3CA_STATUS.SUMMARY",
                                   kras="KRAS_STATUS.SUMMARY",
                                   braf="BRAF_STATUS.SUMMARY",
                                   msi="MSI.status..IHC."),
               continuousFields=list(age="age.at.surgery.date" ),
               censoredFields=list(os=c("OS.days.from.surgery.date","os.event"),
                                   rfs=c("Relapse.free.days.from.surgery.date","rfs.event")))
         })(),
         
           gse13067=(function(){
           tmp <- read.table(synGet(allDatasets$gse13067$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           return (new("phenoObj", data=tmp,
                       discreteFields=list(msi="microsatellite"),
                       continuousFields=list(),
                       censoredFields=list()))
         })(),
         gse13294=(function(){
           tmp <- read.table(synGet(allDatasets$gse13294$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           tmp <- read.table(synGet(allDatasets$gse13294$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           return (new("phenoObj", data=tmp,
                       discreteFields=list(msi="microsatellite"),
                       continuousFields=list(),
                       censoredFields=list()))
           
         })(), 
         gse14333=(function(){
           tmp1 <- read.table(synGet(allDatasets$gse14333$phenoSynId)@filePath, sep="\t",
                              header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           tmp2 <- read.table(synGet("syn2181048")@filePath, sep="\t",
                              header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           tmp <- merge(tmp1, tmp2, by="row.names",all=TRUE)
           rownames(tmp) <- tmp$Row.names
           return (new("phenoObj", data=tmp,
                       discreteFields=list(gender="gender",location="location"),
                       continuousFields=list(),
                       censoredFields=list(rfs=c("dfs_time","dfs_event"))))
           
         })(), 
         gse17536=(function(){
           tmp <- read.table(synGet(allDatasets$gse17536$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           rownames(tmp) <- tmp$sample
           ## has survival time, but missing formatting
           return (new("phenoObj", data=tmp,
                       discreteFields=list(gender="gender",stage="ajcc_stage"),
                       continuousFields=list("age"),
                       censoredFields=list()))
         })(), 
         gse17537=(function(){
           tmp1 <- read.table(synGet(allDatasets$gse17537$phenoSynId)@filePath, sep="\t",
                              header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           tmp2 <- read.table(synGet("syn2178130")@filePath, sep="\t",
                              header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
           tmp <- merge(tmp1, tmp2, by="row.names",all=TRUE)
           rownames(tmp) <- tmp$Row.names
           return (new("phenoObj", data=tmp,
                       discreteFields=list(gender="gender",stage="ajcc_stage"),
                       continuousFields=list(),
                       censoredFields=list(rfs=c("dfs_time","dfs_event"),os=c("overall_time","overall_event"))))
         })(), 
         gse20916=(function(){
           tmp <- read.table(synGet(allDatasets$gse20916$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE,as.is=TRUE)
           tmp <- tmp[!grepl("norma",tmp$tissue),]
           return (new("phenoObj", data=tmp,
                       discreteFields=list(gender="gender"),
                       continuousFields=list(),
                       censoredFields=list()))
         })(), 
         gse2109=(function(){
           tmp <- read.table(synGet(allDatasets$gse2109$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE,as.is=TRUE)
           tmp$stage <- gsub("(\\d).*","\\1", tmp$"Pathological Stage")
           tmp$stage[tmp$stage=="Unknown"] <- NA
           rownames(tmp) <- tmp$sample
           return (new("phenoObj", data=tmp,
                       discreteFields=list(stage="stage"),
                       continuousFields=list(),
                       censoredFields=list()))
         })(), 
         gse23878=(function(){
           tmp <- read.table(synGet(allDatasets$gse23878$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE,as.is=TRUE)
           tmp <- tmp[tmp$sample_type=="tumor",]
           return (new("phenoObj", data=tmp,
                       discreteFields=list(gender="gender"),
                       continuousFields=list(),
                       censoredFields=list()))
         })(), 
         gse37892=(function(){
           tmp <- read.table(synGet(allDatasets$gse37892$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE,as.is=TRUE)
           return (new("phenoObj", data=tmp,
                       discreteFields=list(gender="gender",stage="stage",location="location"),
                       continuousFields=list(age="age"),
                       censoredFields=list()))
         })(), 
      
         gse8671=(function(){
           tmp <- read.table(synGet(allDatasets$gse8671$phenoSynId)@filePath, sep="\t",
                             header=T,na.strings=c("NA","NaN",""),check.names=FALSE,as.is=TRUE)
           tmp <- tmp[tmp$Tissue=="adenoma",]
           rownames(tmp) <- tmp$sample
           return (new("phenoObj", data=tmp,
                       discreteFields=list(location="Location"),
                       continuousFields=list(),
                       censoredFields=list()))
         })(),
         NA)
}


#####
## FUNCTION TO GET EXPRESSION DATA BY PROVIDING DATASET NAME
#####
getExprSet <- function(ds){
  exprSynId <- allDatasets[[ds]]$exprSynId
  phenoSynId <- allDatasets[[ds]]$phenoSynId
  switch(ds,
         agendia_gse42284=toEntrez(exprSynId, NULL, discoverprint_19742_Map, loader=agendia_data_loader, filename="GSE42284_normalized_data_matrix.txt"),
         agendia_ico208=toEntrez(exprSynId, NULL, discoverprint_19742_Map, loader=agendia_data_loader, filename="ICO208_normalized_data.txt"),
         agendia_vhb70=toEntrez(exprSynId, NULL, discoverprint_32627_Map, loader=agendia_data_loader, filename="VHB70_normalized_data.txt"),
         mdanderson=toEntrez(exprSynId, NULL, discoverprint_32627_Map, loader=agendia_data_loader, filename=NULL),
         kfsyscc=toEntrez(exprSynId, NULL, u133plus2Map),
         french=toEntrez(exprSynId, NULL, u133plus2Map),
         amc_ajccii=toEntrez(exprSynId, NULL, u133plus2Map,sep=","),
         nki_az=toEntrez(exprSynId, NULL, u133plus2Map),
         petacc3=toEntrez(exprSynId, NULL, petaccMap),
         tcga_rnaseq=toEntrez(exprSynId, NULL, symbolMap),
         
         gse10961=toEntrez(exprSynId, phenoSynId),
         gse13067=toEntrez(exprSynId, phenoSynId),
         gse13294=toEntrez(exprSynId, phenoSynId),
         gse14333=toEntrez(exprSynId, phenoSynId),
         gse15960=toEntrez(exprSynId, phenoSynId),
         gse17537=toEntrez(exprSynId, phenoSynId),
         gse20916=toEntrez(exprSynId, phenoSynId),
         gse2109=toEntrez(exprSynId, phenoSynId),
         gse23878=toEntrez(exprSynId, phenoSynId),
         gse37892=toEntrez(exprSynId, phenoSynId),
         gse4107=toEntrez(exprSynId, phenoSynId),
         gse4183=toEntrez(exprSynId, phenoSynId),
         gse8671=toEntrez(exprSynId, phenoSynId),
         NA)
}





#####
## FUNCTIONS TO PULL GROUP INFORMATION
#####
cit.changeRange <- function (v, newmin = 0, newmax = 1){
  oldmin <- min(v, na.rm = TRUE)
  oldmax <- max(v, na.rm = TRUE)
  newmin + ((newmax - newmin) * (v - oldmin)/(oldmax - oldmin))
}

groupBHandler <- function(m){ 
  df <- as.data.frame(t(apply(m[,1:6],1,
                              function(z){
                                x<-1-cit.changeRange(z);
                                x/sum(x)}))) #NB: lower value (0)=best class
  colnames(df) <- gsub("distTo(.*)","\\1", colnames(df))
  df
}

groupDHandler <- function(M){
  df <- as.data.frame(t(apply(M,1,
                              function(z){
                                x<-cit.changeRange(z);
                                x/sum(x)})))
  colnames(df) <- colnames(M)
  df
}

groupEHandler <- function(M){
  if("sample_names" %in% colnames(M)){
    rownames(M) <- M$sample_names
    M <- M[,-which(colnames(M) == "sample_names")]
  }
  M[,-which(colnames(M) == "CCS")]
}


getGroupResults <- function(groupId, filter=patientDatasets){
  
  # reimplement using functions below

}
  
getGroupResultId <- function(group, ds){
  parentId <- groupFolders[[group]]
  tmp <- synapseQuery(paste('SELECT id, name FROM entity WHERE parentId=="', parentId, '"', sep=""))
  
  synIds <- lapply(as.list(tmp$entity.name), function(x){ gsub(".*?_(syn.*?)_.*","\\1",x)})
  these <- sapply(synIds, getDatanameForExprSynId)
  
  tmp$entity.id[ which(these == ds) ]
}

getGroupResult <- function(synId, groupId){
  file <- getFileLocation(synGet(synId))
  sep <- ifelse(grepl("\\.csv$", file), ",", "\t")
  pMatrix <- read.table(file, sep=sep, header=T, as.is=T, row.names=1, check.names=FALSE)
  
  pMatrix <- switch(groupId,
                    GroupB = groupBHandler(pMatrix),
                    GroupD = groupDHandler(pMatrix),
                    GroupE = groupEHandler(pMatrix),
                    pMatrix)
  #rownames(pMatrix) <- clean.names(rownames(pMatrix))
  return(pMatrix)
}






