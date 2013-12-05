dataset <- setRefClass("crcDataset", fields=list(exprSynId="character",
                                                 phenoSynId="character"))
phenoObj <- setRefClass("phenoObj", fields=list(data="data.frame",
                                                discreteFields="list",
                                                continuousFields="list",
                                                censoredFields="list"))

coreDatasets <- list(amc_ajccii=dataset(exprSynId="syn2159423",phenoSynId="syn2159427"),
                     tcga_rnaseq=dataset(exprSynId="syn2161141",phenoSynId="syn2165691"),
                     kfsyscc=dataset(exprSynId="syn2169565",phenoSynId="syn2171240"),
                     french=dataset(exprSynId="syn2171434",phenoSynId="syn2171548"),
                     petacc=dataset(exprSynId="syn2175581",phenoSynId="syn2280515"),
                     nki_az=dataset(exprSynId="syn2176657",phenoSynId="syn2176653"),
                     agendia_gse42284=dataset(exprSynId="syn2192792",phenoSynId="syn2192794"),
                     agendia_ico208=dataset(exprSynId="syn2192796",phenoSynId="syn2289240"),
                     agendia_vhb70=dataset(exprSynId="syn2192799",phenoSynId="syn2289239"),
                     mdanderson=dataset(exprSynId="syn2233387",phenoSynId="syn2290781"))

publicDatasets <- list(
                      #gse10961=dataset(exprSynId="syn2177194",phenoSynId="syn2177195"),
                      gse13067=dataset(exprSynId="syn2177888",phenoSynId="syn2177889"),
                      gse13294=dataset(exprSynId="syn2177894",phenoSynId="syn2177895"), 
                      gse14333=dataset(exprSynId="syn2181079",phenoSynId="syn2181006"),
                      #gse15960=dataset(exprSynId="syn2177199",phenoSynId="syn2177200"),
                      gse17536=dataset(exprSynId="syn2178137",phenoSynId="syn2178136"),
                      gse17537=dataset(exprSynId="syn2178128",phenoSynId="syn2178129"),
                      gse20916=dataset(exprSynId="syn2177899",phenoSynId="syn2177898"), 
                      gse2109=dataset(exprSynId="syn2177169",phenoSynId="syn2177168"),
                      gse23878=dataset(exprSynId="syn2177902",phenoSynId="syn2178063"), 
                      gse37892=dataset(exprSynId="syn2178082",phenoSynId="syn2178089"),
                      #gse4107=dataset(exprSynId="syn2177179",phenoSynId="syn2177180"), 
                      gse4183=dataset(exprSynId="syn2177187",phenoSynId="syn2177188"))
                      #gse8671=dataset(exprSynId="syn2181088",phenoSynId="syn2181090"))

getDatanameForExprSynId <- function(synId){
  allDatasets <- c(coreDatasets, publicDatasets)
  idx <- which(sapply(allDatasets, function(ds) { ds$exprSynId==synId } ))
  if(length(idx) == 0){
    warning(paste("Unable to find dataset with expr syn id:", synId))
    return("unknown")
  }else{
    return(names(allDatasets)[idx])
  }
}


##############################
## core data sets
##
corePhenoObjs <- list(
  agendia_gse42284=(function(){
  tmp <- read.table(synGet(coreDatasets$agendia_gse42284$phenoSynId)@filePath, sep=",",
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
  tmp <- read.table(synGet(coreDatasets$agendia_ico208$phenoSynId)@filePath, sep="\t",
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
  tmp <- read.table(synGet(coreDatasets$agendia_vhb70$phenoSynId)@filePath, sep="\t",
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
  tmp <- read.table(synGet(coreDatasets$french$phenoSynId)@filePath, sep="\t",
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
  tmp <- read.table(synGet(coreDatasets$amc_ajccii$phenoSynId)@filePath, sep=",",
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
  tmp <- read.table(synGet(coreDatasets$kfsyscc$phenoSynId)@filePath, sep="\t",
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
  tmp <- read.table(synGet(coreDatasets$petacc$phenoSynId)@filePath, sep="\t",
                    header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
  rownames(tmp) <- tmp$ID
  tmp$os.event <- tmp$os.status == "dead"
  tmp$rfs.event <- tmp$rfs.status == "recurred"
  new("phenoObj",
         data=tmp,
         discreteFields=list(gender="Gender",
                             tstage="Tstage",
                             nstage="Nstage",
                             location="Tumor.location",
                             msi="Microsatellite",
                             kras="KRAS.mut",
                             braf="BRAF.mut"),
         continuousFields=list(age="Age"),
          censoredFields=list(os=c("os.time","os.event"),
                              rfs=c("rfs.time","rfs.event")))
})(),

nki_az=(function(){
  tmp <- read.table(synGet(coreDatasets$nki_az$phenoSynId)@filePath, row.names=1,sep="\t",
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
  tmp <- read.table(synGet(coreDatasets$tcga$phenoSynId)@filePath, row.names=1,sep="\t",
                    header=T,na.strings=c("NA","NaN",""),check.names=FALSE,)
   new("phenoObj", data=tmp,
                    discreteFields=list(gender="gender",
                                        stage="stage"),
                    continuousFields=list(age="age"),
                    censoredFields=list(os=c("osMo","osStat")))
})(),

mdanderson=(function(){
  tmp <- read.table(synGet(coreDatasets$mdanderson$phenoSynId)@filePath,sep="\t",header=TRUE,comment="",as.is=TRUE,
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
})())

##############################
## public data sets
##


publicPhenoObjs <- list(
    gse10961=(function(){
      tmp <- read.table(synGet(publicDatasets$gse10961$phenoSynId)@filePath, sep="\t",
                        header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
      return (new("phenoObj", data=tmp,
          discreteFields=list(),
          continuousFields=list(),
          censoredFields=list()))
     })(),
    gse13067=(function(){
       tmp <- read.table(synGet(publicDatasets$gse13067$phenoSynId)@filePath, sep="\t",
                         header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
       return (new("phenoObj", data=tmp,
                   discreteFields=list(msi="microsatellite"),
                   continuousFields=list(),
                   censoredFields=list()))
     })(),
    gse13294=(function(){
       tmp <- read.table(synGet(publicDatasets$gse13294$phenoSynId)@filePath, sep="\t",
                         header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
       tmp <- read.table(synGet(publicDatasets$gse13294$phenoSynId)@filePath, sep="\t",
                         header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
       return (new("phenoObj", data=tmp,
                   discreteFields=list(msi="microsatellite"),
                   continuousFields=list(),
                   censoredFields=list()))
       
     })(), 
    gse14333=(function(){
       tmp1 <- read.table(synGet(publicDatasets$gse14333$phenoSynId)@filePath, sep="\t",
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
       tmp <- read.table(synGet(publicDatasets$gse17536$phenoSynId)@filePath, sep="\t",
                          header=T,na.strings=c("NA","NaN",""),check.names=FALSE)
       rownames(tmp) <- tmp$sample
       ## has survival time, but missing formatting
       return (new("phenoObj", data=tmp,
                   discreteFields=list(gender="gender",stage="ajcc_stage"),
                   continuousFields=list("age"),
                   censoredFields=list()))
     })(), 
    gse17537=(function(){
      tmp1 <- read.table(synGet(publicDatasets$gse17537$phenoSynId)@filePath, sep="\t",
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
      tmp <- read.table(synGet(publicDatasets$gse20916$phenoSynId)@filePath, sep="\t",
                         header=T,na.strings=c("NA","NaN",""),check.names=FALSE,as.is=TRUE)
      tmp <- tmp[!grepl("norma",tmp$tissue),]
      return (new("phenoObj", data=tmp,
                  discreteFields=list(gender="gender"),
                  continuousFields=list(),
                  censoredFields=list()))
     })(), 
    gse2109=(function(){
      tmp <- read.table(synGet(publicDatasets$gse2109$phenoSynId)@filePath, sep="\t",
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
      tmp <- read.table(synGet(publicDatasets$gse23878$phenoSynId)@filePath, sep="\t",
                        header=T,na.strings=c("NA","NaN",""),check.names=FALSE,as.is=TRUE)
      tmp <- tmp[tmp$sample_type=="tumor",]
      return (new("phenoObj", data=tmp,
                  discreteFields=list(gender="gender"),
                  continuousFields=list(),
                  censoredFields=list()))
     })(), 
    gse37892=(function(){
      tmp <- read.table(synGet(publicDatasets$gse37892$phenoSynId)@filePath, sep="\t",
                        header=T,na.strings=c("NA","NaN",""),check.names=FALSE,as.is=TRUE)
      return (new("phenoObj", data=tmp,
                  discreteFields=list(gender="gender",stage="stage",location="location"),
                  continuousFields=list(age="age"),
                  censoredFields=list()))
     })(), 
    gse4183=(function(){
       tmp <- read.table(synGet(publicDatasets$gse4183$phenoSynId)@filePath, sep="\t",
                         header=T,na.strings=c("NA","NaN",""),check.names=FALSE,as.is=TRUE)
       tmp <- tmp[tmp$sample_type != "normal",]
       return (new("phenoObj", data=tmp,
                   discreteFields=list(),
                   continuousFields=list(),
                   censoredFields=list()))
     })(), 
    gse8671=(function(){
       tmp <- read.table(synGet(publicDatasets$gse8671$phenoSynId)@filePath, sep="\t",
                         header=T,na.strings=c("NA","NaN",""),check.names=FALSE,as.is=TRUE)
       tmp <- tmp[tmp$Tissue=="adenoma",]
       rownames(tmp) <- tmp$sample
       return (new("phenoObj", data=tmp,
                   discreteFields=list(location="Location"),
                   continuousFields=list(),
                   censoredFields=list()))
     })())
