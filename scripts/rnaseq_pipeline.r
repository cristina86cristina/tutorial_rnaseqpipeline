##This is a script with all function needed to run the QC 
#RNAseq pipeline 12-03-2018




#Equaliser - if number of genes differs in the TPM files (different Ensembl versions)
# Based on agilp with a few edits for RNAseq

Equaliser <-
  function(input=file.path(system.file(package="agilp"),"input",""),output=file.path(system.file(package="agilp"),"output","")){
    arrays<-dir(input)
    n<-length(arrays)
    
    i<-1
    name<-paste(input,arrays[i], sep="")
    error<-paste("First file is absent so this data is not reliable")
    if(file.exists(name)){
      data<-read.table(file = name, row.names=1, header=FALSE, sep = "\t", fill = TRUE,  stringsAsFactors=FALSE)
    }else message("First file is absent so this data is not reliable")
    
    
    #########################################################################################################
    #first pass runs though all files and finds common demoninator in file data
    for (i in 2:n) {
      name<-paste(input,arrays[i], sep="")
      error<-paste("Arrays number",arrays[i], "is not present")
      if(file.exists(name)){
        data1<-read.table(file = name, row.names=1, header=FALSE, sep = "\t", fill = TRUE,  stringsAsFactors=FALSE)
        common <- merge(data,data1,by.x="row.names",by.y="row.names")
        data <-common[,-3]
        data<-data.frame(data[,-1])
        rownames(data)<-common[,1]
      }else message("Array number ",arrays[i], " is not present")
    }
    
    #########################################################################################################
    #second pass runs though all files again and merges with common demoninator and then saves with same name
    for (i in 1:n) {
      
      name<-paste(input,arrays[i], sep="")
      if(file.exists(name)){
        data1<-read.table(file = name, row.names=1, header=FALSE, sep = "\t", fill = TRUE,  stringsAsFactors=FALSE)
        common <- merge(data,data1,by.x="row.names",by.y="row.names")
        
        outfile<-data.frame(common[,3])
        rownames(outfile)<-common[,1]
        colnames(outfile)<-colnames(data1)
        
        #Output
        array_name<-arrays[i]
        dg<-paste(output,array_name,sep = "")
        write.table(outfile,dg,sep="\t",col.names=FALSE,row.names=TRUE)
      } else message("Arrays number",arrays[i], "is not present")
      #end of for loop inputing files
    }
    #end of function
  }


###Annotation function - to get gene names


annotation = function(Species,datafile,fileout){  
  
  #Read of the output of SARTools
  data2=read.csv(toString(datafile), header =T) ##change csv if you have that
  
  #annotation
  bmdataset = toString(Species)
  mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset= bmdataset)
  #mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset=bmdataset, host = "dec2016.archive.ensembl.org")
  
  #data2$ensembl_gene_id <- data2$V1; data2$V1 <- NULL
  ann <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), filters="ensembl_gene_id",values=data2$ensembl_gene_id,mart = mart)
  #data2$ensembl_gene_id<-data2$Id
  #data2$Id<-NULL
  data_ann<-merge(data2,ann,by="ensembl_gene_id")
  fileout<- paste("Annotated",datafile,sep="")
  write.csv(data_ann,fileout,row.names=FALSE)
}

#fileNames <- Sys.glob("*.txt")
#for (fileName in fileNames) {
 # annotation_meow("hsapiens_gene_ensembl",fileName) #change here if you have another species
  
#}


##function to remove duplicates 

dedup_genesymbol <- function(filein,fileoutdedup){
  tpm <- read.csv(filein)
  tpm <- tpm %>% select(ensembl_gene_id,external_gene_name,description,gene_biotype,everything())
  tpm$mean <- rowMeans(tpm[,c(5:ncol(tpm))])
  tpm <- tpm[order(tpm$external_gene_name, -abs(tpm$mean) ), ] #sort by id and reverse of abs(value)
  tpm <- tpm[ !duplicated(tpm$external_gene_name), ]  # take the first row within each id
  tpm <- tpm[-ncol(tpm)] # remove the mean col
  write.csv(tpm,fileoutdedup,row.names = FALSE)
}
