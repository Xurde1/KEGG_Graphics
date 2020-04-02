
preparateData <- function(datos, sheetName ) {
  colName <- "UCSC_RefGene_Name"
  Gender_lesional<-read.xlsx(datos, sheet =sheetName, skipEmptyRows=TRUE)# hoja Gender_GSE115797
  Gender_lesional<-na.omit(Gender_lesional, cols= colName) #eliminamos las filas sin gen
  Gender_lesional<-subset(Gender_lesional, Gender_lesional$adj.P.Val<0.05) #filtramos solo los significativos
  Gender_lesional_Norm<-scale(Gender_lesional[, 2:3]) #Normalizamos enre 0 y 1
  Gender_lesional_Norm<-cbind(Gender_lesional[, c(1,4)], Gender_lesional_Norm) # unimos los datos
  Gender_lesional_Norm1<-Gender_lesional_Norm %>% group_by(UCSC_RefGene_Name) %>% summarize(LogFC_mean= mean(logFC))
  Gender_lesional_Norm2<-Gender_lesional_Norm %>% group_by(UCSC_RefGene_Name)%>% summarize(adj.Pval_mean= mean(adj.P.Val))
  Gender_lesional_Norm_filter<-cbind( Gender_lesional_Norm1[, c(1,2)], Gender_lesional_Norm2[, 2]) # unimos los datos
}

mergeData <- function(lesional_Norm_filter, psoriasis_Norm_filter) {
  Gender<-merge(lesional_Norm_filter, psoriasis_Norm_filter, by= "UCSC_RefGene_Name")
  #calculamos la media de los Hombres
  Gender$adj.P.Val.Mean<-(Gender$adj.Pval_mean.x + Gender$adj.Pval_mean.y)/2 
  Gender$LogFC.Mean<-(Gender$LogFC_mean.x+Gender$LogFC_mean.y)/2
  Gender<-Gender[, c(1,6,7)] #filtramos las columnas con las medias
}