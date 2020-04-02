
preparateData <- function(datos, sheetName ) {
  colName <- "UCSC_RefGene_Name"
  Males_lesional<-read.xlsx(datos, sheet =sheetName, skipEmptyRows=TRUE)# hoja Males_GSE115797
  Males_lesional<-na.omit(Males_lesional, cols= colName) #eliminamos las filas sin gen
  Males_lesional<-subset(Males_lesional, Males_lesional$adj.P.Val<0.05) #filtramos solo los significativos
  Males_lesional_Norm<-scale(Males_lesional[, 2:3]) #Normalizamos enre 0 y 1
  Males_lesional_Norm<-cbind(Males_lesional[, c(1,4)], Males_lesional_Norm) # unimos los datos
  Males_lesional_Norm1<-Males_lesional_Norm %>% group_by(UCSC_RefGene_Name) %>% summarize(LogFC_mean= mean(logFC))
  Males_lesional_Norm2<-Males_lesional_Norm %>% group_by(UCSC_RefGene_Name)%>% summarize(adj.Pval_mean= mean(adj.P.Val))
  Males_lesional_Norm_filter<-cbind( Males_lesional_Norm1[, c(1,2)], Males_lesional_Norm2[, 2]) # unimos los datos
}