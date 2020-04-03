
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

generateBoth <- function(Females, Males) {
  Both<-merge(Females, Males, by= "UCSC_RefGene_Name") #comunes hombre y mujeres
  Both$LogFC_Mean<-(Both$LogFC.Mean.x + Both$LogFC.Mean.y)/2
  Both$adj.P.Val_Mean<-(Both$adj.P.Val.Mean.x+Both$adj.P.Val.Mean.y)/2
  Both<-Both[, c(1,6,7)]
}

generateDataGene <- function(){
  library(annotables)
  data_gene <- grch38[grch38[, "biotype"]=="protein_coding", ] #obtenemos los datos de la base de datos, para cada gen del genoma
  data_gene<- distinct(data_gene, entrez, .keep_all = TRUE)# filtramos duplicados
}

getKeggMaps <- function() {
  #obtengp los identificadores map
  Identifier <- as.list(KEGGPATHNAME2ID)
  Description<-names(Identifier)
  maps<-cbind(Description,  paste0(Identifier))
}

getPathways <- function() {
  pathways.list<-keggList("pathway", "hsa")
  head(pathways.list)
  names <- substr(names(pathways.list), 6,13)
  terms<-str_remove(pathways.list, pattern = " - Homo sapiens (human)")
  terms<-gsub('.{23}$', "", terms)
  pathways<-cbind(names, terms)
}
  
getRutas <- function(gene) {
  enrichr2<-enrichr(gene, databases ="KEGG_2019_Human") #enrichment igual a enrichr
  Overlap<-enrichr2$KEGG_2019_Human$Overlap
  #Hacer la tabla tipo enrichR
  Num_genes<-str_count(enrichr2$KEGG_2019_Human$Genes, ";")+1
  Term<-enrichr2$KEGG_2019_Human$Term
  Adjusted.P.value <- enrichr2$KEGG_2019_Human$Adjusted.P.value
  Odds.Ratio<-enrichr2$KEGG_2019_Human$Odds.Ratio
  Combined.Score<-enrichr2$KEGG_2019_Human$Combined.Score
  rutas<-cbind(Term,Adjusted.P.value ,Odds.Ratio ,Combined.Score,  Num_genes ,Overlap)
  rutas<-as.data.frame(rutas, stringsAsFactors = FALSE)
  rutas[, c(2:5)]<-sapply(rutas[, c(2:5)], as.numeric)
  rutas<-subset(rutas, rutas$Adjusted.P.value<0.05)#filtro de rutas significativos
}

getMapsWithIdentifiers <- function(rutas, pathways){
  rutas<-as.data.frame(rutas, stringsAsFactors = FALSE)
  rutas[, c(2:5)]<-sapply(rutas[, c(2:5)], as.numeric)
  pathways<-as.data.frame(pathways, stringsAsFactors = FALSE)
  maps<-as.data.frame(maps)
}

getRutasIdentificadasMap <- function (rutas, maps){
  rutas_identificadas_map<-merge(rutas, maps, by.x= "Term", by.y="Description")
  rutas_identificadas_map[, c(2:5)]<-sapply(rutas_identificadas_map[, c(2:5)], as.numeric)
  rutas_identificadas_map<-rutas_identificadas_map[order(rutas_identificadas_map$Adjusted.P.value),]
}

getRutasIdentificadas <- function (rutas, pathways){
  ################################################## tabla de rutas significativas #########################################################
  rutas_identificadas<-merge(rutas, pathways, by.x= "Term", by.y="terms")
  rutas_identificadas[, c(2:5)]<-sapply(rutas_identificadas[, c(2:5)], as.numeric)
  rutas_identificadas[, -c(1, 2, 5,6,7)]<-round(rutas_identificadas[, -c(1, 2, 5,6,7)],2) # 2 decimales
  rutas_identificadas<-rutas_identificadas[order(rutas_identificadas$Adjusted.P.value),]
}

getDatosRepresentar <- function( Both2){
  enterz_name<-Both2[,3]
  rownames(Both2)<-enterz_name
  datos<-Both2[, c(10)]
  names(datos)<-enterz_name
  datos
}

getRutasRepresentar <- function(rutas_identificadas){
  rutas_representar<-5
  rutas_representar<-head(rutas_identificadas$Identifier)
}

getKeggView <- function(datos, rutas_representar, nombre_graficos ){
  pv.out <- pathview(gene.data = datos, pathway.id = rutas_representar,
                     species = "hsa", out.suffix = nombre_graficos,
                     kegg.native = TRUE)
}

getGraphvizView <- function(datos, rutas_representar, nombre_graficos){
  pv.out <- pathview(gene.data = datos, pathway.id = rutas_representar,
                     species = "hsa", out.suffix = nombre_graficos,
                     kegg.native = FALSE, sign.pos =  "bottomleft")
}

executeKeggGraphic<- function(data_gene, Both, nombre_graficos){
  Both2<-merge(data_gene, Both, by.x="symbol", by.y="UCSC_RefGene_Name" ) #obtenemos los datos de los genes con ID enterz
  gene<- as.vector(Both2$symbol)
  ####################################### combinamos con los identificadores KEGG ############################################################
  #obtengo toda la lista de rutas KEGG con sus identificadores
  maps <- getKeggMaps()
  pathways <-getPathways()
  rutas<-getRutas(gene)
  head(rutas)
  
  #tabla para hacer los graficos, de ella se obtienen los identificadores de rutas KEGG
  maps <- getMapsWithIdentifiers(rutas, pathways)
  rutas_identificadas_map <- getRutasIdentificadasMap(rutas, maps)
  rutas_identificadas <- getRutasIdentificadas(rutas, pathways)
  names(rutas_identificadas)[names(rutas_identificadas) == "names"]<-"Identifier"
  head(rutas_identificadas)
  ###########################################################################################################################################
  ####################################### Generacion de los graficos KEGG ###################################################################
  #load data
  rutas_representar<-getRutasRepresentar(rutas_identificadas)
  datos<-getDatosRepresentar(Both2)
  #KEGG view: gene data only
  pv.out <- getKeggView(datos, rutas_representar, nombre_graficos )
  str(pv.out)
  #Graphviz view: gene data only
  pv.out <- getGraphvizView(datos, rutas_representar,  nombre_graficos)
  #result PDF file in current directory
  head(rutas_identificadas)
}