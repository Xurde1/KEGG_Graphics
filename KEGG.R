  ############### Crear graficos KEGG con datos de metilacion array 450k #########################################################################
  
  ######################################## Preparacion de los datos ########################################################################
  #cargar los datos
  setwd("~/Proyectos/Carmen/KEGG")
  library(openxlsx)
  library(dplyr)
  source("Code/KEGG_Graphics/function.R")
  datos<-"Input/methylation.xlsx" #excel con los datos
  
  Males_lesional_Norm_filter <- preparateData(datos,"Males GSE115797" )
  Males_psoriasis_Norm_filter <- preparateData(datos,"Males GSE63315" )
  Females_lesional_Norm_filter <- preparateData(datos,"Females les GSE115797")
  Females_psoriasis_Norm_filter <- preparateData(datos,"Females GSE63315" )

      #Hacemos el merge
  
  Males<-mergeData(Males_lesional_Norm_filter, Males_psoriasis_Norm_filter)
  Females<-mergeData(Females_lesional_Norm_filter, Females_psoriasis_Norm_filter)
  
  ############## Merge de ambos sexos y los exclusivos de cada sexo
  Both<-generateBoth(Females, Males)
  unique_male<-anti_join(Males, Both, By = "UCSC_RefGene_Name") #exclusivos de hombres
  unique_female<-anti_join(Females, Both, By = "UCSC_RefGene_Name") #exclusivos de mujere
  ########################################### Aqui ya tenemos los datos preparados ##########################################################
  
  ####################################################Enrichment KEGG para el data frame Both ################################################
  data_gene <- generateDataGene()
  
  Both2<-merge(data_gene, Both, by.x="symbol", by.y="UCSC_RefGene_Name" ) #obtenemos los datos de los genes con ID enterz
  gene<- as.vector(Both2$symbol)
  
  library(enrichR)
  library(stringr)
  library(KEGGREST)
  
  enrichr2<-enrichr(gene, databases ="KEGG_2019_Human") #enrichment igual a enrichr
  #Hacer la tabla tipo enrichR
  Num_genes<-str_count(enrichr2$KEGG_2019_Human$Genes, ";")+1
  Term<-enrichr2$KEGG_2019_Human$Term
  Adjusted.P.value <- enrichr2$KEGG_2019_Human$Adjusted.P.value
  Odds.Ratio<-enrichr2$KEGG_2019_Human$Odds.Ratio
  Combined.Score<-enrichr2$KEGG_2019_Human$Combined.Score
  Overlap<-enrichr2$KEGG_2019_Human$Overlap
  
  rutas<-cbind(Term,Adjusted.P.value ,Odds.Ratio ,Combined.Score,  Num_genes ,Overlap)
  rutas<-as.data.frame(rutas, stringsAsFactors = FALSE)
  rutas[, c(2:5)]<-sapply(rutas[, c(2:5)], as.numeric)
  rutas<-subset(rutas, rutas$Adjusted.P.value<0.05)#filtro de rutas significativos
  head(rutas)
  ####################################### combinamos con los identificadores KEGG ############################################################
  #obtengo toda la lista de rutas KEGG con sus identificadores
  maps <- getKeggMaps()
  
  rutas<-as.data.frame(rutas, stringsAsFactors = FALSE)
  rutas[, c(2:5)]<-sapply(rutas[, c(2:5)], as.numeric)
  pathways<-as.data.frame(pathways, stringsAsFactors = FALSE)
  maps<-as.data.frame(maps)
  
  #tabla para hacer los graficos, de ella se obtienen los identificadores de rutas KEGG
  rutas_identificadas_map<-merge(rutas, maps, by.x= "Term", by.y="Description")
  rutas_identificadas_map[, c(2:5)]<-sapply(rutas_identificadas_map[, c(2:5)], as.numeric)
  rutas_identificadas_map<-rutas_identificadas_map[order(rutas_identificadas_map$Adjusted.P.value),]
  
  ################################################## tabla de rutas significativas #########################################################
  rutas_identificadas<-merge(rutas, pathways, by.x= "Term", by.y="terms")
  rutas_identificadas[, c(2:5)]<-sapply(rutas_identificadas[, c(2:5)], as.numeric)
  rutas_identificadas[, -c(1, 2, 5,6,7)]<-round(rutas_identificadas[, -c(1, 2, 5,6,7)],2) # 2 decimales
  rutas_identificadas<-rutas_identificadas[order(rutas_identificadas$Adjusted.P.value),]
  names(rutas_identificadas)[names(rutas_identificadas) == "names"]<-"Identifier"
  head(rutas_identificadas)
  ###########################################################################################################################################
  ####################################### Generacion de los graficos KEGG ###################################################################
  library(pathview)
  
  
  #load data
  rutas_representar<-5
  rutas_representar<-head(rutas_identificadas_map$V2)
  enterz_name<-Both2[,3]
  rownames(Both2)<-enterz_name
  datos<-Both2[, c(10)]
  names(datos)<-enterz_name
  
  setwd("~/KEGG_graficos/")
  #KEGG view: gene data only
  nombre_graficos<-"Both"
  pv.out <- pathview(gene.data = datos, pathway.id = rutas_representar,
                    species = "hsa", out.suffix = nombre_graficos,
                    kegg.native = TRUE)
  str(pv.out)
  
  #result PNG file in current directory
  
  #Graphviz view: gene data only
  position_sig<-"bottomleft"
  pv.out <- pathview(gene.data = datos, pathway.id = rutas_representar,
                     species = "hsa", out.suffix = nombre_graficos,
                     kegg.native = FALSE, sign.pos = position_sig)
  #result PDF file in current directory
  head(rutas_identificadas)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ############################################################################################################################################
  ####################################################Enrichment KEGG para el data frame Male ################################################
  Males2<-merge(data_gene, unique_male, by.x="symbol", by.y="UCSC_RefGene_Name" ) #obtenemos los datos de los genes con ID enterz
  
  gene_males<- as.vector(Males2$symbol)
  
  
  enrichr_males<-enrichr(gene_males, databases ="KEGG_2019_Human") #enrichment igual a enrichr
  #Hacer la tabla tipo enrichR
  Num_genes<-str_count(enrichr_males$KEGG_2019_Human$Genes, ";")+1
  Term<-enrichr_males$KEGG_2019_Human$Term
  Adjusted.P.value <- enrichr_males$KEGG_2019_Human$Adjusted.P.value
  Odds.Ratio<-enrichr_males$KEGG_2019_Human$Odds.Ratio
  Combined.Score<-enrichr_males$KEGG_2019_Human$Combined.Score
  Overlap<-enrichr_males$KEGG_2019_Human$Overlap
  
  rutas_males<-cbind(Term,Adjusted.P.value ,Odds.Ratio ,Combined.Score,  Num_genes ,Overlap)
  rutas_males<-as.data.frame(rutas_males, stringsAsFactors = FALSE)
  rutas_males[, c(2:5)]<-sapply(rutas_males[, c(2:5)], as.numeric)
  rutas_males<-subset(rutas_males, rutas_males$Adjusted.P.value<0.05)#filtro de rutas significativos
  head(rutas_males)
  ####################################### combinamos con los identificadores KEGG ############################################################
  
  #tabla para hacer los graficos, de ella se obtienen los identificadores de rutas KEGG
  rutas_identificadas_map_males<-merge(rutas_males, maps, by.x= "Term", by.y="Description")
  rutas_identificadas_map_males[, c(2:5)]<-sapply(rutas_identificadas_map_males[, c(2:5)], as.numeric)
  rutas_identificadas_map_males<-rutas_identificadas_map_males[order(rutas_identificadas_map_males$Adjusted.P.value),]
  
  ################################################## tabla de rutas significativas #########################################################
  rutas_identificadas_males<-merge(rutas_males, pathways, by.x= "Term", by.y="terms")
  rutas_identificadas_males[, c(2:5)]<-sapply(rutas_identificadas_males[, c(2:5)], as.numeric)
  rutas_identificadas_males[, -c(1, 2, 5,6,7)]<-round(rutas_identificadas_males[, -c(1, 2, 5,6,7)],2) # 2 decimales
  rutas_identificadas_males<-rutas_identificadas_males[order(rutas_identificadas_males$Adjusted.P.value),]
  names(rutas_identificadas_males)[names(rutas_identificadas_males) == "names"]<-"Identifier"
  head(rutas_identificadas_males)
  ###########################################################################################################################################
  ####################################### Generacion de los graficos KEGG ###################################################################
  
  rutas_representar_males<-head(rutas_identificadas_map_males$V2, rutas_representar)
  enterz_name_males<-Males2[,3]
  rownames(Males2)<-enterz_name_males
  datos_males<-Males2[, c(10)]
  names(datos_males)<-enterz_name_males
  
  
  #KEGG view
  nombre_graficos_males<-"Males"
  pv.out <- pathview(gene.data = datos_males, pathway.id = rutas_representar_males,
                     species = "hsa", out.suffix = nombre_graficos_males,
                     kegg.native = TRUE)
  
  
  
  #Graphviz view
  position_sig<-"bottomleft"
  pv.out <- pathview(gene.data = datos_males, pathway.id = rutas_representar_males,
                     species = "hsa", out.suffix = nombre_graficos_males,
                     kegg.native = FALSE, sign.pos = position_sig)
  
  ############################################################################################################################################
  head(rutas_identificadas_males)
  ####################################################Enrichment KEGG para el data frame Female ################################################
  Females2<-merge(data_gene, unique_female, by.x="symbol", by.y="UCSC_RefGene_Name" ) #obtenemos los datos de los genes con ID enterz
  
  gene_females<- as.vector(Females2$symbol)
  
  
  enrichr_females<-enrichr(gene_females, databases ="KEGG_2019_Human") #enrichment igual a enrichr
  #Hacer la tabla tipo enrichR
  Num_genes<-str_count(enrichr_females$KEGG_2019_Human$Genes, ";")+1
  Term<-enrichr_females$KEGG_2019_Human$Term
  Adjusted.P.value <- enrichr_females$KEGG_2019_Human$Adjusted.P.value
  Odds.Ratio<-enrichr_females$KEGG_2019_Human$Odds.Ratio
  Combined.Score<-enrichr_females$KEGG_2019_Human$Combined.Score
  Overlap<-enrichr_females$KEGG_2019_Human$Overlap
  
  rutas_females<-cbind(Term,Adjusted.P.value ,Odds.Ratio ,Combined.Score,  Num_genes ,Overlap)
  rutas_females<-as.data.frame(rutas_females, stringsAsFactors = FALSE)
  rutas_females[, c(2:5)]<-sapply(rutas_females[, c(2:5)], as.numeric)
  rutas_females<-subset(rutas_females, rutas_females$Adjusted.P.value<0.05)#filtro de rutas significativos
  head(rutas_females)
  ####################################### combinamos con los identificadores KEGG ############################################################
  
  #tabla para hacer los graficos, de ella se obtienen los identificadores de rutas KEGG
  rutas_identificadas_map_females<-merge(rutas_females, maps, by.x= "Term", by.y="Description")
  rutas_identificadas_map_females[, c(2:5)]<-sapply(rutas_identificadas_map_females[, c(2:5)], as.numeric)
  rutas_identificadas_map_females<-rutas_identificadas_map_females[order(rutas_identificadas_map_females$Adjusted.P.value),]
  
  ################################################## tabla de rutas significativas #########################################################
  rutas_identificadas_females<-merge(rutas_females, pathways, by.x= "Term", by.y="terms")
  rutas_identificadas_females[, c(2:5)]<-sapply(rutas_identificadas_females[, c(2:5)], as.numeric)
  rutas_identificadas_females[, -c(1, 2, 5,6,7)]<-round(rutas_identificadas_females[, -c(1, 2, 5,6,7)],2) # 2 decimales
  rutas_identificadas_females<-rutas_identificadas_females[order(rutas_identificadas_females$Adjusted.P.value),]
  names(rutas_identificadas_females)[names(rutas_identificadas_females) == "names"]<-"Identifier"
  head(rutas_identificadas_females)
  ###########################################################################################################################################
  ####################################### Generacion de los graficos KEGG ###################################################################
  
  rutas_representar_females<-head(rutas_identificadas_map_females$V2, rutas_representar)
  enterz_name_females<-Females2[,3]
  rownames(Females2)<-enterz_name_females
  datos_females<-Males2[, c(10)]
  names(datos_females)<-enterz_name_females
  
  
  #KEGG view
  nombre_graficos_females<-"Females"
  pv.out <- pathview(gene.data = datos_females, pathway.id = rutas_representar_females,
                     species = "hsa", out.suffix = nombre_graficos_females,
                     kegg.native = TRUE)
  
  
  
  #Graphviz view
  position_sig<-"bottomleft"
  pv.out <- pathview(gene.data = datos_females, pathway.id = rutas_representar_females,
                     species = "hsa", out.suffix = nombre_graficos_females,
                     kegg.native = FALSE, sign.pos = position_sig)
  
  ############################################################################################################################################
  head(rutas_identificadas_females)