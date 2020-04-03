  ############### Crear graficos KEGG con datos de metilacion array 450k #########################################################################
  
  ######################################## Preparacion de los datos ########################################################################
  #cargar los datos
  setwd("~/Proyectos/Carmen/KEGG")
  library(openxlsx)
  library(dplyr)
  library(enrichR)
  library(stringr)
  library(KEGGREST)
  library(KEGG.db)
  library(pathview)
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
  setwd("Results/")
  executeKeggGraphic(data_gene, Both, "Both")
  executeKeggGraphic(data_gene, unique_male, "Males")
  executeKeggGraphic(data_gene, unique_female, "Females")
 