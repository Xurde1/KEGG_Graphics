R Script documentation

KEGG graphics


Data loading:

Data set is store in an Excel file. This file named “Results methylation general .xlsx” has 4 sheet from two data bases. Each sheet contain CpG island ID, Adjusted p-value, LogFold-change and USCS gene identifier (symbol).

Pre-procesing:

Each sheet:

-  was read and discard any row without gene identifier (removing rows without gene association).
-  was filtered by adjusted p-value, only rows with < 0.05 was filtered in
-  a normalization value was calculated for columns adjusted p-vale and LogFold-change with scale function
-  values of different CpG island of the same gene identifier were collapsed using the mean

Merge by gender: 

Both data set  from same gender were join by gene identifier and means value of adjusted p-value and LogFold-change were calculated

Final Dataset 

Both: merge of females and males data sets, obtaining data sharing by both sexs
unique_male : data set of data exclusive for males
unique_female: data set of data exclusive for females

Enrichment

Library annotables provide a complete list of gene annotations.
Merge by gene symbol with each final dataset to obtain KEGG identifier of each gene

library enrichR provide the function to perform an enrichment similar to enrichr web.
A table with significant pathways from KEGG was obtained for each final data set. (unique_female was an empty table => none significant pathway in KEGG)

Library KEGG.db and KEGGREST were used to download all KEGG pathways identifiers and assign to each significant pathway.

KEGG graphics

Library pathview provide a function to generate the final graphics. Taking each gene list and LogFold-change as its value of each final data set and tables for each one with significant pathways and their KEGG identifier, the 5 more significant pathway were represented
