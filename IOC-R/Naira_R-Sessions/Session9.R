### Session application dplyr & ggplot
# Help https://thinkr.fr/pdf/dplyr-french-cheatsheet.pdf
## N'oubliez pas de loader les librairies

## importer le fichier control_allmarkers au format tibble dans une variable nommée markers


## filtrer le cluster 4 dans c4


## ordonner c4 en fonction de pct.2
toto


## creer un vecteur signature sig1 avec les genes : CLEC2D,CCR7,CD3D,CD6,,CD7,CXCR4,IL7R,,LCK,TRAC,TRBC1

## Aide : utilisation de %in% (in = inclus dans, qui sont dans)


## filtrer markers avec les gènes de la signature dans sig1


## grouper sig1 par cluster



## Calculer la taille (en nombre de genes) de chaque cluster de sig1
## aide n() compte le nombre d'elements
## aide 2 summarise



## Faire les etapes precedents filtrer les genes de sig1, grouper par cluster et compage en 1 ligne de commande avec %>%

  
## Calculer la taille en nombre de genes pour chaque cluster dans markers


## creer une tibble avec les 2 tables precedentes  avec un join



## Calculer la moyenne d'expression de la signature dans chaque cluster




## Creer un vecteur sig2 avec les gènes TRBC2,ARID5B,CAMK4,CCL4, CD3E, CD96, CREM, CUTA, DUSP2,GIMAP7,HLA-DRA  HLA-DRB1 ITM2A,JAML, LEPROTL1

## Calculer la moyennes des pct.1 pour la signature sig2 pour les gènes significatifs (pval> 0.0005 et foldchange>2) par cluster

### nombre de gènes significatif par cluster 

## Proportion de genes differentiellement exprimes dans chaque (en pourcentage)



#######
## ggplot
#######


# Introduction avec https://monashbioinformaticsplatform.github.io/r-more/topics/tidyverse.html
# Lignes de commandes pour le fichier fatsqc a faire ensemble
# Help : CheatSheet 
# https://thinkr.fr/pdf/ggplot2-french-cheatsheet.pdf
# Partie RNASeq 













