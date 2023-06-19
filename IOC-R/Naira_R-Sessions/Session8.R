## Révisions sur les fonctions


library(readr)
control_all_markers <- read_csv("control_all_markers.csv")

#Exercice 1.
#fonction qui prend en argument une table et un numéro de cluster et qui retourne
##l'ensemble de gènes du cluster


#Exercice 2a.
## Fonction qui récupére tous les clusters dans lequel ce gène se trouve
##Faire tourner la fonction pour le gène ABHD12

# Exercice 2b. 
## Utiliser la fonction subset pour réaliser la même chose que precedement
##


#Exercice 3 
# Fonction qui récupère le sous tableau pour plusieurs gènes
# Faire tourner la fonction pour vectGenes contenant les gènes IL1B et S100A9
## Aide %in%




ClusterGene <- function(tableau, gene){
return(subset(tableau, select = cluster, subset = gene == gene))}


## Exercice 4. Vérifier la présence de votre gene d'interet est dans une cluster
# fonction qui prend en argument un tableau, un gene et un cluster et qui imprime
##"le gène est présent" ou "pas present"



## Exercice 5. Cas particuliers tester pour le gène CCL2. QUe remarquez vous?



##fonction Grep:permet de trouver les positions des gènes
   
?grep
positionsCCL2 <- grep("CCL2",control_all_markers$gene)

## In regular expression you have special characters :
## ^ means "must begin with
## $ means "must end with"
positionsCCL2 <- grep("^CCL2$",control_all_markers$gene)

# Exercice 6.
## Tester le gene ANXA2 et retrouvez à quel cluster il appartient
#grep ANXA2 lines in the table

##Exercice 7
## Comment modifier un vecteur
# replacer le vecteur des clusters par le character "C" + le numéro du cluster
## sub



# Exercice 8 
# Pour le vect <- c("JANCA","JANA","JAN")
# remplacer seulement le nom du gène JAN par POL

