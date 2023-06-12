## Objectif 1 : savoir se servir du help :)
## Objectif 2 : se familiariser avec les objets Seurat
## Objectif 3 : savoir manipuler les objets
## Objetctif 4: Faire un plot 

## Creer un repertoire pbmc_data et uploader les 3 fichiers de https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
## Ces fichiers correspondent à l'output généré par un séquençage 10X


### Loader les librairies dplyr et Seural

## Importer les data dans la variable pbmc.data avec la fonction Read10X

## Explorer pbmc.data (regarder les data, utiliser dim)

## utiliser la fonction head pour visualiser les differents éléments
## Pour acceder aux elements, basez vous sur l'annotation de l'objet pbmc.data
## Utiliser @ comme dans l'exemple suivant :

## Creer  l'objet Seurat pbmc à partir des données 10X avec la fonction CreateSeuratObject
## Avec les paramètres suivants : 
## comptages : pbmc.data
## nom du projet : "pbmc3k", 
## nombre minimum de cellules : 3
## Nombre minimum de features : 200


## Explorer l'objet pbmc (toujours utiliser la fonction head)
## 

## Après l'analyse des données (différentes étapes de normalisation et clustering)
## Votre Object Seurat s'enrichit avec différents paramètres
## Importer le fichier pbms_Seurat.rds dans pbmc avec la fonction readRDS


## Explorer les metadata ce nouvel objet Seurat 

#Utiliser la petite icône dans l'onglet view pour récupérer le "chemin" de la variable


## La fonction FindAllMarkers permet de trouver les marqueurs pour chaque cluster à partir d'un objet Seurat préanalysé
## Utiliser la fonction pour retrouver les marqueurs de pbmc
# seulement pour les marqueurs positifs
# avec un pourcentage minimal de detection de 25%
# et un cutoff pour le log fold change de 0.25
## Stocker le résultat dans pbmc_markers



# A partir des metadata de pbmc, recupérer dans la variable a la valeur de nFeature_RNA de la 3eme ligne du tableau de metadata

## Combien y a t-il de clusters identifiés ?
#9

## Combien y a t-il de gènes dans chaque cluster ?


## Plot

# On va visualiser les genes qui varient le plus
# Première étape : crée une table qui contient les données d'intéret :

## On peut récupérer les "highly variable features (ici des gènes)" à partir d'un objet seurat avec la fonction HVFInfo
## Utiliser cette fonction pour récupérer la table hvf
## Cette table contient 3 colonnes et les noms de gènes dans le nom des lignes

# Faire un premier plot plot1 avec 
# en abscisse le log de l'expression moyenne de chaque gène et
# en ordonnée la variance standardisée


# Ajouter une colonne variable avec la valeur 
# "yes" si la variance standardisée est supérieur à 1
# "no" sinon


#ARRETER DE METTRE LE NOM DE VARIABLE POUR LA FAIRE
#FAIRE DES HEAD A LA PLACE

# Refaire plot1 en colorant les points en fonction de l'annotation de la colonne variable 

# Sauvegarder le plot1 dans un fichier pdf et/ou jpeg

#Ajouter à cette table une colonne gène avec le nom des gènes

# Récupérer la sous table avec les 10 gènes qui varient le plus dans top10
# (avec des valeurs de variance standardizées les plus grandes)
# Aide : utiliser slice_max


# Ajouter le nom des gènes qui varient le plus au plot précédement réaliser


## Pour une meilleure ecriture du texte, on utilise la librairie ggrepel dont voici un exemple :
## geom_text_repel permet un décalage "intelligent" des noms des variables
## par rapport aux points qui lui sont associés



## Et pour finir des exemples interessants pour les plots :
## http://www.sthda.com/english/wiki/wiki.php?id_contents=7991





