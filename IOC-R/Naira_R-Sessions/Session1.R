## Ceci est mon 1er script R

## Addition 2+3 
2+3 
50*328

## Création de variables
x <- 5
y <- 10.3
z <- "genes"
a <- FALSE
b <- TRUE

# Types de variables
typeof(x)


# Les vecteurs on utilise c()

numeros <- 
x_y <- 

print(x_y)

# exercice :
# Crer un vecteur vect1  avec 3 valeurs numériques
# Creer un vecteur vect2 avec 3 valeurs de type character
# Creer un vecteur vect3 avec les valeurs de vect1 et vect2
# Question : quel est le type des vecteurs



### Valeurs séquentielles

val <- 1:100
?seq
# Récupérer le 5eme element des vecteurs val et val2

# Récupérer les éléments 1,3 et 5 du vecteur val2

## Excercice :
# Créer les vecteurs vect1 et vect2 qui contiennent les valeurs 1,3 et 5 de 2 manières
# Créer le vecteur vect3 qui contient les nombres entiers de 50 à 100
# Creer le vecteur vect4 qui contient les élements 1,3 et 5 du vecteur vect3


## fonctions importantes 

# la taille d'un vecteur
?length



# La valeur NA
vect5 <- c(vect3,NA)
print(vect5)
length(vect5)


# la fonction paste 
## Créer un vecteur vgenes1 de type caractère qui contient les gènes MBNL1, BRCA2

## Créer un vecteur cgenes1 de type caratcter qui contient le numéro des clusters C1 et C2

?paste


## Créer un vectuer de charactères de taille 10 qui aura comme élements Rep1, Rep2, Rep3... Rep10
# s'aider le la fonction rep
# aide : créer un vecteur de taille 10 qui contient "Rep" 10 fois

?rep
## Utiliser Rep pour créer un vecteur avec 1 fois 1, 2 fois 2, 3 fois 3, jusquà 10. 1 2 2 3 3 3 ...


## Calculer la moyenne (mean) la mediane (median) du vecteur v suivant

??random
?runif
v <- runif(30,min= -50, max=100) 

