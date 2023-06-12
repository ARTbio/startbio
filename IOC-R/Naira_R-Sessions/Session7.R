#Listes : https://www.tutorialspoint.com/r/r_lists.htm
 # Create a list containing strings, numbers, vectors and a logical
 # values.
 list_data <- list("Red", "Green", c(21,32,11), TRUE, 51.23, 119.1)
 print(list_data)
 list_data[1]
list_data[[1]]

# Create a list containing a vector, a matrix and a list.
list_data <- list(c("Jan","Feb","Mar"), matrix(c(3,9,5,1,-2,8), nrow = 2),
                  list("green",12.3))

# Give names to the elements in the list.
names(list_data) <- c("1st Quarter", "A_Matrix", "A Inner list")

# Show the list.
print(list_data)

names(list_data) <- c("a", "b", "c", "d", "e", "f")
print(list_data) 
LETTERS
names(list_data) <- LETTERS[1:6]
print(list_data)

# Create a list containing a vector, a matrix and a list.
list_data <- list(c("Jan","Feb","Mar"), matrix(c(3,9,5,1,-2,8), nrow = 2),
                  list("green",12.3))

# Give names to the elements in the list.
names(list_data) <- c("1st Quarter", "A_Matrix", "A Inner list")

# Access the first element of the list.
print(list_data[1])

# Access the thrid element. As it is also a list, all its elements will be print
ed.
print(list_data[3])

# Access the list element using the name of the element.
print(list_data$A_Matrix)


# Create a list containing a vector, a matrix and a list.
list_data <- list(c("Jan","Feb","Mar"), matrix(c(3,9,5,1,-2,8), nrow = 2),
                  list("green",12.3))

# Give names to the elements in the list.
names(list_data) <- c("1st Quarter", "A_Matrix", "A Inner list")

# Add element at the end of the list.
list_data[4] <- "New element"
print(list_data[4])

# Remove the last element.
list_data[4] <- NULL

# Print the 4th Element.
print(list_data[4])

# Update the 3rd Element.
list_data[3] <- "updated element"
print(list_data[3])

# Create two lists.
list1 <- list(1,2,3)
list2 <- list("Sun","Mon","Tue")

# Merge the two lists.
merged.list <- c(list1,list2)

# Print the merged list.
print(merged.list)

# créer 3 vecteurs v1, v2, et v3 de plusieurs tailes, contenant des nombres aléa
toires
v1 <- rnorm(3)
v2 <- rnorm(5)
v3 <- rnorm(7)
l1 <- list(v1,v2,v3)

print(l1)
names(l1) <- c("v1", "v2", "v3")
l1
mean(l1$v1)

lapply(l1, mean)
lapply(l1, length)
unlist(lapply(l1,length))


