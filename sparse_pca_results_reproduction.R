#######
#Elasticnet
#######

#install.packages("elasticnet")
library(elasticnet)
library(FactoMineR)
library(corrplot)
library(ggplot2)
rm(list=objects()); graphics.off()

######
##Variance à la main selon la section 3.5 
######

explained_variance <- function(U, X) {
  U <- U / sqrt(sum(U^2))  # Normalize U
  percent_var <- 100 * adjusted_variance(U, X) / sum(X^2)
  return(percent_var)
}

adjusted_variance <- function(U, X) {
  pcs <- t(U) %*% X
  qr_decomp <- qr(pcs)
  R2 <- qr.R(qr_decomp)^2
  v <- sum(diag(R2))
  return(v)
}
##On peut utiliser var.all et pev déjà intégrés à elasticnet

########
#Data exploration
########

# Charger la bibliothèque pour les opérations sur les matrices (par exemple : la trace d'une matrice).
library(psych)

# Fixer la graine pour des résultats reproductibles.
set.seed(2)

# Importer le jeu de données dans R (jeu de données normalisé ou matrice de corrélation).

data("pitprops")
Data <- pitprops

#Explication de l'article original de Jeffers (1967):

#TOPDIAM - le diamètre supérieur du support en pouces ;
#LENGTH - la longueur du support en pouces ;
#MOIST - la teneur en humidité du support, exprimée en pourcentage du poids sec ;
#TESTSG - la densité spécifique du bois au moment du test ;
#OVENSG - la densité spécifique du bois après séchage au four ;
#RINOTOP - le nombre d'anneaux annuels au sommet du support ;
#RINGBUT - le nombre d'anneaux annuels à la base du support ;
#BowMAx - l'arc maximum en pouces ;
#BOWDIST - la distance du point d'arc maximum depuis le sommet du support en pouces ;
#WHORLS - le nombre de tours de nœuds ;
#CLEAR - la longueur du support dépourvue de nœuds depuis le sommet en pouces ;
#KNOTS - le nombre moyen de nœuds par tour de nœuds ;
#DIAKNOT - le diamètre moyen des nœuds en pouces.

# Nombre de variables dans Data.
nvar <- ncol(Data)
corrplot(Data) #Les variables sont assez corrélées.









########
##Reproduction des résultats
########


# 1) Reproduction de la figure 1 pour des données générées (Soft thresholding)
x <- seq(-5, 5, length.out = 100)

# Definir soft-thresholding 
soft_thresholding <- function(x, delta) {
  return(pmax(abs(x) - delta, 0) * sign(x))
}
# Set delta value
delta <- 1
y <- soft_thresholding(x, delta)
plot(x, y, type = "l", col = "black", lwd = 3, xlab = "x", ylab = "y",
     main = "Soft-thresholding Rule with Delta=1",
     asp = 1, xlim = c(-2, 2), ylim = c(-2, 2))
abline(h = 0, col = "black", lty = 1, lwd = 1)  # Axes
abline(v = 0, col = "black", lty = 1, lwd = 1)
abline(a = 0, b = 1, col = "green", lty = 2, lwd = 1)  # bissectrice
legend("topleft", legend = expression(y == (abs(x) - delta) + Sign(x)),
       col = "black", lty = 1, lwd = 3, bty = "n")


########
#Reproduction des tableaux 1 et 3 à l'aide d'elasticnet
########

out<-spca(Data,K=6,type="Gram",sparse="penalty",trace=TRUE,para=c(0,0,0,0,0,0)) #pour la table 1, on a une simple SVD

out3<-spca(Data,K=6,type="Gram",sparse="penalty",trace=TRUE,para=c(0.06,0.16,0.1,0.5,0.5,0.5))

out2<-spca(Data,K=6,type="Gram",sparse="varnum",trace=TRUE,para=c(7,4,4,1,1,1))


### TABLE 1: PCA normale (voir plus bas pour afficher une jolie table)
out$loadings #les loadings
out$pev #la variance
# la variance cumulée
cumulative_variance <- cumsum(out$pev) 
print(cumulative_variance)
#### Table 2: S
##table out2
names(out2)
##  loadings
out2$loadings
out2$prop.var.explained
out2$var.all #Variance totale
out2$pev #Explained variance
summary(out2)
#### Table 3 et Figure 2: PCA SPARSE SUR PITPROPS
count_non_zero <- function(x) {
  if (is.matrix(x)) {
    # Count non-zero elements each column
    counts <- apply(x, 2, function(col) sum(col != 0))
  } else if (is.vector(x)) {
    # Count non-zero elements vector
    counts <- sum(x != 0)
  } else {
    stop("Input must be a vector or a matrix")
  }
  
  return(counts)
}

##  loadings
out3$loadings
non_zero<-count_non_zero(out1$loadings)
non_zero #la contrainte de sparsity Non zeros
out3$pev #Explained variance
# la variance cumulée
cumulative_variance3 <- cumsum(out3$pev) 
print(cumulative_variance3)
### figure 2: Pev=f(lamda)
# Set up lambda grid
lambda.grid <- seq(0, 3.5, 0.01)

par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))  # 3x2 layout for subplots

# Loop through principal components
for (j in 1:6) {
  # on stocke dans une matrice les variances ajustées
  PEV <- matrix(0, ncol = 6, nrow = length(lambda.grid))
  
  # on fait une boucle
  for (i in seq_along(lambda.grid)) {
    cat(i, "\n")
    out <- spca(pitprops, K = 1, type = "Gram", sparse = "penalty", lambda = 0,
                 para = rep(lambda.grid[i], 6), eps.conv = 1e-4)
    PEV[i, ] <- out$pev[]
  }
  
  # Plot
  plot(plot(lambda.grid, PEVmatrix[,1], type = "l", xlim = c(0, 0.3)) #Pas le même plot: ça décroit trop rapidement;

######
##Exemple synthétique (Table 5)
######


# Number of observations
n <- 10

# Hidden
V1 <- rnorm(n, mean = 0, sd = sqrt(290))
V2 <- rnorm(n, mean = 0, sd = sqrt(300))
V3 <- -0.3 * V1 + 0.925 * V2 + rnorm(n, mean = 0, sd = 1)

# Observable 
X <- data.frame(
  X1 = V1 + rnorm(n, mean = 0, sd = 1),
  X2 = V1 + rnorm(n, mean = 0, sd = 1),
  X3 = V1 + rnorm(n, mean = 0, sd = 1),
  X4 = V1 + rnorm(n, mean = 0, sd = 1),
  X5 = V2 + rnorm(n, mean = 0, sd = 1),
  X6 = V2 + rnorm(n, mean = 0, sd = 1),
  X7 = V2 + rnorm(n, mean = 0, sd = 1),
  X8 = V2 + rnorm(n, mean = 0, sd = 1),
  X9 = V3 + rnorm(n, mean = 0, sd = 1),
  X10 = V3 + rnorm(n, mean = 0, sd = 1)
)


#En prenant une PCA classique, les valeurs ne sont pas tout à fait les mêmes
PCA <- prcomp(X, scale. = TRUE)
PCA$rotation[, 1:3] #comme dans l'article
PCA$sdev[1:3]
PCA$sdev[1:3]/sum(PCA$sdev[1:3])

out5 <- spca(var(X), K = 10, para = rep(5, 10), sparse = "penalty",
            max.iter = 100, lambda = 0,
            use.corr = TRUE)

out5$pev
#### Faire une jolie table


# Load necessary libraries
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}
library(gridExtra)
if (!requireNamespace("grid", quietly = TRUE)) {
  install.packages("grid")
}
library(grid)

# Assuming 'out' is the PCA result
# Replace 'out' with your actual PCA result object

# Load necessary libraries
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}
library(gridExtra)

# Create a data frame with loadings, variance, and cumulative variance
table_df <- data.frame(
  X = c("topdiam", "length", "moist", "testsg", "ovensg", "ringtop", "ringbut", "bowmax", "bowdist", "whorls", "clear", "knots", "diaknot", "Variance (%)", "Cumulative variance (%)"),
  PC1 = round(c(out$loadings[, 1], out$pev[1], cumsum(out$pev)[1]), 2),
  PC2 = round(c(out$loadings[, 2], out$pev[2], cumsum(out$pev)[2]), 2),
  PC3 = round(c(out$loadings[, 3], out$pev[3], cumsum(out$pev)[3]), 2),
  PC4 = round(c(out$loadings[, 4], out$pev[4], cumsum(out$pev)[4]), 2),
  PC5 = round(c(out$loadings[, 5], out$pev[5], cumsum(out$pev)[5]), 2),
  PC6 = round(c(out$loadings[, 6], out$pev[6], cumsum(out$pev)[6]), 2)
)

# Create a tableGrob 
table_grob <- tableGrob(table_df)

# Display the table image
grid.newpage()
grid.draw(table_grob)





##########
#PMA et sparse PCA
##########
#Le code est du package sparsepca est très clair et disponible à l'adresse https://github.com/erichson/spca/blob/master/R/spca.R

