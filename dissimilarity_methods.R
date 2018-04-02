#' Packages needed to run C ++ code from R
#' You need to install the C ++ compiler (gcc) and language (g ++). If it is in linux, 
#' Can be installed from the package manager "synaptic"
#' In R, install.packages("Rcpp") and install.packages("inline")
library(Rcpp)
library(inline)
library(vegan)    # To calculate bray-curtis distance
#library(bioDist) # To calculate spearman distance. It loads when the sparman function is called
library(FD)       # To calculate mahanalobis distance. package ‘bioDist’ is not available (for R version 3.3.2)
library(proxy)    # To calculate cosine distance, 



# ......................... Standardized  euclidean distance (C++ Codes with Rcpp package) ............................
# C++ codes
src <- '
Rcpp::NumericMatrix Data(data);
int nrows = Data.nrow();
int ncolumns = Data.ncol();
Rcpp::NumericMatrix dissimilarity(nrows, nrows);

for (int i = 0; i < nrows; i++) 
{
  for (int j = 0; j < nrows; j++) 
  {
    for(int z = 0; z < ncolumns; ++z)
    {
      dissimilarity(i,j) = dissimilarity(i,j) + pow((Data(i,z) - Data(j,z)),2)/sd(Data(_,z));
    }
  }
}
return dissimilarity;
'
# R function
std_euclidean_distance <- cxxfunction(signature(data = "numeric"), body = src, plugin="Rcpp")



#' Function to calculate the matrix of dissimilarity in a set of individuals.
#' The variables in the data set must be of the numerical type or factors with numerical levels.
#' The diagonal of the matrix is zero, since by definition the distance of an individual with it is zero.
#' Distances: "manhattan", "euclidean", "std euclidean", "mahanalobis", "cosine", "correlation", "spearman", "chebyshev", "Canberra", "bray-Curtis"
#' param: mydata Data.frame containing the individuals on which the comparison was made
#' param: method Character string indicating distance measure for comparing individuals
#'
#' return: matrix contains the value of the distance between each pair of individuals of the original data set
#' Depends: vegan
#'
#' examples:
#' data <- iris[,1:4]
#' #Calculation of the euclidean distance between the first 10 individuals of the iris dataset
#' dissimilarity_iris <- dissimilar(mydata = data[1:10,], method = "euclidean")
#' dissimilarity_iris
dissimilar <- function(mydata = NULL, 
                       method = c("manhattan", "euclidean", "std euclidean", "mahanalobis", "cosine", "correlation", "spearman", "chebyshev", "Canberra", "bray-Curtis")) 
{
  if (is.null(mydata) == TRUE | class(mydata) != "data.frame")
  {
    stop("'mydata' parameter must contain a data.frame with the individuals to compare")  
  }
  
  if (method != "manhattan" && method != "euclidean" && method != "std euclidean" && method != "mahanalobis" && method != "cosine" && method != "correlation" && method != "spearman" && method != "chebyshev" && method != "Canberra" && method != "bray-Curtis") 
  {
    stop("No method found in this name.")
  }
  
  mydata <- as.matrix(mydata)
  colnames(mydata) <- NULL
  m <- nrow(mydata) # no. of rows
  n <- ncol(mydata) # no. of columns

  dissimilarity <- matrix(0,nrow = m, ncol = m)

  # 1 Manhattan Distance
  if (method == "manhattan")
  {
    dissimilarity <- as.matrix(dist(as.matrix(mydata), method = "manhattan", diag = TRUE, upper = TRUE))
  }

  # 2 Euclidean Distance
  if (method == "euclidean")
  {
    # 'dist' function contained in 'stats'package
    dissimilarity <- as.matrix(dist(as.matrix(mydata), method = "euclidean", diag = TRUE, upper = TRUE))
  }

  # 3 standardized Euclidean Distance with a diagonal variance matrix 
  if (method == "std euclidean")
  {
    dissimilarity <- std_euclidean_distance(as.matrix(mydata))
  }

  # 4 Mahanalobis Distance with a diagonal covariance matrix 
  if (method == "mahanalobis")
  {
    # 'mahaldist' function contained in 'FD' package
    dissimilarity <- mahaldis(as.matrix(mydata))
    dissimilarity <- as.matrix(dissimilarity)
  }

  # 5 Cosine Distance
  if (method == "cosine")
  {
    # proxy package used
    # 'dist' function are masked from ‘package:stats’ and 'as.matrix' s masked from ‘package:base’ 
    dissimilarity <- 1 - dist(as.matrix(mydata), method = "cosine", diag = TRUE, upper = TRUE)
    dissimilarity <- as.matrix(dissimilarity)
  }
  
  # 6 Correlation distance
  if (method == "correlation" )
  {
    # proxy package used
    # 'dist' function are masked from ‘package:stats’ and 'as.matrix' s masked from ‘package:base’ 
    dissimilarity <- dist(as.matrix(mydata), method = "correlation", diag = TRUE, upper = TRUE)
    dissimilarity <- as.matrix(dissimilarity)
  } 

  # 7 Spearman distance
  if (method == "spearman") 
  {
    # 'spearman.dist' function contained in 'biioDist' package
    library(bioDist)
    dissimilarity <- spearman.dist(as.matrix(mydata))
    dissimilarity <- as.matrix(dissimilarity)
    detach("package:bioDist", unload=TRUE)
  }  

  # 8 Chebyshev distance
  if (method == "chebyshev")
  {
    # proxy package used
    # 'dist' function are masked from ‘package:stats’ and 'as.matrix' s masked from ‘package:base’ 
    dissimilarity <- dist(as.matrix(mydata), method = "Chebyshev", diag = TRUE, upper = TRUE)
    dissimilarity <- as.matrix(dissimilarity)
  }            

  # 9 Canberra distance   
  if (method == "Canberra")
  {
    dissimilarity <-dist(mydata, method = "canberra", diag = TRUE, upper = TRUE)
    dissimilarity <- as.matrix(dissimilarity)
  }

  #10  Bray-Curtis distance 
  if (method == "bray-Curtis") 
  {
    A <- vegdist(mydata, method = "bray",  diag = TRUE, upper = TRUE)   
    dissimilarity <- as.matrix(A) 
  }             

  return(dissimilarity)
}
