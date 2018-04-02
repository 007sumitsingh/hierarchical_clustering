setup <- function(){
  Sys.setenv(LANGUAGE="en")
  Sys.setlocale(category = "LC_ALL", locale = "English")
  source('installSingle.R')
  
  installSingle('Rcpp')
  installSingle('FD')
  installSingle('RcppArmadillo')
  installSingle('inline')
  installSingle('vegan')
  installSingle('mgcv')
  installSingle('proxy')
  installSingle('pdfCluster')
  installSingle('clue')


  availableFiles <- list.files(pattern="[.]R$", path="./", full.names=TRUE);
  availableFiles <- availableFiles[ !grepl("./setup.R", availableFiles) ];
  x<-sapply(availableFiles, source,echo=FALSE);
}