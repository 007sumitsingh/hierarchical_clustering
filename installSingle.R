installSingle<-function(nameOfPackage)
{
  if(nameOfPackage %in% rownames(installed.packages()) == FALSE) 
  {
    install.packages(nameOfPackage)
  }
  library(nameOfPackage,character.only = TRUE)
}