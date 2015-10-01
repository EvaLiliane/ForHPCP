#!/bin/bash
R CMD BATCH --no-save --no-restore '--args "/home/evaliliane/Documents/PhD/Codes/NewData/CD4Cat_Children2015-04-20.csv" "ChildrenIndFitData.csv" 5 5' MakeIndFits.R IndOutChildren.out
R CMD BATCH --no-save --no-restore '--args "/home/evaliliane/Documents/PhD/Codes/NewData/CD4Cat_Adults2015-04-20.csv" "AdultsIndFitData.csv" 5 5' MakeIndFits.R IndOutAdults.out

# To run this file use:
# To make it executatble:  chmod +x MasterFile.sh
# To actuall run it: ./MasterFile.sh

#R CMD BATCH --no-save --no-restore '--args "/input/CD4Cat_Children2015-04-20.csv" "ChildrenIndFitData.csv" 5 5' MakeIndFits.R IndOutChildren.out
#R CMD BATCH --no-save --no-restore '--args "/input/CD4Cat_Adults2015-04-20.csv" "AdultsIndFitData.csv" 5 5' MakeIndFits.R IndOutAdults.out
