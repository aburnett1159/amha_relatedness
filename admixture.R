#install dependencies
#install.packages(c("ggplot2","gridExtra","label.switching","tidyr","remotes"),repos="https://cloud.r-project.org")
# install pophelper package from GitHub
#remotes::install_github('royfrancis/pophelper')
# Load library
library(pophelper)

setwd("~/Documents/Research/AMHA_Genetics/Analysis/Admixture_May2024_SRER")

#### Analyse cross validation to find optimal number of clusters 

#read cross validation .txt file and make a plot. Save the plot as in your current diretory
#note that cross validation table might have k=10 as the first row! 
#So might need to go in and change that so that the following code actually puts K=1 at K=1. Otherwise will have K=10 at K=1 position.
cv.table <- read.table("cross_validation.txt") #contains cross validation error for each value of K tested
rows.table<- c(1:10) # create a list of number 1 to 10 to add to cv table 
cv.dataframe <- data.frame(rows.table, cv.table[,4]) # create cv data frame with 2 columns only 
plot(cv.dataframe[,1], cv.dataframe[,2], type = "o", main = "Admixture cross-validation error", xlab = "K", ylab = "Cross-validation error")
##look at the lowest number to see likely pop number and/or ancestries -- mixed population K refers more to ancestry 

# Create multi-line barplots from qlist 
# To format plots see http://www.royfrancis.com/pophelper/articles/index.html

# Set population and individual labels
pops <- read.table("pops-for-admixture.txt", stringsAsFactors=FALSE)
inds <- read.delim("inds-for-admixture.txt", header=FALSE, stringsAsFactors=FALSE)
popnames <- c("SRER")

# Admixture plot for K=1 with NO COLORS. It will be saved in working directory as .pdf file
admix_1Q <- readQ(files = "recode_plink.1.Q")
rownames(admix_1Q[[1]]) <- inds$V1

# add individual labels
colnames(pops) <- " "    #Prevents group label from showing up in bottom right corner of graph. 
plotQ(admix_1Q, splab = "K=1", grplab = pops,  splabsize=5, legendtextsize=5, 
      grplabangle=90, grplabsize=2, grplabspacer=0, grplabjust=0, grplabpos = 0, ordergrp = TRUE, 
      useindlab = T, showindlab = T, indlabsize=3, indlabcol = "black", linepos = 0.9, barbordercolour="white",barbordersize=0.05, 
      subsetgrp = popnames,
      exportpath=getwd(), imgtype="pdf", width = 20, height = 6, dpi = 600,
      showlegend = TRUE)




# Admixture plot for K=2 with NO COLORS. It will be saved in working directory as .pdf file
admix_2Q <- readQ(files = "recode_plink.2.Q")
rownames(admix_2Q[[1]]) <- inds$V1

# add individual labels
colnames(pops) <- " "    #Prevents group label from showing up in bottom right corner of graph. 
plotQ(admix_2Q, splab = "K=2", grplab = pops,  splabsize=5, legendtextsize=5, 
      grplabangle=90, grplabsize=2, grplabspacer=0, grplabjust=0, grplabpos = 0, ordergrp = TRUE, 
      useindlab = T, showindlab = T, indlabsize=3, indlabcol = "black", linepos = 0.9, barbordercolour="white",barbordersize=0.05, 
      subsetgrp = popnames,
      exportpath=getwd(), imgtype="pdf", width = 20, height = 6, dpi = 600,
      showlegend = TRUE)


#no clear genetic substructure in my population, pretty admixed 
#so since K=2, but no obvious grouping in substructure (would see blue bars grouped separately from purple bars if really two pops),
#what we're probably seeing is 2 ancestries coming together and mixing together in this location. 

#If you want to see how males & females compare, can add f_ and m_ to ind ID's in "inds-for-admixture.txt" file and rerun 


# K2 and K3 and K4, with colors, indiv. labels and export as pdf (number of colors can be the same number as Ks 

standardcol <- c("#2121D9","#9999FF","#DF0101","#04B404","#FFFB23","#FF9326","#A945FF","#0089B2","#B26314","#610B5E")

#add individual labels to all .Q files to be included

admix_2_3_4Q <- readQ(files = c("recode_plink.2.Q", "recode_plink.3.Q", "recode_plink.4.Q" ))
if(length(unique(sapply(admix_2_3_4Q, nrow)))==1) admix_2_3_4Q <- lapply(admix_2_3_4Q,"rownames<-",inds$V1)
lapply(admix_2_3_4Q,rownames)

align2.3.4 <- alignK(admix_2_3_4Q, type = "across")   #Makes sure clusters are numbered consistently for each value of K.

colnames(pops) <- " "
plotQ(align2.3.4, splab = c("K=2", "K=3", "K=4"),  grplab = pops, splabsize=5, legendtextsize=5,
      grplabangle=90, grplabsize=1.5, grplabspacer=-0.2, grplabjust=0, grplabpos = 0, ordergrp = TRUE,
      useindlab = T, showindlab = T, indlabsize=3, indlabcol = "gray", linepos = 0.9, barbordercolour="white",barbordersize=0.05, 
      subsetgrp = popnames,
      clustercol = standardcol,
      exportpath=getwd(), imgtype="pdf", width = 20, height = 3.2, dpi = 600,
      showlegend = TRUE, imgoutput="join")
