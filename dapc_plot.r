# Discriminant Analysis of Principal Components (DAPC)
# Nice video that explains PCA ==> https://www.youtube.com/watch?v=_UVHneBUBW0
# It was designed to identify and describe clusters of genetically related individuals or describe genetic clusters using a few variables.
# Visual assessment of between-population differentiation.
# DAPC identifies the best number of clusters for the data increasing intergroup variance and reducing intra-group variance. 
# Reference: Jombart T, Devillard S and Balloux, F (2010). Discriminant analysis of principal components: a new method for the analysis of genetically structured populations. BMC Genetics 11: 94.
# DAPC tutorial: https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf

######################################## Load Libraries ########################################

library(adegenet)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

######################################## Set working directory ##################################

setwd("~/Documents/Research/AMHA_Genetics/Analysis/final_pop_output_sex")


#getwd()

######################################## Input data #############################################

#import file as genind object (use .gen file from stacks=>populations output)
infile <- import2genind("populations.snps.gen")


#import population information and add column names. In pops-for-pca.txt, first columns is popualtion and second column is sample name
pops <- read.table("~/Documents/Research/AMHA_Genetics/Analysis/squirrel_popmap_dupsremoved_sex.txt", header=F, stringsAsFactors = TRUE)
colnames(pops) <- c('individual', 'population')

#If you look at infile@pop, it doesn't have the correct populations assigned to each sample
infile@pop

#Define population map by assigning population names from the first column of the table we imported ("pops-for-pca.txt"), which is called pops
infile@pop <- pops$population

#check names of populations again
infile@pop

#check genind file
infile

#' GENIND OBJECT /////////
#'   
#'   // 47 individuals; 21,958 loci; 43,916 alleles; size: 19.4 Mb
#' 
#' // Basic content
#' @tab:  47 x 43978 matrix of allele counts
#' @loc.n.all: number of alleles per locus (range: 2-2)
#' @loc.fac: locus factor for the 43978 columns of @tab
#' @all.names: list of allele names for each locus
#' @ploidy: ploidy of each individual  (range: 2-2)
#' @type:  codom
#' @call: read.genepop(file = file, quiet = quiet)
#' 
#' // Optional content
#' @pop: population of each individual (group size range: 47-47)
#' > 

#Check percentage of complete genotypes per locus
# @loc.fac is a factor indicating which columns in @tab correspond to which marker
# propTyped investigates the structure of missing data in adegenet objects.
# propTyped returns the proportion of available (i.e. non-missing) data per individual/population, locus, or the combination of both in with case the matrix indicates which entity (individual or population) was typed on which locus.
# in propTyped the argument "by" could be "ind","loc", or "both"

infile@loc.fac
locmiss_infile = propTyped(infile, by = "loc")

(length(locmiss_infile))
#21958  ### which is the number of loci

#save the percentage_of_complete_genotypes_per_loci as a csv file in your working directory
write.csv(locmiss_infile, file = "percentage_of_complete_genotypes_per_loci.csv")
#we had a cutoff of -.25 for maf, so all of our percentages are above 75%

########################################################################################################################
#################################### Identify clusters with find.clusters  ###############################################
########################################################################################################################

####### Selecting number of clusters with find.clusters ###########
# Info about find.clusters => https://search.r-project.org/CRAN/refmans/adegenet/html/find.clusters.html

# DAPC requires prior groups to be defined.

# These functions implement the clustering procedure used in Discriminant Analysis of Principal Components (DAPC, Jombart et al. 2010). 
# This procedure consists in running successive K-means with an increasing number of clusters (k), 
# after transforming data using a principal component analysis (PCA). For each model, 
# a statistical measure of goodness of fit (by default, BIC) is computed, which allows to choose the optimal k.

# find.clusters requires selecting # of PCs and # of clusters to be retained. 

# find.clusters will display a graph of cumulative variance. No reason for keeping a small number of components.
# When asked, choose the number PCs to retain (>= 1): respond your number of samples-1   <========= !!!!!!!!!!

# Then, a graph showing BIC (Bayesian Information Criterion) values for increasing values of k will show.
# To identify the optimal number of clusters, k-means is run sequentially with increasing values of k, 
# and different clustering solutions are compared using Bayesian Information Criterion (BIC).
# Ideally, the optimal clustering solution should correspond to the lowest BIC. <========= !!!!!!!!!!

# specify the max number of groups or clusters to evaluate in max.n.clust
k <- 5

clusters <- find.clusters(infile, clust = NULL, n.pca = NULL, n.clust = NULL, method = "kmeans", stat = "BIC", 
                         choose.n.clust = TRUE, criterion = "diffNgroup", 
                         max.n.clust = k, n.iter = 1e5, n.start = 10, scale = FALSE, truenames = TRUE)
#PCs to retain = 46 (n-1)
#will show basically straight positive line, meaning that it is all one pop- lowest number (1) is most likely. 
#But won't let you choose just 1 cluster to plot, so just for plotting purposes choose 2 clusters to retain
clusters

# now we will check how well groups have been retrieved by the procedure. In the graph that will be generated, 
# rows correspond to original grouping (”ori”), and columns correspond to inferred groups (”inf”).

table(pop(infile), clusters$grp)

table.value(table(pop(infile), clusters$grp), col.lab=paste("inf", 1:9),
            row.lab=paste("ori", 1:14))
#if your snps.gen file / pop map were done with files separated by region/population, 
#then it will show your original grouping vs. the inferred grouping
#If you input files as all one population, the original grouping won't be included in the plot
#remember that we forced this into 2, so not necessarily 2 populations just because it grouped some individuals into 2nd group. 

#export the "original_vs_inferred_membership" manually


########## Describing clusters using dapc ############

# ”How many clusters are useful to describe the data?”
# IMPORTANT NOTE: There is no ”true k”, but some values of k are better or more efficient summaries of the data than others.
# clusters are constructed as linear combinations of the original variables (alleles),
# which have the largest between-group variance and the smallest within-group variance.
# dapc provides membership probabilities of each individual for the different groups based on the retained discriminant functions. 

# Run the analysis on the dataset, using the inferred groups stored in clusters$grp

# When the graph of cumulative variance is shown, pick PCs to retain BUT now minimize # of PCs without sacrificing too much info.
# Retaining too many components with respect to the number of individuals can lead to over-fitting. 
# The number of principal components retained has a large effect on the outcome of the data. Here n.pca = NULL (default)
# Notice when the curve plateaus meaning little information is gained by adding PCs after that number.  <========= !!!!!!!!!!
# See the section below for a statistical method called cross- validation as an aid for choosing n.pca (below)

dapc1 = dapc(infile, clusters$grp, n.pca=20)
#n.pca argument derived from cross-validation below 
#can also ommit this argument and choose 

# Then, a barplot of eigenvalues will be displayed, asking for a number of discriminant functions to retain (unless argument n.da is provided).
# For small number of clusters, all eigenvalues can be retained since all discriminant functions can be examined 
# without difficulty. Whenever more (say, tens of) clusters are analysed, 
# it is likely that the first few dimensions will carry more information than the others, 
# and only those can then be retained and interpreted.

# check information in dapc1. For details about the content,read the documentation (?dapc).

### for my project - ideal PCA=20 (based on cross validation test below)
#only shows one discriminate function for bar plot - so just choose 1. 


dapc1

# Basically,ind.coord and grp.coord contain the coordinates of the individuals and of the groups used in scatterplots.
# Contributions of the alleles to each discriminant function are stored in var.contr.

# Check individual assignments: 
# Numerical group assignment can change each time you rerun this,
# so it's good to check that the order of group labels used below for plotting is still correct.

dapc1$grp

# make a basic scatterplot
plot1 = scatter.dapc(dapc1)
plot2 = scatter.dapc(dapc1, label.inds =list(air=0.5, cex=0.05, pch=NA), clab=0, cstar=0, cell=1, legend = TRUE, posi.leg = "topleft", cleg = 0.5 )

#view details of scatter.plot for help
?scatter.dapc

# To look at how each individual assigns to the clusters, we can run commonplot.
# This is a good way to understand how strong the structuring is in your population.
# If individuals assign with high probability to their respective clusters, you know you have arrived at a reasonable solution.

nb.cols <- 14

mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)

compoplot.dapc(dapc1, lab="", txt.leg=c("Cluster1","Cluster2"), n.col=2, cleg=0.8, xlab="individuals", cex.lab=1, col.pal=mycolors)

# Analyse how much percent of genetic variance is explained by each axis
percent = dapc1$eig/sum(dapc1$eig)*100

#Export the probability of membership for each individual
write.table(dapc1$posterior, "DAPC_membership_probability.txt", quote = FALSE)

#check which are the most ’admixed’ individuals.
#for my output, most have 100% membership probability into their clusters
#Let us consider admixed individuals as having no more than 90% of probability of membership in a single cluster

admix_ind <- which(apply(dapc1$posterior,1, function(e) all(e<0.90)))
admix_ind 
#all of ours are above 90%, no admixed individuals.

#########  dapc plot with ggplot #########
#I think the following code is just making scatter.dapc plot prettier
# ggplot cheat sheet ==> https://www.maths.usyd.edu.au/u/UG/SM/STAT3022/r/current/Misc/data-visualization-2.1.pdf

#Create a data.frame containing individual coordinates
ind_coords = as.data.frame(dapc1$ind.coord)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(infile)

# Add a column with the site IDs
ind_coords$Site = infile$pop

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

# Define color palette
cols = mycolors

# Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15) 
)

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = TRUE)+
  
  # centroids ==> remove next part if you don't want Site (population) labels on the plot
  #geom_label(data = centroid, aes(label = Site, fill = Site), size = 2, show.legend = FALSE)+
  
  # coloring
  scale_fill_manual(values = mycolors)+
  scale_colour_manual(values = mycolors)+
  
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("DAPC plot")+
  # custom theme
  ggtheme

# Export plot #could change extension to .png
ggsave("infile_dpca_plot.pdf", width = 12, height = 8, dpi = 600)


########################################################################################################################
##################### Cross validation to find the optimal number of PCs to retain in DAPC #############################
########################################################################################################################

# DAPC tutorial: https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
# Cross-validation for DAPC info: https://rdrr.io/cran/adegenet/man/xvalDapc.html
#a better way to find ideal number of PCs to retain, instead of just looking at it and guestimating.
#rerun DAPC code above based on the results of this
# xvalDapc uses a training and testing set to find the optimal trade-off between selecting too few, and too many PC axes
# NA.method, a character indicating if missing values should be replaced by the mean frequency ("mean"), or left as is ("asis").

set.seed(124)
z = tab(infile, NA.method = "mean")

# xvalDapc_plot (might take a few minutes).
# When xval.plot is TRUE, a scatterplot of the DAPC cross-validation is generated. 
# Export cross_val_plot manually
crossval = xvalDapc(z, infile$pop,  n.pca.max = 300, n.rep = 30, training.set = 0.9, result = "groupMean", xval.plot = TRUE)

# The number of PCs retained in each DAPC varies along the x-axis, and the proportion of successful outcome prediction varies along the y-axis. 
# Individual replicates appear as points, and the density of those points in different regions of the plot is displayed in blue.
# As one might expect (or hope) for an optimization procedure, the results of cross-validation take on an arc-like shape.
#  Predictive success is sub-optimal with both too few and too many retained PCA axes.
crossval

# Number of PCs with best stats (lower score = better)
crossval$`Root Mean Squared Error by Number of PCs of PCA`

#2         4         6         8        10        12        14        16        18        20 
#0.5948522 0.6443166 0.5147057 0.4838452 0.5581231 0.5767204 0.5800016 0.5617855 0.5381241 0.5677540 

crossval$`Number of PCs Achieving Highest Mean Success`
#[1] "20"

crossval$`Number of PCs Achieving Lowest MSE`
#[1] "20"

numPCs = as.numeric(crossval$`Number of PCs Achieving Lowest MSE`)
numPCs
#[1] 8

