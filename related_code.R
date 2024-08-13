
#install.packages("~/related_1.0.tar.gz", repos = NULL, type = "source")
#install.packages("related")
#library(related)
library(terra)
library(adegenet)
library(MASS)
library(FSA)
library(ggplot2)

################################### Calculating Relatedness ######################
setwd("~/Documents/Research/AMHA_Genetics/Analysis/final_pop_output_sex")

input<-readgenotypedata("related_genos3.txt")

##dyadml estimator
#coancestry estimates relatedness
#can add other estimators as arguments to compare
dyadmldf<-coancestry("related_genos3.txt", dyadml=2)

###returns spreadsheet for relatedness values for all possible pairs of individuals 

hist(dyadmldf$dyadml)
#ensure there are no duplicates included in the dataset (0.99 relatedness)

mean.relatedness<-mean(dyadmldf$dyadml)
sd.relatedness<-sd(dyadmldf$dyadml)
se.relatedness<-sd.relatedness/sqrt(length(dyadmldf$dyadml))

#look at whether relatedness between females is higher than relatedness between males
#f.f<-subset(dyadmldf, dyadmldf$ind1.id<500)
f.f<-subset(dyadmldf, dyadmldf$group=="F-F-")
mean.f.relatedness<-mean(f.f$dyadml)
sd.f.relatedness<-sd(f.f$dyadml)
se.f.relatedness<-sd.f.relatedness/sqrt(length(f.f$dyadml))


#m.m<-subset(dyadmldf, dyadmldf$ind1.id>400)
m.m<-subset(dyadmldf, dyadmldf$group=="M-M-")
mean.m.relatedness<-mean(m.m$dyadml)
sd.m.relatedness<-sd(m.m$dyadml)
se..m.relatedness<-sd.m.relatedness/sqrt(length(m.m$dyadml))

#f.m<-subset(dyadmldf, dyadmldf$ind1.id<500 & dyadmldf$ind2.id>500)
f.m<-subset(dyadmldf, dyadmldf$group=="F-M-" | dyadmldf$group=="M-F-")
mean.f.m.relatedness<-mean(f.m$dyadml)
sd.f.m.relatedness<-sd(f.m$dyadml)
se.f.m.relatedness<-sd.f.m.relatedness/sqrt(length(f.m$dyadml))

(mean_r<-c(mean.f.relatedness, mean.m.relatedness, mean.f.m.relatedness))
se<-c(se.f.relatedness, se..m.relatedness, se.f.m.relatedness)
(table<-rbind(mean_r, se))

#seems like two datasets are slightly different because 2nd round ended up with slightly higher estimates for every squirrel

########## Statistical Comparisons of Relatedness between Sex ##############

f.f$group<-"f.f"
f.m$group<-"f.m"
m.m$group<-"m.m"

new<-list(f.f, f.m, m.m)
df<-bind_rows(new)
df$group<-as.factor(df$group)

dyadmldf$group_new[dyadmldf$group=="M-F-"]<- "F-M-"
dyadmldf$group_new[is.na(dyadmldf$group_new)] <- dyadmldf$group[is.na(dyadmldf$group_new)]


#Since distributions are very skewed, probably better to run nonparametric test (but group values are not independent from each other?)
kruskal.test(dyadml~group_new, data=dyadmldf)

#Then follow up with Dunn post-hoc test to see which groups are significantly different
dunnTest(dyadml~group_new, data=dyadmldf, method="bonferroni")

boxplot(dyadml~group_new, data=dyadmldf)

############################# PCA ####################################
#setwd("~/Documents/Research/AMHA_Genetics/Analysis/final_pop_output_region/")
#Infile <- import2genind("populations.snps.gen")


Infile <- import2genind("populations.snps.gen")

x.infile <- tab(Infile, freq =TRUE, NA.method="mean")
pca.infile <- dudi.pca(x.infile, center = TRUE, scale = FALSE, nf=4)
#chose 3 axis, but should choose the number that covers the most amount of variation (where histogram plateaus)
summary(pca.infile)
#calculate percent of variation explained:
eig.perc <- 100*pca.infile$eig/sum(pca.infile$eig)
head(eig.perc)

s.label(pca.infile$li)
s.class(pca.infile$li, fac=pop(Infile), label = c(""), col = funky(15))

#89 = ID 533 (male)
#98 = ID 142 (female)

##results of PCA/axes defined by groups of SNPS/ composition of SNPS, 
#when tested .gen file associated with diff. regions - ours looks to be one population (no clear delineation into two or more groups) 
#also no clear delineation between sexes 

# Can map out SNPS on this PCA plot to see which SNPS are pulling the samples in which direction
#also helpful for visualizing variation within groups if some pops are very clustered and others are very spread out


#to assign labels to populations/regions
#48=3
#40=1
#98=2
#81=4
s.class(pca.infile$li, fac=pop(Infile), label = c("Four", "Two", "Three", "One"), col = funky(15))
#s.class(pca.infile$li, fac=pop(Infile), label = c("F", "M"), col = funky(15))
s.class(pca.infile$li, fac=pop(Infile), label = c("F","M"), col = funky(15))

#?s.class for more info about visualization

#################### Mantel Test ####################################


###when running mantel test....only expect results for females, not males.

setwd("~/Documents/Research/AMHA_Genetics/Analysis/final_pop_output_sex")
#setwd("~/Documents/Research/AMHA_Genetics/Analysis/final_pop_output_SRER")

data<-import2genind("populations.snps.gen")

#data<-import2genind("~/Documents/Research/AMHA_Genetics/Analysis/populations.snps_dispersal.gen")
#data_F<-import2genind("~/Documents/Research/AMHA_Genetics/Analysis/populations.snps_female.gen")
#data_M<-import2genind("~/Documents/Research/AMHA_Genetics/Analysis/populations.snps_male.gen")

#read in long/lat data for squirrel locations
long_lat <- read.csv("~/Documents/Research/AMHA_Genetics/Analysis/long_lat_dupsremoved_sex.csv")
#make sure this is in the same order as the .gen file !
#.gen file changes order based on whether it is divided by population

shared <- propShared(data)
Dps = 1 - shared
head(Dps)
head(data)
#calculate genetic distance
Dgen <- as.dist(Dps)
#calculate geographic distance (reference UTM columns for analysis to be done in meters)
Dgeo<- dist(long_lat[,14:15])
ibd <- mantel.randtest(Dgen,Dgeo,nrepet = 9999)
ibd


plot(ibd)
plot(Dgeo, Dgen)
dist_lm <- lm(as.vector(Dgen) ~ as.vector(Dgeo))
abline(dist_lm, col="red", lty=2)


dens <- kde2d(as.vector(Dgeo), as.vector(Dgen), n=300)
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(dist_lm)

##to look at different sexes, subset by population
#if using a .gen file was generated with a popmap that separated F/M (see populations.snps.gen file from final_output_sex folder), it will include a "pop" section for each 
#Adegenet will interpret the first population in the file as pop1, and the second as pop2
data_F <- data[pop=1, drop=TRUE]

#calculate proportion of shared alleles
shared_F <- propShared(data_F)

#calculate inverse
Dps_F = 1 - shared_F
head(Dps_F)
head(data)
#calculate genetic distance
Dgen_F <- as.dist(Dps_F)

#calculate geographic distance between individuals (m)
#need to exclude lat_long data for males
long_lat_F<-subset(long_lat, long_lat$sex=="F")
Dgeo_F <- dist(long_lat_F[,14:15])
ibd_f <- mantel.randtest(Dgen_F,Dgeo_F,nrepet = 9999)
ibd_f


data_M <- data[pop=2, drop=TRUE]
shared_M <- propShared(data_M)
Dps_M = 1 - shared_M
head(Dps_M)
head(data)
#calculate genetic distance
Dgen_M <- as.dist(Dps_M)
long_lat_M<-subset(long_lat, long_lat$sex=="M")
Dgeo_M <- dist(long_lat_M[,14:15])
ibd_m <- mantel.randtest(Dgen_M,Dgeo_M,nrepet = 9999)
ibd_m


plot(ibd_m)
#observation is right in the center where we would expect - nonsig 
plot(Dgeo_M, Dgen_M)
dist_lm <- lm(as.vector(Dgen_M) ~ as.vector(Dgeo_M))
abline(dist_lm, col="red", lty=2)


dens <- kde2d(as.vector(Dgeo_M), as.vector(Dgen_M), n=300)
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(Dgeo_M, Dgen_M, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(dist_lm)


plot(ibd_f)
#observation is far outside of what is expected-sig
plot(Dgeo_F, Dgen_F)
dist_lm <- lm(as.vector(Dgen_F) ~ as.vector(Dgeo_F))
abline(dist_lm, col="red", lty=2)


dens <- kde2d(as.vector(Dgeo_F), as.vector(Dgen_F), n=300)
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(Dgeo_F, Dgen_F, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(dist_lm)
#have a few clusters of samples, mantel test is best when you have a continuous distribution
#if you have two subpopulations that are geographically separated and not sampled continuously between - 
#can get patterns like this, so need to look at substructures more to see if they are clustering by relatedness or not. 
#can run addmixture test-- shows how many clusters of ancestry is supported



