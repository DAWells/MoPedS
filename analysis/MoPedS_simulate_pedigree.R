#Simulate a pedigree based on the pdata with user specified coefficients for the variables
#Simulate genotypes based on that pedigree
#Fit the pedigree?

setwd("/Users/davidwells/Dropbox/Mongoose/Pedigree/MoPedS/data")
packages<-c("MasterBayes", "dplyr")
lapply(packages,require,character.only=T)

#load phenotypic data created in MoPedS_2.R
pdata<-read.csv("phenotypic_data_oct_2018.csv")

#Error rates
error<-read.csv('rawdata/error rates.csv')
#When the error rate is 0 change it to 0.005 as it is unlikely to be absolutely 0
error$E1[error$E1==0]<-0.005
error$E2[error$E2==0]<-0.005

#microsatellite data
gendata<-read.csv("rawdata/BMRP_Genotype_data_July_2016.csv")

#Use NA for missing data, which are currently 0s
gendata[gendata==0]<-NA

#Calculate the proportion of individuals which are blank at a given locus. This is used for simulating realistic genotypes later.
#High missing proportions are because we stopped genotyping certain loci.
fails<-apply(gendata[,2:ncol(gendata)],2,function(col) mean(is.na(col)))
fails<-matrix(fails,byrow=T,ncol=2)
fails<-apply(fails,1,mean)

#Add the ungenotyped individuals to gendata but with NA for all loci
ungen<-unique(pdata$id[!pdata$id%in%gendata$id])
ungen<-data.frame(id=ungen)

gendata<-merge(gendata,ungen,by="id",all=T)

#Alphabetically reorder the rows of error and columns of gendata so that they match.
gendata<-gendata[,c(1,(1+order(names(gendata)[2:87])))]
error<-error[order(error$locus),]

print("Please check that the names of loci in gendata and error correctly match up")
cbind(names(gendata)[2:87], as.character(rep(error$locus,each=2)))


###############################
#                             #
#      Simulate pedigree      #
#                             #
###############################

#Set offspring to 0 for all those you do not want to attempt to assign parents to.
#This should be done in MoPeds_2
#pdata[format.Date(pdata$date,"%Y")!=2015,]$offspring<-0
#Only look at 2015 as examples for brevity.
#pdata<-pdata[format.Date(pdata$date,"%Y") %in% c(2015,2014),]

!!!!! only Male and Female should be sexes. P must therefore be changed to NA.
#levels(pdata$sex)[c(2,4)]<-NA
levels(pdata$sex)[3]<-NA

offspring<-filter(pdata,offspring==1)

####################
#Restrictions
####################
#Dam restrictions.
#Same pack as offspring at birth.
dres1<-expression(varPed(x="pack", gender="Female", relational="OFFSPRING", restrict="=="))
#Over 6 months old.
dres2<-expression(varPed(x="juv", gender="Female", restrict=F))
#Offspring cannot be its own dam.
dres3<-expression(varPed(x="id", gender="Female", relational="OFFSPRING", restrict="!="))
#Dam must be alive when the pup was born.
dres4<-expression(varPed(x="datesalive", gender="Female", relational="OFFSPRING", restrict="=="))

#Sire restrictions.
#Sires must be over six months old.
sres2<-expression(varPed(x="juv", gender="Male", restrict=F, lag=c(-60,-60)))
#Offspring cannot be its own sire.
sres3<-expression(varPed(x="id", gender="Male", relational="OFFSPRING", restrict="!=", lag=c(-60,-60)))
#Sire must be alive 2 months before the pup was born, this is approximately the gestation period.
sres4<-expression(varPed(x="datesalive", gender="Male", relational="OFFSPRING", restrict="==60", lag=c(-60,-60)))


####################
#    Predictor variables
####################

#has the potential dam given birth in the month the individual was born in or the month on either side?
var1<-expression(varPed(x="givenbirth", gender="Female", lag=c(-30,30)))

#Was the potential sire in the pack that the pup was born into at any point in the three months leading up to birth?
var2<-expression(varPed(x="pack", gender="Male", relational="OFFSPRING", lag=c(-90,0)))

#All age and age squared to have separate effects for each sex
var3<-expression(varPed(x='scaled_age', gender="Female"))
var4<-expression(varPed(x='scaled_age', gender="Male", lag=c(-60,-60)))
var5<-expression(varPed(x="sq_scaled_age", gender="Female"))
var6<-expression(varPed(x="sq_scaled_age", gender="Male", lag=c(-60,-60)))

PdP<-PdataPed(list(dres1,dres2,dres3,dres4, sres2,sres3,sres4, var1, var2, var3, var4, var5, var6), data=pdata, USdam=T, USsire=T, timevar=pdata$datesalive)
#saveRDS(PdP,"PdataPed_object.RDS")

ped_beta<-c(6,4,0.4,0.9,-0.02,-0.03)
simped<-simpedigree(PdP, beta=ped_beta, nUS=c(10,10))

#To include imigrants, you must edit the pedigree
!THIS is trial code to check the scripts work. Do not use in actual analyses!!!

#immigrnats
I <-filter(PdP$data,new_immigrant==T)$id

#Make several immigrants full siblings
simped$ped[simped$ped[,1]%in%I[1:4],2]<-"Fakedam1"
simped$ped[simped$ped[,1]%in%I[1:4],3]<-"Fakesire1"

#Make some immigrants full sibs with real offspring.
#Must check the simped and find an offspring with two assigned parents and use those parents here.
simped$ped[simped$ped[,1]%in%I[5:6],2]<-"BF484"
simped$ped[simped$ped[,1]%in%I[5:6],3]<-"BM500"

#Make some immigrants half sibs with real offspring, they should not be assigned
simped$ped[simped$ped[,1]%in%I[7],2]<-"DF208"
simped$ped[simped$ped[,1]%in%I[8],3]<-"HM216"


#saveRDS(simped,"true_simulated_pedigree.RDS")
################################
#                              #
#      Simulate genotypes      #
#                              #
################################

#When simulating genotypes do remember that 8 loci are no longer genotyped. FS50, MON28, Mon30, AHT130, Ag8, hic1.95, fs41,ss7-1, fs46, fs48.

#Simulate using the complete set of loci.
simgen<-simgenotypes(A=extractA(gendata), E1=error$E1, E2=error$E2, ped=simped$ped, prop.missing=0)

#Simulate using only the current set of 35 loci.
#current_gen<-!names(gendata) %in% c("fs46", "fs46.1","fs48","fs48.1","fs50","fs50.1","AHT130a","AHT130b","Ag8a","Ag8b","Ss7.1a","ss7.1b","fs41", "fs41.1", "hic.1.95", "hic.1.95.1")
#error<-error[current_gen[2:87][rep(c(T,F),43)],]
#simgen<-simgenotypes(A=extractA(gendata[,current_gen]), E1=error$E1, E2=error$E2, ped=simped$ped, prop.missing=0)


#saveRDS(simgen, "simulated_genotypes.RDS")


#Save the pedigree and genotype data so that it is easy to fit into mopeds 3 and 4.
#Convert the genotype to two column format
#include columns for the old loci.
gendata_2col<-data.frame(id=simgen$id)

reformat_genotype <- function(genotype){
	#Convert to string
	genotypeS <- as.character(genotype)
	#split on the /
	split_genotype <- strsplit(genotypeS, "/")
	a1<-split_genotype[[1]][1]
	a2<-split_genotype[[1]][2]
	
	return(c(a1,a2))
}

#Apply this function to each locus
for(i in 1:length(simgen$Gobs)){

	locus <- names(simgen$Gobs[i])
	
	#Evens are allele 1, odds are allele 2
	locus_genotypes <- unlist(lapply(simgen$Gobs[[i]], reformat_genotype))

#Extract the 2 alleles separately and add them to the gendata_2col
	gendata_2col[,substr(locus,1,nchar(locus)-1)]<-locus_genotypes[seq(1,2*length(simgen$id),2)]

	gendata_2col[,locus]<-locus_genotypes[seq(2,2*length(simgen$id),2)]
}

#Save 2 column genotypes.
#write.csv(gendata_2col, "simulated_gendata_2col.csv", row.names=F)


########################
#                      #
#      Fit model       #
#                      #
########################

GdP<-GdataPed(G=simgen$Gobs, id=simgen$id, perlocus=T)

sP<-startPed(estG=F, E1=error$E1, E2=error$E2)
pP<-priorPed(beta=list(mu=c(0,0,0,0,0,0), sigma=diag(c(rep(1+pi^(2/3),6)))))
tP<-tunePed(USsire=0.03, USdam=0.03, beta=0.3)

#Fit pedigree model to simulated data
simmod<-MCMCped(PdP, GdP, sP, tP, pP, write_postP="JOINT", DSapprox=T, jointP=F, verbose=T
,nitt=1300,burnin=300,thin=1)


#saveRDS(simmod,"simmod_190425.RDS")
##################################
#                                #
#         Assess Pedigree        #
#                                #
##################################
#Trace plots
plot(simmod$beta)
plot(simmod$USdam)
plot(simmod$USsire)

#Autocorrelation plots
autocorr.plot(simmod$beta)
autocorr.plot(simmod$USdam)
autocorr.plot(simmod$USsire)

effectiveSize(simmod$beta)
effectiveSize(simmod$USdam)
effectiveSize(simmod$USsire)

#geweke diag
geweke.plot(simmod$beta)
geweke.plot(simmod$USdam)
geweke.plot(simmod$USsire)

#Estimate vs true value
var_names<-c("Female gave birth", "Male in pack", "Female age", "Male age", "Female age^2", "Male age^2")
op<-par(mfrow=c(2,3),mar=c(5,4,1,2))
for(i in 1:6){
	xmin<-min(c(simmod$beta[,i],ped_beta[i]))
	xmax<-max(c(simmod$beta[,i],ped_beta[i]))
	densplot(simmod$beta[,i], xlim=c(xmin,xmax), main=var_names[i])
	abline(v=ped_beta[i], col="red")
}
par(op)

#Assignment probability and accuracy
pup_row_in_simped<-match(modeP(simmod$P)$P[,1],simped$ped[,1])
correct_parents<-modeP(simmod$P)$P[,2] == simped$ped[pup_row_in_simped,2] & modeP(simmod$P)$P[,3] == simped$ped[pup_row_in_simped,3]
summary(correct_parents)

hist(modeP(simmod$P)$prob)
hist(modeP(simmod$P)$prob[correct_parents],breaks=seq(0,1,0.1), main="Assignment probability, errors in red"); hist(modeP(simmod$P)$prob[!correct_parents],add=T,breaks=seq(0,1,0.1), col="red")

##########################
#                        #
#         To do          #
#                        #
##########################

#genotyping error rates
#unsampled parents

#The simped may not be correct if you are doing a multigenerational pedigree, this shorter version was used initially becaues it easy easier to check by hand.

#Should only run if the pdata is already correct, only male and female sexes.

#How do I control the number of unsampled dams and sires for the simulation?, can these change over time?

#The male age terms cause a missing covariate data error. This was because i hadnt specified to look at the correct lag.

#Does this mean that only records that pass the restrictions are examined for variables?

#Errors must match up correctly with loci

#Include a checks section of unit tests before the pedigree is constructed

#check that id is the first column in gendata

#Check prior settings and tune ped

#Assess the pedigree

#Explain why age **2 doesnt work

#!!!!! only Male and Female should be sexes.

#Should probably include a prior for the number of unsampled males and females

#save the pedigree and genotype data so that it is easy to fit into mopeds 3 and 4.
#Convert the genotype to two column format
#include columns for the old loci.