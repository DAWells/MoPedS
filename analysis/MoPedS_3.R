#Asign parentage with MasterBayes

setwd("/Users/davidwells/Dropbox/Mongoose/Pedigree/MoPedS/data")
packages<-c("MasterBayes")
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

#Add the ungenotyped individuals to gendata but with NA for all loci
ungen<-unique(pdata$id[!pdata$id%in%gendata$id])
ungen<-data.frame(id=ungen)

gendata<-merge(gendata,ungen,by="id",all=T)

#Alphabetically reorder the rows of error and columns of gendata so that they match.
gendata<-gendata[,c(1,(1+order(names(gendata)[2:87])))]
error<-error[order(error$locus),]

print("Please check that the names of loci in gendata and error correctly match up")
cbind(names(gendata)[2:87], as.character(rep(error$locus,each=2)))

########################
#                      #
#    Create model      #
#                      #
########################


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

saveRDS(PdP,"PdataPed_object.RDS")
#Join all of these restrictions and predictor variable rules together with the phenotypic data.

#Allow the rate of unsampled dams and sires to vary across time.



#Genotypic data with separate error rates for each locus.
GdP<-GdataPed(gendata, perlocus=T)

#Secify what is to be estimated by the model
#We are not estimating the genotypes and therefore we must provide rates of genotyping errors. As we have not supplied allele frequencies they are automatically calculated using extractA(). 
sP<-startPed(estG=F, E1=error$E1, E2=error$E2)

#Specify the prior
pP<-priorPed(beta=list(mu=c(0,0,0,0,0,0), sigma=diag(c(rep(1+pi^(2/3),6)))))

tP<-tunePed(USsire=0.03, USdam=0.03, beta=0.3)

########################
#                      #
#      Fit model       #
#                      #
########################

#temporary short fit data





########################
#                      #
#   Check model fit    #
#                      #
########################

modeP(modt$P)
'TM478'
dim(modt$P)

########################
#                      #
#        To do         #
#                      #
########################

#Check why the original takes 35024 from P$datesalive to calculate $timevar. I assume that it is to make it into months since records began. 
#Answer: That is correct

#How well do the newly calculated ages match with those originally used because that will impact whether masterbayes correctly uses this information.

#Check that restiction 5 is correctly specified, males must be alive 2 months before the birth, which is the conception month.

#Check that lags and time differences in the restrictions are correct

#check that the unsampled sire and dam rates are correctly allowed to vary through time.

#Use the prior estimated from Jenny's original model

#check that cohorts and timevar are correctly used

#Why is the scaled_age and sq_scaled age of BP373,374 and 376 NA?

#What to do about individuals who were genotyped more than once? Especially if the genotypes are not identical.

#What is the rate of mismatches in the existing pedigree?
#To speed up the pedigree we can pass a mismatch tolerence
#We should certainly return the number of mismatches for each assigned parent. Possibly also produce a plot of the number of mismatches vs that parents probability.

#Errors must match up correctly with loci

#Should probably include a prior for the number of unsampled males and females

#Make the brevity section good, explain how to pick those you actually want to asign.

#Allow the rate of unsampled dams and sires to vary across time.

#By specifiying jointP = FALSE to MCMCped mothers are assigned first and then fathers are assigned which can speed up comupation.