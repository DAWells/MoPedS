#Simulate a pedigree to check that the restrictions are functioning correctly.

#getXlist() returns the potential parents identified within PdataPed objects. Also includes two extra potential parents which I believe indicate an unsampled dam and sire.

packages<-c("MasterBayes", "dplyr")
lapply(packages,require,character.only=T)

#############################
#                           #
#       Generate pdata      #
#                           #
#############################

n<-40
id<-1:n

sex<-sample(c("Female","Male"),n,T)
year<-round(seq(10,14,length.out=40))
indivs<-data.frame(id=id,sex=sex,year)

###########
#Life span
###########

indivs$span<-rpois(n,1)

###########
#Yearly records
###########

#Create a record for each year they're alive.
#the + 1 ensure we have a record for those which do not survive to 1 year old.
pdata<-indivs[rep(rownames(indivs),indivs$span+1),1:3]

#Calculate individuals age in the pdata.
pdata$age<-unlist(sapply(indivs$span,function(span) 0:span))

#Update the year for each year of individuals life.
pdata$year<-pdata$year+pdata$age

###########
#Giving birth
###########

#Did an individual give birth? assign births everywhere then reset all non-potential breeders to 0.
pdata$gave_birth<-rbinom(nrow(pdata),1,0.8)
#Only females over 1 can give birth.
pdata$gave_birth[pdata$sex!="Female" | pdata$age<1]<-0


###########
#Does an individual need to be assigned parents?
###########

#Individuals are offspring to be assigned when their age is 0.  MasterBayes complains if there are no potential parents which is the case for founders.
pdata$offspring<-F
pdata[pdata$age==0 & pdata$year!=10,]$offspring<-T

#############################
#                           #
#       Construct the       #
#      PdataPed object      #
#                           #
#############################

res1<-expression(varPed(x="year",gender="Female",relational="OFFSPRING",restrict="=="))
res2<-expression(varPed(x="year",gender="Male",relational="OFFSPRING",restrict=("==")))
res3<-expression(varPed(x="offspring",restrict=F))
res4<-expression(varPed(x="gave_birth", gender="Female",restrict=1, lag = c(-1,1)))

PdP<-PdataPed(list(res1, res2, res3, res4), data=pdata, USsire=T, USdam=T, timevar=pdata$year)

s_ped<-simpedigree(PdP)


#Function to return the potential parents I expect.
get_potential_parents<-function(row){
	r1<-pdata$year <= row['year'] & pdata$year+1 >= row['year']
	r2<-pdata$year == row['year']
	r3<-pdata$offspring == F
	r4<-pdata$gave_birth == 1
	pdams<-pdata[pdata$sex=="Female" & r1 & r3 & r4,]
	psires<-pdata[pdata$sex=="Male" & r2 & r3,]
	pot_parents<-rbind(pdams,psires)
	names(pot_parents)[1]<-"potential_parent"
	pot_parents$possible_pup<-row['id']
	
	return(pot_parents)	
}

offspring<-filter(pdata,offspring==T)
pot_rents<-apply(offspring,1,get_potential_parents)
pot_parents<-bind_rows(pot_rents)

filter(pot_parents,possible_pup == "28")
pp<-getXlist(PdP)$X$'28'[1:2]
pp
filter(pdata,id==28)
filter(pdata,id%in%pp$dam.id)



pdata[pdata$id==26 & pdata$year==15,]$gave_birth<-0

hold<-
pdata<-pdata[rownames(pdata)!=26.2,]

<-hold

#########################################
#                                       # 
#     Specified very small pedigree     #
#                                       # 
#########################################
#restrictions to take only females which gave_birth 1 year either side of the offspring
pd<-data.frame(id=c("a","a","b","c","d"), sex=as.factor(c(2,2,2,2,1)), year=c(0,1,2,3,2), gave_birth=c(1,1,1,1,0), offspring=c(0,0,0,0,1)); levels(pd$sex)<-c("Male","Female")
res1<-expression(varPed(x="gave_birth",gender="Female",restrict=1,lag=c(-1,1)))
PdP<-PdataPed(list(res1), data=pd,USsire=T,USdam=T,timevar=pd$year)
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$dam.id]
pd


#Restrictions to take only females with a record 1 year before the birth
pd<-data.frame(id=c("a","b","c","d","e"), sex=as.factor(c(1,2,2,2,2)), year=c(1*60,1*60,2*60,3*60,2*60), offspring=c(0,0,0,0,1)); levels(pd$sex)<-c("Male","Female")
res1<-expression(varPed(x="year", gender="Female", relational="OFFSPRING", restrict="==60", lag=c(-60,-60)))
PdP<-PdataPed(list(res1), data=pd,USsire=T,USdam=T,timevar=pd$year)
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$dam.id]
pd


#record from the same year, and gave birth in a -1,1 window
pd<-data.frame(id=c("a","b","c","d","e","b"), sex=as.factor(c(1,2,2,2,2,2)), year=c(0,1,2,3,2,2), offspring=c(0,0,0,0,1,0), gave_birth=c(0,1,1,0,0,0)); levels(pd$sex)<-c("Male","Female")
res1<-expression(varPed(x="year", gender="Female", relational="OFFSPRING", restrict="=="))
res2<-expression(varPed(x="gave_birth", gender="Female", restrict=1, lag=c(-1,1)))
PdP<-PdataPed(list(res1, res2), data=pd,USsire=T,USdam=T,timevar=pd$year)
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$dam.id]
pd

#has a record the same year, gave birth in a -1,1 window but we use gave_birth as a variable not a restrictions
pd<-data.frame(id=c("a","b","c","d","e"), sex=as.factor(c(1,2,2,2,2)), year=c(0,1,1,1,1), offspring=c(0,0,0,0,1), gave_birth=c(0,1,1,0,0)); levels(pd$sex)<-c("Male","Female")
res1<-expression(varPed(x="year", gender="Female", relational="OFFSPRING", restrict="=="))
var1<-expression(varPed(x="gave_birth", gender="Female"))
PdP<-PdataPed(list(var1), data=pd,USsire=T,USdam=T,timevar=pd$year)
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$dam.id]
getXlist(PdP)$X[[1]]$XDus
pd

#How are the variables encoded?
#It seems they are stored in $XDus and centred but not standardised. Im not sure how it works with multiple records of an individual. There seems to be some trouble with "missing covariates" when you try.

pd<-data.frame(id=c("a","b","c","c","e"), sex=as.factor(c(1,2,2,2,2)), year=c(0,1,1,2,2), offspring=c(0,0,0,0,1), gave_birth=c(0,1,1,1,0)); levels(pd$sex)<-c("Male","Female")
res1<-expression(varPed(x="year", gender="Female", relational="OFFSPRING", restrict="==1", lag=c(-1,-1)))
var1<-expression(varPed(x="gave_birth", gender="Female"))
PdP<-PdataPed(list(res1, var1), data=pd,USsire=T,USdam=T,timevar=pd$year)
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$dam.id]
pd

#Dont know why mulitple records is causing problems. Create a new dataframe from scratch, see if that works.
#restricted to those with a record from the same year, the variable is gave birth this year or last.
pd<-data.frame(id=c("a","b","b","c","d","e"), sex=as.factor(c(1,2,2,2,2,1)), year=c(0,0,1,1,1,1), offspring=c(0,0,0,0,0,1), gave_birth=c(0,1,0,0,1,0)); levels(pd$sex)<-c("Male","Female")
res1<-expression(varPed(x="year", gender="Female", relational="OFFSPRING",restrict="=="))
var1<-expression(varPed(x="gave_birth", gender="Female", lag=c(-1,0)))
PdP<-PdataPed(list(res1, var1), data=pd, USsire=T, USdam=T, timevar=pd$year)
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$dam.id]
getXlist(PdP)$X[[1]]$XDus
pd

####################
#Testing the actual moped restrictions

#This gives wrong answers for males because res 5 insists males must have a record the year before the pup, but the other restrictions are only applied to the year OF the pup. Therefore those without records from the year of the pup fail the other restrictions.
pd<-data.frame(id=as.factor(c("m1","m2","m3","m4","f1","f2","f3","p1")), pack=c(1,1,1,1,1,1,2,1), sex=as.factor(c(2,2,2,2,1,1,1,NA)), year=c(0,0,1,1,0,1,1,1), juv=c(F,F,F,F,F,F,F,T), offspring=c(0,0,0,0,0,0,0,1));levels(pd$sex)<-c("Female","Male")
res1<-expression(varPed(x="pack", gender="Female", relational="OFFSPRING", restrict="=="))
res2<-expression(varPed(x="juv", restrict=F))
res3<-expression(varPed(x="id", relational="OFFSPRING", restrict="!="))
res4<-expression(varPed(x="year", gender="Female", relational="OFFSPRING", restrict="=="))
res5<-expression(varPed(x="year", gender="Male", relational="OFFSPRING", restrict="==1", lag=c(-1,-1)))
PdP<-PdataPed(list(res1,res2,res3,res4), data=pd, USdam=T, USsire=T, timevar=pd$year)
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$dam.id]
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$sire.id]
pd


#To correct this, ensure that all male restricts have a lag -1,-1
pd<-data.frame(id=as.factor(c("m1","m2","m3","m3","f1","f2","f2","p1")), pack=c(1,1,1,1,1,1,2,1), sex=as.factor(c(2,2,2,2,1,1,1,NA)), year=c(0,0,0,1,1,0,1,1), juv=c(F,F,T,F,F,F,F,T), offspring=c(0,0,0,0,0,0,0,1));levels(pd$sex)<-c("Female","Male")
fres1<-expression(varPed(x="pack", gender="Female", relational="OFFSPRING", restrict="=="))
fres2<-expression(varPed(x="juv", gender="Female", restrict=F))
fres3<-expression(varPed(x="id", gender="Female", relational="OFFSPRING", restrict="!="))
fres4<-expression(varPed(x="year", gender="Female", relational="OFFSPRING", restrict="=="))

mres2<-expression(varPed(x="juv", gender="Male", restrict=F, lag=c(-1,-1)))
mres3<-expression(varPed(x="id", gender="Male", relational="OFFSPRING", restrict="!=", lag=c(-1,-1)))
mres4<-expression(varPed(x="year", gender="Male", relational="OFFSPRING", restrict="==1", lag=c(-1,-1)))
PdP<-PdataPed(list(fres1,fres2,fres3,fres4,mres2,mres3,mres4), data=pd, USdam=T, USsire=T, timevar=pd$year)
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$dam.id]
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$sire.id]
pd

#######################
#There is a problem translating it to the actual pdata, possibly due to the change in lag?
pd<-data.frame(id=as.factor(c("m1","m2","m3","m3","f1","f2","f2","p1")), pack=c(1,1,1,1,1,1,2,1), sex=as.factor(c(2,2,2,2,1,1,1,NA)), datesalive=60*c(0,0,0,1,1,0,1,1), juv=c(F,F,T,F,F,F,F,T), offspring=c(0,0,0,0,0,0,0,1));levels(pd$sex)<-c("Female","Male")
dres1<-expression(varPed(x="pack", gender="Female", relational="OFFSPRING", restrict="=="))
dres2<-expression(varPed(x="juv", gender="Female", restrict=F))
dres3<-expression(varPed(x="id", gender="Female", relational="OFFSPRING", restrict="!="))
dres4<-expression(varPed(x="datesalive", gender="Female", relational="OFFSPRING", restrict="=="))

sres2<-expression(varPed(x="juv", gender="Male", restrict=F, lag=c(-60,-60)))
sres3<-expression(varPed(x="id", gender="Male", relational="OFFSPRING", restrict="!=", lag=c(-60,-60)))
sres4<-expression(varPed(x="datesalive", gender="Male", relational="OFFSPRING", restrict="==60", lag=c(-60,-60)))
PdP<-PdataPed(list(dres1,dres2,dres3,dres4,sres2,sres3,sres4), data=pd, USdam=T, USsire=T, timevar=pd$datesalive)
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$dam.id]
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$sire.id]
pd


#############################
#Successfully geting pdp constructed using restrictions and varialbes. Im not sure if records that do not pass the restrictions are evaluated for variables. it seems that they do
pd<-data.frame(id=as.factor(c("m1","m2","m3","m3","f1","f2","f2","p1")), pack=c(1,1,1,1,1,1,1,1), sex=as.factor(c(2,2,2,2,1,1,1,NA)), datesalive=60*c(0,0,0,1,1,0,1,1), juv=c(F,F,T,F,F,F,F,T), offspring=c(0,0,0,0,0,0,0,1), gave_birth=c(0,0,0,0,0,1,0,0));levels(pd$sex)<-c("Female","Male")
dres1<-expression(varPed(x="pack", gender="Female", relational="OFFSPRING", restrict="=="))
dres2<-expression(varPed(x="juv", gender="Female", restrict=F))
dres3<-expression(varPed(x="id", gender="Female", relational="OFFSPRING", restrict="!="))
dres4<-expression(varPed(x="datesalive", gender="Female", relational="OFFSPRING", restrict="=="))

sres2<-expression(varPed(x="juv", gender="Male", restrict=F, lag=c(-60,-60)))
sres3<-expression(varPed(x="id", gender="Male", relational="OFFSPRING", restrict="!=", lag=c(-60,-60)))
sres4<-expression(varPed(x="datesalive", gender="Male", relational="OFFSPRING", restrict="==60", lag=c(-60,-60)))

var1<-expression(varPed(x="gave_birth", gender="Female",lag=c(-60,0)))
PdP<-PdataPed(list(dres1,dres2,dres3,dres4,sres2,sres3,sres4, var1), data=pd, USdam=T, USsire=T, timevar=pd$datesalive)
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$dam.id]
getXlist(PdP)$id[getXlist(PdP)$X[[1]]$sire.id]
getXlist(PdP)$X[[1]]$XDus
pd


#############################
#############################

#If you do not specify a time variable it seems to group all of the entries for individuals together so if an individual meets all of the criteria collectively then they are a potential mate.

#But when you do specify a time variable that is the first filtering step so you cannot get records from other datevalues.

#possibly lag and lag_relational are the solutions here

#But will have to make certain that it isn't pooling records again but instead is correctly only accepting potential parents if they meet the criteria at a single time point

#each restriction is separately evaluated to produce a list of parents which meet the criteria.