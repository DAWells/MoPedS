#Generate phenotypic data for pedigree construction

#Create a dataframe of phenotypic data which can be used by MasterBayes to assign parentage
#Columns of the dataframe are:
	#id, the individual's name
	#pack, the social group that individual is in at that point in time
	#datesalive, each record refers to an individual at a given time point. Traditionally, a record is made every 30 days for an individuals life. This time period is refered to as a space. This value is numeric
	#date is the datesalive value in the date format yyyy-mm-dd
	#sex, 
	#givenbirth, did the individual give birth during this space.
	#offspring, is this the first record of the individual and therefore do they need their parentage assigned. Note this may be changed later so that only specified individuals need to be assigned parentage
	#age_space, the rough age in days of the individual at the start of space. It is only rough because it is counted from the start of the space they were born during, not from their date of birth.
	#juv, is the individual a juvenile? i.e. are they less than 6 months old. For the purposes of our pedigree individuals who are juveniles are not considered as possible parents.
	#scaled_age is age_space converted into months and standardised and scaled by the mean and standard deviation. For immigrants the age_space is unknown and so set to the average for their sex at first observations
	#sq_scaled_age is scaled_age squared

setwd("/Users/davidwells/Dropbox/Mongoose/Pedigree/MoPedS/data")
library(dplyr)

#load the life history data prepared in MoPedS_1.R
lhdata<-read.csv("rawdata/lhdata_march_2019.csv")
#Remove some strange records based on Jenny's original code
lhdata<-lhdata[which(lhdata$indiv!="UM PUP"),]
lhdata<-lhdata[which(lhdata$indiv!="UM PUP1"),]
lhdata<-lhdata[which(lhdata$indiv!="UM PUP2"),]
lhdata<-lhdata[which(lhdata$indiv!="UM PUP3"),]
lhdata<-lhdata[which(lhdata$indiv!="UM AD"),]
lhdata<-lhdata[which(lhdata$indiv!="UM F"),]
lhdata<-lhdata[which(lhdata$indiv!="UM SUB"),]
lhdata<-lhdata[which(lhdata$indiv!="UMSUB"),]
lhdata<-lhdata[which(lhdata$indiv!="UNK"),]
#appear
lhdata<-lhdata[which(lhdata$indiv!="UM"),]
lhdata<-lhdata[which(lhdata$indiv!="TF UNK"),]
#microsatellite data
gendata<-read.csv("rawdata/BMRP_Genotype_data_July_2016.csv")

######################
#                    #
#   Generate data    #
#                    #
######################

#Get all recordings of individuals born, must use BORN and START to ensure it is individuals born and not litters
#In the original pedigree Jenny also limited to daten>=35144 (20th March 1996), presumably because the lhdata before then is very sparse.
borns<-filter(lhdata,code=="BORN" & stend=="START")
#Get records of individuals giving birth
births<-filter(lhdata,code=="BIRTH")

########################
#Create a vector of dates roughly from the first born to the last lhdata record, with a defined spacing (30 in the original pedigree). Also include 120 days before the first born because information 3 months before the born is used in the pedigree. Also include a full spacing at the end 
spacing<-30
if(120%%spacing!=0){
	stop("Spacing must fit exactly into 120")}
dates<-seq(min(borns$daten)-120,(max(lhdata$daten)+spacing),spacing)
#Each of these dates refers to itself and the spacing after it

#Get all of the starts/restarts/stops/ends
st<-numeric()
st[which(lhdata$stend=="START")]<-1
st[which(lhdata$stend=="RESTART")]<-1
st[which(lhdata$stend=="STOP")]<-2
st[which(lhdata$stend=="END")]<-2

##starts/restarts
stinds2<-lhdata[which(st==1),]
##stops/ends
endinds<-lhdata[which(st==2),]

#For each start or restart in stinds2 get either the earliest stop/end for the same individual, in the same pack, that is greater than or equal to the date of the start. If there is no recorded end or stop, add a fake end.
add_min_or_fake_end<-function(row){
#This function is to be used within apply()
	id<-row['indiv']
	pack<-row['pack']
	daten<-row['daten']

	same_id<-endinds$indiv == id
	same_pack<-endinds$pack == pack
	after<-endinds$daten>=daten
	
	applicable_ends <- endinds[same_id & same_pack & after,]

#If no applicable ends are found, add a fake end. 155 is an arbitrary value.
#Note that individuals with this fake end cannot have their true age calculated (because they are not dead yet)
	if(nrow(applicable_ends)==0){
		endn<-max(lhdata$daten)+155
	}else{
#Otherwise use the earliest applicable end
		endn<-min(applicable_ends$daten)
	}
	return(endn)
}
#These are the end dates which match up with the rows in stinds2
endn<-apply(stinds2,1,add_min_or_fake_end)
stinds2$endn<-endn

#Note these fake ends are not included in the pdata as they are after the end of dates, which only extends to the final date in the lhdata.
##################

#Which individuals do not have genetic data?
no_gen<-sum(!unique(stinds2$indiv)%in%gendata$id)
print(paste(no_gen, " Mongooses with start or restart records have no genetic data",sep=""))

################
#Create a matrix where each column is a mongoose and each row is a date
#sets the dates which are outside of the start end period to NA
keep_relevant_dates<-function(row){
#This function is to be used within apply()

#which dates in the ith column are before the start?
#+spacing ensures that the actual start date is in the space following the first date which is not set to NA
	b4_start<- dates+spacing < (row['daten'])
#which dates in the ith column are after the end?
	after_end<- dates > row['endn']
#the dates outside of the start end period are set to NA
	relevant_dates<-dates
	relevant_dates[b4_start | after_end]<-NA
	
	return(relevant_dates)	
}



aliveframe<-apply(stinds2, 1, keep_relevant_dates)

#NEW For Loop, use this OR keep_relevant_dates
#In each column makes the dates that are outside of the stinds2 endn period NA
#for(i in 1:ncol(aliveframe)){
#which dates in the ith column are before the start?
#-spacing ensures that the actual start date is in the space following the first date which is not set to NA
#	b4_start<- dates < (stinds2$daten[i]-spacing)
#which dates in the ith column are after the end?
#	after_end<- dates > endn[i]
#the dates outside of the start end period are set to NA
#	aliveframe[b4_start | after_end,i]<-NA
#}

###########################
#Create a series of vectors which will become columns in the pdata frame.

###############
#    Dates    #
###############

#Flatten aliveframe (which sould be a matrix) into a vector
datelist<-as.vector(aliveframe)
#remove the NAs
datesalive<-datelist[!is.na(datelist)]

#add a date column which is the $datesalive column formatted as date
date<-as.Date("1899-12-30")+datesalive

#####################
#     ID & pack     #
#####################

#create a vector of individual mongooses which will line up with the dates such that they match up with the dates that they are alive
id<-rep(stinds2$indiv,each=length(dates))[!is.na(datelist)]
#create a vector of the packs that the individual was in at that date
packs<-rep(stinds2$pack, each=length(dates))[!is.na(datelist)]

#############
#    Sex    #
#############
#create a vector of sexes for these mongooses
#There are some gaps in the sex of individuals
#eg GF175 in the march 2019 release of the life history data
filter(lhdata,indiv=="GF175")
#Her first entry (with the current ordering of lhdata) has sex set as "". This means that if we look up her sex, using match(), R takes the first instance which means her sex would be "".

#create a table of individuals sexes which ignores "" as entries for sex
lhdata_individuals<-filter(lhdata,indiv%in%id)

table_of_sexes<-summarise(group_by(lhdata_individuals,indiv),sexes=paste(unique(sex),collapse=""))


non_FMP<-table_of_sexes[!table_of_sexes$sexes %in% c("F","M","P"),]

if(nrow(non_FMP)>0){
	print(non_FMP)
	non_FMP_warning<-" Individuals have a non F, M, or P sex. This may because their sex is initially entered as P. Please check and correct this."
	stop(paste(nrow(non_FMP),non_FMP_warning))
}

sex<-table_of_sexes$sexes[match(id,table_of_sexes$indiv)]
sex<-as.factor(sex)

#When assigning parentage M and F can cause confusion as it is not clear if they refer to male and female or mother and father
levels(sex)[levels(sex)=="F"]<-"Female"
levels(sex)[levels(sex)=="M"]<-"Male"



###############
#  Offspring  #
###############

#The individuals which we have born records for
bornIDs<-unique(borns$indiv)

#if they are an offspring in need of assigning parentage at this time we will set to 1
offspring<-rep(0,length(id))

which_need_parentage<-function(x){
#function for use in apply
#the relevent individual
	correct_id<-id==x
#the earliest date for that individual
	min_date<-datesalive==min(datesalive[correct_id])

#which records are the first record of an individual with a known birthday?
	wnp<-which(correct_id & min_date)
	if(length(wnp)>1){
		stop("Each individual should only be 'offspring' once")
	}
	return(wnp)
}
#this tells us which records need parentage assigned
need_parentage<-sapply(bornIDs,which_need_parentage)
#this line actually changes these records
offspring[unlist(need_parentage)]<-1

#A for loop to do the same thing, they actually seem to be ~speed
#for (i in 1:length(bornIDs)){
#	indiv<-id==bornIDs[i]
#	min_date<-datesalive==min(datesalive[indiv])
#	which(indiv & min_date)
#	offspring[indiv & min_date]<-1
#}

###############
#  Juvenile   #
###############
#A vector of whether or not individuals are over 6 months old, this is used in the pedigree to determine if they are old enough to have reproduced
#Individuals without born records are assumed to be over 6 months

#birthspaces is taken as the earliest datesalive for that individual, although when spacing is not 30 it will not neccessarily be a month
birthspaces<-as.vector(tapply(datesalive,id,min))[match(id,levels(id))]
#this is set to NA for those without birth records, since we dont know when they were born
birthspaces[!id%in%bornIDs]<-NA

#This age is not exact as it is from the start of the birth space, not from th birth day 
age_space<-datesalive-birthspaces

juv<-factor(levels=c("TRUE","FALSE"))
juv[1:length(age_space)]<-"FALSE"
juv[id%in%bornIDs & age_space <=(30*6)]<-T

###################
#   Given birth   #
###################
givenbirth<-numeric()
givenbirth[1:length(juv)]<-0

for (i in 1:nrow(births)){
	dob <- births$daten[i]
	mother <- births$indiv[i]
	
	space_starts_b4_birth<-datesalive<=dob
	last_space_b4_birth<-(datesalive+spacing)>dob
	correct_id <-id == mother
	
	givenbirth[space_starts_b4_birth & last_space_b4_birth & correct_id]<-1
}

#an apply function to do the same job but it runs very slowly
#gives_birth<-function(row){
#	dob<-row['daten']
#	mother<-row['indiv']
	
#	space_starts_b4_birth<-datesalive<=dob
#	last_space_b4_birth<-(datesalive+spacing)>dob
#	correct_id <-id == mother

#	which(space_starts_b4_birth & last_space_b4_birth & correct_id)
#}

#wgb<-apply(births, 1, gives_birth)
#givenbirth[wgb]<-1

################################
#                              #
#       Bind the vectors       #
#       into a dataframe       #
#                              #
################################

pdata<-data.frame(id,pack=packs,date,datesalive,sex, givenbirth,offspring,age_space,juv)


########################
#                      #
#       Immigrant      #
#                      #
########################

#To create the COLONY input file we need to know which immigrants we would like to assign to sibship groups. Therefore we create a variable declaring which individuals are immigrants to the population.
#Any individual which we do not have a birth date for counts as an immigrant because without that we cannot assign parentage to them using MasterBayes because we cannot compile candidate parents.
pdata$immigrant <- as.numeric(is.na(pdata$age_space))

#For COLONY we only need the first record of an immigrant to get their siblings.

#Get the first datesalive for each individual
first_records<-tapply(pdata$datesalive,pdata$id,min)

get_new_immigrants<-function(row){
	#Get the individual's first record
	first_record<-first_records[row['id']]
	#Is current record the first record?
	current_first <- first_record == row['datesalive']

	#Is this individual an immigrant?
	#The as.locigal and as.numeric just maintain the data type locial.
	immi<-as.logical(as.numeric(row['immigrant']))
	#Individual is a new immigrant if there current record is th first record AND they are an immigrant.
	new_immigrant <- current_first & immi
	
	return(new_immigrant)
}
pdata$new_immigrant<-apply(pdata, 1, get_new_immigrants)


################################
#                              #
#         Center and           #
#       standardise ages       #
#                              #
################################

#$age_space is scaled based on sex. Mathematically this means: (age-mean(age))/sd(age). A separate mean and standard deviation are calculated for males and females.
#Reminder, age is in days
#Imigrants to the population have $age_space = NA so we replace this with the average for that sex. However, for subsequent measures we know that they have continuned to age so that information is included. Note though that this was not the case for Jenny's original pedigree or David's extension 1, imigrant's age was fixed as the average.

av_male_age<-mean(pdata$age_space[pdata$sex=="Male"],na.rm=T)
av_female_age<-mean(pdata$age_space[pdata$sex=="Female"],na.rm=T)

sd_male_age<-sd(pdata$age_space[pdata$sex=="Male"],na.rm=T)
sd_female_age<-sd(pdata$age_space[pdata$sex=="Female"],na.rm=T)


immigrant_age<-function(row){
	indiv<-row['id']
	#immigrants have NA for age_space as we do not know when they were born
	immigrant<-is.na(row['age_space'])
	
	#If they are not an immigrant return their existing age_space
	if(!immigrant){
		return(as.numeric(row['age_space']))
	}else{
		#get all the records for this individual
		all_p_records<-pdata[pdata$id==indiv,]
		#get the date of their first record
		first_date<- min(all_p_records[,'datesalive'])
		#get the date of the relevant record
		record_date <- as.numeric(row['datesalive'])
		
		#If male, return the average male age + the length of time since their first record
		if(row['sex']=="Male"){
			return(as.numeric(av_male_age+(record_date-first_date)))
		}

		#If female, return the average female age + the length of time since their first record
		if(row['sex']=="Female"){
			return(as.numeric(av_female_age+(first_date-record_date)))
		}else{
			#If they have a different sex return NA
			return(NA)
		}
	}
}

age<-apply(pdata,1,immigrant_age)
#convert age into months for compatibility with the original pedigree
age<-age/30

Male<-pdata$sex=='Male'
Female<-pdata$sex=='Female'

pdata$scaled_age[Male]<-scale(age[Male],T,T)
pdata$scaled_age[Female]<-scale(age[Female],T,T)

pdata$sq_scaled_age<-pdata$scaled_age**2


##########################
#                        #
#         Checks         #
#                        #
##########################

#$age_space should be 0 when offspring is 1
if(sum(filter(pdata,offspring==1)$age_space==0) != sum(pdata$offspring)){
	stop("$Offspring should all have $age_space == 0")
}

#There should only be the three sexes
stop("This check for the 3 sexes is untested\nand needs completing\n\n!!!!")

lhdata_individuals<-filter(lhdata,indiv%in%id)

table_of_sexes<-summarise(group_by(lhdata_individuals,indiv),sexes=paste(unique(sex),collapse=""))

non_FMP<-table_of_sexes[!table_of_sexes$sexes %in% c("F","M","P"),]

if(nrow(non_FMP)>0){
	print(non_FMP)
	non_FMP_warning<-" Individuals have a non F, M, or P sex. This may because they have more than one sex in lhdata, e.g. initially entered as P. Please check and correct this."
	stop(paste(nrow(non_FMP),non_FMP_warning))
}


#Do any individuals over 6 months have P as their sex?
#unsexed candidates
unsexed_candidates<-unique(filter(pdata,juv==F & !sex%in%c("Male","Female"))$id)
if(length(unsexed_candidates)>0){
	print(paste(length(unsexed_candidates)," individuals are over six months, and therefore potential parents, but have sex 'P'. This means we cannot assign them as parents."))
}

#Are any individuals born more than once?
if(sum(duplicated(borns$indiv))>0){
	stop("Some offspring are born more than once")
}
born2x<-borns[duplicated(borns$indiv),]$indiv
filter(borns,indiv%in%born2x)

#Have any males given birth?
male_births<-sum(filter(pdata,sex=="Male")$givenbirth)
if(male_births >0){
	stop("Some males have given birth")
}

#Do mothers give birth after their last END?
#the start of all individual's (in the pdata) last recorded space
#Note that there are NAs in last_pdate because it calculates the max for all levels of id but those which are not in id are assigned NA
last_pdate<-tapply(datesalive,id,max)
#all the mothers identified from the lhdata
mothers<-births$indiv
#the dates they gave birth
bds<-births$daten
#Get the last pdata daten for each of the mothers
mums_last_pdate<-last_pdate[match(mothers,levels(id))]
#Did they give birth after the end of their last space?
dead_b4_birth<-bds>mums_last_pdate+spacing
summary(dead_b4_birth)
if(sum(dead_b4_birth,na.rm=T)>0){
	stop(paste(sum(dead_b4_birth,na.rm=T)," mothers are dead before they gave birth. Please check the life history data and correct this."))
}
if(sum(is.na(dead_b4_birth))>0){
	stop(paste(sum(is.na(dead_b4_birth))," mothers do not have final dates in the pdata. Possibly because in the life history they have been given place holder names. If we have no further information on who gave birth then remove it from the life history data."))
}
plot(mums_last_pdate,bds,col=dead_b4_birth+1)

zombie_birth<-cbind(births[,1:5],mums_last_pdate,dead_b4_birth)
zombie_birth$dif<-mums_last_pdate-births$daten
filter(zombie_birth,dead_b4_birth)

if( 0 < sum(zombie_birth$dead_b4_birth,na.rm=T)){
	print(filter(zombie_birth,dead_b4_birth))
	stop("Females giving birth after their final END/STOP in the life history data")
}

#Check that all individuals with sufficient lhdata are in the pdata, 
#Who is not in pdata but has a code which implies they are an individual
#summary(lhdata$code)
mongooses<-filter(lhdata,code%in%c("ABORT","BIRTH","ADIED","DIED","EM","FPREG"))
non_starters<-mongooses[!mongooses$indiv%in%pdata$id,]
if(nrow(non_starters)>0){
	print(non_starters)
	stop(paste(nrow(non_starters),' apparent individuals from lhdata have no start dates, and are therefore ommited from the pdata',sep=""))
}

#Which individuals with start records do not have genetic data?
no_gen<-unique(stinds2[!stinds2$indiv%in%gendata$id,]$indiv)
print(paste(length(no_gen), " Mongooses with start or restart records have no genetic data",sep=""))
#Life histories of all individuals with no genetic data
lh_no_gen<-filter(lhdata,indiv%in%no_gen)
#get the year of the earliest life history record of each individual without genetic data
nogen_start<-summarise(group_by(lh_no_gen,indiv),year=format(min(as.Date(date)),"%Y"))
nogen_start$year<-as.factor(nogen_start$year)
#Bar plot showing how many lifetime records began per year that have not been genotyped
plot(nogen_start$year,ylab="Number of ungenotyped  starts")

#How many of these ungenotyped individuals survive to over 6 months? And are therefore potentially parents.
nogen_6m<-unique(filter(pdata,juv==F & !id%in%gendata$id)$id)
print(paste(length(nogen_6m), " mongooses over 6 months  have no genetic data",sep=""))
lh_nogen_6m<-filter(lhdata,indiv%in%nogen_6m)
#get the year of the first first life history record where they were over 6 months for each individual with no genetic data
nogen_6m_start<-summarise(group_by(lh_nogen_6m,indiv),year=format(min(as.Date(date)),"%Y"))
plot(as.factor(nogen_6m_start$year),ylab="NO.  ungenotyped individuals reaching 6 months")
text(x=18,y=60,labels="Note: only the 1st instance for\n each individual is shown")

#Which individuals that have not been genotyped are still alive?
still_alive<-filter(stinds2, endn==max(lhdata$daten)+155)

plot(nogen_start[nogen_start$indiv%in%still_alive$indiv,]$year, ylab="Number of ungenotyped  starts still alive")


#Do any individuals have genetic data but no life history data?
no_lh<-!gendata$id%in%unique(stinds2$indiv)
print(paste(sum(no_lh)," mongooses have genetic data but no life history data. Individuals without life history data cannot be included in the pedigree.

1f2, 1i4 etc refer to samples genotyped in 2016 that the name had rubbed off of.

The individuals TM or TF401-406 do not show up because the lhdata records them as TP not TM40_ or TF40_. Please correct this so the same individuals have the same name accross datasets.",sep=""))

#Some pups have inconsistent names beause they died at 1 day old and were given TP___ names instead of TM___/TF___. The litter is T1406
c('TF405','TM402','TF403','TF404','TF406','TM401')

gendata$id[no_lh]

######################################
#                                    #
#          Subset only the           #
#        relevant individuals        #
#                                    #
######################################

#For brevity only specify the new offspring to be assigned
!!!make easy to use

pdata[format.Date(pdata$date,"%Y")!=2015,]$offspring<-0
pdata[format.Date(pdata$date,"%Y")!=2015,]$new_immigrant<-F
#Only look at 2015 as examples for brevity.
pdata<-pdata[format.Date(pdata$date,"%Y") %in% c(2015,2014),]


######################
#                    #
#     Save pdata     #
#                    #
######################

#write.csv(pdata,"phenotypic_data_oct_2018.csv",row.names=F)

####################
   #Issues to resolve
####################

#set only the new pups to be assigned otherwise it will take for ever...
#Output a report with all the information about hte data, including how many individuals do not have genotypes
#Do all pdatas for an individual join up correctly?


###########
#         #
#  To do  #
#         #
###########

#There are some errors in the pedigree, some parents are MUCH younger than their "offspring" eg DM165 and P6

# deal with the issue noted in PROBLEM.R

#Report the number of ungenotyped individuals, specifying also the number of ungenotyped offspring, and candidate parents.

#Lots of the scaled ages are NA

#Choose the file names.