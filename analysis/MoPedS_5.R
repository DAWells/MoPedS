#190728
#Extract the COLONY results and combine them with the MasterBayes pedigree.


setwd("/Users/davidwells/Dropbox/Mongoose/Pedigree/MoPedS/data")
packages<-c("MasterBayes", "dplyr", "pedantics")
lapply(packages,require,character.only=T)

#Load the previous pedigree which is to be updated.
previous_ped<-read.csv("rawdata/previous_pedigree.csv")[,c('id','dam','sire')]
#"Fix pedigree", i.e. ensure that all individuals in the pedigree appear in $id, not just as parents.
previous_ped<-fixPedigree(previous_ped)

#Load MasterBayes pedigree output.
mb<-readRDS("simmod_190425.RDS")

#load phenotypic data created in MoPedS_2.R
pdata<-read.csv("phenotypic_data_oct_2018.csv")
#Add to this list of pup id the new_immigrants because COLONY takes them as offspring to assign sibships.
new_immigrants<-as.character(filter(pdata, new_immigrant==T)$id)


#####################################
#                                   #
#   Load the COLONY results files   #
#                                   #
#####################################
#The COLONY output files must be moved to a folder named "COLONY" in the "data" folder of the MoPedS project.

output_name<-readline(prompt="What is the name of the Colony results file?")
colony_test_file

#COLony PATernity
colpat<-read.csv(paste("COLONY/", output_name, ".Paternity", sep=""))

#COLony MATernity
colmat<-read.csv(paste("COLONY/", output_name, ".Maternity", sep=""))

#COLony Parent Pair
colpp<-read.csv(paste("COLONY/", output_name, ".ParentPair", sep=""))

#COLony best CONFIGeration
colconfig<-read.table(paste("COLONY/", output_name, ".BestConfig_Ordered", sep=""), comment.char="", header=T)

#COLony best Full Sibling Dyads
colfsd<-read.csv(paste("COLONY/", output_name, ".FullSibDyad", sep=""))

#COLony best Full Sibling Family
colfsf<-read.table(paste("COLONY/", output_name, ".BestFSFamily", sep=""), comment.char="", header=T)
names(colfsf)[4]<-"Family"

##################################
#                                #
#      #Correct dummy names      #
#                                #
##################################

#The dummy parents named in the COLONY best configuration file must be changed.
#Unsampled sires are labelled with *number and unsampled dams are labelled with #number. These numbers start from 1 for each colony run. The existing pedigree already includes some dummy parents but #1 will not be the same individual. Therefore dummy parents must be given higher numbers to avoid confusion.

#Get existing dummy parents
#PREvious DUMMY DAMS
pre_dummy_dams<-previous_ped$id[substr(as.character(previous_ped$id), 1, 1) == "#"]

#Dummys from the first pedigree extension (2016) were differentiated from the original pedigree by an x, eg "#x46". These cannot be converted to numeric so are ignored.
#Get the highest number dummy dam
max_dummy_dam<-max(as.numeric(substr(pre_dummy_dams,2,999)),na.rm=T)#Ignore the warning about Coersed NAs

#NEW DUMMY DAMS
#Get the dummy dams in the new COLONY results file
unedited_new_dummy_dams<-levels(colconfig$MotherID)[substr(as.character(levels(colconfig$MotherID)), 1, 1) == "#"]

#Function to create new dummy dam levels which carry on from the max_dummy
update_dummy_id<-function(dummy_id, parent){
	if(parent == "Dam"){
		current_index<-as.numeric(substr(dummy_id,2,999))
		new_index<- current_index + max_dummy_dam
		new_dummy_id<- paste('#', new_index,sep="")
	}else if(parent == "Sire"){
		current_index<-as.numeric(substr(dummy_id,2,999))
		new_index<- current_index + max_dummy_sire
		new_dummy_id<- paste('*', new_index,sep="")
	}else{
		stop('parent must be "Dam" or "Sire.')
	}

	return(new_dummy_id)
}

#Replace the current dummy names with new ones so there is no overlap with existing dummy names
levels(colconfig$MotherID)[substr(as.character(levels(colconfig$MotherID)), 1, 1) == "#"]<-update_dummy_id(unedited_new_dummy_dams, "Dam")

#################
# Sires
#################
#PREvious DUMMY SIRES
pre_dummy_sires<-previous_ped$id[substr(as.character(previous_ped$id), 1, 1) =="*"]
max_dummy_sire<-max(as.numeric(substr(pre_dummy_sires,2,999)),na.rm=T)#Ignore the warning about Coersed NAs

unedited_new_dummy_sires<-levels(colconfig$FatherID)[substr(as.character(levels(colconfig$FatherID)), 1, 1) == "*"]
levels(colconfig$FatherID)[substr(as.character(levels(colconfig$FatherID)), 1, 1) == "*"]<-update_dummy_id(unedited_new_dummy_sires, "Sire")

###############################
#                             #
#       Create marginal       #
#       COLONY pedigree       #
#                             #
###############################

#Marginal means assign maternity and paternity based on their probability alone, not on the joint probability.
#PEDigree based on MARginal COLony probabilities
ped_mar_col<-data.frame(id=colconfig$OffspringID, dam=NA, sire=NA, dam_prob=NA, sire_prob=NA)

#Add sires.
ped_mar_col$sire<-colpat$InferredDad1[match(ped_mar_col$id, colpat$OffspringID)]
#Add sires' probability
ped_mar_col$sire_prob<-colpat$ProbDad1[match(ped_mar_col$id, colpat$OffspringID)]

#Add dams.
ped_mar_col$dam<-colmat$InferredMum1[match(ped_mar_col$id, colmat$OffspringID)]
#Add dams' probability.
ped_mar_col$dam_prob<-colmat$ProbMum1[match(ped_mar_col$id, colmat$OffspringID)]

###############################
#                             #
#        Create joint         #
#       COLONY pedigree       #
#                             #
###############################

#Joint means the probability applies to the parent pair rather than either individual.

ped_joint_col<-data.frame(id=colconfig$OffspringID, dam=NA, sire=NA, joint_prob=NA)

#The match function returns the first match. Where there are multiple possible parent pairs COLONY puts them in order of probability. Therefore as long as the colpp dataframe is ordered by probability within offspring match will return the correct pair.
#The file should not need to be edited in that case but to be sure, order by probability within offspring.
colpp<-colpp[order(colpp$OffspringID, -colpp$Probability),]

#Add the most likely joint sire.
ped_joint_col$sire<-colpp$InferredDad[match(ped_joint_col$id, colpp$OffspringID)]
#Add the most likely joint dam.
ped_joint_col$dam<-colpp$InferredMum[match(ped_joint_col$id, colpp$OffspringID)]
#Add the joint probability of the parent pair.
ped_joint_col$joint_prob<-colpp$Probability[match(ped_joint_col$id, colpp$OffspringID)]

###################################
#                                 #
#          Create joint           #
#       MasterBayes pedigree      #
#                                 #
###################################

#modeP gets the details of the most probable parentage assignments. By setting marginal to FALSE we get the probability of both parents being assigned.
ped_joint_mb<-as.data.frame(modeP(mb$P, marginal=F)$P)
names(ped_joint_mb)<-c("id","dam","sire")
ped_joint_mb$joint_prob<-modeP(mb$P, marginal=F)$prob

####################################
#                                  #
#          Create marginal         #
#       MasterBayes pedigree       #
#                                  #
####################################

#By setting marginal to TRUE, modeP returns the probability of the most likely dam and sire independant of the other parent.
ped_mar_mb<-as.data.frame(modeP(mb$P, marginal=T)$P)
names(ped_mar_mb)<-c("id","dam","sire")
ped_mar_mb$dam_prob<-modeP(mb$P,marginal=T)$prob
ped_mar_mb$sire_prob<-modeP(mb$P,marginal=T)$prob.male

##############################
#                            #
#         Joint MB           #
#   probability threshold    #
#                            #
##############################

#The current 4 pedigrees, marginal and joint for COLONY and MasterBayes, include all of the most likely parents. However, some may still be quite low probability. In this section we remove the parents below our probability threshold from the pedigree.

threshold<-0.8
#Start with the MasterBayes joint pedigree.
#NOTE that by starting with the MasterBayes pedigree the non-offspring (founders and immigrants) are initially not included as they were only used in COLONY.
ped<-ped_joint_mb

#Dams and sires with a joint probability below the threshold are set to NA.
low_joint_prob <- ped$joint_prob < threshold
ped[low_joint_prob, c("dam","sire")]<-NA

##############################
#                            #
#      Marginal MB           #
#   probability threshold    #
#                            #
##############################

#Although both parents were not confidently assigned we may be able to assign one parent confidently.
#Add the marginal probabilities.
ped$dam_prob<-ped_mar_mb$dam_prob[match(ped$id, ped_mar_mb$id)]
ped$sire_prob<-ped_mar_mb$sire_prob[match(ped$id, ped_mar_mb$id)]

#If the parent is NA, and their marginal probability is greater than the threshold, assign the marginal parent.

assign_marginal_parents<-function(row){
	#If the individual alread has a dam, keep it.
	if(!is.na(row['dam'])){
		dam<-row['dam']
	#If they do not have an assigned dam and the marginal dam probability exceeds the threshold, set dam to the marginal dam.
	}else if(row['dam_prob'] >= threshold){
		dam<- as.character(ped_mar_mb[ped_mar_mb$id == row['id'],'dam'])
	}else{#Otherwise, no dam is assigned.
		dam<-NA
	}

	#If the individual alread has a sire, keep it.
	if(!is.na(row['sire'])){
		sire<-row['sire']
	#If they do not have an assigned sire and the marginal sire probability exceeds the threshold, set sire to the marginal sire.
	}else if(row['sire_prob'] >= threshold){
		sire<- as.character(ped_mar_mb[ped_mar_mb$id == row['id'],'sire'])
	}else{#Otherwise, no sire is assigned.
		sire<-NA
	}
	
	return(c(dam, sire))
}

ped[,c('dam','sire')]<-t(apply(ped,1, assign_marginal_parents))

#The MasterBayes consensus pedigree, used to check full siblings later.
MBped<-ped

##############################
#                            #
#    Identifying siblings    #
#                            #
##############################

#The sibships should only be added relevant to the masterbayes assignments. Therefore sibships must be assessed before the joint and marginal COLONY parents are added.

#So far niether COLONY pedigree incorporates information about full siblings which is the principal reason we use COLONY, to infer relatedness when parents are unsampled or cannot be assigned, as is the case with immigrants.

#Each family in colfsf is currently a single string of characters, this is not useful. Convert each family member into a separate column.
#Split each Family into a vector of its members
Family<-as.character(colfsf$Family)
Family_list<-strsplit(Family,",",fixed=T)
#Get the size of each full sibship.
fs_size<-sapply(Family_list,length)
#Get the size of the largest full sibship.
max_fs_size<-max(fs_size)
#Make all Family vectors equal length
Family<-lapply(Family_list,function(item){length(item)<-max_fs_size;return(item)})
#Bind them together into a data frame.
Family<-do.call(rbind.data.frame,Family)
names(Family)<-paste("Fullsib",1:max_fs_size,sep="")
#Pair the new Family data frame with their index and probabilities.
colfsf<-cbind(colfsf[,1:3],Family)

#Match each individual with their Full Sibling INDEX.
fsindex<-data.frame(index=rep(colfsf$FullSibshipIndex, fs_size), iprob=rep(colfsf$Prob.Inc.., fs_size), id=unlist(Family_list))


#Checking that COLONY full sibs have the same parents in masterbayes too.

#Full Sibling Index
ped$fsindex<-fsindex$index[match(ped$id, fsindex$id)]
#Confidence in the full sibling group
ped$fs_iprob <- fsindex$iprob[match(ped$id, fsindex$id)]

#Grouping by full sibling index, all offspring should have the same parents, this includes NAs because if the true parent was sampled they should be assigned for all offspring by MasterBayes.
#Immigrants and founders are not assigned parents by MasterBayes so they can be NA, but if all their full siblings have the same parents then they can be given those parents.
#This is only true where the full sibling iprob > threshold
#Note, at this point the immigrants and founders included in the COLONY analysis are not included in the MB analysis. Therefore, all full siblings currently in ped should have exactly the same parent pair.

fs_consistancy <- summarise(group_by(filter(ped, fs_iprob>threshold), fsindex),
	confident_sibship_size = n(),
	fs_iprob = unique(fs_iprob),
	n_MB_dams = length(unique(dam)),
	n_MB_sires = length(unique(sire)),
	#Number of unique MasterBayes assigned parent pairs.
	n_MB_pp = length(unique(paste(dam,sire))))

#Do all offspring have the same parent pair within fullsibships?
#Essentially does MasterBayes agree that they are full siblings?
fs_consistancy$mb_consistant <- fs_consistancy$n_MB_pp == 1

########################
#                      #
#       Add new        #
#     immigrants       #
#                      #
########################

#Add the new immigrants which COLONY tried to assign parents/siblings to.
new_immigrants_df<-data.frame(id=new_immigrants, dam=NA, sire=NA, joint_prob=NA, dam_prob=NA, sire_prob=NA, fsindex=NA, fs_iprob=NA)

ped <- rbind(ped, new_immigrants_df)

##############################
#                            #
#       Joint COLONY         #
#   probability threshold    #
#                            #
##############################

#Add Joint COLONY parents above the threshold

ped$joint_col_prob<-ped_joint_col$joint_prob[match(ped$id, ped_joint_col$id)]

#Function which assigns the joint parent pair from COLONY if both parents are still unassigned and $joint_col_prob exceeds the threshold
assign_joint_col_parents <- function(row){
	#If both parents are NA check the joint_col_probabilty against the threshold
	if(is.na(row['dam']) & is.na(row['sire'])){
		#If probability is not NA and over the threshold
		if(row['joint_col_prob'] > threshold & !is.na(row['joint_col_prob'])){
			#Joint parent pair
			jpp<-ped_joint_col[ped_joint_col$id == as.character(row['id']),]
			dam <- as.character(jpp$dam)
			sire <- as.character(jpp$sire)
		}else{
			dam<-NA
			sire<-NA
			}
	}else{
		dam<-row['dam']
		sire<-row['sire']
	}
	return <- c(dam, sire)
}


ped[,c('dam','sire')]<-t(apply(ped,1,assign_joint_col_parents))

#Those with dummy parents assigned are only given # or *. Here we specify which dummy parent using colconfig.


##############################
#                            #
#      Marginal COLONY       #
#   probability threshold    #
#                            #
##############################

#Add Marginal COLONY parents above the threshold.

#Add the marginal COLONY probabilities to ped.
ped$dam_col_prob <-ped_mar_col$dam_prob[match(ped$id, ped_mar_col$id)]
ped$sire_col_prob <-ped_mar_col$sire_prob[match(ped$id, ped_mar_col$id)]

assign_marginal_parents_col<-function(row){
	#If the individual alread has a dam, keep it.
	if(!is.na(row['dam'])){
		dam<-row['dam']
	#If they do not have an assigned dam and the marginal dam probability exceeds the threshold, set dam to the marginal dam.
	}else if(row['dam_prob'] >= threshold &
	!is.na(row['dam_prob'])){
		dam<- as.character(ped_mar_col[ped_mar_col$id == row['id'],'dam'])
	}else{#Otherwise, no dam is assigned.
		dam<-NA
	}

	#If the individual alread has a sire, keep it.
	if(!is.na(row['sire'])){
		sire<-row['sire']
	#If they do not have an assigned sire and the marginal sire probability exceeds the threshold, set sire to the marginal sire.
	}else if(row['sire_prob'] >= threshold &
	!is.na(row['sire_prob'])){
		sire<- as.character(ped_mar_col[ped_mar_col$id == row['id'],'sire'])
	}else{#Otherwise, no sire is assigned.
		sire<-NA
	}
	
	return(c(dam, sire))
}

ped[,c("dam", "sire")] <- t(apply(ped, 1, assign_marginal_parents_col))

################################
#                              #
#    Assign parentage using    #
#        COLONY sibships       #
#                              #
################################

#It seems that the COLONY results are internally inconsistant and so it is very difficult to automate the identificiation of siblings. Therefore I will try to make it easier but the user will have to decide which individuals to include as siblings.

# COLONY is inconsistant in the sense that full siblings can be identified above the 0.8 threshold but the best assignment may be that they have different parents.


confident_colfsd<-filter(colfsd, Probability>0.8)
#How to assign parentage based on full sibship:
#Parentage is ONLY assigned via full sibling status if MasterBayes assigns all offspring in the sibship (may be more than 2) the same dam and sire. Note that this includes NAs. The only exception is new_immigrants as MasterBayes does not attempt to assign them parents.
new_immigrants

#If full siblings are identified to share dummy parents, look up the dummy parents' name in colconfig
colconfig

#List of individuals with at least one full sibling from COLONY full sibling dyad file (colfsd).
fs_individuals <- c(as.character(confident_colfsd$OffspringID1), as.character(confident_colfsd$OffspringID2))

#Individuals with more than one sibling.
print("When checking that MasterBayes agrees that these individuals are full siblings must check all siblings.")
names(table(fs_individuals)[table(fs_individuals)>1])

#Use this line of code to check the master bayes consensus pedigree for identified full siblings by adding the appropriate sibling names. Note that new_immigrants will not appear as they are not in the MasterBayes pedigree.
filter(MBped, id%in%c(sibling1, sibling2, sibling3...))

#Once the full siblings have been validated and their parentage correctly set, individuals with uninformative dummy parents should be removed. Where dummy parents have been identified through full siblings the parent name should have been changed to a dummy name from colconfig. (Currently this will have to be done by hand).
#This means that those with "#" or "*" only for a parent are uniformative and so we set them to NA.

ped$dam[which(ped$dam == "#")]<-NA
ped$sire[which(ped$sire == "*")]<-NA

########################
#                      #
#     Merge old and    #
#     new pedigrees    #
#                      #
########################

new_ped_name <- paste("new_pedigree_only_", Sys.Date(), ".csv", sep="")

#write.csv(ped, new_ped_name, row.names=F)

combined_ped <- rbind(previous_ped, ped[, c("id", "dam", "sire")])

combined_ped_name <- paste("extended_pedigree_", Sys.Date(), ".csv", sep="")

#write.csv(combined_ped, combined_ped_name, row.names=F)

###############################
#                             #
#        Compare ped to       #
#  a true simulated pedigree  #
#                             #
###############################

#This is only applicable if the pedigree you have made is based on simulated data generated using MoPedS_simulate_pedigree.R

#load the true simulated pedigree
simped <- readRDS("true_simulated_pedigree.RDS")
true_ped<-as.data.frame(simped$ped)
names(true_ped)<-c("id","true_dam","true_sire")

ped_check <- merge(ped[,1:3],true_ped, on="id")
summary(ped_check$dam == ped_check$true_dam)
ped_check[-which(ped_check$dam == ped_check$true_dam),]

summary(ped_check$sire == ped_check$true_sire)
ped_check[-which(ped_check$sire == ped_check$true_sire),]

#####################
#                   #
#       Checks      #
#                   #
#####################

#Check no duplicates in the colony pedigrees.

joint_duplicates<-sum(duplicated(ped_joint_col$id))
marginal_duplicates<-sum(duplicated(ped_mar_col$id))

if(joint_duplicates>0){
	warning(paste(joint_duplicates, " offspring duplicated in the colony joint pedigree 'ped_joint_col'. Returned offspring may not be the most probable."))
}
if(marginal_duplicates>0){
	warning(paste(marginal_duplicates, " offspring duplicated in the colony marginal pedigree 'ped_mar_col'. Returned offspring may not be the most probable."))
}

#Have any indivdiuals been included multiple times in the pedigree?
duplicated_id<-combined_ped$id[which(duplicated(combined_ped$id))]
if(length(duplicated_id) > 0){
	print("Some individuals have been assigned parents more than once in the combined pedigree!")
	print(duplicated_id)
}

#colconfig should be consistant with the definition of full siblings.
#Add full sibling indexes to colconfig
colconfig$fsindex<-fsindex$index[match(colconfig$OffspringID, fsindex$id)]

colconfig$fs_iprob<-fsindex$iprob[match(colconfig$OffspringID, fsindex$id)]

#
col_fs_check <- summarise(group_by(filter(colconfig, fs_iprob>threshold), fsindex),
	fs_iprob = unique(fs_iprob),
	#Number of unique MasterBayes assigned parent pairs.
	n_MB_pp = length(unique(paste(MotherID,FatherID))))

col_fs_check$consistant <- col_fs_check$n_MB_pp == 1

internally_inconsistant_col_fs <- col_fs_check[!col_fs_check$consistant,]

if(nrow(internally_inconsistant_col_fs) > 0){
	warning("Confident sibships from 'colfsf' are not consistent with 'colconfig'.")
	print(internally_inconsistant_col_fs)
}

#No individual should be their own parent
if(length(which(ped$id == ped$dam)) > 0){
	print("Some individuals in the new pedigree are their own dam!")
	print(ped[which(ped$id == ped$dam),c("id","dam")])
}

if(length(which(ped$id == ped$sire)) > 0){
	print("Some individuals in the new pedigree are their own sire!")
	print(ped[which(ped$id == ped$dam),c("id","sire")])
}



#####################
#                   #
#       To do       #
#                   #
#####################

#Perhaps this plot would be more informative if sorted by joint prob?
plot(rep(1:181,3),c(jointprob,marprobd,marprobs),col=rep(1:3,each=181))



#Check that assigned parents were not dead, i.e. they were in the candidate paretns list. They should be automatically due to the exclusion list but i have seen problems.R

#Does COLONY assign full sibs amoungst parents?
#Answer: Yes, because a cut off threshold is not passed to COLONY confident full siblings (according to colfsd) are not necessarily given the same parents in colconfig.

#Use simulations to see if using the probability of joint assignment improves the pedigree construction reliably.

#Check that all accepted full sibs have the same parents.

#Deal with Unsampled parents in the MB pedigree

#Add in the unsampled parents where appropriate.

#Final pedigree should record assignment type as well as that particular probability

#Include the probabilty threshold for full sibships

#Why do not all pups appear in the colony file? Is it because they weren't included to speed up assignments?

#What happens when the colony assignment is NA or #/*
#Answer: if the dummy parent is useful, ie helps identifiy full siblings it should be changed by hand to the more specific dummy name in colconfig. If not it is automatically changed to NA.

#Check how well the colony parent functions work when there is a range of parent types, at time of writing almost all parents are accurately assigned by MasterBayes.

#Possibly change ped so that it contains all of the accepted assignments but in a separate column for each assignment type. Then an apply function could accept the highest priority assignment type MB_joint>MB_marginal>COL_joint>COL_marginal>COL_sibship.

#Dummy parents should only be assigned where inferred through full siblings, otherwise it is not informative or is based on limited information.