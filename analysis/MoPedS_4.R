#Generate the input data for COLONY2. Uses the PdataPed object generated in MoPedS_3.

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


#Alphabetically reorder the rows of error and columns of gendata so that they match.
gendata<-gendata[,c(1,(1+order(names(gendata)[2:87])))]
error<-error[order(error$locus),]

print("Please check that the names of loci in gendata and error correctly match up")
cbind(names(gendata)[2:87], as.character(rep(error$locus,each=2)))

gendataNA<-gendata
gendataNA[gendataNA==0]<-NA
####################################
#                                  #
#      Generate the potential      #
#           parent files           #
#                                  #
####################################

#Create a dataframe of potential parents. Colony takes an list of the excluded parents but this is used to generate the excluded parents and may be useful for checks.

#load the MasterBayes PdataPed object we created in MoPedS_3.R which specifies the potential parents.
PdP<-readRDS("PdataPed_object.RDS")

#Xlist is the design matrix used by MasterBayes and it includes the potential parents of each offspring. These parent IDs are stored as integers which refer to an ID in Xlist$id. Integers greater than the length of Xlist$id return NA and indicates an unsampled parent.
Xlist<-getXlist(PdP)

#Get all the pup ids. The name of each item in X is an integer refering to an individual in Xlist$id, and within that item is all the potential parents of the individual.
pup_id<-Xlist$id[as.numeric(names(Xlist$X))]

#Add to this list of pup id the new_immigrants because COLONY takes them as offspring to assign sibships.
new_immigrants<-as.character(filter(pdata, new_immigrant==T)$id)

pup_id<-c(pup_id, new_immigrants)

#Function to get the id of potential parents
get_potential_parents<-function(offspring,parent){
	dams<-Xlist$id[offspring$restdam.id]
	sires<-Xlist$id[offspring$restsire.id]
	
	if(parent=="dam"){
		return(dams)
	}else if(parent =="sire"){
		return(sires)
	}else{
		stop('parent must be "dam" or "sire".')
	}
}

#################
#Potential DAMS
#################

pdams_list<-lapply(Xlist$X,get_potential_parents,parent="dam")

#Extend pdams_list with NAs by the length of new_immigrants .
pdams_list <- c(pdams_list, as.list(rep(NA, length(new_immigrants))))

#Make all lists of potential dams the same length. To do this we get the highest number of potential dams.
max_dams<-max(unlist(lapply(pdams_list,length)))
#Function to extend each list.
extend_list_dam<-function(list){
	length(list)<-max_dams
	return(list)
}
#Pass this function to lapply
pdams_evenlist<-lapply(pdams_list,extend_list_dam)

#Convert this to a data.frame
pdams_df<-do.call(rbind.data.frame,pdams_evenlist)
names(pdams_df)<-NULL

#Create a dataframe with pup/immigrant id in the first column and subsequent columns filled with potential dams.
pdams<-cbind(pup_id,pdams_df)


#################
#Potential SIRES
#################

psires_list<-lapply(Xlist$X,get_potential_parents,parent="sire")

#Extend psires_list with NAs by the length of new_immigrants .
psires_list <- c(psires_list, as.list(rep(NA, length(new_immigrants))))

#Make all lists of potential sires the same length. To do this we get the highest number of potential sires.
max_sires<-max(unlist(lapply(psires_list,length)))

extend_list_sire<-function(list){
	length(list)<-max_sires
	return(list)
}

#Use the above defined extend_list function.
psires_evenlist<-lapply(psires_list,extend_list_sire)
#Convert this to a data.frame
psires_df<-do.call(rbind.data.frame,psires_evenlist)
names(psires_df)<-NULL

#Add the pup's id.
psires<-cbind(pup_id,psires_df)

############################
#                          #
#       Project name       #
#                          #
############################
project_name<-readline(prompt="Enter a project name, less than 40 characters. NO spaces.")
colony_test
output_name<-readline(prompt="Enter an output file name, less than 40 characters. NO spaces. All colony output files will have this name but different extensions.")
colony_test_file

project_name_string<-paste(project_name, "! Row 1) Project name.")
output_name_string<-paste(output_name, "! Row 2) Output file name.")

############################
#                          #
#    Number of samples     #
#                          #
############################
number_of_loci<-0.5*(ncol(gendata)-1)
random_seed<-round(runif(1,0,9999)) #If the data are informative changing the seed should not affect assignment but using the same seed and data should give the same results.
number_of_loci_string<-paste(number_of_loci, "! Row 4) Number of loci.")
random_seed_string<-paste(random_seed, "! Row 5) Random seed.")

#################################
#                               #
#     Prespecified settings     #
#                               #
#################################
#These are detailed in the colony2 user guide

inputs_6to14_file<-"../analysis/colony/colony_inputs_6-14.txt"
inputs_6to14<-readChar(inputs_6to14_file,file.info(inputs_6to14_file)$size)

inputs_17to23_file<-"../analysis/colony/colony_inputs_17-23.txt"
inputs_17to23<-readChar(inputs_17to23_file,file.info(inputs_17to23_file)$size)

inputs_33to39_file<-"../analysis/colony/colony_inputs_33-39.txt"
inputs_33to39<-readChar(inputs_33to39_file,file.info(inputs_33to39_file)$size)

############################
#                          #
#    Allele frequencies    #
#                          #
############################

#Calculate Allele Frequencies.
AF<-extractA(gendataNA)
#Calculate the number of alleles for each locus
number_of_alleles<-sapply(AF,length)
names(number_of_alleles)<-NULL

number_of_alleles_string<-paste(paste(number_of_alleles, collapse=", "), "! Row 15) Number of alleles per locus.")

#Loci names as comments, "!" is the comment character for colony2.dat, \n creates a new line.
loci_names_comment<-paste("!", names(AF), sep="")

#For each locus, output a row listing the alleles and a row listing their frequencies.
convert_AF_to_string<-function(item1, item2){
	locus<-item2
	alleles<-paste(names(item1), collapse=", ")
	freq<-paste(item1, collapse=", ")
	
	return(c(alleles,freq))
}

#Get a separate string for each locus' alleles and frequencies.
AF_string_vector<-mapply(convert_AF_to_string, AF, loci_names_comment)

############################
#                          #
#      Locus specific      #
#       error rates        #
#                          #
############################

#The names of the loci
marker_names<-error$locus
marker_names_string<-paste(paste(marker_names, collapse=", "), "! Row 24) Marker names.")

#Indicate codominant markers with a 0.
marker_type<-rep(0,length(marker_names))
marker_type_string<-paste(paste(marker_type,collapse=", "), "! Row 25) Marker type. 0/1 = codominant/dominant")

dropout_rate<-error$E1
dropout_rate_string<-paste(paste(dropout_rate, collapse=", "), "! Row 26) dropout rate, E1.")

other_errors<-error$E2
other_errors_string<-paste(paste(other_errors, collapse=", "), "! Row 27) Other errors, E2.")


###############################
#                             #
#     Genotype duplicates     #
#                             #
###############################

#Unlike MasterBayes, COLONY does not allow multiple genotypes of an individual. Create a function to remove duplicated genotypes. How should this be accomplished? Should any mismatches be set to NA? or should heterozygous/homozygous genotypes imply drop out in one? Fill in missing loci and set mismatches to missing?

#Get the IDs of individuals with multiple genotypes.
duplicated_gen_id<-unique(gendata$id[duplicated(gendata$id)])
#Use the IDs to get the duplicate genotypes.
duplicated_gen<-gendata[gendata$id %in% duplicated_gen_id,]
#Order by id so duplicates appear together.
duplicated_gen<-duplicated_gen[order(duplicated_gen$id),]

#NOTE: The generate_consensus_genotype function assumes that the two alleles listed for a locus are listed small, large. This is the case for gene mapper and genemarker.

#Check allele order
allele1<-duplicated_gen[,seq(2, ncol(duplicated_gen), 2)]
allele2<-duplicated_gen[,seq(3, ncol(duplicated_gen), 2)]

if(0 < sum(!(allele1<=allele2))){
	print("The alleles are not ordered small to large within the loci. Therefore the function generate_consensus_genotypes will not function correctly. Please ensure that each individual has only one genotype record for COLONY by creating a new gendata file to feed to MoPedS_4.R to create a COLONY input file.")
}

#Given multiple genotypes of a single individual this function combines the information into a single genotype. Specifically, missing loci are filled from other genotyping attempts, and conflicting loci are labelled as missing.
generate_consensus_genotype<-function(allele_set){
	#Get all the non-missing alleles at one of two spots for a locus.
	non_missing_set<-unique(allele_set[allele_set != 0])
	#If there is only 1 non-missing allele this is inferred.
	if(length(non_missing_set) == 1){
		consensus_allele<-rep(non_missing_set, length(allele_set))
	}else if(length(non_missing_set) != 1){
		#If there is more or less than 1 non-missing allele at this position it is set to missing.
		consensus_allele<-rep(0,length(allele_set))
	}
	
	return(consensus_allele)
}


get_and_merge_duplicate<-function(indiv){
	#Get all duplicate genotypes for an individual
	duplicated_indiv<-gendata[gendata$id == indiv,]
	#generate consensus genotype from these duplicates	
	consensus_gen<-apply(duplicated_indiv[,2:ncol(gendata)], 2, generate_consensus_genotype)[1,]
	
	return(consensus_gen)
}

#This line generates the consensus genotypes for all duplicated individuals.
consensus_genotypes<-t(sapply(duplicated_gen_id,get_and_merge_duplicate))
rownames(consensus_genotypes)<-duplicated_gen_id

#Remove the duplicated individuals from gendata and update the remaining individuals with the consensus genotype.
gendata<-gendata[!duplicated(gendata$id),]
gendata[match(rownames(consensus_genotypes), gendata$id),2:ncol(gendata)]<-consensus_genotypes


##############################
#                            #
#  Offspring and candidates  #
#                            #
##############################

#We use COLONY to assign parentage to offspring and immigrants, in contrast MasterBayes only assigns parentage to offspring. Therefore this "offspring" dataframe must also include the new_immigrants that we wish to assign to sibship groups.
#new_immigrants must also be added to pup_id.

offspring<-filter(pdata,offspring==1 | new_immigrant==T)
offspring_gen<-filter(gendata, id %in% offspring$id)
offspring_gen_string<-apply(offspring_gen, 1, paste, collapse=", ")

number_of_offspring<-nrow(offspring_gen)
#number_of_offspring<-length(Xlist$X)
number_of_offspring_string<-paste(number_of_offspring, "! Row 3) Number of offspring.")

candidate_sires<-filter(pdata,sex=="Male" & juv==F)
candidate_sires_string<-apply(candidate_sires, 1, paste, collapse=", ")

candidate_dams<-filter(pdata,sex=="Female" & juv==F)
candidate_dams_string<-apply(candidate_dams, 1, paste, collapse=", ")

#probability that the true father and mother, respectively, are included in the candidate parents.
true_parent_inclusion_prob<-paste(0.8, 0.8, sep=", ")
true_parent_inclusion_prob_string<-paste(true_parent_inclusion_prob,"! Row 29) Probability that the true father and mother are included in the candidate parents, respectively.")

#Sire and dam genotypes
sire_gen<-filter(gendata, id %in% candidate_sires$id)
dam_gen<-filter(gendata, id %in% candidate_dams$id)

sire_gen_string<-apply(sire_gen,1,paste,collapse=", ")
dam_gen_string<-apply(dam_gen,1,paste,collapse=", ")

#Number of candidate sires and dams identified for MasterBayes.
NO_candidate_sires<-length(unique(candidate_sires$id))
NO_candidate_dams<-length(unique(candidate_dams$id))

#The number of candidate sires and dams with genotypes and therefore passed to COLONY
NO_candidate_sires<-nrow(sire_gen)
NO_candidate_dams<-nrow(dam_gen)

number_of_candidate_parents<-paste(NO_candidate_sires, NO_candidate_dams, sep=", ")
number_of_candidate_parents_string<-paste(number_of_candidate_parents,"! Row 30) Number of candidate fathers and mothers provided to colony, respectively.")

##########################
#                        #
#       Exclusions       #
#                        #
##########################

##########
# Dams
##########

#All offspring have some genotyped individuals which cannot be their dam.
number_of_dam_exclusions<-paste(number_of_offspring, "! Row 43) Number of offspring with known maternal exclusions.")

#Function to return the offspring, number of excluded dams (or sires), and list the excluded dams (or sires).
#This function is to be handed to mapply and will be executed on two vectors/lists. Not that this function only excludes individuals with genotypes as possible parents. We may already know that they cannot be the parent without looking at the genetic data but COLONY only uses the genetic data so they may as well be unsampled.

exclude_parents<-function(item1, item2, parent){
	if(parent == "Dam"){
		#How many exclusions
		number_of_exclusions<-NO_candidate_dams-sum(!is.na(item2))
		#Which dams are excluded
		exclusions<-unique(dam_gen$id[!dam_gen$id %in% item2])
		
	} else if(parent == "Sire"){
		#How many exclusions
		number_of_exclusions<-NO_candidate_sires-sum(!is.na(item2))
		#Which sires are excluded
		exclusions<-unique(sire_gen$id[!sire_gen$id %in% item2])
		
	} else{
		stop('parent must be either "Dam" or "Sire".')
	}
	
	exclusions<-as.character(exclusions)
	
	return(c(item1,number_of_exclusions, exclusions))
}

dam_exclusion_list<-mapply(exclude_parents, pup_id, pdams_list, "Dam")

#Convert each pup_id, exclusion number and exclusion list into a single string.
dam_exclusion_vector<-sapply(dam_exclusion_list,paste,collapse=", ")

#Only pass the exclusion list of genotyped individuals to Colony.
dam_exclusion_of_genotyped_offspring<-dam_exclusion_vector[names(dam_exclusion_vector) %in% offspring_gen$id]

################
#    Sires
################

#Each row of exclusions should be the pup it refers to, the number of excluded sires, and then list those excluded dams.

#All offspring have some genotyped individuals which cannot be their sire.
number_of_sire_exclusions<-paste(number_of_offspring, "! Row 41) Number of offspring with known paternal exclusions.")

sire_exclusion_list<-mapply(exclude_parents, pup_id, psires_list, "Sire")
sire_exclusion_vector<-sapply(sire_exclusion_list,paste,collapse=", ")

#Only pass the exclusion list of genotyped individual's to Colony.
sire_exclusion_of_genotyped_offspring<-sire_exclusion_vector[names(sire_exclusion_vector) %in% offspring_gen$id]


#########################
#                       #
#      Write colony     #
#       input file      #
#                       #
#########################
input_text<-c(
project_name_string,
output_name_string,
number_of_offspring_string,
number_of_loci_string,
random_seed_string,
inputs_6to14,
number_of_alleles_string,
#"! row 16) Alleles and their frequency per locus",
AF_string_vector,
inputs_17to23,
marker_names_string,
marker_type_string,
dropout_rate_string,
other_errors_string,
#"! row 28) Offspring genotypes",
offspring_gen_string,
true_parent_inclusion_prob_string,
number_of_candidate_parents_string,
#"! row 31) Candidate male genotypes",
sire_gen_string,
#"! row 32) Candidate female genotypes",
dam_gen_string,
inputs_33to39,
number_of_sire_exclusions,
#"! Row 42) Excluded paternity, first is the offspring then all of the males which cannot be their sire are listed.",
sire_exclusion_of_genotyped_offspring,
number_of_dam_exclusions,
#"! Row 44) Excluded maternity, first is the offspring then all of the females which cannot be their dam are listed.",
dam_exclusion_of_genotyped_offspring,
"0     ! Row 45) Number of offspring with known excluded Paternal sibships",
"0     ! Row 47) Number of offspring with known excluded Maternal sibships")

input_file_name<-paste(project_name, "_colony_inputs.dat",sep="")
writeLines(input_text,input_file_name)

##########################
#                        #
#     Running colony2    #
#                        #
##########################
#COLONY is stand-alone software presented in:
#Jones, O. and Wang, J. (2009) COLONY: a program for parentage and sibship inference from multilocus genotype data. Molecular Ecology Resources 10: 551 - 555.
#and
#Wang, J. (2004) Sibship reconstruction from genetic data with typing errors. Genetics 166: 1963 - 1979.

#COLONY2 can be downloaded from:
#https://www.zsl.org/science/software/colony

#There is also a R version but I am not sure it is complete:
# rcolony, Owen R. Jones <jones@demogr.mpg.de>

#Download and install COLONY2 from the link above
#Copy the .dat file you have just created using MoPedS_4.R into the folder containing Colony2p.exe.

#Open terminal or command depending on your opporating system and navigate to the folder containing Colony2p.exe and the .dat file. This can be done using cd and the absolute folder location. E.g.
# cd "/Applications/colony2.mac"

#Then, in the terminal window run the start command. The exact start comand depends on how you intend to run it, serial or parallel, and your opporating system. The colony user guide gives specific details and the readme.pdf which accompanied your colony download will have the specific command for your installation.

#On mac
#./colony2s.out IFN:colony_test_colony_inputs.dat
#multicore
#mpiexec -n 4 ./colony2p.out IFN:colony_test_colony_inputs.dat


#########################
#                       #
#         Checks        #
#                       #
#########################
#project and file names must not include spaces " ".

#Check that the number of offspring matches the number of offspring genotype records

#same for dams and sires.

#Ensure that alleles are listed small large as output by genemapper


############################
#
#    To do
#
############################
#The loci names cause errors, but that's not all that is not working...

#Missing loci must be NA when calculating allele frequencies

#set checks before the writeLines

#Make sure exclusions deals with duplicates in the candidate parents.

#offspring genotypes, candidate dam and sire genotypes
#Other file info.

#Convert potential parents into excluded parents

#duplicate genotypes are not allowed and will have to be removed. Ideally missing values are filled in.

#Individuals with only 1 or fewer loci have to be removed.

#Not all offspring, in the pdata exist in the microsatellite data. This must be made clear in the final report so we know how many offspring that we attempted to assign and actually did manage to assign.
#same for the dams and sires

#Check for and report individuals with identical genotypes. Possibly using paste to combine genotypes and groupby to identify >1 unique id.

#only 1 genotype entry is allowed for each individual.

#Report the genotype changes made for the final report.

#Ensure the allele order check is unavoidable.

#Count the number of mismatches between duplicates of the same individual to catch "duplicates" which due to misslabelling are infact different individuals.

#Should individuals in the pdata but without genotypes be provided with missing at all loci or should they not be provided? I think it will not let you submit offspring with no loci but what about parents? We could 

#Should probably remove NA from exclusion lists.

#Reading Jenny's report COLONY was used to assign parentage even to non-offspring. Perhaps I should find some way to include them in the colony input file so that immigrants which are full sibs can be identified.

#In order to get relatedness information about immigrants, all newly genotyped individuals should be entered as offspring for COLONY. How to ensure that this works correctly and parent offspring are not identified as full sibs?

#Increase the number of runs? recomended 5, are the results of all runs combined?

#immigrants must also be added to pup_id.