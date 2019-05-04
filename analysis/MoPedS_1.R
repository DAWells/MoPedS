#How to prepare the microsoft access file for R

#In access, under the external data tab, export the life history file to excel.
#Open the file in excel, check that the very early litter codes have not been converted to dates (this can happen if excel doesnt know they are text).
#Save as .csv

#After the above has been done to create a .csv file execute the following code:
library(dplyr)

#Set your own working directory
setwd("/Users/davidwells/Dropbox/Mongoose/Pedigree/MoPedS/data")

#read in the new life history data
lhdata<-read.csv("rawdata/new_life_history_march_2019_not_r_ready.csv")

#convert the date (currently a factor) to date format
lhdata$DATE<-as.Date(as.character(lhdata$DATE), format="%d/%m/%Y")

#create daten, the date in numeric form
lhdata$daten<- as.numeric(difftime(lhdata$DATE, as.Date("1899-12-30")))


#select the columns you need, I have changed them to lower case because it is easier to type and it fits with my existing code
lhdata<-dplyr::select(lhdata,date=DATE,daten=daten,pack=PACK,indiv=INDIV,sex=SEX,agecategory=AGE.CAT,status=STATUS,stend=START.END,code=CODE,exact=EXACT,lseen=LSEEN,cause=CAUSE, litter=LITTER)

#save lhdata, now reday to use
#write.csv(lhdata,"rawdata/lhdata_march_2019.csv",row.names=F)