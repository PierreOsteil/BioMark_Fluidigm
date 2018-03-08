#Pierre OSTEIL
#v1.0
#2018-03-02

#This code will allow user to read files from FLuidigm software and generate normalised data (delta Ct)


### PACKAGES #############################################################################################################################
library("HTqPCR")

### FUNCTIONS ############################################################################################################################

#################################
### filtering bad quality data ##
#################################

## The column ``Call'' in the sample file contains information about the result of the qPCR reaction.
## Per default, a call of ``Pass'' is translated into ``OK'' in the  \Rcode{featureCategory},
## and ``Fail'' as ``Undetermined''

filter_data <- function(raw){
  
  ## filter empty and dDris annotated samples
  ## set entire column to NA
  exprs(raw)[,c(grep("H2O", colnames(exprs(raw))),grep("-T", colnames(exprs(raw))),grep("-D", colnames(exprs(raw))),grep("water", colnames(exprs(raw))),grep("EMPTY", colnames(exprs(raw))),grep("Empty", colnames(exprs(raw))),grep("-R", colnames(exprs(raw))))]  <- NA
  
  ## ctCall = Failed set to NA
  fail_ind <- which(flag(raw) == "Undetermined")
  exprs(raw)[fail_ind] <- NA
  
  ## remove the RTC and PPC controls and housekeeping genes
  rm_controls <- c("RTC", "PPC", "MGDC")
  
  ##("Actb", "B2m", "Gapdh", "Gusb","Hsp90ab1")
  rm_controls_ind <- which(rownames(exprs(raw)) %in% rm_controls)
  filtered <-  exprs(raw)[-rm_controls_ind,]
  
  ## remove column or row if entire is NA
  filtered_ind_col <- which(colSums(is.na(filtered)) == nrow(filtered))
  if(length(filtered_ind_col) > 0){
    filtered_col <- filtered[,-filtered_ind_col]
  }else{
    filtered_col <- filtered
  }
  
  filtered_ind_row <- which(rowSums(is.na(filtered_col)) == ncol(filtered_col))
  if(length(filtered_ind_row) > 0){
    filtered_row <- filtered_col[-filtered_ind_row,]
  } else{
    filtered_row <- filtered_col
  }
  return(filtered_row)
}



##########################
## Calculate expression ##
##########################

calculate_expression <- function(filtered, LOD){
  exp <- LOD - filtered
  exp[which(exp <= 0)] <- 0.000001
  return(exp)
}

####################
## Remove missing ##
####################

## based on the above assessment have set 40 NA for samples and 23 for genes

remove_missing <- function(exp){
  
  missing <- list()
  
  fail.sample <- which(colSums(is.na(exp))  >= 75)
  if(length(fail.sample) != 0){
    filt.sample <- exp[,-fail.sample]
  }else{
    filt.sample <- exp
  }
  
  ## filter failed genes
  fail.gene <- which(rowSums(is.na(exp))  >= 75)
  if(length(fail.gene) != 0){
    filt.sample2 <- filt.sample[-fail.gene,]
  }else{
    filt.sample2 <- filt.sample
  }
  
  missing <- filt.sample2
  return(missing)
}

####################
## Remove HKGenes ##
####################

remove_HK <- function(exp){
  
  HK <- list()
  
  fail.HK <- which(apply(exp,1,sd, na.rm=T)>= 2) #HK genes SD > 2 CT will be removed
  if(length(fail.HK) != 0){
    filt.HK <- exp[-fail.HK,]
  }else{
    filt.HK <- exp
  }
  
  HK <- filt.HK
  return(HK)
}





## Read in Data to Normalisation ##############################################################################################################################"


## Before to start the analyse:
##  - remove the commas on the two empty rows
##  - save


path <-"C:/Users/Pierre/Desktop/..."
FILE  <- ".csv"

samples_path <- paste(path, FILE , sep="/")
data <- read.csv(samples_path, sep=",", header=F)

### !! NOTE !! ###
## assumes sheets are sorted according to  ID or well name
## need to load in sample sheet with names
## but the easiest way is to load the samples/genes in Fluidigm software


## remove white spaces in sample names
data[,2] <- gsub(" ", "", data[,2], fixed = TRUE)
sample_data <- as.vector(unique((data[-c(1:10),2]))) #extract samples name

## read in CT data
raw_data <- readCtData(files = FILE ,path = path, format="BioMark", n.features= 96, n.data = 96, 
                       samples = sample_data)

dim(raw_data)
rownames(exprs(raw_data))
colnames(exprs(raw_data))

#filtering data
filtered_data <- filter_data(raw_data) #remove controls, water and other things(see function for more details)

#LOD filtering
exp_data <- calculate_expression(filtered_data, 28) #LOD is set 28 but optimal is 24

#remove missing
data_na_filt <- remove_missing(exp_data)

#Normalisation
data_HK  <- data_na_filt[c((nrow(data_na_filt)-4):nrow(data_na_filt)),] # The last 5 genes are the HK genes but not always
data_noHK  <- data_na_filt[c(1:(nrow(data_na_filt)-5)),] #remove HK genes from dataset
sd_HK <- apply(data_HK,1,sd, na.rm=T)
sd_HK
data_HK_filt <- remove_HK(data_HK) #remove HK with high SD. Can be changed in function (see above)

# Average the HK genes
data_HK_mean  <- apply (data_HK_filt, 2, mean, na.rm = T)
sd(data_HK_mean)

# Normalisation
data_hk_norm <- sweep(data_noHK, 2, data_HK_mean) # substarct HK average Ct value to GOI values

#delta CT
data_dCT  <- 2^(data_hk_norm)

#save it!
write.csv(data_dCT, "dCT.csv")



