library(tidyr)
require(reshape)
require(data.table)

A = c(601:604, 608, 610:620, 622:631)

for (i in 1:length(A)){
  pNumber = A[i]
  setwd("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/CODING")
  infile1 <- paste(pNumber,"logfile_manual.txt",sep="_")
  
  data <- read.delim(infile1, header = F)
  data <- separate(data = data, col = V4, into = c("Trial", "rand"), sep = "-")
  data <- separate(data = data, col = Trial, into = c("Trial", "rand2"), sep = "l")
  as.numeric(data$rand2)->data$rand2
  data <- data[order(data$rand2),]
  
  if (length(data$V1) > 46) {
    print(pNumber)
  }
  
  wd1 <-  paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/CODING/English_Finaltest/", pNumber, "_FinaltestEng", sep="")
  setwd(wd1)
  infile2 <- paste(pNumber,"Finaltest.txt",sep="_")
  currentFile <- as.data.frame(read.delim(infile2, stringsAsFactors=FALSE, sep = "\t", header = T, skipNul = TRUE))
  
  for (j in 1:nrow(currentFile)) {
    pos <- which(tolower(as.character(data$rand2)) == tolower(as.character(currentFile$Trial_nr[j])))
    currentFile$RT_new[j] <- data$V5[pos]
  }
  
  for (j in 1:nrow(currentFile)) {
    if(currentFile$RT_new[j]==0){
      currentFile$RT_new[j] <- NA
    }
  }
  
  setwd(wd1)
  write.table(currentFile, infile2, quote = F, row.names = F, col.names = T, sep = "\t")
}

### Read in all data and merge into one data file for Cami and Panthea
data_list <- list()

for (i in 1:length(A)){
  pNumber = A[i]
  wd1 <-  paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/CODING/English_Finaltest/", pNumber, "_FinaltestEng", sep="")
  wd2 <-  paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/CODING/Spanish_Posttest/", pNumber, "_PosttestSpa", sep="")
  infile1 <- paste(pNumber,"Posttest.txt",sep="_")
  infile2 <- paste(pNumber,"Finaltest.txt",sep="_")

  setwd(wd2)
  currentFile <- as.data.frame(read.delim(infile1, stringsAsFactors=FALSE, sep = "\t", header = T, skipNul = TRUE))
  if (length(currentFile[currentFile$Error == 999,]$Error) > 0){
    currentFile[currentFile$Error == 999,]$Error<-1
  }
  
  setwd(wd1)
  currentFile2 <- as.data.frame(read.delim(infile2, stringsAsFactors=FALSE, sep = "\t", header = T, skipNul = TRUE))
  
  ## marking unlearned words as missing values in posttest ##
  for (j in 1:nrow(currentFile)) {
    pos <- which(tolower(as.character(currentFile2$Item )) == tolower(as.character(currentFile$Item[j])))
    if (currentFile$Error[j] == 1) {
      currentFile2$Error[pos] <- NA
      currentFile2$ErrorDetail[pos] <- NA
      currentFile2$VoiceOnset[pos] <- NA
      currentFile2$RT_new[pos] <- NA
  }}
  
  if (length(currentFile2[ifelse(is.na(currentFile2$Error),
                                 1,currentFile2$Error) == 999,]$Error) > 0) {
    currentFile2[ifelse(is.na(currentFile2$Error),
                        1,currentFile2$Error) == 999,]$Error<-1
  }
    
  if (length(currentFile2[ifelse(is.na(currentFile2$Error),
                                 1,currentFile2$Error) == 1,]$VoiceOnset) > 0) {
    currentFile2[ifelse(is.na(currentFile2$Error),
                        1,currentFile2$Error) == 1,]$VoiceOnset <- NA # this excludes words that were produced with errors after interference from RT analysis
  }
  
  if (length(currentFile2[ifelse(is.na(currentFile2$Error),
                                 1,currentFile2$Error) == 1,]$RT_new) > 0) {
    currentFile2[ifelse(is.na(currentFile2$Error),
                        1,currentFile2$Error) == 1,]$RT_new <- NA # this excludes words that were produced with errors after interference from RT analysis
  }
  
  data_list[[i]] <- currentFile2
  
  print(A[i])
  
}

### 07.12.2018 adding pretest corrected RTs and information on article use for pre and posttest
data_list <- list()
## add the article use information to the post dataframe
setwd("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/CODING")
pretest <- as.data.frame(read.delim("PretestArticles.txt", stringsAsFactors=FALSE, sep = "\t", header = T, skipNul = TRUE))
posttest <- as.data.frame(read.delim("FinaltestArticles.txt", stringsAsFactors=FALSE, sep = "\t", header = T, skipNul = TRUE))

for (i in 1:length(A)){
  pNumber = A[i]
  wd1 <-  paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/CODING/English_Finaltest/", pNumber, "_FinaltestEng", sep="")
  wd2 <-  paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/CODING/Spanish_Posttest/", pNumber, "_PosttestSpa", sep="")
  wd3 <-  paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/BACK-UP/", pNumber, "/",pNumber,"_Pretest", sep="")
  wd4 <- "//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/CODING/English_Pretest/"
  infile1 <- paste(pNumber,"Posttest.txt",sep="_")
  infile2 <- paste(pNumber,"Finaltest.txt",sep="_")
  infile3 <- paste(pNumber,"Pretest.txt",sep="_")
  infile4 <- paste(pNumber,"Pretest_logfile_manual.txt",sep="_")
  
  setwd(wd3)
  curPretest <- as.data.frame(read.delim(infile3, stringsAsFactors=FALSE, sep = "\t", header = T, skipNul = TRUE))
  setwd(wd4)
  curNaming <- as.data.frame(read.delim(infile4, stringsAsFactors=FALSE, sep = "\t", header = F))
  
  data <- separate(data = curNaming, col = V4, into = c("Trial", "rand"), sep = "-")
  data <- separate(data = data, col = Trial, into = c("Trial", "rand2"), sep = "l")
  as.numeric(data$rand2)->data$rand2
  data <- data[order(data$rand2),]
  
  for (j in 1:nrow(curPretest)) {
    pos <- which(tolower(as.character(data$rand2)) == tolower(as.character(curPretest$Trial_nr[j])))
    if (length(pos)!=0){
      curPretest$RT_new[j] <- data$V5[pos]} else {
      curPretest$RT_new[j] <- NA
    }
  }
  
  for (l in 1:nrow(curPretest)) {
    if(is.na(curPretest$RT_new[l])==0 && curPretest$RT_new[l] == 0){
      curPretest$RT_new[l] <- NA
    }
  }
  
  curPretest$ArticlesPre <- NA
  if (any(pretest$PP %in% curPretest$Subject_nr[1]) == T){
    for (f in 1:nrow(pretest[pretest$PP==curPretest$Subject_nr[1],])){
      num <- which(tolower(as.character(curPretest$Trial_nr)) == tolower(as.character(pretest$trial[f])))
      curPretest$ArticlesPre[num] <- 1
    }}
  
  curPretest$ArticlesPost <- NA
  if (any(posttest$PP %in% curPretest$Subject_nr[1]) == T){
    for (m in 1:nrow(posttest[posttest$PP==curPretest$Subject_nr[1],])){
      num <- which(tolower(as.character(curPretest$Trial_nr)) == tolower(as.character(posttest$trial[m])))
      curPretest$ArticlesPost[num] <- 1
    }}
  
  setwd(wd2)
  currentFile <- as.data.frame(read.delim(infile1, stringsAsFactors=FALSE, sep = "\t", header = T, skipNul = TRUE))
  if (length(currentFile[currentFile$Error == 999,]$Error) > 0){
    currentFile[currentFile$Error == 999,]$Error<-1
  }
  
  setwd(wd1)
  currentFile2 <- as.data.frame(read.delim(infile2, stringsAsFactors=FALSE, sep = "\t", header = T, skipNul = TRUE))
  
  ## marking unlearned words as missing values in posttest ##
  for (j in 1:nrow(currentFile)) {
    pos <- which(tolower(as.character(currentFile2$Item )) == tolower(as.character(currentFile$Item[j])))
    if (currentFile$Error[j] == 1) {
      currentFile2$Error[pos] <- NA
      currentFile2$ErrorDetail[pos] <- NA
      currentFile2$VoiceOnset[pos] <- NA
      currentFile2$RT_new[pos] <- NA
    }}
  
  if (length(currentFile2[ifelse(is.na(currentFile2$Error),
                                 1,currentFile2$Error) == 999,]$Error) > 0) {
    currentFile2[ifelse(is.na(currentFile2$Error),
                        1,currentFile2$Error) == 999,]$Error<-1
  }
  
  if (length(currentFile2[ifelse(is.na(currentFile2$Error),
                                 1,currentFile2$Error) == 1,]$VoiceOnset) > 0) {
    currentFile2[ifelse(is.na(currentFile2$Error),
                        1,currentFile2$Error) == 1,]$VoiceOnset <- NA # this excludes words that were produced with errors after interference from RT analysis
  }
  
  if (length(currentFile2[ifelse(is.na(currentFile2$Error),
                                 1,currentFile2$Error) == 1,]$RT_new) > 0) {
    currentFile2[ifelse(is.na(currentFile2$Error),
                        1,currentFile2$Error) == 1,]$RT_new <- NA # this excludes words that were produced with errors after interference from RT analysis
  }
  
  
  for (n in 1:nrow(currentFile2)){
    pos <- which(tolower(as.character(curPretest$English_Label)) == tolower(as.character(currentFile2$Item[n])))
    if (length(pos) != 0) {
      currentFile2$RT_pre[n] <- curPretest$RT_new[pos] 
      currentFile2$ArticlesPre[n] <- curPretest$ArticlesPre[pos]
      currentFile2$ArticlesPost[n] <- curPretest$ArticlesPost[pos]
    } else {
      currentFile2$RT_pre[n] <- NA
      currentFile2$ArticlesPre[n] <- NA
      currentFile2$ArticlesPost[n] <- NA
    }
  }
  
  data_list[[i]] <- currentFile2
  
  wd5 <-  paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/CODING/English_Pretest/", pNumber, "_Pretest_subset/", sep="")
  setwd(wd5)
  write.table(curPretest, infile3, quote = F, row.names = F, col.names = T, sep = "\t")
  
  print(A[i])
  
}

post <- rbindlist(data_list)
setwd("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan")
write.table(post, "CompleteDataset_Honours.txt", quote = F, row.names = F, col.names = T, sep = "\t")
