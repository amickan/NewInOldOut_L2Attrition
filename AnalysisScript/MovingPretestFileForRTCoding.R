### Moving files for Experiment 3

A <- c(601:604, 608, 610:620, 622:631)

### for each participant get subset of pretest trials that need RT coding
for (i in 1:length(A)){
  pNumber <- A[i]
  new.parent <- paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/BACK-UP/", pNumber,"/",pNumber,"_Pretest/", sep="")
  new.finaltest <- paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/BACK-UP/", pNumber,"/",pNumber,"_FinaltestEng/", sep="")
  new.folder.org <- paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/CODING/English_Pretest/", pNumber,"_Pretest_subset/", sep="")
  
  if (dir.exists(new.folder.org)){}else{
    dir.create(new.folder.org)}
  
    # find the files that you want
    setwd(new.finaltest)
    infile1 <- paste(pNumber,"Finaltest.txt",sep="_")
    finaltestfile <- as.data.frame(read.delim(infile1, sep = "\t", header = T, skipNul = TRUE))
    
    setwd(new.parent)
    infile1 <- paste(pNumber,"Pretest.txt",sep="_")
    pretest <- as.data.frame(read.delim(infile1, sep = "\t", header = T, skipNul = TRUE))
    files <- as.data.frame(list.files(new.parent, full.names=F, pattern = ".wav"))
    colnames(files) <- "FileName"
  
    for (k in 1:nrow(finaltestfile)){
      num <- which(tolower(as.character(pretest$English_Label)) == finaltestfile$Item[k])
      if (length(num) ==0){
        print(paste("Word did not match:", as.character(finaltestfile$Item[k])))} else {
          trial <- pretest$Trial_nr[num]
          ppfiles <- gsub("PP\\d+.?_Pretest_Trial(\\d+)-001.wav","\\1", files$FileName)
          ind <- which(as.character(ppfiles) == trial)
          file.copy(files$FileName[ind], new.folder.org, overwrite = T, recursive = T, copy.mode = T)
        }
      }
  print(pNumber)
 }

