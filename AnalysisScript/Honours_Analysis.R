### Honours project data analysis ###
require(reshape)
require(data.table)
require(here)

A = c(601:603, 608, 610:620, 622, 623, 625:631)
post <- read.table(here("CompleteDataset_Honours.txt"), header = T)
post<-post[!(post$Subject_nr == 604 | post$Subject_nr == 624),]
post <- droplevels(post)

# setting variables
post$Subject_nr <- as.factor(post$Subject_nr)
post$Condition <- as.factor(post$Condition)
post$Errorfact <- as.factor(post$Error)

# checking RTs
min(post$RT_new, na.rm=T)
hist(post$RT_new)
shapiro.test(post$RT_new) ## data are not normal 
# log-transform RTs
post$RT_new_log <- log(post$RT_new)

# checking coding instances with error code "6"
subset <- post[is.na(post$ErrorDetail)==0 && post$ErrorDetail == 6,] # only one case

# new corrected for reaction time
post$RTdiff <- post$RT_pre - post$RT_new
post$Prelog <- log(post$RT_pre)
post$RTdifflog <- post$Prelog - post$RT_new_log

###### How many people used articles on each test #####
table(post$Subject_nr, post$ArticlesPre) # Pps 604, 618, 624 and 630 used a lot of articles, let's try to exclude them and see what happens to the effect
#post<-post[!(post$Subject_nr== 604 | post$Subject_nr == 618 | post$Subject_nr == 624 | post$Subject_nr == 630),]
#post <- droplevels(post)

# exclude trials in which articles were used from RT analysis 
for (i in 1:nrow(post)){
  if (is.na(post$ArticlesPost[i]) == 0 && is.na(post$ArticlesPre[i]) == 0){
  } else if (is.na(post$ArticlesPost[i]) == 0 && is.na(post$ArticlesPre[i]) == 0) {
  } else if (is.na(post$ArticlesPost[i]) == 0 | is.na(post$ArticlesPre[i]) == 0){
    post$RT_new[i] <- NA
    post$RT_pre[i] <- NA
    post$Prelog[i] <- NA
    post$RT_new_log[i] <- NA
    post$RTdiff[i] <- NA
    post$RTdifflog[i] <- NA
  } 
}
# check how many trials per person we have left
table(post[is.na(post$RTdiff)==0,]$Subject_nr)
table(post[is.na(post$RTdiff)==0,]$Subject_nr, post[is.na(post$RTdiff)==0,]$Condition) # per condition

# percentage of article traisl per person per condition
article <- post[(is.na(post$ArticlesPre)==0 | is.na(post$ArticlesPost)==0),]
article1 <- (table(article$Subject_nr, article$Condition)/23)*100
article2 <- (table(article$Subject_nr)/46)*100


### deleting trials of words that were already known in Spanish before the learning phaser
known <- read.delim("KnownWords.txt")
for (i in 1:nrow(known)){
  pNumber <- known$PP[i]
  num <- which(tolower(post[post$Subject_nr == pNumber,]$Spanish_Label) == tolower(known$Word[i]))
  if (length(num)!= 0 ){
    post[post$Subject_nr == pNumber,]$RT_new[num] <- NA
    post[post$Subject_nr == pNumber,]$RT_pre[num] <- NA
    post[post$Subject_nr == pNumber,]$Prelog[num] <- NA
    post[post$Subject_nr == pNumber,]$RT_new_log[num] <- NA
    post[post$Subject_nr == pNumber,]$RTdiff[num] <- NA
    post[post$Subject_nr == pNumber,]$RTdifflog[num] <- NA
    post[post$Subject_nr == pNumber,]$Error[num] <- NA
  }
}

########## Plots with GGplot ###########
require(plyr)
require(ggplot2)

### Fine-grained error rates ###

# histogram of results 
hist(post$Error)

ddply(post, .(Condition, Subject_nr), 
      summarise, N=length(Error), 
      mean   = mean(Error, na.rm = TRUE), 
      sem = sd(Error, na.rm = TRUE)/sqrt(N)) -> aggregatedError

aggregated_means_error <- ddply(post, .(Condition), 
                                summarise,
                                condition_mean = mean(Error,na.rm = T),
                                condition_sem = sd(Error,na.rm = T)/sqrt(length(Error[!is.na(Error)])))

aggregatedError <- merge(aggregatedError, aggregated_means_error, by = c("Condition"))

aggregated_means_error$condition_mean <- aggregated_means_error$condition_mean*100
aggregated_means_error$condition_sem <- aggregated_means_error$condition_sem*100
aggregatedError$mean <- aggregatedError$mean*100
aggregatedError$sem <- aggregatedError$sem*100
aggregatedError$condition_mean <- aggregatedError$condition_mean*100
aggregatedError$condition_sem <- aggregatedError$condition_sem*100


lineplot <- ggplot(aggregatedError, aes(y = mean, x = Condition, group = Subject_nr))
lineplot + geom_point(color="darkgrey") +
  geom_line(color="darkgrey") +
  geom_point(aes(y = condition_mean,
                 color = Condition), color="black") +
  geom_text(aes(label=Subject_nr)) +
  geom_line(aes(y = condition_mean,color="red")) +
  geom_errorbar(aes(ymin=condition_mean-condition_sem,
                    ymax=condition_mean+condition_sem,
                    color = "red",
                    na.rm = T),
                width = 0.5) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20)) + 
  scale_x_discrete(labels=c("Interference", "No interference"), breaks = 1:2, expand = c(0.1,0.1)) +
  ylab("Percentage correctly recalled words in Spanish") +
  # scale_color_manual(guide=F, "Frequency Condition", values=c("dodgerblue4","firebrick"),labels=c("High","Low")) +
  theme_bw()

barplot <- ggplot(aggregated_means_error, aes(y = condition_mean, x = Condition, fill = Condition))
barplot + geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=condition_mean-condition_sem,
                    ymax=condition_mean+condition_sem),
                width = 0.5, position=position_dodge(0.9)) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20)) + 
  coord_cartesian(ylim=c(0,5)) +
  scale_x_discrete(labels=c("Interference", "No Interference"), breaks = 1:2, expand = c(0.1,0.1)) +
  ylab("Percentage incorrectly recalled words in English") +
  scale_fill_grey(labels=c("Interference","No Interference")) +
  theme_bw()

# Accuracy instead of error rates

post$Accuracy <- 1-post$Error
ddply(post, .(Condition, Subject_nr), 
      summarise, N=length(Accuracy), 
      mean   = mean(Accuracy, na.rm = TRUE), 
      sem = sd(Accuracy, na.rm = TRUE)/sqrt(N)) -> aggregatedAcc

aggregated_means_acc <- ddply(post, .(Condition), 
                                summarise,
                                condition_mean = mean(Accuracy,na.rm = T),
                                condition_sem = sd(Accuracy,na.rm = T)/sqrt(length(Accuracy[!is.na(Accuracy)])))

aggregatedAcc <- merge(aggregatedAcc, aggregated_means_acc, by = c("Condition"))

aggregated_means_acc$condition_mean <- aggregated_means_acc$condition_mean*100
aggregated_means_acc$condition_sem <- aggregated_means_acc$condition_sem*100
aggregatedAcc$mean <- aggregatedAcc$mean*100
aggregatedAcc$sem <- aggregatedAcc$sem*100
aggregatedAcc$condition_mean <- aggregatedAcc$condition_mean*100
aggregatedAcc$condition_sem <- aggregatedAcc$condition_sem*100


barplot <- ggplot(aggregated_means_acc, aes(y = condition_mean, x = Condition, fill = Condition))
barplot + geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=condition_mean-condition_sem,
                    ymax=condition_mean+condition_sem),
                width = 0.5, position=position_dodge(0.9)) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20)) + 
  coord_cartesian(ylim=c(90,100)) +
  scale_x_discrete(labels=c("Interference", "No Interference"), breaks = 1:2, expand = c(0.1,0.1)) +
  ylab("Percentage correctly recalled words in English") +
  scale_fill_grey(labels=c("Interference","No Interference")) +
  theme_bw()


#### Plot for RTs ###
ddply(post, .(Condition, Subject_nr), 
      summarise, N=length(RT_new), 
      mean   = mean(RT_new, na.rm = TRUE), 
      sem = sd(RT_new, na.rm = TRUE)/sqrt(N)) -> aggregatedrt

aggregated_means_rt<- ddply(post, .(Condition), 
                                summarise,
                                condition_mean = mean(RT_new,na.rm = T),
                                condition_sem = sd(RT_new,na.rm = T)/sqrt(length(RT_new[!is.na(RT_new)])))

aggregatedrt <- merge(aggregatedrt, aggregated_means_rt, by = c("Condition"))


lineplot <- ggplot(aggregatedrt, aes(y = mean, x = Condition, group = Subject_nr))
lineplot + geom_point(color="darkgrey") +
  geom_line(color="darkgrey") +
  geom_point(aes(y = condition_mean,
                 color = Condition), color="black") +
  geom_text(aes(label=Subject_nr)) +
  geom_line(aes(y = condition_mean,color="red")) +
  geom_errorbar(aes(ymin=condition_mean-condition_sem,
                    ymax=condition_mean+condition_sem,
                    color = "red",
                    na.rm = T),
                width = 0.5) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20)) + 
  scale_x_discrete(labels=c("Interference", "No interference"), breaks = 1:2, expand = c(0.1,0.1)) +
  ylab("Naming latencies in ms") +
  # scale_color_manual(guide=F, "Frequency Condition", values=c("dodgerblue4","firebrick"),labels=c("High","Low")) +
  theme_bw()

barplot <- ggplot(aggregated_means_rt, aes(y = condition_mean, x = Condition, fill = Condition))
barplot + geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=condition_mean-condition_sem,
                    ymax=condition_mean+condition_sem),
                width = 0.5, position=position_dodge(0.9)) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20)) + 
  scale_x_discrete(labels=c("Interference", "No Interference"), breaks = 1:2, expand = c(0.1,0.1)) +
  ylab("Naming latencies in ms") +
  scale_fill_grey(labels=c("Interference","No Interference")) +
  theme_bw()

#### Plot for RTs difference ###
ddply(post, .(Condition, Subject_nr), 
      summarise, N=length(RTdiff), 
      mean   = mean(RTdiff, na.rm = TRUE), 
      sem = sd(RTdiff, na.rm = TRUE)/sqrt(N)) -> aggregatedrtdiff

aggregated_means_rtdiff<- ddply(post, .(Condition), 
                            summarise,
                            condition_mean = mean(RTdiff,na.rm = T),
                            condition_sem = sd(RTdiff,na.rm = T)/sqrt(length(RTdiff[!is.na(RTdiff)])))

aggregatedrtdiff <- merge(aggregatedrtdiff, aggregated_means_rtdiff, by = c("Condition"))


lineplot <- ggplot(aggregatedrtdiff, aes(y = mean, x = Condition, group = Subject_nr))
lineplot + geom_point(color="darkgrey") +
  geom_line(color="darkgrey") +
  geom_point(aes(y = condition_mean,
                 color = Condition), color="black") +
  geom_text(aes(label=Subject_nr)) +
  geom_line(aes(y = condition_mean,color="red")) +
  geom_errorbar(aes(ymin=condition_mean-condition_sem,
                    ymax=condition_mean+condition_sem,
                    color = "red",
                    na.rm = T),
                width = 0.5) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20)) + 
  scale_x_discrete(labels=c("Interference", "No interference"), breaks = 1:2, expand = c(0.1,0.1)) +
  ylab("Naming latencies in ms") +
  # scale_color_manual(guide=F, "Frequency Condition", values=c("dodgerblue4","firebrick"),labels=c("High","Low")) +
  theme_bw()

barplot <- ggplot(aggregated_means_rtdiff, aes(y = condition_mean, x = Condition, fill = Condition))
barplot + geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=condition_mean-condition_sem,
                    ymax=condition_mean+condition_sem),
                width = 0.5, position=position_dodge(0.9)) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20)) + 
  scale_x_discrete(labels=c("Interference", "No Interference"), breaks = 1:2, expand = c(0.1,0.1)) +
  ylab("Speed up in naming latencies from English pre- to posttest (in ms)") +
  scale_fill_grey(labels=c("Interference","No Interference")) +
  theme_bw()



###### Stats on behavioral results ######

require(lme4)
require(lmerTest)
require(lmtest)

# setting contrasts to the mean of each condition 
contrasts(post$Condition) <- c(-0.5,0.5)
# turning my factors into numerical factors reflecting a dummy coding 
post$ConditionN <- (-(as.numeric(post$Condition)-2))-0.5


###### Accuracy after interference #####

## Full model with maximal random effects structure
modelfull <- glmer(Error ~ ConditionN + (1|Item) + (1+ConditionN|Subject_nr), family = "binomial", control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)), data = post)
summary(modelfull)
# the model converges with the maximal justifyable random effects structure, and none of the random effects are highly correlated with each other, so we leave it this complex 
# no comparisons needed, you report the beta weights from this model in a table in your paper

## Simple Anova for accuracy
## Arcsine transformed data Anova
anova <- aov(Error ~ Condition, data = post)
summary(anova)


## Model reporting - fullest model above
# It is best to report Chi-square p-values for each of the effects serpately 
# First let's take out the main effect for Condition (-Condition below in the code)
modelCondition<- glmer(Error ~ ConditionN -ConditionN + (1|Item) + (1+ConditionN|Subject_nr), family = binomial, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)), data = post)
anova(modelfull, modelCondition)
# The chi-suare p-value from the Anova table is the p-value for the main effect of Condition. This p-value is slightly higher than the one from the model output itself because the distribution against which it is calcualted is different (chi-square vs z-distribution)
#IMPORTANT: the intercept in these models is always the grand mean: the effect over all conditions: mean over the mean of each cell. cells being: Interference condition for Block 1, Interference Block 2, No interference Block 1, No interference Block 2
# So now it is not correct anymore what you say in your methods section: the intercept DOES NOT reflect the no interference condition any longer, it represents the mean of both conditions over both blocks!!! 
# The p-values you get out of these comparisons are what you report in the paper and in the table along with the estimates.  

###### Modelling for RTs #####
# simple Anova for RTs (log-transformed)
anova_rt <- aov(RT_new_log ~ Condition, data = post)
summary(anova_rt)

## Full model on log transformed data 
# Full model with maximum random effects structure 
# We take the log of the reaction times because the distribution is very non-normal, and we subtract 2000ms because that's the lowest value there is currently (due to 2s delay), log transform works better if there are values close to 0 and between 0-1
modelRT2full <- lmer(RTdifflog ~ ConditionN + (1|Item) + (1+ConditionN|Subject_nr), control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)),data = post)
summary(modelRT2full)

## Model reporting - fullest model above
# Same as above
# First let's take out the main effect for Condition (-Condition below in the code)
modelRT2Condition <- lmer(RTdifflog ~ ConditionN - ConditionN + (1|Item) + (1+ConditionN|Subject_nr), control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)),data = post)
anova(modelRT2full, modelRT2Condition)
# IMPORTANT: the intercept in these models is always the grand mean: the effect over all conditions: mean over the mean of each cell. cells being: Interference condition for Block 1, Interference Block 2, No interference Block 1, No interference Block 2
# So now it is not correct anymore what you say in your methods section: the intercept DOES NOT reflect the no interference condition any longer, it represents the mean of both conditions over both blocks!!! 
# The p-values you get out of these comparisons are what you report in the paper and in the table along with the estimates. 



### Correlations with LexTale score, interference performance 

#LexTale
library(tidyr)
A = c(601:603, 608, 610:620, 622, 623, 625:631)
data_list <- list()
for (i in 1:length(A)){
  pNumber = A[i]
  wd <- paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/BACK-UP/", pNumber, "/", sep="")
  setwd(wd)
  infile1 <- paste(pNumber,"score_LexTale.txt",sep="_")
  
  currentFile <- as.data.frame(read.delim(infile1, sep = "\t", header = T, skipNul = TRUE))
  
  data_list[[i]] <- currentFile
  
  print(A[i])
}
lextale <- rbindlist(data_list)

lextalescore <- lextale[-seq(1, nrow(lextale), 2),]
lextalescore <- separate(data = lextalescore, col = Number.of.correctly.identified.words..35, into = c("text", "score"), sep = "\\: ")
lextalescore[,2] <- as.numeric(unlist(lextalescore[,2]))

Pnames <- A
lextalescore$text <- Pnames

#### Interference Phase: Spanish adaptive learning and posttest ####
#### Adaptive learning task ####
# count number of exposure and learning success after the frist two rounds of this test
#B <- c(601:603, 608, 610:620, 622:631)# no adaptive picture naming fuile available for pp 604
data_list <- list()
for (i in 1:length(A)){
  pNumber = A[i]
  #setwd(file.path("//cnas.ru.nl/Wrkgrp/L2-Attrition-Mickan/RESULTS_EXP1/", pNumber))
  wd <- paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/BACK-UP/", pNumber, "/", pNumber, "_PicNaming", sep="")
  setwd(wd)
  infile1 <- paste(pNumber,"LearnPicNaming_A.txt",sep="_")
  
  currentFile <- as.data.frame(read.delim(infile1, stringsAsFactors=FALSE, sep = "\t", header = T, skipNul = TRUE))
  
  if (length(currentFile[currentFile$Error == 999,]$Error) > 0){
    currentFile[currentFile$Error == 999,]$Error<-1
  }
  
  data_list[[i]] <- currentFile
  
  print(A[i])
}
adap <- rbindlist(data_list)
blocks <- data.frame(tapply(adap$Block_nr, adap$Subject_nr,max))

adaptive <- adap[adap$Block_nr <3,] # keep only the second block
#adaptive <- adaptive[adaptive$Block_nr >1,]
successAdap <- 1-tapply(adaptive$Error, adaptive$Subject_nr,mean)

#### Exposure per item/pp ####
exposures<-data.frame(table(adap$Item, adap$Subject_nr))
exposures <- exposures[exposures$Freq != 0,]
exposures$Freq <- exposures$Freq + 8
expavg <- data.frame(tapply(exposures$Freq, exposures$Var2, mean))

######## Forgetting score ########
# difference between error rates in interference and no interfernce condition
forgetting <- data.frame(tapply(post$Error, list(post$Subject_nr, post$Condition), mean, na.rm = T))
forgetting$Difference <- forgetting$X2 - forgetting$X1
forgetting2 <- data.frame(tapply(post$RTdiff, list(post$Subject_nr, post$Condition), mean, na.rm = T))
forgetting2$Difference <- forgetting2$X1 - forgetting2$X2
forgetting$ForgettingRT <- forgetting2$Difference
forgetting$Interference_RT <- forgetting2$X1
forgetting$NoInterference_RT <- forgetting2$X2
colnames(forgetting) <- c("Interference_Error", "NoInterference_Error","Difference_Error", "Difference_RT", "Interference_RT","NoInterference_RT")

### Mean learning success at posttest in Spanish
data_list <- list()
for (i in 1:length(A)){
  pNumber = A[i]
  #setwd(file.path("//cnas.ru.nl/Wrkgrp/L2-Attrition-Mickan/RESULTS_EXP1/", pNumber))
  wd <- paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/BACK-UP/", pNumber, "/", pNumber, "_PosttestSpa", sep="")
  setwd(wd)
  infile1 <- paste(pNumber,"Posttest.txt",sep="_")
  
  currentFile <- as.data.frame(read.delim(infile1, stringsAsFactors=FALSE, sep = "\t", header = T, skipNul = TRUE))
  
  if (length(currentFile[currentFile$Error == 999,]$Error) > 0){
    currentFile[currentFile$Error == 999,]$Error<-1
  }
  
  data_list[[i]] <- currentFile
  
  print(A[i])
}
posttestSpanish <- rbindlist(data_list)
m1 <- 100-(tapply(posttestSpanish$Error, posttestSpanish$Subject_nr, mean)*100)
m2 <- tapply(posttestSpanish$VoiceOnset, posttestSpanish$Subject_nr, mean)

### Correlations #### 
correlations <- matrix(nrow = length(A), ncol = 7)
for (i in 1:length(A)) {
  pNumber = A[i]
  correlations[i,1] <- pNumber
  correlations[i,2] <- forgetting[i,3]                                    # Forgetting score error rate
  correlations[i,3] <- forgetting[i,4]                                    # Forgetting score RT
  correlations[i,4] <- m1[[i]]                                              # Learning success Spanish posttest
  #correlations[i,5] <- blocks[[i,1]]                                         # Number of blocks in adaptive learning
  #correlations[i,6] <- successAdap[[i]]                                    # Percent learned after second adaptive round
  #correlations[i,7] <- expavg[[i,1]]                                       # Average exposures to items (minimum 8)
  num <- which(tolower(as.character(rownames(pretest)))== pNumber)
  correlations[i,5] <- pretest[[num,1]]                                # Pretest percent known among first 101 words
  #correlations[i,15] <- LBQ$SRmean[i]                                      # Self-ratings average Spanish  
  #correlations[i,19] <- LBQ$SpanishExposureLengthMonth[i]                  # Spanish exposure in month
  #correlations[i,20] <- LBQ$FreqUseTotal[i]                                # Amount of time spent with spanish per week
  correlations[i,6] <- m2[[i]]                                     # Mean RT during Spanish Posttest
  num <- which(tolower(as.character(rownames(lextalescore)))== pNumber)
  correlations[i,7] <- lextalescore[[i,2]]                                   # Lextale score
}
#as.data.frame(correlations)->correlations
#colnames(correlations) <- c("Pnumber","Forgetting_Error","Forgetting_RT","MeanLearnSpa","PretestScore", "MeanRTSpanishPost","LexTale")

require(Hmisc)
corrtable <- round(rcorr(correlations),2)

cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

pvalues <- cor.test.p(corrtable)

for (i in 2:length(correlations)) {
  for (j in 2:length(correlations)) {
    if (pvalues[i,j]<0.05){
    } else {
      pvalues[i,j] = NA
      corrtable[i,j]=NA
    }
  }}

pvalues <- pvalues[-1,-1]
corrtable <- corrtable[-1,-1]
corrtable[corrtable==1] <- NA
pvalues[pvalues==0]<- NA

### Language background questionnaire ###
LBQ <- read.delim(here("LBQ_Honours.txt"))
# subset to only those participants that are included in the analyses
LBQfin <- LBQ[(LBQ$Subj %in% A),]

# age 
mean(LBQfin$Age)
sd(LBQfin$Age)

# AoA
mean(LBQfin$EnglishAoA)

# LoE
mean(LBQfin$EnglishLengthExposure)

# Frequency of use 
mean(LBQfin$EnglishFreqSpeakingMin)
mean(LBQfin$EnglishFreqListeningMin)
mean(LBQfin$EnglishFreqReadingMin)
mean(LBQfin$EnglishFreqWritingMin)

# Proficiency 
mean(LBQfin$EnglishProficiencySpeaking)
mean(LBQfin$EnglishProficiencyListening)
mean(LBQfin$EnglishProficiencyReading)
mean(LBQfin$EnglishProficiencyWriting)

#### Pretest English performance ####
A <- c(601:604, 608, 610:620, 622:631)
data_list <- list()
for (i in 1:length(A)){
  pNumber <- A[i]
  wd <- paste("//cnas.ru.nl/wrkgrp/STD-Honours-Mickan/BACK-UP/", pNumber,"/",pNumber,"_Pretest/", sep="")
  setwd(wd)
  infile1 <- paste(pNumber,"Pretest.txt",sep="_")
  pretest <- as.data.frame(read.delim(infile1, sep = "\t", header = T, skipNul = TRUE))
  
  data_list[[i]] <- pretest
  
  print(A[i])
}

pre <- rbindlist(data_list)

pre$Subject_nr <- as.factor(pre$Subject_nr)
pre$Condition <- as.factor(pre$Condition)
pre$Trial_nr <- as.factor(pre$Trial_nr)
pre$English_Label <- as.factor(pre$English_Label)

mean(tapply(pre$Unknown, pre$Subject_nr, mean))
pretest <- as.data.frame(tapply(pre$Unknown, pre$Subject_nr, mean))
cor.test(pretest$`tapply(pre$Unknown, pre$Subject_nr, mean)`, lextalescore$score)

for (i in 1:nrow(LBQ)){
  LBQ$MeanProfEng[i] <- sum(LBQ$EnglishProficiencyListening[i], LBQ$EnglishProficiencyReading[i], LBQ$EnglishProficiencySpeaking[i], LBQ$EnglishProficiencyWriting[i])/4}
LBQfin <- LBQ[(LBQ$Subj %in% A),]
cor.test(pretest$`tapply(pre$Unknown, pre$Subject_nr, mean)`, LBQfin$EnglishAoA)

