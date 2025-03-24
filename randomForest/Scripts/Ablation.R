library(randomForest)
library(missForest)
library(pROC)
library(ggplot2)
library(MASS)
library(caret)
library(data.table)

StatusFile = read.csv("Status_Mario_Brockman.csv", sep = ",", header = T)

health_outcome = StatusFile$HIVStatus[1:91]
health_outcome[which(health_outcome == 2)] = 1
health_outcome[which(health_outcome == 3)] = 1
health_outcome[which(health_outcome == 4)] = 1

features_all.imp = read.csv("features_withMissForestImputation_IR_INR_together_jul25_63Features.csv")
length(features_all.imp[1,])


'%!in%' <- function(x,y)!('%in%'(x,y))
x = features_all.imp[1:91,2:64]
y = health_outcome



n_individuals = length(y)
k = 5
Limit = 1000
auc=vector()
varIMPORT.cumAVG.i = vector()
varIMPORT.cumAVG.total = vector()
testSetCount = numeric(n_individuals)
CorrectClassification_counter = numeric(n_individuals)
INCorrectClassification_counter = numeric(n_individuals)
accuracy_vector_generate = vector()
sensitivity = vector()
specificity = vector()
precision = vector()
f1_score = vector()
n_individuals.training = 40
set.seed(42)
for(i in seq(Limit))
{
  #20:20 randomly sampled HIV-:HIV+
  KeepHIVIndices = sample(24:90, size = 20)
  KeepHIVNegativeIndices = sample(1:23, size = 20)
  indexes = c(KeepHIVNegativeIndices, KeepHIVIndices)
  y.downSamp.TRAIN = health_outcome[indexes]
  x.downSamp.TRAIN = x[indexes,]
  
  #For test: Keep all 3 remaining HIV- indiciduals. Sample 12 remaining HIV+ individuals
  seqHIV.Neg.Test = which(seq(1,23,1) %!in%  KeepHIVNegativeIndices)
  seqHIV.Pos.Test = which(seq(24,90,1) %!in%  KeepHIVIndices)
  HIV.Neg.Test.finalSeq = sample(seq(24,90,1)[seqHIV.Pos.Test], size = 12)
  indexes.TEST = c(seqHIV.Neg.Test, HIV.Neg.Test.finalSeq)
  #Check to make sure indices in TEST are not in TRAIN
  #which(indexes %in% indexes.TEST)
  #which(HIV.Neg.Test.finalSeq %in% KeepHIVIndices)
  y.TEST = health_outcome[indexes.TEST]
  x.TEST = x[indexes.TEST,]
  
  # Train the Random Forest model
  rf_model <- randomForest(x = x.downSamp.TRAIN, y = as.factor(y.downSamp.TRAIN), ntree = 64)
  
  #Keep track of a counter of how many times each ID is put into the test group
  testSetCount[indexes.TEST] = testSetCount[indexes.TEST] + 1
  #Predict for Accuracy
  predictions <- predict(rf_model, newdata = x.TEST)
  accuracy <- sum(predictions == y.TEST) / length(y.TEST)
  accuracy_vector_generate[i] = accuracy
  #Predict for probabilities
  predictions.prob = predict(rf_model,  newdata = x.TEST, type='prob')
  auc[i]=roc(as.factor(y.TEST), predictions.prob[,1], quiet = TRUE)$auc
  
  # Create a confusion matrix
  conf_matrix <- confusionMatrix(predictions, factor(y.TEST))
  
  # Calculate and store sensitivity and specificity and precision
  sensitivity[i] <- conf_matrix$byClass['Sensitivity']
  specificity[i] <- conf_matrix$byClass['Specificity']
  precision[i] <- conf_matrix$byClass['Pos Pred Value']
  f1_score[i] <- 2 * ((precision[i] * sensitivity[i]) / (precision[i] + sensitivity[i]))
  
  SeqCorrect = which(predictions == y.TEST)
  CorrectClassification_IDS = as.numeric(labels(predictions[SeqCorrect]))
  INCorrectClassification_IDS = as.numeric(labels(predictions[-SeqCorrect]))
  CorrectClassification_counter[CorrectClassification_IDS] = CorrectClassification_counter[CorrectClassification_IDS] + 1
  INCorrectClassification_counter[INCorrectClassification_IDS] = INCorrectClassification_counter[INCorrectClassification_IDS] + 1
  
  if(i == 1){
    varIMPORT.cumAVG = varImp(rf_model)
  }else{
    varIMPORT.cumAVG = (varImp(rf_model) + (i-1)*varIMPORT.cumAVG)/(i)
  }
  
  if(i %% 100 == 0){
    print(i)
  }    
}


########################################################
#####   Forward Ablation
########################################################

FeatureImportance = (as.data.frame(varIMPORT.cumAVG$Overall))


n_individuals = length(y)
#k = 5
Limit = 1000
######Reset these after jth loop ##########
auc=vector()
varIMPORT.cumAVG.i = vector()
varIMPORT.cumAVG.total = vector()
testSetCount = numeric(n_individuals)
CorrectClassification_counter = numeric(n_individuals)
INCorrectClassification_counter = numeric(n_individuals)
accuracy_vector_generate = vector()
sensitivity = vector()
specificity = vector()
precision = vector()
f1_score = vector()
#############################################

###### Data storage after each jth loop (keep all these) ##########
nFeatures = length(x[1,])
n.columns = nFeatures - 1

#varIMPORT.cumAVG.i.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = length(features_all.imp[1,])))
#varIMPORT.cumAVG.total.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = length(features_all.imp[1,])))
TREES = 200
#testSetCount.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = n_individuals))
CorrectClassification_counter.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = n_individuals))
INCorrectClassification_counter.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = n_individuals))

ProbHIVpos.MEAN = data.frame(matrix(ncol = n.columns, nrow = n_individuals))
ProbHIVpos.SD  = data.frame(matrix(ncol = n.columns, nrow = n_individuals))


testSetCount.forwardAblation = data.frame(matrix(ncol = n.columns, nrow = n_individuals))

#Columns will all be the jth iteration while rows the inner run loop summary evaluation
#Save to column j-1 in the for loop. 
auc.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
accuracy_vector_generate.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
sensitivity.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
specificity.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
precision.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
f1_score.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = Limit))

set.seed(42)
nFeatures = length(x[1,])
for(j in 2:nFeatures){#Outer loop for feature manipulation

  if(j == 2){
    keep = SeqOrder[c(1,2)]
    x.forwardAb = x[,keep]
  }else{
    keep = SeqOrder[1:j]
    x.forwardAb = x[,keep]
  }
  #For each particular j need to store all the model summary measures in a data.frame
  for(i in seq(Limit)){#Inner loop for running through the RF algorithm 
    
    #20:20 randomly sampled HIV-:HIV+
    KeepHIVIndices = sample(24:91, size = 20)
    KeepHIVNegativeIndices = sample(1:23, size = 20)
    indexes = c(KeepHIVNegativeIndices, KeepHIVIndices)
    y.downSamp.TRAIN = health_outcome[indexes]
    x.downSamp.TRAIN = x.forwardAb[indexes,]
    
    #For test: Keep all 3 remaining HIV- individuals. Sample 12 remaining HIV+ individuals
    seqHIV.Neg.Test = which(seq(1,23,1) %!in%  KeepHIVNegativeIndices)
    seqHIV.Pos.Test = which(seq(24,91,1) %!in%  KeepHIVIndices)
    HIV.Neg.Test.finalSeq = sample(seq(24,91,1)[seqHIV.Pos.Test], size = 12)
    
    indexes.TEST = c(seqHIV.Neg.Test, HIV.Neg.Test.finalSeq)
    y.TEST = health_outcome[indexes.TEST]
    x.TEST = x.forwardAb[indexes.TEST,]
    
    #Keep track of a counter of how many times each ID is put into the test group
    testSetCount[indexes.TEST] = testSetCount[indexes.TEST] + 1
    
    # Train the Random Forest model
    rf_model <- randomForest(x = x.downSamp.TRAIN, y = as.factor(y.downSamp.TRAIN), ntree = TREES)
    
    #Keep track of a counter of how many times each ID is put into the test group
    #testSetCount[indexes.TEST] = testSetCount[indexes.TEST] + 1
    #Predict for Accuracy
    predictions <- predict(rf_model, newdata = x.TEST)
    accuracy <- sum(predictions == y.TEST) / length(y.TEST)
    accuracy_vector_generate[i] = accuracy
    #Predict for probabilities
    predictions.prob = predict(rf_model,  newdata = x.TEST, type='prob')
    auc[i]=roc(as.factor(y.TEST), predictions.prob[,1], quiet = TRUE)$auc
    
    # Create a confusion matrix
    conf_matrix <- confusionMatrix(predictions, factor(y.TEST))
    
    # Calculate and store sensitivity and specificity and precision
    sensitivity[i] <- conf_matrix$byClass['Sensitivity']
    specificity[i] <- conf_matrix$byClass['Specificity']
    precision[i] <- conf_matrix$byClass['Pos Pred Value']
    f1_score[i] <- 2 * ((precision[i] * sensitivity[i]) / (precision[i] + sensitivity[i]))
    
    SeqCorrect = which(predictions == y.TEST)
    CorrectClassification_IDS = as.numeric(labels(predictions[SeqCorrect]))
    INCorrectClassification_IDS = as.numeric(labels(predictions[-SeqCorrect]))
    CorrectClassification_counter[CorrectClassification_IDS] = CorrectClassification_counter[CorrectClassification_IDS] + 1
    INCorrectClassification_counter[INCorrectClassification_IDS] = INCorrectClassification_counter[INCorrectClassification_IDS] + 1
    
    Individual.Predictions[i,indexes.TEST] = as.numeric(predictions.prob[,2])
    
    
    
    if(i %% 250 == 0){
      print(c(j,i))
    }    
    
  }#End i
  #Individual.Predictions[,90]
  column_means <- apply(Individual.Predictions, 2, function(x) mean(x, na.rm = TRUE))
  column_sd <- apply( Individual.Predictions, 2, function(x) sd(x, na.rm = TRUE))

  ProbHIVpos.MEAN[,(j-1)] = as.numeric(column_means)
  ProbHIVpos.SD[,(j-1)] = as.numeric(column_sd)
  #ProbHIVpos.MEAN[39,1:60]
  #Store jth data
  #testSetCount.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = n_individuals))
  CorrectClassification_counter.forwardAbilation[,(j-1)] = CorrectClassification_counter
  INCorrectClassification_counter.forwardAbilation[,(j-1)] = INCorrectClassification_counter
  
  testSetCount.forwardAblation[,(j-1)] = testSetCount
  auc.forwardAbilation[,(j-1)] = auc
  accuracy_vector_generate.forwardAbilation[,(j-1)] = accuracy_vector_generate
  sensitivity.forwardAbilation[,(j-1)] = sensitivity
  specificity.forwardAbilation[,(j-1)] = specificity
  precision.forwardAbilation[,(j-1)] = precision
  f1_score.forwardAbilation[,(j-1)] = f1_score
  
  #Wipe temporary objects for next loop
  auc=vector()
  varIMPORT.cumAVG.i = vector()
  varIMPORT.cumAVG.total = vector()
  testSetCount = numeric(n_individuals)
  CorrectClassification_counter = numeric(n_individuals)
  INCorrectClassification_counter = numeric(n_individuals)
  accuracy_vector_generate = vector()
  sensitivity = vector()
  specificity = vector()
  precision = vector()
  f1_score = vector()
  Individual.Predictions = data.frame(matrix(ncol = length(y), nrow = Limit))
}#End j


save(CorrectClassification_counter.forwardAbilation, INCorrectClassification_counter.forwardAbilation, auc.forwardAbilation, 
     accuracy_vector_generate.forwardAbilation,sensitivity.forwardAbilation,  specificity.forwardAbilation, 
     precision.forwardAbilation, f1_score.forwardAbilation, ProbHIVpos.MEAN,ProbHIVpos.SD, file = "forwardAblation.RData")

########################################################
#####   Reverse Ablation
########################################################


#####
n_individuals = length(y)
#k = 5
Limit = 1000
######Reset these after jth loop ##########
auc=vector()
varIMPORT.cumAVG.i = vector()
varIMPORT.cumAVG.total = vector()
testSetCount = numeric(n_individuals)
CorrectClassification_counter = numeric(n_individuals)
INCorrectClassification_counter = numeric(n_individuals)
accuracy_vector_generate = vector()
sensitivity = vector()
specificity = vector()
precision = vector()
f1_score = vector()
#############################################

###### Data storage after each jth loop (keep all these) ##########
n.columns = nFeatures - 1

#varIMPORT.cumAVG.i.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = length(features_all.imp[1,])))
#varIMPORT.cumAVG.total.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = length(features_all.imp[1,])))

#testSetCount.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = n_individuals))
CorrectClassification_counter.reverseAbilation = data.frame(matrix(ncol = n.columns, nrow = n_individuals))
INCorrectClassification_counter.reverseAbilation = data.frame(matrix(ncol = n.columns, nrow = n_individuals))

#Columns will all be the jth iteration while rows the inner run loop summary evaluation
#Save to column j-1 in the for loop. 
auc.reverseAbilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
accuracy_vector_generate.reverseAbilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
sensitivity.reverseAbilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
specificity.reverseAbilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
precision.reverseAbilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
f1_score.reverseAbilation = data.frame(matrix(ncol = n.columns, nrow = Limit))

set.seed(42)

nFeatures = length(x[1,])

for(j in 1:(nFeatures-2)){#Outer loop for feature manipulation
  
  if(j == 1){
    x.reverseAb = x[,-SeqOrder[1]]
  }else{
    x.reverseAb = x[,-SeqOrder[1:j]]
  }
  #For each particular j need to store all the model summary measures in a data.frame
  for(i in seq(Limit)){#Inner loop for running through the RF algorithm 
    
    #20:20 randomly sampled HIV-:HIV+
    KeepHIVIndices = sample(24:91, size = 20)
    KeepHIVNegativeIndices = sample(1:23, size = 20)
    indexes = c(KeepHIVNegativeIndices, KeepHIVIndices)
    y.downSamp.TRAIN = health_outcome[indexes]
    x.downSamp.TRAIN = x.reverseAb[indexes,]
    
    #For test: Keep all 3 remaining HIV- individuals. Sample 12 remaining HIV+ individuals
    seqHIV.Neg.Test = which(seq(1,23,1) %!in%  KeepHIVNegativeIndices)
    seqHIV.Pos.Test = which(seq(24,91,1) %!in%  KeepHIVIndices)
    HIV.Neg.Test.finalSeq = sample(seq(24,91,1)[seqHIV.Pos.Test], size = 12)
    
    indexes.TEST = c(seqHIV.Neg.Test, HIV.Neg.Test.finalSeq)
    y.TEST = health_outcome[indexes.TEST]
    x.TEST = x.reverseAb[indexes.TEST,]
    
    #Keep track of a counter of how many times each ID is put into the test group
    #testSetCount[indexes.TEST] = testSetCount[indexes.TEST] + 1
    
    # Train the Random Forest model
    rf_model <- randomForest(x = x.downSamp.TRAIN, y = as.factor(y.downSamp.TRAIN), ntree = TREES)
    
    #Keep track of a counter of how many times each ID is put into the test group
    testSetCount[indexes.TEST] = testSetCount[indexes.TEST] + 1
    #Predict for Accuracy
    predictions <- predict(rf_model, newdata = x.TEST)
    accuracy <- sum(predictions == y.TEST) / length(y.TEST)
    accuracy_vector_generate[i] = accuracy
    #Predict for probabilities
    predictions.prob = predict(rf_model,  newdata = x.TEST, type='prob')
    auc[i]=roc(as.factor(y.TEST), predictions.prob[,1], quiet = TRUE)$auc
    
    # Create a confusion matrix
    conf_matrix <- confusionMatrix(predictions, factor(y.TEST))
    
    # Calculate and store sensitivity and specificity and precision
    sensitivity[i] <- conf_matrix$byClass['Sensitivity']
    specificity[i] <- conf_matrix$byClass['Specificity']
    precision[i] <- conf_matrix$byClass['Pos Pred Value']
    f1_score[i] <- 2 * ((precision[i] * sensitivity[i]) / (precision[i] + sensitivity[i]))
    
    SeqCorrect = which(predictions == y.TEST)
    CorrectClassification_IDS = as.numeric(labels(predictions[SeqCorrect]))
    INCorrectClassification_IDS = as.numeric(labels(predictions[-SeqCorrect]))
    CorrectClassification_counter[CorrectClassification_IDS] = CorrectClassification_counter[CorrectClassification_IDS] + 1
    INCorrectClassification_counter[INCorrectClassification_IDS] = INCorrectClassification_counter[INCorrectClassification_IDS] + 1
    
    if(i %% 250 == 0){
      print(c(j,i))
    }    
    
  }#End i
  
  #Store jth data
  #testSetCount.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = n_individuals))
  CorrectClassification_counter.reverseAbilation[,(j)] = CorrectClassification_counter
  INCorrectClassification_counter.reverseAbilation[,(j)] = INCorrectClassification_counter
  auc.reverseAbilation[,(j)] = auc
  accuracy_vector_generate.reverseAbilation[,(j)] = accuracy_vector_generate
  sensitivity.reverseAbilation[,(j)] = sensitivity
  specificity.reverseAbilation[,(j)] = specificity
  precision.reverseAbilation[,(j)] = precision
  f1_score.reverseAbilation[,(j)] = f1_score
  
  #Wipe temporary objects for next loop
  auc=vector()
  varIMPORT.cumAVG.i = vector()
  varIMPORT.cumAVG.total = vector()
  testSetCount = numeric(n_individuals)
  CorrectClassification_counter = numeric(n_individuals)
  INCorrectClassification_counter = numeric(n_individuals)
  accuracy_vector_generate = vector()
  sensitivity = vector()
  specificity = vector()
  precision = vector()
  f1_score = vector()
}#End j


save(CorrectClassification_counter.reverseAbilation, INCorrectClassification_counter.reverseAbilation, auc.reverseAbilation,
     accuracy_vector_generate.reverseAbilation, sensitivity.reverseAbilation, specificity.reverseAbilation, 
     precision.reverseAbilation, f1_score.reverseAbilation, file = "ReverseAblation.RData")

###############################################################
####### Summarize ablation analysis and make plots #######
###############################################################

#Get all the medians for a plot
reverse.median_f1 = numeric(61)
reverse.median_auc = numeric(61)

for(j in 1:61){

  reverse.median_f1[j] = mean(na.omit(f1_score.reverseAbilation[,j]))
  reverse.median_auc[j] = mean(na.omit(auc.reverseAbilation[,j]))

}


#Get all the medians for a plot
forward.median_f1 = numeric(62)
forward.median_auc = numeric(62)

for(j in 1:62){
  
  forward.median_f1[j] = mean(na.omit(f1_score.forwardAbilation[,j]))
  foreward.median_auc[j] = mean(na.omit(auc.forwardAbilation[,j]))
  
}

forward_index = seq(1, 62, 1)
reverse_index = seq(1, 61, 1)

quartz()
plot(forward_index, forward.median_f1, pch = 15, col = "#44C0DE", ylim = c(0, 1), ylab = "mean proportion", yaxt = "n", xaxt = "n", xlab = "Ablation index", cex = 1.7)
points(forward_index, foreward.median_auc, pch = 0, col = "#44C0DE", cex = 1.7)
points(reverse_index, reverse.median_f1, pch = 17, col = "#F6B541", cex = 1.7)
points(reverse_index, reverse.median_auc, pch = 2, col = "#F6B541", cex = 1.7)

axis(2,las= 1, tck=0.02,mgp=c(3, .5, 0), cex.axis = 1.5)
axis(1,las= 1, tck=0.02,mgp=c(3, .5, 0), cex.axis = 1.5)

####### Plot trajectories of incorrect classification as a function of forward ablation procedure


INCorrectClassification_counter.forwardAbilation[,1]/testSetCount.forwardAblation[,1]

quartz()
ID = seq(1,91,1)
par(mar=c(11,4,4,4))
barplot(INCorrectClassification_counter.forwardAbilation[,60]/testSetCount.forwardAblation[,60],names.arg = ID, col="#69b3a2", las=2, ylim = c(0, 1.1), ylab = "Probability of incorrect classification")
box(bty = "o")
#So I want to plot the trajectories as a function of index for each person. 
index = seq(1,62,1) #Always xAxis
#loop over all the 91 individuals? 

probTrajectories.forwardAblation = data.frame(matrix(ncol = 91, nrow = 61))


quartz()
#Set up the plots to add points to it
plot(as.numeric(INCorrectClassification_counter.forwardAbilation[39,1:62]/testSetCount.forwardAblation[39,1:62]), ylim = c(0,1), xlim = c(1,62))
for(j in 1:62){#loop over index
  for(i in 1:91){
    probTrajectories.forwardAblation[j,i] = INCorrectClassification_counter.forwardAbilation[i,j]/testSetCount.forwardAblation[i,j]
    
  }
}
probTrajectories.forwardAblation[index,39]
col = "#4573AD"
quartz()
plot(index, probTrajectories.forwardAblation[index,39], ylim = c(0,1), xlim = c(1,62), ylab = "", xlab = "", col = "#B43B2B", pch = 20, xaxt = "n", yaxt = "n")
axis(2,las= 1, tck=0.02,mgp=c(3, .5, 0), cex.axis = 1.5)
axis(1,las= 1, tck=0.02,mgp=c(3, .5, 0), cex.axis = 1.5)
lines(index, probTrajectories.forwardAblation[index,39], col = rgb(180/255, 59/255, 43/255, 0.5))
#abline(h = 0.5, lty = 2)
for(i in 1:91){
  if(i <= 23){
    points(index, probTrajectories.forwardAblation[index,i], ylim = c(0,1), xlim = c(1,62), ylab = "", xlab = "", col = "#4573AD", pch = 15)
    lines(index,  probTrajectories.forwardAblation[index,i], col = rgb(69/255, 115/255, 173/255, 0.5))
  }
  if(i > 23){
    points(index, probTrajectories.forwardAblation[index,i], ylim = c(0,1), xlim = c(1,62), ylab = "", xlab = "", col = "#B43B2B", pch = 20)
    lines(index,  probTrajectories.forwardAblation[index,i], col = rgb(180/255, 59/255, 43/255, 0.5))
  }
}

probTrajectories.forwardAblation[1,1:91]/probTrajectories.forwardAblation[60,1:91]
lines(index,  probTrajectories.forwardAblation[index,6], col = rgb(69/255, 115/255, 173/255, 0.2), lwd = 15) #ID 6
lines(index,  probTrajectories.forwardAblation[index,12], col = rgb(69/255, 115/255, 173/255, 0.2), lwd = 15) #ID 12

lines(index,  probTrajectories.forwardAblation[index,39], col = rgb(180/255, 59/255, 43/255, 0.2), lwd = 15) #ID 39
lines(index,  probTrajectories.forwardAblation[index,64], col = rgb(180/255, 59/255, 43/255, 0.2), lwd = 15) #ID 64
lines(index,  probTrajectories.forwardAblation[index,76], col = rgb(180/255, 59/255, 43/255, 0.2), lwd = 15) #ID 76



ProbHIVpos.MEAN[,(j-1)] = column_means
ProbHIVpos.SD[,(j-1)] = column_sd
ProbHIVpos.SD[index, 39]
quartz()
plot(index, ProbHIVpos.MEAN[39,index], ylim = c(0,1), xlim = c(1,62), ylab = "", xlab = "", col = "white", pch = 20, xaxt = "n", yaxt = "n")
axis(2,las= 1, tck=0.02,mgp=c(3, .7, 0), cex.axis = 2.2)
axis(1,las= 1, tck=0.02,mgp=c(3, .7, 0), cex.axis = 2.2)
lines(index, ProbHIVpos.MEAN[39, index], col = "white")
#abline(h = 0.5, lty = 2)
for(i in 1:91){
  if(i <= 23){
    points(index, ProbHIVpos.MEAN[i,index], ylim = c(0,1), xlim = c(1,62), ylab = "", xlab = "", col = rgb(69/255, 115/255, 173/255, 0.2), pch = 15)
    lines(index,  ProbHIVpos.MEAN[i, index], col = rgb(69/255, 115/255, 173/255, 0.2))
  }
  if(i > 23){
    points(index, ProbHIVpos.MEAN[i, index], ylim = c(0,1), xlim = c(1,62), ylab = "", xlab = "", col = rgb(180/255, 59/255, 43/255, 0.2), pch = 20)
    lines(index,  ProbHIVpos.MEAN[i, index], col = rgb(180/255, 59/255, 43/255, 0.2))
  }
}
abline(h = 0.5, lty = 2, lwd= 3.0)

lines(index,  ProbHIVpos.MEAN[6, index], col = rgb(69/255, 115/255, 173/255, 0.8), lwd = 15) #ID 6
#lines(index,  ProbHIVpos.MEAN[12, index], col = rgb(69/255, 115/255, 173/255, 0.8), lwd = 15) #ID 12

lines(index,  ProbHIVpos.MEAN[39, index], col = rgb(180/255, 59/255, 43/255, 0.8), lwd = 15) #ID 39
lines(index,  ProbHIVpos.MEAN[72, index], col = rgb(180/255, 59/255, 43/255, 0.8), lwd = 15) #ID 72

lines(index,  ProbHIVpos.MEAN[64, index], col = rgb(180/255, 59/255, 43/255, 0.8), lwd = 15) #ID 64
lines(index,  ProbHIVpos.MEAN[76, index], col = rgb(180/255, 59/255, 43/255, 0.8), lwd = 15) #ID 76
lines(index,  ProbHIVpos.MEAN[72, index], col = rgb(180/255, 59/255, 43/255, 0.8), lwd = 15) #ID 72
lines(index,  ProbHIVpos.MEAN[69, index], col = rgb(180/255, 59/255, 43/255, 0.8), lwd = 15) #ID 69
lines(index,  ProbHIVpos.MEAN[88, index], col = rgb(180/255, 59/255, 43/255, 0.8), lwd = 15) #ID 88
lines(index,  ProbHIVpos.MEAN[87, index], col = rgb(180/255, 59/255, 43/255, 0.8), lwd = 15) #ID 87

abline(v = 9, lwd = 2.0)



####################################################################
################# Jessica's Suggestion - use the colour code to compute the AUCs of the respective data type 
####################################################################
### If we just used cytokines, versus blood, saliva data how meaningful are the classifications? 

#IgG serum features
#IgG Saliva features
#IgA Saliva features
#Cytokine features
#Neutralization/ACE2 displacement

Limit = 1000
######Reset these after jth loop ##########
TREES = 200
Ntrain = 17
Ntest_HIVpos = 18 
PID = sample(seq(2*Ntrain)) 
n_individuals = length(y)
auc=vector()
varIMPORT.cumAVG.i = vector()
varIMPORT.cumAVG.total = vector()
testSetCount = numeric(n_individuals)
CorrectClassification_counter = numeric(n_individuals)
INCorrectClassification_counter = numeric(n_individuals)
accuracy_vector_generate = vector()
sensitivity = vector()
specificity = vector()
precision = vector()
f1_score = vector()
#############################################

seq.IgG_serum = seq(1, 28, 1)
seq.IgG_saliva = seq(29,38, 1)
seq.IgA_saliva = seq(39, 48, 1)
seq.Cytokine = seq(49, 57, 1)
seq.neut_ACE2 = seq(58, 63, 1)

x = features_all.imp[1:91,2:64]
y = health_outcome


###### Data storage after each jth loop (keep all these) ##########
n.columns = 6

#varIMPORT.cumAVG.i.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = length(features_all.imp[1,])))
#varIMPORT.cumAVG.total.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = length(features_all.imp[1,])))

#testSetCount.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = n_individuals))

CorrectClassification_counter.Abilation = data.frame(matrix(ncol = n.columns, nrow = n_individuals))
INCorrectClassification_counter.Abilation = data.frame(matrix(ncol = n.columns, nrow = n_individuals))

#Columns will all be the jth iteration while rows the inner run loop summary evaluation
#Save to column j-1 in the for loop. 
auc.Abilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
accuracy_vector_generate.Abilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
sensitivity.Abilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
specificity.Abilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
precision.Abilation = data.frame(matrix(ncol = n.columns, nrow = Limit))
f1_score.Abilation = data.frame(matrix(ncol = n.columns, nrow = Limit))

# ############################### k-fold cross-validation function
kfoldcrossvalidation <- function(k, features, health_outcome, pid, validation) {
  ev.pred = matrix(NA, nrow = k, ncol = nrow(validation))
  w = matrix(NA, nrow = k, ncol = ncol(features))
  pred = rep(NA, length(health_outcome))
  folds <- cut(pid, breaks = k, labels = FALSE) 
  
  for (i in 1:k) {
    test_indices <- which(folds == i, arr.ind = TRUE)
    test_data <- features[test_indices, ]
    train_data <- features[-test_indices, ]
    train_labels <- health_outcome[-test_indices]
    rf_model <- randomForest(x = train_data, y = as.factor(train_labels))
    w[i, ] = importance(rf_model)
    pred[test_indices] = predict(rf_model, newdata = test_data, type = 'prob')[, 2]
    ev.pred[i, ] = predict(rf_model, newdata = validation, type = 'prob')[, 2]
  }
  
  return(list(pred = pred, weights = w, ev.pred = ev.pred))
}

set.seed(42)

for(j in 1:6){#Outer loop for feature manipulation
  
  if(j == 1){# all features
    x = features_all.imp[1:91,2:64]
    y = health_outcome
  }
  if(j == 2){# IgG serum
    x = features_all.imp[1:91,2:64]
    y = health_outcome
    
    x = x[, seq.IgG_serum]
  }
  if(j == 3){# IgG saliva
    x = features_all.imp[1:91,2:64]
    y = health_outcome
    
    x = x[, seq.IgG_saliva]
    
  }
  if(j == 4){# IgA saliva
    x = features_all.imp[1:91,2:64]
    y = health_outcome
    
    x = x[, seq.IgA_saliva]
    
  }
  if(j == 5){# Cytokines
    x = features_all.imp[1:91,2:64]
    y = health_outcome
    
    x = x[, seq.Cytokine]
    
  }
  if(j == 6){# Neut/ACE2
    x = features_all.imp[1:91,2:64]
    y = health_outcome
    
    x = x[, seq.neut_ACE2]
    
  }
  

  #For each particular j need to store all the model summary measures in a data.frame
  for(i in seq(Limit)){#Inner loop for running through the RF algorithm 
    #20:20 randomly sampled HIV-:HIV+
    KeepHIVIndices = sample(24:91, size = Ntrain)
    KeepHIVNegativeIndices = sample(1:23, size = Ntrain)
    indexes = c(KeepHIVNegativeIndices, KeepHIVIndices)
    y.downSamp.TRAIN = health_outcome[indexes]
    x.downSamp.TRAIN = x[indexes,]
    
    #For test: Keep all 3 remaining HIV- individuals. Sample 12 remaining HIV+ individuals
    seqHIV.Neg.Test = which(seq(1,23,1) %!in%  KeepHIVNegativeIndices)
    seqHIV.Pos.Test = which(seq(24,91,1) %!in%  KeepHIVIndices)
    HIV.Neg.Test.finalSeq = sample(seq(24,91,1)[seqHIV.Pos.Test], size = 12)
    
    indexes.TEST = c(seqHIV.Neg.Test, HIV.Neg.Test.finalSeq)
    y.TEST = health_outcome[indexes.TEST]
    x.TEST = x[indexes.TEST,]
    
    #Keep track of a counter of how many times each ID is put into the test group
    #testSetCount[indexes.TEST] = testSetCount[indexes.TEST] + 1
    
    # Train the Random Forest model
    rf_model <- randomForest(x = x.downSamp.TRAIN, y = as.factor(y.downSamp.TRAIN), ntree = TREES)
    
    #Keep track of a counter of how many times each ID is put into the test group
    testSetCount[indexes.TEST] = testSetCount[indexes.TEST] + 1
    #Predict for Accuracy
    predictions <- predict(rf_model, newdata = x.TEST)
    accuracy <- sum(predictions == y.TEST) / length(y.TEST)
    accuracy_vector_generate[i] = accuracy
    #Predict for probabilities
    predictions.prob = predict(rf_model,  newdata = x.TEST, type='prob')
    #auc[i]=roc(as.factor(y.TEST), predictions.prob[,1], quiet = TRUE)$auc
    
    # Create a confusion matrix
    conf_matrix <- confusionMatrix(predictions, factor(y.TEST))
    
    # Calculate and store sensitivity and specificity and precision
    sensitivity[i] <- conf_matrix$byClass['Sensitivity']
    specificity[i] <- conf_matrix$byClass['Specificity']
    precision[i] <- conf_matrix$byClass['Pos Pred Value']
    f1_score[i] <- 2 * ((precision[i] * sensitivity[i]) / (precision[i] + sensitivity[i]))
    
    SeqCorrect = which(predictions == y.TEST)
    CorrectClassification_IDS = as.numeric(labels(predictions[SeqCorrect]))
    INCorrectClassification_IDS = as.numeric(labels(predictions[-SeqCorrect]))
    CorrectClassification_counter[CorrectClassification_IDS] = CorrectClassification_counter[CorrectClassification_IDS] + 1
    INCorrectClassification_counter[INCorrectClassification_IDS] = INCorrectClassification_counter[INCorrectClassification_IDS] + 1
    
    
    out = kfoldcrossvalidation(5, x.downSamp.TRAIN , y.downSamp.TRAIN, PID, x.TEST) #Actual Data
    
    auc[i]=suppressMessages(auc(roc(y.TEST, colMeans(out$ev.pred)), quiet = TRUE))
   # auc[i]=auc(roc(y.downSamp.TRAIN, out$pred))
    
    if(i %% 250 == 0){
      print(c(j,i))
    }    
    
  }#End i
  
  #Store jth data
  #testSetCount.forwardAbilation = data.frame(matrix(ncol = n.columns, nrow = n_individuals))
  CorrectClassification_counter.Abilation[,(j)] = CorrectClassification_counter
  INCorrectClassification_counter.Abilation[,(j)] = INCorrectClassification_counter
  auc.Abilation[,(j)] = auc
  accuracy_vector_generate.Abilation[,(j)] = accuracy_vector_generate
  sensitivity.Abilation[,(j)] = sensitivity
  specificity.Abilation[,(j)] = specificity
  precision.Abilation[,(j)] = precision
  f1_score.Abilation[,(j)] = f1_score
  
  #Wipe temporary objects for next loop
  auc=vector()
  varIMPORT.cumAVG.i = vector()
  varIMPORT.cumAVG.total = vector()
  testSetCount = numeric(n_individuals)
  CorrectClassification_counter = numeric(n_individuals)
  INCorrectClassification_counter = numeric(n_individuals)
  accuracy_vector_generate = vector()
  sensitivity = vector()
  specificity = vector()
  precision = vector()
  f1_score = vector()
}#End j


save(CorrectClassification_counter.Abilation, INCorrectClassification_counter.Abilation, auc.Abilation,
     accuracy_vector_generate.Abilation, sensitivity.Abilation, specificity.Abilation, 
     precision.Abilation, f1_score.Abilation, file = "Ablation_byDataType_withKfold.RData")

sd(na.omit(auc.Abilation[,2])) #IgG serum
na.omit(auc.Abilation[,3]) # IgG saliva
na.omit(auc.Abilation[,4]) # IgA saliva
na.omit(auc.Abilation[,5]) # Cytokines
sd(na.omit(auc.Abilation[,6])) #ACE2/Neutralization

seq.IgG_serum = seq(1, 28, 1)
seq.IgG_saliva = seq(29,38, 1)
seq.IgA_saliva = seq(39, 48, 1)
seq.Cytokine = seq(49, 57, 1)
seq.neut_ACE2 = seq(58, 63, 1)

length(seq.neut_ACE2)

barplot.colours = c(rep("#DC3A7A", 8), rep("#9E3934", 10),rep("#ED220D", 10), rep("skyblue", 10),rep("#3675B8", 10), rep("#F8B959", 9), rep("#94C652", 6))

mean(na.omit(auc.Abilation[,5]))


library(vioplot)
quartz()
vioplot(na.omit(auc.Abilation[,1]), #all features
        na.omit(auc.Abilation[,2]), #IgG serum
        na.omit(auc.Abilation[,3]), # IgG saliva
        na.omit(auc.Abilation[,4]), # IgA saliva
        na.omit(auc.Abilation[,5]), # Cytokines
        na.omit(auc.Abilation[,6]), #ACE2/Neutralization
        main = "", 
        at = seq(1,6,1), 
        names = c("All features \n (n = 63)", "IgG Serum \n (n = 28)", "IgG Saliva \n (n = 10)", "IgA Saliva \n (n = 10)", "Cytokines  \n (n = 9)", "ACE2/Neut. \n (n = 9)"), 
        col = c("#757476", "#9E3934", "skyblue", "#3675B8", "#F8B959", "#94C652"),
        border = c("#141414", "#D46D61", "#7E9EB2", "#70A0F2", "#AE8E59", "#D7FB92"),
        horizontal = FALSE, 
        notch = FALSE, 
        yaxt = "n",
        boxlwd = 1,
        lwd = 1, 
        las = 2, 
        cex.axis = 1.0,
        ylab = "AUC",
        xlab = "",
        ylim = c(0,1)
)
axis(2,las= 1, tck=0.02,mgp=c(3, .5, 0), cex.axis = 1.9)


quartz()
vioplot(na.omit(f1_score.Abilation[,1]), #all features
        na.omit(f1_score.Abilation[,2]), #IgG serum
        na.omit(f1_score.Abilation[,3]), # IgG saliva
        na.omit(f1_score.Abilation[,4]), # IgA saliva
        na.omit(f1_score.Abilation[,5]), # Cytokines
        na.omit(f1_score.Abilation[,6]), #ACE2/Neutralization
        main = "F1 Score by feature type", 
        at = seq(1,6,1), 
        names = c("All features \n (n = 63)", "IgG Serum \n (n = 28)", "IgG Saliva \n (n = 10)", "IgA Saliva \n (n = 10)", "Cytokines  \n (n = 9)", "ACE2/Neut. \n (n = 9)"), 
        col = c("#757476", "#9E3934", "skyblue", "#3675B8", "#F8B959", "#94C652"),
        border = c("#141414", "#D46D61", "#7E9EB2", "#70A0F2", "#AE8E59", "#D7FB92"),
        horizontal = FALSE, 
        notch = FALSE, 
        yaxt = "n",
        boxlwd = 1,
        lwd = 4, 
        las = 2, 
        cex.axis = 1.0,
        ylab = "F1 Score",
        xlab = "",
        ylim = c(0,1)
)
axis(2,las= 1, tck=0.02,mgp=c(3, .5, 0), cex.axis = 1.5)


####################################
# Compute mean ROC curve and plot 1 STD around the mean

####################################

kfoldcrossvalidation <- function(k, features, health_outcome, pid)
{
  w=matrix(NA, nrow = k, ncol = ncol(features))
  pred=vector()
  folds <- cut(pid, breaks=k, labels=FALSE) 
  for (i in 1:k) {
    test_indices <- which(folds == i, arr.ind=TRUE)
    test_data <- features[test_indices, ]
    train_data <- features[-test_indices, ]
    train_labels <- health_outcome[-test_indices]
    rf_model <- randomForest(x = train_data, y = as.factor(train_labels))
    w[i,]=importance(rf_model)
    pred[test_indices] <- predict(rf_model, newdata = test_data, type='prob')[,1]
  }
  return(list(pred, w))
}

kfoldcrossvalidation_externalValidation <- function(k, features, health_outcome, pid, featuresEV)
{
  w=matrix(NA, nrow = k, ncol = ncol(features))
  pred=vector()
  pred.externalValidation = vector()
  folds <- cut(pid, breaks=k, labels=FALSE) 
  for (i in 1:k) {
    test_indices <- which(folds == i, arr.ind=TRUE)
    test_data <- features[test_indices, ]
    train_data <- features[-test_indices, ]
    train_labels <- health_outcome[-test_indices]
    rf_model <- randomForest(x = train_data, y = as.factor(train_labels))
    w[i,]=importance(rf_model)
    #k-fold prediction and then prediction against external validation data set
    pred[test_indices] <- predict(rf_model, newdata = test_data, type='prob')[,1]
    pred.externalValidation = predict(rf_model, newdata = featuresEV, type = 'prob')[,1]
  }
  return(list(pred, pred.externalValidation, w))
}



StatusFile = read.csv("Status_Mario_Brockman.csv", sep = ",", header = T)

health_outcome = StatusFile$HIVStatus[1:91]
health_outcome[which(health_outcome == 2)] = 1
health_outcome[which(health_outcome == 3)] = 1
health_outcome[which(health_outcome == 4)] = 1

features_all.imp = read.csv("features_withMissForestImputation_IR_INR_together_jul25_63Features.csv")
length(features_all.imp[1,])


'%!in%' <- function(x,y)!('%in%'(x,y))
x = features_all.imp[1:91,2:64]
y = health_outcome

keep = SeqOrder[59:63]
x.forwardAb = x[,keep]
x = x.forwardAb

x = features_all.imp[1:91,2:64]
y = health_outcome

x = x[, seq.IgG_serum]


x = features_all.imp[1:91,2:64]
y = health_outcome

x = x[, seq.neut_ACE2]


x = features_all.imp[1:91,2:64]
y = health_outcome

x = x[, seq.Cytokine]


x = features_all.imp[1:91,2:64]
y = health_outcome

x = x[, seq.IgA_saliva]



############ synthetic Data Generate below ############
#x = features_all.imp.MVsynthetic
#x = features_all.imp.COPULAsynthetic
############
#Compute mean ROC curve for all 63 features, and do it again for the top 9, which apparently is nearly equivalent. 
n_individuals = length(y)

Limit = 1000

#Above is initialized with NA, and a cell will only be filled if that individual was in the testing set 
auc=vector()
auc.KV = vector()
roc_data <- list()
roc_data.KV = list()
for(i in seq(Limit))
{
  
  #20:20 randomly sampled HIV-:HIV+
  KeepHIVIndices = sample(24:91, size = 20)
  KeepHIVNegativeIndices = sample(1:23, size = 10)
  indexes = c(KeepHIVNegativeIndices, KeepHIVIndices)
  y.downSamp.TRAIN = health_outcome[indexes]
  x.downSamp.TRAIN = x[indexes,]
  
  #For test: Keep all 3 remaining HIV- indiciduals. Sample 12 remaining HIV+ individuals
  seqHIV.Neg.Test = which(seq(1,23,1) %!in%  KeepHIVNegativeIndices)
  seqHIV.Pos.Test = which(seq(24,91,1) %!in%  KeepHIVIndices)
  HIV.Pos.Test.finalSeq = sample(seq(24,91,1)[seqHIV.Pos.Test], size = 40)

  indexes.TEST = c(seqHIV.Neg.Test, HIV.Pos.Test.finalSeq)
  #Check to make sure indices in TEST are not in TRAIN
  #which(indexes %in% indexes.TEST)
  #which(HIV.Neg.Test.finalSeq %in% KeepHIVIndices)
  y.TEST = health_outcome[indexes.TEST]
  x.TEST = x[indexes.TEST,]
  
  # Train the Random Forest model
  rf_model <- randomForest(x = x.downSamp.TRAIN, y = as.factor(y.downSamp.TRAIN), ntree = 100)
  predictions.prob = predict(rf_model,  newdata = x.TEST, type='prob')
  auc[i]=roc(as.factor(y.TEST), predictions.prob[,1], quiet = TRUE)$auc
  
  #plot(roc(as.factor(y.TEST), predictions.prob[,1], quiet = TRUE))

  # Compute ROC curve for this iteration
  roc_iter <- roc(as.factor(y.TEST), predictions.prob[,1],  thresholds = seq(0, 1, length.out = 100), quiet = TRUE)
  # Compute ROC curve for this iteration with 100 thresholds
  #roc_iter <- roc(test_labels, rf_predictions, thresholds = seq(0, 1, length.out = 100), quiet = TRUE)
  
  #plot( 1-roc_iter$specificities,roc_iter$sensitivities, type = "l")
  
  # Store FPR and TPR values, including (0, 0) and (1, 1)
  #fpr <- c(0, roc_iter$specificities, 1)
  #tpr <- c(0, roc_iter$sensitivities, 1)
  
  roc_data[[i]] <- data.frame(fpr = 1-roc_iter$specificities, tpr = roc_iter$sensitivities)
  
  
  ######With k-fold CV ###
  PID = sample(seq(20+10)) 
  
  out =  kfoldcrossvalidation_externalValidation(5, x.downSamp.TRAIN , y.downSamp.TRAIN, PID, x.TEST)
  

  #out=kfoldcrossvalidation(5, x.downSamp.TRAIN , y.downSamp.TRAIN, PID) #Old KCV function
  #auc.KV[i]=roc(as.factor(y.downSamp.TRAIN), out[[1]], quiet = TRUE)$auc
  auc.KV[i] = roc(as.factor(y.TEST), out[[2]], quiet = TRUE)$auc
  roc_iter.KV = roc(as.factor(y.TEST), out[[2]],  thresholds = seq(0, 1, length.out = 100), quiet = TRUE)
  #roc_iter.KV = roc(as.factor(y.downSamp.TRAIN), out[[1]],  thresholds = seq(0, 1, length.out = 100), quiet = TRUE)
  
  # Extract and store FPR and TPR values, explicitly adding (0, 0) and (1, 1)
  fpr <- c(0, roc_iter.KV$specificities, 1)
  tpr <- c(0, roc_iter.KV$sensitivities, 1)
  
  
  #roc_data.KV[[i]] <- data.frame(fpr = 1-fpr, tpr = tpr)
  roc_data.KV[[i]] = data.frame(fpr = rev(1 - roc_iter.KV$specificities), tpr = rev(roc_iter.KV$sensitivities))

  ########################

  if(i %% 100 == 0){
    print(i)
  }   
  
}
#Might be to do the ROC curves manually with all the thresholds. When we do prediction, command 
quartz()
vioplot(auc, ylab = "AUC")
vioplot(auc.KV, ylab = "AUC with KV")

#What does it mean that HIV+ IRs have a much better cytokine response? More prone to inflamation? Challanged with SARS-CoV-2, better outcomes? 

# Interpolate all ROC curves to a common set of FPR points (including 0 and 1)
fpr_common <- seq(0, 1, length.out = 1000)
tpr_all <- sapply(roc_data.KV, function(roc_curve) {
  approx(roc_curve$fpr, roc_curve$tpr, xout = fpr_common, rule = 2, ties = mean)$y  # Linear interpolation
})

# Calculate the mean and standard deviation of TPR at each common FPR point
tpr_mean <- rowMeans(tpr_all, na.rm = TRUE)
tpr_sd <- apply(tpr_all, 1, sd, na.rm = TRUE)

# Ensure the mean TPR curve includes (0, 0) and (1, 1)
tpr_mean[1] <- 0  # Make sure it starts at (0, 0)
tpr_mean[length(tpr_mean)] <- 1  # Make sure it ends at (1, 1)

# Ensure the std dev curve also follows these rules
tpr_sd[1] <- 0  # Std dev starts at 0 at (0, 0)
tpr_sd[length(tpr_sd)] <- 0  # Std dev is 0 at (1, 1)

# Clip the shaded region to ensure it stays between 0 and 1
tpr_upper <- pmin(tpr_mean + tpr_sd, 1)  # Clip values above 1 to 1
tpr_lower <- pmax(tpr_mean - tpr_sd, 0)  # Clip values below 0 to 0

# Plot all ROC curves
quartz()
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), type = "n", xlab = "False Positive Rate (FPR)", 
     ylab = "True Positive Rate (TPR)", main = "ROC Curves (100 Iterations)")
for (i in 1:50) {
  lines(roc_data.KV[[i]]$fpr, roc_data.KV[[i]]$tpr, col = rgb(0.1, 0.1, 0.8, 0.2), lwd = 1)
}
abline(a = 0, b = 1, lty = 2, col = "gray")  # Random guessing line

# Plot the mean ROC curve
lines(fpr_common, tpr_mean, col = "red", lwd = 2)


# Plot the shaded region representing ±1 std dev around the mean ROC curve
polygon(c(fpr_common, rev(fpr_common)),
        c(tpr_mean + tpr_sd, rev(tpr_mean - tpr_sd)),
        col = rgb(0.1, 0.1, 0.8, 0.2), border = NA) 


quartz()

plot(fpr_common, tpr_mean, col = "#924038", lwd = 2, type = "l", ylim = c(0,1), xlim = c(0,1))

# Plot the shaded region representing ±1 std dev around the mean ROC curve
polygon(c(fpr_common, rev(fpr_common)),
        c(tpr_upper, rev(tpr_lower)),
        col = rgb(205/255, 88/255, 79/255, 0.2), border = NA) 

abline(a = 0, b = 1, lty = 2, col = "black")  # 


IgASaliva_fpr_common = fpr_common
IgASaliva_tpr_mean = tpr_mean
IgASaliva_std_upper = tpr_upper
IgASaliva_std_lower = tpr_lower

Cytokine_fpr_common = fpr_common
Cytokine_tpr_mean = tpr_mean
Cytokine_std_upper = tpr_upper
Cytokine_std_lower = tpr_lower


Neut_fpr_common = fpr_common
Neut_tpr_mean = tpr_mean
Neut_std_upper = tpr_upper
Neut_std_lower = tpr_lower

serum_fpr_common = fpr_common
serum_tpr_mean = tpr_mean
serum_std_upper = tpr_upper
serum_std_lower = tpr_lower

quartz()
plot(serum_fpr_common, serum_tpr_mean, col = "#924038", lwd = 5, type = "l", ylim = c(0,1), xlim = c(0,1), xaxt = "n", yaxt = "n")
axis(2,las= 1, tck=0.02,mgp=c(3, .5, 0), cex.axis = 1.9)
axis(1,las= 1, tck=0.02,mgp=c(3, .5, 0), cex.axis = 1.9)

# Plot the shaded region representing ±1 std dev around the mean ROC curve
polygon(c(serum_fpr_common, rev(serum_fpr_common)),
        c(serum_std_upper, rev(serum_std_lower)),
        col = rgb(205/255, 88/255, 79/255, 0.2), border = NA) 


# lines(Neut_fpr_common, Neut_tpr_mean, col = "#9EC562", lwd = 2, type = "l", ylim = c(0,1), xlim = c(0,1))
# polygon(c(Neut_fpr_common, rev(Neut_fpr_common)),
#         c(Neut_std_upper, rev(Neut_std_lower)),
#         col = rgb(158/255, 197/255, 98/255, 0.2), border = NA) 


lines(IgASaliva_fpr_common, IgASaliva_tpr_mean, col = "#4674B3", lwd = 5, type = "l", ylim = c(0,1), xlim = c(0,1))
polygon(c(IgASaliva_fpr_common, rev(IgASaliva_fpr_common)),
        c(IgASaliva_std_upper, rev(IgASaliva_std_lower)),
        col = rgb(70/255, 116/255, 79/255, 0.2), border = NA) 

lines(Cytokine_fpr_common, Cytokine_tpr_mean, col = "#EEBC69", lwd = 5, type = "l", ylim = c(0,1), xlim = c(0,1))
polygon(c(Cytokine_fpr_common, rev(Cytokine_fpr_common)),
        c(Cytokine_std_upper, rev(Cytokine_std_lower)),
        col = rgb(238/255, 188/255, 105/255, 0.2), border = NA) 

abline(a = 0, b = 1, lty = 2, col = "black", lwd = 5)  # 


save(serum_fpr_common, serum_tpr_mean,serum_std_upper ,serum_std_lower, file = "meanROC_serumData.RData")

save(Neut_fpr_common, Neut_tpr_mean,Neut_std_upper ,Neut_std_lower, file = "meanROC_NeutData.RData")

save(Cytokine_fpr_common, Cytokine_tpr_mean, Cytokine_std_upper ,Cytokine_std_lower, file = "meanROC_CytokineData.RData")

save(IgASaliva_fpr_common, IgASaliva_tpr_mean, IgASaliva_std_upper ,IgASaliva_std_lower, file = "meanROC_IgAData.RData")

i <- 115  # Change this to any index you want to plot

# Extract FPR and TPR from the i-th ROC curve
fpr_single <- roc_data.KV[[i]]$fpr
tpr_single <- roc_data.KV[[i]]$tpr
quartz()
# Plot the single ROC curve
plot(1- fpr_single, tpr_single, type = "l", col = "red", lwd = 2,
     xlab = "False Positive Rate (1 - Sp)", ylab = "True Positive Rate (Se)", 
     main = paste("ROC Curve for Iteration", i),
     xlim = c(0, 1), ylim = c(0, 1))

for(i in 1:100){
  fpr_single <- roc_data.KV[[i]]$fpr
  tpr_single <- roc_data.KV[[i]]$tpr
  lines(1- fpr_single, tpr_single, col = "black")
}
# Optionally, you can add a diagonal line to represent random guessing
abline(a = 0, b = 1, lty = 2, col = "gray")  # Random guess line (diagonal)
roc_data.KV[[1]]$fpr

