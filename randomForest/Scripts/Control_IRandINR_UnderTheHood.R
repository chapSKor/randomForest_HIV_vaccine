library(randomForest)
library(missForest)
library(pROC)
library(ggplot2)
library(MASS)
library(caret)
library(vioplot)

StatusFile = read.csv("Status_Mario_Brockman.csv", sep = ",", header = T)

health_outcome = StatusFile$HIVStatus[1:91]
health_outcome[which(health_outcome == 2)] = 1
health_outcome[which(health_outcome == 3)] = 1
health_outcome[which(health_outcome == 4)] = 1

#features_all.imp = read.csv("features_all.imp_June17.csv")
#features_all.imp = read.csv("features_withMissForestImputation_IR_INR_separate_July9.csv")
#features_all.imp = read.csv("features_withMissForestImputation_IR_INR_together.csv")
#features_all.imp = read.csv("features_withMissForestImputation_IR_INR_together_jul11_64Features.csv")
#features_all.imp = read.csv("features_withMissForestImputation_IR_INR_together_jul17_64Features.csv")
#features_all.imp = read.csv("features_withMissForestImputation_IR_INR_together_jul18_64Features.csv")
features_all.imp = read.csv("features_withMissForestImputation_IR_INR_together_jul25_63Features.csv")


'%!in%' <- function(x,y)!('%in%'(x,y))
x = features_all.imp[1:91,2:64]
y = health_outcome

############ synthetic Data Generate below ############
#x = features_all.imp.MVsynthetic
#x = features_all.imp.COPULAsynthetic
############

n_individuals = length(y)

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
#Columns are the IDs and rows are the probabilities after each iteration. 
Individual.Predictions = data.frame(matrix(ncol = length(y), nrow = Limit))
#Above is initialized with NA, and a cell will only be filled if that individual was in the testing set 

for(i in seq(Limit))
{
  #20:20 randomly sampled HIV-:HIV+
  KeepHIVIndices = sample(24:91, size = 20)
  KeepHIVNegativeIndices = sample(1:23, size = 20)
  indexes = c(KeepHIVNegativeIndices, KeepHIVIndices)
  y.downSamp.TRAIN = health_outcome[indexes]
  x.downSamp.TRAIN = x[indexes,]
  
  #For test: Keep all 3 remaining HIV- indiciduals. Sample 12 remaining HIV+ individuals
  seqHIV.Neg.Test = which(seq(1,23,1) %!in%  KeepHIVNegativeIndices)
  seqHIV.Pos.Test = which(seq(24,91,1) %!in%  KeepHIVIndices)
  HIV.Neg.Test.finalSeq = sample(seq(24,91,1)[seqHIV.Pos.Test], size = 12)
  indexes.TEST = c(seqHIV.Neg.Test, HIV.Neg.Test.finalSeq)
  #Check to make sure indices in TEST are not in TRAIN
  #which(indexes %in% indexes.TEST)
  #which(HIV.Neg.Test.finalSeq %in% KeepHIVIndices)
  y.TEST = health_outcome[indexes.TEST]
  x.TEST = x[indexes.TEST,]
  
  # Train the Random Forest model
  rf_model <- randomForest(x = x.downSamp.TRAIN, y = as.factor(y.downSamp.TRAIN), ntree = 80)
  
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
  #Store probabilities
  
  #Save probabilities they are classified as HIV+, note that the probabiity of HIV- is 1 - prob(HIV+)
  Individual.Predictions[i,indexes.TEST] = as.numeric(predictions.prob[,2])
  
  
  if(i == 1){
    varIMPORT.cumAVG = varImp(rf_model)
  }else{
    varIMPORT.cumAVG = (varImp(rf_model) + (i-1)*varIMPORT.cumAVG)/(i)
  }
  
  if(i %% 100 == 0){
    print(i)
  }    
}

###### Save #####
write.csv(INCorrectClassification_counter, "INCorrectClassifications_downSampled_control_IRandINR_july18.csv")
write.csv(CorrectClassification_counter, "CorrectClassifications_downSampled_control_IRandINR_july18.csv")
write.csv(varIMPORT.cumAVG$Overall, "ImportanceFactor_downSampled_control_IRandINR_july18.csv")
write.csv(testSetCount, "totalTestCount_downSampled_control_IRandINR_july18.csv")

write.csv(sensitivity, "sensitivity_control_IRandINR_july18.csv")
write.csv(specificity, "specificity_control_IRandINR_july18.csv")
write.csv( precision, " precision_control_IRandINR_july18.csv")
write.csv(f1_score, "f1_score_control_IRandINR_july18.csv")
write.csv(auc, "auc_control_IRandINR_july18.csv")

####################
###### Plots #####
####################

variables.Names <- c('spikeProduction_D1D2', 'spikeDecay_D1D2', 'spikeProduction_D3', 'spikeDecay_D3', 'RBDProduction_D1D2', 'RBDDecay_D1D2', 'RBDProduction_D3', 'RBDDecay_D3',
                     'V1_blood_IgGspike', 
                     'V4_blood_IgGspike', 
                     'V4a_blood_IgGspike', 
                     'V6_blood_IgGspike',
                     'V8_blood_IgGspike', 
                     'V8a_blood_IgGspike', 
                     'V8b_blood_IgGspike', 
                     'V9_blood_IgGspike', 
                     'V10_blood_IgGspike', 
                     'V11_blood_IgGspike', 
                     'V1_blood_IgGRBD',
                     'V4_blood_IgGRBD',
                     'V4a_blood_IgGRBD',
                     'V6_blood_IgGRBD',
                     'V8_blood_IgGRBD',
                     'V8a_blood_IgGRBD',
                     'V8b_blood_IgGRBD',
                     'V9_blood_IgGRBD',
                     'V10_blood_IgGRBD',
                     'V11_blood_IgGRBD',
                     'V4_Saliva_IgGspike', 'V4_Saliva_IgGRBD', 
                     'V5_Saliva_IgGspike', 'V5_Saliva_IgGRBD',
                     'V8_Saliva_IgGspike', 'V8_Saliva_IgGRBD', 
                     'V8b_Saliva_IgGspike', 'V8b_Saliva_IgGRBD', 
                     'V9_Saliva_IgGspike', 'V9_Saliva_IgGRBD', 
                     'V4_Saliva_IgAspike', 'V4_Saliva_IgARBD', 
                     'V5_Saliva_IgAspike', 'V5_Saliva_IgARBD', 
                     'V8_Saliva_IgAspike', 'V8_Saliva_IgARBD', 
                     'V8b_Saliva_IgAspike', 'V8b_Saliva_IgARBD', 
                     'V9_Saliva_IgAspike', 'V9_Saliva_IgARBD', 
                     'V8_IFNg', 'V9_IFNg','V8_IL2', 'V9_IL2','V9Dual', 'V8Dual', 'IFNG_production','Il2_production', 'RATIO_CD4CD8', 
                     'V8Neut', 'V9Neut',
                     'V7_ACE2', 'V8_ACE2', 'V8b_ACE2','V9_ACE2')


quartz()
colors = c(rep("#757476", 23), rep("#00A1FF", length(seq(24,66,1))), rep("#F8B959", length(seq(67,91,1))))
barplot.colours = c(rep("firebrick", 28), rep("skyblue", 20), rep("#F8B959", 9), rep("#757476", 6))
barplot.colours = c(rep("#DC3A7A", 8), rep("#9E3934", 10),rep("#ED220D", 10), rep("skyblue", 10),rep("#3675B8", 10), rep("#F8B959", 9), rep("#94C652", 6))
par(mar=c(11,4,4,4))
barplot(varIMPORT.cumAVG$Overall,names.arg = variables.Names, col=barplot.colours, las=2, ylim = c(0, 6), ylab = "Measure of Importance")
box(bty = "o")
legend("topleft", 
       legend = c("IgG model rates", "Serum - IgG spike", "Serum - IgG RBD", "Saliva - IgG spike/RBD", "Saliva - IgA spike/RBD", "Cytokines", "Neutralization/ACE2"), 
       fill = c("#DC3A7A", "#9E3934", "#ED220D", "skyblue", "#3675B8", "#F8B959", "#94C652"), cex = 1.8)

#Correct classification barplot
quartz()
ID = seq(1,91,1)
par(mar=c(11,4,4,4))
barplot(CorrectClassification_counter/testSetCount, names.arg = ID, col="#69b3a2", las=2, ylim = c(0, 1.1), ylab = "Probability of correct classification")
box(bty = "o")

#Incorrect classification barplot
quartz()
ID = seq(1,91,1)
par(mar=c(11,4,4,4))
barplot(INCorrectClassification_counter/testSetCount,names.arg = ID, col="#69b3a2", las=2, ylim = c(0, 1.1), ylab = "Probability of incorrect classification")
box(bty = "o")

quartz()
par(mfrow = c(2, 2))
boxplot(sensitivity, ylab = "sensitivity") #you can replace x, and y with your real data values (x~ features, y~outcome or response vector) with appropriate x, y axis labels
boxplot(specificity, ylab = "specificity")
boxplot(f1_score, ylab = "F1 Score")
boxplot(auc, ylab = "AUC")
mean(auc)
(mean(na.omit(f1_score)))
mean(specificity)
mean(sensitivity)
#Basic auc plot
dev.off()

quartz()
vioplot(auc, ylab = "AUC")
vioplot(f1_score, ylab = "F1 Score")
##In large synthetic data limit does one of the clinical features for {IR,INR} emerge as vary important. This could lead to a prediction for the most clinically important feature in distinguishing these two classes.

##Violin plots of the individual probabilties 
#Note these are the probabilities the individuals are determined to be HIV+ SARS-CoV-2 vaccine responders. 
#1:23 are HIV-, 24:66 are IR, 67:90 are INR
length(seq(24,66,1))
colors = c(rep("#757476", 23), rep("#00A1FF", length(seq(24,66,1))), rep("#F8B959", length(seq(67,91,1))))
colorsBORDER = c(rep("#141414", 23), rep("#22548C", length(seq(24,66,1))), rep("#B18545", length(seq(67,91,1))))
quartz()
vioplot(na.omit(Individual.Predictions[,1]), 
        na.omit(Individual.Predictions[,2]), 
        na.omit(Individual.Predictions[,3]), 
        na.omit(Individual.Predictions[,4]), 
        na.omit(Individual.Predictions[,5]), 
        na.omit(Individual.Predictions[,6]), 
        na.omit(Individual.Predictions[,7]), 
        na.omit(Individual.Predictions[,8]), 
        na.omit(Individual.Predictions[,9]), 
        na.omit(Individual.Predictions[,10]), 
        na.omit(Individual.Predictions[,11]), 
        na.omit(Individual.Predictions[,12]), 
        na.omit(Individual.Predictions[,13]), 
        na.omit(Individual.Predictions[,14]), 
        na.omit(Individual.Predictions[,15]), 
        na.omit(Individual.Predictions[,16]), 
        na.omit(Individual.Predictions[,17]), 
        na.omit(Individual.Predictions[,18]), 
        na.omit(Individual.Predictions[,19]), 
        na.omit(Individual.Predictions[,20]), 
        na.omit(Individual.Predictions[,21]), 
        na.omit(Individual.Predictions[,22]), 
        na.omit(Individual.Predictions[,23]), 
        na.omit(Individual.Predictions[,24]), 
        na.omit(Individual.Predictions[,25]), 
        na.omit(Individual.Predictions[,26]), 
        na.omit(Individual.Predictions[,27]), 
        na.omit(Individual.Predictions[,28]), 
        na.omit(Individual.Predictions[,29]), 
        na.omit(Individual.Predictions[,30]), 
        na.omit(Individual.Predictions[,31]), 
        na.omit(Individual.Predictions[,32]), 
        na.omit(Individual.Predictions[,33]), 
        na.omit(Individual.Predictions[,34]), 
        na.omit(Individual.Predictions[,35]), 
        na.omit(Individual.Predictions[,36]), 
        na.omit(Individual.Predictions[,37]), 
        na.omit(Individual.Predictions[,38]), 
        na.omit(Individual.Predictions[,39]), 
        na.omit(Individual.Predictions[,40]), 
        na.omit(Individual.Predictions[,41]), 
        na.omit(Individual.Predictions[,42]), 
        na.omit(Individual.Predictions[,43]), 
        na.omit(Individual.Predictions[,44]), 
        na.omit(Individual.Predictions[,45]), 
        na.omit(Individual.Predictions[,46]), 
        na.omit(Individual.Predictions[,47]), 
        na.omit(Individual.Predictions[,48]), 
        na.omit(Individual.Predictions[,49]), 
        na.omit(Individual.Predictions[,50]), 
        na.omit(Individual.Predictions[,51]), 
        na.omit(Individual.Predictions[,52]), 
        na.omit(Individual.Predictions[,53]), 
        na.omit(Individual.Predictions[,54]), 
        na.omit(Individual.Predictions[,55]), 
        na.omit(Individual.Predictions[,56]), 
        na.omit(Individual.Predictions[,57]), 
        na.omit(Individual.Predictions[,58]), 
        na.omit(Individual.Predictions[,59]), 
        na.omit(Individual.Predictions[,60]), 
        na.omit(Individual.Predictions[,61]), 
        na.omit(Individual.Predictions[,62]), 
        na.omit(Individual.Predictions[,63]), 
        na.omit(Individual.Predictions[,64]), 
        na.omit(Individual.Predictions[,65]), 
        na.omit(Individual.Predictions[,66]), 
        na.omit(Individual.Predictions[,67]), 
        na.omit(Individual.Predictions[,68]), 
        na.omit(Individual.Predictions[,69]), 
        na.omit(Individual.Predictions[,70]), 
        na.omit(Individual.Predictions[,71]), 
        na.omit(Individual.Predictions[,72]), 
        na.omit(Individual.Predictions[,73]), 
        na.omit(Individual.Predictions[,74]), 
        na.omit(Individual.Predictions[,75]), 
        na.omit(Individual.Predictions[,76]), 
        na.omit(Individual.Predictions[,77]), 
        na.omit(Individual.Predictions[,78]), 
        na.omit(Individual.Predictions[,79]), 
        na.omit(Individual.Predictions[,80]), 
        na.omit(Individual.Predictions[,81]), 
        na.omit(Individual.Predictions[,82]), 
        na.omit(Individual.Predictions[,83]), 
        na.omit(Individual.Predictions[,84]), 
        na.omit(Individual.Predictions[,85]), 
        na.omit(Individual.Predictions[,86]), 
        na.omit(Individual.Predictions[,87]), 
        na.omit(Individual.Predictions[,88]), 
        na.omit(Individual.Predictions[,89]), 
        na.omit(Individual.Predictions[,90]), 
        na.omit(Individual.Predictions[,91]), 
        main = "Prob predict HIV+",
        ylim = c(0,1), 
        col = colors,
        border =colorsBORDER,  
        lwd = 2, 
        cex.axis = 2.0
)
abline(h = 0.5, col = "red", lty = 2, lwd = 3)

#Below red line means classification is HIV-, above red line means classification is HIV+ 

####################
###### Look under the hood #####
####################

#Stop the RF cost at a high AUC RF outcome. Then use the lime package to "look under the hood" of the predictions and feature importance for each individual of that iteration. 
#Can check training and predicting separately. 

library(lime)

#Stop when the AUC is high 
n_individuals = length(y)
Limit = 1000
auc=vector()

for(i in seq(Limit))
{
  #20:20 randomly sampled HIV-:HIV+
  KeepHIVIndices = sample(24:91, size = 20)
  KeepHIVNegativeIndices = sample(1:23, size = 20)
  indexes = c(KeepHIVNegativeIndices, KeepHIVIndices)
  y.downSamp.TRAIN = y[indexes]
  x.downSamp.TRAIN = x[indexes,]
  
  #For test: Keep all 3 remaining HIV- indiciduals. Sample 12 remaining HIV+ individuals
  seqHIV.Neg.Test = which(seq(1,23,1) %!in%  KeepHIVNegativeIndices)
  seqHIV.Pos.Test = which(seq(24,91,1) %!in%  KeepHIVIndices)
  HIV.Neg.Test.finalSeq = sample(seq(24,91,1)[seqHIV.Pos.Test], size = 12)
  indexes.TEST = c(seqHIV.Neg.Test, HIV.Neg.Test.finalSeq)
  
  #Check to make sure indices in TEST are not in TRAIN
  #which(indexes %in% indexes.TEST)
  #which(HIV.Neg.Test.finalSeq %in% KeepHIVIndices)
  y.TEST = y[indexes.TEST]
  x.TEST = x[indexes.TEST,]
  
  # Train the Random Forest model
  rf_model <- randomForest(x = x.downSamp.TRAIN, y = as.factor(y.downSamp.TRAIN), ntree = 64)
  
  #Keep track of a counter of how many times each ID is put into the test group
  testSetCount[indexes.TEST] = testSetCount[indexes.TEST] + 1
  #Predict for Accuracy
  predictions <- predict(rf_model, newdata = x.TEST)
  accuracy <- sum(predictions == y.TEST) / length(y.TEST)
  #Predict for probabilities
  predictions.prob = predict(rf_model,  newdata = x.TEST, type='prob')
  auc[i]=roc(as.factor(y.TEST), predictions.prob[,1], quiet = TRUE)$auc
  if(auc[i] > 0.98){
    break #Take a high auc RF model example
  }

}

y.downSamp.TRAIN.FACTOR = as.factor(y.downSamp.TRAIN)
x.downSamp.TRAIN$outcome = y.downSamp.TRAIN

model_caret <- train(outcome ~ ., data = as.data.frame(x.downSamp.TRAIN), method = "rf", tuneLength = 4)
# Prepare the explainer with the model wrapped by caret
explainer <- lime(as.data.frame(x.downSamp.TRAIN[, -ncol(x.downSamp.TRAIN)]), model = model_caret, 
                  bin_continuous = TRUE, predict_function = predict_proba, type = "classification")

# Plot the explanation
# Plot for Testing Data set
LIM =15
explanation_2 <- explain(x.TEST[1:5,], explainer, n_labels = 1, n_features = 63)
quartz()
plot_explanations(explanation_2)
plot_features(explanation_2)


#Plot for training Data set
explanation <- explain(as.data.frame(x.downSamp.TRAIN[, -ncol(x.downSamp.TRAIN)]), explainer, n_labels = 1, n_features = 63)
quartz()
plot_explanations(explanation, cex.axis = 0.5)


explanation_2 <- explain(as.data.frame(x.downSamp.TRAIN[25, -ncol(x.downSamp.TRAIN)]), explainer, n_labels = 1, n_features = 63)
quartz()
plot_features(explanation_2)

# Assuming `plot_explanations` returns a ggplot object
plot <- plot_explanations(explanation)

# Customize y-axis text size
plot + theme(axis.text.y = element_text(size = 4))

#ID 40 for HIV+ and ID 12 for control seem interesting. 
#OKay I want to run the training script, but every time ID 39 gets put into the test group I want to save their test landscape. 
#Maybe just keep running the script until there are 100 ID 39 tests ?

seqID39 = 3
#Only run model_cart and the explainer if ID 39 is randomly chosen. 
x.TEST[seqID39,]

explanation_2 <- explain(x.TEST, explainer, n_features = 63)

quartz()
plot_features(explanation_2)
#Cases can be used to isolate features weights corresponding to a specific ID. They are in the same order as x.TEST.
explanation_2$case
#For each person I need to make a data column that is N = 63 rows with m columns, where each column is a particular ith iteration. 


library(lime)
##Similar code as above, however, I want to first check to see if ID 39 is in the test group, if they are, then I want to keep going and then save their explanation. If they are not, then I 
#want to go to the next iteraction. Stop when we have 10 ID 39 explanation results. 
predict_proba <- function(model, newdata, type = 'prob') {
  predictions <- predict(model, newdata, type = type)
  data.frame(Healthy = predictions[, "0"], Unhealthy = predictions[, "1"])
}
#rwos are iterations, columns are the features values
underTheHood_ID39 = data.frame(matrix(ncol = 63, nrow = 20))
underTheHood_ID43 = data.frame(matrix(ncol = 63, nrow = 20))

colnames(underTheHood_ID43) <- c('spikeProduction_D1D2', 'spikeDecay_D1D2', 'spikeProduction_D3', 'spikeDecay_D3', 'RBDProduction_D1D2', 'RBDDecay_D1D2', 'RBDProduction_D3', 'RBDDecay_D3',
                            'V1_blood_IgGspike', 
                            'V4_blood_IgGspike', 
                            'V4a_blood_IgGspike', 
                            'V6_blood_IgGspike',
                            'V8_blood_IgGspike', 
                            'V8a_blood_IgGspike', 
                            'V8b_blood_IgGspike', 
                            'V9_blood_IgGspike', 
                            'V10_blood_IgGspike', 
                            'V11_blood_IgGspike', 
                            'V1_blood_IgGRBD',
                            'V4_blood_IgGRBD',
                            'V4a_blood_IgGRBD',
                            'V6_blood_IgGRBD',
                            'V8_blood_IgGRBD',
                            'V8a_blood_IgGRBD',
                            'V8b_blood_IgGRBD',
                            'V9_blood_IgGRBD',
                            'V10_blood_IgGRBD',
                            'V11_blood_IgGRBD',
                            'V4_Saliva_IgGspike', 'V4_Saliva_IgGRBD', 
                            'V5_Saliva_IgGspike', 'V5_Saliva_IgGRBD',
                            'V8_Saliva_IgGspike', 'V8_Saliva_IgGRBD', 
                            'V8b_Saliva_IgGspike', 'V8b_Saliva_IgGRBD', 
                            'V9_Saliva_IgGspike', 'V9_Saliva_IgGRBD', 
                            'V4_Saliva_IgAspike', 'V4_Saliva_IgARBD', 
                            'V5_Saliva_IgAspike', 'V5_Saliva_IgARBD', 
                            'V8_Saliva_IgAspike', 'V8_Saliva_IgARBD', 
                            'V8b_Saliva_IgAspike', 'V8b_Saliva_IgARBD', 
                            'V9_Saliva_IgAspike', 'V9_Saliva_IgARBD', 
                            'V8_IFNg', 'V9_IFNg','V8_IL2', 'V9_IL2','V9Dual', 'V8Dual', 'IFNG_production','Il2_production', 'RATIO_CD4CD8', 
                            'V8Neut', 'V9Neut',
                            'V7_ACE2', 'V8_ACE2', 'V8b_ACE2','V9_ACE2' 
)

n_individuals = length(y)
ID39_predictions = numeric(20)
ID43_predictions = numeric(20)
Limit = 1000
auc=vector()
counter = 1
for(i in seq(Limit))
{
  #20:20 randomly sampled HIV-:HIV+
  KeepHIVIndices = sample(24:91, size = 20)
  KeepHIVNegativeIndices = sample(1:23, size = 20)
  indexes = c(KeepHIVNegativeIndices, KeepHIVIndices)
  y.downSamp.TRAIN = y[indexes]
  x.downSamp.TRAIN = x[indexes,]
  
  #For test: Keep all 3 remaining HIV- indiciduals. Sample 12 remaining HIV+ individuals
  seqHIV.Neg.Test = which(seq(1,23,1) %!in%  KeepHIVNegativeIndices)
  seqHIV.Pos.Test = which(seq(24,91,1) %!in%  KeepHIVIndices)
  HIV.Neg.Test.finalSeq = sample(seq(24,91,1)[seqHIV.Pos.Test], size = 12)
  indexes.TEST = c(seqHIV.Neg.Test, HIV.Neg.Test.finalSeq)
  #Check to make sure indices in TEST are not in TRAIN
  #which(indexes %in% indexes.TEST)
  #which(HIV.Neg.Test.finalSeq %in% KeepHIVIndices)
  y.TEST = y[indexes.TEST]
  x.TEST = x[indexes.TEST,]
  
  temp = which(indexes.TEST == 43)
  check = length(temp)
  if(check == 0){
    next
  }else{
    
  #if this is selected then ID 39 corresponds on the temp index. 
  # Train the Random Forest model
  rf_model <- randomForest(x = x.downSamp.TRAIN, y = as.factor(y.downSamp.TRAIN), ntree = 100)
  
  y.downSamp.TRAIN.FACTOR = as.factor(y.downSamp.TRAIN)
  x.downSamp.TRAIN$outcome = y.downSamp.TRAIN
  model_caret <- train(outcome ~ ., data = as.data.frame(x.downSamp.TRAIN), method = "rf", tuneLength = 4)
  explainer <- lime(as.data.frame(x.downSamp.TRAIN[, -ncol(x.downSamp.TRAIN)]), model = model_caret, 
                    bin_continuous = TRUE, predict_function = predict_proba, type = "classification")
  
  #Only run model_cart and the explainer if ID 39 is randomly chosen. 
  explanation_2 <- explain(x.TEST[temp,], explainer, n_features = 63)
  
  #underTheHood_ID39[counter,] = explanation_2$feature_weight
  underTheHood_ID43[counter,] = explanation_2$feature_weight
  #Keep track of a counter of how many times each ID is put into the test group
  #testSetCount[indexes.TEST] = testSetCount[indexes.TEST] + 1
  #Predict for Accuracy
  predictions <- predict(model_caret, newdata = x.TEST)
  predictions.prob <- predict(rf_model, newdata = x.TEST, type = 'prob')
  #Save probility of being assigned HIV+ (over 0.5 means its a 1)

  #ID39_predictions[counter] = predictions.prob[,2][temp]
  ID43_predictions[counter] = predictions.prob[,2][temp]
  
  #accuracy <- sum(predictions == y.TEST) / length(y.TEST)
  #Predict for probabilities
  auc[counter]=roc(as.factor(y.TEST), predictions.prob[,1], quiet = TRUE)$auc
  
  counter = counter + 1
  print(counter)
  if(counter == 21){
    break
  }
  }
}


feature.labels <- c('spikeProduction_D1D2', 'spikeDecay_D1D2', 'spikeProduction_D3', 'spikeDecay_D3', 'RBDProduction_D1D2', 'RBDDecay_D1D2', 'RBDProduction_D3', 'RBDDecay_D3',
                                 'V1_blood_IgGspike', 
                                 'V4_blood_IgGspike', 
                                 'V4a_blood_IgGspike', 
                                 'V6_blood_IgGspike',
                                 'V8_blood_IgGspike', 
                                 'V8a_blood_IgGspike', 
                                 'V8b_blood_IgGspike', 
                                 'V9_blood_IgGspike', 
                                 'V10_blood_IgGspike', 
                                 'V11_blood_IgGspike', 
                                 'V1_blood_IgGRBD',
                                 'V4_blood_IgGRBD',
                                 'V4a_blood_IgGRBD',
                                 'V6_blood_IgGRBD',
                                 'V8_blood_IgGRBD',
                                 'V8a_blood_IgGRBD',
                                 'V8b_blood_IgGRBD',
                                 'V9_blood_IgGRBD',
                                 'V10_blood_IgGRBD',
                                 'V11_blood_IgGRBD',
                                 'V4_Saliva_IgGspike', 'V4_Saliva_IgGRBD', 
                                 'V5_Saliva_IgGspike', 'V5_Saliva_IgGRBD',
                                 'V8_Saliva_IgGspike', 'V8_Saliva_IgGRBD', 
                                 'V8b_Saliva_IgGspike', 'V8b_Saliva_IgGRBD', 
                                 'V9_Saliva_IgGspike', 'V9_Saliva_IgGRBD', 
                                 'V4_Saliva_IgAspike', 'V4_Saliva_IgARBD', 
                                 'V5_Saliva_IgAspike', 'V5_Saliva_IgARBD', 
                                 'V8_Saliva_IgAspike', 'V8_Saliva_IgARBD', 
                                 'V8b_Saliva_IgAspike', 'V8b_Saliva_IgARBD', 
                                 'V9_Saliva_IgAspike', 'V9_Saliva_IgARBD', 
                                 'V8_IFNg', 'V9_IFNg','V8_IL2', 'V9_IL2','V9Dual', 'V8Dual', 'IFNG_production','Il2_production', 'RATIO_CD4CD8', 
                                 'V8Neut', 'V9Neut',
                                 'V7_ACE2', 'V8_ACE2', 'V8b_ACE2','V9_ACE2' 
)
feature.labels[57]

length(seq(24,66,1))
colors = c(rep("#DC3A7A", 8), rep("#9E3934", 10),rep("#ED220D", 10), rep("skyblue", 10),rep("#3675B8", 10), rep("#F8B959", 9), rep("#94C652", 6))
colorsBORDER = c(rep("#141414", 23), rep("#22548C", length(seq(24,66,1))), rep("#B18545", length(seq(67,91,1))))
data.plot = underTheHood_ID43
quartz()
vioplot(na.omit(data.plot[,1]), 
        na.omit(data.plot[,2]), 
        na.omit(data.plot[,3]), 
        na.omit(data.plot[,4]), 
        na.omit(data.plot[,5]), 
        na.omit(data.plot[,6]), 
        na.omit(data.plot[,7]), 
        na.omit(data.plot[,8]), 
        na.omit(data.plot[,9]), 
        na.omit(data.plot[,10]), 
        na.omit(data.plot[,11]), 
        na.omit(data.plot[,12]), 
        na.omit(data.plot[,13]), 
        na.omit(data.plot[,14]), 
        na.omit(data.plot[,15]), 
        na.omit(data.plot[,16]), 
        na.omit(data.plot[,17]), 
        na.omit(data.plot[,18]), 
        na.omit(data.plot[,19]), 
        na.omit(data.plot[,20]), 
        na.omit(data.plot[,21]), 
        na.omit(data.plot[,22]), 
        na.omit(data.plot[,23]), 
        na.omit(data.plot[,24]), 
        na.omit(data.plot[,25]), 
        na.omit(data.plot[,26]), 
        na.omit(data.plot[,27]), 
        na.omit(data.plot[,28]), 
        na.omit(data.plot[,29]), 
        na.omit(data.plot[,30]), 
        na.omit(data.plot[,31]), 
        na.omit(data.plot[,32]), 
        na.omit(data.plot[,33]), 
        na.omit(data.plot[,34]), 
        na.omit(data.plot[,35]), 
        na.omit(data.plot[,36]), 
        na.omit(data.plot[,37]), 
        na.omit(data.plot[,38]), 
        na.omit(data.plot[,39]), 
        na.omit(data.plot[,40]), 
        na.omit(data.plot[,41]), 
        na.omit(data.plot[,42]), 
        na.omit(data.plot[,43]), 
        na.omit(data.plot[,44]), 
        na.omit(data.plot[,45]), 
        na.omit(data.plot[,46]), 
        na.omit(data.plot[,47]), 
        na.omit(data.plot[,48]), 
        na.omit(data.plot[,49]), 
        na.omit(data.plot[,50]), 
        na.omit(data.plot[,51]), 
        na.omit(data.plot[,52]), 
        na.omit(data.plot[,53]), 
        na.omit(data.plot[,54]), 
        na.omit(data.plot[,55]), 
        na.omit(data.plot[,56]), 
        na.omit(data.plot[,57]), 
        na.omit(data.plot[,58]), 
        na.omit(data.plot[,59]), 
        na.omit(data.plot[,60]), 
        na.omit(data.plot[,61]), 
        na.omit(data.plot[,62]), 
        na.omit(data.plot[,63]), 
        main = "ID 43 feature weight",
        ylim = c(min(data.plot),max(data.plot)), 
        lwd = 2, 
        col = colors,
        cex.axis = 1.0, 
        names = feature.labels, 
        las = 2
)

#Okay now go through and do this for every individual. So far it looks like the wrongly classed individuals tend to have cytokines that look like the opposing class values. 

# Define the dimensions
num_ids <- 91
num_features <- 63
num_iterations <- 50

# Create the 3-dimensional array with the specified dimensions
# Initialize the array with NA or 0 or any other value you want
my_array <- array(NA, dim = c(num_ids, num_features, num_iterations))

# Optionally, you can name the dimensions for better clarity
dimnames(my_array) <- list(
  ID = paste0("", 1:num_ids),
  Feature = feature.labels,
  Iteration = paste0("i_", 1:num_iterations)
)

# Print the structure of the array to verify
str(my_array)

#my_array[id, feature, iteration]
#Save the values of the weights for those IDs who are in the testing group. 


Limit = 1000
auc=vector()
counter = 1
testSetCount = numeric(n_individuals)
Individual.Predictions = data.frame(matrix(ncol = length(y), nrow = num_iterations))
for(i in 1:Limit)
{
  #20:20 randomly sampled HIV-:HIV+
  KeepHIVIndices = sample(24:91, size = 20)
  KeepHIVNegativeIndices = sample(1:23, size = 20)
  indexes = c(KeepHIVNegativeIndices, KeepHIVIndices)
  y.downSamp.TRAIN = y[indexes]
  x.downSamp.TRAIN = x[indexes,]
  
  #For test: Keep all 3 remaining HIV- indiciduals. Sample 12 remaining HIV+ individuals
  seqHIV.Neg.Test = which(seq(1,23,1) %!in%  KeepHIVNegativeIndices)
  seqHIV.Pos.Test = which(seq(24,91,1) %!in%  KeepHIVIndices)
  HIV.Neg.Test.finalSeq = sample(seq(24,91,1)[seqHIV.Pos.Test], size = 12)
  indexes.TEST = c(seqHIV.Neg.Test, HIV.Neg.Test.finalSeq)
  #Check to make sure indices in TEST are not in TRAIN
  #which(indexes %in% indexes.TEST)
  #which(HIV.Neg.Test.finalSeq %in% KeepHIVIndices)
  y.TEST = y[indexes.TEST]
  x.TEST = x[indexes.TEST,]
  testSetCount[indexes.TEST] = testSetCount[indexes.TEST] + 1

    #if this is selected then ID 39 corresponds on the temp index. 
    # Train the Random Forest model
    rf_model <- randomForest(x = x.downSamp.TRAIN, y = as.factor(y.downSamp.TRAIN), ntree = 100)
    
    y.downSamp.TRAIN.FACTOR = as.factor(y.downSamp.TRAIN)
    x.downSamp.TRAIN$outcome = y.downSamp.TRAIN
    model_caret <- train(outcome ~ ., data = as.data.frame(x.downSamp.TRAIN), method = "rf", tuneLength = 4)
    explainer <- lime(as.data.frame(x.downSamp.TRAIN[, -ncol(x.downSamp.TRAIN)]), model = model_caret, 
                      bin_continuous = TRUE, predict_function = predict_proba, type = "classification")

    explanation_2 <- explain(x.TEST, explainer, n_features = 63)
  
    
    #Save to Array
    for(j in indexes.TEST){
      savesSEQ = which(explanation_2$case == j)
      index.to.save = testSetCount[j]
      if(index.to.save > num_iterations){
        next #Only save up to num_iterations for each ID
      }
      my_array[j, ,index.to.save] = explanation_2$feature_weight[savesSEQ]

    }
    if(min(testSetCount) > num_iterations){
      break
    }

    #Keep track of a counter of how many times each ID is put into the test group
    #testSetCount[indexes.TEST] = testSetCount[indexes.TEST] + 1
    #Predict for Accuracy
    predictions <- predict(model_caret, newdata = x.TEST)
    predictions.prob <- predict(rf_model, newdata = x.TEST, type = 'prob')
    #Save probility of being assigned HIV+ (over 0.5 means its a 1)
    
    #Need to save all the predictions for each iteraction
    Individual.Predictions[i,indexes.TEST] = as.numeric(predictions.prob[,2])
    
    #accuracy <- sum(predictions == y.TEST) / length(y.TEST)
    #Predict for probabilities
    auc[i]=roc(as.factor(y.TEST), predictions.prob[,1], quiet = TRUE)$auc
    

  print(c("iteration=",i))
}
save(my_array, file = "FeatureWeights_allID_50Iterations.RData")
#my_array[id, feature, iteration]
my_array[1, ,1:20]
#Seq of HIV- who is always guessed correctly
seqHIVNEG_correct = c(1, 2, 4, 5, 9, 13, 14 , 15, 16 , 17, 18, 19, 20, 21, 22, 23)
length(seqHIVNEG_correct)
seqHIVNEG_correct_feature_means = numeric(63)
seqHIVNEG_correct_feature_SDs = numeric(63)
for(i in 1:63){
seqHIVNEG_correct_feature_means[i] = mean(na.omit(my_array[seqHIVNEG_correct, i,1:20]))
seqHIVNEG_correct_feature_SDs[i] = sd(na.omit(my_array[seqHIVNEG_correct, i,1:20]))
}
#Seq of HIV- who tend to be guessed incorrectly
seqHIVNEG_INcorrect = c(6, 7, 11, 12)
#Seq of HIV+ who is always guessed correctly 
seqHIVPOS_correct = c(25, 26, 28, 29, 30, 31, 32, 33, 34, 36, 37, 38, 41, 42, 43, 44, 46, 47, 49, 50, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 63, 65, 66, 68, 70, 71, 73, 74, 78, 79, 80, 81, 82, 83, 84, 85, 86, 90, 91)
length(seqHIVPOS_correct)
#Seq of HIV+ who is usually guessed incorrectly 
seqHIVPOS_INcorrect = c(39, 40, 51, 67, 72, 75, 77, 87)
length(seqHIVPOS_INcorrect)

library(vioplot)
colors = c(rep("#DC3A7A", 8), rep("#9E3934", 10),rep("#ED220D", 10), rep("skyblue", 10),rep("#3675B8", 10), rep("#F8B959", 9), rep("#94C652", 6))
colorsBORDER = c(rep("#141414", 23), rep("#22548C", length(seq(24,66,1))), rep("#B18545", length(seq(67,91,1))))
IDsequence = seqHIVPOS_correct
quartz()
par(mar = c(10, 3, 2, 2)) # Set the margin on all sides to 2
vioplot(c(na.omit(my_array[IDsequence, 1,1:50])),
        c(na.omit(my_array[IDsequence, 2,1:50])),
        c(na.omit(my_array[IDsequence, 3,1:50])),
        c(na.omit(my_array[IDsequence, 4,1:50])),
        c(na.omit(my_array[IDsequence, 5,1:50])),
        c(na.omit(my_array[IDsequence, 6,1:50])),
        c(na.omit(my_array[IDsequence, 7,1:50])),
        c(na.omit(my_array[IDsequence, 8,1:50])),
        c(na.omit(my_array[IDsequence, 9,1:50])),
        c(na.omit(my_array[IDsequence, 10,1:50])),
        c(na.omit(my_array[IDsequence, 11,1:50])),
        c(na.omit(my_array[IDsequence, 12,1:50])),
        c(na.omit(my_array[IDsequence, 13,1:50])),
        c(na.omit(my_array[IDsequence, 14,1:50])),
        c(na.omit(my_array[IDsequence, 15,1:50])),
        c(na.omit(my_array[IDsequence, 16,1:50])),
        c(na.omit(my_array[IDsequence, 17,1:50])),
        c(na.omit(my_array[IDsequence, 18,1:50])),
        c(na.omit(my_array[IDsequence, 19,1:50])),
        c(na.omit(my_array[IDsequence, 20,1:50])),
        c(na.omit(my_array[IDsequence, 21,1:50])),
        c(na.omit(my_array[IDsequence, 22,1:50])),
        c(na.omit(my_array[IDsequence, 23,1:50])),
        c(na.omit(my_array[IDsequence, 24,1:50])),
        c(na.omit(my_array[IDsequence, 25,1:50])),
        c(na.omit(my_array[IDsequence, 26,1:50])),
        c(na.omit(my_array[IDsequence, 27,1:50])),
        c(na.omit(my_array[IDsequence, 28,1:50])),
        c(na.omit(my_array[IDsequence, 29,1:50])),
        c(na.omit(my_array[IDsequence, 30,1:50])),
        c(na.omit(my_array[IDsequence, 31,1:50])),
        c(na.omit(my_array[IDsequence, 32,1:50])),
        c(na.omit(my_array[IDsequence, 33,1:50])),
        c(na.omit(my_array[IDsequence, 34,1:50])),
        c(na.omit(my_array[IDsequence, 35,1:50])),
        c(na.omit(my_array[IDsequence, 36,1:50])),
        c(na.omit(my_array[IDsequence, 37,1:50])),
        c(na.omit(my_array[IDsequence, 38,1:50])),
        c(na.omit(my_array[IDsequence, 39,1:50])),
        c(na.omit(my_array[IDsequence, 40,1:50])),
        c(na.omit(my_array[IDsequence, 41,1:50])),
        c(na.omit(my_array[IDsequence, 42,1:50])),
        c(na.omit(my_array[IDsequence, 43,1:50])),
        c(na.omit(my_array[IDsequence, 44,1:50])),
        c(na.omit(my_array[IDsequence, 45,1:50])),
        c(na.omit(my_array[IDsequence, 46,1:50])),
        c(na.omit(my_array[IDsequence, 47,1:50])),
        c(na.omit(my_array[IDsequence, 48,1:50])),
        c(na.omit(my_array[IDsequence, 49,1:50])),
        c(na.omit(my_array[IDsequence, 50,1:50])),
        c(na.omit(my_array[IDsequence, 51,1:50])),
        c(na.omit(my_array[IDsequence, 52,1:50])),
        c(na.omit(my_array[IDsequence, 53,1:50])),
        c(na.omit(my_array[IDsequence, 54,1:50])),
        c(na.omit(my_array[IDsequence, 55,1:50])),
        c(na.omit(my_array[IDsequence, 56,1:50])),
        c(na.omit(my_array[IDsequence, 57,1:50])),
        c(na.omit(my_array[IDsequence, 58,1:50])),
        c(na.omit(my_array[IDsequence, 59,1:50])),
        c(na.omit(my_array[IDsequence, 60,1:50])),
        c(na.omit(my_array[IDsequence, 61,1:50])),
        c(na.omit(my_array[IDsequence, 62,1:50])),
        c(na.omit(my_array[IDsequence, 63,1:50])),
        main = "HIV+ ID 39 ",
        ylim = c(-0.35,0.35), 
        lwd = 2, 
        cex.axis = 1.0, 
        col = colors,
        names = feature.labels, 
        las = 2
)
#Idea: Ablation based on parameters whose feature weight is ~0 but who have a large dispersion around that median feature weight of ~0
