
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

###############################################################
#### Gaussian Mixture Model : good copy 
###############################################################
#GMM
library(MASS)
library(mclust)

# Fit a Gaussian Mixture Model to the data and generate synthetic data

gmm_synthetic_data <- function(data, num_samples) {
  # Fit the GMM model to the data
  gmm_model <- Mclust(data)
  
  # Extract model parameters
  means <- gmm_model$parameters$mean
  variances <- gmm_model$parameters$variance$sigma
  probs <- gmm_model$parameters$pro
  
  # Get column names to preserve them
  col_names <- colnames(data)
  
  # Function to regenerate a sample if it contains negative values
  regenerate_if_negative <- function(label) {
    sample <- mvrnorm(1, means[, label], variances[,, label])
    while(any(sample < 0)) {  # Check for negative values
      sample <- mvrnorm(1, means[, label], variances[,, label])  # Regenerate if negative value found
    }
    return(sample)
  }
  
  # Sample component labels based on the mixture probabilities
  component_labels <- sample(1:length(probs), size = num_samples, replace = TRUE, prob = probs)
  
  # Generate synthetic data by sampling from the corresponding Gaussian distributions with negative value handling
  synthetic_data <- t(sapply(component_labels, function(label) {
    regenerate_if_negative(label)
  }))
  
  # Convert to a data frame and set the column names to preserve the original feature names
  synthetic_data <- as.data.frame(synthetic_data)
  colnames(synthetic_data) <- col_names
  
  return(synthetic_data)
}

gmm_synthetic_data_2 <- function(data, num_samples, max_iter = 100000) {
  # Fit the GMM model to the data
  gmm_model <- Mclust(data)
  
  # Check if GMM fitting was successful
  if (is.null(gmm_model$parameters)) {
    stop("GMM model fitting failed. Please check your input data.")
  }
  
  # Extract model parameters
  means <- gmm_model$parameters$mean
  variances <- gmm_model$parameters$variance$sigma
  probs <- gmm_model$parameters$pro
  
  # Get column names to preserve them
  col_names <- colnames(data)
  
  # Function to regenerate a sample if it contains negative values
  regenerate_if_negative <- function(label, max_iter) {
    iter <- 0
    while (TRUE) {
      # Try generating the sample with the current iteration count
      sample <- mvrnorm(1, means[, label], variances[,, label])
      iter <- iter + 1
      
      # Check if the sample is valid (i.e., all non-negative values)
      if (all(sample >= 0)) {
        return(sample)
      }
      
      # If max_iter is reached, reset seed and iteration counter
      if (iter >= max_iter) {
        cat("Max iterations reached. Choosing a new seed...\n")
        new_seed <- sample(1:1e6, 1)
        set.seed(new_seed)
        cat("New seed set to:", new_seed, "\n")
        iter <- 0  # Reset iteration counter
      }
    }
  }
  
  # Sample component labels based on the mixture probabilities
  component_labels <- sample(1:length(probs), size = num_samples, replace = TRUE, prob = probs)
  
  # Generate synthetic data by sampling from the corresponding Gaussian distributions with negative value handling
  synthetic_data <- t(sapply(component_labels, function(label) {
    regenerate_if_negative(label, max_iter)
  }))
  
  # Convert to a data frame and set the column names to preserve the original feature names
  synthetic_data <- as.data.frame(synthetic_data)
  colnames(synthetic_data) <- col_names
  
  return(synthetic_data)
}


# Example usage
set.seed(49)
#input_data <- matrix(runif(90 * 60), nrow = 90, ncol = 60)  # Example input matrix
input_data = x 
num_samples = length(x[,1])
synthetic_data_gmm <- gmm_synthetic_data(input_data, num_samples)

##Generate synthetic data for HIV- 

set.seed(3)
input_data_hivNEG = x[1:23,]
num_samples_NEG = 23
synthetic_data_gmm_hivNEG <- gmm_synthetic_data_2(input_data_hivNEG, num_samples_NEG)


##Generate synthetic data for HIV+
input_data_hivPOS = x[24:91,]
num_samples_POS = 68
synthetic_data_gmm_hivPOS <- gmm_synthetic_data(input_data_hivPOS, num_samples_POS)

###Next need to combine the two synthetic arrays back to a single feature matrix
GMM_synthetic = data.frame(matrix(ncol = length(x[1,]), nrow = 91))
GMM_synthetic[1:23,] = synthetic_data_gmm_hivNEG
GMM_synthetic[24:91,] = synthetic_data_gmm_hivPOS
col_names <- colnames(x)
colnames(GMM_synthetic) <- col_names

write.csv(GMM_synthetic, "GMM_synthetic_23HIVneg_68HIVpos.csv")

#Generate lots of synethic data via the GMM approach. 

# Main loop to generate synthetic data in batches
generate_synthetic_samples <- function(data, num_samples_per_batch, n_iterations) {
  all_synthetic_data <- data.frame()
  
  for (i in 1:n_iterations) {

    set.seed(i)  # Set seed to match the loop index
    
    
    # Generate synthetic samples
    synthetic_batch <- gmm_synthetic_data_2(data, num_samples_per_batch)
    
    # Add the batch to the combined data frame
    all_synthetic_data <- rbind(all_synthetic_data, synthetic_batch)
  }
  
  return(all_synthetic_data)
}

input_data_hivNEG = x[1:23,]
num_samples_NEG = 23
combined_synthetic_data <- generate_synthetic_samples(input_data_hivNEG, num_samples_per_batch = num_samples_NEG, n_iterations = 10)
write.csv(combined_synthetic_data, "GMM_synthetic_23HIVneg_230Individuals.csv")

input_data_hivPOS = x[24:91,]
num_samples_POS = 68
combined_synthetic_data_POS <- generate_synthetic_samples(input_data_hivPOS, num_samples_per_batch = num_samples_POS, n_iterations = 10)
write.csv(combined_synthetic_data_POS, "GMM_synthetic_23HIVPOS_680Individuals.csv")


GMM_synthetic_all = data.frame(matrix(ncol = length(x[1,]), nrow = 680+230))
GMM_synthetic_all[1:230,] = combined_synthetic_data
GMM_synthetic_all[231:(230+680),] = combined_synthetic_data_POS
col_names <- colnames(x)
colnames(GMM_synthetic_all) <- col_names

write.csv(GMM_synthetic_all, "GMM_synthetic_230HIVneg_680HIVpos.csv")
save(GMM_synthetic_all, file = "GMM_synthetic_10x.RData")
#comebine and then save. 


###############################################################
############ Multivariate normal
###############################################################

library(MASS)
library(dplyr)

generate_positive_data_MVN <- function(n, mu, Sigma) {
  synthetic_data <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  # Check and replace negative values
  for (i in 1:ncol(synthetic_data)) {
    while (any(synthetic_data[,i] < 0)) {
      # Resample only the negative entries
      negatives <- synthetic_data[,i] < 0
      synthetic_data[negatives,i] <- MASS::mvrnorm(n = sum(negatives), mu = mu[i], Sigma = Sigma[i, i])
    }
  }
  return(synthetic_data)
}

control_data <- x[1:23,]
HIVposData <- x[24:91,]

# Calculate mean and covariance matrix
control_mean <- colMeans(control_data)
control_cov <- cov(control_data)

HIV_mean <- colMeans(HIVposData)
HIV_cov <- cov(HIVposData)

# Number of synthetic samples to generate
num_samples_control <- nrow(control_data)
num_samples_HIV <- nrow(HIVposData)


synthetic_data_control <- generate_positive_data_MVN(num_samples_control, control_mean, control_cov)
synthetic_data_control_df <- as.data.frame(synthetic_data_control)

synthetic_data_HIV <- generate_positive_data_MVN(num_samples_HIV, HIV_mean, HIV_cov)
synthetic_data_HIV_df <- as.data.frame(synthetic_data_HIV)

#Stitch the two feature matrices together
#First 23 are Control the remaining are HIV+
features_all.imp.MVsynthetic = data.frame(matrix(ncol = ncol(x), nrow = nrow(x)))
features_all.imp.MVsynthetic[1:num_samples_control,] = synthetic_data_control_df
features_all.imp.MVsynthetic[24:91,] = synthetic_data_HIV_df
MVN_synthetic = features_all.imp.MVsynthetic
colnames(MVN_synthetic) <- col_names

write.csv(MVN_synthetic, "MVN_synthetic_23HIVneg_68HIVpos.csv")


###############################################################
############ KNN  : good copy 
###############################################################

library(FNN)
knn_interpolation <- function(data, num_samples, k) {
  # Find the k-nearest neighbors for each data point
  nbrs <- get.knn(data, k = k)
  
  synthetic_data <- matrix(NA, nrow = num_samples, ncol = ncol(data))
  for (i in 1:num_samples) {
    # Randomly choose a point and one of its neighbors
    idx <- sample(1:nrow(data), 1)
    neighbor_idx <- sample(nbrs$nn.index[idx, ], 1)
    
    # Generate a new sample by interpolating between the chosen point and its neighbor
    alpha <- runif(1)  # Random interpolation factor between 0 and 1
    synthetic_data[i, ] <- data[idx, ] + alpha * (data[neighbor_idx, ] - data[idx, ])
  }
  
  return(synthetic_data)
}

# Example usage
set.seed(123)
input_data <- matrix(runif(90 * 60), nrow = 90, ncol = 60)  # Example input matrix
num_samples <- 100  # Number of synthetic samples to generate

synthetic_data_knn <- knn_interpolation(input_data, num_samples)


##Generate synthetic data for HIV- 
input_data_hivNEG = x[1:23,]
input_data_hivNEG = as.matrix(input_data_hivNEG)
num_samples = 23
synthetic_data_knn_hivNEG_k3 <- knn_interpolation(input_data_hivNEG, num_samples, k = 3)
synthetic_data_knn_hivNEG_koptimal <- knn_interpolation(input_data_hivNEG, num_samples, k = round(sqrt(91)))
synthetic_data_knn_hivNEG_k20 <- knn_interpolation(input_data_hivNEG, num_samples, k = 20)


##Generate synthetic data for HIV+
input_data_hivPOS = x[24:91,]
input_data_hivPOS = as.matrix(input_data_hivPOS)
num_samples = 68
#synthetic_data_knn_hivPOS <- knn_interpolation(input_data_hivPOS, num_samples)

synthetic_data_knn_hivPOS_k3 <- knn_interpolation(input_data_hivPOS, num_samples, k = 3)
synthetic_data_knn_hivPOS_koptimal <- knn_interpolation(input_data_hivPOS, num_samples, k = round(sqrt(91)))
synthetic_data_knn_hivPOS_k20 <- knn_interpolation(input_data_hivPOS, num_samples, k = 20)


###Next need to combine the two synthetic arrays back to a single feature matrix
KNN_synthetic_k3 = data.frame(matrix(ncol = length(x[1,]), nrow = 91))
KNN_synthetic_k3[1:23,] = synthetic_data_knn_hivNEG_k3
KNN_synthetic_k3[24:91,] = synthetic_data_knn_hivPOS_k3

KNN_synthetic_koptimal = data.frame(matrix(ncol = length(x[1,]), nrow = 91))
KNN_synthetic_koptimal[1:23,] = synthetic_data_knn_hivNEG_koptimal
KNN_synthetic_koptimal[24:91,] = synthetic_data_knn_hivPOS_koptimal

KNN_synthetic_k20 = data.frame(matrix(ncol = length(x[1,]), nrow = 91))
KNN_synthetic_k20[1:23,] = synthetic_data_knn_hivNEG_k20
KNN_synthetic_k20[24:91,] = synthetic_data_knn_hivPOS_k20


col_names <- colnames(x)
colnames(KNN_synthetic_k3) <- col_names
colnames(KNN_synthetic_koptimal) <- col_names
colnames(KNN_synthetic_k20) <- col_names


write.csv(KNN_synthetic_k3, "KNN_synthetic_k3_23HIVneg_68HIVpos.csv")
write.csv(KNN_synthetic_koptimal, "KNN_synthetic_koptimal_23HIVneg_68HIVpos.csv")
write.csv(KNN_synthetic_k20, "KNN_synthetic_k20_23HIVneg_68HIVpos.csv")


###############################################################
############ SMOTE  : good copy 
###############################################################

library(smotefamily)

generate_synthetic_data_equal_size <- function(data, labels, K = 5) {
  # Convert data to a data frame
  data <- as.data.frame(data)
  
  # Separate the minority and majority classes
  minority_data <- data[labels == 0, ]
  majority_data <- data[labels == 1, ]
  
  minority_labels <- labels[labels == 0]
  majority_labels <- labels[labels == 1]
  
  # Calculate the number of synthetic samples needed for each class to match the original size
  num_synthetic_minority <- length(minority_labels)
  num_synthetic_majority <- length(majority_labels)
  
  # Apply SMOTE to generate synthetic samples for the minority class
  smote_minority_result <- SMOTE(X = minority_data, target = factor(minority_labels), K = K, dup_size = ceiling(num_synthetic_minority / nrow(minority_data)))
  
  # Apply SMOTE to generate synthetic samples for the majority class
  smote_majority_result <- SMOTE(X = majority_data, target = factor(majority_labels), K = K, dup_size = ceiling(num_synthetic_majority / nrow(majority_data)))
  
  # Extract the synthetic data
  synthetic_minority_data <- smote_minority_result$syn_data[, -ncol(smote_minority_result$syn_data)]  # Remove the label column
  synthetic_majority_data <- smote_majority_result$syn_data[, -ncol(smote_majority_result$syn_data)]  # Remove the label column
  
  synthetic_minority_labels <- rep(0, nrow(synthetic_minority_data))
  synthetic_majority_labels <- rep(1, nrow(synthetic_majority_data))
  
  # Combine the synthetic data
  synthetic_data <- rbind(synthetic_minority_data, synthetic_majority_data)
  synthetic_labels <- c(synthetic_minority_labels, synthetic_majority_labels)
  
  # Randomly shuffle the synthetic data to mix classes
  #set.seed(123)  # For reproducibility
  #shuffle_indices <- sample(1:nrow(synthetic_data))
  
  #synthetic_data <- synthetic_data[shuffle_indices, ]
  #synthetic_labels <- synthetic_labels[shuffle_indices]
  
  return(list(data = synthetic_data, labels = synthetic_labels))
}

# Example usage
set.seed(123)
input_data <- matrix(runif(91 * 63), nrow = 91, ncol = 63)  # Example input matrix with 91 samples and 63 features
input_labels <- c(rep(0, 23), rep(1, 68))  # 23 class 0s and 68 class 1s

# Generate synthetic data of the same size and class distribution as the original data
result <- generate_synthetic_data_equal_size(input_data, input_labels)

# Extract the synthetic data and labels
synthetic_data <- result$data
synthetic_labels <- result$labels

# Check the class distribution
table(synthetic_labels)

input_data = x 
input_labels = y
synthetic_data_SMOTE <- generate_synthetic_data_equal_size(input_data, input_labels)

synthetic_data_SMOTE$data
synthetic_data_SMOTE$labels

write.csv(synthetic_data_SMOTE$data, "SMOTE_synthetic_K5_23HIVneg_68HIVpos.csv")

##Generate synthetic data for HIV- 
#input_data_hivNEG = x[1:23,]
#input_labels = y[1:23]
#num_samples = 23
#synthetic_data_SMOTE_hivNEG <- gmm_synthetic_data(input_data_hivNEG, input_labels)


##Generate synthetic data for HIV+
#input_data_hivPOS = x[24:91,]
#input_labels = y[]
#num_samples = 68
#synthetic_data_SMOTE_hivPOS <- gmm_synthetic_data(input_data_hivPOS, input_labels)

###############################################################
#### ROSE : over sampling technique
###############################################################

library(ROSE)
x = features_all.imp[1:91,2:64]
y = health_outcome

# Simulate the original dataset with 91 individuals and 63 features
set.seed(2024)
data <- data.frame(matrix(rnorm(91 * 63), nrow = 91, ncol = 63))
data$class <- factor(c(rep(0, 23), rep(1, 68)))  # Class 0 has 23 individuals, Class 1 has 68 individuals

data = x
data$class = y
# Generate synthetic data using ROSE with oversampling
synthetic_data <- ROSE(class ~ ., data = data, N = 2 * nrow(data))$data

# Function to resample only negative values in the synthetic data
resample_negative_values <- function(sampled_row, class_data) {
  for (j in 1:(length(sampled_row) - 1)) {  # Loop over all columns except the class column
    while (sampled_row[j] < 0) {
      # Resample the negative value only from the corresponding column in the class_data
      sampled_row[j] <- class_data[sample(1:nrow(class_data), 1), j]
    }
  }
  return(sampled_row)
}

# Separate the synthetic data by class
synthetic_class_0 <- synthetic_data[synthetic_data$class == 0, ]
synthetic_class_1 <- synthetic_data[synthetic_data$class == 1, ]

# Apply the resampling of negative values to each row, ensuring to pass the entire data matrix
final_synthetic_class_0 <- t(apply(synthetic_class_0, 1, resample_negative_values, class_data = synthetic_class_0))
final_synthetic_class_1 <- t(apply(synthetic_class_1, 1, resample_negative_values, class_data = synthetic_class_1))

# Combine the sampled synthetic data
final_synthetic_data <- rbind(final_synthetic_class_0, final_synthetic_class_1)

# Convert back to data.frame and set correct column names
final_synthetic_data <- as.data.frame(final_synthetic_data)
colnames(final_synthetic_data) <- colnames(data)

# Ensure the class column is a factor
final_synthetic_data$class <- factor(final_synthetic_data$class, levels = c(0, 1))

# Randomly sample 23 class 0 and 68 class 1 individuals, and organize them in the desired order
random_number <- sample(1:10000, 1)
set.seed(random_number)  # Set seed for reproducibility
sampled_class_0 <- final_synthetic_data[final_synthetic_data$class == 0, ]
sampled_class_1 <- final_synthetic_data[final_synthetic_data$class == 1, ]

sampled_class_0 <- sampled_class_0[sample(1:nrow(sampled_class_0), 23), ]
sampled_class_1 <- sampled_class_1[sample(1:nrow(sampled_class_1), 68), ]

# Combine them into a new data frame with the desired order
final_sampled_data <- rbind(sampled_class_0, sampled_class_1)

# Verify the structure and distribution of the new data frame
print(dim(final_sampled_data))  # Should be 91 rows and 64 columns (including the class column)
print(table(final_sampled_data$class))  # 

# The final_synthetic_data object now contains the generated data with preserved class balance
# and all non-negative values.
ROSE_synthetic = final_sampled_data[,-64]
x.ROSE = ROSE_synthetic


pca <- prcomp(ROSE_synthetic, scale. = TRUE)
pc <-predict(pca,ROSE_synthetic)


quartz()
plot(pc[,1], pc[,2], pch = 20, cex = 1.5, yaxt = "n", xaxt = "n", )
axis(1,las= 1, tck=0.02,mgp=c(3, .8, 0), cex.axis = 1.8)
axis(2,las= 1, tck=0.02,mgp=c(3, .5, 0), cex.axis = 1.5)
points(pc[24:91,1], pc[24:91,2], col = "#C42D20", pch = 15, cex = 1.5)
points(pc[1:23,1], pc[1:23,2], col = "#3474B2", pch = 20, cex = 2.0)

dataEllipse(pc[24:91,1], pc[24:91,2], col = "#C42D20", lwd = 2, add = TRUE, levels = c(0.68), center.pch = NULL)
dataEllipse(pc[1:23,1], pc[1:23,2], col = "#3474B2", lwd = 2, add = TRUE, levels = c(0.68), center.pch = NULL)
###############################################################
############ Save all
###############################################################

save(x, y, GMM_synthetic, MVN_synthetic, KNN_synthetic_k3, KNN_synthetic_koptimal,KNN_synthetic_k20,synthetic_data_SMOTE, x.ROSE, file = "allSyntheticData.RData")



###############################################################
#### OLD CODE: just for example
###############################################################
#GMM

library(mclust)

# Fit a Gaussian Mixture Model to the data and generate synthetic data
#gmm_synthetic_data <- function(data, num_samples) {
# Fit the GMM model to the data
# gmm_model <- Mclust(data)

# Generate synthetic data from the GMM model
#  synthetic_data <- rGMM(num_samples, gmm_model$parameters)

# return(synthetic_data)
#}

gmm_synthetic_data <- function(data, num_samples) {
  # Fit the GMM model to the data
  gmm_model <- Mclust(data)
  
  # Extract model parameters
  means <- gmm_model$parameters$mean
  variances <- gmm_model$parameters$variance$sigma
  probs <- gmm_model$parameters$pro
  
  # Sample component labels based on the mixture probabilities
  component_labels <- sample(1:length(probs), size = num_samples, replace = TRUE, prob = probs)
  
  # Generate synthetic data by sampling from the corresponding Gaussian distributions
  synthetic_data <- sapply(component_labels, function(label) {
    mvrnorm(1, means[, label], variances[,, label])
  })
  
  return(t(synthetic_data))
}

# Example usage
set.seed(123)
input_data <- matrix(runif(90 * 60), nrow = 90, ncol = 60)  # Example input matrix
num_samples <- 100  # Number of synthetic samples to generate

synthetic_data_gmm <- gmm_synthetic_data(input_data, num_samples)


################################################################
############ KNN 
################################################################
#install.packages("FNN")
library(FNN)
knn_interpolation <- function(data, num_samples, k = 5) {
  # Find the k-nearest neighbors for each data point
  nbrs <- get.knn(data, k = k)
  
  synthetic_data <- matrix(NA, nrow = num_samples, ncol = ncol(data))
  for (i in 1:num_samples) {
    # Randomly choose a point and one of its neighbors
    idx <- sample(1:nrow(data), 1)
    neighbor_idx <- sample(nbrs$nn.index[idx, ], 1)
    
    # Generate a new sample by interpolating between the chosen point and its neighbor
    alpha <- runif(1)  # Random interpolation factor between 0 and 1
    synthetic_data[i, ] <- data[idx, ] + alpha * (data[neighbor_idx, ] - data[idx, ])
  }
  
  return(synthetic_data)
}

# Example usage
set.seed(123)
input_data <- matrix(runif(90 * 60), nrow = 90, ncol = 60)  # Example input matrix
num_samples <- 100  # Number of synthetic samples to generate

synthetic_data_knn <- knn_interpolation(input_data, num_samples)


################################################################
############ SMOTE
################################################################
install.packages("smotefamily")
library(DMwR2)
library(smotefamily)



smote_synthetic_data <- function(data, num_samples) {
  # Ensure the data is a data frame
  data <- as.data.frame(data)
  
  # Assign random class labels (0 or 1) to the data
  set.seed(123)  # Ensure reproducibility
  random_labels <- factor(sample(0:1, nrow(data), replace = TRUE))
  
  # Apply SMOTE to generate synthetic samples
  smote_result <- SMOTE(X = data, target = random_labels, K = 5, dup_size = ceiling(num_samples / nrow(data)))
  
  # Extract the synthetic data (without the labels)
  synthetic_data <- smote_result$syn_data[, -ncol(smote_result$syn_data)]  # Remove the label column
  
  return(as.matrix(synthetic_data))
}

# Example usage
set.seed(123)
input_data <- matrix(runif(90 * 60), nrow = 90, ncol = 60)  # Example input matrix
num_samples <- 100  # Number of synthetic samples to generate

synthetic_data_smote <- smote_synthetic_data(input_data, num_samples)



library(smotefamily)

generate_synthetic_data_equal_size <- function(data, labels, K = 5) {
  # Convert data to a data frame
  data <- as.data.frame(data)
  
  # Separate the minority and majority classes
  minority_data <- data[labels == 0, ]
  majority_data <- data[labels == 1, ]
  
  minority_labels <- labels[labels == 0]
  majority_labels <- labels[labels == 1]
  
  # Calculate the number of synthetic samples needed for each class to match the original size
  num_synthetic_minority <- length(minority_labels)
  num_synthetic_majority <- length(majority_labels)
  
  # Apply SMOTE to generate synthetic samples for the minority class
  smote_minority_result <- SMOTE(X = minority_data, target = factor(minority_labels), K = K, dup_size = ceiling(num_synthetic_minority / nrow(minority_data)))
  
  # Apply SMOTE to generate synthetic samples for the majority class
  smote_majority_result <- SMOTE(X = majority_data, target = factor(majority_labels), K = K, dup_size = ceiling(num_synthetic_majority / nrow(majority_data)))
  
  # Extract the synthetic data
  synthetic_minority_data <- smote_minority_result$syn_data[, -ncol(smote_minority_result$syn_data)]  # Remove the label column
  synthetic_majority_data <- smote_majority_result$syn_data[, -ncol(smote_majority_result$syn_data)]  # Remove the label column
  
  synthetic_minority_labels <- rep(0, nrow(synthetic_minority_data))
  synthetic_majority_labels <- rep(1, nrow(synthetic_majority_data))
  
  # Combine the synthetic data
  synthetic_data <- rbind(synthetic_minority_data, synthetic_majority_data)
  synthetic_labels <- c(synthetic_minority_labels, synthetic_majority_labels)
  
  # Randomly shuffle the synthetic data to mix classes
  set.seed(123)  # For reproducibility
  shuffle_indices <- sample(1:nrow(synthetic_data))
  
  synthetic_data <- synthetic_data[shuffle_indices, ]
  synthetic_labels <- synthetic_labels[shuffle_indices]
  
  return(list(data = synthetic_data, labels = synthetic_labels))
}

# Example usage
set.seed(123)
input_data <- matrix(runif(91 * 63), nrow = 91, ncol = 63)  # Example input matrix with 91 samples and 63 features
input_labels <- c(rep(0, 23), rep(1, 68))  # 23 class 0s and 68 class 1s

# Generate synthetic data of the same size and class distribution as the original data
result <- generate_synthetic_data_equal_size(input_data, input_labels)

# Extract the synthetic data and labels
synthetic_data <- result$data
synthetic_labels <- result$labels

# Check the class distribution
table(synthetic_labels)

# 
# ROSE (Random Over-Sampling Examples) and SMOTE (Synthetic Minority Over-sampling Technique) are both techniques used to address class imbalance in datasets by generating synthetic data. However, they approach this problem in different ways, leading to differences in the nature of the synthetic data they produce and their impact on downstream tasks like classification.
# 
# 1. Core Methodology:
#   ROSE:
#   
#   Approach: ROSE generates synthetic examples by randomly perturbing existing instances in the feature space. It does so by adding a random noise component to the data, which is sampled from a kernel density estimate around the existing data points.
# Resulting Data: The synthetic data points generated by ROSE are close to the original data points but can exhibit more variability due to the added noise. This can lead to the creation of synthetic examples that are not directly interpolations between two existing points but rather lie in the vicinity of a single data point.
# SMOTE:
#   
#   Approach: SMOTE generates synthetic examples by interpolating between existing instances of the minority class. It selects a random minority class instance and then identifies its k-nearest neighbors in the feature space. A synthetic point is generated along the line segment between the selected instance and one of its neighbors.
# Resulting Data: The synthetic data points generated by SMOTE are direct interpolations between existing minority class instances. This tends to create new instances that are linearly distributed between the original data points, preserving the local structure of the minority class.
# 2. Impact on Local Structure:
#   ROSE:
#   
#   Preservation of Local Structure: ROSE may not preserve the local structure as well as SMOTE because it introduces random noise. This noise can cause synthetic data points to drift slightly away from the original clusters, potentially leading to less precise representation of local patterns.
# Variability: ROSE can introduce more variability into the synthetic data, which might be beneficial in certain scenarios (e.g., avoiding overfitting) but could also dilute the distinction between classes if the noise is too large.
# SMOTE:
#   
#   Preservation of Local Structure: SMOTE is designed to preserve local structures by creating synthetic data points that lie along the line segments between existing data points. This tends to maintain the original class distribution and the relationships between points within the minority class.
# Risk of Overfitting: By generating synthetic points only along the lines between existing minority instances, SMOTE might lead to overfitting, particularly in scenarios where the minority class instances are very close together or form tight clusters.
# 3. Applicability and Use Cases:
#   ROSE:
#   
#   General Use: ROSE is more general-purpose and can be applied to both binary and multiclass problems. It is particularly useful when there's a need to generate synthetic data that introduces some variability, which can help in generalizing the model.
# Balancing Flexibility: ROSE allows for flexible oversampling by adjusting the amount of noise and controlling how synthetic data points are generated. This makes it adaptable to a wide range of datasets and imbalance scenarios.
# SMOTE:
# 
# Targeted Use: SMOTE is specifically designed for binary or multiclass problems where the focus is on increasing the number of minority class instances. It's widely used when the goal is to enhance the representation of minority classes without altering the majority class.
# Cluster Preservation: SMOTE is particularly effective when it's important to maintain the integrity of minority class clusters in the feature space, as it doesn't introduce noise but rather interpolates directly between points.
# 4. Real-World Consequence of Differences:
#   ROSE:
#   
#   Potential Issue: If the added noise in ROSE is not carefully controlled, it could lead to the creation of synthetic data points that do not represent the true nature of the minority class, potentially leading to less accurate models.
# Benefit: The added variability can help prevent models from overfitting to specific patterns in the minority class, making them more robust in certain scenarios.
# SMOTE:
#   
#   Potential Issue: SMOTE might lead to overfitting, especially when the synthetic data is too closely tied to the existing minority instances, potentially causing the model to perform poorly on new, unseen data.
# Benefit: By preserving local structures and relationships between minority class instances, SMOTE can produce more realistic synthetic data that accurately reflects the underlying distribution of the minority class.
# Summary:
#   ROSE introduces random perturbations to generate synthetic data, which can increase variability but may not preserve local structures as effectively as SMOTE.
# SMOTE generates synthetic data by interpolating between existing minority class instances, which preserves local structure but might lead to overfitting if the synthetic data is too similar to the original data.
# The choice between ROSE and SMOTE depends on the specific requirements of the task, such as the importance of local structure preservation versus the need to introduce variability to prevent overfitting.