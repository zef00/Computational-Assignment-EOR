rm(list=ls())
library(ggplot2)
library(MASS)
library(tidyverse)
library(zoo)

#Datasets
ScanRecords <- read.csv("ScanRecords.csv")
ScanRecords$Date <- as.Date(ScanRecords$Date)
Type1 <- dplyr::filter(ScanRecords, PatientType=="Type 1") 
Type2 <- dplyr::filter(ScanRecords, PatientType=="Type 2") 


#Duration Scatterplot 
ggplot(ScanRecords, aes(x=1:618, y=Duration, color = PatientType))+
  geom_point() +
  geom_hline(yintercept = mean(Type1$Duration), color = "red") +
  geom_hline(yintercept = mean(Type2$Duration), color = "blue")

#Duration Histogram
ggplot(Type1) + 
  geom_histogram(aes(x = Duration))
ggplot(Type2) + 
  geom_histogram(aes(x = Duration))

#Duration qqplots
qqnorm(Type1$Duration) #This suggests a normal distribution, which corresponds to the information given
qqnorm(Type2$Duration) #This one is bit weird 

#Duration Normal Distribution PDF plot
pdf1 <- dnorm(Type1$Duration, mean=mean(Type1$Duration), sd = sd(Type1$Duration))
plot(Type1$Duration, pdf1)
pdf2 <- dnorm(Type2$Duration, mean=mean(Type2$Duration), sd = sd(Type2$Duration))
plot(Type2$Duration, pdf2)

#Bootstrapping the mean duration for type 2 
X <- Type2$Duration
B <- 1000 

n <- length(X)                                          # Sample size
X.bar <- mean(X)                                        # Sample mean of X
St.Dev <- sd(X)                                         # Standard deviation of X

X.star_mean <- rep(NA, n)
X.star_sd <- rep(NA,n)

for (b in 1:B) {
  J <- sample.int(n, size = n, replace = TRUE)            # Draw the indices J
  X.star <- X[J]                                      # Draw the bootstrap sample
  X.star_mean[b] <- mean(X.star) 
  X.star_sd[b] <- sd(X.star)
}

mean_duration_type2 <- mean(X.star_mean) #these are the variables of interest to us
SD_duration_type2 <- mean(X.star_sd)

# Monte Carlo Study for the type 2 bootstrap mean (this can take long)
N = 2000
B = 500 
bias = rep(NA, N)
for (i in 1:N){
  mc_type2_duration <- rnorm(length(Type2$Duration),
                             mean = mean_duration_type2,
                             sd = SD_duration_type2)
  X.star_mean <- rep(NA, length(Type2$Duration))
  for (b in 1:B){
    J <- sample.int(length(Type2$Duration), size = length(Type2$Duration), replace = TRUE)
    X.star <- mc_type2_duration[J]                                      # Draw the bootstrap sample
    X.star_mean[b] <- mean(X.star) 
  }
  MCbootstrap_mean <- mean(X.star_mean)
  bias[i] <- MCbootstrap_mean - mean(mc_type2_duration)
}
avg_bias = mean(bias) #very low so should be good 


## BOOTSTRAP patient 1 arrival times with interarival times

# Calculate inter-arrival times
Type1 <- Type1 %>%
  group_by(Date) %>%
  mutate(
    InterArrivalTime = Time - lag(Time),  # Standard inter-arrival time within the same day
    FirstOfDay = row_number() == 1        # Flag for the first record of the day
  ) %>%
  ungroup()

# Handle first times of each day
Type1 <- Type1 %>%
  mutate(
    InterArrivalTime = ifelse(
      FirstOfDay,
      # If first record of the day, calculate time difference with lagged time + 9 hours (17,00 - 8,00 = 9,00)
      Time - lag(Time) + 9,
      InterArrivalTime
    )
  )

# Remove remaining NA values (the first time of the first day)
InterArrivalTimes <- Type1$InterArrivalTime[!is.na(Type1$InterArrivalTime)]

# Calculate rate parameter for exponential distribution
lambda <- 1 / mean(InterArrivalTimes)

# Perform bootstrapping based on exponential distribution
B <- 10000  # Number of bootstrap samples
bootstrap_means <- numeric(B)

for (b in seq_len(B)) {
  # Generate bootstrap sample from exponential distribution
  bootstrap_sample <- rexp(length(InterArrivalTimes), rate = lambda)
  bootstrap_means[b] <- mean(bootstrap_sample)
}

# Calculate bootstrapped mean and 95% confidence interval
bootstrapped_mean <- mean(bootstrap_means)
ci <- quantile(bootstrap_means, probs = c(0.025, 0.975))

# Print results
cat("Bootstrapped Mean Inter-Arrival Time (Exponential):", bootstrapped_mean, "\n")
cat("95% Confidence Interval for Mean Inter-Arrival Time (Exponential):", ci, "\n")

# Plot the bootstrap distribution
ggplot(data.frame(bootstrap_means = bootstrap_means), aes(x = bootstrap_means)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_vline(xintercept = bootstrapped_mean, color = "red", linetype = "dashed") +
  geom_vline(xintercept = ci[1], color = "blue", linetype = "dotted") +
  geom_vline(xintercept = ci[2], color = "blue", linetype = "dotted") +
  labs(title = "Bootstrap Distribution of Mean Inter-Arrival Time (Exponential)",
       subtitle = paste("95% CI:", round(ci[1], 3), "-", round(ci[2], 3)),
       x = "Mean Inter-Arrival Time", y = "Density") +
  theme_minimal()