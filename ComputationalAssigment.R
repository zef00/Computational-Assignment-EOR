library(ggplot2)
library(MASS)
library(tidyverse)
library(zoo)
library(dplyr)

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


##IMPORTANT METRICS FOR DISCRETE EVENT SIMULATION (TYPE 1 PATIENT SCAN DURATION)##
# Calculate mean and standard deviation
mean_duration <- mean(Type1$Duration)
std_dev <- sd(Type1$Duration)
n <- length(Type1$Duration)
std_error <- std_dev / sqrt(n)


# Quantiles
duration_quantiles <- quantile(Type1$Duration, probs = c(0.025, 0.975))

# Create a summary table
duration_type1_summary <- tibble(
  Metric = c(
    "Sample Mean",
    "Standard Deviation",
    "Quantile (2.5%)",
    "Quantile (97.5%)"
  ),
  Value = c(
    mean_duration,
    std_dev,
    duration_quantiles[1],
    duration_quantiles[2]
  )
)

# Print the summary table
print(duration_type1_summary)


ggplot(Type1, aes(x = Duration)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black", alpha = 0.7) +
  geom_vline(xintercept = duration_quantiles[1], color = "blue", linetype = "dashed", size = 1) +
  geom_vline(xintercept = duration_quantiles[2], color = "blue", linetype = "dashed", size = 1) +
  labs(
    title = "Histogram of Scan Durations (Type 1)",
    x = "Duration",
    y = "Frequency"
  ) +
  annotate(
    "text",
    x = duration_quantiles[1],
    y = max(table(cut(Type1$Duration, breaks = 20))),
    label = paste0("2.5% Quantile: ", round(duration_quantiles[1], 2)),
    hjust = -0.1,
    color = "blue"
  ) +
  annotate(
    "text",
    x = duration_quantiles[2],
    y = max(table(cut(Type1$Duration, breaks = 20))),
    label = paste0("97.5% Quantile: ", round(duration_quantiles[2], 2)),
    hjust = -0.1,
    color = "blue"
  )


### INTERARRIVAL TIMES FOR TYPE 1###
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

# Calculate statistical metrics
mean_interarrival <- mean(InterArrivalTimes)

# Rate parameter for exponential distribution
lambda_rate <- 1 / mean_interarrival


ggplot(Type1, aes(x = InterArrivalTime)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black", alpha = 0.7) +
  labs(
    title = "Histogram of Interarrivaltimes (Type 1)",
    x = "Interarrivaltimes",
    y = "Frequency"
  ) 




#Bootstrapping the mean duration for type 2 
X <- Type2$Duration
B <- 1000 

n <- length(X)                                          # Sample size
X.bar <- mean(X)                                        # Sample mean of X
St.Dev <- sd(X)                                         # Standard deviation of X

X.star_mean <- rep(NA, n)
X.star_sd <- rep(NA,n)
X.star_low <- rep(NA, n)      # 0.025 bound
X.star_high <- rep(NA, n)     # 0.975 bound
for (b in 1:B) {
  J <- sample.int(n, size = n, replace = TRUE)            # Draw the indices J
  X.star <- X[J]                                      # Draw the bootstrap sample
  X.star_mean[b] <- mean(X.star) 
  X.star_sd[b] <- sd(X.star)
  X.star_high[b] <- quantile(X.star, probs = 0.975)
  X.star_low[b] <- quantile(X.star, probs = 0.025)
}

mean_duration_type2 <- mean(X.star_mean) #these are the variables of interest to us
SD_duration_type2 <- mean(X.star_sd)
duration2_quantiles <- c(mean(X.star_low), mean(X.star_high))

duration_type2_summary <- tibble(
  Metric = c(
    "Sample Mean",
    "95% Confidence Interval (Lower)",
    "95% Confidence Interval (Upper)",
    "Standard Deviation",
    "Quantile (2.5%)",
    "Quantile (97.5%)"
  ),
  Value = c(
    mean_duration_type2,
    mean_duration_type2-1.96*SD_duration_type2/sqrt(n),
    mean_duration_type2+1.96*SD_duration_type2/sqrt(n),
    SD_duration_type2,
    duration2_quantiles[1],
    duration2_quantiles[2]
  )
)

# Monte Carlo Study for the type 2 bootstrap mean (this can take long)
N = 2000
B = 500 
bias_mean = rep(NA, N)
bias_sd = rep(NA, N)
for (i in 1:N){
  mc_type2_duration <- rnorm(length(Type2$Duration),
                             mean = mean_duration_type2,
                             sd = SD_duration_type2)
  mc_mean <- mean(mc_type2_duration)
  mc_sd <- sd(mc_type2_duration)
  X.star_mean <- rep(NA, length(Type2$Duration))
  X.star_sd <- rep(NA, length(Type2$Duration))
  for (b in 1:B){
    J <- sample.int(length(Type2$Duration), size = length(Type2$Duration), replace = TRUE)
    X.star <- mc_type2_duration[J]                                      # Draw the bootstrap sample
    X.star_mean[b] <- mean(X.star) 
    X.star_sd[b] <- sd(X.star)
  }
  MCbootstrap_mean <- mean(X.star_mean)
  MCbootstrap_sd <- mean(X.star_sd)
  bias_mean[i] <- (MCbootstrap_mean - mc_mean)^2
  bias_sd[i] <- (MCbootstrap_sd - mc_sd)^2

}
avg_bias_mean = mean(bias_mean) #very low so should be good 
avg_bias_sd = mean(bias_sd)

## BOOTSTRAP patient 2 arrival times with interarival times

# Calculate inter-arrival times
Type2 <- Type2 %>%
  group_by(Date) %>%
  mutate(
    InterArrivalTime = Time - lag(Time),  # Standard inter-arrival time within the same day
    FirstOfDay = row_number() == 1        # Flag for the first record of the day
  ) %>%
  ungroup()

# Handle first times of each day
Type2 <- Type2 %>%
  mutate(
    InterArrivalTime = ifelse(
      FirstOfDay,
      # If first record of the day, calculate time difference with lagged time + 9 hours (17,00 - 8,00 = 9,00)
      Time - lag(Time) + 9,
      InterArrivalTime
    )
  )

# Remove remaining NA values (the first time of the first day)
InterArrivalTimes <- Type2$InterArrivalTime[!is.na(Type2$InterArrivalTime)]

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


# Bootstrap duration for all the patients together
X <- ScanRecords$Duration
B <- 1000 

n <- length(X)                                          # Sample size
X.bar <- mean(X)                                        # Sample mean of X
St.Dev <- sd(X)                                         # Standard deviation of X

X.star_mean <- rep(NA, n)
X.star_sd <- rep(NA,n)
X.star_low <- rep(NA, n)      # 0.025 bound
X.star_high <- rep(NA, n)     # 0.975 bound
for (b in 1:B) {
  J <- sample.int(n, size = n, replace = TRUE)            # Draw the indices J
  X.star <- X[J]                                      # Draw the bootstrap sample
  X.star_mean[b] <- mean(X.star) 
  X.star_sd[b] <- sd(X.star)
  X.star_high[b] <- quantile(X.star, probs = 0.75)
  X.star_low[b] <- quantile(X.star, probs = 0.06)
}

mean_duration_all <- mean(X.star_mean) #these are the variables of interest to us
SD_duration_all <- mean(X.star_sd)
duration_all_quantiles <- c(mean(X.star_low), mean(X.star_high))

duration_all_summary <- tibble(
  Metric = c(
    "Sample Mean",
    "95% Confidence Interval (Lower)",
    "95% Confidence Interval (Upper)",
    "Standard Deviation",
    "Quantile (60%)",
    "Quantile (75%)"
  ),
  Value = c(
    mean_duration_all,
    mean_duration_all-1.96*SD_duration_all/sqrt(n),
    mean_duration_all+1.96*SD_duration_all/sqrt(n),
    SD_duration_all,
    duration_all_quantiles[1],
    duration_all_quantiles[2]
  )
)

# Monte Carlo Study for the duration bootstrap mean and standard deviation (this can take long)
N = 2000
B = 500 
bias_mean = rep(NA, N)
bias_sd = rep(NA, N)
for (i in 1:N){
  mc_type_duration <- rnorm(length(ScanRecords$Duration),
                             mean = mean_duration_all,
                             sd = SD_duration_all)
  mc_mean <- mean(mc_type_duration)
  mc_sd <- sd(mc_type_duration)
  X.star_mean <- rep(NA, length(ScanRecords$Duration))
  X.star_sd <- rep(NA, length(ScanRecords$Duration))
  for (b in 1:B){
    J <- sample.int(length(ScanRecords$Duration), size = length(ScanRecords$Duration), replace = TRUE)
    X.star <- mc_type_duration[J]                                      # Draw the bootstrap sample
    X.star_mean[b] <- mean(X.star) 
    X.star_sd[b] <- sd(X.star)
  }
  MCbootstrap_mean <- mean(X.star_mean)
  MCbootstrap_sd <- mean(X.star_sd)
  bias_mean[i] <- (MCbootstrap_mean - mc_mean)^2
  bias_sd[i] <- (MCbootstrap_sd - mc_sd)^2
  
}
avg_bias_mean = mean(bias_mean) #very low so should be good 
avg_bias_sd = mean(bias_sd)


## BOOTSTRAP both patients interarrival times 

# Calculate inter-arrival times for all patients
ScanRecords <- ScanRecords %>%
  arrange(Date, Time) %>%  # Ensure records are sorted by date and time
  group_by(Date) %>%
  mutate(
    InterArrivalTime = Time - lag(Time),  # Standard inter-arrival time within the same day
    FirstOfDay = row_number() == 1        # Flag for the first record of the day
  ) %>%
  ungroup()

# Handle first times of each day
ScanRecords <- ScanRecords %>%
  mutate(
    InterArrivalTime = ifelse(
      FirstOfDay,
      # If first record of the day, calculate time difference with lagged time + 9 hours (17,00 - 8,00 = 9,00)
      Time - lag(Time) + 9,
      InterArrivalTime
    )
  )

# Remove remaining NA values (the first time of the first day)
InterArrivalTimes <- ScanRecords$InterArrivalTime[!is.na(ScanRecords$InterArrivalTime)]

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

# Monte Carlo Simulation for Inter-Arrival Times (All Patients)
N <- 2000  # Number of Monte Carlo iterations
B <- 500  # Number of bootstrap samples in each Monte Carlo iteration

# Initialize vectors to store biases
bias_mean <- numeric(N)
bias_sd <- numeric(N)

# Perform Monte Carlo simulation
for (i in seq_len(N)) {
  # Generate synthetic data from exponential distribution
  mc_interarrival_times <- rexp(length(InterArrivalTimes), rate = lambda)
  
  # Calculate true mean and standard deviation of synthetic data
  true_mean <- mean(mc_interarrival_times)
  true_sd <- sd(mc_interarrival_times)
  
  # Bootstrap the synthetic data
  bootstrap_means <- numeric(B)
  bootstrap_sds <- numeric(B)
  
  for (b in seq_len(B)) {
    # Generate bootstrap sample
    bootstrap_sample <- sample(mc_interarrival_times, replace = TRUE)
    
    # Calculate mean and standard deviation of bootstrap sample
    bootstrap_means[b] <- mean(bootstrap_sample)
    bootstrap_sds[b] <- sd(bootstrap_sample)
  }
  
  # Calculate average bootstrap mean and standard deviation
  bootstrap_mean <- mean(bootstrap_means)
  bootstrap_sd <- mean(bootstrap_sds)
  
  # Calculate squared bias for mean and standard deviation
  bias_mean[i] <- (bootstrap_mean - true_mean)^2
  bias_sd[i] <- (bootstrap_sd - true_sd)^2
}

# Calculate average bias for mean and standard deviation
avg_bias_mean <- mean(bias_mean) # again very low I guess
avg_bias_sd <- mean(bias_sd)

# Print results
cat("Monte Carlo Simulation Results:\n")
cat("Average Squared Bias of Mean:", avg_bias_mean, "\n")
cat("Average Squared Bias of Standard Deviation:", avg_bias_sd, "\n")



