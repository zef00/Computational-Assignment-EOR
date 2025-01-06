library(ggplot2)
library(MASS)
library(tidyverse)
library(zoo)
library(dplyr)
library(vcd)

#Datasets
ScanRecords <- read.csv("ScanRecords.csv")
ScanRecords$Date <- as.Date(ScanRecords$Date)
Type1 <- dplyr::filter(ScanRecords, PatientType=="Type 1") 
Type2 <- dplyr::filter(ScanRecords, PatientType=="Type 2") 

#Descriptives 
ScanRecords_grouped <- ScanRecords %>%
  group_by(PatientType) %>%
  summarise(count = length(PatientType))

ggplot(ScanRecords_grouped, aes(x = PatientType, y = count, fill = PatientType)) +
  geom_bar(stat = "identity") +
  labs(x = "Type", y = "Counts") +
  theme_minimal()

Type1_grouped <- Type1 %>% 
  group_by(Date) %>%
  summarise(count = length(Date))

ggplot(Type1_grouped)+
  geom_histogram(aes(x=count))

Type2_grouped <- Type2 %>% 
  group_by(Date) %>%
  summarise(count = length(Date))

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
qqline(Type1$Duration, col = "red", lwd = 2)
qqnorm(Type2$Duration) #This one is bit weird 
qqline(Type2$Duration, col = "red", lwd = 2)

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
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  labs(
    title = "Histogram of Interarrival times (Type 1)",
    x = "Interarrival times",
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

# plot histogram of interarrival times of patienttype 2 
ggplot(Type2, aes(x = InterArrivalTime)) +
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "darkblue", alpha = 0.7) +
  labs(
    title = "Histogram of Interarrival times (Type 2)",
    x = "Interarrival times",
    y = "Frequency"
  ) 

# Check of we can use poisson/exponential distribution for bootstrap (seems to NOT be the case)
gf <- goodfit(InterArrivalTimes,type= "poisson",method= "ML")
summary(gf)
ks.test(InterArrivalTimes, "pexp", rate =  1 / mean(InterArrivalTimes))

# Testing if we can use normal distribution for bootstrap (seems like it)
qqnorm(InterArrivalTimes) 
qqline(InterArrivalTimes, col = "red", lwd = 2)

# Perform bootstrapping based on exponential distribution
B <- 10000  # Number of bootstrap samples
bootstrap_means <- numeric(B)
bootstrap_sd <- numeric(B)

for (b in seq_len(B)) {
  # Generate bootstrap sample from exponential distribution
  bootstrap_sample <- rnorm(length(InterArrivalTimes), mean = mean(InterArrivalTimes), sd = sd(InterArrivalTimes))
  bootstrap_means[b] <- mean(bootstrap_sample)
  bootstrap_sd[b] <- sd(bootstrap_sample)
}

# Calculate bootstrapped mean and 95% confidence interval
bootstrapped_mean <- mean(bootstrap_means)
bootstrapped_sd <- mean(bootstrap_sd)
ci <- quantile(bootstrap_means, probs = c(0.025, 0.975))

# Print results
cat("Bootstrapped Mean Inter-Arrival Time:", bootstrapped_mean, "\n")
cat("Bootstrapped SD Inter-Arrival Time:", bootstrapped_sd, "\n")
cat("95% Confidence Interval for Mean Inter-Arrival Time (Exponential):", ci, "\n")

# Plot the bootstrap distribution
ggplot(data.frame(bootstrap_means = bootstrap_means), aes(x = bootstrap_means)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  geom_vline(xintercept = bootstrapped_mean, color = "red", linetype = "dashed") +
  geom_vline(xintercept = ci[1], color = "blue", linetype = "dotted") +
  geom_vline(xintercept = ci[2], color = "blue", linetype = "dotted") +
  labs(title = "Bootstrap Distribution of Mean Inter-Arrival Time",
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

# Plot the Interarrival time 
ggplot(ScanRecords, aes(x = InterArrivalTime)) +
  geom_histogram(binwidth = 0.1, fill = "violet", color = "purple", alpha = 0.7) +
  labs(
    title = "Histogram of Interarrival times (All Types)",
    x = "Interarrival times",
    y = "Frequency"
  ) 

#Check if we can use parametric bootstrap (doubtful as p-value = 0.06)
ScanRecords_grouped <- ScanRecords %>%
  group_by(Date) %>%
  summarise(count = length(PatientType))

gf <- goodfit(ScanRecords_grouped$count,type= "poisson",method= "ML")
summary(gf)
plot(gf,main="Count data vs Poisson distribution")

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


################################################################################
# PART II
# Old Policy = Having a facility for each type
# New Policy = Having uniform facilities
################################################################################

# Distributions PART I
# TYPE 1: Duration
# sample_mean <- 0.43266080
# standard_deviation <- 0.09777424
# quantile_2.5 <- 0.24122472
# quantile_97.5 <- 0.63157661

# TYPE 1: Inter-arrivals (Exponential)
# bootstrapped_mean <- 1.83370524886

# TYPE 2: Duration
# Non-parametric bootstrap sample used for simulation (as defined earlier)

# TYPE 2: Inter-arrivals (Exponential)
# bootstrapped_mean <- 0.8667086
# confidence_95_lower <- 0.7595675
# confidence_95_upper <- 0.9816138

# COMBINED: Duration
# sample_mean <- 0.5242738
# confidence_95_lower <- 0.5100280
# confidence_95_upper <- 0.5385197
# standard_deviation <- 0.1806868
# quantile_60 <- 0.3051035
# quantile_75 <- 0.6118824

# COMBINED: Inter-arrivals (Exponential)
# bootstrapped_mean <- 0.3345074
# confidence_95_lower <- 0.3086699
# confidence_95_upper <- 0.3613692

######################### RESULTS FROM PART I ##################################
set.seed(777)
# TYPE 1
# Duration
mu_type1 <- 0.43266080
sig_type1 <- 0.09777424
var_type1 <- sig_type1**2
# Inter-arrivals
lambda_type1 <- 1.83370524886
mean_type1 <- 1 / lambda_type1

# TYPE 2
# Duration
# Non-parametric bootstrap sample is used for Type 2 durations
bootstrap_durations_type2 <- replicate(1000, {
  sample(Type2$Duration, size = nrow(Type2), replace = TRUE)
})

# Inter-arrivals
lambda_type2 <- 0.8667086
mean_type2 <- 1 / lambda_type2

# INITIALIZING THE (FIXED) TIME SLOTS
z_score <- qnorm(0.95) # One-tailed 95th percentile score, i.e., P(Z <= z) = 0.95
duration_95_type1 <- mu_type1 + z_score * sig_type1 # z_score = (X - mu_type1) / sig_type1
duration_95_type2 <- max(bootstrap_durations_type2)  # Use bootstrap quantiles
duration_95_type1_min <- duration_95_type1 * 60 # 35.6 is approximately 35 min
duration_95_type2_min <- duration_95_type2 * 60 # Example duration quantile
time_slot_type1 <- 35 / 60
time_slot_type2 <- 60 / 60

# DISCRETE EVENT SIMULATION
# Policy: both
# Parameters
num_facilities <- 2
sim_period <- 23 # Approximately one month
start_time <- 8.00 # Every day starts at 8.00h
end_time <- 17.00 # And after 17.00h considered 'overtime'

# Policy: both
# Generating patient calls (input: certain day of the 30-day month -> t+1)
generate_calls <- function(lambda, day, patient_type) {
  calls <- c()
  time <- start_time
  while (time < end_time) { # Call center closes at 17.00h
    inter_arrival <- rexp(1, rate = lambda) # Time between arrivals \sim Exp(\lambda)
    time <- time + inter_arrival # Update time of the day variable
    if (time < end_time) { # Add call if happens before closing hour
      calls <- c(calls, time)
    }
  }
  data.frame(
    Day = day,
    CallTime = calls,
    PatientType = patient_type
  )
}

# Policy: both
# Combining daily calls from type I and type II (input: certain day of the 30-day month)
generate_daily_calls <- function(day) {
  calls_type1 <- generate_calls(lambda_type1, day, "Type1")
  calls_type2 <- generate_calls(lambda_type2, day, "Type2")
  calls <- rbind(calls_type1, calls_type2) # rbind() function adds type II calls below type I calls
  return(calls)
}

# Policy: Old
# Assigning patients to facilities and time slots
schedule_appointments_old <- function(calls) {
  # Order patients w.r.t. their calling times
  calls <- calls[order(calls$CallTime), ]
  # Initialize schedules for each facility
  schedule <- data.frame( # Initializing a data frame for a daily facility schedule
    PatientID = seq_len(nrow(calls)),
    PatientType = calls$PatientType,
    CallDay = calls$Day,
    ScheduledDay = calls$Day + 1, # Patients have to be scheduled the next day
    ScheduledStart = NA,
    Facility = NA
  )
  # Split patients by type
  idx_type1 <- which(schedule$PatientType == "Type1")
  idx_type2 <- which(schedule$PatientType == "Type2")
  # Assign facilities
  schedule$Facility[idx_type1] <- 1
  schedule$Facility[idx_type2] <- 2
  # For facility 1, schedule Type1 patients sequentially
  if (length(idx_type1) > 0) {
    start_time_vector <- rep(start_time, length(idx_type1))
    for (i in seq_along(idx_type1)) {
      if (i > 1) {
        start_time_vector[i] <- start_time_vector[i-1] + time_slot_type1
      }
      schedule$ScheduledStart[idx_type1[i]] <- start_time_vector[i]
    }
  }
  # For facility 2, schedule Type2 patients sequentially
  if (length(idx_type2) > 0) {
    start_time_vector <- rep(start_time, length(idx_type2))
    for (i in seq_along(idx_type2)) {
      if (i > 1) {
        start_time_vector[i] <- start_time_vector[i-1] + time_slot_type2
      }
      schedule$ScheduledStart[idx_type2[i]] <- start_time_vector[i]
    }
  }
  schedule
}

# Policy: New
# Assigning patients to facilities and time slots
schedule_appointments_new <- function(calls) {
  calls <- calls[order(calls$CallTime), ]
  schedule <- data.frame(
    PatientID = seq_len(nrow(calls)),
    PatientType = calls$PatientType,
    CallDay = calls$Day,
    ScheduledDay = calls$Day + 1,
    ScheduledStart = NA,
    Facility = NA
  )
  
  # Track next available time for each facility
  facility_next_time <- rep(start_time, num_facilities)
  
  for (i in seq_len(nrow(schedule))) {
    ptype <- schedule$PatientType[i]
    # Decide slot length
    slot_len <- if (ptype == "Type1") time_slot_type1 else time_slot_type2
    
    # Find facility that is available first
    f <- which.min(facility_next_time)  # facility index
    schedule$Facility[i] <- f
    schedule$ScheduledStart[i] <- facility_next_time[f]
    
    # Update that facility's next available time
    facility_next_time[f] <- facility_next_time[f] + slot_len
  }
  
  schedule
}

# The rest of the simulation metrics, day simulation, and policy analysis functions remain unchanged.

# Policy: both
# Simulate scanning process
simulate_day <- function(schedule) {
  # Facility availability tracking
  facility_end_time <- rep(start_time, num_facilities)
  
  # Add columns
  schedule$ActualDuration <- NA
  schedule$ScanStart <- NA
  schedule$ScanEnd <- NA
  
  # Order patients by scheduled start time
  schedule <- schedule[order(schedule$ScheduledStart), ]
  
  # Bootstrap index for Type 2 duration
  bootstrap_index <- sample(1:ncol(bootstrap_durations_type2), size = nrow(schedule), replace = TRUE)
  
  for (i in seq_len(nrow(schedule))) {
    f <- schedule$Facility[i]
    ptype <- schedule$PatientType[i]
    
    # Generate random actual duration
    if (ptype == "Type1") {
      duration <- rnorm(1, mean = mu_type1, sd = sig_type1)  # Keep Type 1 as normal distribution
    } else if (ptype == "Type2") {
      # Use the bootstrap sample for Type 2 durations
      duration <- bootstrap_durations_type2[sample(1:nrow(bootstrap_durations_type2), 1), bootstrap_index[i]]
    }
    
    duration <- max(duration, 0)  # Ensure no negative durations
    
    # Actual scan start is the max of scheduled start or facility availability
    scan_start <- max(schedule$ScheduledStart[i], facility_end_time[f])
    schedule$ScanStart[i] <- scan_start
    schedule$ActualDuration[i] <- duration
    schedule$ScanEnd[i] <- scan_start + duration
    
    # Update facility availability
    facility_end_time[f] <- schedule$ScanEnd[i]
  }
  
  # Return in original patient order
  schedule[order(schedule$PatientID), ]
}


# Policy: both
# Our valuation variables
calculate_metrics <- function(schedule) {
  # If schedule is empty, return zeros and logical defaults
  if (nrow(schedule) == 0) {
    return(list(
      Overtime = 0,
      Utilization = 0,
      AvgWaiting = 0,
      FinishedOnTime = TRUE
    ))
  }
  facility_finish <- numeric(num_facilities)
  for (f in seq_len(num_facilities)) {
    # Rows for this facility
    idx <- which(schedule$Facility == f)
    if (length(idx) == 0) {
      # No patients on this facility => finishing time = 0
      facility_finish[f] <- start_time  # or 0
    } else {
      facility_finish[f] <- max(schedule$ScanEnd[idx])
    }
  }
  raw_overtime <- facility_finish - end_time
  day_overtime <- pmax(0, raw_overtime)  # vector of length num_facilities
  
  total_scan_time <- sum(schedule$ActualDuration)
  total_capacity <- (end_time - start_time) * num_facilities
  utilization <- total_scan_time / total_capacity
  
  waiting_times <- schedule$ScanStart - schedule$ScheduledStart
  avg_wait <- mean(waiting_times, na.rm = TRUE)
  
  # Determine if all scans finished on time
  finished_on_time <- all(facility_finish <= end_time)
  
  # Return list of metrics
  list(
    Overtime = sum(day_overtime),
    Utilization = utilization,
    AvgWaiting = avg_wait,
    FinishedOnTime = finished_on_time
  )
}

# Policy: both
# SIMULATING A DAY
simulate_one_day <- function(day, policy = c("Old", "New")) {
  calls_today <- generate_daily_calls(day)
  
  # If no calls, return zero metrics
  if (nrow(calls_today) == 0) {
    return(list(
      Overtime = 0, Utilization = 0,
      AvgWaiting = 0, FinishedOnTime = TRUE
    ))
  }
  
  if (policy == "Old") {
    sched <- schedule_appointments_old(calls_today)
  } else {
    sched <- schedule_appointments_new(calls_today)
  }
  
  # Now simulate
  sched_sim <- simulate_day(sched)
  
  # Return metrics
  calculate_metrics(sched_sim)
}

# SIMULATING MULPTIPLE DAYS
run_simulation <- function(num_days = 23, policy = c("Old", "New"), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Collect daily metrics
  daily_results <- data.frame(
    Day = integer(),
    Overtime = numeric(),
    Utilization = numeric(),
    AvgWaiting = numeric(),
    FinishedOnTime = logical(),
    stringsAsFactors = FALSE
  )
  
  for (d in seq_len(num_days)) {
    day_metrics <- simulate_one_day(d, policy = policy)
    daily_results[d, ] <- c(
      d,
      day_metrics$Overtime,
      day_metrics$Utilization,
      day_metrics$AvgWaiting,
      day_metrics$FinishedOnTime
    )
  }
  daily_results
}

# Example: single run
results_old <- run_simulation(num_days = 23, policy = "Old", seed = 777)
results_new <- run_simulation(num_days = 23, policy = "New", seed = 777)

# Summarize
summary_old <- data.frame(
  Overtime = mean(results_old$Overtime),
  Utilization = mean(results_old$Utilization),
  AvgWaiting = mean(results_old$AvgWaiting),
  FinishedOnTimePct = mean(results_old$FinishedOnTime) * 100
)
summary_new <- data.frame(
  Overtime = mean(results_new$Overtime),
  Utilization = mean(results_new$Utilization),
  AvgWaiting = mean(results_new$AvgWaiting),
  FinishedOnTimePct = mean(results_new$FinishedOnTime) * 100
)

summary_old
summary_new

# SIMULATING MULTIPLE DAYS AND REPLICATIONS
run_replications <- function(num_days, policy, R) {
  metrics_list <- list()
  for (r in 1:R) {
    # Each replication can have a different seed or be random
    daily_res <- run_simulation(num_days, policy)
    # Store the average over the daily run
    metrics_list[[r]] <- c(
      Overtime = mean(daily_res$Overtime),
      Utilization = mean(daily_res$Utilization),
      AvgWaiting = mean(daily_res$AvgWaiting),
      FinishedOnTimePct = mean(daily_res$FinishedOnTime)*100
    )
  }
  # Combine
  mat <- do.call(rbind, metrics_list)
  # Return data frame of replication results
  data.frame(mat)
}

# Example usage:
rep_old <- run_replications(23, policy = "Old", R = 4)
rep_new <- run_replications(23, policy = "New", R = 4)

# Summaries
apply(rep_old, 2, mean)
apply(rep_new, 2, mean)

set.seed(777)
# Code snippet the look at a day
policy = "Old"
inspect_one_day <- function(day, policy) {
  calls_today <- generate_daily_calls(day)
  if (policy == "Old") {
    sched <- schedule_appointments_old(calls_today)
  } else {
    sched <- schedule_appointments_new(calls_today)
  }
  simulate_day(sched)
}
schedule_day_5 <- inspect_one_day(day = 1, policy)
print(schedule_day_5)

set.seed(777)
# Code snippet to check the zero call/arrival day(s)
inspect_all_days <- function(num_days = 23, policy = c("Old", "New")) {
  for (d in seq_len(num_days)) {
    # Generate calls for day d
    calls_today <- generate_daily_calls(d)
    
    # Print how many calls were generated
    cat("\n--- Day", d, "---\n")
    cat("Number of calls: ", nrow(calls_today), "\n")
    if (nrow(calls_today) > 0) {
      # Schedule according to the chosen policy
      if (policy == "Old") {
        sched <- schedule_appointments_old(calls_today)
      } else {
        sched <- schedule_appointments_new(calls_today)
      }
      
      # Simulate scanning
      sched_sim <- simulate_day(sched)
      
      # Print the final schedule (after simulating actual durations)
      cat("Final Schedule:\n")
      print(sched_sim)
      
    } else {
      cat("No arrivals on day", d, "!\n")
    }
  }
}

# Let's say we want to see all days (1..30) for the Old policy:
inspect_all_days(num_days = 23, policy = "Old")

