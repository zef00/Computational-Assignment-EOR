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
    "Sample Mean of scanduration Type 1",
    "Standard Deviation of scanduration Type 1",
    "Quantile (2.5%) of scanduration Type 1",
    "Quantile (97.5%) of scanduration Type 1"
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
InterArrivalTimes1 <- Type1$InterArrivalTime[!is.na(Type1$InterArrivalTime)]

# Calculate statistical metrics
mean_interarrival <- mean(InterArrivalTimes1)

# Rate parameter for exponential distribution
lambda_rate <- 1 / mean_interarrival

interarrivaltime1_summary <- tibble(
  Metric = c(
    "Mean Interarrival time Type 1",
    "lambda"
  ),
  Value = c(
    mean_interarrival,
    lambda_rate
  )
)

ggplot(Type1, aes(x = InterArrivalTime)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "darkgreen", alpha = 0.7) +
  labs(
    title = "Histogram of Interarrival times (Type 1)",
    x = "Interarrival times",
    y = "Frequency"
  ) 


#Bootstrapping the mean duration for type 2 
X <- Type2$Duration
B <- 2000 

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
    "Mean scanduration Type 2",
    "95% Confidence Interval (Lower) of mean scanduration Type 2",
    "95% Confidence Interval (Upper) of mean scanduration Type 2",
    "Standard Deviation of scanduration Type 2",
    "Quantile (2.5%) of scanduration Type 2",
    "Quantile (97.5%) of scanduration Type 2"
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

# Perform bootstrapping based on normal distribution
B <- 10000  # Number of bootstrap samples
bootstrap_means <- numeric(B)
bootstrap_sd <- numeric(B)

for (b in seq_len(B)) {
  # Generate bootstrap sample from normal distribution
  bootstrap_sample <- rnorm(length(InterArrivalTimes), mean = mean(InterArrivalTimes), sd = sd(InterArrivalTimes))
  bootstrap_means[b] <- mean(bootstrap_sample)
  bootstrap_sd[b] <- sd(bootstrap_sample)
}

# Calculate bootstrapped mean and 95% confidence interval
bootstrapped_mean <- mean(bootstrap_means)
bootstrapped_sd <- mean(bootstrap_sd)
ci <- quantile(bootstrap_means, probs = c(0.025, 0.975))

interarrivaltimes2_summary <- tibble(
  Metric = c(
    "Mean Interarrival time Type 2",
    "95% Confidence Interval (Lower) of mean Interarrival time Type 2",
    "95% Confidence Interval (Upper) of mean Interarrival time Type 2",
    "Standard Deviation of Interarrival time Type 2"
  ),
  Value = c(
    bootstrapped_mean,
    bootstrapped_mean -1.96*bootstrapped_sd/sqrt(length(InterArrivalTimes)),
    bootstrapped_mean +1.96*bootstrapped_sd/sqrt(length(InterArrivalTimes)),
    bootstrapped_sd
  )
)

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


# Overview of the parameter tables 
duration_type1_summary
interarrivaltime1_summary
duration_type2_summary
interarrivaltimes2_summary

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
# sample_mean <- 0.6689911
# confidence_95_lower <- 0.6453036
# confidence_95_upper <- 0.6926785
# standard_deviation <- 0.1868363
# quantile_2.5 <- 0.3634360
# quantile_97.5 <- 1.0652090

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
mu_type1 <- 0.43266080
sig_type1 <- 0.09777424
var_type1 <- sig_type1^2

# Inter-arrivals (Type1)
lambda_type1 <- 1.83370524886  # rate
mean_type1   <- 1 / lambda_type1

# TYPE 2
# mu_type2 <- 0.6689911
# sig_type2 <- 0.1868363
# var_type2 <- sig_type2^2

# Inter-arrivals (Type2)
# lambda_type2 <- 0.8667086 # exp. distr.
# mean_type2   <- 1 / lambda_type2
mean_int_type2 <- 0.8663962 # normal distr.
sig_int_type2 <- 0.3103057
var_int_type2 <- sig_int_type2^2

# INITIALIZING THE (FIXED) TIME SLOTS
z_score <- qnorm(0.95) # One-tailed 95th percentile => P(Z <= z) = 0.95
duration_95_type1 <- mu_type1 + z_score * sig_type1
duration_95_type2 <- mu_type2 + z_score * sig_type2
duration_95_type1_min <- duration_95_type1 * 60
duration_95_type2_min <- duration_95_type2 * 60
time_slot_type1 <- 29/60
time_slot_type2 <- 48/60

# DISCRETE EVENT SIMULATION
# Parameters
num_facilities <- 2
sim_period <- 23 # ~ one month without weekend days
start_time <- 8.00  # working day start
end_time   <- 17.00 # working day end

############################## PART II #########################################
# Generating patient calls
generate_calls_type1 <- function(lambda, day, patient_type) {
  calls <- c()
  time <- start_time
  while (time < end_time) {
    inter_arrival <- rexp(1, rate = lambda)
    time <- time + inter_arrival
    if (time < end_time) {
      calls <- c(calls, time)
    }
  }
  data.frame(
    Day         = day,
    CallTime    = calls,
    PatientType = patient_type
  )
}

generate_calls_type2 <- function(mean_int_type2, sig_int_type2, day, patient_type) {
  calls <- c()
  time <- start_time
  while (time < end_time) {
    inter_arrival <- sample(InterArrivalTimes, 1, replace = TRUE)
    time <- time + inter_arrival
    if (time < end_time) {
      calls <- c(calls, time)
    }
  }
  data.frame(
    Day         = day,
    CallTime    = calls,
    PatientType = patient_type
  )
}

# Combine Type1 + Type2 calls
generate_daily_calls <- function(day) {
  calls_type1 <- generate_calls_type1(lambda_type1, day, "Type1")
  calls_type2 <- generate_calls_type2(mean_int_type2, sig_int_type2, day, "Type2")
  rbind(calls_type1, calls_type2)
}

# Scheduling (Old Policy)
schedule_appointments_old <- function(calls) {
  if (nrow(calls) == 0) return(data.frame()) # the low prob. case of no calls
  
  calls <- calls[order(calls$CallTime), ] # order the calls according to time of calling
  # creating the schedule with the information columns below
  schedule <- data.frame(
    PatientID = seq_len(nrow(calls)), # ID given 1, 2, 3, ...
    PatientType = calls$PatientType,
    CallDay     = calls$Day,
    ScheduledDay= calls$Day + 1, # have to be scheduled next day
    ScheduledStart = NA,
    Facility      = NA
  )
  
  idx_type1 <- which(schedule$PatientType == "Type1") # group of type 1 patients
  idx_type2 <- which(schedule$PatientType == "Type2") # group of type 2 patients
  
  schedule$Facility[idx_type1] <- 1 # assigns to the right facility (old policy)
  schedule$Facility[idx_type2] <- 2 # ''
  
  # allocates the type 1 patients to timeslots according to the chosen timeslot length
  if (length(idx_type1) > 0) {
    start_time_vector <- rep(start_time, length(idx_type1))
    for (i in seq_along(idx_type1)) {
      if (i > 1) {
        start_time_vector[i] <- start_time_vector[i - 1] + time_slot_type1
      }
      schedule$ScheduledStart[idx_type1[i]] <- start_time_vector[i]
    }
  }
  # allocates the type 2 patients to timeslots according to the chosen timeslot length
  if (length(idx_type2) > 0) {
    start_time_vector <- rep(start_time, length(idx_type2))
    for (i in seq_along(idx_type2)) {
      if (i > 1) {
        start_time_vector[i] <- start_time_vector[i - 1] + time_slot_type2
      }
      schedule$ScheduledStart[idx_type2[i]] <- start_time_vector[i]
    }
  }
  schedule
}

# Scheduling (New Policy)
schedule_appointments_new <- function(calls) {
  if (nrow(calls) == 0) return(data.frame()) # the low prob. case of no calls
  
  calls <- calls[order(calls$CallTime), ] # order the calls according to time of calling
  # creating the schedule with the information columns below
  schedule <- data.frame(
    PatientID     = seq_len(nrow(calls)),
    PatientType   = calls$PatientType,
    CallDay       = calls$Day,
    ScheduledDay  = calls$Day + 1,
    ScheduledStart= NA,
    Facility      = NA
  )
  
  facility_next_time <- rep(start_time, num_facilities)
  
  for (i in seq_len(nrow(schedule))) {
    ptype <- schedule$PatientType[i]
    slot_len <- if (ptype == "Type1") time_slot_type1 else time_slot_type2
    f <- which.min(facility_next_time)
    schedule$Facility[i] <- f
    schedule$ScheduledStart[i] <- facility_next_time[f]
    facility_next_time[f] <- facility_next_time[f] + slot_len
  }
  schedule
}

# Simulate a day's scanning
simulate_day <- function(schedule) {
  if (nrow(schedule) == 0) return(schedule) # the low prob. case of no calls
  
  facility_end_time <- rep(start_time, num_facilities) # since identical facilities end time can be treated as twice as long
  
  schedule$ActualDuration <- NA
  schedule$ScanStart <- NA
  schedule$ScanEnd   <- NA
  schedule <- schedule[order(schedule$ScheduledStart), ]
  
  for (i in seq_len(nrow(schedule))) {
    f <- schedule$Facility[i]
    ptype <- schedule$PatientType[i]
    if (ptype == "Type1") {
      duration <- rnorm(1, mean = mu_type1, sd = sig_type1)
    } else {
      duration <- sample(Type2$Duration, size = 1, replace = TRUE)
      # duration <- rnorm(1, mean = mu_type2, sd = sig_type2)
    }
    duration <- max(duration, 0)
    scan_start <- max(schedule$ScheduledStart[i], facility_end_time[f])
    schedule$ScanStart[i] <- scan_start
    schedule$ActualDuration[i] <- duration
    schedule$ScanEnd[i] <- scan_start + duration
    
    facility_end_time[f] <- schedule$ScanEnd[i]
  }
  schedule[order(schedule$PatientID), ]
}

# This function calculates per-facility, per-day metrics
calculate_per_facility_metrics <- function(day_schedule) {
  # If no patients, return 0-rows
  if (nrow(day_schedule) == 0) {
    return(data.frame(
      Facility             = integer(0),
      Day                  = integer(0),
      OvertimeHours        = numeric(0),
      OvertimeCount        = integer(0),  # 1 if overtime > 0, else 0
      OvertimeCountExtreme = integer(0),  # 1 if finishing > 19 (chosen)
      ScanTime            = numeric(0),
      Capacity            = numeric(0),   # capacity = 9 hours for each facility
      TotalWait           = numeric(0),
      NumPatients         = integer(0),
      WaitExtremes        = integer(0)
    ))
  }
  # We want a row per facility
  out_list <- list()
  facilities <- unique(day_schedule$Facility)
  for (f in facilities) {
    sched_f <- subset(day_schedule, Facility == f)
    if (nrow(sched_f) == 0) {
      # Means no patients for facility f => 0 usage
      out_list[[f]] <- data.frame(
        Facility             = f,
        Day                  = sched_f$ScheduledDay[1], # or NA or the day
        OvertimeHours        = 0,
        OvertimeCount        = 0,
        OvertimeCountExtreme = 0,
        ScanTime             = 0,
        Capacity             = (end_time - start_time), 
        TotalWait            = 0,
        NumPatients          = 0,
        WaitExtremes         = 0
      )
    } else {
      # Compute finishing times
      last_finish <- max(sched_f$ScanEnd)
      raw_overtime <- last_finish - 17
      overtime_hrs <- max(0, raw_overtime)
      # 1 if overtime>0
      ot_count <- if (overtime_hrs>0) 1 else 0
      # 1 if finishing > 19
      ot_extreme <- if (last_finish>19) 1 else 0
      
      # sum of scan times => used for utilization
      scan_time <- sum(sched_f$ActualDuration)
      capacity  <- (end_time - start_time)  # 9 hours for the day, per facility
      
      # waiting times
      wtimes <- sched_f$ScanStart - sched_f$ScheduledStart
      total_wait  <- sum(wtimes)
      n_patients  <- nrow(sched_f)
      w_extremes  <- sum(wtimes > (30/60)) # halfuurke
      
      out_list[[f]] <- data.frame(
        Facility             = f,
        Day                  = sched_f$ScheduledDay[1],
        OvertimeHours        = overtime_hrs,
        OvertimeCount        = ot_count,
        OvertimeCountExtreme = ot_extreme,
        ScanTime             = scan_time,
        Capacity             = capacity,
        TotalWait            = total_wait,
        NumPatients          = n_patients,
        WaitExtremes         = w_extremes
      )
    }
  }
  do.call(rbind, out_list)
}

# simulate_one_day => schedule and metrics per facility
simulate_one_day <- function(day, policy = c("Old","New")) {
  calls_today <- generate_daily_calls(day)
  if (nrow(calls_today) == 0) {
    # Return 2 facility-rows with zero usage, or 0 rows if you prefer
    # We'll produce 2 rows so the structure is consistent
    df_empty <- data.frame(
      Facility             = 1:2,
      Day                  = day,
      OvertimeHours        = 0,
      OvertimeCount        = 0,
      OvertimeCountExtreme = 0,
      ScanTime             = 0,
      Capacity             = (end_time - start_time),
      TotalWait            = 0,
      NumPatients          = 0,
      WaitExtremes         = 0
    )
    return(df_empty)
  }
  
  if (policy == "Old") {
    sched <- schedule_appointments_old(calls_today)
  } else {
    sched <- schedule_appointments_new(calls_today)
  }
  sched_sim <- simulate_day(sched)
  
  # Return a data.frame with 1 row per facility
  calculate_per_facility_metrics(sched_sim)
}

# run_simulation => combine daily facility-level rows
run_simulation <- function(num_days = 23, policy=c("Old","New"), seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # We'll store each day's facility stats in a big data.frame
  all_stats <- data.frame()
  
  for (d in seq_len(num_days)) {
    day_stats <- simulate_one_day(d, policy)
    all_stats <- rbind(all_stats, day_stats)
  }
  all_stats
}

# Summarize the entire simulation (single run) => single final metrics per facility
summarize_simulation <- function(sim_data) {
  # sim_data has columns: Facility, Day, OvertimeHours, OvertimeCount,
  #                      OvertimeCountExtreme, ScanTime, Capacity,
  #                      TotalWait, NumPatients, WaitExtremes
  
  # We'll produce a summary for each facility
  facilities <- unique(sim_data$Facility)
  
  results <- list()
  
  for (f in facilities) {
    df_f <- subset(sim_data, Facility == f)
    n_days <- length(unique(df_f$Day))  # how many days in the simulation
    
    # OvertimeCount = how many days had overtime>0 => sum of daily OvertimeCount
    OvertimeCount <- sum(df_f$OvertimeCount)
    
    # OvertimeAverage = average overtime hours per day => sum(OvertimeHours)/n_days
    OvertimeAverage <- sum(df_f$OvertimeHours) / n_days
    
    # OvertimeCountExtremes = how many days finishing time>19 => sum daily OvertimeCountExtreme
    OvertimeCountExtremes <- sum(df_f$OvertimeCountExtreme)
    
    # Utilization => total ScanTime / (Capacity*n_days)
    total_scan <- sum(df_f$ScanTime)
    total_capacity <- sum(df_f$Capacity)*1  # = 9 * n_days if consistent
    Utilization <- total_scan / (9 * n_days)  # or total_scan / total_capacity
    
    # WaitingtimeAverage => total waiting / total patients
    total_wait   <- sum(df_f$TotalWait)
    total_pats   <- sum(df_f$NumPatients)
    if (total_pats>0) {
      WaitingtimeAverage <- total_wait / total_pats
    } else {
      WaitingtimeAverage <- 0
    }
    
    # WaitingtimeExtremes => sum of daily WaitExtremes ( # of patients waiting>20m )
    WaitingtimeExtremes <- sum(df_f$WaitExtremes)
    
    results[[paste0("Facility",f)]] <- data.frame(
      Facility = f,
      OvertimeCount         = OvertimeCount,
      OvertimeAverage       = OvertimeAverage,
      OvertimeCountExtremes = OvertimeCountExtremes,
      Utilization           = Utilization,
      WaitingtimeAverage    = WaitingtimeAverage,
      WaitingtimeExtremes   = WaitingtimeExtremes
    )
  }
  # Combine rows for each facility
  final_df <- do.call(rbind, results)
  rownames(final_df) <- NULL
  final_df
}

#--------------------------------------------------------------------------------
# Example single run usage
results_old <- run_simulation(num_days = 23, policy = "Old", seed = 777)
results_new <- run_simulation(num_days = 23, policy = "New", seed = 777)

summary_old <- summarize_simulation(results_old)
summary_new <- summarize_simulation(results_new)

cat("\n--- OLD POLICY (Single Run) ---\n")
print(summary_old)
cat("\n--- NEW POLICY (Single Run) ---\n")
print(summary_new)

#--------------------------------------------------------------------------------
# MULTIPLE REPLICATIONS
run_replications <- function(num_days, policy, R) {
  # We'll run R times, each time storing the summary for each facility
  # Then we can average the results across replications if we like
  rep_list <- list()
  
  for (r in seq_len(R)) {
    sim_data <- run_simulation(num_days, policy)
    sum_data <- summarize_simulation(sim_data)
    # We'll add a column for replication index
    sum_data$Replication <- r
    rep_list[[r]] <- sum_data
  }
  # Combine all replication-level summaries
  big_df <- do.call(rbind, rep_list)
  rownames(big_df) <- NULL
  big_df
}

# Example: 4 replications
rep_old <- run_replications(23, "Old", R = 12)
rep_new <- run_replications(23, "New", R = 12)

cat("\n--- OLD POLICY (Multiple Replications) ---\n")
print(rep_old)
cat("\n--- NEW POLICY (Multiple Replications) ---\n")
print(rep_new)

# If you want to see the average over all replications for each facility:
aggregate(. ~ Facility, data=rep_old[ , -ncol(rep_old)], FUN=mean)
aggregate(. ~ Facility, data=rep_new[ , -ncol(rep_new)], FUN=mean)
# (That excludes Replication col, or use a more advanced aggregator)

#--------------------------------------------------------------------------------
# Inspect Single Day
set.seed(777)
policy = "Old"
inspect_one_day <- function(day, policy) {
  sched <- simulate_one_day(day, policy)
  cat("\n--- Sim Stats for Day", day, "Policy:", policy, "---\n")
  print(sched)
}
inspect_one_day(day = 1, policy)

#--------------------------------------------------------------------------------
# Print Daily Schedules
set.seed(777)
inspect_all_days <- function(num_days = 23, policy = c("Old","New")) {
  for (d in seq_len(num_days)) {
    cat("\n--- Day", d, "---\n")
    day_stats <- simulate_one_day(d, policy)
    # day_stats is facility-level info, not the patient schedule
    print(day_stats)
  }
}
# Example usage
inspect_all_days(23, "Old")
