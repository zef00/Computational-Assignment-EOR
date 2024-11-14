
library(ggplot2)
library(lubridate)
library(tidyverse)
library(zoo)

#Datasets
ScanRecords <- read.csv("ScanRecords.csv")
ScanRecords$Date <- as.Date(ScanRecords$Date)
Type1 <- dplyr::filter(ScanRecords, PatientType=="Type 1") #Ik moest dplyr:: toevoegen want hij had een conflict met een ander package
Type2 <- dplyr::filter(ScanRecords, PatientType=="Type 2") #Ik moest dplyr:: toevoegen want hij had een conflict met een ander package

#Duration Scatterplot 
ggplot(ScanRecords, aes(x=1:618, y=Duration, color = PatientType))+
  geom_point() +
  geom_hline(yintercept = mean(Type1$Duration), color = "red") +
  geom_hline(yintercept = mean(Type2$Duration), color = "blue")

#Duration Normal Distribution PDF plot
pdf1 <- dnorm(Type1$Duration, mean=mean(Type1$Duration), sd = sd(Type1$Duration))
plot(Type1$Duration, pdf1)
pdf2 <- dnorm(Type2$Duration, mean=mean(Type2$Duration), sd = sd(Type2$Duration))
plot(Type2$Duration, pdf2)

#Arrival time/Poisson Process 
# For now this only gives the PMF and lambda for both patienttypes together, we
# should still use bootstrap to estimate the lambda and do it also per patientype. 

# Step 1: Create a timestamp by combining Date and Time, then round to hourly intervals
ScanRecords <- ScanRecords %>%
  mutate(timestamp = as.POSIXct(Date) + dhours(Time),   # Combine Date and Time
         hour = floor_date(timestamp, "hour"))          # Round to hourly intervals

# Step 2: Count the number of arrivals per hour per day
hourly_counts <- ScanRecords %>%
  group_by(Date, PatientType, hour) %>%
  summarise(count = n()) %>%
  ungroup()

# Step 3: Calculate the mean arrival rate (lambda) per hour
lambda <- mean(hourly_counts$count)

# Step 4: Create the observed PMF by calculating the frequency of each count value
observed_pmf <- hourly_counts %>%
  count(count) %>%
  mutate(prob = n / sum(n))

# Step 5: Generate the theoretical Poisson PMF based on the mean arrival rate
# Generate the counts we observed (unique arrival values)
max_count <- max(hourly_counts$count)
theoretical_pmf <- data.frame(
  count = 0:max_count,
  prob = dpois(0:max_count, lambda)
)

# Step 6: Plot the observed PMF and overlay the theoretical Poisson PMF
ggplot() +
  geom_bar(data = observed_pmf, aes(x = count, y = prob), stat = "identity", fill = "skyblue", alpha = 0.6, width = 0.8) +
  geom_line(data = theoretical_pmf, aes(x = count, y = prob), color = "red", size = 1.2) +
  geom_point(data = theoretical_pmf, aes(x = count, y = prob), color = "red", size = 2) +
  labs(x = "Number of Arrivals per Hour", y = "Probability", title = "Probability Mass Function of Arrival Times per Hour") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = 0:max_count)


# ARRIVAL TIMES TYPE 1 patient
# Count daily arrivals for Type 1
daily_arrivals_type1 <- Type1 %>%
  group_by(Date) %>%
  summarise(daily_count = n())

# Calculate mean arrival rate (lambda) for Poisson distribution
lambda_type1 <- mean(daily_arrivals_type1$daily_count)

# Parametric Bootstrap to estimate uncertainty for lambda
bootstrap_lambda <- replicate(1000, {
  sample_counts <- rpois(n = nrow(daily_arrivals_type1), lambda = lambda_type1)
  mean(sample_counts)
})

# Create a data frame for the bootstrapped lambdas
bootstrap_lambda_df <- data.frame(lambda = bootstrap_lambda)

# Calculate the mean of the bootstrap lambdas
mean_lambda <- mean(bootstrap_lambda)

# Calculate 95% confidence interval for bootstrap lambda values
lambda_ci <- quantile(bootstrap_lambda_df$lambda, probs = c(0.025, 0.975))

# Plot the distribution of bootstrap lambda values with mean and confidence interval lines
ggplot(bootstrap_lambda_df, aes(x = lambda)) +
  geom_density(fill = "lightblue", alpha = 0.5) +             # Density plot of bootstrap lambdas
  geom_vline(xintercept = mean_lambda, color = "red",         # Line at mean lambda
             linetype = "dashed", size = 1) +
  geom_vline(xintercept = lambda_ci[1], color = "blue",       # Lower bound of CI
             linetype = "dotted", size = 1) +
  geom_vline(xintercept = lambda_ci[2], color = "blue",       # Upper bound of CI
             linetype = "dotted", size = 1) +
  labs(title = "Distribution of Bootstrapped Lambda Estimates for Type 1 Daily Arrivals",
       subtitle = paste("95% CI:", round(lambda_ci[1], 2), "-", round(lambda_ci[2], 2)),
       x = "Lambda (Average Daily Arrivals)", y = "Density") +
  theme_minimal()

#dit is testtestest
#dfsdklfjsd
#ftyhhgfngfhb
