
library(ggplot2)
library(lubridate)
library(tidyverse)
library(zoo)

#Datasets
ScanRecords <- read.csv("ScanRecords.csv")
ScanRecords$Date <- as.Date(ScanRecords$Date)
Type1 <- filter(ScanRecords, PatientType=="Type 1")
Type2 <- filter(ScanRecords, PatientType=="Type 2")

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

  
