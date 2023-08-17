library(dplyr)             # Load the dplyr package
library(tidyr)
library(purrr)
library(data.table)
library(geosphere) #nearest Neighbourt
library(stringdist) #cluster 
library(lubridate)
library(mgcv)

# List all files in the folder
AZ_files <- list.files(path = "Rain_Data/Arizona", pattern = "\\.csv$", full.names = TRUE)
NM_files <- list.files(path = "Rain_Data/NewMexico", pattern = "\\.csv$", full.names = TRUE)

# Read and bind all CSV files
AZ_data <- AZ_files %>%
  map_dfr(read.csv)

NM_data <- NM_files %>%
  map_dfr(read.csv)

#combined and format date
SW_raindata <- rbind(AZ_data, NM_data)
SW_raindata$DATE <- as.Date(SW_raindata$DATE)

#REad ANPP Data
ANPP <- read.csv("data/ANPP_2020-2022.csv")
ANPP$Forecast <- as.Date(ANPP$Forecast)
ANPP <- split(ANPP, ANPP$Forecast)

Data <- function(N) {

#Set Dates
start_date <- as.Date(names(ANPP)[N]) - days(6)
end_date <- as.Date(names(ANPP)[N])

# Filter only rows within the specific date range
ANPP_Spring <- do.call(rbind, ANPP) %>%
  filter(Forecast >= start_date, Forecast <= end_date)

# Filter only rows within the specific date range
SW_rain <- SW_raindata %>%
  filter(DATE >= start_date, DATE <= end_date) %>%
  select(1:6, PRCP)

sort(unique(SW_rain$DATE), decreasing = TRUE)

# Create biweekly intervals using the cut function
SW_rain <- as.data.table(SW_rain)
SW_rain[, biweek_end_date := cut(DATE, breaks = "2 weeks")]

# Summarize precipitation values for each biweek using data.table
SW_rain_biweek <- SW_rain[, .(total_PRCP = sum(PRCP)),
                                by = .(STATION, biweek_end_date, NAME, LATITUDE, LONGITUDE, ELEVATION)] %>%
                                 mutate(week_number = week(biweek_end_date))



##### NEAREST NEIGHBOUR ####

# Function to find the nearest neighbor in df2 for each row in df1
find_nearest_neighbor <- function(lat, long, df) {
  distances <- distGeo(matrix(c(long, lat), ncol = 2), matrix(c(df$long, df$lat), ncol = 2))
  nearest_index <- which.min(distances)
  return(df[nearest_index, ])
}

# Create an empty list to store the results
result_list <- list()

# Iterate over each row in df1 and find the nearest neighbor in df2
for (i in 1:nrow(SW_rain_biweek)) {
nearest_neighbor <- find_nearest_neighbor(SW_rain_biweek$LATITUDE[i], SW_rain_biweek$LONGITUDE[i], ANPP_Spring)
  result_list[[i]] <- cbind(SW_rain_biweek[i, ], nearest_neighbor[3:4],NPP_predict_avg = nearest_neighbor$NPP_predict_avg)
}

# Combine the results into a new data frame
merged_df <- do.call(rbind, result_list)

# Split names into two columns using the comma separator
split_names <- strsplit(merged_df$NAME, ",")

# Create new columns for First_Name and Last_Name
merged_df$County <- sapply(split_names, "[", 1)
merged_df$State <- sapply(split_names, "[", 2)

# Calculate the string distance matrix
dist_matrix <- stringdist::stringdistmatrix(merged_df$County, merged_df$County)

# Define a threshold for similarity
similarity_threshold <- 6 # You can adjust this value based on your data

# Perform hierarchical clustering
hc <- hclust(as.dist(dist_matrix))
clusters <- cutree(hc, h = similarity_threshold)

# Add the cluster information to the data frame
merged_df$cluster <- clusters

# Group by cluster and summarize the data
Cluster_summary <- merged_df %>%
  group_by(cluster) %>%
  summarise(
    Forecast = Forecast,
    mean_PRCP = mean(total_PRCP, na.rm = TRUE),
    avg_NPP_predict_avg = mean(NPP_predict_avg, na.rm = TRUE),
    LATITUDE = names(which.max(table(LATITUDE))),
    LONGITUDE = names(which.max(table(LONGITUDE))),
    County = names(which.max(table(County))),
    State = names(which.max(table(State))),
    week_number = names(which.max(table(week_number)))
  )

return(Cluster_summary)

}

# Initialize an empty data frame to store the results
Prep_Anpp <- data.frame()

# Loop to generate data frames and append to the combined data frame
for (i in 1:length(ANPP)) {
  result_df <- Data(i)
  Prep_Anpp <- rbind(Prep_Anpp, result_df)
}

Prep_Anpp <- Prep_Anpp %>%
  mutate_at(vars(mean_PRCP, avg_NPP_predict_avg),
            ~ ifelse(is.nan(.), 0, .)) %>%
  rename(y = LATITUDE, x = LONGITUDE)

# Convert y and x to numeric 
Prep_Anpp$y <- as.numeric(Prep_Anpp$y)
Prep_Anpp$x <- as.numeric(Prep_Anpp$x)

#Add the rain ratio with nearest neighbour

for (i in 1:nrow(Prep_Anpp)) {
  nearest_neighbor <- find_nearest_neighbor(Prep_Anpp$y[i], Prep_Anpp$x[i], ANPP)
    ppt_ratio <- nearest_neighbor$pptRatioSummerWinter
  Prep_Anpp$pptRatioSummerWinter[i] <- ppt_ratio
}

#############################RUN THE MODELS ###################

#### All Data
GAM_All <- gam(avg_NPP_predict_avg ~ s(x, y) + s(mean_PRCP), data = Prep_Anpp, method = "REML")
GAM_AllPlot <-plot(GAM, scheme = 2, pages = 1)

summary(GAM_All)
### Year 2020
Prep_Anpp_2020 <- Prep_Anpp %>%
  filter(year(Forecast) == 2020) 

GAM2020 <- gam(avg_NPP_predict_avg ~ s(x, y) + s(mean_PRCP), data = Prep_Anpp_2020, method = "REML")
GAM2020_plot <-plot(GAM2020, scheme = 2, pages = 1)

summary(GAM2020)


### Year 2021
Prep_Anpp_2021 <- Prep_Anpp %>%
  filter(year(Forecast) == 2021)

GAM2021 <- gam(avg_NPP_predict_avg ~ s(x, y) + s(mean_PRCP), data = Prep_Anpp_2021, method = "REML")
GAM2021_plot <-plot(GAM2021, scheme = 2, pages = 1)

summary(GAM2021)


### Year 2022
Prep_Anpp_2022 <- Prep_Anpp %>%
  filter(year(Forecast) == 2022) 

GAM2022 <- gam(avg_NPP_predict_avg ~ s(x, y) + s(mean_PRCP), data = Prep_Anpp_2022, method = "REML")
GAM2022_plot <- plot(GAM2022, scheme = 2, pages = 1)

summary(GAM2022)

###Spliting by rain zones ####

### Monsson
Prep_Anpp_Mon <- Prep_Anpp %>%
  filter(pptRatioSummerWinter >= 1.2 & pptRatioSummerWinter <= 3)

GAM_mon <- gam(avg_NPP_predict_avg ~ s(x, y) + s(mean_PRCP), data = Prep_Anpp_Mon, method = "REML")
GAM_monPlot <- plot(GAM_mon, scheme = 2, pages = 2)

summary(GAM_mon)

### Winter
Prep_Anpp_Win <- Prep_Anpp %>%
  filter(pptRatioSummerWinter >= 0.3 & pptRatioSummerWinter <= 0.8)

GAM_win <- gam(avg_NPP_predict_avg ~ s(x, y) + s(mean_PRCP), data = Prep_Anpp_Win, method = "REML")
GAM_winPlot <- plot(GAM_win, scheme = 2, pages = 1)

summary(GAM_win)

### Transition
Prep_Anpp_Tran <- Prep_Anpp %>%
  filter(pptRatioSummerWinter >= 0.8 & pptRatioSummerWinter <= 1.2)

GAM_Tran <- gam(avg_NPP_predict_avg ~ s(x, y) + s(mean_PRCP), data = Prep_Anpp_Tran, method = "REML")
GAM_TranPlot <- plot(GAM_Tran, scheme = 2, pages = 1)

summary(GAM_Tran)
#