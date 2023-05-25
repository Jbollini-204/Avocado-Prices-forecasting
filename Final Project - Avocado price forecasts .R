### Final Assignment 
### Time Series and Forecasting Analysis
### Josefina Bollini
library(forecast)
library(xts)
library(tseries)
library(lubridate)
library(TSA)
library(dplyr)
library(ggplot2)
library(prophet)
library(tibble)

# Read the data - weekly avocado prices in Chicago
Avocado_prices = read.csv('/Users/Jbollini/Desktop/Time Series/Final Project/avocado-updated-2020.csv')
Avocado_prices <- Avocado_prices[,c(1,2,11,13)]
Avocado_prices <- Avocado_prices[Avocado_prices$geography  == "Chicago", ] 

# Make sure I ahve one observation by week by type
Avocado_prices  <- aggregate(Avocado_prices$average_price,
                              by=list(date = Avocado_prices$date,article = Avocado_prices$type), 
                              FUN = mean)

# EDA
Average_prices_by_year <- aggregate(Avocado_prices$x,
                            by=list(Year = year(Avocado_prices$date),article = Avocado_prices$article), 
                            FUN = function(x) c(mean = mean(x), median = median(x), sd = sd(x)))


# there are 2 types of avocados (Organic and conventional)
# check if the series are uniform or not - 
## The store might be closed some days of the week or special holidays

uniform_ts <- function(data.in,type) {
  data.in <- data.in[data.in$article == type, ]
  data.in$date <- as.Date(data.in$date, format = "%Y-%m-%d")
  date_range <- seq(as.Date("2015-01-04"), as.Date("2020-11-29"), by = "week")
  new_xts <- xts(rep(NA, length(date_range)), order.by = date_range)
  
  x_a_avgp <- aggregate(data.in$x, by = list(date=data.in$date), FUN = sum)
  x_a_avgp_xts <- xts(x_a_avgp$x, order.by = data.in$date)
  
  alldates_x_a_avgp_xts <- merge(x_a_avgp_xts, new_xts, all = TRUE)
  alldates_x_a_avgp_xts <-  alldates_x_a_avgp_xts[ , 1]
  # Linear interpolaiton
  ts_imputed_li <- na.approx(alldates_x_a_avgp_xts)
  print(plot(ts_imputed_li))
  
  # Spline 
  #ts_inputed_spline <- na.spline(alldates_x_a_avgp_xts)
  colnames(ts_imputed_li)[colnames(ts_imputed_li) == "x_a_sum_xts"] <- type
  return(ts_imputed_li)
}

par(mfrow=c(2,1))
organic_xts = uniform_ts(Avocado_prices,'organic')
conventional_xts = uniform_ts(Avocado_prices,'conventional')

### Split in train and test 
### Forecast for 2020
split_train_test <- function(data.in){
  train = data.in["2015-01-04/2019-12-29"]
  test = data.in["2019-12-29/2020-11-29"]
  return(list(train,test))
}

organic_xts.split = split_train_test(organic_xts)
conventional_xts.split = split_train_test(conventional_xts)

# STATIONARITY TESTS
stationarity_checks <- function(data) {
  # ADF test
  print('ADF test - H0 : non stationarity')
  adf_result <- adf.test(data)
  print(adf_result)
  
  print('KPSS test - H0 : stationarity')
  kpss_result <- kpss.test(data)
  print(kpss_result)
  
  # Automatic differencing order determination
  diff_order <- ndiffs(data)
  print(diff_order)
}

stationarity_checks(organic_xts)
stationarity_checks(conventional_xts)

explore_tests <- function(data.in,avo_type) {
  data.in %>% mstl() %>%
    autoplot()
  p1 <- periodogram(data.in)
  acf(data.in, main=avo_type)
  pacf(data.in, main=avo_type)
  
  Box.test(data.in, lag = 1, type = "Ljung-Box")
}
par(mfrow=c(2,3))
explore_tests(organic_xts,'Organic') 
explore_tests(conventional_xts,'Conventional') 


decomposition_ts <- function(data.in,type_avo) {
  data_ts <- as.ts(data.in)
  p1 <- autoplot(data_ts) +
  ylab(type_avo) + xlab("Weeks") 

p2 <- autoplot(window(data_ts, end=5)) +
  ylab("Baguette sales") + xlab("Years") 
gridExtra::grid.arrange(p1,p2)

data_ts %>% mstl() %>%
  autoplot() + xlab("Week")
}

decomposition_ts(organic_xts.split[[1]],'organic')
decomposition_ts(conventional_xts.split[[1]],'conventional')
performance_metrics_all <- function(actual_values,predicted_values){
  MAPE <- mean(abs((actual_values - predicted_values) / actual_values)) * 100
  SMAPE <- mean(2 * abs(actual_values - predicted_values) / (abs(actual_values) + abs(predicted_values))) * 100
  RMSE <- sqrt(mean((actual_values - predicted_values)^2))
  
  # Print the calculated metrics
  print(paste("MAPE:", MAPE))
  print(paste("SMAPE:", SMAPE))
  print(paste("RMSE:", RMSE))
  
}

Test_models <- function(data.in,data.in.test,type_art) {
  ts_data <- ts(data.in, frequency = 52) 
  actuals <- ts(data.in.test, frequency = 52) 
  # auto Arima
  print('Auto Arima Results')
  arima_model <- auto.arima(ts_data,seasonal = TRUE,nmodels = 500)
  print(summary(arima_model))
  forecast_arima <-forecast(arima_model,h=50)
  plot(forecast_arima)
  forecast::checkresiduals(arima_model)

  performance_metrics_all(c(actuals[,1]),c(forecast_arima$mean))
  
  # Holt Winters
  print('Holt Winters Results')
  hw_model <- HoltWinters(ts_data)    
  print(summary(hw_model))
  forecast_hw <- forecast(hw_model, h = 50)  # Forecast 12 periods ahead
  plot(forecast_hw) 
  performance_metrics_all(c(actuals[,1]),c(forecast_hw$mean))
  
  #Arfima
  print('Arfima Results')
  arfima_model <- forecast::arfima(ts_data)
  print(summary(arfima_model))
  forecast::checkresiduals(arfima_model)
  forecast_arfima<-forecast(arfima_model,h=50)
  plot(forecast_arfima)
  performance_metrics_all(c(actuals[,1]),c(forecast_arfima$mean))

}  

# Prophet Model

transform_to_df <-  function(data.input){
  df <- data.frame(ds = index(data.input), y = coredata(data.input))
  colnames(df) <- c("ds", "y")
  return(df)
}

prophet_model <- function(data.in,data.in.test,type){
  prophet_model <- prophet(data.in)
  
  # Generate future dates for forecasting
  future_dates <- make_future_dataframe(prophet_model, periods = 50, freq = "week") 
  # Perform the forecast
  forecast_data <- predict(prophet_model, future_dates)
  print(plot(prophet_model, forecast_data))
  prophet_plot_components(prophet_model, forecast_data)

  predicted_values <- forecast_data$yhat
  actual_values <- data.in.test$y
  
  performance_metrics_all(actual_values,predicted_values)
}
par(mfrow=c(1,1))
Test_models(organic_xts.split[[1]],organic_xts.split[[2]],'Organic')  
prophet_model(transform_to_df(organic_xts.split[[1]]),
              transform_to_df(organic_xts.split[[2]]),
              'organic')
Test_models(conventional_xts.split[[1]],organic_xts.split[[2]],'Conventional')
prophet_model(transform_to_df(conventional_xts.split[[1]]),
              transform_to_df(conventional_xts.split[[2]]),
              'Conventional')



