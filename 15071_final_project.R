library(dplyr)
library(jsonlite)
library(rvest)
library(randomForest)
library(caret)
library(nnet)
library(zoo)
library(data.table)

library(quanteda)
library(tm)


# DATA PREPROCESSING ------------------------------------------------------
preprocess_tweets <- function(tweets) {
  # Convert columns to proper types
  tweets$user_created <-
    strptime(tweets$user_created, "%Y-%m-%d %H:%M:%S")
  tweets$user_verified <- as.logical(tweets$user_verified)
  tweets$date <-
    strptime(tweets$date, "%Y-%m-%d %H:%M:%S")
  tweets$hashtags <-
    lapply(strsplit(tweets$hashtags, "[][']|,\\s*"), function(x)
      x[nzchar(x)])
  tweets$source <- as.factor(tweets$source)
  tweets$is_retweet <- as.logical(tweets$is_retweet)

  # Match user location string to coordinates
  geocode <- function(location) {
    coordinates <- data.frame()
    
    src_url <- "https://nominatim.openstreetmap.org/search?q="
    requests <-
      paste0(src_url, gsub("#|&", "", location), "&format=geojson")
    
    for (i in 1:length(requests)) {
      response <- tryCatch({
        res <- read_html(requests[i]) %>%
          html_node("p") %>%
          html_text() %>%
          fromJSON()
        if (length(res$features$geometry) == 0)
          NA
        else
          res$features$geometry$coordinates[[1]]
      },
      error = function (c)
        NA,
      warning = function (c)
        NA)
      
      coordinates <- rbind(coordinates, response)
    }
    
    names(coordinates) <- c("user_longitude", "user_latitude")
    return(coordinates)
  }
  coordinates <- geocode(tweets$user_location) # takes a long time!
  tweets <- cbind(tweets, coordinates)
  
  # Save tweets dataframe to file
  saveRDS(tweets, file = "tweets.rds")
  return(tweets)
}
# tweets <- read.csv("vaccination_all_tweets.csv")
# tweets <- preprocess_tweets(tweets)
tweets <- readRDS("tweets.rds") # load preprocessed data to save time

preprocess_stock <-
  function(stock,
           volatility,
           name = "stock",
           lags = 7,
           lag_vars = c('Close.Last', 'Volume', 'Volatility', 'Spread')) {
    # Convert columns to proper types
  stock$Date <- as.Date(strptime(stock$Date, "%m/%d/%Y"))
  stock$Close.Last <- as.numeric(gsub("\\$", "", stock$Close.Last))
  stock$Open <- as.numeric(gsub("\\$", "", stock$Open))
  stock$High <- as.numeric(gsub("\\$", "", stock$High))
  stock$Low <- as.numeric(gsub("\\$", "", stock$Low))
  
  stock <- stock[order(stock$Date),]
  
  volatility$Date <- as.Date(strptime(volatility$Date, "%m/%d/%y"))
  stock <- merge(stock, volatility, by.x = "Date", by.y = "Date")
  colnames(stock)[ncol(stock)] <- 'Volatility'
  
  # Add other potential outcome variables
  stock$Spread <- stock$High - stock$Low
  stock$Return <-
    stock$Close.Last / shift(stock$Close.Last, n = 1, type = 'lag')
  
  # Add lagged variables
  for (k in 1:lags) {
    for (var in lag_vars) {
      Lag <- shift(stock[,var], n = k, type = 'lag')
      stock <- cbind(stock, Lag)
      colnames(stock)[ncol(stock)] <-
        c(paste(var, k, sep = "."))
    }
  }
  
  # Save stock dataframe to file
  saveRDS(stock, file = paste(name, ".rds", sep = ""))
  return(stock)
}
# volatility <- read.csv("historical-volatility.csv")
# pfizer <- read.csv("PfizerHistoricalData.csv")
# pfizer <-
#   preprocess_stock(pfizer, volatility[, c('Date', 'Pfizer')], name = "pfizer")
# moderna <- read.csv("ModernaHistoricalData.csv")
# moderna <-
#   preprocess_stock(moderna, volatility[, c('Date', 'Moderna')], name = "moderna")
pfizer <- readRDS("pfizer.rds") # load preprocessed data to save time
moderna <- readRDS("moderna.rds") # load preprocessed data to save time

merge_data <- function(tweets, stock, hashtag, verified = T) {
  # Filter tweets by hashtag
  company_tweets <-
    dplyr::filter(tweets, sapply(hashtags, function(list)
      hashtag %in% list))
  # Combine tweet text and stats by day
  if (verified) {
    company_tweets <- subset(company_tweets, user_verified)
  }
  daily_tweets <-
    company_tweets %>% group_by(Date = as.Date(date)) %>% summarize(
      Text = paste(text, collapse = " "),
      User.Followers = sum(user_followers),
      User.Friends = sum(user_friends),
      User.Favourites = sum(user_favourites),
      Retweets = sum(retweets),
      Favorites = sum(favorites),
      Tweets = n()
    )
  daily_tweets$Nchar <- nchar(daily_tweets$Text)
  
  # Merge tweet and stock data
  merged_data <- merge(daily_tweets, stock, by.x = "Date", by.y = "Date", all.y = T)
  merged_data <- subset(merged_data, Date >= as.Date('2020-12-14'))
  return(merged_data)
}
pfizer <- merge_data(tweets, pfizer, "PfizerBioNTech")
moderna <- merge_data(tweets, moderna, "Moderna")


# VISUALIZATIONS ----------------------------------------------------------
library(reshape2)
library(lattice)
# Appendix 1
pfizer.mm <-
  melt(subset(
    pfizer,
    select = c('Date', 'Close.Last', 'Volume', 'Volatility', 'Spread')
  ), id.var = "Date")
xyplot(
  value ~ Date | variable,
  data = pfizer.mm,
  type = "l",
  scales = list(y = list(relation = "free", rot = 0)),
  ylab = "Value",
  main = "Pfizer Stock Data",
  layout = c(1, 4)
)
moderna.mm <-
  melt(subset(
    moderna,
    select = c('Date', 'Close.Last', 'Volume', 'Volatility', 'Spread')
  ), id.var = "Date")
xyplot(
  value ~ Date | variable,
  data = moderna.mm,
  type = "l",
  scales = list(y = list(relation = "free", rot=0)),
  ylab = "Value",
  main = "Moderna Stock Data",
  layout = c(1, 4)
)


# RUN MODELS --------------------------------------------------------------
data <- pfizer # set data to use
split <- 0.7
train <- data[1:as.integer(nrow(data) * split), ]
test <- data[(as.integer(nrow(data) * split) + 1):nrow(data), ]

# Function to generate models for a given formula
run_models <- function(outcome, formula, train_set = train) {
  results <- c()
  # Linear regression model
  lm.cv <- train(
    formula,
    data = train_set,
    method = "lm",
    trControl = trainControl(method = "cv", number = 10),
    preProcess = "medianImpute",
    na.action = na.pass
  )
  lm.model <- lm.cv$finalModel
  
  # Random forest model
  rf.cv = train(
    formula,
    data = train_set,
    method = "rf",
    trControl = trainControl(method = "cv", number = 10),
    tuneGrid = data.frame(mtry = seq(1, 10, 1)),
    preProcess = "medianImpute",
    na.action = na.pass
  )
  rf.model <- rf.cv$finalModel
  
  # CART model
  tree.cv = train(
    formula,
    data = train_set,
    method = "rpart",
    trControl = trainControl(method = "cv", number = 10),
    tuneGrid = expand.grid(cp = seq(0.001, 0.01, .0005)),
    preProcess = "medianImpute",
    na.action = na.pass
  )
  tree.model <- tree.cv$finalModel
  
  return(list("lm" = lm.model, "rf" = rf.model, "cart" = tree.model))
}
# Function to plot comparison of model predictions and actual data
plot_comparison <- function(models, outcome, test_set = test, title = "") {
  test_set <- na.exclude(test_set)
  pred.lm <- predict(models$lm, test_set)
  pred.rf <- predict(models$rf, test_set)
  pred.cart <- predict(models$cart, test_set)

  ylim <-
    c(
      min(test_set[, outcome], pred.lm, pred.rf, pred.cart),
      max(test_set[, outcome], pred.lm, pred.rf, pred.cart)
    )
  plot(
    test_set$Date,
    test_set[, outcome],
    type = 'l',
    xlab = "Date",
    ylab = outcome,
    ylim = ylim,
    main = title
  )
  lines(test_set$Date, pred.lm, col = 'blue')
  lines(test_set$Date, pred.rf, col = 'red')
  lines(test_set$Date, pred.cart, col = 'green')
  legend(
    "topright",
    legend = c("Actual", "Linear Regression", "Random Forest", "CART"),
    col = c("black", "blue", "red", "green"),
    lty = 1,
    cex = 0.8
  )
}
# Function to calculate evaluation metrics on a given model
calculate_metrics <-
  function(model,
           train_set,
           test_set,
           outcome = 'Close.Last') {
    pred.train <- predict(model, train_set)
    r2 <-
      1 - sum((pred.train - train_set[, outcome]) ^ 2) / sum((mean(train_set[, outcome]) - train_set[, outcome]) ^ 2)
    aic <- nrow(train_set) * log(sum((pred.train - train_set[, outcome]) ^ 2) / nrow(train_set))
    pred.test <- predict(model, test_set)
    osr2 <-
      1 - sum((pred.test - test_set[, outcome]) ^ 2) / sum((mean(train_set[, outcome]) - test_set[, outcome]) ^ 2)
    metrics <- 
      c("R2" = r2, "OSR2" = osr2, "AIC" = aic)
    if (outcome == 'Close.Last') {
      return <- get_return(pred.test, test_set)
      metrics <- c(metrics, "Return" = return)
    }
    return(metrics)
  }
get_return <- function(pred, test_set) {
  pred <- c(diff(pred), NA)
  buy <- c(as.numeric(pred > 0))
  dailyreturn <- (test_set$Return - 1) * buy
  dailyreturn <- dailyreturn[1:(length(dailyreturn) - 1)]
  cash <- 1000000
  for (i in 1:length(dailyreturn)) {
    cash <- cash * (as.numeric((dailyreturn[i] + 1)))
  }
  return(cash / 1000000 * 100 - 100)
}
# Function to summarize evaluation metrics on all models
metrics_comparison <- function(models, outcome, train_set = train, test_set = test) {
  train_set <- na.exclude(train_set)
  test_set <- na.exclude(test_set)
  
  results <- rbind(
    "lm" = calculate_metrics(models$lm,
                             train_set,
                             test_set,
                             outcome = outcome),
    "rf" = calculate_metrics(models$rf,
                             train_set,
                             test_set,
                             outcome = outcome),
    "cart" = calculate_metrics(models$cart,
                               train_set,
                               test_set,
                               outcome = outcome)
  )
  
  return(results)
}

# Find optimal number of lags
trainSubset <- train[, c(10, 17, 3:9)]
testSubset <- test[, c(10, 17, 3:9)]
R2 <- c()
OSR2 <- c()
AIC <- c()
Return <- c()
for (k in 0:7) {
  if (k > 0) {
    trainSubset <- train[, c(10, 17, 3:9, 18:(18 + 4 * k - 1))]
    testSubset <- test[, c(10, 17, 3:9, 18:(18 + 4 * k - 1))]
  }
  models <- run_models('Close.Last', Close.Last ~ ., train_set = na.exclude(trainSubset))
  results <- metrics_comparison(models, 'Close.Last', train_set = na.exclude(trainSubset), test_set = na.exclude(testSubset))
  R2 <- cbind(R2, results[, 'R2'])
  OSR2 <- cbind(OSR2, results[, 'OSR2'])
  AIC <- cbind(AIC, results[, 'AIC'] + 2 * (k + 7 + 1))
  Return <- cbind(Return, results[, 'Return'])
}
colnames(R2) <- 0:7
colnames(OSR2) <- 0:7
colnames(AIC) <- 0:7
plot_statistic <- function(statistic, ylab, title, pos) {
  plot(0:7,
       statistic[1,],
       type = 'l',
       col = 'blue',
       ylim = c(min(statistic), max(statistic)),
       xlab = "Number of Lags",
       ylab = ylab,
       main = title)
  lines(0:7, statistic[2, ], col = 'red')
  lines(0:7, statistic[3, ], col = 'green')
  legend(
    4.5, pos,
    legend = c("Linear Regression", "Random Forest", "CART"),
    col = c("blue", "red", "green"),
    lty = 1,
    cex = 0.9
  )
}
# Appendix 2
plot_statistic(R2, ylab = "R^2", title = "R^2 for Different Autoregression Models", pos = 0.8)
plot_statistic(OSR2, ylab = "OSR^2", title = "OSR^2 for Different Autoregression Models", pos = 0.2)
plot_statistic(AIC, ylab = "AIC", title = "AIC for Different Autoregression Models", pos = -200)
plot_statistic(Return, ylab = "Return", title = "Return for Different Autoregression Models", pos = 4)

# Try models with just Tweet data
close.last.tweets.models <- run_models(
  'Close.Last',
  Close.Last ~ User.Followers + User.Friends + User.Favourites + Retweets + Favorites + Tweets + Nchar
)
close.last.tweets <- metrics_comparison(close.last.tweets.models, 'Close.Last')
volatility.tweets.models <- run_models(
  'Volatility',
  Volatility ~ User.Followers + User.Friends + User.Favourites + Retweets + Favorites + Tweets + Nchar
)
volume.tweets.models <- run_models(
  'Volume',
  Volume ~ User.Followers + User.Friends + User.Favourites + Retweets + Favorites + Tweets + Nchar
)
spread.tweets.models <- run_models(
  'Spread',
  Spread ~ User.Followers + User.Friends + User.Favourites + Retweets + Favorites + Tweets + Nchar
)

# Try models with just lagged stock data
close.last.stock.models <- run_models(
  'Close.Last',
  Close.Last ~ Close.Last.1 + Volume.1 + Volatility.1 + Spread.1 + Close.Last.2 + Volume.2 + Volatility.2 + Spread.2
)
close.last.stock <- metrics_comparison(close.last.stock.models, 'Close.Last')
volatility.stock.models <- run_models(
  'Volatility',
  Volatility ~ Close.Last.1 + Volume.1 + Volatility.1 + Spread.1 + Close.Last.2 + Volume.2 + Volatility.2 + Spread.2
)
volume.stock.models <- run_models(
  'Volume',
  Volume ~ Close.Last.1 + Volume.1 + Volatility.1 + Spread.1 + Close.Last.2 + Volume.2 + Volatility.2 + Spread.2
)
spread.stock.models <- run_models(
  'Spread',
  Spread ~ Close.Last.1 + Volume.1 + Volatility.1 + Spread.1 + Close.Last.2 + Volume.2 + Volatility.2 + Spread.2
)

# Try models with both Tweet and stock data
close.last.both.models <- run_models(
  'Close.Last',
  Close.Last ~ User.Followers + User.Friends + User.Favourites + Retweets + Favorites + Tweets + Nchar + 
    Close.Last.1 + Volume.1 + Volatility.1 + Spread.1 + Close.Last.2 + Volume.2 + Volatility.2 + Spread.2
)
close.last.both <- metrics_comparison(close.last.both.models, 'Close.Last')
volatility.both.models <- run_models(
  'Volatility',
  Volatility ~ User.Followers + User.Friends + User.Favourites + Retweets + Favorites + Tweets + Nchar +
    Close.Last.1 + Volume.1 + Volatility.1 + Spread.1 + Close.Last.2 + Volume.2 + Volatility.2 + Spread.2
)
volatility.both <- metrics_comparison(volatility.both.models, 'Volatility')
volume.both.models <- run_models(
  'Volume',
  Volume ~ User.Followers + User.Friends + User.Favourites + Retweets + Favorites + Tweets + Nchar +
    Close.Last.1 + Volume.1 + Volatility.1 + Spread.1 + Close.Last.2 + Volume.2 + Volatility.2 + Spread.2
)
volume.both <- metrics_comparison(volume.both.models, 'Volume')
spread.both.models <- run_models(
  'Spread',
  Spread ~ User.Followers + User.Friends + User.Favourites + Retweets + Favorites + Tweets + Nchar +
    Close.Last.1 + Volume.1 + Volatility.1 + Spread.1 + Close.Last.2 + Volume.2 + Volatility.2 + Spread.2
)
spread.both <- metrics_comparison(spread.both.models, 'Spread')

# Appendix 3
close.last.return <- rbind(close.last.tweets[, 'Return'], close.last.stock[, 'Return'], close.last.both[, 'Return'])
rownames(close.last.return) <- c("Tweet data", "Stock data", "Combined data")
colnames(close.last.return) <- c("Linear Regression", "Random Forest", "CART")
md <-
  barplot(
    close.last.return,
    beside = T,
    legend = T,
    ylab = "Return (%)",
    main = "Model Return with Given Data",
    las = 1,
    col = c('grey85', 'grey65', 'grey45')
  )
text(md, close.last.return - 2, labels = paste(as.character(round(close.last.return, 2)), '%', sep = ''))

# Appendix 4
plot_comparison(close.last.both.models, 'Close.Last', title = "Closing Price Prediction Comparison")
plot_comparison(volatility.both.models, 'Volatility', title = "Volatility Prediction Comparison")
plot_comparison(volume.both.models, 'Volume', title = "Volume Prediction Comparison")
plot_comparison(spread.both.models, 'Spread', title = "Spread Prediction Comparison")

# Appendix 5
library(RColorBrewer)
osr2.both <-
  rbind(close.last.both[, 'OSR2'],
        volatility.both[, 'OSR2'],
        volume.both[, 'OSR2'],
        spread.both[, 'OSR2'])
rownames(osr2.both) <- c("Close.Last", "Volatility", "Volume", "Spread")
colnames(osr2.both) <- c("Linear Regression", "Random Forest", "CART")
barplot(
  osr2.both,
  beside = T,
  legend = T,
  ylab = "Out-of-Sample R^2",
  main = "Out-of-Sample Prediction Performance",
  las = 1,
  ylim = c(-1, 1),
  col = brewer.pal(4, "Set2")
)

# Appendix 6
library(rpart.plot)
close.last.variable.impact <-
  t(rbind(
    "Linear Model P-Value" = summary(close.last.both.models$lm)$coefficients[-1, 4],
    "Random Forest Importance" = t(close.last.both.models$rf$importance)
  ))
colnames(close.last.variable.impact) <- c("Linear Model P-Value", "Random Forest Importance")
prp(close.last.both.models$cart, type = 2, varlen = 0)
volatility.variable.impact <-
  t(rbind(
    "Linear Model P-Value" = summary(volatility.both.models$lm)$coefficients[-1, 4],
    "Random Forest Importance" = t(volatility.both.models$rf$importance)
  ))
colnames(volatility.variable.impact) <- c("Linear Model P-Value", "Random Forest Importance")
prp(volatility.both.models$cart, type = 2, varlen = 0)
volume.variable.impact <-
  t(rbind(
    "Linear Model P-Value" = summary(volume.both.models$lm)$coefficients[-1, 4],
    "Random Forest Importance" = t(volume.both.models$rf$importance)
  ))
colnames(volume.variable.impact) <- c("Linear Model P-Value", "Random Forest Importance")
prp(volume.both.models$cart, type = 2, varlen = 0)
spread.variable.impact <-
  t(rbind(
    "Linear Model P-Value" = summary(spread.both.models$lm)$coefficients[-1, 4],
    "Random Forest Importance" = t(spread.both.models$rf$importance)
  ))
colnames(spread.variable.impact) <- c("Linear Model P-Value", "Random Forest Importance")
prp(spread.both.models$cart, type = 2, varlen = 0)
