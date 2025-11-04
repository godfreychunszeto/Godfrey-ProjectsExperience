###############################################################################
# Load required libraries
library(quantmod)
library(ggplot2)
library(doParallel)

# Acquisition of computer cores
num_cores <- detectCores() - 1  # Set aside a core for other work
cl <- makeCluster(num_cores)
# Register for Parallel calculation
registerDoParallel(cl)

###############################################################################
# Functions:

# 1.Functions to calculate sigma of S&P500
# 20-year volatility
calculate_sigma_20y <- function(start_date, end_date) {
  sp500_data_20y <- readRDS("sp500_data_20y.rds")
  log_returns <- diff(log(Ad(sp500_data_20y)))
  sigma <- as.numeric(sd(log_returns, na.rm = TRUE) * sqrt(252)) #annual volatility
  return(as.numeric(sigma))
}
# 2-year volatility
calculate_sigma_2y <- function(start_date, end_date) {
  sp500_data_2y <- readRDS("sp500_data_2y.rds")
  log_returns <- diff(log(Ad(sp500_data_2y)))
  sigma <- as.numeric(sd(log_returns, na.rm = TRUE) * sqrt(252)) #annual volatility
  return(as.numeric(sigma))
}

# 2.Function to calculate historical annualised return（mu) of S&P500 
SP500_return_expectation <- function(start_date, end_date) {
  # Load S&P500 historical data
  sp500_data_20y <- readRDS("sp500_data_20y.rds")
  # Catch adjusted closed price from historical data，because no dividend data of S&P500 has been free provided externally
  # Thus, we use adjusted closed price to avoid calculating the dividend in all simulations
  adj_close <- Ad(sp500_data_20y)
  # Calculate the log rate of return
  log_returns <- diff(log(adj_close))
  # Calculate the total rate of return
  total_return <- (as.numeric(last(adj_close)) / as.numeric(first(adj_close))) - 1
  # Calculation of annualized rate of return
  T_days <- as.numeric(difftime(end_date, start_date, units = "days")) / 365  # Convert to number of years
  annualized_return <- (1 + total_return)^(1 / T_days) - 1
  return(annualized_return)
}

# 3.Function to generate standard normal variables (mean = 0) for all simulation to use the same
# Consistency and reproducibility of all simulation paths and steps
random <- function(paths, steps, seed = 4){
  set.seed(seed)
  all_Z <- rnorm(paths * steps)
}

# 4.Function for up-and-out European call option (Barrier option)
up_and_out_call <- function(S0, K, B_out, r, sigma, T = 1, n_paths = 20000, n_steps = 252) {
  dt <- T / n_steps
  payoffs <- numeric(n_paths)
  # Calling random function to generate standard normal variables
  all_Z <- random(n_paths, n_steps)
  dim(all_Z) <- c(n_steps, n_paths)
  
  payoffs <- foreach(i = 1:n_paths, .combine = c) %dopar% {
    S <- numeric(n_steps + 1)
    S[1] <- S0
    knocked_out <- FALSE
    
    for (j in 2:(n_steps + 1)) {
      Z <- all_Z[j-1, i]
      S[j] <- S[j - 1] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
      if (S[j] >= B_out) {
        knocked_out <- TRUE
        return(0)
      }
    }
    if (!knocked_out){
      return(max(S[n_steps + 1] - K, 0))
    }
  }
  price <- exp(-r * T) * mean(payoffs)
  return(price)
}

# 5.Function for down-and-in European put option (Barrier option)
down_and_in_put <- function(S0, K, B_in, r, sigma, T, n_paths = 20000, n_steps = 252) {
  dt <- T / n_steps
  payoffs <- numeric(n_paths)
  # Calling random function to generate standard normal variables
  all_Z <- random(n_paths, n_steps)
  dim(all_Z) <- c(n_steps, n_paths)
  
  payoffs <- foreach(i = 1:n_paths, .combine = c) %dopar% {
    S <- numeric(n_steps + 1)
    S[1] <- S0
    knock_in <- FALSE
    
    for (j in 2:(n_steps + 1)) {
      Z <- all_Z[j - 1, i]
      S[j] <- S[j - 1] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
      if (S[j] <= B_in) {
        knock_in <- TRUE
      } 
    }
    if (knock_in) {
      return(max(K - S[n_steps + 1], 0))
    } else {
      return(0)
    }
  }
  price <- exp(-r * T) * mean(payoffs)
  return(price)
} 

###############################################################################
# Structured Product Design:

# Download start date adjusted close price of S&P500 as Strike price
sp500_data_20y <- readRDS("sp500_data_20y.rds")
price_data <- if(file.exists("SP500_start_price.rds")) {
  readRDS("SP500_start_price.rds")
} else {
  readRDS("SP500_latest_price.rds")
}
cat("Current value of S&P500: $", price_data, "\n")

# Parameters (Very important !!!!!!!!!!!!!!!!!!!!!!!)
# Three scenarios: No recession / Slight recession / Strong recession
is_recession <- TRUE  # !!! TRUE if you think American economy will go into recession (S&P500 index will go down)
# !!! FALSE if it won't
Strong <- TRUE # !!! TRUE if you think the American economy will go into strong recession
# !!! FALSE if it is only slight recession
I <- 1000000 # Initial investment
Time <- 1 # The product period is 1 year
n_trials <- 20000 # Paths of the simulation
n_steps <- 252 * Time # Steps from t=0 to t=T (252 trading day)
dt <- Time/n_steps # Time interval between steps
r <- 0.0432 # Fed fund rate（risk-free rate)
r_1year <- 0.0392 # 1 year ZCB yield on 2025-04-04
P_ZCB <- exp(-r_1year * Time) # 1 year zero coupon bond
S0 <- start_price # Start price of S&P500 index in the product
K <- S0 # At the money options (K = S_0)
B_in <- 0.9 * K # Knock in with 90% of Strike price, due to betting on a decline of S&P500 index
start_date_option <- '2023-04-09' # Start date to measure the volatility for 1 year barrier option 
end_date_option <- '2025-04-09' # End date to measure the volatility for 1 year barrier option 
start_date <- "2005-04-09" # Start date to measure the expected return and volatility for simulate index trend
end_date <- "2025-04-09" # End date to measure the expected return and volatility for simulate index trend

# Calculate the S&P500 expected rate of return (20 years)
return_mu <- SP500_return_expectation(start_date, end_date)
cat("S&P500 annualized return from", start_date, "to", end_date, "is:", return_mu * 100, "%\n")

# Volatility of S&P500
# 2 years volatility
sigma <- calculate_sigma_2y(start_date_option, end_date_option)
cat("S&P500 volatility from historical data is (2 years) :", sigma, "\n")
# 20 years volatility
sigma_20yrs <- calculate_sigma_20y(start_date, end_date)
cat("S&P500 volatility from historical data is (20 years) :", sigma_20yrs, "\n")

Initial_ZCB <- I * P_ZCB # Initial investment in 1 year ZCB
residual_invest <- I - Initial_ZCB
cat("The knock-in price of barrier call is :", B_in, "\n")
DIP_price <- down_and_in_put(S0, K, B_in, r, sigma, Time, n_trials, n_steps)
cat("Down_and_In Put simulated price is :", DIP_price, "\n")
N <- floor(residual_invest / DIP_price)  # Number of down-and-in put we buy, short same number call
unused_cash <- residual_invest - N * DIP_price
index_down <- (K - B_in) / K * 100
cat("The put option will knock in if the index fall", index_down, "%", "\n")

B_out <- K + (0.20 * I) / N # Barrier to cap call maximum loss at 20% of Initial investment
cat("The knock-out price of barrier call is :", B_out, "\n")
UOC_price <- up_and_out_call(S0, K, B_out, r, sigma, Time, n_trials, n_steps)
cat("Up_and_Out Call simulated price is :", UOC_price, "\n")
index_up <- (B_out - K) / K
cat("The call option will knock out if the index raise", index_up * 100, "%", "\n")

stopCluster(cl)

# Participation rate
participation_rate <- (N * (DIP_price + UOC_price)) / residual_invest
cat("The participation rate is :", participation_rate * 100, "%", "\n")

# Fixed coupon rate
Real_ZCB_0 <- Initial_ZCB + N * UOC_price + unused_cash # Real investment in 1 year ZCB pay at t=0
Real_ZCB_T <- Real_ZCB_0 * exp(r_1year) # Amount of real investment in 1 year ZCB received at t=T
Fixed_coupon <- (Real_ZCB_T - I) / I
cat("Fixed coupon of the structure product is :", Fixed_coupon * 100 , "%", "\n")

# Maximum loss & return of the product
max_loss <- (Real_ZCB_T - index_up * K * N - DIP_price * N + UOC_price * N - I) / I
cat("The maximum loss of the product at the maturity is :", max_loss * 100, "%", "\n")
max_return <- (Real_ZCB_T + (K - 0) * N - DIP_price * N + UOC_price * N - I) / I
cat("The maximum return of the product at the maturity is :", max_return * 100, "%", "\n")

# Minimum return at the maturity (if S&P500 index < knock-in price at time T)
return_below_knockin <- (Real_ZCB_T + (K - B_in) * N - DIP_price * N + UOC_price * N - I) / I 
cat("The minimun return at the maturity when the S&P500 index is below the knock-in price :", return_below_knockin * 100, "%", "\n")

# Minimum return at the maturity (if S&P500 index <= strike price at time T)
return_same_as_strike <- (Real_ZCB_T - DIP_price * N + UOC_price * N - I) / I
cat("The minimum return at the maturity when S&P500 index is less than or equal to strike price :", return_same_as_strike * 100, "%", "\n")

# Fixed return at the maturity (if S&P500 index > knock-out price at time T)
return_above_knockout <- (Real_ZCB_T  - DIP_price * N + UOC_price * N - I) / I 
cat("The fixed return at the maturity when the S&P500 index is higher than the knock-out price :", return_above_knockout * 100, "%", "\n")

# Simulate the S&P500 index by different economic conditions you setted in the first and second rows of the parameters
return_rates <- numeric(n_trials)  # Final rate of return for each path
final_ST <- numeric(n_trials)      # Price at the maturity of each path
zone1 <- character(n_trials)       # Type 1 of potential return interval (depends on index level at the maturity)
zone2 <- character(n_trials)       # Type 2 of potential return interval (depends on potential return at the maturity)

# Calling random function to generate standard normal variables
all_Z <- random(n_trials, n_steps)
dim(all_Z) <- c(n_steps, n_trials)

# Impact of recession and S&P500 index downtrend on parameters
# If slightly recession，-0.5 * annualised return (mu) and 1.5 times of sigma from 20 years historical data
if (is_recession) {
  return_mu <- -0.5 * return_mu   # In a recession, mu is negative, simulating a market decline
  sigma_20yrs <- 1.5 * sigma_20yrs  # Increase volatility wo model heightened market uncertainty
  # If strong recession，-1 * annualised return (mu) and 2 times of sigma from 20 years historical data
  if (Strong){
    return_mu <- 2 * return_mu
    sigma_20yrs <- 4/3 * sigma_20yrs
  }
}

# Simulate each trials of S&P500 index
for (i in 1:n_trials) {
  S <- numeric(n_steps + 1)
  S[1] <- S0
  knock_in <- FALSE
  knocked_out <- FALSE
  
  for (j in 2:(n_steps + 1)) {
    Z <- all_Z[j - 1, i]
    S[j] <- S[j - 1] * exp((return_mu - 0.5 * sigma_20yrs^2) * dt + sigma_20yrs * sqrt(dt) * Z)
    # Determinate knock-in (Put void)
    if (S[j] <= B_in) knock_in <- TRUE
    # Determinate knock-out (Call void)
    if (S[j] >= B_out) knocked_out <- TRUE
  }
  
  # Final price using full simulation steps (not affected by knock-out)
  ST <- S[n_steps + 1]
  
  # Payoff of two barrier options
  put_payoff <- if (knock_in) max(K - ST, 0) * N else 0
  call_loss  <- if (!knocked_out) max(ST - K, 0) * N else 0
  
  # Potential return of structure product = Put payoff - Call payoff (due to sell a call) - LONG Put paid fee
  potential_return <- (put_payoff - call_loss - DIP_price * N ) * 100
  # If you want calculate and plot total return (potential return + fixed coupon) of structure product, 
  # Please use total return = Put payoff - Call payoff (due to sell a call) - LONG Put paid fee  + SHORT Call received fee:
  # Use code below replace original potential_return above:
  # total_return <- (Real_ZCB_T + put_payoff - call_loss - DIP_price * N - I) * 100
  
  # Potential return on structured products
  return_rates[i] <- potential_return / I
  # If want to calculate total return, please change the above row to the next row:
  # return_rates[i] <- total_return / I
  
  # Interval division (by S_T price)
  if (ST < B_in) {
    zone1[i] <- "S_T < B_in"
  } else if (ST < K) {
    zone1[i] <- "B_in ≤ S_T < K"
  } else if (ST < B_out) {
    zone1[i] <- "K ≤ S_T < B_out"
  } else {
    zone1[i] <- "S_T ≥ B_out"
  }
  
  # Interval diversion (by potential return)
  if (return_rates[i] > 20) {
    zone2[i] <- "higher than 20%"
  } else if (return_rates[i] > 4) {
    zone2[i] <- "4% ~ 20%"
  } else if (return_rates[i] > 0) {
    zone2[i] <- "0% ~ 4%"
  } else if (return_rates[i] > -2) {
    zone2[i] <- "0% ~ -2%"
  } else {
    zone2[i] <- "lower than -2%"
  }
}

# Combine results into a dataframe
results <- data.frame(Return = return_rates, Zone1 = zone1, Zone2 = zone2)

# Aggregate the average potential return and probability about the S_T interval
summary_zone1 <- aggregate(Return ~ Zone1, data = results, FUN = function(x) c(mean = mean(x), prob = length(x)/n_trials))
summary_zone1 <- do.call(data.frame, summary_zone1)
names(summary_zone1) <- c("Zone_1", "Mean_Potential_Return", "Probability")

# Aggregate the average potential return and probability about the potential return interval
summary_zone2 <- aggregate(Return ~ Zone2, data = results, FUN = function(x) c(mean = mean(x), prob = length(x)/n_trials))
summary_zone2 <- do.call(data.frame, summary_zone2)
names(summary_zone2) <- c("Zone_2", "Mean_Potential_Return", "Probability")

# Print results
print(summary_zone1)
print(summary_zone2)

###############################################################################
# Plot diagrams:

# Histogram of S_T price interval (Zone1) 
ggplot(results, aes(x = Zone1, fill = Zone1)) +
  geom_bar() +
  labs(title = "Distribution of S&P500 Index at S_T (20k trials)",
       x = "Interval Based on S_T", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")

# Histogram of potential return interval (Zone2) 
ggplot(results, aes(x = Zone2, fill = Zone2)) +
  geom_bar() +
  labs(title = "Distribution of Structured Product Potential Return (20k trials)",
       x = "Interval Based on Return", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")

###############################################################################
