# Load necessary libraries
library(parallel)
library(doParallel)
library(foreach)
library(doRNG)
library(MASS)
library(gamlss.data)
library(ggplot2)
library(dplyr)
library(asp23hmc)
library(lmls)
library(microbenchmark)
library(cli)
library(zoo)


# To recreate simulation, please enter your own output directory here:

# Define a function to unregister the parallel backend
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar



# Parameters
N <- c(50, 150, 300)
k <- c(2, 3)
rX <- c(0, 0.3, 0.7)
n_simulations <- 1000

# Initialize a list for storing simulation results
combin <- vector("list", length(N) * length(k) * length(rX))
names(combin) <- NULL

count <- 1

# Set a random seed for reproducibility
set.seed(420)

# Nested loops for simulations
for (Ni in N) {
  for (ki in k) {
    for (rx in rX) {
      # True parameter vectors
      beta_vec <- c(0.8, seq(length.out = ki, from = 1, to = 1.4))
      gamma_vec <- c(0.08, seq(length.out = ki, from = 0.06, to = 0.3))
      Jk <- matrix(data = 1, nrow = ki, ncol = ki)

      # Correlation matrices SigmaX
      SigmaX <- rx * Jk + (1 - rx) * diag(ki)

      # Loop for running simulations
      for (simulation in 1:n_simulations) {
        # Generate Design matrices
        Xs_withoutput_intercept <- mvrnorm(n = Ni, mu = rep(0, times = ki), Sigma = SigmaX)
        Xs <- cbind(1, Xs_withoutput_intercept)

        # Generate location and scale
        mu <- Xs %*% beta_vec
        sigma <- exp(Xs %*% gamma_vec)

        # Error vector
        error <- rnorm(n = Ni, mean = 0, sd = sigma)

        # Calculate y
        y = mu + error

        # Store simulation results
        combin[[count]][[simulation]] <- list(N = Ni, k = ki, rX = rx, beta_vec = beta_vec, gamma_vec = gamma_vec, Xs_withoutput_intercept = Xs_withoutput_intercept,
                                              Jk = Jk, SigmaX = SigmaX, Xs = Xs, mu = mu, sigma = sigma, y = y)
      }
      count <- count + 1
    }
  }
}

# Initialize output lists
output <- output_lmls <- vector("list", length(N) * length(k) * length(rX))

# Set up parallel processing
parallel::detectCores()
n.cores <- parallel::detectCores() - 1

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

output <- output_lmls <- vector("list", length(N) * length(k) * length(rX))

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParWorkers()

# Estiamte the coefficients with HMC for all combinations simulations
for (i in 1:(length(N) * length(k) * length(rX))) {

  output[[i]] <- list()  # Initialize each element as a list
  registerDoRNG(420)
  x<- foreach(sim = 1:n_simulations, .combine='rbind', .packages = c('HMC4', 'cli')) %dopar% {
    y <- combin[[i]][[sim]]$y
    x <- combin[[i]][[sim]]$Xs_withoutput_intercept
    z <- x
    res <- hmc(y ~ x, ~ z)
    unlist(res$coefficients)
  }
  output[[i]] <- x
}

# Estiamte the coefficients with HMC for all combinations simulations (baseline model)
for (i in 1:(length(N) * length(k) * length(rX))) {

  output_lmls[[i]] <- list()  # Initialize each element as a list
  registerDoRNG(420)
  x<- foreach(sim = 1:n_simulations, .combine='rbind', .packages = c('lmls')) %dopar% {
    y <- combin[[i]][[sim]]$y
    x <- combin[[i]][[sim]]$Xs_withoutput_intercept
    z <- x
    res <- lmls(y ~ x, ~ x, light = FALSE)

    unlist(res$coefficients)
  }
  output_lmls[[i]] <- x
}

# Stop parallel cluster
parallel::stopCluster(cl = my.cluster)

# Compute Monte Carlo Standard Errors (MCSE)
n_simulations

# Initialize lists for storing results
theta_hats <- theta_bars <- biases <- true_params <- MCSE <- list()
theta_hats_lmls <- theta_bars_lmls <- biases_lmls <- true_params_lmls <- MCSE_lmls <- list()

# Process simulation results
for (i in 1:(length(N) * length(k) * length(rX))){
  theta_hats[[i]] <- theta_hats_lmls[[i]] <- matrix(data = NA, nrow = n_simulations, ncol = dim(output[[i]])[2])
  for (sim in 1:n_simulations){
    theta_hats[[i]][sim,]  <- output[[i]][sim,]
    theta_hats_lmls[[i]][sim,]  <-  output_lmls[[i]][sim,]
  }
  theta_bars[[i]] <- colSums(theta_hats[[i]])/n_simulations
  theta_bars_lmls[[i]] <- colSums(theta_hats_lmls[[i]])/n_simulations
  true_params[[i]] <- c(combin[[i]][[1]]$beta_vec, combin[[i]][[1]]$gamma_vec)
  biases[[i]] <- theta_bars[[i]] - true_params[[i]]
  biases_lmls[[i]] <- theta_bars_lmls[[i]] - true_params[[i]]
}

# Check if we accidentally used wrong true parameters
for (i in 1:(length(N) * length(k) * length(rX))){
  print(c(combin[[i]][[1]]$beta_vec, combin[[i]][[1]]$gamma_vec))
}

# Compute Monte Carlo Standard Errors (MCSE)
MCSE

# Compute Monte Carlo Standard Errors for the baseline model
MCSE_lmls

# Define the output directory
output_directory <- "C:/Users/Soenke/OneDrive/Desktop/SoSe2023/Advanced Statistical Programming/hamiltonian-mcmc"

# Create the output directory if it doesn't exist
if (!file.exists(output_directory)) {
  dir.create(output_directory)
}

# Loop to save results to CSV files
for (i in 1:(length(N) * length(k) * length(rX))) {
  file_name <- paste0("matrix_", i, ".csv")

  write.csv(theta_hats[[i]], file.path(output_directory, paste("theta_hats",file_name)))
  write.csv(theta_bars[[i]], file.path(output_directory, paste("theta_bars",file_name)))
  write.csv(true_params[[i]], file.path(output_directory, paste("true_params",file_name)))
  write.csv(biases[[i]], file.path(output_directory, paste("biases",file_name)))
  write.csv(MCSE[[i]], file.path(output_directory, paste("MCSE",file_name)))

  write.csv(theta_hats_lmls[[i]], file.path(output_directory, paste("theta_hats_lmls",file_name)))
  write.csv(theta_bars_lmls[[i]], file.path(output_directory, paste("theta_bars_lmls",file_name)))
  write.csv(true_params[[i]], file.path(output_directory, paste("true_params__lmls",file_name)))
  write.csv(biases_lmls[[i]], file.path(output_directory, paste("biases_lmls",file_name)))
  write.csv(MCSE_lmls[[i]], file.path(output_directory, paste("MCSE_lmls",file_name)))
}

# Set the directory where the data files are located
data_dir <- "C:/Users/Soenke/OneDrive/Desktop/SoSe2023/Advanced Statistical Programming/hamiltonian-mcmc/matrizen"

# Create empty lists to store the data
theta_hats <- theta_bars <- true_params <- biases <- MCSE <- list()
theta_hats_lmls <- theta_bars_lmls <- true_params_lmls <- biases_lmls <- MCSE_lmls <- list()

plot_biases <- list()
plot_biases_lmls <- list()

# Loop to read in data from CSV files
for (i in 1:(length(N) * length(k) * length(rX))) {
  file_name <- paste0("matrix_", i, ".csv")

  theta_hats[[i]] <- read.csv(file.path(data_dir, paste("theta_hats", file_name)))
  theta_bars[[i]] <- read.csv(file.path(data_dir, paste("theta_bars", file_name)))
  true_params[[i]] <- read.csv(file.path(data_dir, paste("true_params", file_name)))
  biases[[i]] <- read.csv(file.path(data_dir, paste("biases", file_name)))
  MCSE[[i]] <- read.csv(file.path(data_dir, paste("MCSE", file_name)))
  theta_hats_lmls[[i]] <- read.csv(file.path(data_dir, paste("theta_hats_lmls", file_name)))
  theta_bars_lmls[[i]] <- read.csv(file.path(data_dir, paste("theta_bars_lmls", file_name)))
  true_params[[i]] <- read.csv(file.path(data_dir, paste("true_params__lmls", file_name)))

  theta_hats[[i]] <- theta_hats[[i]][,-1]
  theta_bars[[i]] <- theta_bars[[i]][, -1]
  true_params[[i]] <- true_params[[i]][, -1]
  biases[[i]] <- biases[[i]][, -1]
  MCSE[[i]] <- MCSE[[i]][,-1]
  theta_hats_lmls[[i]] <- theta_hats_lmls[[i]][,-1]
  theta_bars_lmls[[i]] <- theta_bars_lmls[[i]][,-1]
  biases_lmls[[i]] <- biases_lmls[[i]][,-1]
  MCSE_lmls[[i]] <- MCSE_lmls[[i]][,-1]

  plot_biases[[i]] <- sweep(theta_hats[[i]], 2, true_params[[i]])
  plot_biases_lmls[[i]]  <- sweep(theta_hats_lmls[[i]], 2, true_params[[i]])
}

# Loop to create and display comparison plots
for (i in 1:(length(N) * length(k) * length(rX))) {
  plot_list <- list()

  biases <- plot_biases[[i]]
  biases_lmls <- plot_biases_lmls[[i]]

  num_parameters <- ncol(biases)

  for (j in 1:num_parameters) {
    data_df <- data.frame(Bias = biases[, j], Method = "HMC", Combination = paste("Combination", i), Parameter =  ifelse(j <= num_parameters/2, paste("Beta Hat", j-1), paste("Gamma Hat", j-num_parameters/2-1)))
    data_df_lmls <- data.frame(Bias = biases_lmls[, j], Method = "Baseline", Combination = paste("Combination", i), Parameter =  ifelse(j <= num_parameters/2, paste("Beta Hat", j-1), paste("Gamma Hat", j-num_parameters/2-1)))

    combined_data <- rbind(data_df, data_df_lmls)

    p <- ggplot(combined_data, aes(x = Bias, fill = Method)) +
      geom_histogram(binwidth = 0.1, position = "dodge") +
      facet_wrap(~ Parameter) +
      labs(title = paste("Bias Comparison - Combination", i)) +
      theme_minimal()

    plot_list[[j]] <- p
  }

  pdf(file.path(output_directory, paste("Bias_Plots_Combination_", i, ".pdf")))
  for (j in 1:num_parameters) {
    print(plot_list[[j]])
  }
  dev.off()
}
