# Standard plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.hmc <- function(m, predictor = c('location', 'scale'), type = c('trace','hist','density','acf', 'roll_avg'), exclude_warmup = TRUE){

  # Match arguments
  type <- match.arg(type)
  predictor <- match.arg(predictor)

  # Define type based on input and accept only one parameter
  predictor <- match.arg(predictor)

  # Get number of coefficients
  num_coef <- length(m$coefficients[[predictor]])
  var_names <- names(coef(output)[[predictor]])

  # Define matrix
  mat <- m$hmc[[predictor]]
  mat_pruned <- m$hmc[[predictor]][-(1:(m$hmc$end_warmup+1)),]


  #Define plot grid
  par(mfrow = c(ceiling(num_coef/2),2))

  if (type == 'trace') {
    plot_trace(m, num_coef, var_names, mat)

  } else if (type == 'hist') {
    if (exclude_warmup) {
      plot_hist(m, num_coef, var_names, mat_pruned)
    } else{
      plot_hist(m, num_coef, var_names, mat)
    }

  } else if (type == 'density') {
    if (exclude_warmup) {
      plot_density(m, predictor, num_coef, var_names, mat_pruned)
    } else{
      plot_density(m, predictor, num_coef, var_names, mat)
    }
  } else if (type == 'acf') {
    if (exclude_warmup) {
      plot_acf(num_coef, var_names, mat_pruned)
    } else{
      plot_acf(num_coef, var_names, mat)
    }
  } else {
    # Case for rolling average
    plot_roll_mean(m, predictor, num_coef, var_names, mat)
  }
}

plot_trace <- function(m, num_coef, var_names, data){
  for (i in 1:num_coef) {
    x_values <- 0:(length(data[,i]) - 1)  # Create x-values starting from 0
    plot(x_values, data[,i], type = 'l', xlab = "Iteration", ylab = 'value', main = var_names[i],
         xlim = c(0, length(data[,i]) - 1))  # Set x-axis limits to start at 0
    abline(v = m$hmc$total_adapt, col = "red", lty = 2)  # Add a vertical red line at num_adapt
    abline(v = m$hmc$end_warmup, col = "blue", lty = 3)  # Add a vertical blue dotted line at end_warm_up
    abline(h = data[1, i], col = "black", lty = 3)  # Add a horizontal dashed line for initialization
  }

  legend("bottomright", bty = "n", inset=c(0,0),
         legend = c("End Adaptation", "End Warmup", "Initialization"),
         col = c("red", "blue", "black"), lty = c(2, 2, 3), cex = 0.6)
}

plot_hist <- function(m, num_coef, var_names, data){
  for (i in 1:num_coef) {
    hist(data[, i], xlab = var_names[i], ylab = 'Frequency', main = var_names[i])
  }
}

plot_density <- function(m, predictor, num_coef, var_names, data){
  for (i in 1:num_coef){
    plot(density(data[, i]), main = var_names[i])
    # Add a vertical line for the mean
    abline(v = m$coefficients[[predictor]][i], col = "red")}
}

plot_acf <- function(num_coef, var_names, data){
  for (i in 1:num_coef) {
    acf(data[,i], main = var_names[i])
  }
}


plot_roll_mean <- function(m, predictor, num_coef, var_names, data){
  for(i in 1:num_coef){
    # Compute variables based on predictor
    if(predictor == 'location'){
      rolling_average <- m$hmc$rolling_average_location[-1, i, drop = FALSE]
      rolling_average_without_warm_up <- m$hmc$rolling_average_location_without_warmup[-1, i, drop = FALSE]

      # Combine both rolling averages and initial estimation values to determine y-axis limits
      combined_data <- c(rolling_average, rolling_average_without_warm_up, data[1, i])
    }else{
      rolling_average <- m$hmc$rolling_average_scale[-1, i, drop = FALSE]
      rolling_average_without_warm_up <- m$hmc$rolling_average_scale_without_warmup[-1, i, drop = FALSE]

      # Combine both rolling averages and initial estimation values to determine y-axis limits
      combined_data <- c(rolling_average, rolling_average_without_warm_up, data[1, i])
    }

    y_min <- min(combined_data, na.rm = TRUE)
    y_max <- max(combined_data, na.rm = TRUE)

    # Set y-axis limits explicitly to ensure the horizontal dashed line is visible
    ylim <- c(y_min - 0.01 * abs(y_min), y_max + 0.01 * abs(y_max))

    matplot(1:m$hmc$num_samples, cbind(rolling_average, rolling_average_without_warm_up),
            type = "l", xlab = "Iteration", ylab = "Rolling Average",
            main = var_names[i], ylim = ylim, col = c("blue", "darkgreen"), lty=c(1))

    abline(v = m$hmc$total_adapt, col = "red", lty = 2)  # Add a vertical red line at num_adapt
    abline(v = m$hmc$end_warmup, col = "blue", lty = 3)  # Add a vertical blue dotted line at end_warm_up
    abline(h = data[1, i], col = "black", lty = 3)  # Add a horizontal dashed line for initialization
  }
  legend("bottomright", bty = "n", inset = 0.01,
         legend = c("End Adaptation", "End Warmup", "Initialization", "Rolling Average (All samples)", "Rolling Average (Without Warmup)"),
         col = c("red", "blue", "black", "blue", "darkgreen"), lty = c(2, 3, 3, 1, 1), cex = 0.6) # Create a legend
}
