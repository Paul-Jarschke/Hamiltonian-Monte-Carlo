#' Extract Coefficients from an HMC Model
#' This function extracts coefficients from an HMC model object.
#'
#' @param object An HMC model object.
#' @param predictor A character vector specifying the predictor(s) for which coefficients are extracted.
#'   Options include "location" and "scale." Default is c("location", "scale"), which returns both location and scale coefficients.
#' @param ... Additional arguments (not used).
#'
#' @return A list of coefficients corresponding to the specified predictor(s).
#'
#' @rdname hmc-methods
#'
#' @export
#'
coef.hmc <- function(object, predictor = c("location", "scale"), ...) {
  # Match the specified predictor(s) to the available options
  predictor <- match.arg(predictor, several.ok = TRUE)

  if (length(predictor) == 1) {
    # Return coefficients for a single predictor
    object$coefficients[[predictor]]
  } else {
    # Return coefficients for multiple predictors
    object$coefficients[predictor]
  }
}


#' Extract Log-Likelihood from an HMC Model
#' This function extracts the log-likelihood from an HMC model object.
#'
#' @param object An HMC model object.
#' @param ... Additional arguments (not used).
#'
#' @importFrom stats nobs
#'
#' @export
#'
logLik.hmc <- function(object, ...) {
  # Calculate the log-likelihood using the loglik method for the HMC model object
  out <- loglik(object)

  # Set attributes for degrees of freedom and number of observations
  attr(out, "df") <- object$df
  attr(out, "nobs") <- nobs(object)

  # Set the class of the result to "logLik"
  class(out) <- "logLik"

  # Return the log-likelihood object
  out
}



#' Plot Method for an HMC Model
#'
#' This function generates various types of plots for an HMC model object.
#'
#' @param x An HMC model object.
#' @param predictor A character vector specifying the predictor for which plots are generated.
#'   Options include "location" and "scale." Default is c("location", "scale"), which generates plots for both location and scale predictors.
#' @param type A character vector specifying the type of plots to generate. Options include "trace," "hist," "density," "acf," and "roll_avg."
#' @param exclude_warmup A logical value indicating whether to exclude warm-up iterations when generating plots. Default is TRUE.
#' @param ... Additional arguments (not used).
#'
#' @export
#'
plot.hmc <- function(x, predictor = c('location', 'scale'), type = c('trace','hist','density','acf', 'roll_avg'), exclude_warmup = TRUE, ...){
  # Match arguments
  type <- match.arg(type)
  predictor <- match.arg(predictor)

  # Get number of coefficients and their names
  num_coef <- length(x$coefficients[[predictor]])
  var_names <- names(coef(x)[[predictor]])

  # Define matrix
  mat <- x$hmc[[predictor]]
  mat_pruned <- x$hmc[[predictor]][-(1:(x$hmc$end_warmup+1)),]

  # Define plot grid
  par(mfrow = c(ceiling(num_coef/2), 2))

  if (type == 'trace') {
    # Generate trace plots
    plot_trace(x, num_coef, var_names, mat)
  } else if (type == 'hist') {
    if (exclude_warmup) {
      # Generate histogram plots (excluding warm-up)
      plot_hist(x, num_coef, var_names, mat_pruned)
    } else {
      # Generate histogram plots (including warm-up)
      plot_hist(x, num_coef, var_names, mat)
    }
  } else if (type == 'density') {
    if (exclude_warmup) {
      # Generate density plots (excluding warm-up)
      plot_density(x, predictor, num_coef, var_names, mat_pruned)
    } else {
      # Generate density plots (including warm-up)
      plot_density(x, predictor, num_coef, var_names, mat)
    }
  } else if (type == 'acf') {
    if (exclude_warmup) {
      # Generate autocorrelation plots (excluding warm-up)
      plot_acf(num_coef, var_names, mat_pruned)
    } else {
      # Generate autocorrelation plots (including warm-up)
      plot_acf(num_coef, var_names, mat)
    }
  } else {
    # Case for rolling average plots
    plot_roll_mean(x, predictor, num_coef, var_names, mat)
  }
}




#' Print Method for an HMC Model
#'
#' This function provides a customized print method for an HMC model object.
#'
#' @param x An HMC model object.
#' @param digits The number of digits to be used for printing numeric values. Default is max(3, getOption("digits") - 3).
#' @param ... Additional arguments (not used).
#'
#' @importFrom stats coef
#'
#' @export
#'
print.hmc <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  # Print the call
  cat(
    "\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )

  pretty <- function(x) {
    # Define a custom printing function
    print.default(format(x, digits = digits), print.gap = 2, quote = FALSE)
  }

  if (length(coef(x, "location"))) {
    # Print location coefficients with identity link
    cat("Location coefficients (identity link):\n")
    pretty(coef(x, "location"))
  } else {
    cat("No location coefficients\n")
  }

  cat("\n")

  if (length(coef(x, "scale"))) {
    # Print scale coefficients with log link
    cat("Scale coefficients (log link):\n")
    pretty(coef(x, "scale"))
  } else {
    cat("No scale coefficients\n")
  }

  cat("\n")

  # Make the function invisible and return the original object
  invisible(x)
}


#' Summary Method for an HMC Model
#'
#' This function provides a summary of an HMC model object.
#'
#' @param object An HMC model object.
#' @param type A character string specifying the type of summary to generate. Default is 'hmc'.
#' @param digits The number of digits to be used for printing numeric values. Default is max(3, getOption("digits") - 3).
#' @param ... Additional arguments (not used).
#'
#' @return A summary of the HMC model object.
#'
#' @export
summary.hmc <- function(object, type = 'hmc',
                        digits = max(3, getOption("digits") - 3), ...) {

  # Define a function to extract coefficient matrices based on predictor
  coefmat <- function(m, predictor) coefmat_samples(m, predictor, type)

  cat(
    "\nCall:\n",
    paste(deparse(object$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )

  if (length(object$coefficients$location)) {
    cat("Location coefficients (identity link):\n")
    printCoefmat(coefmat(object, "location"), digits = digits, signif.stars = FALSE, ...)
  } else {
    cat("No location coefficients\n")
  }

  cat("\n")

  if (length(object$coefficients$scale)) {
    cat("Scale coefficients (log link):\n")
    printCoefmat(coefmat(object, "scale"), digits = digits, signif.stars = FALSE, ...)
  } else {
    cat("No scale coefficients\n")
  }

  cat("\n")

  pretty <- function(x, y) {
    cat(x, ": ", format(signif(y), getOption("digits")), "\n", sep = "")
  }

  cat("Acceptance ratios:\n")
  print(object$acceptance_ratios)
  cat("\n")
  cat("Hyperparameter:\n")
  print(object$hyperparameter)
  cat("\n")

  pretty("Residual degrees of freedom", df.residual(object))
  pretty("Log-likelihood", loglik(object))
  pretty("AIC", AIC(logLik.hmc(object)))
  pretty("BIC", BIC(logLik.hmc(object)))

  cat("\n")

  invisible(object)
}

# method helpers ---------------------------------------------------------------

coefmat <- function(m, predictor) {
  m <- tidy(m, predictor)
  out <-as.matrix(m[,(2:5)])
  colnames(out) <- c("Mean", "2.5%", "50%", "97.5%")
  rownames(out) <- m$term
  out
}

#' @importFrom stats quantile

coefmat_samples <- function(m, predictor, type) {
  samples <- m[[type]][[predictor]]

  coefmat <- apply(samples, 2, function(x) {
    c(Mean = mean(x), quantile(x, c(0.025, 0.5, 0.975)))
  })

  t(coefmat)
}

# summary helpers --------------------------------------------------------------

acceptance_ratios <- function(m, num_adapt, num_samples){
  df <- data.frame(
             Adaptation = c(sum(m$hmc$accepted_beta[1:num_adapt+1])/(num_adapt),
                            sum(m$hmc$accepted_gamma[1:num_adapt+1])/(num_adapt)),
             Tuned = c(sum(m$hmc$accepted_beta[(num_adapt+1):(num_samples)+1])/(num_samples-num_adapt),
                       sum(m$hmc$accepted_gamma[(num_adapt+1):(num_samples)+1])/(num_samples-num_adapt)),
             Overall = c(sum(m$hmc$accepted_beta[1:num_samples+1])/(num_samples),
                         sum(m$hmc$accepted_gamma[1:num_samples+1])/(num_samples)),
             row.names = c('location','scale'))
  df
}


hyperparameter <- function(epsilon_beta, epsilon_gamma, L_beta, L_gamma, LAMBDA){
    df <- data.frame(
      epsilon = c(epsilon_beta, epsilon_gamma),
      L_tuned = c(round(LAMBDA/epsilon_beta), round(LAMBDA/epsilon_gamma)),
      L_realized = c(L_beta, L_gamma),
      row.names = c('location','scale'))
    df
}


tidy <- function(x, predictor = c("location", "scale"), ...) {
  predictor <- match.arg(predictor, several.ok = TRUE)

  if (length(predictor) > 1) {
    out <- lapply(predictor, function(p) tidy(x, p))
    out <- do.call(rbind, out)
  } else if(predictor == 'location') {
    out <- data.frame(
      term  = colnames(x$x),
      Mean  = coef.hmc(x, predictor),
      '2.5%' = apply(x$hmc[[predictor]], 2, quantile, probs = 0.025),
      '50%' = apply(x$hmc[[predictor]], 2, quantile, probs = 0.5),
      '97.5%'  = apply(x$hmc[[predictor]], 2, quantile, probs = 0.975),
      row.names = NULL)
  }else{
    out <- data.frame(
      term  = colnames(x$z),
      Mean  = coef.hmc(x, predictor),
      '2.5%' = apply(x$hmc[[predictor]], 2, quantile, probs = 0.025),
      '50%' = apply(x$hmc[[predictor]], 2, quantile, probs = 0.5),
      '97.5%'  = apply(x$hmc[[predictor]], 2, quantile, probs = 0.975),
      row.names = NULL)
  }

  out
}

