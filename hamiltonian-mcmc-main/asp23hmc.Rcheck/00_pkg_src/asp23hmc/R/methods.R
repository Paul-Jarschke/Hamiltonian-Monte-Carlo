#' Methods for hmc
#'
#' A couple of methods for location-scale regression models from the [hmc()]
#' function are provided.
#'
#' @param object A location-scale regression model from the [hmc()] function.
#' @param predictor The predictor to work on. Either `"location"` or `"scale"`
#'                  or both. If both, a list with the names `"location"` and
#'                  `"scale"` is returned.
#' @param ... Currently ignored.
#' @param newdata A data frame (or list or environment) with the covariate
#'                values at which the predictions are computed. If `NULL`, the
#'                predictions at the original data are returned.
#' @param type Used by `predict()` and `residuals()`:
#'             \itemize{
#'               \item For `predict()`, `"link"` or `"response"`. If `"link"`
#'                     (default), \eqn{\mu} and log(\eqn{\sigma}) are returned.
#'                     If `"response"`, \eqn{\mu} and \eqn{\sigma}
#'                     are returned.
#'               \item For `residuals()`, `"deviance"`, `"pearson"` or
#'                     `"response"`. If `"deviance"` (default) or `"pearson"`,
#'                     (\eqn{y - \mu}) / \eqn{\sigma} is returned.
#'                     If `"response"`, \eqn{y - \mu} is returned.
#'             }
#'
#' @return
#'
#' A numeric vector for `residuals()`. For the other methods, a numeric vector
#' if the argument `predictor` is either `"location"` or `"scale"`, or a list
#' with the names `location` and `scale` if it is both.
#'
#' @name hmc-methods

NULL

#' @rdname hmc-methods
#' @export

coef.hmc <- function(object, predictor = c("location", "scale"), ...) {
  predictor <- match.arg(predictor, several.ok = TRUE)

  if (length(predictor) == 1) {
    object$coefficients[[predictor]]
  } else {
    object$coefficients[predictor]
  }
}

#' @importFrom stats resid
#' @export

deviance.hmc <- function(object, ...) {
  sum(resid(object, "deviance")^2)
}

#' @rdname hmc-methods
#' @export

fitted.hmc <- function(object, predictor = c("location", "scale"), ...) {
  predictor <- match.arg(predictor, several.ok = TRUE)

  if (length(predictor) == 1) {
    object$fitted.values[[predictor]]
  } else {
    object$fitted.values[predictor]
  }
}

#' @importFrom stats nobs
#' @export

logLik.hmc <- function(object, ...) {
  out <- loglik(object)

  attr(out, "df") <- object$df
  attr(out, "nobs") <- nobs(object)
  class(out) <- "logLik"

  out
}

#' @importFrom graphics abline plot
#' @importFrom stats fitted qnorm resid
#' @export

plot.hmc <- function(x,
                      xlab = "Fitted values",
                      ylab = "Deviance residuals",
                      ...) {
  plot(
    x = fitted.hmc(x, "location"),
    y = resid(x, "deviance"),
    xlab = xlab,
    ylab = ylab,
    ...
  )

  abline(h = qnorm(c(0.025, 0.975)), lty = "dashed")
  abline(h = 0)

  invisible(x)
}

#' @rdname hmc-methods
#' @importFrom stats as.formula coef fitted model.matrix predict update
#' @export

predict.hmc <- function(object,
                         newdata = NULL,
                         predictor = c("location", "scale"),
                         type = c("link", "response"),
                         ...) {
  predictor <- match.arg(predictor, several.ok = TRUE)
  type <- match.arg(type)

  if (length(predictor) > 1) {
    out <- lapply(predictor, function(p) {
      predict(object, newdata, p, type)
    })

    names(out) <- predictor

    return(out)
  }

  if (is.null(newdata)) {
    out <- fitted.hmc(object, predictor)

    if (predictor == "scale" && type == "link") {
      out <- log(out)
    }

    return(out)
  }

  arg <- object$call[[predictor]]
  formula <- update(as.formula(arg), NULL ~ .)
  x <- model.matrix(formula, newdata)

  beta <- coef(object, predictor)
  out <- drop(x %*% beta)

  if (predictor == "scale" && type == "response") {
    out <- exp(out)
  }

  out
}

#' @importFrom stats coef
#' @export

print.hmc <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat(
    "\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )

  pretty <- function(x) {
    print.default(format(x, digits = digits), print.gap = 2, quote = FALSE)
  }

  if (length(coef(x, "location"))) {
    cat("Location coefficients (identity link):\n")
    pretty(coef(x, "location"))
  } else {
    cat("No location coefficients\n")
  }

  cat("\n")

  if (length(coef(x, "scale"))) {
    cat("Scale coefficients (log link):\n")
    pretty(coef(x, "scale"))
  } else {
    cat("No scale coefficients\n")
  }

  cat("\n")

  invisible(x)
}

#' @importFrom stats qqline qqnorm resid
#' @export

qqnorm.hmc <- function(y,
                        xlab = "Theoretical quantiles",
                        ylab = "Deviance residuals",
                        ...) {
  qqnorm(resid(y, "deviance"), xlab = xlab, ylab = ylab, ...)
  qqline(resid(y, "deviance"))

  invisible(y)
}

#' @rdname hmc-methods
#' @importFrom stats fitted
#' @export

residuals.hmc <- function(object,
                           type = c("deviance", "pearson", "response"),
                           ...) {
  type <- match.arg(type)
  out <- object$residuals

  if (type != "response") {
    out <- out / fitted.hmc(object, "scale")
  }

  out
}

#' @importFrom stats fitted nobs rnorm runif
#' @export

simulate.hmc <- function(object, nsim = 1, seed = NULL, ...) {
  # RNG handling taken from simulate.lm():
  # https://github.com/wch/r-source/blob/master/src/library/stats/R/lm.R

  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    runif(1)  # initialize the RNG if necessary
  }

  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  n <- nobs(object) * nsim
  out <- rnorm(n, fitted.hmc(object, "location"), fitted.hmc(object, "scale"))
  out <- matrix(out, nrow = nobs(object), ncol = nsim)

  out <- as.data.frame(out)
  names(out) <- paste("sim", seq_len(nsim), sep = "_")
  attr(out, "seed") <- RNGstate

  out
}

#' Summary for hmc
#'
#' Prints a summary for location-scale regression models from the [hmc()]
#' function.
#'
#' @param object A location-scale regression model from the [hmc()] function.
#' @param type Either `"ml"` or `"boot"` or `"mcmc"`:
#'             \itemize{
#'               \item If `"ml"`, the maximum likelihood estimates and the
#'                     asymptotic standard errors are shown.
#'               \item If `"boot"`, the bootstrap estimates and confidence
#'                     intervals are shown.
#'               \item If `"mcmc"`, the Markov chain Monte Carlo (MCMC)
#'                     estimates and credible intervals are shown.
#'             }
#' @param digits The number of digits to print.
#' @param ... Passed on to [printCoefmat()].
#'
#' @return
#'
#' The (unmodified and invisible) `hmc` S3 object, see [hmc()].
#'
#' @importFrom stats AIC BIC coef df.residual printCoefmat resid
#' @export

summary.hmc <- function(object,
                         type = c("ml", "boot", "mcmc", "hmc"),
                         digits = max(3, getOption("digits") - 3),
                         ...) {
  type <- match.arg(type)

  if (type == "hmc") {
    coefmat <- function(m, predictor) coefmat_samples(m, predictor, type)

    if (is.null(object[[type]]$location) || is.null(object[[type]]$scale)) {
      stop("Model does not include samples, run ", type, "() first")
    }
  }

  cat(
    "\nCall:\n",
    paste(deparse(object$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )

 # if (df.residual(object) > 5) {
  #  cat("Deviance residuals:\n")
  #  print.hmc(summary(resid(object, "deviance"), digits = digits))
  #  cat("\n")
  #}

  if (length(object$coefficients$location)) {
    cat("Location coefficients (identity link):\n")
    printCoefmat(coefmat(object, "location"), digits = digits, signif.stars = FALSE, ...)
  } else {
    cat("No location coefficients\n")
  }

  cat("\n")

  if (length(object$coefficients$scale)) {
    cat("Scale coefficients (log link):\n")
    printCoefmat(coefmat(object, "scale"), digits = digits, signif.stars = FALSE,...)
  } else {
    cat("No scale coefficients\n")
  }

  cat("\n")

  pretty <- function(x, y) {
    cat(x, ": ", format(signif(y), getOption("digits")), "\n", sep = "")
  }

  #cat("\n")

  cat("Acceptance ratios:\n")
  print(acceptance_ratios.hmc(object))

  cat("\n")


  pretty("Residual degrees of freedom", df.residual(object))
  pretty("Log-likelihood", loglik(object))
  pretty("AIC", AIC(logLik.hmc(object)))
  pretty("BIC", BIC(logLik.hmc(object)))

  cat("\n")

  invisible(object)
}

#' @rdname hmc-methods
#' @export

vcov.hmc <- function(object, predictor = c("location", "scale"), ...) {
  predictor <- match.arg(predictor, several.ok = TRUE)

  if (length(predictor) == 1) {
    object$vcov[[predictor]]
  } else {
    object$vcov[predictor]
  }
}




coefmat <- function(m, predictor) {
  m <- tidy.hmc(m, predictor)
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


tidy.hmc <- function(x, predictor = c("location", "scale"), ...) {
  predictor <- match.arg(predictor, several.ok = TRUE)

  if (length(predictor) > 1) {
    out <- lapply(predictor, function(p) tidy(x, p))
    out <- do.call(rbind, out)
  } else {
    out <- data.frame(
      term      = colnames(m$x),
      Mean  = coef.hmc(x, predictor),
      '2.5%' = apply(x$hmc[[predictor]], 2, quantile, probs = 0.025),
      '50%' = apply(x$hmc[[predictor]], 2, quantile, probs = 0.5),
      '97.5%'  = apply(x$hmc[[predictor]], 2, quantile, probs = 0.975),
      row.names = NULL
    )
  }

  out
}

acceptance_ratios.hmc <- function(m){
  df <- data.frame(
             Adaptation = c(m$hmc$acceptance_rate_beta_adapt_phase,
                          m$hmc$acceptance_rate_gamma_adapt_phase),
             Tuned = c(m$hmc$acceptance_rate_beta_tuned,
                       m$hmc$acceptance_rate_gamma_tuned),
             Overall = c(m$hmc$acceptance_rate_beta,
                         m$hmc$acceptance_rate_gamma),
             row.names = c('location','scale'))
  df
}

hyperparameter_df <- function(epsilon_beta, epsilon_gamma, L_beta, L_gamma, LAMBDA){
    df <- data.frame(
      epsilon_tuned = c(epsilon_beta, epsilon_gamma),
      L_tuned = c(round(LAMBDA/epsilon_beta), round(LAMBDA/epsilon_gamma)),
      L_realized = c(L_beta, L_gamma)
      )
    df
}




vcov.hmc <- function(object, predictor = c("location", "scale"), ...) {
  predictor <- match.arg(predictor, several.ok = TRUE)

  if (length(predictor) == 1) {
    object$vcov[[predictor]]
  } else {
    object$vcov[predictor]
  }
}




logLik.hmc <- function(object, ...) {
  out <- loglik(object)

  attr(out, "df") <- object$df
  attr(out, "nobs") <- nobs(object)
  class(out) <- "logLik"

  out
}

