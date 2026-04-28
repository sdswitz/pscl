anova.zeroinfl <- function(object, ..., test = c("Chisq", "none", "LRT")) {
  if (is.null(test)) test <- "Chisq"
  else if (isFALSE(test)) test <- "none"
  else test <- match.arg(test)
  if (test == "LRT") test <- "Chisq"
  models <- list(object, ...)
  if (length(models) < 2L) {
    stop("anova.zeroinfl() currently supports only multi-model comparisons: anova(m1, m2, ...).")
  }
  
  # Basic class check
  bad <- !vapply(models, inherits, logical(1), what = "zeroinfl")
  if (any(bad)) stop("All models must be of class 'zeroinfl'.")

  responses <- vapply(models, function(x) paste(deparse(formula(x)[[2L]]), collapse = "\n"), character(1))
  if (!all(responses == responses[1L])) {
    stop("Models must have the same response.")
  }
  
  # Extract logLik objects (pscl provides logLik.zeroinfl with df and nobs attrs)
  ll <- lapply(models, logLik)
  df <- vapply(ll, function(x) as.integer(attr(x, "df")), integer(1))
  n  <- vapply(ll, function(x) as.integer(attr(x, "nobs")), integer(1))
  
  if (length(unique(n)) != 1L) {
    stop("Models were not fit to the same number of observations (nobs differs).")
  }
  
  # Order by increasing df (small -> large)
  ord <- order(df)
  if (!identical(ord, seq_along(models))) {
    warning("Models reordered by increasing number of parameters for comparison.")
    models <- models[ord]
    ll <- ll[ord]
    df <- df[ord]
  }
  
  # Check nestedness: smaller model terms must be a subset of larger model terms
  for (i in seq_along(models)[-1]) {
    small <- models[[i - 1]]
    large <- models[[i]]
    
    count_nested <- all(attr(small$terms$count, "term.labels") %in%
                          attr(large$terms$count, "term.labels"))
    zero_nested  <- all(attr(small$terms$zero, "term.labels") %in%
                          attr(large$terms$zero, "term.labels"))
    
    if (!count_nested || !zero_nested) {
      warning("Models ", i - 1, " and ", i,
              " may not be nested; LRT results may be invalid.")
    }
  }
  
  loglik <- vapply(ll, as.numeric, numeric(1))
  
  # Residual df = nobs - #params; deviance = -2 * logLik
  nobs    <- n[1]
  resdf   <- nobs - df
  resdev  <- -2 * loglik
  
  # Successive differences (NA for first row, like anova.glm)
  df_diff  <- c(NA_integer_, diff(df))
  dev_diff <- c(NA_real_, -diff(resdev))  # reduction in deviance (positive = improvement)
  
  pval <- rep(NA_real_, length(models))
  if (test == "Chisq") {
    for (i in seq_along(models)[-1]) {
      if (df_diff[i] > 0 && dev_diff[i] >= 0) {
        pval[i] <- pchisq(dev_diff[i], df = df_diff[i], lower.tail = FALSE)
      } else if (!is.na(dev_diff[i]) && dev_diff[i] < 0) {
        warning("Negative LRT statistic for models ", i - 1, " vs ", i,
                "; models may not be properly nested.")
      }
    }
  }
  
  tab <- data.frame(
    `Resid. Df`  = resdf,
    `Resid. Dev` = resdev,
    Df           = df_diff,
    Deviance     = dev_diff,
    `Pr(>Chi)`   = pval,
    check.names  = FALSE
  )
  rownames(tab) <- seq_along(models)
  
  # Build heading with model formulas, expanding "." to actual terms
  formulas <- vapply(models, function(m) {
    resp <- all.vars(m$formula)[1]
    count_terms <- attr(m$terms$count, "term.labels")
    zero_terms  <- attr(m$terms$zero, "term.labels")
    count_rhs <- if (length(count_terms)) paste(count_terms, collapse = " + ") else "1"
    zero_rhs  <- if (length(zero_terms)) paste(zero_terms, collapse = " + ") else "1"
    paste0(resp, " ~ ", count_rhs, " | ", zero_rhs)
  }, character(1))
  model_lines <- paste0("Model ", seq_along(models), ": ", formulas)
  heading <- paste0(
    "Likelihood Ratio Tests of Zero-Inflated Count Regression Models\n\n",
    paste(model_lines, collapse = "\n"), "\n"
  )
  
  structure(tab, heading = heading, class = c("anova", "data.frame"))
}

anova.hurdle <- function(object, ..., test = c("Chisq", "none", "LRT")) {
  if (is.null(test)) test <- "Chisq"
  else if (isFALSE(test)) test <- "none"
  else test <- match.arg(test)
  if (test == "LRT") test <- "Chisq"
  models <- list(object, ...)
  if (length(models) < 2L) {
    stop("anova.hurdle() currently supports only multi-model comparisons: anova(m1, m2, ...).")
  }

  # Basic class check
  bad <- !vapply(models, inherits, logical(1), what = "hurdle")
  if (any(bad)) stop("All models must be of class 'hurdle'.")

  responses <- vapply(models, function(x) paste(deparse(formula(x)[[2L]]), collapse = "\n"), character(1))
  if (!all(responses == responses[1L])) {
    stop("Models must have the same response.")
  }

  # Extract logLik objects (pscl provides logLik.hurdle with df and nobs attrs)
  ll <- lapply(models, logLik)
  df <- vapply(ll, function(x) as.integer(attr(x, "df")), integer(1))
  n  <- vapply(ll, function(x) as.integer(attr(x, "nobs")), integer(1))

  if (length(unique(n)) != 1L) {
    stop("Models were not fit to the same number of observations (nobs differs).")
  }

  # Order by increasing df (small -> large)
  ord <- order(df)
  if (!identical(ord, seq_along(models))) {
    warning("Models reordered by increasing number of parameters for comparison.")
    models <- models[ord]
    ll <- ll[ord]
    df <- df[ord]
  }

  # Check nestedness: smaller model terms must be a subset of larger model terms
  for (i in seq_along(models)[-1]) {
    small <- models[[i - 1]]
    large <- models[[i]]

    count_nested <- all(attr(small$terms$count, "term.labels") %in%
                          attr(large$terms$count, "term.labels"))
    zero_nested  <- all(attr(small$terms$zero, "term.labels") %in%
                          attr(large$terms$zero, "term.labels"))

    if (!count_nested || !zero_nested) {
      warning("Models ", i - 1, " and ", i,
              " may not be nested; LRT results may be invalid.")
    }
  }

  loglik <- vapply(ll, as.numeric, numeric(1))

  # Residual df = nobs - #params; deviance = -2 * logLik
  nobs    <- n[1]
  resdf   <- nobs - df
  resdev  <- -2 * loglik

  # Successive differences (NA for first row, like anova.glm)
  df_diff  <- c(NA_integer_, diff(df))
  dev_diff <- c(NA_real_, -diff(resdev))  # reduction in deviance (positive = improvement)

  pval <- rep(NA_real_, length(models))
  if (test == "Chisq") {
    for (i in seq_along(models)[-1]) {
      if (df_diff[i] > 0 && dev_diff[i] >= 0) {
        pval[i] <- pchisq(dev_diff[i], df = df_diff[i], lower.tail = FALSE)
      } else if (!is.na(dev_diff[i]) && dev_diff[i] < 0) {
        warning("Negative LRT statistic for models ", i - 1, " vs ", i,
                "; models may not be properly nested.")
      }
    }
  }

  tab <- data.frame(
    `Resid. Df`  = resdf,
    `Resid. Dev` = resdev,
    Df           = df_diff,
    Deviance     = dev_diff,
    `Pr(>Chi)`   = pval,
    check.names  = FALSE
  )
  rownames(tab) <- seq_along(models)

  # Build heading with model formulas, expanding "." to actual terms
  formulas <- vapply(models, function(m) {
    resp <- all.vars(m$formula)[1]
    count_terms <- attr(m$terms$count, "term.labels")
    zero_terms  <- attr(m$terms$zero, "term.labels")
    count_rhs <- if (length(count_terms)) paste(count_terms, collapse = " + ") else "1"
    zero_rhs  <- if (length(zero_terms)) paste(zero_terms, collapse = " + ") else "1"
    paste0(resp, " ~ ", count_rhs, " | ", zero_rhs)
  }, character(1))
  model_lines <- paste0("Model ", seq_along(models), ": ", formulas)
  heading <- paste0(
    "Likelihood Ratio Tests of Hurdle Count Regression Models\n\n",
    paste(model_lines, collapse = "\n"), "\n"
  )

  structure(tab, heading = heading, class = c("anova", "data.frame"))
}
