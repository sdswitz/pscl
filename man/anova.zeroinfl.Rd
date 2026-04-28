\name{anova.zeroinfl}
\alias{anova.zeroinfl}
\alias{anova.hurdle}

\title{Likelihood Ratio Tests for Zero-Inflated and Hurdle Count Regression Models}

\description{
  Compute likelihood ratio test tables comparing two or more nested
  zero-inflated (\code{\link{zeroinfl}}) or hurdle (\code{\link{hurdle}})
  count regression models.
}

\usage{
\method{anova}{zeroinfl}(object, \dots, test = c("Chisq", "none", "LRT"))
\method{anova}{hurdle}(object, \dots, test = c("Chisq", "none", "LRT"))
}

\arguments{
  \item{object, \dots}{objects of class \code{"zeroinfl"} (for
    \code{anova.zeroinfl}) or \code{"hurdle"} (for \code{anova.hurdle}),
    typically the result of calls to \code{\link{zeroinfl}} or
    \code{\link{hurdle}}. At least two models must be supplied.}
  \item{test}{a character string specifying the test statistic to be
    computed. Either \code{"Chisq"} (the default) or \code{"LRT"} for a
    likelihood ratio chi-squared test, or \code{"none"} or \code{FALSE}
    to suppress the test. \code{NULL} is treated as the default.}
}

\details{
  Unlike \code{\link[stats]{anova.glm}}, these methods only support
  multi-model comparisons: at least two fitted model objects must be
  supplied. A sequential analysis of deviance table for a single fit is
  not produced because zero-inflated and hurdle models have two linear
  predictors (count and zero components), and there is no single
  canonical term ordering.

  When multiple models are supplied, the table has a row for the
  residual degrees of freedom and deviance (\eqn{-2 \log L}) for each
  model. For all but the first model, the change in degrees of freedom
  and deviance is also reported. The models are automatically reordered
  by increasing number of estimated parameters; a warning is issued if
  reordering occurs.

  This only makes statistical sense if the models are nested, i.e., if
  the terms of the smaller model's count and zero-inflation components
  are each a subset of those of the larger model. A warning is issued
  if the supplied models do not appear to be nested in this sense, or
  if the resulting likelihood ratio statistic is negative.

  All models must have the same response and be fit to the same number
  of observations; an error is raised otherwise. This may be a problem
  if there are missing values and R's default of \code{na.action =
  na.omit} is used.

  With \code{test = "Chisq"} or \code{test = "LRT"}, p-values are computed from the
  chi-squared distribution with degrees of freedom equal to the
  difference in the number of estimated parameters between successive
  models.
}

\value{
  An object of class \code{"anova"} inheriting from class
  \code{"data.frame"}, with columns \code{"Resid. Df"},
  \code{"Resid. Dev"}, \code{"Df"}, \code{"Deviance"}, and (if
  \code{test = "Chisq"}) \code{"Pr(>Chi)"}.
}

\section{Warning}{
  The comparison between two or more models will only be valid if they
  are fitted to the same dataset and are properly nested. Nestedness is
  checked only at the level of term labels in the count and
  zero-inflation components; more subtle forms of non-nestedness are
  not detected.
}

\references{
Cameron, A. Colin and Pravin K. Trivedi. 1998. \emph{Regression Analysis of Count
Data.} New York: Cambridge University Press.

Zeileis, Achim, Christian Kleiber and Simon Jackman 2008.
\dQuote{Regression Models for Count Data in R.}
\emph{Journal of Statistical Software}, \bold{27}(8).
URL \url{https://www.jstatsoft.org/v27/i08/}.
}

\seealso{\code{\link{zeroinfl}}, \code{\link{hurdle}},
  \code{\link[stats]{anova.glm}}, \code{\link[stats]{logLik}}
}

\examples{
## data
data("bioChemists", package = "pscl")

## zero-inflated Poisson models
fm_zip0 <- zeroinfl(art ~ 1     | 1, data = bioChemists)
fm_zip1 <- zeroinfl(art ~ fem   | 1, data = bioChemists)
fm_zip2 <- zeroinfl(art ~ .     | 1, data = bioChemists)
fm_zip3 <- zeroinfl(art ~ .     | ., data = bioChemists)

anova(fm_zip0, fm_zip1, fm_zip2, fm_zip3)

## hurdle models
fm_h0 <- hurdle(art ~ 1   | 1, data = bioChemists)
fm_h1 <- hurdle(art ~ fem | 1, data = bioChemists)
fm_h2 <- hurdle(art ~ .   | ., data = bioChemists)

anova(fm_h0, fm_h1, fm_h2)

## suppress p-values
anova(fm_h0, fm_h1, fm_h2, test = "none")
}

\keyword{regression}
\keyword{htest}
