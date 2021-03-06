% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.gbmars}
\alias{predict.gbmars}
\title{Predict function for class "gbmars" (Depreciated)}
\usage{
\method{predict}{gbmars}(model, Xnew, burn = 1, thin = 1)
}
\arguments{
\item{model}{an object of class "gbmars" created by the gbmars() function.}

\item{Xnew}{a matrix of predictor variables with ncol(Xnew) == ncol(X)}

\item{burn}{burn-in}

\item{thin}{thinning}
}
\description{
This function generates posterior predictions using predictors Xnew and a model of class "gbmars"
}
\details{
add stuff here
}
\examples{
foo <- 1 + 1

}
