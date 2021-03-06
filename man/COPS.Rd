% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COPS.R
\name{COPS}
\alias{COPS}
\title{The COPS procedure}
\usage{
COPS(X, Y, h_seq = NULL, n_ss = 20, kappa = 1)
}
\arguments{
\item{X}{design points, should be standardized in (0, 1)}

\item{Y}{responses}

\item{h_seq}{candidates for the bandwidth}

\item{n_ss}{times of sample-splitting for bandwidth selection}

\item{kappa}{controlling the size of removed neighborhood of detected jumps}
}
\value{
a list contains the number, locations and sizes of detected jumps
}
\description{
The COPS procedure
}
\examples{
NA
}
