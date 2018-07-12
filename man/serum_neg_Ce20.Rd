\name{serum_neg_Ce20}
\alias{serum_neg_Ce20}
\docType{data}
\title{
Peak list example of MSMS1 obtained from a pool of serum samples.
}
\description{
Data frame containing the MSMS1 peak list of a pool of serum samples.
}
\usage{data("serum_neg_Ce20")}
\format{
  A data frame with 38174 observations on the following 3 variables.
  \describe{
    \item{\code{m.z}}{a numeric vector}
    \item{\code{RT}}{a numeric vector}
    \item{\code{int}}{a numeric vector}
  }
}
\details{
Q-TOF 6550 from Agilent.
MS/MS DIA acquisition.
ESI -.
Collision energy: 20.
Preprocessed with enviPick(peak-picking) and xcms/CAMERA (isotopes removal).
}
\keyword{datasets}
