\name{serum_pos_fullMS}
\alias{serum_pos_fullMS}
\docType{data}
\title{
Peak list example of MS1 obtained from a pool of serum samples.
}
\description{
Data frame containing the MS1 peak list of a pool of serum samples.
}
\usage{data("serum_pos_fullMS")}
\format{
  A data frame with 1177 observations on the following 3 variables.
  \describe{
    \item{\code{m.z}}{a numeric vector}
    \item{\code{RT}}{a numeric vector}
    \item{\code{int}}{a numeric vector}
  }
}
\details{
Q-TOF 6550 from Agilent.
MS/MS DIA acquisition.
ESI +.
Collision energy: 0.
Preprocessed with enviPick(peak-picking) and xcms/CAMERA (isotopes removal).
}
\keyword{datasets}
