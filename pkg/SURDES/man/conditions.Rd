\name{conditions}

\Rdversion{1.1}

\alias{conditions}

\docType{data}

\title{Additional conditions}

\description{
  A matrix including the different conditions (one row for each condition),
  to be applied after the p-median procedure. Depending on the aim of the
  survey we may need to include other aspects different to environment and
  geographic distances. In the example \code{slope} is the slope of the species
  area curve, the higher the slope the better the sampling effort. \code{rare}
  is an index of rarity of the localities and \code{usab} an ordinal variable
  classifying the localities according to the amount of open lanscape, ranging
  from 1 to 6, 6 being the class with the highest amount of open landscape.
}

\usage{data(conditions)}

\format{
  The format is:
  num [1:3, 1:108] 1 16 1 1 17.7 ...
  - attr(*, "dimnames")=List of 2
  ..$ : chr [1:3] "slope" "rare" "usab"
  ..$ : NULL
}

\source{
   Hortal, J. & Lobo, J.M. (2005) An ED-based protocol for the optimal sampling
   of diversity. Biodiversity and Conservation 14: 2913-2947.
}

\keyword{datasets}

