\name{SURDES-package}

\Rdversion{1.1}

\alias{SURDES-package}

\alias{SURDES}

\docType{package}

\title{
  Functions for survey design based on environmental distance (ED) site 
  allocation criterion.
}

\description{
  A method for selecting sampling points in a given territory using an
  environmental distance site allocation criterion which minimises the
  variation uncovered by the selected sites between all available locations
  using a p-median procedure.
}

\details{
  \tabular{ll}{
    Package: \tab SURDES\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2010-02-02\cr
    License: \tab GLU-2\cr
    LazyLoad: \tab yes\cr
  }
  
  Two different functions are included: 

  1- \code{alloc()}

  in which the allocation procedure is based on one matrix combining
  environmental and geographic distance matrices. It is recommended to multiply
  distance and environmental matrices before running the analysis (see
  example).
     

  2- \code{allocTM()}

  which is based on a sequential procedure that first allocates points in the
  environmental distance matrix and then in the geographical distance matrix.
  Usually the function gives less weight to geographical distance than
  \code{alloc()}.
}

\author{Nagore Garcia Medina & Bernardo Garcia Carreras
  Maintainer: Nagore Garcia Medina<ngmedina@gmail.com>
}

\references{
  Church, R.L. & Sorensen, P. (1994) Integrating Normative location models into
  GIS: problems and prospects  with the p-median model. Technical Report,
  NGCIA.

  Church, R.L. (2002) Geographical information systems and location science.
  Computers and Operation Research 29: 541-562.

  Faith, D.P. & Walker, P.A. (1996) Environmental diversity: on the best
  possible use of surrogate data for assessing the relative biodiversity of
  sets of areas.  Biodiversity and Conservation 5: 399-415.

  Hortal, J., Araujo, M.B. & Lobo, J.M. (2009) Testing the effectiveness of
  discrete and continuous environmental diversity as a surrogate for species
  diversity. Ecological Indicators 9: 139-149. 
   
  Hortal, J. & Lobo, J.M. (2005) An ED-based protocol for the optimal sampling
  of diversity. Biodiversity and Conservation 14: 2913-2947.
}

\keyword{ package }
\keyword{ p-median }
\keyword{ allocation }
\keyword{ procedure }

\examples{
  # see examples in the help files of the main functions: \code{alloc()} and
  # \code{allocTM()}.
}

