\name{allocTM}
\Rdversion{1.1}
\alias{allocTM}

\title{Site allocation criterion that uses two initial matrices      
}

\description{
  A stratified survey design that selects a set of sampling localities from a
  universe of available sites using an allocation procedure in which
  environmental and geographical distances are assumed to be surrogates for
  diversity variations. Environmental and geographical distances are considered
  separately.

  The procedure selects the sampling points with a set of iterative rules. First 
  step: maximizises the amount of environmental coverage using a p-median 
  allocation procedure (for further explanations see Faith 1996, Hortal and Lobo
  2005, Hortal et al. 2009). Second step: Takes the sampling localities selected
  in the previous step and uses a p-median allocation procedure to maximize the 
  amount of geographical distance covered. Next steps: Takes the sampling 
  localities selected after second step and uses a set of rules or conditions 
  defined by the user to select final localities. Conditions should be related 
  to the prioritisation in the selection procedure, for example: conservation 
  status or distance to roads. For each rule a vector of values and the type of 
  criteria should be defined (see definition of criteria and conditions). 
}

\usage{
  allocTM(dist, m, vini=rep(0, dim(dist)[2]), vtarget=rep(1, dim(dist)[2]), 
  	criteria, sdint=rep(1, length(criteria)), conditions, iter)
}

\arguments{
  \item{dist}{
    Geographical distance matrix, a dissimilarity matrix containing all pairwise
    distances between the available sites in the region usually calculated with 
    \code{daisy()} from package cluster or similar function.
  }

  \item{m}{
    Environmental distance matrix, a dissimilarity matrix containing all 
    pairwise distances between the available sites in the region, usually 
    calculated with \code{daisy()} from package cluster or similar function.
  }
  
  \item{vini}{
    See definition in \code{alloc()}.
  }

  \item{vtarget}{
    See definition in \code{alloc()}.
  }

  \item{criteria}{
    See definition in \code{alloc()}.
  }

  \item{sdint}{
    See definition in \code{alloc()}.
  }

  \item{conditions}{
    See definition in \code{alloc()}.
  }

  \item{iter}{
    See definition in \code{alloc()}.
  }
}

\details{
  In each iteration, the procedure first selects a set of sample sites
  minimizing the environmental distance between the non selected sites
  (p-median procedure), then from the sites selected in the first step selects
  a second set of sites minimizing the geographical distance between the non
  selected sites (p-median procedure), and in the following steps, criteria
  defined by the user in \code{conditions} are used iteratively to restrict the
  number of selected sites. If more than one point remain selected after
  applying all the conditions, the function selects one at random. Therefore,
  in each step just one sample site is selected. The next iteration will start
  with all the points minus the selected point in the previous iteration.
}

\value{
  See definition in \code{alloc()}.
}

\references{
  Church, R.L. & Sorensen, P. (1994) Integrating Normative location models into
  GIS: problems and prospects  with the p-median model. Technical Report,
  NGCIA.

  Church, R.L. (2002) Geographical information systems and location science.
  Computers and Operation Research 29: 541-562.

  Faith, D.P. & Walker, P.A. (1996) Environmental diversity: on the best 
  possible use of surrogate data for assessing the relative biodiversity of 
  sets of areas. Biodiversity and Conservation 5: 399-415.

  Hortal, J., Araujo, M.B. & Lobo, J.M. (2009) Testing the effectiveness of 
  discrete and continuous environmental diversity as a surrogate for species 
  diversity. Ecological Indicators 9: 139-149. 
   
  Hortal, J. & Lobo, J.M. (2005) An ED-based protocol for the optimal sampling
  of diversity. Biodiversity and Conservation 14: 2913-2947.
}

\author{Nagore Garcia Medina & Bernardo Garcia Carreras}

\examples{
  # load the environmental and spatial matrices
  data(env)
  data(geogr)

  # load a matrix containing the conditions, the initial vector and the target
  # points
  data(conditions)
  
  # define vector of initial points  
  data(vini)

  # define the criteria to apply to the conditions: maximise conservation
  # (qualitative variable), maximise area of the site (qualitative variable)
  # and minimize distance to roads (quantitative variable)
  data(criteria)

  result <- allocTM(m=env, dist=geogr, vini=vini, criteria=criteria, 
  	conditions=conditions, iter=20)

  uncov <- rowSums(result$pmmatrix)
  plot(1:length(uncov),uncov)
}

\keyword{ function }  
\keyword{ p-median }
\keyword{ survey }
\keyword{ design }

