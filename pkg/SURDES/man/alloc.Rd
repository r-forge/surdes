\name{alloc}
\Rdversion{1.1}
\alias{alloc}

\title{Site allocation procedure using one matrix  
}
\description{A stratified survey design that selects a set of sampling 
  localities from a universe of available sites using an allocation procedure in
  which environmental and geographical distances are assumed to be surrogates 
  for diversity variations. Environmental and geographical distances are 
  combined prior to the calculations.

  Selects the sampling points with a set of iterative rules. First step: 
  Maximizises both the amount of spatio-environmental coverage using a p-median 
  allocation procedure. Next steps: uses a set of rules or conditions defined by
  the user. Conditions should be related to the prioritization in the selection
  procedure (for example: conservation status or distance to roads). For each
  rule a vector of values and the type of criteria should be defined (see 
  definition of criteria and conditions). 
}

\usage{
  alloc(mdist, vini=rep(0, dim(mdist)[2]), vtarget=rep(1, dim(mdist)[2]), 
  	criteria, sdint=rep(1, length(criteria)), conditions, iter)
}

\arguments{

  \item{mdist}{
    A dissimilarity matrix coming usually from the combination of a matrix of 
    environmental distances and a matrix of geographic distances. The maximum 
    value of the matrices need to be approximately equal, so if one of the 
    matrices has very high values, usually the geographic distances matrix, 
    normalization is recommended to restrict the variation within 0 and 1. Once 
    the matrices are normalized multiplying the two matrices is enough to get an
    even representation of selected sample sites. Missing values are not 
    allowed.
  }

  \item{vini}{
    A binary vector (0,1) with the same length as the number of columns in 
    \code{mdist}. Should be a set of previously well surveyed localities, these
    points will be used as starting points for the iteration procedure.

    If \code{vini} is not provided the function will create a vector with all
    values are equal to "0". Including initial values is recommended, if there
    are no previously well known sites the best option is to create a random
    vector with one or two positive values. The starting points or previously
    well known places can affect the result of the selection. Caution is
    required when the initial vector is not provided.
  }

  \item{vtarget}{
    A binary vector where values are equal to "1" are the sample sites
    fulfilling a set of minimum conditions defined by the user. For example,
    locations with good conservation status or an area above a threshold level.

    The length of the vector has to be the same as the number of columns in
    \code{mdist}.  If \code{vtarget} is not provided by the user the function
    will compute a vector with all values equal to "1" values, and therefore
    all the sites are assumed to be adequate in the selection procedure.
  }

  \item{criteria}{
    A vector of character strings specifying the criteria that will be applied 
    to the conditions vector. Currently available options are \code{minimum}, 
    \code{maximum}, \code{rangemax} and \code{rangemin}.

    \code{minimum} selects sample sites with the minimum value within the conditions 
    vector and is intended to be used when the variable in the conditions vector
    is qualitative (i.e. ordinal).

    \code{maximum} selects sample sites with the maximum value within the conditions 
    vector and is intended to be used when the variable in the conditions vector
    is qualitative (i.e. ordinal).

    \code{rangemin} selects the sample sites with a value between the minimum
    (min) and one standard deviation (min + sd) within the conditions vector
    and is intended to be used when the variable in conditions vector is
    quantitative. The function allows changes in the number of standard
    deviations added to the minimum value. By default the function will add 1
    standard deviation, but it can be changed to add 0.5, 2, or other user
    defined value (see item "n").

    \code{rangemax} selects the sample sites with a value between the maximum
    (max) and one standard deviation (max - sd) of the values in the conditions
    vector. It is intended to be used when the variable in conditions vector is
    quantitative. The function allows changes in the number of standard
    deviations subtracted to the minimum value. By default the package will
    subtract one standard deviation, but it can be changed to subtract 0.5, 2
    or other user defined value (see item "n").
    
    For example, if the first condition (the first row in the conditions matrix)
    is conservation status in a quantitative scale from 1 to 6 and we are 
    interested in well conserved sites we should chose \code{maximum} as the
    first criterion. 

    The first criterion in the vector \code{criteria} will be applied to the 
    condition in the first row in \code{conditions} vector and so on.
  }
  
  \item{sdint}{
    A vector of numerical values with length equal to the number of criteria. 

    By default the function defines a vector with all values equal to "1". So
    one standard deviation will added or subtracted when applying the
    continuous criteria (\code{rangemax} or \code{rangemin}). 

    To change the number of standard deviations added (in the case of
    \code{rangemin} criteria) or subtracted (in the case of \code{rangemax}
    criteria) the values of the vector should be changed.

    The value for ordinal criteria in the \code{sdint} vector needs to be equal
    to "1".  In an example with one ordinal criterion, and two continuous
    criteria, if the user needs the third criteria to add 2 standard
    deviations, the vector should look like this: 

    \code{sdint <- c(1, 1, 2)}

    Note that the first "1" in the vector applies to the ordinal criteria. 

    See Example 2 in the examples Section for further explanations. 

    Changes in the values of \code{sdint} are recommended when the variability
    in the data in conditions is very high or very low.   
  }

  \item{conditions}{
    A matrix with the same number of columns as \code{mdist}.

    Rules or conditions are defined by the user. A typical example is the level
    of conservation of each site.
  }

  \item{iter}{
    Number of iterations.
        
    Never bigger than the number of available sites, which are all the sites
    minus the previously well sampled sites (\code{vini=1}) and the sites that
    are not adequate (\code{vtarget=0}).
    
    Note that the procedure selects one site per iteration, so the number of
    iterations and the number of selected points will be the same.
  }
}

\details{
  In each iteration, the procedure first selects a set of sample sites
  minimizing the total distance between the non selected sites (p-median
  procedure), and then, in the next steps, criteria defined by the user in
  \code{conditions} are used iteratively to restrict the number of selected
  sites. If more than one point remain selected after applying all the
  conditions, the function selects one at random. Therefore, in each step just
  one sample site is selected. The next iteration will start with all the
  points minus the selected point in the previous iteration and so on.
}

\value{
  The function prints the criterion used to select the sampling point in each
  iteration. By default the names of the criteria are\code{crit1},
  \code{crit2}, etc.

  \code{crit1} represents the first criterion included in the \code{criteria}
  vector, \code{crit2} the second and so on.

  The function also returns a list including two matrices: \code{pmmatrix},
  \code{selmatrix} (see below).
  
  \item{selmatrix}{
    \code{selmatrix} is a matrix of the selected sampling points after applying
    the allocation procedure, rows being each one of the iterations and columns
    the sampling sites in the same order as entered by the user. 
  }

  \item{pmmatrix}{
    \code{pmmatrix} is a matrix of the distances between the selected sites and
    the non-selected sites in each step, rows being each one of the iterations
    and columns the sampling sites in the same order as entered by the user. 

    The value of each cell is the distance between each site to the closest
    previously selected site. Note that selected sites in each step have "0"
    value in this matrix. 

    Adding the values of the rows (see function \code{rowSum()} and examples
    below) in this matrix gives the amount of uncovered variability in each
    step. It is very useful to represent the decrease of uncovered variability
    in each iteration. See examples and Hortal et al. (2005) for an extended
    explanation.
  }
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

\author{
  Nagore Garcia Medina & Bernardo Garcia Carreras
}

\examples{
  # load the environmental and spatial matrices
  data(env)
  data(geogr)

  # multiply environmental and geographic matrices
  mdist <- env * geogr
       
  # load a matrix containing the conditions, the initial vector and the target 
  # points
  data(conditions)
  
  # define vector of initial points  
  data(vini)

  # define the criteria to apply to the conditions: maximise conservation
  # (qualitative variable), maximise area of the site (qualitative variable)
  # and minimise distance to roads (quantitative variable)
  data(criteria)

  result <- alloc(mdist=mdist, vini=vini, criteria=criteria, 
  		conditions=conditions, iter=20)

  uncov <- rowSums(result$pmmatrix)
  plot(1:length(uncov), uncov)

  # Example 2 
  # With the same data as Example 1, but with a change in the number of
  # standard deviations added in the first criteria, note that the first
  # condition is slope of the species area curve, and the criteria to apply is
  # \code{rangemax}, in this case the range of selection for the condition will
  # be = max-0.1*sd
  result <- alloc(mdist=mdist, vini=vini, criteria=criteria, sdint=c(0.1,1,1), 
  		conditions=conditions, iter=20)

  uncov <- rowSums(result$pmmatrix)
  plot(1:length(uncov), uncov)
}

\keyword{ function }
\keyword{ p-median }
\keyword{ survey }
\keyword{ design }

