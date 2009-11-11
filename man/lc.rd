\name{lc}
\alias{lc}
\alias{lc.0}
\alias{lc.1}

\title{Templates for pseudo-species according to Turboveg layer informations.}

\description{Differentiate species according to layer (tree, shrub, juvenile etc.) or other species-plot informations from Turboveg. 
The data.frames lc.0 and lc.1 are templates for layer aggregation/differentiation.}

\details{
Column \code{layer} point to the Turboveg layer specification (see Turboveg Help) and column \code{comb} defines the aggregation.

\code{lc.1} = Default layer combination. Differenciates only tree layers, juveniles and all the rest.
\tabular{rr}{
layer \tab  comb \cr
 0 \tab 0 \cr
 1 \tab 1 \cr
 2 \tab 1 \cr
 3 \tab 1 \cr
 4 \tab 0 \cr
 5 \tab 0 \cr
 6 \tab 0 \cr
 7 \tab 2 \cr
 8 \tab 2 \cr
 9 \tab 0 \cr
}

\code{lc.0} = Use every layer differentiation from 0 to 9 in Turboveg database as pseudo-species.
\tabular{rr}{
layer \tab comb \cr
 0 \tab 0 \cr
 1 \tab 1 \cr
 2 \tab 2 \cr
 3 \tab 3 \cr
 4 \tab 4 \cr
 5 \tab 5 \cr
 6 \tab 6 \cr
 7 \tab 7 \cr
 8 \tab 8 \cr
 9 \tab 9 \cr
}

}

\seealso{\link{tv.veg}
}
\author{Florian Jansen
\email{jansen@uni-greifswald.de}
}

\keyword{data}