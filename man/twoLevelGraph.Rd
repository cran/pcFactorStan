% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate.R
\name{twoLevelGraph}
\alias{twoLevelGraph}
\title{Create an edge list with a random two level connectivity}
\usage{
twoLevelGraph(name, N, shape1 = 0.8, shape2 = 0.5)
}
\arguments{
\item{name}{vector of vertex names}

\item{N}{number of comparisons}

\item{shape1}{beta distribution parameter for first edge}

\item{shape2}{beta distribution parameter for second edge}
}
\value{
An undirected graph represented as a data frame with each row describing an edge.
}
\description{
Initially, edges are added from the first vertex to all the other
vertices. Thereafter, the first vertex is drawn from a Beta(shape1,
1.0) distribution and the second vertex is drawn from a
Beta(shape2, 1.0) distribution. The idea is that the edges will
tend to connect a small subset of vertices from the top of the tree
to leaf vertices. These vertex connections are similar to the pairs
that you might observe in an elimination tournament. The selected
vertices are sorted so it doesn't matter whether \code{shape1 >
shape2} or \code{shape1 < shape2}.
}
\examples{
twoLevelGraph(letters[1:5], 20)
}
\seealso{
Other graph generators: 
\code{\link{roundRobinGraph}()}
}
\concept{graph generators}
