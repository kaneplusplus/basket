
\name{cluster_baskets}
\alias{cluster_baskets}
\title{Get the Cluster Members of MEM Models}
\usage{
cluster_baskets(x)
}
\arguments{
\item{x}{either an MEM model or basket object.}
}
\value{
A named list is returned where the name is the cluster name and
each element of the list comprises character vector describing the 
baskets assigned to each cluster.
}
\description{
Object returned by the `mem_mcmc()` and `mem_exact()` include
information about the subgroups in the trials and the cluster 
composed of sets of similar subgroups. This function returns the name
of each subgroup in a cluster.
}
\examples{
\dontrun{
# Create an MEM analysis of the Vemurafenib trial data.

data(vemu_wide)

mem_analysis <- mem_exact(vemu_wide$responders, 
                          vemu_wide$evaluable,
                          vemu_wide$baskets)

# Get the subgroups in the clusters.
cluster_baskets(mem_analysis)
}}
