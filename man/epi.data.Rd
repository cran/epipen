\name{epi.data}
\alias{epi.data}
\concept{data generation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate data under two-locus GLM model
}
\description{
Generates data following the two-locus GLM model with the specified model parameters. The phenotype may be generated based on existing genotypes, or genotype data can be generated.
}
\usage{
epi.data(betas, model, distn, x = NULL, z = NULL,
        n = NULL, ncase = NULL, ncontrol = NULL,
        pa = NULL, pb = NULL, ve = NA, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betas}{
The 9 regression parameters for the named \code{model}. See \code{\link{epipen-package}} for details of the parameters of the available models. 
}
\item{model}{
The model used to define the regression parameters. Can be one of: "NOIA_S", "NOIA_F", "G2A", "F2", "Finf", "unweight", or "genotype". The details on each of these models are provided in \code{\link{epipen-package}}.
}
\item{distn}{
either "normal", "logit", or "probit"; specifying the appropriate distribution for the GLM model. If "normal", then the phenotype is continuous with normally distributed errors. If "logit", then the phenotype is dichotomous based on an underlying continuous variable with logistically distributed errors. If "probit", then the phenotype is dichotomous based on an underlying continuous variable with normally distributed errors.
}
  \item{x}{
(optional) a numeric vector of genotype data at the first locus, coded 0,1,2,NA for the number of reference alleles. If omitted, genotype data will be randomly generated.
}
  \item{z}{
(optional) data for the second locus, as \code{x}
}
  \item{n}{
sample size for the generated data. Ignored with a warning if data is provided for \code{x} and/or \code{z}.
}
  \item{ncase}{
(optional) number of affected individuals to include in the generated data. Only applicable for \code{distn="logit"} or \code{"probit"}. Ignored with a warning if data is provided for \code{x} and/or \code{z}.
}
  \item{ncontrol}{
(optional) number of unaffected individuals to include in the generated data. Only applicable for \code{distn="logit"} or \code{"probit"}. Ignored with a warning if data is provided for \code{x} and/or \code{z}.
}
  \item{pa}{
allele or genotype frequencies for the first locus. Can be specified as either a single numeric value (0<pa<1) for the allele frequency, or a vector of three values giving the genotype frequencies for the three possible biallelic genotypes (should sum to 1). If the allele frequency is specified, Hardy-Weinberg equilibrium is assumed. Ignored with a warning if \code{x} is given.
}
  \item{pb}{
allele or genotype frequencies for the second locus, as \code{pa}. Ignored with a warning if \code{z} is given.
}
  \item{ve}{
the variance of the error term in the normal model. If not specified, it will be set so that the total phenotypic variance is 1. If \code{distn="logit"} or \code{"probit"}, this argument is ignored (with a warning).
}
  \item{verbose}{
Logical, indicating whether to provide additional information about the model. Default is FALSE.
}
}
\details{
Data may be generated based on existing genotype data (coded 0,1,2 or NA), or new genotype data may be generated from \eqn{Bin(2,p)}. If existing genotype data is provided for \code{x} and \code{z} then the sample size \code{n} and allele frequencies \code{pa} and \code{pb} are ignored. If no data is provided then the sample size and allele frequencies must be given so that new genotype data can be generated. For the logit or probit model, the sample size can optionally be specified to include a specific number of cases (\code{y=1}) and controls (\code{y=0}), instead of a simple random sample of size \code{n} from the population.

Optionally, by setting \code{verbose=TRUE} expected characteristics of the generated data will be provided in the same format as \code{\link{epi.convert}}. If existing genotype data was used to generate the phenotype, these expected values are computed based on the observed allele frequencies in the input data.

The two loci for the genotype data are assumed to be biallelic and unlinked.
}
\value{
If \code{verbose=FALSE}, the function returns a data frame with columns \code{y}, \code{x}, and \code{z}, where \code{y} is the generated phenotype data and \code{x}, and \code{z} are either the provided genotype data or new generated genotypes.

If \code{verbose=TRUE}, more complete information is returned. The output will include:
\item{Full_Data}{data frame containing the full generated data, including latent values where applicable. Specifically,

\code{y} - the generated phenotype

\code{ystar} - (logit or probit model only) the latent continuous variable dichotomized to \code{y}

\code{XB} - the value of the linear model minus the error term

\code{yprob} - (logit or probit model only) the probability \code{y=1} conditional on the genotypes

\code{x} - genotype data for the first locus

\code{z} - genotype data for the second locus
}
\item{Model_Info}{Expected values and variances for the model. If existing genotype data was used for \code{x} or \code{z}, model information is computed based on the observed allele frequency in the data. If a fixed number of cases and controls were generated, this information reflects the population from which the individuals were sampled, not expected sample statistics. See \code{\link{epi.convert}} for details on the extra provided information.
}
}
\references{
Walters, R.K., Laurin, C., and Lubke, G.H. (in prep). EpiPen: An R package to investigate two-locus epistatic models.
}
\author{
Raymond Walters
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{epipen-package}} for model details; \code{\link{epi.design}} for fitting two-locus models to generated data
}
\examples{
# generate data from existing genotypes
x <- rbinom(100,2,.5)
z <- rbinom(100,2,.3)
betas <- c(0,0,0,0,0,0,.25,0,.5)
dat.M <- epi.data(betas=betas, model="Finf", distn="normal", x=x, z=z)
head(dat.M)

# data for logit model without existing genotypes
# include extra output
dat2.M <- epi.data(betas=betas, model="Finf", distn="logit", 
        n=100, pa=.5, pb=.3, verbose=TRUE)
head(dat2.M)
dat2.M$Model_Info
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
\keyword{ models }
