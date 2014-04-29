\name{epi.power}
\alias{epi.power}
\concept{power analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power analysis for two-locus GLM models
}
\description{
Performs power analysis for one or more parameters in a two-locus epistatic model. This power analysis is based on either Cohen's \eqn{f^2} multiple regression effect size (for continuous outcomes), or the asymptotic power of likelihood ratio tests (for dichotomous outcomes).
}
\usage{
epi.power(truebetas, model, distn, 
        full = rep(TRUE,9), reduced, N, 
        alpha = .05, pa = NULL, pb = NULL, 
        ve = NA, plotit = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{truebetas}{
Numeric vector of the 9 regression parameters for the saturated GLM model named in \code{model}.  See \code{\link{epipen-package}} for details of the parameters of the available models.
}
  \item{model}{
The model used to define the regression parameters. Can be one of: "NOIA_S", "NOIA_F", "G2A", "F2", "Finf", "unweight", or "genotype". The details on each of these models are provided in \code{\link{epipen-package}}.
}
  \item{distn}{
either "normal", "logit", or "probit"; specifying the appropriate distribution for the GLM model. If "normal", then the phenotype is continuous with normally distributed errors. If "logit", then the phenotype is dichotomous based on an underlying continuous variable with logistically distributed errors. If "probit", then the phenotype is dichotomous based on an underlying continuous variable with normally distributed errors.
}
  \item{full}{
Logical vector indicating which of the 9 parameters (from \code{truebetas}) are included in the full model for the desired hypothesis test. Values of TRUE indicate the parameter is included in the full model, and FALSE indicates the parameter is not included.
}
  \item{reduced}{
Logical vector indicating which of the 9 parameters are included in the reduced model for the desired hypothesis test. The included parameters must be a subset of \code{full}. 
}
  \item{N}{
a numeric vector of which sample sizes to use for estimating power
}
  \item{alpha}{
significance level for testing power
}
  \item{pa}{
(optional) allele or genotype frequencies for the first locus. Can be specified as either a single numeric value (0<pa<1) for the allele frequency, or a vector of three values giving the genotype frequencies for the three possible biallelic genotypes (should sum to 1). If the allele frequency is specified, Hardy-Weinberg equilibrium is assumed. If omitted, defaults to allele frequency of .5 with a warning.
}
  \item{pb}{
(optional) allele or genotype frequency for the second locus. See \code{pa}.
}
  \item{ve}{
variance of the error term for the normally distributed model
}
  \item{plotit}{
logical, if \code{TRUE} a plot of the estimated power at each sample size will be produced
}
}
\details{
Power is computed for testing the difference between the \code{full} and \code{reduced} models based on a same of size \code{N} drawn from the population defined by the \code{truebetas} and \code{model}. 

For continuous outcomes (\code{distn="normal"}), the power is computed based on Cohen's \eqn{f^2} effect size for the difference in variance explained by the two models (see \code{\link{pwr.f2.test}}. For dichotomous outcomes (\code{distn="logit"} or \code{"probit"}), the asymptotic power of the likelihood ratio test is used instead (see \code{\link[asypow]{asypow.power}}).

Allele/genotype frequencies \code{pa} and \code{pb} are used to determine the expected joint genotype frequencies for the sample with size \code{N}. This is necessary to get the expected parameter estimates for the reduced model (and the full model if the saturated model is not used). 

}
\value{
A length(N) by 2 matrix of the estimated power at each given sample size. If \code{plotit=TRUE}, a line plot of the power as a function of sample size is also produced.
}
\references{
Walters, R.K., Laurin, C., and Lubke, G.H. (2014). EpiPen: An R package to investigate two-locus epistatic models. \emph{Twin Research and Human Genetics}. doi:10.1017/thg.2014.25
}
\author{
Raymond Walters, Charles Laurin
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
### Example test of interaction using variance explained model
# true effects: additive effect of locus 1, 
# additive x additive interaction, and 
# dominant x dominant interaction
varex <- c(.003,0,0,0,.003,0,0,.003)

# get NOIA_S model parameters for those effects
betas <- epi.convert(modelin="varex", params=varex, modelout="NOIA_S", 
        pa=.5, pb=.3, distn = "normal")

# compute power for 4 df test of all interaction effects
epi.power(truebetas=betas$NOIA_S, model="NOIA_S", distn="normal", 
        reduced=c(rep(TRUE,4),rep(FALSE,5)), 
        N=c(1250,2500,5000,10000,15000,20000,25000,27500,30000), 
        pa=.5, pb=.3, plotit = FALSE)

# compute power for 1 df test of dominant x dominant interaction
epi.power(truebetas=betas$NOIA_S, model="NOIA_S", distn="normal", 
        reduced=c(rep(TRUE,8),rep(FALSE,1)), 
        N=c(1250,2500,5000,10000,15000,20000,25000,27500,30000), 
        pa=.5, pb=.3, plotit = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
