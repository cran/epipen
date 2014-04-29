\name{epi.convert}
\alias{epi.convert}
\concept{convert models}
\concept{variance decomposition}
\title{
Conversions between models for two-locus interactions
}
\description{
Two-locus epistatic models can be parametrized in a number of ways, including GLM regression models, penetrance or log odds models, or variance component models. This function provides easy conversion between these models in order to allow comparison of these different parametrizations.
}
\usage{
epi.convert(modelin, params, modelout = c(modelin,"NOIA_S"),
        pa = NULL, pb = NULL, distn, mu = NA, ve = NA, 
        signs = rep(1,8), totaleff = NA, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelin}{
The model being converted from. Can be one of "penetrance", "logodds", "varex", "NOIA_S", "NOIA_F", "G2A", "F2", "Finf", "unweight", or "genotype". The details on each of these models are provided in \code{\link{epipen-package}}.
}
  \item{params}{
Complete parameters for the saturated model named in \code{modelin}. For \code{modelin} = \code{"penetrance"} or \code{"logodds"}, should be a 3x3 numeric matrix containing the conditional penetrance or log odds for each pair of genotypes, with the columns corresponding to the genotype at locus A and the rows corresponding to locus B. For \code{modelin="varex"}, should be a numeric vector containing the variance explained by each of the 8 possible variance components. For all other models, \code{param} should be the vector of 9 regression coefficients for the saturated model. See \code{\link{epipen-package}} for more details on each model.
}
  \item{modelout}{
A vector naming the model(s) to be converted to. Can include one or more of "penetrance", "logodds", "varex", "NOIA_S", "NOIA_F", "G2A", "F2", "Finf", "unweight", or "genotype".  The details on each of these models are provided in \code{\link{epipen-package}}.
}
  \item{pa}{
(optional) allele or genotype frequencies for the first locus. Can be specified as either a single numeric value (0<pa<1) for the allele frequency, or a vector of three values giving the genotype frequencies for the three possible biallelic genotypes (should sum to 1). If the allele frequency is specified, Hardy-Weinberg equilibrium is assumed. If omitted, defaults to allele frequency of .5 with a warning.
}
  \item{pb}{
(optional) allele or genotype frequency for the second locus. See \code{pa}.
}
  \item{distn}{
The distribution used for the link function for GLM models. For a continuous phenotype, should be "normal". For dichotomous phenotypes, can be either "logit" or "probit". 
}
  \item{mu}{
(optional) Value of the mean parameter under the NOIA statistical model. Only used for \code{modelin="varex"}, due to the mean being undetermined for the variance component model. If omitted, default value of \code{mu} is selected to set the expected value of the full linear model to zero.
}
  \item{ve}{
(optional) Error variance for the linear model. Only used for either \code{modelin} or \code{modelout="varex"}, and only when \code{distn="normal"}. If omitted for those models, default value of \code{ve} is selected so that the expected total variance of the phenotype equals 1.
}
  \item{signs}{
(optional) Sign of the square root to use for each variance component in \code{params} for \code{modelin="varex"}. Necessary because variance components do not identify the sign of each effect. Default is to use all positive roots.
}
  \item{totaleff}{
(optional) Total proportion of variance explained, including the effects of the two modelled loci as well as any other covariates. Only used for for \code{modelin="varex"}. Useful for computing model parameters that account for other covariates, or to model phenotypes whose variance is not equal to 1. If omitted, default is \code{totaleff = sum(params)}.
}
  \item{verbose}{
Logical value, indicating whether to include extra information on the coverted model, such as marginal and conditional expected values (see output values). If \code{FALSE}, only the model parameters for each \code{modelout} are returned.
}
}
\details{
The input model for two-locus interactions (\code{modelin}) is converted to parameters for one or more alternative models (\code{modelout}). Details on the available models are provided in the package overview (\code{\link{epipen-package}}). The two loci are assumed to be biallelic and unlinked, and are assumed to be in Hardy-Weinberg equilibrium unless otherwise specified with \code{pa} and \code{pb}.

The provided allele or genotype frequencies affect parameters for variance explained (\code{"varex"}), NOIA statistical (\code{"NOIA_S"}), NOIA functional (\code{NOIA_F}) and G2A models, as well as marginal expected values and variances. Other models do not depend on allele frequencies. 

For dichotomous outcomes (\code{distn="logit"} or \code{"probit"}), the variance decomposition is based on the variance of the underlying latent variable \eqn{y^*}{y*} in the GLM framework. This usage is consistent with the pseudo-R^2 for logistic regression introduced by McKelvey & Zavoina (1975). Specifically, the variable \eqn{y^*}{y*} is assumed to be continuous with errors following either a standard logistic distribution (\code{distn="logit"}) or a standard normal distribution (\code{distn="probit"}); the dichotomous phenotype \eqn{y} is created by thresholding this variable.
}
\value{
A list of models, with each list item containing the vector of coefficient values for the specified \code{modelout}. If \code{verbose=TRUE}, then one additional list item is added with the following elements:
\item{Expected_Values}{a list with elements:

\code{Conditional} - a 3x3 matrix of the expected value of the phenotype by genotype (\code{distn="normal"} only)

\code{Penetrances} - a 3x3 matrix of the expected probability of the phenotype by genotype (\code{distn="logit"} or \code{"probit"} only)

\code{MarginA} - length 3 vector of the marginal expected value conditional on the genotype of the first locus

\code{MarginB} - length 3 vector of the marginal expected value conditional on the genotype of the second locus

\code{Unconditional} - overall expected value of the phenotype (\code{distn="normal"} only)

\code{Prevalence} - overall expected probability of the phenotype (\code{distn="logit"} or \code{"probit"} only)
}
\item{Var_y}{expected variance of the phenotype (\code{distn="normal"} only)}

\item{Var_ystar}{expected variance of underlying latent variable \eqn{y^*}{y*} for the phenotype (\code{distn="logit"} or \code{"probit"} only)}

\item{error_variance}{variance of the error term for the phenotype}
}
\references{

Neuman, R.J., and Rice, J.P. (1992). Two-locus models of disease. \emph{Genetic Epidemiology}, 9(5): 347-365.

McKelvey, R., and Zavoina, W. (1975). A Statistical model for the analysis of ordinal level dependent variables. \emph{Journal of Mathematical Sociology}, 4: 103-120.

Walters, R.K., Laurin, C., and Lubke, G.H. (2014). EpiPen: An R package to investigate two-locus epistatic models. \emph{Twin Research and Human Genetics}. doi:10.1017/thg.2014.25
}
\author{
Raymond Walters
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{epipen-package}} for details on the included models
}
\examples{
### convert penetrance model
# using "modifying effect" model (Neuman & Rice, 1992)
pen.M <- matrix(c(.01,.01,.01,
                .01,.01,.02,
                .02,.02,.02),3,3,byrow=TRUE)
epi.convert(modelin="penetrance", params=pen.M, 
        modelout = c("NOIA_S","NOIA_F"),
        pa = .5, pb = .3, distn="logit")

### get variance components from GLM parameters
betas <- c(1,0,0,0,0,0,.25,0,.5)
epi.convert(modelin="Finf", params=betas, modelout="varex", 
        pa=.5, pb=.3, distn = "normal")

### get model from variance components
# additive effect of locus 1, 
# additive x additive interaction, and 
# dominant x dominant interaction
varex <- c(.003,0,0,0,.003,0,0,.003)
epi.convert(modelin="varex", params=varex, modelout="NOIA_S", 
        pa=.5, pb=.3, distn = "normal")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
