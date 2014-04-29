\name{epi.design}
\alias{epi.design}
\concept{model fitting}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create design matrix for two-locus GLM models
}
\description{
The two-locus GLM model is defined based on data for two unlinked biallelic loci. In order to fit the model, it is necessary to convert the data for the two loci into appropriate terms for the nine corresponding model parameters. This function creates that \eqn{N \times 9}{N by 9} design matrix from the data.
}
\usage{
epi.design(x, z, model=NULL, pa=NULL, pb=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the numeric vector of data at the first locus, coded 0,1,2 for the number of reference alleles. NAs are allowed.
}
  \item{z}{
data for the second locus, as \code{x}
}
\item{model}{
One of "NOIA_S", "NOIA_F", "G2A", "F2", "Finf", "unweight", or "genotype" specifying the desired regression model. See \code{\link{epipen-package}} for details of these models.
}
\item{pa}{
(optional) true allele or genotype frequencies for the first locus. Used for models "NOIA_S", "NOIA_F", and "G2A". Can be specified as either a single numeric value (0<pa<1) for the allele frequency, or a vector of three values giving the genotype frequencies for the three possible biallelic genotypes (should sum to 1). If the allele frequency is specified, Hardy-Weinberg equilibrium is assumed. If omitted, values are estimated from the data.
}
\item{pb}{
(optional) true allele or genotype frequencies for the second locus, as pa above
}
}
\value{
An \eqn{N \times 9}{N by 9} design matrix, with columns corresponding to the parameters of the two-locus GLM model, where \eqn{N} is the length of the data vectors \code{x} and \code{z}. Specifically, the columns of the output contain the appropriate dummy variables for the selected \code{model}. Columns of the matrix are labeled with the appropriate variable names.
}
\references{
Walters, R.K., Laurin, C., and Lubke, G.H. (2014). EpiPen: An R package to investigate two-locus epistatic models. \emph{Twin Research and Human Genetics}. doi:10.1017/thg.2014.25
}
\author{
Raymond Walters
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{epi.data}} to generate SNP and phenotype data; \code{\link{epipen-package}} for details on available GLM models
}
\examples{
# Generate random data for 2 loci
x <- rbinom(100,2,.5)
z <- rbinom(100,2,.3)

# Create design matrix
dat.M <- epi.design(x,z,model="NOIA_S")
head(dat.M)
dim(dat.M)

# fit model using design matrix
u <- .2*x + .2*(x*z)
y <- rbinom(100,1,u)
fit <- glm(y~.,data=cbind(y=y,dat.M),family=binomial("logit"))
summary(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
