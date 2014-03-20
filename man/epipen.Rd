\name{epipen-package}
\alias{epipen-package}
\alias{epipen}
\docType{package}
\title{
EpiPen: Epistatic and penetrance models for two-locus genetic interactions
}
\description{
A small set of utilities for manipulating two-locus epistatic models. Includes methods for converting between common models (penetrances, log odds, variance components, GLM regression parameters), data generation, and power analysis.
}
\details{
A wide variety of models exist for genetic interactions between two biallelic loci. This package provides a set of utilities for easily working with this array of possible models. Specifically, four basic functions are offered:

\itemize{
\item \code{\link{epi.convert}}, for converting between different two-locus models
\item \code{\link{epi.data}}, for generating data under each model
\item \code{\link{epi.design}}, for creating the design matrix to fit each model
\item \code{\link{epi.power}}, for power analysis of testing the parameters in each model
}

See the documentation for the details of using each function. Additional information on this package can be found in Walters et al. (in preparation).
}

\section{Models}{
The models that are currently implemented in this package are listed below, including a description of their respective parameters. The saturated forms of all models have equivalent fit, but can be expected to differ when parameters are constrained. Additional information on these models can be found in the listed references.

\subsection{Penetrances (\code{penetrance})}{
A model for dichotomous phenotypes that specifies the conditional probability of the phenotype (i.e. \eqn{P[y=1|x,z]}) for each set of genotypes for the two loci (\eqn{x,z}). In practice, subsets of these probabilities may be constrained to be equal to create a desired pattern of effects (e.g. see Neuman and Rice, 1992). In \pkg{epipen}, these models are specified by providing the 3 x 3 table of conditional probabilities, with columns conditional on the first locus and rows conditional on the second locus. See \code{\link{epi.convert}} for an example that utilizes this kind of model.
}

\subsection{Log Odds (\code{logodds})}{
Similar to the penetrance model, this model for dichotomous phenotypes specifies the log odds of the phenotype conditional on the genotypes (i.e. \eqn{ln(P[y=1|x,z]/(1-P[y=1|x,z]))}). Like penetrance models, these log odds models are specified in \pkg{epipen} by providing the 3 x 3 table of conditional log odds, with columns conditional on the first locus and rows conditional on the second locus.
}

\subsection{Variance Components (\code{varex})}{
Variance component models focus on decomposing the total genetic effect of the two loci into additive, dominant, and interaction components (e.g. see Falconer and Mackay, 1996). For two biallelic loci this decomposition yields eight distinct variance components. In \pkg{epipen}, decomposition of the genetic variance is performed based on the NOIA model (Alvarez-Castro and Carlborg, 2007) in order to insure that the variance components are orthogonal even when Hardy-Weinberg equilibrium does not hold. As such, the variance component model is specified using the proportion of phenotypic variance explained by each of the following genetic effects: 
\enumerate{
\item additive effect of locus A, 
\item additive effect of locus B, 
\item dominance effect of A, 
\item dominance effect B, 
\item additive by additive interaction,
\item A dominance by B additive interaction, 
\item A additive by B dominance interaction, and 
\item dominance by dominance interaction
}
}

\subsection{GLM Models}{
The remaining models all follow a generalized linear model (GLM) framework, specifying a set of variables based on the genotypes at each loci that can be used for regression. Broadly, the GLM framework defines some inverse link function \eqn{g()} between the phenotype and a linear model of the genotypes based on the distribution of the phenotypes. For the designed set of variables \eqn{(x_1, x_2, \ldots, x_J)} and phenotype \eqn{(y)},

\eqn{g(y|x_1, \ldots, x_J) = \Sigma \beta_j x_j}{g(y|x_1, \ldots, x_J) = sum(beta_j * x_j)}

The appropriate distribution is specified in \pkg{epipen} using the \code{distn} argument. For continuous phenotypes, \code{distn="normal"} sets \eqn{g(y)=E(y)} resulting in standard linear regression. For dichotomous phenotypes, \code{distn="logit"} defines \eqn{g(y)=ln[P(y)/(1-P(y))]} leading to logistic regression, and \code{distn="probit"} defines \eqn{g(y)=\Phi^{-1}(P(y))}{g(y)=qnorm(P(y))} for probit regression. See \code{\link{epi.convert}} for info on how variance components are computed for these dichotomous models by defining a latent variable \eqn{y^*=g(y)+\epsilon}{y*=g(y)+error} with some known error distribution.

Within this framework, the variables \eqn{x_j} based the genotypes for the two loci can be specified in a number of ways. A saturated model for two biallelic loci will require regression coefficients for nine such variables. A number of different parametrizations of these variables have been suggested based on various theoretical frameworks. The variants implemented in the \pkg{epipen} package are described below.

\subsection{Natural and Orthogonal Interactions - Statistical (\code{NOIA_S})}{
The NOIA Statistical model is one of two formulations of the NOIA model proposed by Alvarez-Castro and Carlborg (2007). The Statistical formulation focuses on ensuring the defined variables for the model are orthogonal even when Hardy-Weinberg equilibrium does not hold for either locus. This is done by adjusting the variables for the genotype frequencies of each locus. Full formulas are presented in Alvarez-Castro and Carlborg (2007). The coefficients for the resulting variables correspond to the additive, dominance, and interaction effects. Specifically, the parameters for \pkg{epipen} are: 
\enumerate{
\item "mu", the mean
\item "alpha_a", the average additive effect of locus A
\item "delta_a", the dominance deviation of locus A
\item "alpha_b", the average additive effect of locus B
\item "delta_b", the dominance deviation of locus B
\item "alphaalpha_ab", the additive by additive interaction
\item "deltaalpha_ab", the interaction of the dominance effect of locus A and the additive effect of locus B
\item "alphadelta_ab", the interaction of the additive effect of locus A and the dominance of locus B
\item "deltadelta_ab", the dominance by dominance interaction
}
}

\subsection{Natural and Orthogonal Interactions - Functional (\code{NOIA_F})}{
The NOIA Functional formulation follows a similar development as the NOIA Statistical, but with a focus on ensuring the parameter estimates correspond to the conventional genotypic values \eqn{a} and \eqn{d} (e.g. Falconer and Mackay, 1996). The full development is presented by Alvarez-Castro and Carlborg (2007), and only maintains orthogonal effects under certain conditions. The coefficients for the variables defined by the NOIA Functional model and implemented in \pkg{epipen} correspond to these genotype values as follows: 
\enumerate{
\item "R", the reference level
\item "a_a", the additive value of locus A
\item "d_a", the dominance deviation of locus A
\item "a_b", the additive value of locus B
\item "d_b", the dominance deviation of locus B
\item "aa_ab", the additive by additive interaction
\item "da_ab", the dominance (locus A) by additive (locus B) interaction
\item "ad_ab", the additive (locus A) by dominance (locus B) interaction
\item "dd_ab",  the dominance by dominance interaction
}
}

\subsection{General Two-Allele (\code{G2A})}{
The General Two-Allele model, proposed by Zeng et al. (2005), is designed to be orthogonal when Hardy-Weinberg equilibrium holds for all allele frequencies. This can be viewed as a special case of the NOIA Statistical model, and allows for a simpler correction for allele frequencies. Full details can be found in Zeng et al. (2005), with discussion of its relation to the NOIA model in Alvarez-Castro and Carlborg (2007). The resulting variables for the General Two-Allele model as implemented in \pkg{epipen} are as follows: 
\enumerate{
\item "mu", the mean
\item "a_a", the average additive effect of locus A
\item "d_a", the dominance deviation of locus A
\item "a_b", the average additive effect of locus B
\item "d_b", the dominance deviation of locus B
\item "aa_ab", the additive by additive interaction
\item "da_ab", the interaction of the dominance effect of locus A and the additive effect of locus B
\item "ad_ab", the interaction of the additive effect of locus A and the dominance of locus B
\item "dd_ab", the dominance by dominance interaction
}
Note that although these values are average effects as in the NOIA Functional model, the parameters are labelled based on "a" and "d" instead of "alpha" and "delta" in order to maintain consistency with Zeng et al. (2005).
}

\subsection{\eqn{F_2} (\code{F2})}{
The \eqn{F_2} model (Anderson and Kempthorne, 1954) is a more conventional model that does not adjust for allele frequencies, but still yields orthogonal effects in the special case of biallelic loci with allele frequencies of 0.5 (i.e. an \eqn{F_2} population). For a single locus, the effects defined by coding the (0,1,2) genotype as (-1, 0, 1) for the additive effect and as (-0.5, 0.5, -0.5) for the dominance effect. Interaction effects are then defined by the products of the corresponding effects for the two loci. Following the pattern of the previous models, the \pkg{epipen} package implements the \eqn{F_2} model with the following variables: 
\enumerate{
\item "mu", the mean of the F_2 population
\item "a_a", the additive value of locus A
\item "d_a", the dominance deviation of locus A
\item "a_b", the additive value of locus B
\item "d_b",  the dominance deviation of locus B
\item "aa_ab", the additive by additive interaction
\item "da_ab", the dominance (locus A) by additive (locus B) interaction
\item "ad_ab", the additive (locus A) by dominance (locus B) interaction
\item "dd_ab",  the dominance by dominance interaction
}
}

\subsection{\eqn{F_\infty}{F_infinity} (\code{Finf})}{
The \eqn{F_\infty}{F_infinity} model (Hayman and Mather, 1955) is the traditional genetic model with additive effects coded (-1, 0, 1) for the genotypes at a single locus and dominance effects coded (0,1,0). As in the \eqn{F_2} model, interaction effects are then defined by cross-products of these variables for the two loci. The resulting effects are only orthogonal under unusual circumstances. In \pkg{epipen}, the variables for the \eqn{F_\infty}{F_infinity} model follow the same naming scheme as the \eqn{F_2} model above.
}

\subsection{Unweighted (\code{unweight})}{
The Unweighted model, proposed by Cheverud and Routman (1995), uses the same (-1,0,1) coding for the additive effects at a single locus but adjusts the variable for dominance effects to (-1/3, 2/3, -1/3) coding. This change is designed to make the mean variable an unweighted average of the expected phenotype value for the nine possible genotype combinations of the two loci. Users considering using this model may wish to consider the discussion in Zeng et al. (2005) regarding other controversial claims related to this model. In the \pkg{epipen} package, the variables are again labelled following the names used for the \eqn{F_2} and \eqn{F_\infty}{F_infinity} models.
}

\subsection{Genotype (\code{genotype})}{
Cordell (2009; Supplementary Box S1) proposes a somewhat different model to allow simpler parametrization of the expected value of the phenotype for each combination of genotypes. The resulting parameters are more focused on differentiating cells of the 3x3 matrix of possible genotypes rather than indicating specific genetic effects, but the parameters still distinguish between univariate effects of locus A (\eqn{\gamma}{gamma}) and locus B (\eqn{\beta}{beta}) and interaction effect (\eqn{i}). Cordell (2009) refers to this model as a general "genotype" model. Following the naming scheme in Cordell (2009), the \pkg{epipen} package specifies the variables for the Genotype model as follows: 
\enumerate{
\item "mu", the mean of the homozygote reference cell
\item "gam_1", the effect of first minor allele at locus A
\item "gam_2", the effect of two minor alleles at locus A
\item "beta_1", the effect of first minor allele at locus B
\item "beta_2", the effect of two minor alleles at locus B
\item "i_11", the interaction deviation for one minor allele at each locus
\item "i_12", the interaction deviation for two minor alleles at locus A and one minor allele at locus B
\item "i_21", the interaction deviation for one minor allele at locus A and two minor allele at locus B
\item "i_22", the interaction deviation for two minor alleles at both loci
}
}
}
}
\author{
Raymond Walters, Charles Laurin, Gitta Lubke
}
\references{

Alvarez-Castro, J.M., and Carlborg, O. (2007). A unified model for functional and statistical epistasis and its application in quantitative trait loci analysis. \emph{Genetics}, 176: 1151-1167.

Anderson, V.L., and Kempthorne, O. (1954). A model for the study of quantitative inheritance. \emph{Genetics}, 39: 883-898.

Cheverud, J.M. and Routman, E.J. (1995). Epistasis and its contribution to genetic variance components. \emph{Genetics}, 139: 1455-1461.

Cordell, H.J. (2009). Detecting gene-gene interactions that underlie human diseases. \emph{Nature Reviews Genetics}, 10: 392-404.

Falconer, D.S. and Mackay, T.F.C. (1996). \emph{Introduction to Quantitative Genetics}. 4th edition. Pearson: Harlow, UK.

Hayman, B.I. and Mather, K. (1955). The description of genetic interactions in continuous variation. \emph{Biometrics}, 11: 69-82.

McKelvey, R., and Zavoina, W. (1975). A Statistical model for the analysis of ordinal level dependent variables. \emph{Journal of Mathematical Sociology}, 4: 103-120.

Neuman, R.J., and Rice, J.P. (1992). Two-locus models of disease. \emph{Genetic Epidemiology}, 9(5): 347-365.

Walters, R.K., Laurin, C., and Lubke, G.H. (in preparation). EpiPen: An R package to investigate two-locus epistatic models.

Zeng, Z.-B., Wang, T., and Zou, W. (2005) Modeling quantitative trait loci and interpretation of models. \emph{Genetics}, 169: 1711-1725.
}
\keyword{ package }
