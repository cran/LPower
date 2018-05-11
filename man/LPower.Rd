\name{LPower}
\alias{LPower}

\title{
Calculates either the power, sample size, or detectable effect for a longitudinal study with a repeated measures design.
}
\description{
Provide two of three parameters (power, sample size, detectable effect) and it supplies the third in a design with repeated measures.  It requires the design matrix, and the variance covariance matrix of the repeated measures. It can also take into account of attrition at each of the time points in the model.
}
\usage{
LPower(rx_effect = NULL, sample_size = NULL, allocationRatio = c(1, 1), power = NULL,
          contrast=c(rep(0,dim(xMatrix[[1]])[2]-1),1), xMatrix, vMatrix,
          attritionRates = 0, alpha = 0.025,
          simulate=FALSE,nsims=1000,betas=c(rep(0,dim(xMatrix[[1]])[2]-1),1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rx_effect}{
The size of the effect to be detected. Set to Null if this is the parameter to be calculated.
}
  \item{sample_size}{
The total sample size of the study.
}
  \item{allocationRatio}{
The allocation ratio, the allocation to each arm in the study or to each group of patients that have a distinct design matrix.
}
  \item{power}{
The desired power. Null if the power is to be computed.
}
  \item{contrast}{
The contrast to be estimated, the default value, which is valid if the xMatrix paramter is a list is c(0,...1). That is the last parameter is the effect of interest.
}
  \item{xMatrix}{
A list of matricies giving the regression coefficients for each patient group, note that the all must have the same dimensions.
}
  \item{vMatrix}{
A list of variance covariance matricies for each patient group. A single matrix also will work if all patients groups have the same variance covariance matrix. All must have the same dimensions.
}
  \item{attritionRates}{
A vector which is the rate of attrition between each visit. Attrition is considered to be exponential between visits.
}
  \item{alpha}{
The significance level
}
  \item{simulate}{
  Logical, indicating that you also want to run a simulation to calculate the power given the calculated sample size or detectable difference.
  }
  \item{nsims}{Number of simulations to use}
  \item{betas}{Coefficient value for simulations. Note that betas is rescaled so that the value of the contrast is \code{rx_effect}. The code is
  \code{betas=rep(rx_effect/(matrix(betas,1,m[2])\%*\%mcontrast),m[2])*betas}
  }
}

\value{
A vector giving the detectable difference, sample size and power.
}

\note{
The code to analyse the model used in the simulation is something like
\code{
    mod2 = nlme::gls(y~X1+X2+X3, correlation = corSymm(form = ~visits|subject),
        weights = varIdent(form =~1|visits), na.action = na.omit, data = df)}
with test statistic.
    \code{sum(contrast * mod2$coefficients)/sqrt(t(mcontrast) %*%
        mod2$varBeta \%*\% mcontrast)}
}

\references{
Diggle, P., Liang, K.Y. and Zeger, S.L., 1994. Longitudinal data analysis. New York: Oxford University Press, 5, pp 59, ISBN:9780198522843
}
\author{
David A. Schoenfeld
}


\seealso{
 \code{\link{randomEffectsMatrix}}, \code{\link{randomSlopesMatrix}}
}
\examples{
#This would be what would be used for an analysis of covariance assuming
#a correlation of 0.3 and a standard deviation of 5.46.
LPower(sample_size=60,power=.8,
     xMatrix=list(matrix(c(1,1,0,1,0,0),2,3),matrix(c(1,1,0,1,0,1),2,3)),
     vMatrix=5.46^2*matrix(c(1,0.3,0.3,1),2,2),attritionRates=0.1)
}
\keyword{~design }% use one of  RShowDoc("KEYWORDS")

