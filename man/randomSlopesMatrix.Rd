\name{randomSlopesMatrix}
\alias{randomSlopesMatrix}

\title{
Creates the xMatrix and zMatrix, and attrition rates for a two treatment clinical trial analyzed using the random slopes model.
}
\description{
In the random slopes model each patient has a linear trajectory over time with a random intercept and slope. The intercepts are assumed to be the same for each of two treatment groups and the treatment effect is measured by the difference in average slopes.
}
\usage{
randomSlopesMatrix(visit, vs, sigma2, dropPerMonth,baselineTreatment=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{visit}{
A vector of visit times or a list of two visit time vectors if the treatments have different visit times.
}
  \item{vs}{
The variance covariance matrix of the intercept and slope random effects.
}
  \item{sigma2}{
The error variance.
}
  \item{dropPerMonth}{
Either a single number which is the attrition rate or a vector of attrition rates for each visit. Note this would have length one less than the number of visits since the attrition after the last visit would not be used.
}
  \item{baselineTreatment}{A logical indicating whether their treatment is in the model as a main effect. In a random slopes model the effect or treatment is measured by the treatment-time interaction.

  }
}
\details{
This calculates the matrices for the random slopes model
\eqn{y_{t}=\mu+\beta_1 t+ \beta_2 t*I(\rm{rx}=1)+u+b t+\sigma^2 \epsilon},
where \eqn{u,b,\epsilon} are random variables. Note that a treatment main effect is not included in the model by default, because in a randomized study the treatments should be the same at the baseline visit. This practice may vary.
}
\value{
A list of \code{xMatrix,vMatrix,attritionRates} for input into \code{\link{LPower}}
}

\references{
Q Yi and T. Panzarella. Estimating sample size for tests on
trends across repeated measurements with missing data based on
the interaction term in a mixed model. Control Clin Trials,
23(5):481–96, 2002.

}
\author{
David A. Schoenfeld
}

\seealso{\code{\link{LPower}}
}
\examples{
vars=randomSlopesMatrix(list(c(0,1,2,4,6),c(0,1,2,4,6)),
                matrix(c(31.8,.8527,.8527,.6687),2,2),2.7085,.02)
LPower(sample_size=40,power=.8,xMatrix=vars$xMatrix,
vMatrix=vars$vMatrix,attritionRates=vars$attritionRates)
}
\keyword{ ~design }% use one of  RShowDoc("KEYWORDS")

