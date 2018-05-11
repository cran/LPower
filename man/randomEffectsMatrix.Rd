\name{randomEffectsMatrix}
\alias{randomEffectsMatrix}
\title{
Calculates the variance covariance matrix for a multivariate normal vector when there are random effects.
}
\description{
Computes the variance covariance matrix of an \eqn{m} vector which results from a random effects model.
}
\usage{
randomEffectsMatrix(zMatrix, vs, sigma2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zMatrix}{
An \eqn{m X p} design matrix which specifies how \eqn{p} random variables with zero mean are linearly related to the m-dimensional normal vector.
}
  \item{vs}{
  The \eqn{p X p} variance covariance matrix of the random effects,
}
  \item{sigma2}{
The error variance.
}
}
\details{
We assume that \eqn{y_{t}=\mu_t+\Sigma \gamma_j z_{t,j}+\sigma^2 \epsilon},
where \eqn{\gamma_j} are random variables with mean \eqn{0} and and variance covariance \code{vs}, and \eqn{z} is \code{zMatrix}, \eqn{\epsilon} is a standard normal random variable.
The \code{zMatrix} could be a list of matricies
}
\value{
Either a single variance covariance matrix or a list of them if zMatrix is a list.
}
\author{
David A. Schoenfeld
}
\seealso{\code{\link{LPower}},\code{\link{randomSlopesMatrix}}
}
\examples{
 #Creates random variance covariance matrix for random follow up model
 #where baseline is random among patients and all follow up have a compound symetry structure
 #from a common random effect
vars=randomEffectsMatrix(cbind(rep(1,5),matrix(c(0,rep(1,4)),5,1)),
            matrix(c(31.8,.8527,.8527,.6687),2,2),2.7085)
LPower(sample_size=40,power=.8,
          xMatrix=list(cbind(1,c(0,rep(1,4)),0),cbind(1,c(0,rep(1,4)),c(0,rep(1,4)))),vMatrix=vars)
#Creates random variance covariance matrix for random slopes model
vars=randomEffectsMatrix(cbind(rep(1,5),0:4),
            matrix(c(31.8,.8527,.8527,.6687),2,2),2.7085)
LPower(sample_size=40,power=.8,
       xMatrix=list(cbind(1,0:4,0),cbind(1,0:4,0:4)),vMatrix=vars)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }% use one of  RShowDoc("KEYWORDS")

