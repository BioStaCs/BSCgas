\name{lift.chart}
\alias{lift.chart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lift Chart
}
\description{
Model evaluation: lift chart
}
\usage{
lift.chart(prediction, labels, dec.point = 0.5, nclass = 10, colorize = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prediction}{
  predicted result (double)
%%     ~~Describe \code{prediction} here~~
}
  \item{labels}{
  real target lable (binary)
%%     ~~Describe \code{labels} here~~
}
  \item{dec.point}{
  cut off point, default=.5
%%     ~~Describe \code{dec.point} here~~
}
  \item{nclass}{
  Classifed in to n class, default=10
%%     ~~Describe \code{nclass} here~~
}
  \item{colorize}{
  whether use color, default=F
%%     ~~Describe \code{colorize} here~~
}
}
\details{
Lift charts
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Picture
}
\references{
http://sweb.uky.edu/~yya234
%% ~put references to the literature/web site here ~
}
\author{
Yifan Yang
%%  ~~who you are~~
}
\note{
Lift Charts
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
AUC
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (prediction, labels, dec.point = 0.5, nclass = 10, colorize = T) 
{
    if (nclass < 5) {
        print("nclass should be larger than (or equals) 5")
        nclass = 5
    }
    n <- length(labels)
    if (length(prediction) != n) {
        print("lengths of prediction and labels are different.")
        return(F)
    }
    if (n < nclass) {
        print("length of prediction should be larger than N-class.")
        return(F)
    }
    if (n\%\%nclass == 0) {
        col <- rainbow(nclass)
        n.perclass <- n/nclass
        pred.rate <- rep(0, nclass)
        tmp <- as.integer(sort(prediction, decreasing = T) > 
            dec.point)
        for (i in 1:nclass) {
            pred.rate[]
        }
    }
    else {
        col <- rainbow(nclass) + 1
        print("Extra 1 new class")
        n.perclass <- as.integer(n/nclass)
        n.res <- n\%\%nclass
        pred.rate <- rep(0, nclass)
    }
    return(pred.rate)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ lift charts }
\keyword{ model evaluation }% __ONLY ONE__ keyword per line
