\name{svm.prehandle}
\alias{svm.prehandle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepair for SVM(e1071/libsvm) 
}
\description{
SVM matrix scale
}
\usage{
svm.prehandle(x,att.biostacs=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{input matrix/data.frame.
  }
  \item{att.biostacs}{Trace Information. Attention.
  }
}
\details{For e1071/libsvm
}
\value{
\item{mat}{return matrix.}
\item{scale}{scale for each feature.}
\item{factor.nu}{Index of Factors. Need attention.}
%% ...
}
\references{
e1071, libsvm
%% ~put references to the literature/web site here ~
}
\author{
Yifan Yang
%%  ~~who you are~~
}
\note{
PREHANDLE
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
mat<-data.frame(X=1:10,Y=factor(sample(1:5,10,replace=T)),Z=sample(letters[1:3],10,replace=T))
example<-svm.prehandle(mat,att.biostacs=T)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ auc }

