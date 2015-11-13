\name{auc}
\alias{auc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate AUC
}
\description{
Calculate AUC with or without missing
}
\usage{
auc(pred, lable, draw = F, option = "all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pred}{prediction results (real value)
%%     ~~Describe \code{pred} here~~
}
  \item{lable}{target value (binary)
%%     ~~Describe \code{lable} here~~
}
  \item{draw}{If TRUE, then draw a ROC curve
%%     ~~Describe \code{draw} here~~
}
  \item{option}{Missing Value handling
%%     ~~Describe \code{option} here~~
}
}
\details{This function is a wrapper of ROCR
%%  ~~ If necessary, more details than the description above ~~
}
\value{
return value is area under curve.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
ROCR
%% ~put references to the literature/web site here ~
}
\author{
Yifan Yang
%%  ~~who you are~~
}
\note{

%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (pred, lable, draw = F, option = "all") 
{
    if (option == "all") {
        pred.roc <- prediction(pred, lable)
        perf <- performance(pred.roc, "tpr", "fpr")
        if (draw) 
            plot(perf, col = 2, main = "GLM")
        perf.auc <- performance(pred.roc, "auc")
        perf.auc.areas <- slot(perf.auc, "y.values")
        curve.area <- mean(unlist(perf.auc.areas))
    }
    else {
        nomissing <- !is.na(lable)
        pred <- pred[nomissing]
        lable <- lable[nomissing]
        pred.roc <- prediction(pred, lable)
        perf <- performance(pred.roc, "tpr", "fpr")
        if (draw) 
            plot(perf, col = 2, main = "GLM")
        perf.auc <- performance(pred.roc, "auc")
        perf.auc.areas <- slot(perf.auc, "y.values")
        curve.area <- mean(unlist(perf.auc.areas))
    }
    return(curve.area)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ auc }

