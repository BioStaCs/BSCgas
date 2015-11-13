\name{glmauc.greed}
\alias{glmauc.greed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GAS
%%  ~~function to do ... ~~
}
\description{
A novel Greedy AUC-maximize Strategy
}
\usage{
glmauc.greed(data.train, data.test = NA, varlist, terms, MaxVar = 8, auc.inf0 = rep(0, MaxVar), varlist.in = NULL, trace.info = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.train}{Train Data which is needed to be a data-frame.
%%     ~~Describe \code{data.train} here~~
}
  \item{data.test}{If NA, the output is trained data, otherwise, one should offer test data, which is of the same form of train data. 
%%     ~~Describe \code{data.test} here~~
}
  \item{varlist}{Features that are used in this model.
%%     ~~Describe \code{varlist} here~~
}
  \item{terms}{Features that are used in this model.
%%     ~~Describe \code{terms} here~~
}
  \item{MaxVar}{Max number of features that are allowed.
%%     ~~Describe \code{MaxVar} here~~
}
  \item{auc.inf0}{trace info.
%%     ~~Describe \code{auc.inf0} here~~
}
  \item{varlist.in}{trace info.
%%     ~~Describe \code{varlist.in} here~~
}
  \item{trace.info}{trace info.
%%     ~~Describe \code{trace.info} here~~
}
}
\details{GAS
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{term } Terms predicted.
\item{model} Final predicted model.
\item{auc.train } AUC of train data. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
http://sweb.uky.edu/~yya234/wordpress
%% ~put references to the literature/web site here ~
}
\author{
Yifan Yang
%%  ~~who you are~~
}
\note{
The license is BSD
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
GAS
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data.train, data.test = NA, varlist, terms, MaxVar = 8, 
    auc.inf0 = rep(0, MaxVar), varlist.in = NULL, trace.info = T) 
{
    tic <- proc.time()
    L <- length(varlist)
    trace.matrix <- matrix(0, 8, L)
    ii = 1
    while (ii <= MaxVar) {
        auc.list <- rep(0, length(varlist))
        for (i in varlist) {
            model.term <- paste.my(terms[c(i, varlist.in)])
            model.parse <- paste("glm(lable~", model.term, ",data=data.train,family=binomial('probit'))", 
                sep = " ")
            model.tmp <- eval(parse(text = model.parse))
            pred <- logit(predict.glm(model.tmp, data.train[, 
                -16], type = "link"))
            auc.list[i] <- auc(pred, data.train$lable)
            if (trace.info) 
                cat(ii, "-th Turn, testing var+=", terms[i], 
                  "    AUC=", auc.list[i], "\n")
        }
        auc.inf0[ii] <- max(auc.list)
        if (ii > 1) {
            if (auc.inf0[ii] > auc.inf0[ii - 1]) {
                varlist.in <- c(varlist.in, varlist[which.max(auc.list)])
                trace.matrix[ii, ] <- sort(auc.list, index.return = T)$ix
                varlist <- varlist[!(varlist == which.max(auc.list))]
                ii = ii + 1
            }
            else {
                if (trace.info) 
                  cat(ii, "-th Turn, testing var --\n")
                ii = ii - 1
                pop.ind <- trace.matrix[ii, 1]
                trace.matrix[ii, ] = c(trace.matrix[ii, 2:L], 
                  0)
                if (trace.matrix[ii, 1] == 0) 
                  break
                varlist.in <- c(varlist.in[1:ii], trace.matrix[ii, 
                  1])
                varlist <- c(varlist, pop.ind)
                ii = ii + 1
            }
        }
        else {
            varlist.in <- c(varlist.in, varlist[which.max(auc.list)])
            trace.matrix[ii, 1:(L + 1 - ii)] <- sort(auc.list, 
                index.return = T)$ix
            varlist <- varlist[!(varlist == which.max(auc.list))]
            ii = ii + 1
        }
    }
    model.term <- paste.my(terms[varlist.in])
    model.parse <- paste("glm(lable~", model.term, ",data=data.train,family=binomial('probit'))", 
        sep = " ")
    model.final <- eval(parse(text = model.parse))
    if (!is.na(data.test)) {
        pred <- logit(predict.glm(model.final, data.test[, -16], 
            type = "link"))
        auc.test <- auc(pred, data.test$lable, T)
        cat("AUC on test set is", auc.inf0[MaxVar], "\nAUC on test set is:", 
            auc.test, "\n")
        toc <- proc.time()
        toc - tic
        return(list(term = model.term, model = model.final, prediction = pred, 
            auc.train = auc.inf0, auc.test = auc.test))
    }
    else {
        return(list(term = model.term, model = model.final, prediction = pred, 
            auc.train = auc.inf0))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GAS }
\keyword{ AUC }% __ONLY ONE__ keyword per line
