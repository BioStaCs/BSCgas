###############
#Copyright: BioStaCs Group;
#Lisence:   GPL v2
###############


###performance evaluation
#1. lift.chart
#2. ROC curve
#3. confusion Matrix

#fusion matrix

require(e1071)
require(compiler)
library(c060)


Identify.biostacs<-function(X){ ## check vegnette for detail!!!
  X.real<-model.matrix(~.,data=X)[,-1]
  return (X.real);
}

##simple newton raphson method
nf.biostacs<-function(func,dfunc,x){
  if (abs(dfunc(x))<10*.Machine$double.eps){
    return (x);
  }else{
    return(x-func(x)/dfunc(x));
  }
}

#mean, check vignette for detail!!!
mode.biostacs<-function(vec){
  re<-table(vec)
  re<-(names(re)[re==max(re)])
  l<-length(re);
  if(l>1) {
    re<-(re[sample(l,1)]);
  }
  return(re);
}

confusion.matrix<-function (
  prediction, labels, test = T, draw = F, family = "binomial", 
  sep.point = 0.5) {
  pre <- as.integer(prediction >= sep.point)
  TP <- sum((pre == 1) * (labels == 1))
  FP <- sum((pre == 1) * (labels == 0))
  TN <- sum((pre == 0) * (labels == 0))
  FN <- sum((pre == 0) * (labels == 1))
  mat <- matrix(c(TP, FP, FN, TN), byrow = T, 2, 2)
  colnames(mat)=c("Real +","Real -")
  rownames(mat)=c("Pred +","Pred -")
  if (draw) {
    par(mfrow = c(1, 2))
    plot(1, type = "n", xaxt = "n", yaxt = "n", xlab = " ", 
         ylab = " ", main = "Fusion Table")
    heatmap(mat)
  }
  if (test) {
    
    cat('\n--------------------------------------------------\n')
    cat('The Seprate point is:',sep.point,'\n');
    cat('--------------------------------------------------\n')
    acc = (TP + TN)/sum(sum(mat))
    
    re <- chisq.test(mat, simulate.p.value = TRUE, B = 1000)
    print(re)
    if (re$p.value>0.05){cat('We don\' have enough evidence to reject the classifier is random. We need a better model.\n\n')}else{
    cat('We have strong evidence to say the classifier is not random. It is a potential choice.\n\n');  
    }
    MCC = (TP * TN - FN * FP)/sqrt(prod(c((TP + FP), (TP + 
      FN), (TN + FP), (TN + FN))))
    cat("\nMatthews correlation coefficient:\t", MCC, "\n");
    cat('Sn=',TP/(TP+FN),'\n');
    cat('Sp=',TN/(TN+FP),'\n');
    cat("ACC=\n", acc, "\n")
  }
  return(mat)
}

confusion.matrixOptSep<-function (
  prediction, labels
  ) {
  ff=Vectorize(function(sep.point,...){
    pre <- as.integer(prediction >= sep.point)
    TP <- sum((pre == 1) * (labels == 1))
    FP <- sum((pre == 1) * (labels == 0))
    TN <- sum((pre == 0) * (labels == 0))
    FN <- sum((pre == 0) * (labels == 1))
    MCC = (TP * TN - FN * FP)/sqrt(prod(c((TP + FP), (TP + 
                                                        FN), (TN + FP), (TN + FN))))
    return(MCC)
  })
  sep.points=(1:99)/100;
  sep.points[which.max(ff(sep.points))]
}


lift.chart<-function(prediction,labels,dec.point=.5,nclass=10,colorize=T){
  if (nclass<5) {
    print('nclass should be larger than (or equals) 5');
    nclass=5;
  }
  
  n<-length(labels);
  if (length(prediction)!=n) {print('lengths of prediction and labels are different.');return (F);}
  if (n<nclass)  {print('length of prediction should be larger than N-class.');return (F);}
  
  if (n%%nclass==0) {
    col<-rainbow(nclass);
    n.perclass<-n/nclass;
    pred.rate <- rep(0,nclass);
    tmp<-as.integer(sort(prediction,decreasing = T)>dec.point);
    for(i in 1:nclass){pred.rate[i]};
  }else{
    col<-rainbow(nclass+1);
    print('Extra 1 new class');
    n.perclass<-as.integer(n/nclass);
    n.res<-n%%nclass;
    pred.rate <- rep(0,nclass);
  }
  return(pred.rate);   
}

###############
# lars: Efron
###############

lars.biostacs <-function(x, y, type = c("lasso", "lar", "forward.stagewise","stepwise"), trace = FALSE,
           normalize=TRUE, intercept=TRUE, Gram, 
           eps = .Machine$double.eps,  max.steps, use.Gram = TRUE)
  {
    call <- match.call()
    type <- match.arg(type)
    TYPE <- switch(type,
                   lasso = "LASSO",
                   lar = "LAR",
                   forward.stagewise = "Forward Stagewise",
                   stepwise = "Forward Stepwise")
    if(trace)
      cat(paste(TYPE, "sequence\n"))
    
    nm <- dim(x)
    n <- nm[1]
    m <- nm[2]
    im <- inactive <- seq(m)
    one <- rep(1, n)
    vn <- dimnames(x)[[2]]  
    ### Center x and y, and scale x, and save the means and sds
    if(intercept){
      meanx <- drop(one %*% x)/n
      x <- scale(x, meanx, FALSE)	# centers x
      mu <- mean(y)
      y <- drop(y - mu)
    }
    else {
      meanx <- rep(0,m)
      mu <- 0
      y <- drop(y)
    }
    if(normalize){
      normx <- sqrt(drop(one %*% (x^2)))
      nosignal<-normx/sqrt(n) < eps
      if(any(nosignal))# ignore variables with too small a variance
      {
        ignores<-im[nosignal]
        inactive<-im[-ignores]
        normx[nosignal]<-eps*sqrt(n)
        if(trace)
          cat("LARS Step 0 :\t", sum(nosignal), "Variables with Variance < eps; dropped for good\n")	#
      }
      else ignores <- NULL #singularities; augmented later as well
      names(normx) <- NULL
      x <- scale(x, FALSE, normx)	# scales x
    }
    else {
      normx <- rep(1,m)
      ignores <- NULL
    }
    if(use.Gram & missing(Gram)) {
      if(m > 500 && n < m)
        cat("There are more than 500 variables and n<m;\nYou may wish to restart and set use.Gram=FALSE\n"
        )
      if(trace)
        cat("Computing X'X .....\n")
      Gram <- t(x) %*% x	#Time saving
    }
    Cvec <- drop(t(y) %*% x)
    ssy <- sum(y^2)	### Some initializations
    residuals <- y
    if(missing(max.steps))
      max.steps <- 8*min(m, n-intercept)
    beta <- matrix(0, max.steps + 1, m)	# beta starts at 0
    lambda=double(max.steps)
    Gamrat <- NULL
    arc.length <- NULL
    R2 <- 1
    RSS <- ssy
    first.in <- integer(m)
    active <- NULL	# maintains active set
    actions <- as.list(seq(max.steps))	
    # a signed index list to show what comes in and out
    drops <- FALSE	# to do with type=="lasso" or "forward.stagewise"
    Sign <- NULL	# Keeps the sign of the terms in the model
    R <- NULL	###
    ### Now the main loop over moves
    ###
    k <- 0
    while((k < max.steps) & (length(active) < min(m - length(ignores),n-intercept)) )
    {
      action <- NULL
      C <- Cvec[inactive]	#
      ### identify the largest nonactive gradient
      Cmax <- max(abs(C))
      if(Cmax<eps*100){ # the 100 is there as a safety net
        if(trace)cat("Max |corr| = 0; exiting...\n")
        break
      }
      k <- k + 1
      lambda[k]=Cmax
      ### Check if we are in a DROP situation
      if(!any(drops)) {
        new <- abs(C) >= Cmax - eps
        C <- C[!new]	# for later
        new <- inactive[new]	# Get index numbers
        ### We keep the choleski R  of X[,active] (in the order they enter)
        for(inew in new) {
          if(use.Gram) {
            R <- updateR(Gram[inew, inew], R, drop(Gram[
              inew, active]), Gram = TRUE,eps=eps)
          }
          else {
            R <- updateR(x[, inew], R, x[, active], Gram
                         = FALSE,eps=eps)
          }
          if(attr(R, "rank") == length(active)) {
            ##singularity; back out
            nR <- seq(length(active))
            R <- R[nR, nR, drop = FALSE]
            attr(R, "rank") <- length(active)
            ignores <- c(ignores, inew)
            action <- c(action,  - inew)
            if(trace)
              cat("LARS Step", k, ":\t Variable", inew, 
                  "\tcollinear; dropped for good\n")	#
          }
          else {
            if(first.in[inew] == 0)
              first.in[inew] <- k
            active <- c(active, inew)
            Sign <- c(Sign, sign(Cvec[inew]))
            action <- c(action, inew)
            if(trace)
              cat("LARS Step", k, ":\t Variable", inew, 
                  "\tadded\n")	#
          }
        }
      }
      else action <-  - dropid
      Gi1 <- backsolve(R, backsolvet(R, Sign))	
      ### Now we have to do the forward.stagewise dance
      ### This is equivalent to NNLS
      dropouts<-NULL
      if(type == "forward.stagewise") {
        directions <- Gi1 * Sign
        if(!all(directions > 0)) {
          if(use.Gram) {
            nnls.object <- nnls.lars(active, Sign, R, 
                                     directions, Gram[active, active], trace = 
                                       trace, use.Gram = TRUE,eps=eps)
          }
          else {
            nnls.object <- nnls.lars(active, Sign, R, 
                                     directions, x[, active], trace = trace, 
                                     use.Gram = FALSE,eps=eps)
          }
          positive <- nnls.object$positive
          dropouts <-active[-positive]
          action <- c(action, -dropouts)
          active <- nnls.object$active
          Sign <- Sign[positive]
          Gi1 <- nnls.object$beta[positive] * Sign
          R <- nnls.object$R
          C <- Cvec[ - c(active, ignores)]
        }
      }
      A <- 1/sqrt(sum(Gi1 * Sign))
      w <- A * Gi1	# note that w has the right signs
      if(!use.Gram) u <- drop(x[, active, drop = FALSE] %*% w)	###
      ### Now we see how far we go along this direction before the
      ### next competitor arrives. There are several cases
      ###
      ### If the active set is all of x, go all the way
      if( (length(active) >=  min(n-intercept, m - length(ignores) ) )|type=="stepwise") {
        gamhat <- Cmax/A
      }
      else {
        if(use.Gram) {
          a <- drop(w %*% Gram[active,  - c(active,ignores), drop = FALSE])
        }
        else {
          a <- drop(u %*% x[,  - c(active, ignores), drop=FALSE])
        }
        gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))	
        ### Any dropouts will have gam=0, which are ignored here
        gamhat <- min(gam[gam > eps], Cmax/A)	
      }
      if(type == "lasso") {
        dropid <- NULL
        b1 <- beta[k, active]	# beta starts at 0
        z1 <-  - b1/w
        zmin <- min(z1[z1 > eps], gamhat)
        if(zmin < gamhat) {
          gamhat <- zmin
          drops <- z1 == zmin
        }
        else drops <- FALSE
      }
      beta[k + 1,  ] <- beta[k,  ]
      beta[k + 1, active] <- beta[k + 1, active] + gamhat * w
      if(use.Gram) {
        Cvec <- Cvec - gamhat * Gram[, active, drop = FALSE] %*% w
      }
      else {
        residuals <- residuals - gamhat * u
        Cvec <- drop(t(residuals) %*% x)
      }
      Gamrat <- c(Gamrat, gamhat/(Cmax/A))
      arc.length <- c(arc.length, gamhat)	
      ### Check if we have to drop any guys
      if(type == "lasso" && any(drops)) {
        dropid <- seq(drops)[drops]	
        #turns the TRUE, FALSE vector into numbers
        for(id in rev(dropid)) {
          if(trace)
            cat("Lasso Step", k+1, ":\t Variable", active[
              id], "\tdropped\n")
          R <- downdateR(R, id)
        }
        dropid <- active[drops]	# indices from 1:m
        beta[k+1,dropid]<-0  # added to make sure dropped coef is zero
        active <- active[!drops]
        Sign <- Sign[!drops]
      }
      if(!is.null(vn))
        names(action) <- vn[abs(action)]
      actions[[k]] <- action
      inactive <- im[ - c(active, ignores)]
      if(type=="stepwise")Sign=Sign*0
    }
    beta <- beta[seq(k + 1), ,drop=FALSE ]	#
    lambda=lambda[seq(k)]
    dimnames(beta) <- list(paste(0:k), vn)	### Now compute RSS and R2
    if(trace)
      cat("Computing residuals, RSS etc .....\n")
    residuals <- y - x %*% t(beta)
    beta <- scale(beta, FALSE, normx)
    RSS <- apply(residuals^2, 2, sum)
    R2 <- 1 - RSS/RSS[1]
    actions=actions[seq(k)]
    netdf=sapply(actions,function(x)sum(sign(x)))
    df=cumsum(netdf)### This takes into account drops
    if(intercept)df=c(Intercept=1,df+1)
    else df=c(Null=0,df)
    rss.big=rev(RSS)[1]
    df.big=n-rev(df)[1]
    if(rss.big<eps|df.big<eps)sigma2=NaN
    else
      sigma2=rss.big/df.big
    Cp <- RSS/sigma2 - n + 2 * df
    attr(Cp,"sigma2")=sigma2
    attr(Cp,"n")=n
    object <- list(call = call, type = TYPE, df=df, lambda=lambda,R2 = R2, RSS = RSS, Cp = Cp, 
                   actions = actions[seq(k)], entry = first.in, Gamrat = Gamrat, 
                   arc.length = arc.length, Gram = if(use.Gram) Gram else NULL, 
                   beta = beta, mu = mu, normx = normx, meanx = meanx)
    class(object) <- "lars"
    object
  }



###############
#GAS algorithom
###############

#Author: Yifan Yang@uky

mode.avg<-function(vec){
  re<-table(vec)
  re<-(names(re)[re==max(re)])
  l<-length(re);
  if(l>1) {
    re<-(re[sample(l,1)]);
  }
  return(re);
}

logit<-function(x)(1/(1+exp(-x)));


liftcharts.biostacs<-function (pred, lable,ADD=F,caption=' ') 
{
  pred.roc <- prediction(pred, lable)
  perf <- performance(pred.roc, "lift", "rpp")
  plot(perf,colorize=T, main = caption,add=ADD)
  
}


auc.biostacs<-function (pred, lable, draw = F, option = "all",col=2,ADD=F,caption='GLM') 
{
  if (option == "all") {
    pred.roc <- prediction(pred, lable)
    perf <- performance(pred.roc, "tpr", "fpr")
    if (draw) 
      plot(perf, col = col, main = caption,add=ADD)
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
      plot(perf, col = col, main = caption)
    perf.auc <- performance(pred.roc, "auc")
    perf.auc.areas <- slot(perf.auc, "y.values")
    curve.area <- mean(unlist(perf.auc.areas))
  }
  return(curve.area);
}

outer.ch1<-function(x,y){
  crossing<-function(x,y,ch=':'){
    paste(x,y,sep=ch)
  }
  I=length(x);
  J=length(y);
  re<-matrix(' ',I,J)
  for (i in 1:I){
    for(j in 1:J){
      if (x[i]==y[j]){
        tmp=x[i];
      }else{
        tmp=crossing(x[i],y[j]);
      }
      re[i,j]=tmp;
    }
  }
  return (re);
}

paste.my1<-function(x,sep='+'){
  L=length(x)
  if (L==1) {
    re=x;
  }else{
    re=x[1];
    for (i in 2:L){
      re=paste(re,x[i],sep=sep);
    }
  }
  return (re);
}
outer.ch<-cmpfun(outer.ch1);
paste.my<-cmpfun(paste.my1);

glmauc.greed<-function(data.train,data.test=NA,varlist,terms,MaxVar=8,auc.inf0=rep(0,MaxVar),varlist.in=NULL,trace.info=T){
  
  lable.nu<-names(data.train);lable.nu<-(1:length(lable.nu))[lable.nu=='lable']
  tic<-proc.time();
  L<-length(varlist);
  trace.matrix<-matrix(0,MaxVar,L);
  ii=1;
  while (ii <=MaxVar){
    auc.list<-rep(0,length(varlist));
    for (i in varlist){
      model.term<-paste.my(terms[c(i,varlist.in)]);
      model.parse<-paste('glm(lable~',model.term,",data=data.train,family=binomial('probit'))",sep=' ');
      model.tmp<-eval(parse(text=model.parse));
      pred<-logit(predict.glm(model.tmp,data.train[,-lable.nu],type='link'));
      auc.list[i]<-auc.biostacs(pred,data.train$lable);
      if (trace.info) cat(ii,'-th Turn, testing var+=',terms[i],'    AUC=',auc.list[i],'\n');
    }
    auc.list[is.na(auc.list)]=0;
    auc.inf0[ii]<-max(auc.list);
    if (ii>1){
      if (auc.inf0[ii]>auc.inf0[ii-1]){
        varlist.in <- c(varlist.in,which.max(auc.list));
        sort(auc.list,index.return = T)$ix->tmp.ix
        trace.matrix[ii,1:length(tmp.ix)] <- tmp.ix
        varlist <- varlist[!(varlist==which.max(auc.list))];
        ii=ii+1;
      }else{
        if (trace.info) cat(ii,'-th Turn, testing var --\n');
        ii=ii-1;
        pop.ind<-trace.matrix[ii,1];
        trace.matrix[ii,]=c(trace.matrix[ii,2:L],0);
        if (trace.matrix[ii,1]==0) break;
        varlist.in <- c(varlist.in[1:ii],trace.matrix[ii,1]);
        varlist <- c(varlist,pop.ind);
        ii=ii+1;
      }
    }else{
      varlist.in <- c(varlist.in,which.max(auc.list));
      trace.matrix[ii,1:(L+1-ii)] <- sort(auc.list,index.return = T)$ix 
      varlist <- varlist[!(varlist==which.max(auc.list))];
      ii=ii+1;
    }
  }
  model.term<-paste.my(terms[varlist.in]);
  model.parse<-paste('glm(lable~',model.term,",data=data.train,family=binomial('probit'))",sep=' ');
  model.final<-eval(parse(text=model.parse));
  
  if (!is.na(data.test)){
    pred<-logit(predict.glm(model.final,data.test[,-lable.nu],type='link'));  
    #auc.test<-auc.biostacs(pred,data.test$lable,T);
    #cat('AUC on test set is',auc.inf0[MaxVar],'\nAUC on test set is:', auc.test,'\n');
    toc<-proc.time();
    cat('\nRunning time:',toc-tic,'\n');
    return (list(term=model.term,
                 model=model.final,
                 prediction=pred,
                 auc.train=auc.inf0
                 #,auc.test=auc.test
    )); 
  }else{
    return(list(term=model.term,
                model=model.final,
                auc.train=auc.inf0));
  }
}

normalize.matrix.biostacs<-function(x){##nomalized to range near [-1,1]
  n<-dim(x)[1];
  f.L2<-function(x){x/sqrt(sum(x^2)/n)}
  tmp<-apply(x,2,f.L2);
  return(list(mat=tmp,scale=apply(x,2,function(x){sqrt(sum(x^2)/n)})));
}

svm.prehandle<-function(x,att.biostacs=F,keep=T){
  if(!keep){
    ISMATRIX<-is.matrix(x);
    if (ISMATRIX){
      re<-normalize.matrix.biostacs(x);
      FACTOR.biostacs<-NA;
    } else{
      type.x<-do.call(cbind,lapply(x[1,],class));
      tmp.x<-x;
      FACTOR.biostacs<-(1:length(x[1,]))[ type.x=='factor'];
      for (i in FACTOR.biostacs){
        cat(i,'th feature may cause error!!!Transform Factor into numeric!!! Using expand method here.\n');
        tmp.x[,i]<-as.numeric(tmp.x[,i]);  
      }   
      re<-normalize.matrix.biostacs(tmp.x);
    }
    if (att.biostacs){
      cat('Usage: return value include the scale parameters and train matrix!');
    }
    push.list.biostacs<-function(re,addname,add){
      re[[length(re)+1]]=add;
      names(re)[length(re)]<-addname;
      return(re);
    }
    rerere=push.list.biostacs(re,'factor.nu',FACTOR.biostacs)
  }else{
    rerere=list(mat=x,scale=rep(1,dim(x)[1]));
  }

  return(rerere);
}

svm.tune.biostacs<-function(X,lable,control=list(gamma = 3^(-10:1), 
                                       cost = seq(0.5,4,0.5),
                                       tunecontrol = tune.control(sampling = "fix")
                                       )
                  ){
  dimdim<-dim(X)[1]
  # svm requires tuning
  re<-svm.prehandle(X);
  data.train=as.data.frame(re$mat);
  scale<-re$scale;
  x.svm.tune <- tune(svm, lable~., data = cbind(data.train,lable),
                     ranges = control)
  ##default cv fold =10 in tunecontrols
  gc();
  weights=table(as.factor(lable));
  weights[2] -> wtp;
  weights[2] <- weights[1];
  weights[1] <- wtp;  
  
  x.svm <- svm(lable~.,weights=weights, data =cbind(data.train,lable), cost= x.svm.tune$best.parameters$cost, gamma= x.svm.tune$best.parameters$gamma, probability = TRUE);
  return (list(model=x.svm,tune=x.svm.tune,scale=scale));  
}

svm.prediction<-function(X,model){
  x.svm=model$model;
  x.svm.prob <- predict(x.svm,type="prob",newdata=X,probability = TRUE);
return (x.svm.prob);
}

glm.tuning.biostacs<-function(){

#Define Model Controls

MultiControl <- trainControl(workers = 2, #2 cores

  method = 'repeatedcv',

 number = 10, #10 Folds

 repeats = 25, #25 Repeats

 classProbs = TRUE,

 returnResamp = "all",

 summaryFunction = twoClassSummary,#Use 2 class-summary function to get AUC

 computeFunction = mclapply) #Use the parallel apply function

}
