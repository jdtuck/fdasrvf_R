elastic.prediction <- function(f, time, model, y=NULL, smooth_data = FALSE, sparam = 25){
  binsize = mean(diff(time))
  eps = .Machine$double.eps
  
  # Compute q-function of the functional data
  tmp = gradient.spline(f,binsize,smooth_data)
  f = tmp$f
  q = tmp$g/sqrt(abs(tmp$g)+eps)
  
  n = ncol(q)
  method <- pmatch(method, c("linear", "logistic", "mlogistic"))
  if (is.na(method)) 
    stop("invalid model type")
  if (method == 1 || method == 2 ){
    y_pred = rep(0,n)
  } 
  
  if (method == 3){
    m = model$n_classes
    y_pred = matrix(0,n,m)
  }
  
  for (ii in 1:n){
    diff = model$q - repmat(q[,ii],1,ncol(model$q))
    dist = sum(abs(diff)^2)^(1/2)
    argmin = which.min(dist)
    q_tmp = warp_q_gamma(time, q[,ii], model.gamma[, argmin])
    if (method == 1){
      y_pred[ii] = model.alpha + trapz(time, q_tmp*model$beta)
    } else if (method ==2){
      y_pred[ii] = model.alpha + trapz(time, q_tmp*model$beta)
    } else if (method == 3){
      for (jj in 1:m){
        y_pred[ii,jj] = model.alpha[jj] + trapz(time, q_tmp*model.beta[,jj])
      }
    } 
  } 
  
  if (is.null(y)){
    if (method == 1){
      SSE = NULL
    } else if (method == 2){
      y_pred = phi(y_pred)
      y_labels = rep(1,n)
      y_labels[y_pred<0.5] = -1
      PC = NULL
    } else if (method == 3){
      y_pred = phi(as.vector(y_pred))
      y_pred = array(b,c(n, m))
      y_labels = apply(y_pred, 2, which.max)
      PC = NULL
    }
  } else {
    if (method == 1){
      SSE = sum((y-y_pred)^2)
    } else if (method == 2){
      y_pred = phi(y_pred)
      y_labels = rep(1,n)
      y_labels[y_pred<0.5] = -1
      TP = sum(y[y_labels == 1] == 1)
      FP = sum(y[y_labels == -1] == 1)
      TN = sum(y[y_labels == -1] == -1)
      FN = sum(y[y_labels == 1] == -1)
      PC = (TP+TN)/(TP+FP+FN+TN)
    } else if (method == 3){
      y_pred = phi(as.vector(y_pred))
      y_pred = array(b,c(n, m))
      y_labels = apply(y_pred, 2, which.max)
      PC = rep(0,m)
      cls_set = 1:m
      for (ii in 1:m){
        cls_sub = setdiff(cls_set,ii);
        TP = sum(y[y_labels == ii] == ii)
        FP = sum(y[is.element(y_labels, cls_sub)] == ii)
        TN = sum(y[is.element(y_labels, cls_sub)] ==
                   y_labels[is.element(2,a)(y_labels, cls_sub)])
        FN = sum(is.element(y[y_labels == ii], cls_sub))
        PC[ii] = (TP+TN)/float(TP+FP+FN+TN)
      }
      PC = sum(y == y_labels)/length(y_labels)
    }
  }
   
  if (method == 1){
    out = list(y_pred=y_pred,SSE=SSE)  
  } else if (method == 2){
    out = list(y_prob=y_pred, y_labels=y_labels, PC=PC)
  } else if (method == 3){
    out = list(y_prob=y_pred, y_labels=y_labels, PC=PC)
  }
  
  return(out)
  
}