#' Elastic Prediction for functional multinomial logistic PCR Model
#'
#' This function performs prediction from an elastic multinomial logistic fPCR regression model
#'  with phase-variability
#'
#' @param object Object of class inheriting from "elastic.pcr.regression"
#' @param newdata An optional matrix in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param y An optional vector of labels to calculate PC. If omitted, PC is NULL
#' @param ... additional arguments affecting the predictions produced
#' @return Returns a list containing
#' \item{y_pred}{predicted probabilities of the class of newdata}
#' \item{y_labels}{class labels of newdata}
#' \item{PC}{probability of classification per class}
#' \item{PC.comb}{total probability of classification}
#' @keywords srvf alignment regression
#' @references J. D. Tucker, J. R. Lewis, and A. Srivastava, “Elastic
#'  Functional Principal Component Regression,” Statistical Analysis and Data
#'  Mining, 10.1002/sam.11399, 2018.
#' @export
predict.mlpcr <- function(object, newdata=NULL, y=NULL, ...){
    m <- ncol(object$Y)
    if (is.null(newdata)){
        n <- nrow(object$pca$coef)
        y_pred <- matrix(0,n,m)
        for (ii in 1:n){
            for (jj in 1:m){
                y_pred[ii,jj] <- object$alpha[jj]+sum(object$pca$coef[ii,] * object$b[,jj])
            }
        }

        y_pred <- y_pred/apply(apply(y_pred,2,abs),1,sum)
        y_pred <- phi(c(y_pred))
        y_pred <- matrix(y_pred,ncol=m)
        y_labels <- apply(y_pred, 1, which.min)
    } else {
        n <- ncol(newdata)
        M <- nrow(newdata)
        y_pred <- matrix(0,n,m)

        # Project newdata onto basis
        time <- object$warp_data$time
        lambda <- object$warp_data$lambda
        omethod <- object$warp_data$omethod
        q <- f_to_srvf(newdata, time)
        mq <- rowMeans(object$warp_data$qn)
        fn <- matrix(0,M,n)
        qn <- matrix(0,M,n)
        gam <- matrix(0,M,n)
        for (ii in 1:n){
            gam[, ii] <- optimum.reparam(mq,time,q[,ii],time,lambda,omethod)
            fn[, ii] <- stats::approx(time,newdata[,ii],xout=(time[length(time)]-time[1])*gam[, ii] +
                               time[1])$y
            qn[, ii] <- f_to_srvf(fn[, ii], time)
        }

        m_new <- sign(fn[object$pca$id,])*sqrt(abs(fn[object$pca$id,]))  # scaled version
        qn1 <- rbind(qn,m_new)

        pca.method1 <- object$pca.method
        pca.method <- pmatch(pca.method1, c("combined", "vert", "horiz")) # 1 - combined, 2 - vert, 3 - horiz
        if (is.na(pca.method))
            stop("invalid method selection")

        U <- object$pca$U
        no <- ncol(object$pca$U)
        if (pca.method == 1){
            C <- object$pca$C
            TT <- length(time)
            mu_g <- object$pca$mu_g
            mu_psi <- object$pca$mu_psi
            vec <- matrix(0,M,n)
            psi <- matrix(0,TT,n)
            binsize <- mean(diff(time))
            for (i in 1:n){
                psi[,i] <- sqrt(gradient(gam[,i],binsize))
            }

            for (i in 1:n){
                vec[,i] <- inv_exp_map(mu_psi, psi[,i])
            }

            g <- rbind(qn1,C*vec)
            a <- matrix(0,n,no)
            for (i in 1:n){
                for (j in 1:no){
                    a[i,j] <- (g[,i]-mu_g)%*%U[,j]
                }
            }
        }

        if (pca.method == 2){
            a <- matrix(0,n,no)
            for (k in 1:no){
                for (i in 1:n){
                    a[i,k] <- (qn1[,i]-object$pca$mqn)%*%U[,k]
                }
            }
        }

        if (pca.method == 3){
            a <- matrix(0,n,no)
            mu_psi <- object$pca$mu
            vec <- matrix(0,M,n)
            TT <- length(time)
            psi <- matrix(0,TT,n)
            binsize <- mean(diff(time))
            for (i in 1:n){
                psi[,i] <- sqrt(gradient(gam[,i],binsize))
            }

            for (i in 1:n){
                vec[,i] <- inv_exp_map(mu_psi, psi[,i])
            }

            vm <- rowMeans(object$pca$vec)

            for (k in 1:no){
                for (i in 1:n){
                    a[i,k] <- sum((vec[,i]-vm)*U[,k])
                }
            }
        }

        for (ii in 1:n){
            for (jj in 1:m){
                y_pred[ii,jj] <- object$alpha[jj]+sum(a[ii,] * object$b[,jj])
            }
        }

        y_pred <- y_pred/apply(apply(y_pred,2,abs),1,sum)
        y_pred <- phi(c(y_pred))
        y_pred <- matrix(y_pred,ncol=m)
        y_labels <- apply(y_pred, 1, which.min)
    }

    if (missing(newdata)){
        PC <- rep(0,m)
        cls_set <- 1:m
        for (ii in 1:m){
            cls_sub <- setdiff(cls_set,ii);
            TP <- sum(object$y[y_labels == ii] == ii)
            FP <- sum(object$y[is.element(y_labels, cls_sub)] == ii)
            TN <- sum(object$y[is.element(y_labels, cls_sub)] ==
                          y_labels[is.element(y_labels, cls_sub)])
            FN <- sum(is.element(object$y[y_labels == ii], cls_sub))
            PC[ii] <- (TP+TN)/(TP+FP+FN+TN)
        }
        PC.comb <- sum(object$y == y_labels)/length(y_labels)
    } else if (missing(y)){
        PC <- NULL
        PC.comb <- NULL
    } else {
        PC <- rep(0,m)
        cls_set <- 1:m
        for (ii in 1:m){
            cls_sub <- setdiff(cls_set,ii);
            TP <- sum(y[y_labels == ii] == ii)
            FP <- sum(y[is.element(y_labels, cls_sub)] == ii)
            TN <- sum(y[is.element(y_labels, cls_sub)] ==
                          y_labels[is.element(y_labels, cls_sub)])
            FN <- sum(is.element(y[y_labels == ii], cls_sub))
            PC[ii] <- (TP+TN)/(TP+FP+FN+TN)
        }
        PC.comb <- sum(y == y_labels)/length(y_labels)
    }
    out <- list(y_pred=y_pred, y_labels=y_labels, PC=PC, PC.comb=PC.comb)

    return(out)

}
