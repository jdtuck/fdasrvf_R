#' Elastic Prediction for functional logisitc PCR Model
#'
#' This function performs prediction from an elastic logistic fPCR regression model
#'  with phase-variability
#'
#' @param object Object of class inheriting from "elastic.pcr.regression"
#' @param newdata An optional matrix in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param y An optional vector of labels to calculate PC. If omitted, PC is NULL
#' @param ... additional arguments affecting the predictions produced
#' @return Returns a list containing
#' \item{y_pred}{predicted probabilites of the class of newdata}
#' \item{y_labels}{class labels of newdata}
#' \item{PC}{probability of classification}
#' @keywords srvf alignment regression
#' @references J. D. Tucker, J. R. Lewis, and A. Srivastava, “Elastic
#'  Functional Principal Component Regression,” Statistical Analysis and Data
#'  Mining, 10.1002/sam.11399, 2018.
#' @export
predict.lpcr <- function(object, newdata=NULL, y=NULL, ...){
    if (is.null(newdata)){
        n <- nrow(object$pca$coef)
        y_pred <- rep(0,n)
        for (ii in 1:n){
            y_pred[ii] <- object$alpha + sum(object$pca$coef[ii,] * object$b)
        }
        y_pred <- phi(y_pred)
        y_labels <- rep(1,n)
        y_labels[y_pred<0.5] <- -1
    } else {
        n <- ncol(newdata)
        M <- nrow(newdata)
        y_pred <- rep(0,n)

        # Project newdata onto basis
        time <- object$warp_data$time
        lambda <- object$warp_data$lambda
        omethod <- object$warp_data$omethod
        q <- f_to_srvf(newdata, time)
        mq <- object$warp_data$mqn
        fn <- matrix(0,M,n)
        qn <- matrix(0,M,n)
        gam <- matrix(0,M,n)
        for (ii in 1:n){
            gam[, ii] <- optimum.reparam(mq,time,q[,ii],time,lambda,omethod)
            fn[, ii] <- approx(time,newdata[,ii],xout=(time[length(time)]-time[1])*gam[, ii] +
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
            y_pred[ii] <- object$alpha + sum(a[ii,] * object$b)
        }

        y_pred <- phi(y_pred)
        y_labels <- rep(1,n)
        y_labels[y_pred<0.5] <- -1
    }

    if (missing(newdata)){
        TP <- sum(object$y[y_labels == 1] == 1)
        FP <- sum(object$y[y_labels == -1] == 1)
        TN <- sum(object$y[y_labels == -1] == -1)
        FN <- sum(object$y[y_labels == 1] == -1)
        PC <- (TP+TN)/(TP+FP+FN+TN)
    } else if (missing(y)){
        PC <- NULL
    } else {
        TP <- sum(y[y_labels == 1] == 1)
        FP <- sum(y[y_labels == -1] == 1)
        TN <- sum(y[y_labels == -1] == -1)
        FN <- sum(y[y_labels == 1] == -1)
        PC <- (TP+TN)/(TP+FP+FN+TN)
    }
    out <- list(y_pred=y_pred, y_labels=y_labels, PC=PC)

    return(out)

}
