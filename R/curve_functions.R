calculatecentroid <- function(beta){
    n = nrow(beta)
    T1 = ncol(beta)

    betadot = matrix(0, n, T1)
    for (i in 1:n){
        betadot[i,] = gradient(beta[i,], 1.0/(T1-1))
    }

    normbetadot = rep(0, T1)
    integrand = matrix(0, n, T1)
    for (i in 1:T1){
        normbetadot[i] = norm(betadot[,i], 2)
        integrand[,i] = beta[,i] * normbetadot[i]
    }
    scale = trapz(seq(0,1,length.out=T1), normbetadot)
    centroid = trapz(seq(0,1,length.out=T1), integrand, 2)/scale

    return(centroid)
}


innerprod_q2 <- function(q1, q2){
    T1 = ncol(q1)
    val = sum(sum(q1*q2))/T1

    return(val)
}


find_best_rotation <- function(q1, q2){
    eps = .Machine$double.eps
    n = nrow(q1)
    T1 = ncol(q1)
    A = q1%*%t(q2)
    out = svd(A)
    s = out$d
    U = out$u
    V = out$v
    if (abs(det(U)*det(V)-1) < (10*eps)){
        S = diag(1,n)
    } else {
        S = diag(1,n)
        S[,n] = -S[,n]
    }
    R = U*S*t(V)
    q2new = R%*%q2

    return(list(q2new=q2new, R=R))
}


calculate_variance <- function(beta){
    n = nrow(beta)
    T1 = ncol(beta)
    betadot = matrix(0, n, T1)
    for (i in 1:n){
        betadot[i,] = gradient(beta[i,], 1.0/(T1-1))
    }
    normbetadot = rep(0,T)
    centroid = calculatecentroid(beta)
    integrand = array(0, c(n,n,T1))
    time = seq(0,1,length.out=T1)
    for (i in 1:T1){
        normbetadot[i] = norm(betadot[,i],2)
        a1 = beta[,i] - centroid
        integrand[,,i] = a1 %*% t(a1) * normbetadot[i]
    }
    l = trapz(t, normbetadot)
    variance = trapz(t, integrand, 3)
    varaince = variance / l

    return(variance)
}
