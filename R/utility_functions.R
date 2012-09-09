cumtrapz <- function(x,y){
	m = length(y)
	
	dt = diff(x)/2
	z = c(0, cumsum(dt*(y[1:(m-1)] + y[2:m])))
	
	dim(z) = c(m,1)
	
	return(z)
}

trapzs <- function(x,y){
	M = length(y)
	out = sum(diff(x)*(y[-M]+y[-1])/2)
	return(out)
}

trapz <- function(x,y){
	M = nrow(y)
	N = ncol(y)
	out = rep(0,N)
	for (i in 1:N){
		out[i] = sum(diff(x)*(y[-M,i]+y[-1,i])/2)
	}
	return(out)
}

simpson <- function(x,y){
	M = nrow(y)
	N = ncol(y)
	
	# use  trapz if M < 3
	if (M < 3){
		out = trapz(x,y)
	}else{
		out = rep(0,N)
		dx = diff(x)
		dx1 = dx[1:(length(dx)-1)]
		dx2 = dx[2:length(dx)]
		alpha = (dx1+dx2)/dx1/6
		a0 = alpha*(2*dx1-dx2)
		a1 = alpha*(dx1+dx2)^2/dx2
		a2 = alpha*dx1/dx2*(2*dx2-dx1)
		for (i in 1:N){
			require("matrixcalc")
			out[i] = sum(a0[seq(1,length(a0),2)]*y[seq(1,M-2,2),i] + a1[seq(1,length(a1),2)]*y[seq(2,M-1,2),i]+a2[seq(1,length(a2),2)]*y[seq(3,M,2),i])
			if (M %% 2 == 0){
				A = vandermonde.matrix(x[(length(x)-2):length(x)],3)
				C = solve(A[,3:1],y[(length(y)-2):length(y)])
				out[i] = out[i] + C[1]*(x[length(x)]^3-x[(length(x)-1)]^3)/3 + C[2]*(x[length(x)]^3-x[(length(x)-1)]^2)/2 + C[3]*dx[length(dx)]
			}
		}
	}
	return(out)
}

simpson_s <- function(x,y){
	M = length(y)
	if (M < 3){
		out = trapz(x,y)
	}else{
		dx = diff(x)
		dx1 = dx[1:(length(dx)-1)]
		dx2 = dx[2:length(dx)]
		alpha = (dx1+dx2)/dx1/6
		a0 = alpha*(2*dx1-dx2)
		a1 = alpha*(dx1+dx2)^2/dx2
		a2 = alpha*dx1/dx2*(2*dx2-dx1)
		
		out = sum(a0[seq(1,length(a0),2)]*y[seq(1,M-2,2)] + a1[seq(1,length(a1),2)]*y[seq(2,M-1,2)]+a2[seq(1,length(a2),2)]*y[seq(3,M,2)])
		if (M %% 2 == 0){
			require("matrixcalc")
			A = vandermonde.matrix(x[(length(x)-2):length(x)],3)
			C = solve(A[,3:1],y[(length(y)-2):length(y)])
			out = out + C[1]*(x[length(x)]^3-x[(length(x)-1)]^3)/3 + C[2]*(x[length(x)]^2-x[(length(x)-1)]^2)/2 + C[3]*dx[length(dx)]
		}
	}
	return(out)
}

cumtraps <- function(x,y){
	M = length(y)
	dx = diff(x)
	dx1 = dx[1:(length(dx)-1)]
	dx2 = dx[2:length(dx)]
	alpha = (dx1+dx2)/dx1/6
	a0 = alpha*(2*dx1-dx2)
	a1 = alpha*(dx1+dx2)^2/dx2
	a2 = alpha*dx1/dx2*(2*dx2-dx1)
	
	require("matrixcalc")
	A = vandermonde.matrix(x[1:3],3)
	C = solve(A[,3:1],y[1:3])
	z = rep(0,M)
	z[2] = C[1]*(x[2]^3-x[1]^3)/3 + C[2]*(x[2]^2-x[1]^2)/2 + C[3]*dx[1]
	z[seq(3,length(z),2)] = cumsum(a0[seq(1,length(a0),2)]*y[seq(1,M-2,2)] + a1[seq(1,length(a1),2)]*y[seq(2,M-1,2)] + a2[seq(1,length(a1),2)]*y[seq(3,M,2)])
	z[seq(4,length(z),2)] = cumsum(a0[seq(2,length(a0),2)]*y[seq(2,M-2,2)] + a1[seq(2,length(a1),2)]*y[seq(3,M-1,2)] + a2[seq(2,length(a1),2)]*y[seq(4,M,2)])+z[2]
	
	
	return(z)
}

pvecnorm <-function(v,p){
	sum(abs(v)^p)^(1/p)
}

pvecnorm2 <-function(dt,x){
	sqrt(sum(abs(x)*abs(x))*dt)
}

gradient <- function(f,binsize){
	n = nrow(f)
	if (is.null(n)){
		f = as.vector(f)
		n = length(f)
		g = rep(0,n)
		h = binsize*(1:n)
		# Take forward differences on left and right edges
		g[1] = (f[2] - f[1])/(h[2]-h[1])
		g[n] = (f[n] - f[(n-1)])/(h[length(h)]-h[(length(h)-1)])
		
		# Take centered differences on interior points
		h = h[3:n]-h[1:(n-2)]
		g[2:(n-1)] = (f[3:n]-f[1:(n-2)])/h[1]
	}else {
		f = as.matrix(f)
		p = ncol(f)
		g = matrix(0,n,p)
		h = binsize*(1:n)
		# Take forward differences on left and right edges
		g[1,] = (f[2,] - f[1,])/(h[2]-h[1])
		g[n,] = (f[n,] - f[(n-1),])/(h[length(h)]-h[(length(h)-1)])
		
		# Take centered differences on interior points
		h = h[3:n]-h[1:(n-2)]
		g[2:(n-1),] = (f[3:n,]-f[1:(n-2),])/h[rep(1,p)]
	}
	
	return(g)
}

SqrtMeanInverse <- function(gam){
	n = nrow(gam)
	T1 = ncol(gam)
	dt = 1/(T1-1)
	psi = matrix(0,n,T1-1)
	for (i in 1:n){
		psi[i,] = sqrt(diff(gam[i,])/dt)
	}
	
	# Find direction
	mu = psi[1,]
	tt = 1
	maxiter = 5
	lvm = rep(0,1,maxiter)
	vec = matrix(0,n,T1-1)
	for (iter in 1:maxiter){
		for (i in 1:n){
			v = psi[i,] - mu
			len = acos(sum(mu*psi[i,])*dt)
			if (len > 0.0001){
				vec[i,] = (len/sin(len))*(psi[i,] - cos(len)*mu)
			}else{
				vec[i,] = rep(0,1,T1-1)
			}	
		}
		vm = colMeans(vec)
		lvm[iter] = sqrt(sum(vm*vm)*dt)
		if (lvm[iter] > 0.0001){
			mu = cos(tt*lvm[iter])*mu + (sin(tt*lvm[iter])/lvm[iter])*vm
		}
		
	}
	
	gam_mu = c(0,cumsum(mu*mu))/T1
	gam_mu = (gam_mu - min(gam_mu))/(max(gam_mu)-min(gam_mu))
	gamI = invertGamma(gam_mu)
	return(gamI)
}

invertGamma <- function(gam){
	N = length(gam)
	x = (0:(N-1))/(N-1)
	tmp = approx(gam,x,xout=x)
	gamI = tmp$y
	gamI[N] = 1
	gamI = gamI/gamI[N]
	return(gamI)
}

randomGamma <- function(gam,num){
	out = SqrtMean(gam)
	mu = out$mu
	psi = out$psi
	vec = out$vec
	
	K = cov(t(vec))
	out = svd(K)
	s = out$d
	U = out$u
	n = 5
	TT = nrow(vec) + 1
	vm = rowMeans(vec)
	
	rgam = matrix(0,num,TT)
	for (k in 1:num){
		a = rnorm(n)
		v = rep(0,length(vm))
		for (i in 1:n){
			v = v + a[i]*sqrt(s[i])*U[,i]
		}
		vn = pvecnorm(v,2)/sqrt(TT)
		psi = cos(vn)*mu + sin(vn)*v/vn
		tmp = rep(0,TT)
		tmp[2:TT] = cumsum(psi*psi)
		rgam[k,] = tmp/TT
	}
	return(rgam)
}

f_K_fold <- function(Nobs,K=5){
	rs <- runif(Nobs)
	id <- seq(Nobs)[order(rs)]
	k <- as.integer(Nobs*seq(1,K-1)/K)
	k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
	k[,1] <- k[,1]+1
	l <- lapply(seq.int(K),function(x,k,d) 
		list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
				 test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
	return(l)
}

cov_samp <- function(x,y=NULL){
	x = scale(x,scale=F)
	N = dim(x)[1]
	if (length(y) == 0){
		sigma = 1/N * t(x) %*% x
	}else{
		y = scale(y,scale=F)
		sigma = 1/N * t(x) %*% y
	}
	
	return(sigma)
}

diffop <- function(n,binsize = 1){
	m = matrix(0,nrow=n,ncol=n)
	diag(m[-1,]) <- 1
	diag(m) <- -2
	diag(m[,-1]) <- 1
	m = t(m) %*% m
	m[1,1] = 6
	m[n,n] = 6
	m = m/(binsize^4)
	return(m)
}

geigen <- function (Amat, Bmat, Cmat) 
{
	Bdim <- dim(Bmat)
	Cdim <- dim(Cmat)
	if (Bdim[1] != Bdim[2]) 
		stop("BMAT is not square")
	if (Cdim[1] != Cdim[2]) 
		stop("CMAT is not square")
	p <- Bdim[1]
	q <- Cdim[1]
	s <- min(c(p, q))
	if (max(abs(Bmat - t(Bmat)))/max(abs(Bmat)) > 1e-10) 
		stop("BMAT not symmetric.")
	if (max(abs(Cmat - t(Cmat)))/max(abs(Cmat)) > 1e-10) 
		stop("CMAT not symmetric.")
	Bmat <- (Bmat + t(Bmat))/2
	Cmat <- (Cmat + t(Cmat))/2
	Bfac <- chol(Bmat)
	Cfac <- chol(Cmat)
	Bfacinv <- solve(Bfac)
	Cfacinv <- solve(Cfac)
	Dmat <- t(Bfacinv) %*% Amat %*% Cfacinv
	if (p >= q) {
		result <- svd2(Dmat)
		values <- result$d
		Lmat <- Bfacinv %*% result$u
		Mmat <- Cfacinv %*% result$v
	}
	else {
		result <- svd2(t(Dmat))
		values <- result$d
		Lmat <- Bfacinv %*% result$v
		Mmat <- Cfacinv %*% result$u
	}
	geigenlist <- list(values, Lmat, Mmat)
	names(geigenlist) <- c("values", "Lmat", "Mmat")
	return(geigenlist)
}

svd2 <- function (x, nu = min(n, p), nv = min(n, p), LINPACK = FALSE) 
{
	dx <- dim(x)
	n <- dx[1]
	p <- dx[2]
	svd.x <- try(svd(x, nu, nv, LINPACK))
	if (class(svd.x) == "try-error") {
		nNA <- sum(is.na(x))
		nInf <- sum(abs(x) == Inf)
		if ((nNA > 0) || (nInf > 0)) {
			msg <- paste("sum(is.na(x)) = ", nNA, "; sum(abs(x)==Inf) = ", 
									 nInf, ".  'x stored in .svd.x.NA.Inf'", sep = "")
			assign(".svd.x.NA.Inf", x, envir = .GlobalEnv)
			stop(msg)
		}
		attr(x, "n") <- n
		attr(x, "p") <- p
		attr(x, "LINPACK") <- LINPACK
		.x2 <- c(".svd.LAPACK.error.matrix", ".svd.LINPACK.error.matrix")
		.x <- .x2[1 + LINPACK]
		assign(.x, x, envir = .GlobalEnv)
		msg <- paste("svd failed using LINPACK = ", LINPACK, 
								 " with n = ", n, " and p = ", p, ";  x stored in '", 
								 .x, "'", sep = "")
		warning(msg)
		svd.x <- try(svd(x, nu, nv, !LINPACK))
		if (class(svd.x) == "try-error") {
			.xc <- .x2[1 + (!LINPACK)]
			assign(.xc, x, envir = .GlobalEnv)
			stop("svd also failed using LINPACK = ", !LINPACK, 
					 ";  x stored in '", .xc, "'")
		}
	}
	svd.x
}

cbcolorscale <- function(n){
	require(ggplot2)
	colors = c("#56B4E9","#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
	if (n <= length(colors)){
		cbgColourPalette <- scale_colour_manual(values=colors[1:n])
	}else{
		colors2 = rep(colors,round(n/length(colors)))
		cbgColourPalette <- scale_colour_manual(values=sample(colors2))
	}
	return(cbgColourPalette)
}

cumtrapzmid <- function(x,y,c){
	a = length(x)
	mid = round(a/2)
	
	# case < mid
	fn = rep(0,a)
	tmpx = x[seq(mid-1,1,-1)]
	tmpy = y[seq(mid-1,1,-1)]
	tmp = c + cumtrapz(tmpx,tmpy)
	fn[1:mid-1] = rev(tmp)
	
	# case >= mid
	fn[mid:a] = c + cumtrapz(x[mid:a],y[mid:a])
	
	return(fn)
}