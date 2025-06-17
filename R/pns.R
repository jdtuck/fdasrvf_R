# PNS  The Principal Nested Spheres code (PNS) for spheres and shapes has
# been written by Kwang-Rae Kim, and builds closely on the original matlab
# code for PNS by Sungkyu Jung, code from I. Dryden (shapes) with numerical changes
# for warping functions
pcscore2sphere3 <- function (n.pc, X.hat, Xs, Tan, V) {
  d = nrow(Tan)
  n = ncol(Tan)
  W = matrix(NA, d, n)
  for (i in 1:n) {
    W[, i] =  acos(sum(Xs[i, ] * X.hat)) * Tan[, i] / sqrt(sum(Tan[, i]^2))
  }
  lambda = matrix(NA, n, d)
  for (i in 1:n) {
    for (j in 1:n.pc) {
      lambda[i, j] = sum(W[, i] * V[, j])
    }
  }
  U = matrix(0, n, d)
  for (i in 1:n) {
    for (j in 1:n.pc) {
      U[i, ] = U[i, ] + lambda[i, j] * V[, j]
    }
  }
  S.star = matrix(NA, n, n.pc + 1)
  for (i in 1:n) {
    U.norm = sqrt(sum(U[i, ]^2))
    S.star[i, ] = c(cos(U.norm), sin(U.norm) / U.norm * lambda[i, 1:n.pc])
  }
  return(S.star)
}

Enorm <- function(a)
{
  return(sqrt(sum(diag(t(
    a
  ) %*% a))))
}


fastpns <- function (x,
                     n.pc = "Full",
                     sphere.type = "seq.test",
                     mean.type = "Frechet",
                     alpha = 0.1,
                     R = 100,
                     nlast.small.sphere = 1,
                     output = TRUE) {
  n <- dim(x)[2]
  pdim <- dim(x)[1]
  if (n.pc == "Full") {
    n.pc = min(c(pdim - 1 , n - 1))
  }
  Xs <- t(x)
  for (i in 1:n) {
    Xs[i, ] <- Xs[i, ] / Enorm(Xs[i, ])
  }
  muhat <- apply(Xs, 2, mean)
  muhat <- muhat / Enorm(muhat)
  TT <- Xs
  for (i in 1:n) {
    TT[i, ] <- Xs[i, ] - sum(Xs[i, ] * muhat) * muhat
  }
  pca <- stats::prcomp(TT)
  pcapercent <- sum(pca$sdev[1:n.pc]^2 / sum(pca$sdev^2))
  cat(c("Initial PNS subsphere dimension", n.pc + 1, "\n"))
  cat(c(
    "Percentage of variability in PNS sequence",
    round(pcapercent *         100, 2),
    "\n"
  ))
  TT <- t(TT)
  ans <- pcscore2sphere3(n.pc, muhat, Xs, TT, pca$rotation)
  Xssubsphere <- t(ans)
  out <- pns(
    (Xssubsphere),
    sphere.type = sphere.type,
    mean.type = mean.type,
    alpha = alpha,
    R = R,
    nlast.small.sphere = nlast.small.sphere,
    output = output
  )
  out$percent <- out$percent * pcapercent
  cat(
    c(
      "Percent explained by 1st three PNS scores out of total variability:",
      "\n",
      round(out$percent[1:3], 2),
      "\n"
    )
  )
  out$spheredata <- (Xssubsphere)
  out$pca <- pca
  out$muhat <- muhat
  out$n.pc <- n.pc
  out
}

fastPNSe2s <- function(res , pns) {
  out <- pns
  GG <- PNSe2s(res , out$PNS)
  n <- dim(GG)[2]
  muhat <- pns$muhat
  n.pc <- pns$n.pc### now work out the PC scores for the original high-dimensional coordinates
  s <- acos(GG[1, ])
  HH <- diag(s / sin(s)) %*% t(GG[2:(n.pc + 1), ])
  ones <- rep(1, times = n)#Preferred approx back on sphere (it is unit size) #This is exact if n.pc = "Full"
  approx1 <- t(GG[2:(n.pc + 1), ]) %*% t(out$pca$rotation[, 1:(n.pc)]) +
    diag(cos(s)) %*% ones %*% t(muhat) / Enorm(muhat)
  approx1
}

fastpns_biplot <- function(pns, varnames) {
  pns1 <- pns
  nd <- dim(pns$resmat)[1] + 1
  ndhigh <- length(pns$muhat)
  grDevices::palette(grDevices::rainbow(min(ndhigh, 1024)))
  res1 <- cbind(c((20:(-20)) / 10 * stats::sd(pns1$resmat[1, ])) , matrix(0, 41, nd -
                                                                            2))
  if (nd > 3) {
    res2 <- cbind(cbind(matrix(0, 41, 1) , c((20:(
      -20
    )) / 10 * stats::sd(pns1$resmat[2, ]))) , matrix(0, 41, nd - 3))
  }
  if (nd <= 3) {
    res2 <- cbind(cbind(matrix(0, 41, 1) , c((20:(
      -20
    )) / 10 * stats::sd(pns1$resmat[2, ]))))
  }
  mshape <- fastPNSe2s(t(res1) * 0 , pns1)
  aa1 <- fastPNSe2s(t(res1) , pns1) - mshape
  aa2 <- fastPNSe2s(t(res2) , pns1) - mshape
  nl <- dim(aa1)[1]
  aa1 <- t(aa1)
  aa2 <- t(aa2)
  plot(
    aa1[1, ],
    aa2[1, ],
    xlim = c(min(aa1), max(aa1)) ,
    type = "n",
    col = 2,
    ylim = c(min(aa2), max(aa2)) ,
    xlab = "PNS1",
    ylab = "PNS2"
  )
  for (i in 1:(ndhigh)) {
    graphics::lines(aa1[i, ], aa2[i, ], col = i)
    graphics::arrows(aa1[i, 2], aa2[i, 2], aa1[i, 1], aa2[i, 1], col = i)
    graphics::text(aa1[i, 1], aa2[i, 1], varnames[i], col = i, cex = 1)
  }
  graphics::title("fast PNS biplot")
  grDevices::palette("default")
}

pns = function(x,
               sphere.type = "seq.test",
               mean.type = "Frechet",
               alpha = 0.1,
               R = 100,
               nlast.small.sphere = 1,
               output = TRUE)
{
  n = ncol(x)
  k = nrow(x)
  PNS = list()

  if (abs(sum(apply(x^2, 2, sum)) - n) > 1e-8)
  {
    stop("Error: Each column of x should be a unit vector, ||x[ , i]|| = 1.")
  }
  svd.x = svd(x, nu = nrow(x))
  uu = svd.x$u
  maxd = which(svd.x$d < 1e-15)[1]
  if (is.na(maxd) | k > n)
  {
    maxd = min(k, n) + 1
  }
  nullspdim = k - maxd + 1
  d = k - 1
  if (output) {
    cat("Message from pns() : dataset is on ", d, "-sphere. \n", sep = "")
  }
  if (nullspdim > 0)
  {
    if (output) {
      cat(" .. found null space of dimension ",
          nullspdim,
          ", to be trivially reduced. \n",
          sep = "")
    }
  }

  if (d == 2) {
    PNS$spherePNS <- t(x)
  }

  resmat = matrix(NA, d, n)
  orthaxis = list()
  orthaxis[[d - 1]] = NA
  dist = rep(NA, d - 1)
  pvalues = matrix(NA, d - 1, 2)
  ratio = rep(NA, d - 1)
  currentSphere = x
  if (nullspdim > 0)
  {
    for (i in 1:nullspdim)
    {
      oaxis = uu[, ncol(uu) - i + 1]
      r = pi / 2
      pvalues[i, ] = c(NaN, NaN)
      res = acos(t(oaxis) %*% currentSphere) - r
      orthaxis[[i]] = oaxis
      dist[i] = r
      resmat[i, ] = res
      NestedSphere = rotMat(oaxis) %*% currentSphere
      currentSphere = NestedSphere[1:(k - i), ] /
        repmat.pns(matrix(sqrt(1 - NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), k - i, 1)
      uu = rotMat(oaxis) %*% uu
      tmp.uu = 1 - uu[nrow(uu), ]^2
      tmp.uu[tmp.uu<0] = 0
      sqrt.uu = sqrt(tmp.uu)
      uu = uu[1:(k - i), ] / repmat.pns(matrix(sqrt.uu, nrow = 1), k - i, 1)
      if (output) {
        cat(d - i + 1,
            "-sphere to ",
            d - i,
            "-sphere, by ",
            "NULL space \n",
            sep = "")
      }
    }
  }
  if (sphere.type == "seq.test")
  {
    if (output) {
      cat(" .. sequential tests with significance level ",
          alpha,
          "\n",
          sep = "")
    }
    isIsotropic = FALSE
    for (i in (nullspdim + 1):(d - 1))
    {
      if (!isIsotropic)
      {
        sp = getSubSphere(x = currentSphere, geodesic = "small")
        center.s = sp$center
        r.s = sp$r
        resSMALL = acos(t(center.s) %*% currentSphere) - r.s
        sp = getSubSphere(x = currentSphere, geodesic = "great")
        center.g = sp$center
        r.g = sp$r
        resGREAT = acos(t(center.g) %*% currentSphere) - r.g
        pval1 = LRTpval(resGREAT, resSMALL, n)
        pvalues[i, 1] = pval1
        if (pval1 > alpha)
        {
          center = center.g
          r = r.g
          pvalues[i, 2] = NA
          if (output) {
            cat(
              d - i + 1,
              "-sphere to ",
              d - i,
              "-sphere, by GREAT sphere, p(LRT) = ",
              pval1,
              "\n",
              sep = ""
            )
          }
        } else {
          pval2 = vMFtest(currentSphere, R)
          pvalues[i, 2] = pval2
          if (pval2 > alpha)
          {
            center = center.g
            r = r.g
            if (output) {
              cat(
                d - i + 1,
                "-sphere to ",
                d - i,
                "-sphere, by GREAT sphere, p(LRT) = ",
                pval1,
                ", p(vMF) = ",
                pval2,
                "\n",
                sep = ""
              )
            }
            isIsotropic = TRUE
          } else {
            center = center.s
            r = r.s
            if (output) {
              cat(
                d - i + 1,
                "-sphere to ",
                d - i,
                "-sphere, by SMALL sphere, p(LRT) = ",
                pval1,
                ", p(vMF) = ",
                pval2,
                "\n",
                sep = ""
              )
            }
          }
        }
      } else if (isIsotropic) {
        sp = getSubSphere(x = currentSphere, geodesic = "great")
        center = sp$center
        r = sp$r
        if (output) {
          cat(
            d - i + 1,
            "-sphere to ",
            d - i,
            "-sphere, by GREAT sphere, restricted by testing vMF distn",
            "\n",
            sep = ""
          )
        }
        pvalues[i, 1] = NA
        pvalues[i, 2] = NA
      }
      res = acos(t(center) %*% currentSphere) - r
      orthaxis[[i]] = center
      dist[i] = r
      resmat[i, ] = res
      cur.proj = project.subsphere(x = currentSphere,
                                   center = center,
                                   r = r)
      NestedSphere = rotMat(center) %*% currentSphere
      currentSphere = NestedSphere[1:(k - i), ] /
        repmat.pns(matrix(sqrt(1 - NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), k - i, 1)

      if (nrow(currentSphere) == 3)
      {
        PNS$spherePNS = t(currentSphere)
      }

      if (nrow(currentSphere) == 2)
      {
        PNS$circlePNS = t(cur.proj)
      }

    }
  } else if (sphere.type == "BIC") {
    if (output) {
      cat(" .. with BIC \n")
    }
    for (i in (nullspdim + 1):(d - 1))
    {
      sp = getSubSphere(x = currentSphere, geodesic = "small")
      center.s = sp$center
      r.s = sp$r
      resSMALL = acos(t(center.s) %*% currentSphere) - r.s
      sp = getSubSphere(x = currentSphere, geodesic = "great")
      center.g = sp$center
      r.g = sp$r
      resGREAT = acos(t(center.g) %*% currentSphere) - r.g
      BICsmall = n * log(mean(resSMALL^2)) + (d - i + 1 + 1) * log(n)
      BICgreat = n * log(mean(resGREAT^2)) + (d - i + 1) * log(n)
      if (output) {
        cat("BICsm: ", BICsmall, ", BICgr: ", BICgreat, "\n", sep = "")
      }
      if (BICsmall > BICgreat)
      {
        center = center.g
        r = r.g
        if (output) {
          cat(d - i + 1,
              "-sphere to ",
              d - i,
              "-sphere, by ",
              "GREAT sphere, BIC \n",
              sep = "")
        }
      } else {
        center = center.s
        r = r.s
        if (output) {
          cat(d - i + 1,
              "-sphere to ",
              d - i,
              "-sphere, by ",
              "SMALL sphere, BIC \n",
              sep = "")
        }
      }
      res = acos(t(center) %*% currentSphere) - r
      orthaxis[[i]] = center
      dist[i] = r
      resmat[i, ] = res

      cur.proj = project.subsphere(x = currentSphere,
                                   center = center,
                                   r = r)

      NestedSphere = rotMat(center) %*% currentSphere
      currentSphere = NestedSphere[1:(k - i), ] /
        repmat.pns(matrix(sqrt(1 - NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), k - i, 1)

      if (nrow(currentSphere) == 3)
      {
        PNS$spherePNS = t(currentSphere)
      }

      if (nrow(currentSphere) == 2)
      {
        PNS$circlePNS = t(cur.proj)
      }


    }
  } else if (sphere.type == "small" | sphere.type == "great") {
    pvalues = NaN
    for (i in (nullspdim + 1):(d - 1))
    {
      sp = getSubSphere(x = currentSphere, geodesic = sphere.type)
      center = sp$center
      r = sp$r
      res = acos(t(center) %*% currentSphere) - r
      orthaxis[[i]] = center
      dist[i] = r
      resmat[i, ] = res

      cur.proj = project.subsphere(x = currentSphere,
                                   center = center,
                                   r = r)

      NestedSphere = rotMat(center) %*% currentSphere
      currentSphere = NestedSphere[1:(k - i), ] /
        repmat.pns(matrix(sqrt(1 - NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), k - i, 1)

      if (nrow(currentSphere) == 3)
      {
        PNS$spherePNS = t(currentSphere)
      }

      if (nrow(currentSphere) == 2)
      {
        PNS$circlePNS = t(cur.proj)
      }

    }
  } else if (sphere.type == "bi.sphere") {
    if (nlast.small.sphere < 0)
    {
      cat("!!! Error from pns(): \n")
      cat("!!! nlast.small.sphere should be >= 0. \n")
      return(NULL)
    }
    mx = (d - 1) - nullspdim
    if (nlast.small.sphere > mx)
    {
      cat("!!! Error from pns(): \n")
      cat("!!! nlast.small.sphere should be <= ",
          mx,
          " for this data. \n",
          sep = "")
      return(NULL)
    }
    pvalues = NaN
    if (nlast.small.sphere != mx)
    {
      for (i in (nullspdim + 1):(d - 1 - nlast.small.sphere))
      {
        sp = getSubSphere(x = currentSphere, geodesic = "great")
        center = sp$center
        r = sp$r
        res = acos(t(center) %*% currentSphere) - r
        orthaxis[[i]] = center
        dist[i] = r
        resmat[i, ] = res

        cur.proj = project.subsphere(x = currentSphere,
                                     center = center,
                                     r = r)

        NestedSphere = rotMat(center) %*% currentSphere
        currentSphere = NestedSphere[1:(k - i), ] /
          repmat.pns(matrix(sqrt(1 - NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), k - i, 1)

        if (nrow(currentSphere) == 3)
        {
          PNS$spherePNS = t(currentSphere)
        }

        if (nrow(currentSphere) == 2)
        {
          PNS$circlePNS = t(cur.proj)
        }

      }
    }
    if (nlast.small.sphere != 0)
    {
      for (i in (d - nlast.small.sphere):(d - 1))
      {
        sp = getSubSphere(x = currentSphere, geodesic = "small")
        center = sp$center
        r = sp$r
        res = acos(t(center) %*% currentSphere) - r
        orthaxis[[i]] = center
        dist[i] = r
        resmat[i, ] = res

        cur.proj = project.subsphere(x = currentSphere,
                                     center = center,
                                     r = r)

        NestedSphere = rotMat(center) %*% currentSphere
        currentSphere = NestedSphere[1:(k - i), ] /
          repmat.pns(matrix(sqrt(1 - NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), k - i, 1)


        if (nrow(currentSphere) == 3)
        {
          PNS$spherePNS = t(currentSphere)
        }

        if (nrow(currentSphere) == 2)
        {
          PNS$circlePNS = t(cur.proj)
        }

      }
    }
  } else {
    print("!!! Error from pns():")
    print("!!! sphere.type must be 'seq.test', 'small', 'great', 'BIC', or 'bi.sphere'")
    print("!!!   Terminating execution     ")
    return(NULL)
  }
  S1toRadian = atan2(currentSphere[2, ], currentSphere[1, ])
  meantheta = geodmeanS1(S1toRadian, mean.type = mean.type)$geodmean
  orthaxis[[d]] = meantheta
  resmat[d, ] = mod(S1toRadian - meantheta + pi, 2 * pi) - pi

  radii = 1
  for (i in 1:(d - 1))
  {
    radii = c(radii, prod(sin(dist[1:i])))
  }
  resmat = flipud0(repmat.pns(matrix(radii, ncol = 1), 1, n) * resmat)

  if (d > 1) {
    yy <- orthaxis[[d - 1]]
    xx <- c(-yy[2], yy[1] , yy[3])
    c1 <- Enorm(c(xx[1], xx[2], xx[3]) - c(-PNS$circlePNS[1, 2], PNS$circlePNS[1, 1], PNS$circlePNS[1, 3]))
    costheta <- 1 - c1^2 / 2
    angle <- (1:201) / (200) * 2 * pi
    centre <- xx * costheta
    A <- xx - centre
    B <- diag(3) - A %*% t(A) / Enorm(A) ** 2
    bv <- eigen(B)$vectors
    b1 <- bv[, 1]
    b2 <- bv[, 2]
    cc <- sin(acos(costheta)) * (cos(angle) %*% t(b1) + sin(angle) %*% t(b2)) + rep(1, times =
                                                                                      201) %*% t(centre)

  }
  PNS$scores = t(resmat)
  PNS$radii = radii
  PNS$pnscircle <- cbind(cbind(cc[, 2], -cc[, 1]) , cc[, 3])
  PNS$orthaxis = orthaxis
  PNS$dist = dist
  PNS$pvalues = pvalues
  PNS$ratio = ratio
  PNS$basisu = NULL
  PNS$mean = c(PNSe2s(matrix(0, d, 1), PNS))

  if (sphere.type == "seq.test")
  {
    PNS$sphere.type = "seq.test"
  } else if (sphere.type == "small") {
    PNS$sphere.type = "small"
  } else if (sphere.type == "great") {
    PNS$sphere.type = "great"
  } else if (sphere.type == "BIC") {
    PNS$sphere.type = "BIC"
  } else if (sphere.type == "bi.sphere") {
    PNS$sphere.type = "bi.sphere"
  }
  varPNS = apply(abs(resmat)^2, 1, sum) / n
  total = sum(varPNS)
  propPNS = varPNS / total * 100
  return(list(
    resmat = resmat,
    PNS = PNS,
    percent = propPNS
  ))
}

pns_biplot <- function(pns, varnames = rownames(q)) {
  pns1 <- pns
  nd <- dim(pns$resmat)[1] + 1
  grDevices::palette(grDevices::rainbow(nd))
  res1 <- cbind(c((20:(-20)) / 10 * stats::sd(pns1$resmat[1, ])) , matrix(0, 41, nd - 2))
  if (nd > 3) {
    res2 <- cbind(cbind(matrix(0, 41, 1) , c((20:(
      -20
    )) / 10 * stats::sd(pns1$resmat[2, ]))) , matrix(0, 41, nd - 3))
  } else{
    res2 <- cbind(cbind(matrix(0, 41, 1) , c((20:(
      -20
    )) / 10 * stats::sd(pns1$resmat[2, ]))))
  }
  aa1 <- PNSe2s(t(res1) , pns1$PNS) - pns1$PNS$mean
  aa2 <- PNSe2s(t(res2) , pns1$PNS) - pns1$PNS$mean
  plot(
    aa1[1, ],
    aa2[1, ],
    xlim = c(min(aa1), max(aa1)) ,
    type = "n",
    col = 2,
    ylim = c(min(aa2), max(aa2)) ,
    xlab = "PNS1",
    ylab = "PNS2"
  )
  for (i in 1:(nd)) {
    graphics::lines(aa1[i, ], aa2[i, ], col = i)
    graphics::arrows(aa1[i, 2], aa2[i, 2], aa1[i, 1], aa2[i, 1], col = i)
    graphics::text(aa1[i, 1], aa2[i, 1], varnames[i], col = i, cex = 1)
  }
  graphics::title("PNS biplot")
  grDevices::palette("default")
}


rotMat = function(b, a = NULL, alpha = NULL)
{
  if (is.matrix(b))
  {
    if (min(dim(b)) == 1)
    {
      b = c(b)
    } else {
      stop("Error: b should be a unit vector.")
    }
  }
  d = length(b)
  b = b / norm(b, type = "2")

  if (is.null(a) & is.null(alpha))
  {
    a = c(rep(0, d - 1), 1)
    alpha = acos(sum(a * b))
  } else if (!is.null(a) & is.null(alpha)) {
    alpha = acos(sum(a * b))
  } else if (is.null(a) & !is.null(alpha)) {
    a = c(rep(0, d - 1), 1)
  }
  if (abs(sum(a * b) - 1) < 1e-15)
  {
    rot = diag(d)
    return(rot)
  }
  if (abs(sum(a * b) + 1) < 1e-15)
  {
    rot = -diag(d)
    return(rot)
  }
  c = b - a * sum(a * b)
  c = c / norm(c, type = "2")
  A = a %*% t(c) - c %*% t(a)
  rot = diag(d) + sin(alpha) * A + (cos(alpha) - 1) * (a %*% t(a) + c %*% t(c))
  return(rot)
}


ExpNPd = function(x)
{
  if (is.vector(x))
  {
    x = as.matrix(x)
  }
  d = nrow(x)
  nv = sqrt(apply(x^2, 2, sum))
  Exppx = rbind(matrix(rep(sin(nv) / nv, d), nrow = d, byrow = T) * x, cos(nv))
  Exppx[, nv < 1e-16] = repmat.pns(matrix(c(rep(0, d), 1)), 1, sum(nv < 1e-16))
  return(Exppx)
}


LogNPd = function(x)
{
  n = ncol(x)
  d = nrow(x)
  scale = acos(x[d, ]) / sqrt(1 - x[d, ]^2)
  scale[is.nan(scale)] = 1
  Logpx = repmat.pns(t(scale), d - 1, 1) * x[-d, ]
  return(Logpx)
}


objfn = function(center, r, x)
{
  return(mean((acos(t(
    center
  ) %*% x) - r)^2))
}

project.subsphere = function(x, center, r)
{
  n = ncol(x)
  d = nrow(x)
  x.proj = matrix(NA, d, n)
  for (i in 1:n)
  {
    rho = acos(sum(x[, i] * center))
    x.proj[, i] = (sin(r) * x[, i] + sin(rho - r) * center) / sin(rho)
  }
  return(x.proj)
}

getSubSphere = function(x, geodesic = "small")
{
  svd.x = svd(x)
  initialCenter = svd.x$u[, ncol(svd.x$u)]
  c0 = initialCenter
  TOL = 1e-10
  cnt = 0
  err = 1
  n = ncol(x)
  d = nrow(x)
  Gnow = 1e+10
  while (err > TOL)
  {
    c0 = c0 / norm(c0, type = "2")
    rot = rotMat(c0)
    TpX = LogNPd(rot %*% x)
    fit = sphereFit(x = TpX,
                    initialCenter = rep(0, d - 1),
                    geodesic = geodesic)
    newCenterTp = fit$center
    r = fit$r
    if (r > pi)
    {
      r = pi / 2
      svd.TpX = svd(TpX)
      newCenterTp = svd.TpX$u[, ncol(svd.TpX$u)] * pi / 2
    }
    newCenter = ExpNPd(newCenterTp)
    center = solve(rot, newCenter)
    Gnext = objfn(center, r, x)
    err = abs(Gnow - Gnext)
    Gnow = Gnext
    c0 = center
    cnt = cnt + 1
    if (cnt > 30)
    {
      break
    }
  }
  i1save = list()
  i1save$Gnow = Gnow
  i1save$center = center
  i1save$r = r
  U = stats::princomp(t(x))$loadings[, ]
  initialCenter = U[, ncol(U)]
  c0 = initialCenter
  TOL = 1e-10
  cnt = 0
  err = 1
  n = ncol(x)
  d = nrow(x)
  Gnow = 1e+10

  while (err > TOL)
  {
    c0 = c0 / norm(c0, type = "2")
    rot = rotMat(c0)
    TpX = LogNPd(rot %*% x)
    fit = sphereFit(x = TpX,
                    initialCenter = rep(0, d - 1),
                    geodesic = geodesic)
    newCenterTp = fit$center
    r = fit$r
    if (r > pi)
    {
      r = pi / 2
      svd.TpX = svd(TpX)
      newCenterTp = svd.TpX$u[, ncol(svd.TpX$u)] * pi / 2
    }
    newCenter = ExpNPd(newCenterTp)
    center = solve(rot, newCenter)
    Gnext = objfn(center, r, x)
    err = abs(Gnow - Gnext)
    Gnow = Gnext
    c0 = center
    cnt = cnt + 1
    if (cnt > 30)
    {
      break
    }
  }
  if (i1save$Gnow == min(Gnow, i1save$Gnow))
  {
    center = i1save$center
    r = i1save$r
  }
  if (r > pi / 2)
  {
    center = -center
    r = pi - r
  }
  return(list(center = c(center), r = r))
}


LRTpval = function(resGREAT, resSMALL, n)
{
  chi2 = max(n * log(sum(resGREAT^2) / sum(resSMALL^2)), 0)
  return(stats::pchisq(
    q = chi2,
    df = 1,
    lower.tail = FALSE
  ))
}


vMFtest = function(x, R = 100)
{
  d = nrow(x)
  n = ncol(x)
  sumx = apply(x, 1, sum)
  rbar = norm(sumx, "2") / n
  muMLE = sumx / norm(sumx, "2")
  kappaMLE = (rbar * d - rbar^3) / (1 - rbar^2)
  sp = getSubSphere(x = x, geodesic = "small")
  center.s = sp$center
  r.s = sp$r
  radialdistances = acos(t(center.s) %*% x)
  xi_sample = mean(radialdistances) / stats::sd(radialdistances)
  xi_vec = rep(0, R)
  for (r in 1:R)
  {
    rdata = randvonMisesFisherm(d, n, kappaMLE)
    sp = getSubSphere(x = rdata, geodesic = "small")
    center.s = sp$center
    r.s = sp$r
    radialdistances = acos(t(center.s) %*% rdata)
    xi_vec[r] = mean(radialdistances) / stats::sd(radialdistances)
  }
  pvalue = mean(xi_vec > xi_sample)
  return(pvalue)
}


geodmeanS1 = function(theta, mean.type = "Frechet")
{
  n = length(theta)
  if (mean.type == "Frechet") {
    #kk candidates angles
    kk <- 1000
    meancandi = mod(mean(theta) + 2 * pi * (0:(kk - 1)) / kk, 2 * pi)
    theta = mod(theta, 2 * pi)
    geodvar = rep(0, kk)
    for (i in 1:kk)
    {
      v = meancandi[i]
      dist2 = apply(cbind((theta - v)^2, (theta - v + 2 * pi)^2, (v - theta + 2 * pi)^2), 1, min)
      geodvar[i] = sum(dist2)
    }
    m = min(geodvar)
    ind = which.min(geodvar)
    geodmean = mod(meancandi[ind], 2 * pi)
    geodvar = geodvar[ind] / n
  }
  if (mean.type == "Fisher") {
    mm <- atan2(mean(sin(theta)), mean(cos(theta)))
    geodmean <- mod(mm, 2 * pi)
    geodvar <- 1 - sqrt(mean(sin(theta)) ** 2 + mean(cos(theta)) ** 2)
  }
  return(list(geodmean = geodmean, geodvar = geodvar))
}


PNSe2s = function(resmat, PNS)
{
  dm = nrow(resmat)
  n = ncol(resmat)
  NSOrthaxis = rev(PNS$orthaxis[1:(dm - 1)])
  NSradius = flipud0(matrix(PNS$dist, ncol = 1))
  geodmean = PNS$orthaxis[[dm]]
  res = resmat / repmat.pns(flipud0(matrix(PNS$radii, ncol = 1)), 1, n)
  T = t(rotMat(NSOrthaxis[[1]])) %*%
    rbind(repmat.pns(sin(NSradius[1] + matrix(res[2, ], nrow = 1)), 2, 1) *
            rbind(cos(geodmean + res[1, ]), sin(geodmean + res[1, ])), cos(NSradius[1] + res[2, ]))
  if (dm > 2)
  {
    for (i in 1:(dm - 2))
    {
      T = t(rotMat(NSOrthaxis[[i + 1]])) %*%
        rbind(repmat.pns(sin(NSradius[i + 1] + matrix(res[i + 2, ], nrow = 1)), 2 + i, 1) * T, cos(NSradius[i + 1] + res[i + 2, ]))
    }
  }
  if (!is.null(PNS$basisu))
  {
    T = PNS$basisu %*% T
  }
  return(T)
}


PNSs2e = function(spheredata, PNS)
{
  if (nrow(spheredata) != length(PNS$mean))
  {
    cat(" Error from PNSs2e() \n")
    cat(" Dimensions of the sphere and PNS decomposition do not match")
    return(NULL)
  }
  if (!is.null(PNS$basisu))
  {
    spheredata = t(PNS$basisu) %*% spheredata
  }
  kk = nrow(spheredata)
  n = ncol(spheredata)
  Res = matrix(0, kk - 1, n)
  currentSphere = spheredata
  for (i in 1:(kk - 2))
  {
    v = PNS$orthaxis[[i]]
    r = PNS$dist[i]
    res = acos(t(v) %*% currentSphere) - r
    Res[i, ] = res
    NestedSphere = rotMat(v) %*% currentSphere
    currentSphere = as.matrix(NestedSphere[1:(kk - i), ]) /
      repmat.pns(matrix(sqrt(1 - NestedSphere[nrow(NestedSphere), ]^2), nrow = 1), kk - i, 1)
  }
  S1toRadian = atan2(currentSphere[2, ], currentSphere[1, ])
  devS1 = mod(S1toRadian - rev(PNS$orthaxis)[[1]] + pi, 2 * pi) - pi
  Res[kk - 1, ] = devS1
  EuclidData = flipud0(repmat.pns(PNS$radii, 1, n) * Res)
  return(EuclidData)
}


randvonMisesFisherm = function(m, n, kappa, mu = NULL)
{
  if (is.null(mu))
  {
    muflag = FALSE
  } else {
    muflag = TRUE
  }
  if (m < 2)
  {
    print("Message from randvonMisesFisherm(): dimension m must be > 2")
    print("Message from randvonMisesFisherm(): Set m to be 2")
    m = 2
  }
  if (kappa < 0)
  {
    print("Message from randvonMisesFisherm(): kappa must be >= 0")
    print("Message from randvonMisesFisherm(): Set kappa to be 0")
    kappa = 0
  }
  b = (-2 * kappa + sqrt(4 * kappa^2 + (m - 1)^2)) / (m - 1)
  x0 = (1 - b) / (1 + b)
  c = kappa * x0 + (m - 1) * log(1 - x0^2)
  nnow = n
  w = c()
  while (TRUE)
  {
    ntrial = max(round(nnow * 1.2), nnow + 10)
    Z = stats::rbeta(n = ntrial,
                     shape1 = (m - 1) / 2,
                     shape2 = (m - 1) / 2)
    U = stats::runif(ntrial)
    W = (1 - (1 + b) * Z) / (1 - (1 - b) * Z)

    indicator = kappa * W + (m - 1) * log(1 - x0 * W) - c >= log(U)
    if (sum(indicator) >= nnow)
    {
      w1 = W[indicator]
      w = c(w, w1[1:nnow])
      break
    } else {
      w = c(w, W[indicator])
      nnow = nnow - sum(indicator)
    }
  }
  V = UNIFORMdirections(m - 1, n)
  X = rbind(repmat.pns(sqrt(1 - matrix(w, nrow = 1)^2), m - 1, 1) * V, matrix(w, nrow = 1))
  if (muflag)
  {
    mu = mu / norm(mu, "2")
    X = t(rotMat(mu)) %*% X
  }
  return(X)
}


UNIFORMdirections = function(m, n)
{
  V = matrix(0, m, n)
  nr = matrix(stats::rnorm(m * n), nrow = m)
  for (i in 1:n)
  {
    while (TRUE)
    {
      ni = sum(nr[, i]^2)
      if (ni < 1e-10)
      {
        nr[, i] = stats::rnorm(m)
      } else {
        V[, i] = nr[, i] / sqrt(ni)
        break
      }
    }
  }
  return(V)
}

mod = function(x, y)
{
  return(x %% y)
}


repmat.pns = function(x, m, n)
{
  return(kronecker(matrix(1, m, n), x))
}


flipud0 = function(x)
{
  return(apply(x, 2, rev))
}


sphere.obj = function(center, x, is.greatCircle)
{
  di = sqrt(apply((x - repmat.pns(
    matrix(center, ncol = 1), 1, ncol(x)
  ))^2, 2, sum))
  if (is.greatCircle)
  {
    r = pi / 2
  } else {
    r = mean(di)
  }
  sum((di - r)^2)
}


sphere.res = function(center, x, is.greatCircle)
{
  center = c(center)
  xmc = x - center
  di = sqrt(apply(xmc^2, 2, sum))
  if (is.greatCircle)
  {
    r = pi / 2
  } else {
    r = mean(di)
  }
  (di - r)
}


sphere.jac = function(center, x, is.greatCircle)
{
  center = c(center)
  xmc = x - center
  di = sqrt(apply(xmc^2, 2, sum))
  di.vj = -xmc / repmat.pns(matrix(di, nrow = 1), length(center), 1)
  if (is.greatCircle)
  {
    c(t(di.vj))
  } else {
    r.vj = apply(di.vj, 1, mean)
    c(t(di.vj - repmat.pns(matrix(r.vj, ncol = 1), 1, ncol(x))))
  }
}


sphereFit = function(x,
                     initialCenter = NULL,
                     geodesic = "small")
{
  if (is.null(initialCenter))
  {
    initialCenter = apply(x, 1, mean)
  }
  op = minpack.lm::nls.lm(
    par = initialCenter,
    fn = sphere.res,
    jac = sphere.jac,
    x = x,
    is.greatCircle = ifelse(geodesic == "great", TRUE, FALSE),
    control = minpack.lm::nls.lm.control(maxiter = 1000)
  )
  center = stats::coef(op)
  di = sqrt(apply((x - repmat.pns(
    matrix(center, ncol = 1), 1, ncol(x)
  ))^2, 2, sum))
  if (geodesic == "great")
  {
    r = pi / 2
  } else {
    r = mean(di)
  }

  list(center = center, r = r)
}
