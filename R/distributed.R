lmm.profile <- function(V, s2,             # var components para
                        pooled = F, reml = T,
                        Y, X, Z, id.site, weights=NULL,  # if pooled ipd is available
                        SiXYZ,             # if pooled ipd is not available, use summary stats
                        rcpp = F){
  ## lmm profile likelihood w.r.t. variance components
  if(pooled == T){
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    pz <- ncol(Z)
  }else{
    id.site.uniq <- names(SiXYZ)
    px <- ncol(SiXYZ[[1]]$SiX)
    pz <- ncol(SiXYZ[[1]]$SiXZ)
  }

  # allows assuming the same residual var across sites
  if(length(s2) == 1) {
    s2 <- rep(s2, length(id.site.uniq))
  }

  if(rcpp == F){
    lpterm1 = lpterm2 = remlterm = 0
    remlterm.i = rep(NA, length(id.site.uniq))

    bterm1 <- matrix(0, px, px)
    bterm2 <- rep(0, px)
    Vinv <- solve(V)
    Wi <- ui <- varui <- varui_post <- list()  # save this for each subject for BLUP

    for(ii in seq_along(id.site.uniq)){
      si <- id.site.uniq[ii]
      if(pooled == T){
        SiXYZ <- lmm.get.summary(Y, X, Z, weights, id.site)
      } #else{
      SiX  <- SiXYZ[[si]]$SiX
      SiXZ <- SiXYZ[[si]]$SiXZ
      SiXY <- SiXYZ[[si]]$SiXY
      SiZ <- SiXYZ[[si]]$SiZ
      SiZY <- SiXYZ[[si]]$SiZY
      SiY  <- SiXYZ[[si]]$SiY
      ni <- SiXYZ[[si]]$ni
      # }

      s2i <- s2[ii]  # sigma_i^2 #随机项方差

      tmp <- log(det(diag(1, pz) + SiZ %*% V / s2i))

      if(is.na(tmp)) cat(diag(V), '...', s2i, '...', V[1,], '\n')
      # logdet <- ni * log(s2i) + log(det(diag(1, pz) + SiZ %*% V / s2i))
      logdet <- ni * log(s2i) + log(det(diag(1, pz) + SiZ %*% V / s2i))

      # log(max(1e-14, det(diag(1,pz)+SiZ%*%V/s2i)))
      lpterm1 <- lpterm1 + logdet

      Wi[[ii]] <- solve(s2i * Vinv + SiZ)
      bterm1 <- bterm1 + (SiX - SiXZ %*% Wi[[ii]] %*% t(SiXZ)) / s2i
      bterm2 <- bterm2 + (SiXY - SiXZ %*% Wi[[ii]] %*% SiZY) / s2i
      lpterm2 <- lpterm2 + (SiY - t(SiZY) %*% Wi[[ii]] %*% SiZY) / s2i

      if(reml == T){
        # tmp = log(det((SiX-SiXZ%*%Wi[[ii]]%*%t(SiXZ))/s2i))
        # if(is.na(tmp)) cat(Wi[[ii]], '...', s2i, '\n')
        remlterm.i[ii] <- log(abs(det((SiX - SiXZ %*% Wi[[ii]] %*% t(SiXZ))/s2i))) # abs() to make sure positivity
        # remlterm <- remlterm + log(det((SiX - SiXZ %*% Wi[[ii]] %*% t(SiXZ)) / s2i))
      }
    }

    b <- solve(bterm1, bterm2)
    if(reml == T){
      remlterm = sum(remlterm.i[is.finite(remlterm.i)])
      lp <- - (lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b + remlterm) / 2
    }else{
      lp <- - (lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b) / 2
    }

    lk <- -(lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b) / 2 # To be used in calculating mAIC

    # Loop to estimate ui
    # ui = V Z_i^{T} Sigma_{i}^{-1} (Y_i - X_i \beta)
    for(ii in seq_along(id.site.uniq)){ # re-call this loop
      si <- id.site.uniq[ii]
      if(pooled == T){
        Xi <- X[id.site == si, ]
        Zi <- Z[id.site == si, ]
        Yi <- Y[id.site == si]
        SiX  <- t(Xi) %*% Xi
        SiXZ <- t(Xi) %*% Zi
        SiXY <- t(Xi) %*% Yi
        SiZ  <- t(Zi) %*% Zi
        SiZY <- t(Zi) %*% Yi
        SiY  <- sum(Yi ^ 2)
        ni <- sum(id.site == si)
      }else{
        SiX  <- SiXYZ[[si]]$SiX
        SiXZ <- SiXYZ[[si]]$SiXZ
        SiXY <- SiXYZ[[si]]$SiXY
        SiZ <- SiXYZ[[si]]$SiZ
        SiZY <- SiXYZ[[si]]$SiZY
        SiY <- SiXYZ[[si]]$SiY
        ni <- SiXYZ[[si]]$ni
      }

      s2i <- s2[ii]  # sigma_i^2
      uiterm1 <- (SiZY - SiZ %*% Wi[[ii]] %*% SiZY) / s2i
      uiterm2 <-  ((t(SiXZ) - SiZ %*% Wi[[ii]] %*% t(SiXZ)) / s2i) %*% as.matrix(b)
      ui[[ii]] <- V %*% as.numeric(uiterm1 - uiterm2)

      vterm1 <- V %*% ((SiZ - SiZ %*% Wi[[ii]] %*% SiZ) / s2i) %*% V
      vterm2 <- V %*% ((t(SiXZ) - SiZ %*% Wi[[ii]] %*% t(SiXZ)) / s2i)
      # varui[[ii]] <- V - vterm1 + (vterm2 %*% t(vterm2) / lpterm1)
      # varui_post[[ii]] <- vterm1 - (vterm2 %*% t(vterm2) / lpterm1)
      # 20210112:
      varui[[ii]] <- V - vterm1 + (vterm2 %*% solve(bterm1, t(vterm2)))
      varui_post[[ii]] <- vterm1 - (vterm2 %*% solve(bterm1, t(vterm2)))
    }


    res <- list(lp = lp, b = b, ui = ui, lk = lk, varui = varui, varui_post = varui_post,
                allterms = list(lpterm1 = lpterm1,
                                lpterm2 = lpterm2,
                                remlterm = remlterm,
                                remlterm.i = remlterm.i,
                                bterm1 = bterm1,
                                bterm2 = bterm2))
  } else{ ## CHECK THIS
    res <- NULL # LMM_Profile(SiXYZ = SiXYZ, V = V, s2 = s2, reml = reml)
  }
  # SiY - t(SiZY)%*%Wi%*%SiZY - 2*(t(SiXY)-t(SiZY)%*%Wi%*%t(SiXZ))%*%b + t(b)%*%(SiX-SiXZ%*%Wi%*%t(SiXZ))%*%b
  return(res)
}


lmm.profile0 <- function(V,                # var components para, = original V / s2    # s2,
                         pooled = F, reml = T,
                         Y, X, Z, id.site, weights=NULL,  # if pooled ipd is available
                         SiXYZ,             # if pooled ipd is not available, use summary stats
                         rcpp = F){
  ## further profile out the residual var s2, used if common.s2=T (deemed as usual situation)
  if(pooled == T){
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    pz <- ncol(Z)
  }else{
    id.site.uniq <- names(SiXYZ)
    px <- ncol(SiXYZ[[1]]$SiX)
    pz <- ncol(SiXYZ[[1]]$SiXZ)
  }

  if(rcpp == F){
    lpterm1 = lpterm2 = remlterm = 0
    bterm1 <- matrix(0, px, px)  # sum_i Xi' \Sigma_i^-1 Xi
    bterm2 <- rep(0, px)         # sum_i Xi' \Sigma_i^-1 Yi
    Vinv <- solve(V)
    Wi <- ui <- varui <- varui_post <- list()  # save this for each subject for BLUP

    N <- 0
    for(ii in seq_along(id.site.uniq)){
      si <- id.site.uniq[ii]
      if(pooled == T){
        SiXYZ <- lmm.get.summary(Y, X, Z, weights, id.site)
      } # else{
      SiX  <- SiXYZ[[si]]$SiX
      SiXZ <- SiXYZ[[si]]$SiXZ
      SiXY <- SiXYZ[[si]]$SiXY
      SiZ <- SiXYZ[[si]]$SiZ
      SiZY <- SiXYZ[[si]]$SiZY
      SiY  <- SiXYZ[[si]]$SiY
      ni <- SiXYZ[[si]]$ni
      # }
      N <- N + ni

      # tmp <- log(det(diag(1, pz) + SiZ %*% V))  # improve
      logdet <- log(det(diag(1, pz) + SiZ %*% V))   # -log(det(diag(wti))) omitted as wti are the fixed weights
      lpterm1 <- lpterm1 + logdet

      Wi[[ii]] <- solve(Vinv + SiZ)
      bterm1 <- bterm1 + (SiX - SiXZ %*% Wi[[ii]] %*% t(SiXZ)) # / s2i
      bterm2 <- bterm2 + (SiXY - SiXZ %*% Wi[[ii]] %*% SiZY) # / s2i
      lpterm2 <- lpterm2 + (SiY - t(SiZY) %*% Wi[[ii]] %*% SiZY) #/ s2i
    }

    b <- solve(bterm1, bterm2) #beta
    # quadratic term: = (Yi-Xi*b)'\Sigma_i^-1 (Yi-Xi*b)
    qterm <- as.numeric(lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b )
    if(reml == T){
      remlterm = log(det(bterm1))
      s2 <- qterm / (N-px)
      # lp <- - (lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b + remlterm) / 2
      lp <- - (lpterm1 + (1 + log(qterm * 2 * pi / (N - px))) * (N - px)) / 2     # Bates2015JSS eq(42)
    }else{
      s2 <- qterm / N
      lp <- - (lpterm1 + (1 + log(qterm * 2 * pi / N)) * N) / 2               # Bates2015JSS eq(35)
    }

    # lk <- -(lpterm1 + lpterm2 - 2 * sum(bterm2 * b) + t(b) %*% bterm1 %*% b) / 2 # To be used in calculating mAIC
    lk <- - (lpterm1 + (1+log(qterm*2*pi/N))*N) / 2

    # Loop to estimate ui
    # ui = V Z_i^{T} Sigma_{i}^{-1} (Y_i - X_i \beta)
    for(ii in seq_along(id.site.uniq)){ # re-call this loop
      si <- id.site.uniq[ii]

      SiX  <- SiXYZ[[si]]$SiX
      SiXZ <- SiXYZ[[si]]$SiXZ
      SiXY <- SiXYZ[[si]]$SiXY
      SiZ <- SiXYZ[[si]]$SiZ
      SiZY <- SiXYZ[[si]]$SiZY
      SiY <- SiXYZ[[si]]$SiY
      ni <- SiXYZ[[si]]$ni

      s2i <- s2   # [ii]  # sigma_i^2
      uiterm1 <- (SiZY - SiZ %*% Wi[[ii]] %*% SiZY) / s2i
      uiterm2 <-  ((t(SiXZ) - SiZ %*% Wi[[ii]] %*% t(SiXZ)) / s2i) %*% as.matrix(b)
      ui[[ii]] <- V %*% as.numeric(uiterm1 - uiterm2) * s2i #

      vterm1 <- V %*% ((SiZ - SiZ %*% Wi[[ii]] %*% SiZ) / s2i) %*% V * (s2i ^ 2)    #
      vterm2 <- V %*% ((t(SiXZ) - SiZ %*% Wi[[ii]] %*% t(SiXZ)) / s2i)  * s2i   #
      varui[[ii]] <- (V * s2i) - vterm1 + (vterm2 %*% solve(bterm1, t(vterm2)))     # 20210111
      varui_post[[ii]] <- vterm1 - (vterm2 %*% solve(bterm1, t(vterm2)))            # 20210111
    }

    res <- list(lp = lp, b = b,
                s2 = s2,
                ui = ui, lk = lk, varui = varui, varui_post = varui_post,
                allterms = list(lpterm1 = lpterm1,
                                lpterm2 = lpterm2,
                                qterm = qterm,
                                remlterm = remlterm,
                                # remlterm.i = remlterm.i,
                                bterm1 = bterm1,
                                bterm2 = bterm2))
  } else{ ## CHECK THIS
    res <- NULL # LMM_Profile(SiXYZ = SiXYZ, V = V, s2 = s2, reml = reml)
  }
  return(res)
}  #除随机项方差分量后的目标函数



lmm.fit <- function(Y = NULL, X = NULL, Z = NULL, id.site = NULL, weights = NULL,
                    pooled = F, reml = T,
                    common.s2 = T,      # common residual var across sites
                    SiXYZ = list(),
                    corstr = 'independence', # 'exchangeable', 'ar1', 'unstructured'),
                    mypar.init = NULL,
                    hessian = F,
                    verbose){
  ## fit lmm distributed
  if(pooled == T){
    id.site.uniq <- unique(id.site)
    px <- ncol(X)
    pz <- ncol(Z)
    K <- length(id.site.uniq)
    SiXYZ <- lmm.get.summary(Y, X, Z, weights, id.site)
  }else{
    id.site.uniq <- names(SiXYZ)
    px <- ncol(SiXYZ[[1]]$SiX)
    pz <- ncol(SiXYZ[[1]]$SiXZ)
    K <- length(SiXYZ)
  }

  ## 20200629: now further profile out s2 if common.s2=T
  if(common.s2 == T){
    ns <- 1 # number of s2 para
    fn <- function(mypar){
      if(corstr == 'independence'){
        V <- diag(mypar[1 : pz], pz)
        # V <- diag(exp(mypar[1 : pz]), pz)
        # s2 <- exp(mypar[-c(1 : pz)])
      }else if(corstr == 'exchangeable'){
        V = diag(sqrt(mypar[1 : pz])) %*% (matrix(mypar[pz + 1], pz, pz) + diag(1 - mypar[pz + 1], pz)) %*% diag(sqrt(mypar[1 : pz]))
        # s2 <- mypar[-c(1 : (pz + 1))]
      }else if(corstr == 'unstructured'){
        V = matrix(0, pz, pz)
        diag(V) <- mypar[1 : pz]
        V[lower.tri(V)] <- V[upper.tri(V)] <- mypar[(pz + 1) : (pz * (pz + 1) / 2)]
        # s2 <- mypar[-c(1 : (pz*(pz+1)/2))]
      }
      return(-lmm.profile0(V, pooled=F, reml, Y, X, Z, id.site, weights, SiXYZ)$lp)
    }

    if(is.null(mypar.init)){
      if(corstr == 'independence'){
        mypar.init <- rep(0.5, pz)
        # mypar.init <- log(c(rep(0.5, pz)))
      }else if(corstr == 'exchangeable'){
        mypar.init <- c(rep(0.5, pz), 0.1 )
      }else if(corstr == 'unstructured'){
        mypar.init <- c(rep(0.5, pz), rep(0.1, pz * (pz - 1) / 2) )
      }
      cat('default mypar.init (var comp) = ', mypar.init, '\n')
    }

    # res <- optim(mypar.init, fn, hessian = hessian)
    res <- minqa::bobyqa(mypar.init, fn, lower=rep(1e-6, pz), control=list(maxfun=1e5))
    # res <- nloptwrap(mypar.init, fn, lower=rep(1e-6, pz), upper = rep(1e6, pz), control=list(maxfun=1e5))

    mypar <- res$par
    if(corstr == 'independence'){
      V <- diag(mypar[1 : pz], pz)
      # V <- diag(exp(mypar[1 : pz]), pz)
      # s2 <- exp(mypar[- c(1 : pz)])
    }else if(corstr == 'exchangeable'){
      V <- diag(sqrt(mypar[1 : pz])) %*% (matrix(mypar[pz + 1], pz, pz) + diag(1 - mypar[pz + 1], pz)) %*% diag(sqrt(mypar[1 : pz]))
      # s2 <- mypar[- c(1 : (pz + 1))]
    }else if(corstr == 'unstructured'){
      V <- matrix(0, pz, pz)
      diag(V) <- mypar[1 : pz]
      V[lower.tri(V)] <- V[upper.tri(V)] <- mypar[(pz + 1) : (pz * (pz + 1) / 2)]
      # s2 <- mypar[-c(1 : (pz * (pz + 1) / 2))]
      # error('corstr=="unstructured" not yet implemented')
    }

    res.profile <- lmm.profile0(V = V, pooled=F, reml, Y, X, Z, id.site, weights, SiXYZ)
    s2 <- res.profile$s2
    V <- V * s2             # scale back
  }else{  # if common.s2=F, can't profile out s2 vector
    ns <- K
    fn <- function(mypar){
      if(corstr == 'independence'){
        V <- diag(mypar[1 : pz], pz)
        s2 <- exp(mypar[-c(1 : pz)])
      }else if(corstr == 'exchangeable'){
        V = diag(sqrt(mypar[1 : pz])) %*% (matrix(mypar[pz + 1], pz, pz) + diag(1 - mypar[pz + 1], pz)) %*% diag(sqrt(mypar[1 : pz]))
        s2 <- mypar[-c(1 : (pz + 1))]
      }else if(corstr == 'unstructured'){
        V = matrix(0, pz, pz)
        diag(V) <- mypar[1 : pz]
        V[lower.tri(V)] <- V[upper.tri(V)] <- mypar[(pz + 1) : (pz * (pz + 1) / 2)]
        s2 <- mypar[-c(1 : (pz*(pz+1)/2))]
        # error('corstr=="unstructured" not yet implemented')
      }
      return(-lmm.profile(V, s2, pooled=F, reml, Y, X, Z, id.site, weights, SiXYZ)$lp)
    }

    if(is.null(mypar.init)){
      if(corstr == 'independence'){
        mypar.init <- c(rep(0.5, pz), rep(0.5, ns))
      }else if(corstr == 'exchangeable'){
        mypar.init <- c(rep(0.5, pz), 0.1, rep(0.5, ns))
      }else if(corstr == 'unstructured'){
        mypar.init <- c(rep(0.5, pz), rep(0.1, pz * (pz - 1) / 2), rep(0.5, ns))
      }
      cat('default mypar.init (var comp) = ', mypar.init, '\n')
    }

    # res <- optim(mypar.init, fn, hessian = hessian)
    res <- minqa::bobyqa(mypar.init, fn, lower=rep(1e-6, length(mypar.init)), control=list(maxfun=1e5))
    # res <- nloptwrap(mypar.init, fn, lower=rep(1e-6, length(mypar.init)),
    #                  upper =rep(1e6, length(mypar.init)),  control=list(maxfun=1e5))

    mypar <- res$par
    if(corstr == 'independence'){
      V <- diag(mypar[1 : pz], pz)
      s2 <- mypar[- c(1 : pz)]
    }else if(corstr == 'exchangeable'){
      V <- diag(sqrt(mypar[1 : pz])) %*% (matrix(mypar[pz + 1], pz, pz) + diag(1 - mypar[pz + 1], pz)) %*% diag(sqrt(mypar[1 : pz]))
      s2 <- mypar[- c(1 : (pz + 1))]
    }else if(corstr == 'unstructured'){
      V <- matrix(0, pz, pz)
      diag(V) <- mypar[1 : pz]
      V[lower.tri(V)] <- V[upper.tri(V)] <- mypar[(pz + 1) : (pz * (pz + 1) / 2)]
      s2 <- mypar[-c(1 : (pz * (pz + 1) / 2))]
      # error('corstr=="unstructured" not yet implemented')
    }

    res.profile <- lmm.profile(V = V, s2 = s2, pooled, reml, Y, X, Z, id.site, SiXYZ)
  }




  # # ## New added
  # # ## Inference (Wald test statistic)
  # # vd <- diag(solve(res.profile$allterms$bterm1))
  # # if(common.s2==T)  vd <- diag(solve(res.profile$allterms$bterm1 / s2)); vd_matrix <-  solve(res.profile$allterms$bterm1 / s2) # scale back
  # # wald <- res.profile$b / sqrt(vd)
  # #
  # # ## 95% CI for fixed effects
  # # lb <- res.profile$b -  1.96 * sqrt(vd)
  # # ub <- res.profile$b +  1.96 * sqrt(vd)
  # #
  # # ## Marginal AIC: OKAY to use if the main interest is to model fixed population effects with a reasonable correlation structur (Kneib & Greven (2010))
  # # ###边际AIC：如果主要兴趣是用合理的相关性结构对固定种群效应进行建模，那么可以使用（Knib&Greven（2010））
  # # mAIC <- 2 * res.profile$lk + 2 * (px + (length(mypar) - ns))
  #
  # ## Conditional AIC: Vaida & Blanchard (2005) expression details in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2572765/pdf/nihms62635.pdf
  # ## However a better approach should be in Kneib & Steven (2010) Appendix:B https://arxiv.org/pdf/1803.05664.pdf
  # ## assuming V and sigma^2 are (UN)KNOWN
  # # if(pooled == T & common.s2 == T){
  # #   c2 <- diag(s2, ncol(V)) * V
  # #
  # #   #if(common.s2 == T){
  # #   #  c2 <- bdiag(lapply(seq_len(K), function(kk) diag(s2[1], ncol(V)) * V))  #
  # #   #} else if(common.s2 == F){
  # #   #  c2 <- bdiag(lapply(seq_len(K), function(kk) diag(s2[kk], ncol(V)) * V)) # block diagonal
  # #   #}
  # #
  # #   LLinv <- solve(c2, diag(1, ncol(c2)))
  # #   c3 <- eigen(LLinv)
  # #   Lambda <- c3$vectors %*% diag(sqrt(c3$values))
  # #
  # #   c11 <- cbind(X, Z)
  # #   c12 <- cbind(matrix(0, ncol = px, nrow = nrow(Lambda)), Lambda)
  # #   # dim(c11); dim(c12)
  # #   M <- rbind(c11, c12)
  # #
  # #   MtM <- solve(t(M) %*% M)
  # #   # H1 <- c11 %*% MtM %*% t(c11)
  # #   # trace <- sum(diag(H1))
  # #   trace <- sum(diag(MtM %*% t(c11) %*% c11))
  # #   cAIC_vb <- 2 * res.profile$lk + 2 * trace
  # # } else{
  # cAIC_vb = NULL
  # # }

  ## Prediction
  uihat <- as.matrix(do.call(rbind, lapply(seq_len(K), function(kk) {

    ll <- length(which(id.site == id.site.uniq[kk]))
    dm <- matrix(1, nrow = ll, ncol = pz)

    sweep(dm, MARGIN = 2, res.profile$ui[[kk]], `*`)

  })))

  uihat_subj <- as.matrix(do.call(rbind, lapply(seq_len(K), function(kk) {
    t(res.profile$ui[[kk]])
  })))

  if(pooled == T){
    Yhat <- X %*% as.matrix(res.profile$b) # population-level
    Yihat <- Yhat +  rowSums(Z * uihat) # subject-level
  } else {
    Yhat <- Yihat <- NULL
  }

  return(list(b = res.profile$b,
              b.sd = sqrt(vd),     # sd of fixed effect est
              b.matrix = vd_matrix, #beta covariance
              s2 = s2,
              # wald = wald,   # Wald-test statistic
              # lb = lb,       # lower-bound
              # ub = ub,       # uppper-bound
              XvX = res.profile$allterms$bterm1, # X^{T}V^{-1}X
              ui = res.profile$ui, # BLUP of random effects
              uiM = uihat_subj,
              varui = res.profile$varui,  # Variance (based on prediction error)
              varui_post = res.profile$varui_post, # posterior
              Yhat = Yhat,         # population-level prediction (WOT random effects) XB
              Yihat = Yihat,       # subject-specific prediction XB + Zu
              mAIC = mAIC,
              cAIC_vb = cAIC_vb,
              V = V,
              s2 = s2,
              res = res, res.profile = res.profile))
}






#' @keywords internal
lmm.get.summary <- function(Y = NULL, X = NULL, Z = NULL,
                            weights = NULL, id.site = NULL){
  ## get summary stats from each site for distributed lmm
  # 20201203: incorporate weight (wt) in LMM for distributed PQL, all the summary stats are Xi^TWiX_i, Xi^TWiZ_i, Xi^TWiY_i, etc

  if(is.null(weights)) weights <- rep(1, length(Y))
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  id.site <- as.character(id.site)
  id.site.uniq <- unique(id.site)
  px <- ncol(X)
  pz <- ncol(Z)

  SiXYZ <- list()
  for(ii in seq_along(id.site.uniq)){
    si = id.site.uniq[ii]
    wti = weights[id.site == si]
    Xi <- X[id.site == si, ]
    Zi <- Z[id.site == si, ]
    Yi <- Y[id.site == si]
    # if(any(apply(Xi[,-1], 2, function(a)length(unique(a)))==1))
    #   warning(paste0('singular X in site #', ii, ' detected!'))
    # if(any(apply(Zi[,-1], 2, function(a)length(unique(a)))==1))
    #   warning(paste0('singular Z in site #', ii, ' detected!'))

    SiX  = t(Xi*wti) %*% Xi
    SiXZ = t(Xi*wti) %*% Zi
    SiXY = t(Xi*wti) %*% Yi
    SiZ  = t(Zi*wti) %*% Zi
    SiZY = t(Zi*wti) %*% Yi
    SiY  = sum(Yi ^ 2 *wti)
    ni <- sum(id.site == si)
    SiXYZ[[si]] <- list(SiX  = SiX, SiXZ = SiXZ, SiXY = SiXY,
                        SiZ  = SiZ, SiZY = SiZY, SiY  = SiY, ni = ni)
  }

  return(SiXYZ)
}
