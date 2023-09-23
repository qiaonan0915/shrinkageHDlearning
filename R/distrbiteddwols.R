distrbiteddwols <- function(XX, Y, weights, site){ #XX <- X
  id.site <- as.character(site)
  id.site.uniq <- unique(id.site)
  SiX_s <- matrix(0, nrow = ncol(XX)+1, ncol = ncol(XX)+1)
  SiXY_s <- matrix(0, nrow = ncol(XX)+1, ncol = 1)
  for(ii in seq_along(id.site.uniq)){
    si = id.site.uniq[ii]
    wti = weights[id.site == si]
    XXi <- XX[id.site == si, ]
    Yi <- Y[id.site == si]
    SiX  = t(cbind(1,XXi)*wti) %*% cbind(1,XXi)
    SiXY = t(cbind(1,XXi)*wti) %*% Yi
    ni <- sum(id.site == si)
    SiX_s  <- SiX_s + SiX
    SiXY_s <- SiXY_s + SiXY
    # SiXYZ[[si]] <- list(SiX  = SiX, SiXZ = SiXZ, SiXY = SiXY,  SiZ  = SiZ, SiZY = SiZY, SiY  = SiY, ni = ni)
  }
  beta <- solve(SiX_s, SiXY_s)
  return(beta)
}
