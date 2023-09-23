emexoutcome <- function(R, a, preda, W){ #W is propensity score
  v_fenzi <- sum(R*(a == preda)/W)
  v_fenmu <- sum((a == preda)/W)
  v <- v_fenzi/v_fenmu
  return(v)
}
