# ---- numerical helpers ----
trapz_cumint <- function(t, y) {
  # cumulative integral via trapezoid rule
  n <- length(t)
  out <- numeric(n)
  dt <- diff(t)
  out[-1] <- cumsum(dt * (y[-1] + y[-n]) / 2)
  out
}

sample_discrete <- function(x, w) {
  # sample x with probabilities proportional to w (w must be >=0, not all zero)
  w <- pmax(w, 0)
  s <- sum(w)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  sample(x, size = 1, prob = w / s)
}
 

draw_params_input <- function(parameters, n=1) {
  if(parameters$fixed_param ==FALSE){     
   L =  list(
      T_P = parameters$T_P_star * exp(rnorm(n, sd = parameters$T_P_eta)),
      T_U = parameters$T_U_star * exp(rnorm(n, sd = parameters$T_U_eta)),
      V_P = parameters$V_P_star * exp(rnorm(n, sd = parameters$V_P_eta)),
      m   = parameters$m
    ) 
  }
  
  if(parameters$fixed_param ==TRUE){     
    L =  list(
      T_P = parameters$T_P_star * exp(rnorm(n, sd = 0)),
      T_U = parameters$T_U_star * exp(rnorm(n, sd = 0)),
      V_P = parameters$V_P_star * exp(rnorm(n, sd = 0)),
      m   = parameters$m
    ) 
  }
  
  return(L)
}


compute_q_and_P <- function(T_P, T_U, V_P, m, t_grid) {
  vl <- model_gen_triangle(T_P = T_P, T_U = T_U, V_P = V_P, t = t_grid)
  Vt <- vl$V
  
  ht <- m * sqrt(pmax(Vt, 0))
  Ht <- trapz_cumint(vl$t, ht)
  
  qt <- ht * exp(-Ht)
  
  H_end <- Ht[length(Ht)]
  P <- 1 - exp(-H_end)  # equals integral_0^T qt dt (approx), if T is large enough
  
  list(vl = vl, h = ht, H = Ht, q = qt, P = P)
}



draw_time_conditional <- function(t, q, P) {
  # conditional on a transmission event
  if (!is.finite(P) || P <= 0) return(NA_real_)
  sample_discrete(t, q)  # sampling proportional to q(t); normalization cancels out
}

draw_time_weighted <- function(t, q, P) {
  # unconditional mixture: event happens with prob P, then t ~ q/P
  if (!is.finite(P) || P <= 0) return(list(event = 0L, t = NA_real_))
  event <- rbinom(1, 1, prob = pmin(pmax(P, 0), 1))
  if (event == 0) return(list(event = 0L, t = NA_real_))
  list(event = 1L, t = sample_discrete(t, q))
}

simulate_transmissions_input <- function(parameters,
                                         N = 2000,
                                         t_grid = seq(0, 35, by = 0.05),
                                         keep_curves = FALSE) {
  
  results <- vector("list", N)
  curves  <- if (keep_curves) vector("list", N) else NULL
  
  for (i in seq_len(N)) {
    params <- draw_params_input(parameters)
    
    qp <- compute_q_and_P(
      T_P = params$T_P,
      T_U = params$T_U,
      V_P = params$V_P,
      m   = params$m,
      t_grid = t_grid
    )
    
    t_cond <- draw_time_conditional(qp$vl$t, qp$q, qp$P)
    wdraw  <- draw_time_weighted(qp$vl$t, qp$q, qp$P)
    
    results[[i]] <- data.frame(
      sim = i,
      T_P = params$T_P,
      T_U = params$T_U,
      V_P = params$V_P,
      m   = params$m,
      P   = qp$P,
      t_conditional = t_cond,
      event_weighted = wdraw$event,
      t_weighted = wdraw$t
    )
    
    if (keep_curves) {
      curves[[i]] <- transform(qp$vl,
                               sim = i,
                               h = qp$h,
                               H = qp$H,
                               q = qp$q
      )
    }
  }
  
  out <- list(summary = do.call(rbind, results))
  if (keep_curves) out$curves <- do.call(rbind, curves)
  out
}
 