#' Approximate Gradient (Jacobian for scalar functions)
#' @param f A vector-valued function f(x)
#' @param x A numeric vector at which to evaluate the gradient
#' @return Numeric vector of partial derivatives
numeric_gradient <- function(f, x, h1 = 1e-8) {
  fx <- f(x)
  grad <- numeric(length(fx))
  for (i in seq_along(fx)) {
    h <- numeric(length(x))
    h[i] <- h1
    grad[i] <- (f(x + h)[i] - fx[i]) / h1
  }
  return(grad)
}

#' Jacobian Matrix using Finite Differences
#' @param f A vector-valued function f(x)
#' @param x A numeric vector
#' @return Jacobian matrix of f at x
jacobian <- function(f, x, h1 = 1e-8) {
  fx <- f(x)
  J <- matrix(0, nrow = length(fx), ncol = length(x))
  for (j in seq_along(x)) {
    h <- numeric(length(x))
    h[j] <- h1
    J[, j] <- (f(x + h) - fx) / h1
  }
  return(J)
}

#' Solve System of Nonlinear Equations using Newton-Raphson Method
#' @param f A vector-valued function f(x)
#' @param x0 Initial guess (numeric vector)
#' @param tol Convergence tolerance (default 1e-6)
#' @param max_iter Maximum number of iterations (default 100)
#' @return Solution vector
slfn <- function(f, x0, tol = 1e-6, max_iter = 100) {
  for (i in 1:max_iter) {
    J <- jacobian(f, x0)
    fx <- f(x0)
    delta <- tryCatch(solve(J, fx), error = function(e) stop("Jacobian is singular or ill-conditioned."))
    x1 <- x0 - delta
    if (sqrt(sum((x1 - x0)^2)) < tol) {
      return(x1)
    }
    x0 <- x1
  }
  stop("Did not converge within the maximum number of iterations.")
}

lm_binom <- function(y, x, trials = 2) {
  x <- as.matrix(cbind(1, x))  # Add intercept
  n <- length(y)
  
  # Vectorized score function (gradient of log-likelihood)
  score_fn <- function(theta) {
    eta <- x %*% theta
    mu <- plogis(eta)  # safer and faster than exp / (1 + exp)
    t(x) %*% (y - trials * mu)
  }
  
  # Solve for MLEs using Newton-Raphson
  theta <- slfn(score_fn, rep(0, ncol(x)))
  
  # Predict probabilities and compute log-likelihood
  eta <- x %*% theta
  p <- plogis(eta)
  loglik <- sum(dbinom(y, size = trials, prob = p, log = TRUE))
  
  list(estimate = theta, prob = p, likelihood = loglik)
}



binomial_lrt <- function(gt, y, cv = NULL) {
  x_full<-cbind(y,cv)
  x1<-x_full[,-2]       #Discarding second trait
  x2<-x_full[,-1]       #Discarding first trait
  
  m1 <- lm_binom(gt, x1)
  m2 <- lm_binom(gt, x2)
  m3 <- lm_binom(gt, x_full)
  
  # Pick the stronger null
  if (m1$likelihood > m2$likelihood) {
    lrt <- -2 * (m1$likelihood - m3$likelihood)
    p_g <- m1$prob
  } else {
    lrt <- -2 * (m2$likelihood - m3$likelihood)
    p_g <- m2$prob
  }
  
  list(lrt = lrt, null_prob = p_g, null1 = m1$estimate, null2 = m2$estimate, full = m3$estimate)
}


hoprr <- function(gt, y, cv = NULL, n_perm = 1000, seed = 123) {
  set.seed(seed)
  sm_size <- length(gt)
  test <- binomial_lrt(gt, y, cv)
  lrt_obs <- test$lrt
  p_g <- test$null_prob
  
  # Initialize
  lrt_perm <- c()
  total_attempts <- 0
  
  # Keep collecting until n_perm valid permutations
  while (length(lrt_perm) < n_perm) {
    n_remaining <- n_perm - length(lrt_perm)
    batch_size <- min(500, n_remaining * 2)  # Oversample to account for failures
    
    lrt_batch <- numeric(batch_size)
    
    for (i in 1:batch_size) {
      total_attempts <- total_attempts + 1
      gt_perm <- rbinom(sm_size, size = 2, prob = p_g)
      lrt_batch[i] <- tryCatch({
        binomial_lrt(gt_perm, y, cv)$lrt
      }, error = function(e) NA)
    }
    
    # Keep only successful results
    lrt_perm <- c(lrt_perm, lrt_batch[!is.na(lrt_batch)])
  }
  
  # Trim to exactly n_perm
  lrt_perm <- lrt_perm[1:n_perm]
  p_val <- mean(lrt_perm > lrt_obs)
  decision <- ifelse(p_val <= 0.05, "Significant", "Insignificant")
  
  list(
    sample_size = sm_size,
    est = test$full,
    lrt = lrt_obs,
    p_value = p_val,
    decision = decision,
    valid_permutations = n_perm,
    total_attempted = total_attempts
  )
}

plafiu <- function(gt, y, cv = NULL){
  
  fwd_adj1 <- summary(lm(y[,1] ~ gt + cbind(y[,2], cv)))
  fwd_adj2 <- summary(lm(y[,2] ~ gt + cbind(y[,1], cv)))
  
  fwd_adj_e1 <- fwd_adj1$coefficients[2,1]
  fwd_adj_e2 <- fwd_adj2$coefficients[2,1]
  fwd_adj_s1 <- fwd_adj1$coefficients[2,2]
  fwd_adj_s2 <- fwd_adj2$coefficients[2,2]
  fwd_adj_p1 <- fwd_adj1$coefficients[2,4]
  fwd_adj_p2 <- fwd_adj2$coefficients[2,4]
  
  p_val <- max(fwd_adj_p1, fwd_adj_p2)
  list(effect1 = fwd_adj_e1,
       effect2 = fwd_adj_e2,
       sd1 = fwd_adj_s1,
       sd1 = fwd_adj_s2,
       p1 = fwd_adj_p1,
       p2 = fwd_adj_p2,
       maxp = p_val)
}


##### Marginal Test for Pleiotropy with Bonferroni Correction ##########################################################################
marg_plei <- function(gt, y, cv){
  x1_marg <- cbind(y[,-2], cv)       #Discarding second trait
  x2_marg <- cbind(y[,-1], cv)       #Discarding first trait
  
  m_cv <- lm_binom(gt, cv)
  m1 <- lm_binom(gt, x1_marg)
  m2 <- lm_binom(gt, x2_marg)
  lrt_marg1 <- -2 * (m_cv$likelihood - m1$likelihood)
  lrt_marg2 <- -2 * (m_cv$likelihood - m2$likelihood)
  p_marg1 <- pchisq(lrt_marg1, df = 1, lower.tail = F)
  p_marg2 <- pchisq(lrt_marg2, df = 1, lower.tail = F)
  if(p_marg1 < 0.025 & p_marg2 < 0.025){D_marg <- "Significant"}else{D_marg <- "Insignificant"}
  return(list(marginal1 = lrt_marg1, marginal2 = lrt_marg2, p1 = p_marg1, p2 = p_marg2, est1 = m1$estimate, est2 = m2$estimate, margD = D_marg))
}

################################## Toy Data ########################################################################################################################
set.seed(123)
R1 <- 0^2
R2 <- 0^2
R <- 0.3^2
p1<-0.2              # MAF of Trait Locus #
p2<-0.2              # MAF of Marker Locus #
v1<-25               # Total variance of Trait 1 #
v2<-25               # Total variance of Trait 2 #
q1<-1-p1             # Major Allele Frequency of Trait Locus #
sm_size<-500         # Sample Size #

compute_stats <- function(a) {
  R_G_on_Y1_givenY2 <- abs(sqrt(R1))
  R_G_on_Y2_givenY1 <- abs(sqrt(R2))
  R_Y2_on_Y1_givenG <- abs(sqrt(R))
  varG <- 2*p1*q1
  varY1 <- v1
  varY2 <- v2
  # residual variances (so marginal var fixed)
  var_e1 <- varY1 - a[1]^2 * varG
  var_e2 <- varY2 - ( (a[3] * a[1] + a[2])^2 * varG + a[3]^2 * var_e1 )
  # covariances
  cov_Y1_T <- a[1] * varG
  cov_Y2_T <- (a[3] * a[1] + a[2]) * varG
  cov_Y1_Y2 <- a[3] * varY1 + a[2] * cov_Y1_T
  # Compute correlations needed for partial R^2 of G ~ Y1 | Y2 and G ~ Y2 | Y1
  cov_G_Y1 <- cov_Y1_T
  cov_G_Y2 <- cov_Y2_T
  # conditional cov/vars (conditioning on T)
  cov_Y1Y2_givenT <- cov_Y1_Y2 - (cov_Y1_T * cov_Y2_T) / varG
  var_Y1_givenT  <- varY1 - (cov_Y1_T^2) / varG
  var_Y2_givenT  <- varY2 - (cov_Y2_T^2) / varG
  # partial correlation r_{Y1,Y2 | T}
  r_partial <- cov_Y1Y2_givenT / sqrt(var_Y1_givenT * var_Y2_givenT)
  # partial R^2s: two directional partial R^2 measures
  # R2_{Y1 | Y2, T}  = (corr(Y1_resid, Y2_resid))^2  -- same as r_partial^2
  r1 <- r_partial^2
  # For R2_{Y2 | Y1, T} the same numeric value (partial correlation squared is symmetric).
  # But if you prefer directional partials for "unique contribution to G", compute increments in regressions for G:
  # correlations
  rho_Y1G <- cov_G_Y1 / sqrt(varY1 * varG)
  rho_Y2G <- cov_G_Y2 / sqrt(varY2 * varG)
  rho_Y1Y2 <- cov_Y1_Y2 / sqrt(varY1 * varY2)
  # partial correlation of G with Y1 given Y2:
  r_GY1_givenY2 <- (rho_Y1G - rho_Y1Y2 * rho_Y2G) / sqrt((1 - rho_Y2G^2)*(1 - rho_Y1Y2^2))
  r_GY2_givenY1 <- (rho_Y2G - rho_Y1Y2 * rho_Y1G) / sqrt((1 - rho_Y1G^2)*(1 - rho_Y1Y2^2))
  # partial R2s for directional interpretation
  c((R_G_on_Y1_givenY2 - r_GY1_givenY2), (R_G_on_Y2_givenY1 - r_GY2_givenY1), (R_Y2_on_Y1_givenG - r_partial))
}

xseq <- seq(-3, 3, length.out = 5)
yseq <- seq(-3, 3, length.out = 5)
zseq <- seq(-3, 3, length.out = 5)

# Evaluate residual norm
grid <- expand.grid(x1 = xseq, x2 = yseq, x3 = zseq)
grid$resid <- apply(grid, 1, function(v) sum(compute_stats(v)^2))

# Best initial guess
x0 <- as.numeric(grid[which.min(grid$resid), 1:3])
slv <- slfn(compute_stats, x0)
slv[which(slv < 1e-10)] <- 0

alpha1<-slv[1]             # Heritability of Trait 1 #
beta2<-slv[2]             # Heritability of Trait 2 horizontally #
beta1<-slv[3]            # Heritability for mediation to Trait 2 #
ld<-0.95*min(p1*(1-p2),(1-p1)*p2)    # LD b/w Trait and Marker Locus #
pvln1<-0.1           # Dichotomising Prevalance of Trait 1 #
pvln2<-0.2           # Dichotomising Prevalance of Trait 1 #
ve1<-v1-(alpha1^2)*((p1^2+q1^2)-(p1^2-q1^2)^2)
ve2<-v2-((beta1*alpha1+beta2)^2)*((p1^2+q1^2)-(p1^2-q1^2)^2)-(beta1^2)*ve1
a1<-c(-alpha1,0,alpha1)
b2<-c(-beta2,0,beta2)

#simulation#
A<-vector("numeric",sm_size)
x<-vector("numeric",sm_size)
y<-matrix(nrow = sm_size,ncol = 2)
y1_bin<-vector("numeric",sm_size)
y2_bin<-vector("numeric",sm_size)
M<-vector("numeric",sm_size)


for (i in 1:sm_size) {
  env<-vector("numeric")
  A[i]<-rbinom(1,2,p1)         #Simulation of Trait Locus#
  env[1]<-sqrt(ve1/2)*(rchisq(1,1))-sqrt(ve1/2)
  env[2]<-sqrt(ve2/2)*(rchisq(1,1))-sqrt(ve2/2)
  if(A[i]==0){
    y[i,1]<-a1[1]+env[1]                         #Simulation of Trait 1#
    y[i,2]<-beta1*y[i,1]+b2[1]+env[2]            #Simulation of Trait 2#
  }else if(A[i]==1){
    y[i,1]<-a1[2]+env[1]
    y[i,2]<-beta1*y[i,1]+b2[2]+env[2]
  }else{
    y[i,1]<-a1[3]+env[1]
    y[i,2]<-beta1*y[i,1]+b2[3]+env[2]
  }
  
  #Dichotomising The Traits#
  ds1<-function(C){
    (p1^2)*(1-pnorm((C-a1[3])/sqrt(ve1)))+(2*p1*(1-p1))*(1-pnorm((C-a1[2])/sqrt(ve1)))+((1-p1)^2)*(1-pnorm((C-a1[1])/sqrt(ve1)))-pvln1
  }
  c1<-slfn(ds1,0)
  if(y[i,1]>c1){y1_bin[i]<-1}else{y1_bin[i]<-0}
  ds2<-function(C){
    (p1^2)*(1-pnorm((C-beta1*a1[3]-b2[3])/sqrt(beta1^2*ve1+ve2)))+(2*p1*(1-p1))*(1-pnorm((C-beta1*a1[2]-b2[2])/sqrt(beta1^2*ve1+ve2)))+((1-p1)^2)*(1-pnorm((C-beta1*a1[1]-b2[1])/sqrt(beta1^2*ve1+ve2)))-pvln2
  }
  c2<-slfn(ds2,0)
  if(y[i,1]>c1){y1_bin[i]<-1}else{y1_bin[i]<-0}
  if(y[i,2]>c2){y2_bin[i]<-1}else{y2_bin[i]<-0}
  
  #Simulation of Marker Locus#
  if(A[i]==2){
    M_r<-rmultinom(1,1,c(((-ld+p1*(1-p2))/p1)^2,2*((ld+p1*p2)/p1)*((-ld+p1*(1-p2))/p1),((ld+p1*p2)/p1)^2))
    for(j in 1:3){if(M_r[j]==1){M[i]<-j-1}}
  }else if(A[i]==1){
    M_r<-rmultinom(1,1,c(((p1*(1-p2)-ld)*(ld+(1-p1)*(1-p2)))/(p1*(1-p1)),(((ld+p1*p2)*(ld+(1-p1)*(1-p2)))+(((1-p1)*p2-ld)*(p1*(1-p2)-ld)))/(p1*(1-p1)),((ld+p1*p2)*((1-p1)*p2-ld))/(p1*(1-p1))))
    for(j in 1:3){if(M_r[j]==1){M[i]<-j-1}}              
  }else{
    M_r<-rmultinom(1,1,c(((ld+(1-p1)*(1-p2))/(1-p1))^2,2*((ld+(1-p1)*(1-p2))/(1-p1))*(((1-p1)*p2-ld)/(1-p1)),(((1-p1)*p2-ld)/(1-p1))^2))
    for(j in 1:3){if(M_r[j]==1){M[i]<-j-1}} 
  }
}
Y<-cbind(A,M,y[,1],y[,2])

############################################### Implementation ###########################################################################################################################################################

hoprr(Y[,2], cbind(Y[,3], Y[,4]))
plafiu(Y[,2], cbind(Y[,3], Y[,4]))
