library(pseval)
library(dplyr)
library(purrr)

expit <- function(x){
  exp(x)/(1 + exp(x))
}

generate_data <- function(n, beta = 0, gamma, alpha = c(1, 1, 1, 1), sig.s = .1, mu.s1 = 1){

  Z <- rbinom(n, 1, .5)
  X <- rnorm(n, sd = sqrt(0.5))  ## pre-treatment covariate

  ## generate S conditional on Z

  S.0 <- 0 # rnorm(n, sd = sig.s)
  S.1 <- mu.s1 + gamma * X + rnorm(n, sd = sig.s)

  #risk.obs <- (-1 - 0.0 * S.1 - 0 * Z - beta * (S.1 - S.0) * Z)
  risk.0 <- alpha[1] + alpha[3] * S.1 + beta * X
  risk.1 <- alpha[1] + alpha[2] + (alpha[3] + alpha[4]) * S.1 + beta * X

  Y.0 <- rpois(n, exp(risk.0))
  Y.1 <- rpois(n, exp(risk.1))
  Y.obs <- ifelse(Z == 1, Y.1, Y.0)

  S.obs <- ifelse(Z == 1, S.1, S.0)

  ## CPV measure noisy S.1 at the end of the study for placebo subjects non-event

  CPV <- S.1 + rnorm(n, sd = .1)
  CPV[Z == 1 | Y.obs == 1] <- NA

  ## BSM measure noisy S.0 at start

  BSM <- S.0 + rnorm(n, sd = .1)
  S.1[Z == 0] <- NA
  S.0[Z == 1] <- NA
  S.obs <- ifelse(Z == 1, S.1, S.0)

  data.frame(Z, BIP = X, CPV, BSM, S.obs = S.obs, Y.obs)

}

generate_cep <- function(alpha, S.1) {
  
  
  lam.0 <- exp(alpha[1] + alpha[3] * S.1)
  lam.1 <- exp(alpha[1] + alpha[2] + (alpha[3] + alpha[4]) * S.1)
  
  risk.0 <- ppois(0, lambda = lam.0, lower.tail = FALSE)
  risk.1 <- ppois(0, lambda = lam.1, lower.tail = FALSE)
  
  #1 - risk.1 / risk.0
  #log(log(1 - risk.0) / log(1 - risk.1))
  ## should be the same thing
  log(lam.0 / lam.1)
  
  
}

#set.seed(1123)
#sample(1:12, 12)
dex <- c( 10, 12,  9,  5,  7,  2,  8,  4, 11,  6,  3,  1)

get_sizes <- function(){

  as.integer(seq.int(100, 1000, length.out = 12))[dex]

}



## constants

## parameters
# set.seed(1430)
## 

# mu.list <- runif(12, -2, 2)
mu.list <- c(0.055, -0.046, -1.692, -1.957, -1.575, -1.959,
             -1.351, -1.925, 1.192, -1.549, 1.024, -0.112)

# a1.list <- runif(12, -.2, .75)
#a1.list <- c(-0.532, -0.763, 1.508, 1.055, 0.421, 0.559,
#             -1.704, 0.331, -1.604, -1.955, -1.104, -1.744)

a1.list <- c(0.149, 0.094, 0.633, 0.526, 0.375, 0.408, 
             -0.13, 0.354, -0.106, -0.189, 0.013, -0.139)

# a2.list <- runif(12, -.35, -0.05)
#a2.list <- c(.048, 1.388, 1.807, 1.154, 1.084, 1.758,
#             1.035, 1.705, 1.894, 1.246, 1.489, 1.915)

a2.list <- c(-0.336, -0.234, -0.108, -0.304, -0.325, -0.123, 
             -0.339, -0.139, -0.082, -0.276, -0.203, -0.076)

a1 <- 0
a4 <- -0.35
a3 <- -0.25
g1 <- 0.95

parameter_list <- function(scenario) {
  
  alpha.list <- vector(mode = "list", length = 12)

  if(scenario == 1) {
    
    alpha.list[1:6] <- list(c(a1, 0, a3, a4))
    alpha.list[7:12] <- list(c(a1, 0, a3, 0))
    
    gamma.list <- rnorm(12, mean = g1, sd = .25)
    sig.list <- runif(12, sqrt(0.25), sqrt(0.5))
    mu.list <- rep(1, 12)
    
  } else if(scenario == 2) {
    
    alpha.list <- vector(mode = "list", length = 12)
    alpha.list[1:6] <- list(c(a1, 0, -a3, -a4))
    alpha.list[7:12] <- list(c(a1, 0, a3, a4))
    
    gamma.list <- rnorm(12, mean = g1, sd = .25)
    sig.list <- runif(12, sqrt(0.25), sqrt(0.5))
    mu.list <- rep(1, 12)
    
  } else if(scenario == 3) { 
    
    alpha.list <- lapply(1:12, function(i) c(a1, a1.list[i], a3, a4))
    
    gamma.list <- rnorm(12, mean = g1, sd = .25)
    sig.list <- runif(12, sqrt(0.25), sqrt(0.5))
    mu.list <- rep(1, 12)
    
  } else if(scenario == 4) {
    
    alpha.list <- lapply(1:12, function(i) c(a1, a2.list[i], a3, a4))
    
    gamma.list <- rnorm(12, mean = g1, sd = .25)
    sig.list <- runif(12, sqrt(0.25), sqrt(0.5))
    mu.list <- rep(1, 12)
    
  } else if(scenario == 5) {
    
    alpha.list <- lapply(1:12, function(i) c(a1, 0, a3 + .25, a4 - 0.5))
    
    gamma.list <- rnorm(12, mean = g1, sd = .1)
    sig.list <- runif(12, sqrt(0.15), sqrt(0.35))
    mu.list <- rep(1, 12)
    
  } else if(scenario == 6) {
    
    alpha.list <- lapply(1:12, function(i) c(a1, 0, a3, a4))
    
    gamma.list <- rnorm(12, mean = g1, sd = .25)
    sig.list <- runif(12, sqrt(0.25), sqrt(0.5))
    mu.list <- c(0.055, -0.046, -1.692, -1.957, -1.575, -1.959,
                 -1.351, -1.925, 1.192, -1.549, 1.024, -0.112)
    
  }
  
  return(list(alpha.list = alpha.list, gamma.list = gamma.list, 
              sig.list = sig.list, mu.list = mu.list))
  
}


simulate_scen <- function(scenario, n.trials = 12, samp.sizes = get_sizes()){


  outdata <- vector(mode = "list", length = length(samp.sizes))
  i <- 1

  pars <- parameter_list(scenario)

  for(j in samp.sizes){

    outdata[[i]] <- cbind(trial = i,
                          generate_data(j, gamma = pars$gamma.list[i],
                                        alpha = pars$alpha.list[[i]],
                                        sig.s = pars$sig.list[i], mu.s1 = pars$mu.list[i]))
    

    i <- i + 1
  }

  data_frame(trial = 1:n.trials, data = outdata)

}

true_cep <- function(scenario, n.trials = 12, samp.sizes = get_sizes(), Sobs) {
  
  outdata <- vector(mode = "list", length = length(samp.sizes))
  i <- 1
  
  pars <- parameter_list(scenario)
  
  for(i in 1:length(samp.sizes)){
    
    outdata[[i]] <- data.frame(S.1 = Sobs[[i]], 
                          VE = generate_cep(pars$alpha.list[[i]], Sobs[[i]]))
    
  }
  
  data_frame(trial = 1:n.trials, data = outdata, pars = pars$alpha.list)
  
}


risk_poisson2 <- function(model = Y ~ S.1 * Z, D = 5000 ){
  
  arglist <- as.list(match.call())
  rval <- function(psdesign){
    
    expanded <- pseval:::expand_augdata(model, psdesign, D = D)
    
    trtmat <- expanded$noimp
    Y.trt <- expanded$noimp.Y
    
    untrt.expand <- expanded$imp
    Y.untrt <- expanded$imp.Y
    
    attr(trtmat, "hasoffset") <- attr(expanded, "hasoffset")
    
    likelihood <- function(beta){
      
      
      if(attr(trtmat, "hasoffset")){
        beta.o <- c(beta, 1)
      } else{
        beta.o <- beta
      }
      
      lambda <- exp(trtmat %*% beta.o)
      
      trtlike <- dpois(Y.trt, lambda, log = TRUE)
      
      if(!is.null(untrt.expand) & !is.null(Y.untrt)){
        
        lambda.untrt <- exp(untrt.expand %*% beta.o)
        untrted <- matrix(dpois(Y.untrt, lambda.untrt, log = TRUE), nrow = D, byrow = TRUE)
      } else untrted <- matrix(1)
      
      -1 * (sum(trtlike) + sum(colMeans(untrted)))
      
    }
    
    psdesign$risk.function <- function(data, beta, t = 0){ # P(Y <= t | S, Z)
      
      lambda <- as.vector(exp(model.matrix(model[-2], data) %*% beta))
      
      ppois(t, lambda, lower.tail = FALSE)
      
    }
    
    psdesign$likelihood <- likelihood
    psdesign$risk.model <- list(model = "poisson", args = arglist )
    psdesign$nparam <- ncol(trtmat) + ifelse(attr(trtmat, "hasoffset"), -1, 0)
    psdesign$param.names <- colnames(model.matrix(model, psdesign$augdata))
    
    psdesign
    
  }
  ## return a likelihood closure
  
  class(rval) <- c("ps", "riskmodel")
  rval
  
}



analyze_bytrial <- function(trials, n.samps = 200){

  trials2 <- trials %>% mutate(

    pseval = map(data, ~ {
      psdesign(., Z = Z, Y = Y.obs, S = S.obs, BIP = BIP) +
        integrate_semiparametric(S.1 ~ BIP, S.1 ~ 1) +
        risk_poisson2(D = 50) + ps_estimate() + ps_bootstrap(n.boots = 100, progress.bar = FALSE)

    }))

  trials2 
}


analyze_overall <- function(trials, n.samps = 200){

  fit <- psdesign(do.call(rbind, trials$data), Z = Z, Y = Y.obs, S = S.obs, BIP = BIP) +
    integrate_semiparametric(S.1 ~ BIP, S.1 ~ 1) +
    risk_poisson2(D = 50) + ps_estimate() + ps_bootstrap(n.boots = 100, progress.bar = FALSE)

  list(pseval = fit)

}


predict_leaveoneout <- function(trials, n.samps = 200){

  trials$leftout <- lapply(1:nrow(trials), function(i){

    do.call(rbind, trials$data[-i])

  })

  trials2 <- trials %>% mutate(

    pseval = map(leftout, ~ {
      psdesign(., Z = Z, Y = Y.obs, S = S.obs, BIP = BIP) +
        integrate_semiparametric(S.1 ~ BIP, S.1 ~ 1) +
        risk_poisson2(D = 50) + ps_estimate() + ps_bootstrap(n.boots = 100, progress.bar = FALSE)

    }))


  trials2

}


simulate_one <- function(scenario, output, seed){

  set.seed(seed)

  data <- simulate_scen(scenario)
  bytrial <- analyze_bytrial(data)
  overall <- analyze_overall(data)
  leftout <- predict_leaveoneout(data)

  save(bytrial, file = gsub("Rfile", "Rfile-bytrial", output))
  save(overall, file = gsub("Rfile", "Rfile-overall", output))
  save(leftout, file = gsub("Rfile", "Rfile-leftout", output))

}

moar_bootstraps <- function(scenario, output, seed){

  load(file = gsub("Rboot", "Rfile-bytrial", output))
  load(file = gsub("Rboot", "Rfile-overall", output))
  load(file = gsub("Rboot", "Rfile-leftout", output))

  set.seed(seed)

  for(i in 1:nrow(bytrial)){

    tmpby <- bytrial$pseval[[i]]
    tmpby <- tmpby + ps_bootstrap(n.boots = 100, progress.bar = FALSE)
    bytrial$pseval[[i]]$bootstraps <- rbind(bytrial$pseval[[i]]$bootstraps,
                                            tmpby$bootstraps)

    tmpleft <- leftout$pseval[[i]]
    tmpleft <- tmpleft + ps_bootstrap(n.boots = 100, progress.bar = FALSE)
    leftout$pseval[[i]]$bootstraps <- rbind(leftout$pseval[[i]]$bootstraps,
                                            tmpleft$bootstraps)

  }

  tmpover <- overall$pseval
  tmpover <- tmpover + ps_bootstrap(n.boots = 100, progress.bar = FALSE)
  overall$pseval$bootstraps <- rbind(overall$pseval$bootstraps,
                                     tmpover$bootstraps)


  save(bytrial, file = gsub("Rboot", "Rfile-bytrial", output))
  save(overall, file = gsub("Rboot", "Rfile-overall", output))
  save(leftout, file = gsub("Rboot", "Rfile-leftout", output))

}



