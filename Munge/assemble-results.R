#source("../Run/simulate-data.R")

#library(purrr)
#library(pseval)

calc_ve <- function(betaco, riskfun, Sobs) {
  
  dat1 <- data.frame(S.1 = Sobs, Z = 1)
  dat0 <- data.frame(S.1 = Sobs, Z = 0)
  
  beta <- unlist(betaco[1, -length(betaco)])
  
  R1 <- riskfun(dat1, beta)
  R0 <- riskfun(dat0, beta)
  
  log(log(1 - R0) / log(1 - R1))
  #1 - R1 / R0
  
  
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


calc_rj_boot <- function(psevalj, Sobs) {
  
  ves <- psevalj$bootstraps %>% 
    by_row(calc_ve, riskfun = psevalj$risk.function, Sobs = Sobs)
  
  ves$.out
  
}

Deetz.tilde <- function(psbyj, psleftj, Sobs) {
  
  S.topred <- psbyj$augdata$S.1
  S.topred <- sort(S.topred[!is.na(S.topred)])
  
  byboot <- calc_rj_boot(psbyj, S.topred)
  leftboot <- calc_rj_boot(psleftj, S.topred)
  
  #prod <- cross_d(list(by = byboot, left = leftboot))
  #Dij <- abs(prod$by - prod$left)
  
  Dij <- abs(unlist(byboot) - unlist(leftboot))
  
  Dij #sample(Dij, 5000) 
  
}

calc_D.hat.bias <- function(bytrial, leftout, scenario) {
  
  thisS <- map(bytrial$data, function(d) d$S.obs[!is.na(d$S.obs)])
  truecep <- true_cep(scenario, Sobs = thisS)
  
  estcep <- map2(leftout$pseval, thisS, function(d, x){
    
    data.frame(S.1 = x, estcep = rowMeans(do.call(cbind,calc_rj_boot(d, x))))
                 #calc_ve(as.data.frame(t(c(d$estimates$par, 0))), d$risk.function, Sobs = x))
  })
  
  byestcep <- map2(bytrial$pseval, thisS, function(d, x){
    
    data.frame(S.1 = x, estcep = rowMeans(do.call(cbind, calc_rj_boot(d, x))))
    #calc_ve(as.data.frame(t(c(d$estimates$par, 0))), d$risk.function, Sobs = x))
  })
  
  
  Deesj.hats <- unlist(map2(truecep$data, estcep, .f = function(x, y){ abs(x$VE - y$estcep)}))
  Deesj.tildes <- unlist(map2(byestcep, estcep, .f = function(x, y){ abs(x$estcep - y$estcep)}))
  mean(Deesj.hats - Deesj.tildes)
  
}


Deetz.diff <- function(psbyj, psleftj, Sobs) {
  
  S.topred <- psbyj$augdata$S.1
  S.topred <- sort(S.topred[!is.na(S.topred)])
  
  byboot <- calc_rj_boot(psbyj, S.topred)
  leftboot <- calc_rj_boot(psleftj, S.topred)
  
  #prod <- cross_d(list(by = byboot, left = leftboot))
  #Dij <- abs(prod$by - prod$left)
  
  Dij <- abs(unlist(byboot) - unlist(leftboot))
  
  ## calc D0j by marginalizing over the surrogate
  
  leftboot0 <- lapply(leftboot, function(llb) rep(mean(llb), length(llb)))
  #prod0 <- cross_d(list(by = byboot, left0 = leftboot0))
  
  Dij0 <- abs(unlist(byboot) - unlist(leftboot0))
  
  Dij - Dij0
  
}

CEP.ci <- function(psbyj, psleftj, psoverall, Sobs, D.tilde) {
  
  S.topred <- psbyj$augdata$S.1
  S.topred <- sort(S.topred[!is.na(S.topred)])
  
  ve.J.o <- calc_risk(psleftj, newdata = S.topred,
                      contrast = function(R0, R1) log(log(1 - R0) / log(1 - R1)))[, c("S.1", "Y", "Y.lower.CL.0.95", "Y.upper.CL.0.95")]
  ve.J.o$Y.lower.pred <- ve.J.o$Y.lower.CL.0.95 - D.tilde
  ve.J.o$Y.upper.pred <- ve.J.o$Y.upper.CL.0.95 + D.tilde
  
  ve.J.b <- calc_risk(psbyj, newdata = S.topred,
                      contrast = function(R0, R1) log(log(1 - R0) / log(1 - R1)))[, c("S.1", "Y", "Y.lower.CL.0.95", "Y.upper.CL.0.95")]
  #ve.J.b$Y.lower.pred <- ve.J.b$Y.lower.CL.0.95 #- D.hat
  #ve.J.b$Y.upper.pred <- ve.J.b$Y.upper.CL.0.95 #+ D.hat
  
  ## model based confidence band
  
  parboots <- psbyj$bootstraps[, 3:4]
  cv <- cov(parboots)
  v.z <- cv[1, 1]
  v.sz <- cv[2, 2]
  c.zsz <- cv[1, 2]
  
  se.line <- sqrt(v.z + S.topred^2 * v.sz + S.topred * c.zsz)
  zed <- abs(qnorm(.05 / length(S.topred)))
  ve.J.b$Y.lower.model <- ve.J.b$Y - zed * se.line
  ve.J.b$Y.upper.model <- ve.J.b$Y + zed * se.line
  
  parbootso <- psleftj$bootstraps[, 3:4]
  cvo <- cov(parbootso)
  v.zo <- cv[1, 1]
  v.szo <- cv[2, 2]
  c.zszo <- cv[1, 2]
  
  seo.line <- sqrt(v.zo + S.topred^2 * v.szo + S.topred * c.zszo)
  
  ve.J.o$Y.lower.pred <- ve.J.o$Y - zed * seo.line - D.tilde
  ve.J.o$Y.upper.pred <- ve.J.o$Y + zed * seo.line + D.tilde
  
  data.frame(ve.by = ve.J.b, ve.over = ve.J.o)
  
}


CEP.coverage <- function(scenario, bytrial, leftout, psoverall, Sobs) {
  
  ## calculate D.tilde
  
  Djays <- unlist(map2(bytrial$pseval, leftout$pseval, .f = Deetz.tilde, Sobs))
  D.tilde.95 <- quantile(Djays, 0.975)
  D.mean <- mean(Djays)
  
  Sobsin <- lapply(bytrial$pseval, function(pf){ 
    tmp <- pf$augdata$S.1
    sort(tmp[!is.na(tmp)])
  })
  
  CEP.ests.95 <- map2(bytrial$pseval, leftout$pseval, .f = CEP.ci, psoverall, Sobs, D.tilde.95)
  CEP.ests.mean <- map2(bytrial$pseval, leftout$pseval, .f = CEP.ci, psoverall, Sobs, D.mean)
  CEP.true <- true_cep(scenario = scenario, Sobs = Sobsin)
  CEP.true$intervals.95 <- CEP.ests.95
  CEP.true$intervals.mean <- CEP.ests.mean
  
  
  ## check coverage
  
  cover.by <- map2_lgl(CEP.true$data, CEP.true$intervals.95, .f = function(data, intervals) {
    all(data$VE < intervals$ve.by.Y.upper.pred & 
          data$VE > intervals$ve.by.Y.lower.pred) 
    })
  
  cover.over <- map2_lgl(CEP.true$data, CEP.true$intervals.95, .f = function(data, intervals) {
    all(data$VE < intervals$ve.over.Y.upper.pred & 
          data$VE > intervals$ve.over.Y.lower.pred) 
  })
  
  
  cover.by.mean <- map2_lgl(CEP.true$data, CEP.true$intervals.mean, .f = function(data, intervals) {
    all(data$VE < intervals$ve.by.Y.upper.model & 
          data$VE > intervals$ve.by.Y.lower.model) 
  })
  
  cover.over.mean <- map2_lgl(CEP.true$data, CEP.true$intervals.mean, .f = function(data, intervals) {
    all(data$VE < intervals$ve.over.Y.upper.pred & 
          data$VE > intervals$ve.over.Y.lower.pred) 
  })
  
  
  d.hat.bias <- calc_D.hat.bias(bytrial, leftout, scenario)
  
  c(cover.by.trial = mean(cover.by), cover.overall = all(cover.over), 
    cover.by.trial.mean = mean(cover.by.mean), cover.overall.mean = all(cover.over.mean), d.hat.bias)
  
  
}


calc_DjD0_pee <- function(bytrial, leftout) {
  
  stopifnot(nrow(bytrial) == nrow(leftout))
  
  thisS <- unlist(map(bytrial$data, function(d) d$S.obs[!is.na(d$S.obs)]))
  thisS <- sort(sample(thisS, 300))[-c(1:50, 250:300)]
  Deesj.boots <- map2(bytrial$pseval, leftout$pseval, .f = Deetz.diff, Sobs = thisS)
  
  Djs.boot <- do.call(cbind, Deesj.boots)
  mean(rowMeans(Djs.boot) > 0) # p.value for PIGS value
  
}


test_one <- function(psby) {
  
  ests <- psby$estimates$par
  ses <- sqrt(diag(cov(psby$bootstraps)))
  ses <- ses[-length(ses)]
  
  p.wem <- 2 * pnorm(-abs(ests[4] / ses[4])) ## 2 sided test for WEM
  
  p.issacs <- pnorm(ests[4] / ses[4], lower.tail = TRUE)  ## 1 sided test for alpha4 < 0
  p.harm <- pnorm(ests[3] / ses[3], lower.tail = TRUE)
  harm.est <- ests[3]
  harm.width <- ses[3] * 2 * 1.96
      ## test that alpha1 > 0, i.e., treatment is harmful at S(1) = 0, should *fail to reject*
  
  data.frame(p.wem = p.wem, p.issacs = p.issacs, harm.est = harm.est, harm.width = harm.width, p.harm = p.harm)
  
}


run_tests <- function(bytrial) {
  
  testsbytrial <- map_df(bytrial$pseval, .f = test_one)
  ## check for harm first
  
  rej.issacs <- all(testsbytrial$p.issacs < .05)
  rej.harm <- all(testsbytrial$p.noharm < 0.05)
  
  c(rej.0wem.dir = all(testsbytrial$p.wem  < 0.05), 
    rej.1aascs = rej.issacs,
    harm.mean = mean(testsbytrial$harm.est), 
    harm.width.mean = mean(testsbytrial$harm.width), 
    test.harm = rej.harm)
  
  
}


