# Reparameterized Theta


dir <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(dir,'/Final_Pipeline.R'))


#Helper functions skip to line 354

#Pairwise interaction I_ij
get_I_ij <- function(data, prop, theta=1){
  if (is.numeric(prop)) {
    prop <- data[,prop]
  }
  
  nSpecies = length(prop)
  dat <- data %>% select(all_of(prop))
  
  tibble(do.call(cbind, combn(colnames(dat), 2, FUN= function(x) 
    list(setNames(data.frame((2*((nSpecies^2)*dat[,x[1]]*dat[,x[2]])^theta)/(nSpecies*(nSpecies - 1))), 
                  paste(x, collapse="_")) ))))
}

#Average interaction I_ij
get_I_ij_E_AV <- function(data, prop, what, theta = 1){
  if (!(what=='E' | what=='AV')){
    stop('Please specify whether you want E or AV term')
  }
  
  AV <- rowSums(get_I_ij(data, prop, theta))
  
  if (what=='AV'){
    return(tibble('AV'=AV))
  }
  else {
    nSpecies = length(prop)
    normE <- 2 * nSpecies/(nSpecies - 1)
    return(tibble('E'=normE*AV))
  }
}

#Additive species interaction I_ij
get_I_ij_ADD <- function(data, prop, theta=1){
  nSpecies = length(prop)
  # Get normal _add terms assuming theta =1 first and then scale them 
  # with scaling factor to get I_ij terms
  
  if(any(grepl('_add', colnames(data)))){
    data <- data[, !grepl('_add', colnames(data))]
  }
  add_terms <- (DI_data_ADD(prop = prop, data = data, theta=theta))$ADD
  
  #Scaling the terms to get the I_ij variants
  add_terms <- (2*add_terms*(nSpecies^2)^theta)/(nSpecies*(nSpecies - 1))
  return(add_terms)
}

#Functional group interaction I_ij
get_I_ij_FG <- function(data, prop, FG, theta=1){
  nSpecies = length(prop)
  
  FG_terms <- (DI_data_FG(prop = prop, data = data, FG=FG, theta=theta))$FG
  
  #For maintaining consistency with DImodels package
  colnames(FG_terms) <- paste0('FG_', colnames(FG_terms))
  #Scaling the terms to get the I_ij variants
  FG_terms <- (2*FG_terms*(nSpecies^2)^theta)/(nSpecies*(nSpecies - 1))
  return(FG_terms)
}

#Wrapper function to produce all types of I_ij terms
DI_data_I_ij <- function(data, prop, 
                        what=c('AV', 'E', 'FG', 'ADD', 'FULL'), 
                        theta=1, FG=NULL){
  
  if (('FG' %in% what) & is.null(FG)){
    stop('Please specify functional groupings for FG model')
  }
  
  interactions <- list()
  
  if ("E" %in% what) {
    interactions[[1]] <- get_I_ij_E_AV(data, prop, theta, what='E')
  }
  if ("AV" %in% what) {
    interactions[[2]] <- get_I_ij_E_AV(data, prop, theta, what='AV')
  }
  if ("FG" %in% what) {
    interactions[[3]] <- as_tibble(get_I_ij_FG(data, prop, FG, theta))
  }
  if ("ADD" %in% what) {
    interactions[[4]] <- as_tibble(get_I_ij_ADD(data, prop, theta))
  }
  if ("FULL" %in% what) {
    interactions[[5]] <- get_I_ij(prop = prop, data = data, theta = theta)
  } 
  interactions <- bind_cols(interactions)
  common_cols <- intersect(colnames(data), colnames(interactions))
  
  if(length(common_cols)>0){
    warnings('Columns with names same as interaction terms were present,
             these columns will be updated')
    data %>% select(-all_of(common_cols))
  }
  
  result <- bind_cols(data,interactions)
  return(result)
}



# Main function to estimate theta, added one extra paramter which asks for upper
# boundary for searching theta. If not specified, defaults to 1.5
DI_theta <- function (obj, DImodel, FGnames, prop, nSpecies, family, lower_boundary=0.01, upper_boundary=1.5) 
{
  if (missing(FGnames)) {
    FGnames <- NULL
  }

  mm <- model.matrix(obj)
  
  int_terms <- switch(EXPR = DImodel, AV = grep("AV", 
                                                colnames(mm)),
                      E = grep("E", colnames(mm)), 
                      FG = grep("FG", colnames(mm)),
                      ADD = grep("_add", colnames(mm)), 
                      FULL = get_int_terms_FULL(mm_names = colnames(mm), prop_names = prop))
  upper_boundary <- upper_boundary
  theta_info <- get_theta_info(prop=prop, lower_boundary = lower_boundary, upper_boundary = upper_boundary, 
                               DImodel = DImodel, obj = obj, family = family, int_terms = int_terms, 
                               nSpecies = nSpecies, FGnames = FGnames)
  # Theta estimation is done and finished here
  theta_hat <- theta_info$theta_hat
  profile_loglik <- theta_info$profile_loglik
  #plot(profile_loglik$grid, profile_loglik$prof)
  
  if ((upper_boundary - theta_hat) < 0.01) {
    warning("Theta has reached the upper boundary, this may indicate lack of convergence and/or a problem with non-significant diversity effect.")
  }
  
  # This part of code is there just to maintain similarity with DImodels and 
  # incase we wish to compare the models generated with the new higher theta.
  if (DImodel %in% c("E", "AV")) {
    data_theta_E_AV <- obj$model
    
    data_theta_E_AV[,'AV'] <- get_I_ij_E_AV(data_theta_E_AV, prop, theta = theta_hat, 'AV')
    data_theta_E_AV[,'E'] <- get_I_ij_E_AV(data_theta_E_AV, prop, theta = theta_hat, 'E')
    
    mod_theta <- glm(formula(obj), family = family, data = data_theta_E_AV)
  }
  else if (DImodel == "FG") {
    data_theta_FG <- obj$model
    data_theta <- data.frame(mm)
    #data_theta[1:nSpecies] <- data_theta[1:nSpecies]^theta_hat
    new_FG <- get_I_ij_FG(prop = prop, FG = FGnames, 
                          data = data_theta, theta = theta_hat)
    FG_ <- new_FG
    FG_cols_in_the_data <- grep("FG_", colnames(data_theta_FG))
    if (length(FG_cols_in_the_data) > 0) {
      j <- 1
      for (i in FG_cols_in_the_data) {
        data_theta_FG[, i] <- FG_[, j]
        j <- j + 1
      }
    }
    old_formula <- formula(obj)
    new_formula <- paste0(old_formula[2], " ~ ", old_formula[3])
    mod_theta <- glm(as.formula(new_formula), family = family, 
                     data = data_theta_FG)
  }
  else if (DImodel == "ADD") {
    data_theta <- data.frame(mm, check.names = FALSE)
    new_ADD <- get_I_ij_ADD(prop = prop, data = data_theta, 
                            theta = theta_hat)
    ADD_cols_in_the_data <- int_terms
    if (length(ADD_cols_in_the_data) > 0) {
      j <- 1
      for (i in ADD_cols_in_the_data) {
        data_theta[, i] <- new_ADD[, j]
        j <- j + 1
      }
    }
    old_formula <- formula(obj)
    new_formula <- paste0(old_formula[2], " ~ ", old_formula[3])
    data_theta_ADD <- data_theta
    data_theta_ADD$y <- obj$model[,1]
    names(data_theta_ADD)[length(names(data_theta_ADD))] <- paste(old_formula[2])
    mod_theta <- glm(as.formula(new_formula), family = family, 
                     data = data_theta_ADD)
  }
  else {
    mm[, int_terms] <- as.matrix(get_I_ij(as.data.frame(mm), prop = prop, theta=theta_hat))
    colnames(mm) <- gsub("`", "", colnames(mm))
    names_mm <- paste0("`", colnames(mm)[int_terms], 
                       "`")
    resp_name <- paste(formula(obj))[2]
    ndata <- data.frame(obj$model[, resp_name], mm, check.names = FALSE)
    names(ndata)[1] <- resp_name
    new_formula_theta <- as.formula(paste(resp_name, "~", 
                                          "0+", paste(colnames(mm)[-int_terms], collapse = "+"), 
                                          "+", paste(names_mm, collapse = "+")))
    mod_theta <- glm(formula = new_formula_theta, family = family, 
                     data = ndata)
  }
  mod_theta$coefficients <- c(mod_theta$coefficients, theta = theta_hat)
  mod_theta$df.residual <- mod_theta$df.residual - 1
  mod_theta$profile_loglik <- profile_loglik
  mod_theta$aic <- AIC2(mod_theta)
  # return(list('model'=mod_theta,
  #             'Likelihood' = profile_loglik))
  return(mod_theta)
}



# This is the function that does the profile likelihood 
# Copied as is from DImodels
get_theta_info <- function (prop, lower_boundary, upper_boundary, DImodel, obj, family, int_terms, nSpecies, 
                            FGnames) 
{
  theta_hat <- optimize(proflik_theta, interval = c(lower_boundary, upper_boundary), 
                        maximum = TRUE, DImodel = DImodel, prop=prop, obj = obj, family = family, 
                        int_terms = int_terms, nSpecies = nSpecies, FGnames = FGnames)$maximum
  theta_grid <- seq(lower_boundary, upper_boundary + 1, length = 100)
  proflik_theta_vec <- Vectorize(proflik_theta, "theta")
  profile_loglik <- proflik_theta_vec(theta = theta_grid,prop=prop, obj = obj, 
                                      family = family, int_terms = int_terms, DImodel = DImodel, 
                                      nSpecies = nSpecies, FGnames = FGnames)
  return(list(theta_hat = theta_hat, profile_loglik = data.frame(prof = profile_loglik, 
                                                                 grid = theta_grid)))
}

# The actual function that is optimized to get theta value
# Copied as is from DImodels
proflik_theta <- function (theta, obj, prop, family, int_terms, DImodel, nSpecies, FGnames) 
{
  mm <- model.matrix(obj)
  if (DImodel %in% c("E", "AV")) {
    data_theta_E_AV <- obj$model
    #data_theta <- data.frame(mm) %>% select(species)
    #data_theta[1:nSpecies] <- data_theta[1:nSpecies]^theta
    
    data_theta_E_AV[,'AV'] <- get_I_ij_E_AV(data_theta_E_AV, prop, theta = theta, 'AV')
    data_theta_E_AV[,'E'] <- get_I_ij_E_AV(data_theta_E_AV, prop, theta = theta, 'E')
    
    
    fitted_model_theta <- glm(formula(obj), family = family, 
                              data = data_theta_E_AV)
    n <- nrow(fitted_model_theta$data)
    p <- length(fitted_model_theta$coef) + 1
    mu_hat <- fitted(fitted_model_theta)
    sigma_hat <- sqrt(sum(fitted_model_theta$residuals^2)/(fitted_model_theta$df.residual) * 
                        (n - p)/n)
    llik <- sum(dnorm(fitted_model_theta$y, mu_hat, sigma_hat, 
                      log = TRUE))
    #print(paste('Theta:', theta,'LogLik:',llik))
  }
  else if (DImodel == "FG") {
    data_theta_FG <- obj$model
    data_theta <- data.frame(mm)
    #data_theta[1:nSpecies] <- data_theta[1:nSpecies]^theta
    new_FG <- get_I_ij_FG(prop = prop, FG = FGnames, 
                          data = data_theta, theta=theta)
    FG_ <- new_FG
    FG_cols_in_the_data <- grep("FG_", colnames(data_theta_FG))
    if (length(FG_cols_in_the_data) > 0) {
      if (length(FG_cols_in_the_data) != ncol(FG_)) {
        stop("please rename variables beginning with 'FG_'")
      }
      j <- 1
      for (i in FG_cols_in_the_data) {
        data_theta_FG[, i] <- FG_[, j]
        j <- j + 1
      }
    }
    old_formula <- formula(obj)
    new_formula <- paste0(old_formula[2], " ~ ", old_formula[3])
    fitted_model_theta <- glm(as.formula(new_formula), family = family, 
                              data = data_theta_FG)
    n <- nrow(fitted_model_theta$data)
    p <- length(fitted_model_theta$coef) + 1
    mu_hat <- fitted(fitted_model_theta)
    sigma_hat <- sqrt(sum(fitted_model_theta$residuals^2)/(fitted_model_theta$df.residual) * 
                        (n - p)/n)
    llik <- sum(dnorm(fitted_model_theta$y, mu_hat, sigma_hat, 
                      log = TRUE))
  }
  else if (DImodel == "ADD") {
    data_theta <- data.frame(mm, check.names = FALSE)
    new_ADD <- get_I_ij_ADD(prop = prop, data = data_theta, 
                                 theta = theta)
    ADD_cols_in_the_data <- int_terms
    if (length(ADD_cols_in_the_data) > 0) {
      j <- 1
      for (i in ADD_cols_in_the_data) {
        data_theta[, i] <- new_ADD[, j]
        j <- j + 1
      }
    }
    old_formula <- formula(obj)
    new_formula <- paste0(old_formula[2], " ~ ", old_formula[3])
    data_theta_ADD <- data_theta
    data_theta_ADD$y <- obj$model[,1]
    names(data_theta_ADD)[length(names(data_theta_ADD))] <- paste(old_formula[2])
    fitted_model_theta <- glm(as.formula(new_formula), family = family, 
                              data = data_theta_ADD)
    n <- nrow(fitted_model_theta$data)
    p <- length(fitted_model_theta$coef) + 1
    mu_hat <- fitted(fitted_model_theta)
    sigma_hat <- sqrt(sum(fitted_model_theta$residuals^2)/(fitted_model_theta$df.residual) * 
                        (n - p)/n)
    llik <- sum(dnorm(fitted_model_theta$y, mu_hat, sigma_hat, 
                      log = TRUE))
  }
  else {
    mm[, int_terms] <- as.matrix(get_I_ij(as.data.frame(mm), prop = prop, theta=theta))
    ndata <- obj$model
    ndata$mm <- mm
    fitted_model_theta <- glm(update.formula(formula(obj), 
                                             . ~ 0 + mm), family = family, data = ndata)
    n <- nrow(fitted_model_theta$data)
    p <- length(fitted_model_theta$coef) + 1
    mu_hat <- fitted(fitted_model_theta)
    sigma_hat <- sqrt(sum(fitted_model_theta$residuals^2)/(fitted_model_theta$df.residual) * 
                        (n - p)/n)
    llik <- sum(dnorm(fitted_model_theta$y, mu_hat, sigma_hat, 
                      log = TRUE))
  }
  #print(paste('Theta:', theta,'Likelihood:',llik))
  return(llik)
}


get_int_terms_FULL <- function (mm_names, prop_names) 
{
  mm_names <- gsub("`", "", mm_names)
  all_pairwise_ints <- which(lengths(regmatches(mm_names, gregexpr("_", 
                                                                   mm_names))) == 1)
  prop_match <- lapply(strsplit(mm_names, "_"), function(x) {
    #print(x)
    x %in% prop_names
  })
  prop_pairwise <- which(unlist(lapply(prop_match, sum)) == 
                           2)
  int_terms <- intersect(all_pairwise_ints, prop_pairwise)
  return(int_terms)
}

AIC2 <-  function (obj) 
{
  n <- length(na.omit(obj$y))
  p <- length(na.omit(obj$coef))
  mu_hat <- fitted(obj)
  sigma_hat <- sqrt(sum(obj$residuals^2)/(obj$df.residual) * 
                      (n - p)/n)
  ll <- sum(dnorm(obj$y, mu_hat, sigma_hat, log = TRUE))
  np <- p + 1
  aic <- -2 * ll + 2 * np
  return(aic)
}

# 4 species example theta not equal 1

#Creating dataset of all interaction terms
data(sim2)
prop <- paste0('p',1:4)
FG <- c('G','G','H','H')
dummy2 <- DI_data_I_ij(sim2, prop, FG=FG)

# AvG model
mod_AV <- lm(response ~ p1 + p2 + p3 + p4 + 
               AV + 0, data=dummy2) 
mod_AV_theta <- DI_theta(mod_AV, DImodel = 'AV', prop=prop, nSpecies = length(prop), family='gaussian')
mod_AV_theta_DI <- DI(y='response', prop=prop, DImodel = 'AV',FG=FG, data=sim2, estimate_theta=T)
mod_AV_theta
mod_AV_theta_DI

# ADD model
mod_ADD <- lm(response ~ p1 + p2 + p3 + p4 + 
                        p1_add + p2_add + 
                        p3_add + p4_add + 0, data=dummy2) 
mod_ADD_theta <- DI_theta(mod_ADD, DImodel = 'ADD', prop=prop, nSpecies = length(prop), family='gaussian')
mod_ADD_theta_DI <- DI(y='response', prop=prop, DImodel = 'ADD',FG=FG, data=sim2, estimate_theta=T)
mod_ADD_theta
mod_ADD_theta_DI

# FG model
mod_FG <- lm(response ~ p1 + p2 + p3 + p4 + 
                FG_bfg_G_H + FG_wfg_G +
               FG_wfg_H + 0, data=dummy2) 
mod_FG_theta <- DI_theta(mod_FG, DImodel = 'FG', prop=prop, FG=FG,nSpecies = length(prop), family='gaussian')
mod_FG_theta_DI <- DI(y='response', prop=prop, DImodel = 'FG',FG=FG, data=sim2, estimate_theta=T)
mod_FG_theta
mod_FG_theta_DI

# FULL model
mod_FULL <- lm(response ~ p1 + p2 + p3 + p4 + 
                         p1_p2 + p1_p3 + p1_p4 + 
                         p2_p3 + p2_p4 + p3_p4 +0, data=dummy2) 
mod_FULL_theta <- DI_theta(mod_FULL, DImodel = 'FULL', prop=prop, nSpecies = length(prop), family='gaussian')
mod_FULL_theta_DI <- DI(y='response', prop=prop, DImodel = 'FULL',FG=FG, data=sim2, estimate_theta=T)
mod_FULL_theta
mod_FULL_theta_DI

# Note: Theta stays same across both parameterizations 
# Note: In the new parameterization the interaction coefficients 
# scaled by a factor of s(s-1)/2s^2*theta where s is number of species

scale_factor <- function(s, theta){
  (s*(s-1))/(2*(s^(2*theta)))
}
mod_AV_theta_DI$coefficients['AV']*scale_factor(s=4, theta=mod_AV_theta_DI$coefficients['theta'])
mod_AV_theta$coefficients['AV']

# Note Predictions are unaffected too
all.equal(predict(mod_AV_theta_DI),predict(mod_AV_theta))




#9 species example

data(sim5)

#Creating dataset of all interaction terms
prop <- paste0('p',1:9)
FG <- c('G','G','G','G','G','H','H','L','L')
dummy3 <- DI_data_I_ij(sim5, prop, FG=FG)

#Creating partial formula for lm object
iden <- paste(prop, collapse = ' + ')

# AvG model
ints <- 'AV'
formula <- as.formula(paste('response ~ 0 +', iden, ints, sep =' + '))

mod_AV <- lm(formula, data=dummy3) 
mod_AV_theta <- DI_theta(mod_AV, DImodel = 'AV', prop=prop, nSpecies = length(prop), family='gaussian')
mod_AV_theta_DI <- DI(y='response', prop=prop, DImodel = 'AV',FG=FG, data=sim5, estimate_theta=T)
mod_AV_theta
mod_AV_theta_DI

# ADD model
ints <- paste(dummy3 %>% select(contains('_add')) %>% colnames(.), collapse = ' + ')
formula <- as.formula(paste('response ~ 0 +', iden, ints, sep =' + '))

mod_ADD <- lm(formula, data=dummy3) 
mod_ADD_theta <- DI_theta(mod_ADD, DImodel = 'ADD', FG=FG, prop=prop, nSpecies = length(prop), family='gaussian')
mod_ADD_theta_DI <- DI(y='response', prop=prop, DImodel = 'ADD',FG=FG, data=sim5, estimate_theta=T)
mod_ADD_theta
mod_ADD_theta_DI

# FG model
ints <- paste(dummy3 %>% select(contains('FG_')) %>% colnames(.), collapse = ' + ')
formula <- as.formula(paste('response ~ 0 +', iden, ints, sep =' + '))

mod_FG <- lm(formula, data=dummy3) 
mod_FG_theta <- DI_theta(mod_FG, DImodel = 'FG', FG=FG, prop=prop, nSpecies = length(prop), family='gaussian')
mod_FG_theta_DI <- DI(y='response', prop=prop, DImodel = 'FG',FG=FG, data=sim5, estimate_theta=T)
mod_FG_theta
mod_FG_theta_DI

# FULL model
ints <- paste(dummy3 %>% select(contains('_') & !(contains('FG_') | contains('_add'))) %>% colnames(.), collapse = ' + ')
formula <- as.formula(paste('response ~ 0 +', iden, ints, sep =' + '))

mod_FULL <- lm(formula, data=dummy3) 
mod_FULL_theta <- DI_theta(mod_FULL, DImodel = 'FULL', prop=prop, nSpecies = length(prop), family='gaussian')
mod_FULL_theta_DI <- DI(y='response', prop=prop, DImodel = 'FULL',FG=FG, data=sim5, estimate_theta=T)
mod_FULL_theta
mod_FULL_theta_DI

#6 species example

data(sim4)

#Creating dataset of all interaction terms
prop <- paste0('p',1:6)
FG <- c('G','G','H','H','L','L')
dummy4 <- DI_data_I_ij(sim4, prop, FG=FG)

#Creating partial formula for lm object
iden <- paste(prop, collapse = ' + ')

# AvG model
ints <- 'AV'
formula <- as.formula(paste('response ~ 0 +', iden, ints, sep =' + '))

mod_AV <- lm(formula, data=dummy4) 
mod_AV_theta <- DI_theta(mod_AV, DImodel = 'AV', prop=prop, nSpecies = length(prop), family='gaussian')
mod_AV_theta_DI <- DI(y='response', prop=prop, DImodel = 'AV',FG=FG, data=sim4, estimate_theta=T)
mod_AV_theta
mod_AV_theta_DI

# ADD model
ints <- paste(dummy4 %>% select(contains('_add')) %>% colnames(.), collapse = ' + ')
formula <- as.formula(paste('response ~ 0 +', iden, ints, sep =' + '))

mod_ADD <- lm(formula, data=dummy4) 
mod_ADD_theta <- DI_theta(mod_ADD, DImodel = 'ADD', FG=FG, prop=prop, nSpecies = length(prop), family='gaussian')
mod_ADD_theta_DI <- DI(y='response', prop=prop, DImodel = 'ADD',FG=FG, data=sim4, estimate_theta=T)
mod_ADD_theta
mod_ADD_theta_DI

# FG model
ints <- paste(dummy4 %>% select(contains('FG_')) %>% colnames(.), collapse = ' + ')
formula <- as.formula(paste('response ~ 0 +', iden, ints, sep =' + '))

mod_FG <- lm(formula, data=dummy4) 
mod_FG_theta <- DI_theta(mod_FG, DImodel = 'FG', FG=FG, prop=prop, nSpecies = length(prop), family='gaussian')
mod_FG_theta_DI <- DI(y='response', prop=prop, DImodel = 'FG',FG=FG, data=sim4, estimate_theta=T)
mod_FG_theta
mod_FG_theta_DI

# FULL model
ints <- paste(dummy4 %>% select(contains('_') & !(contains('FG_') | contains('_add'))) %>% colnames(.), collapse = ' + ')
formula <- as.formula(paste('response ~ 0 +', iden, ints, sep =' + '))

mod_FULL <- lm(formula, data=dummy4) 
mod_FULL_theta <- DI_theta(mod_FULL, DImodel = 'FULL', prop=prop, nSpecies = length(prop), family='gaussian')
mod_FULL_theta_DI <- DI(y='response', prop=prop, DImodel = 'FULL',FG=FG, data=sim4, estimate_theta=T)
mod_FULL_theta
mod_FULL_theta_DI


# Using non-linear estimation to estimate theta
estimate_theta_NLS_AV<- function(data, prop, ij=T){
  iden <- prop
  
  if (ij){
    ints <- data %>% select(contains('_') & !(contains('FG_') | contains('_add'))) %>% colnames(.)
  } else {
    ints <- grep(':',colnames(data), value =T)
    new_ints <- str_replace(ints, ':','_')
    data <- data %>% rename_with(~ new_ints[which(ints == .x)], .cols = ints)
    ints <- new_ints
  }
  
  formula <- as.formula(paste('response ~ 0', paste(iden, collapse=' + '),
                              'AV', sep ='+'))
  
  mod_AV <- lm(formula, data=data) 
  svs <- as.vector(c(coef(mod_AV), 1))
  names(svs) <- c(paste("b", 1:length(iden), sep=""),"dav", "theta")
  
  theta_idens <- paste("b", 1:length(iden), sep="", paste("*p",1:length(iden),sep=""))
  theta_ints <- paste("dav","*", sep="", paste(ints, 'theta', sep='^'))	
  
  form <- as.formula(paste('response ~ 0', paste(theta_idens, collapse=' + '),
                           paste(theta_ints, collapse=' + '), sep ='+'))
  
  model <- nls(form, data=data, start=svs, trace=F)
  return(model)
}

estimate_theta_NLS_FULL<- function(data, prop, ij=T){
  iden <- prop
  
  if (ij){
    ints <- data %>% select(contains('_') & !(contains('FG_') | contains('_add'))) %>% colnames(.)
  } else {
    ints <- grep(':',colnames(data), value =T)
    new_ints <- str_replace(ints, ':','_')
    data <- data %>% rename_with(~ new_ints[which(ints == .x)], .cols = ints)
    ints <- new_ints
  }
  
  formula <- as.formula(paste('response ~ 0', paste(iden, collapse=' + '),
                              paste(ints, collapse=' + '), sep ='+'))
  
  mod_FULL <- lm(formula, data=data) 
  svs <- as.vector(c(coef(mod_FULL), 1))
  names(svs) <- c(paste("b", 1:length(iden), sep=""),paste("d", 1:length(ints), sep=""), "theta")
  
  theta_idens <- paste("b", 1:length(iden), sep="", paste("*p",1:length(iden),sep=""))
  theta_ints <- paste("d", 1:length(ints),"*", sep="", paste(ints, 'theta', sep='^'))	
  
  form <- as.formula(paste('response ~ 0', paste(theta_idens, collapse=' + '),
                           paste(theta_ints, collapse=' + '), sep ='+'))
  
  model <- nls(form, data=data, start=svs, trace=F)
  return(model)
}

species4_NLS <- estimate_theta_NLS(sim2, prop=paste0('p',1:4), FG=c('G','G','H','H'))
species4_NLS$AV_I
species4_NLS$AV

species6_NLS <- estimate_theta_NLS(sim4, prop=paste0('p',1:6), FG=c('G','G','H','H','L','L'))
species6_NLS$AV_I
species6_NLS$AV

species9_NLS <- estimate_theta_NLS(sim5, prop=paste0('p',1:9), FG=c('G','G','G','G','G','H','H','L','L'))
species9_NLS$AV_I
species9_NLS$AV

# It can be seen that the theta coefficent estimate is same when
# estimating with non-linear least squares too and that the interaction
# coefficients are scaled by a factor of s(s-1)/2s^2*theta where s is 
# number of species
