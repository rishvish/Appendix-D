
# Main function to estimate theta, added one extra paramter which asks for upper
# boundary for searching theta. If not specified, defaults to 1.5
DI_theta <- function (obj, DImodel, FGnames, prop, nSpecies, family, lower_boundary=0.01, upper_boundary=1.5) 
{
  if (missing(FGnames)) {
    FGnames <- NULL
  }
  mm <- model.matrix(obj)
  int_terms <- switch(EXPR = DImodel, AV = grep("AV", 
                                                colnames(mm)), E = grep("E", colnames(mm)), FG = grep("FG", 
                                                                                                      colnames(mm)), ADD = grep("_add", colnames(mm)), 
                      FULL = get_int_terms_FULL(mm_names = colnames(mm), prop_names = prop))
  upper_boundary <-upper_boundary
  theta_info <- get_theta_info(lower_boundary = lower_boundary, upper_boundary = upper_boundary, 
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
    data_theta_E_AV <- obj$data
    data_theta <- data.frame(mm)
    # Had a question on why this is being done?
    data_theta[1:nSpecies] <- data_theta[1:nSpecies]^theta_hat
    new_E_AV <- DI_data_E_AV_internal(prop = 1:nSpecies, 
                                      data = data_theta)
    data_theta_E_AV$E <- new_E_AV$E
    data_theta_E_AV$AV <- new_E_AV$AV
    mod_theta <- glm(formula(obj), family = family, data = data_theta_E_AV)
  }
  else if (DImodel == "FG") {
    data_theta_FG <- obj$data
    data_theta <- data.frame(mm)
    data_theta[1:nSpecies] <- data_theta[1:nSpecies]^theta_hat
    new_FG <- DI_data_FG_internal(prop = 1:nSpecies, FG = FGnames, 
                                  data = data_theta)
    FG_ <- new_FG$FG
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
    new_ADD <- DI_data_ADD_theta(prop = 1:nSpecies, data = data_theta, 
                                 theta = theta_hat)
    ADD_cols_in_the_data <- int_terms
    if (length(ADD_cols_in_the_data) > 0) {
      j <- 1
      for (i in ADD_cols_in_the_data) {
        data_theta[, i] <- new_ADD$ADD_theta[, j]
        j <- j + 1
      }
    }
    old_formula <- formula(obj)
    new_formula <- paste0(old_formula[2], " ~ ", old_formula[3])
    data_theta_ADD <- data_theta
    data_theta_ADD$y <- obj$y
    names(data_theta_ADD)[length(names(data_theta_ADD))] <- paste(old_formula[2])
    mod_theta <- glm(as.formula(new_formula), family = family, 
                     data = data_theta_ADD)
  }
  else {
    mm[, int_terms] <- (mm[, int_terms])^theta_hat
    colnames(mm) <- gsub("`", "", colnames(mm))
    names_mm <- paste0("`", colnames(mm)[int_terms], 
                       "`")
    resp_name <- paste(formula(obj))[2]
    ndata <- data.frame(obj$data[, resp_name], mm, check.names = FALSE)
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
  return(mod_theta)
}



# This is the function that does the profile likelihood 
# Copied as is from DImodels
get_theta_info <- function (lower_boundary, upper_boundary, DImodel, obj, family, int_terms, nSpecies, 
                            FGnames) 
{
  optimum <- optimize(proflik_theta, interval = c(lower_boundary, upper_boundary), 
                        maximum = TRUE, DImodel = DImodel, obj = obj, family = family, 
                        int_terms = int_terms, nSpecies = nSpecies, FGnames = FGnames)
  theta_hat <- optimum$maximum
  theta_grid <- seq(lower_boundary, upper_boundary + 1, length = 100)
  proflik_theta_vec <- Vectorize(proflik_theta, "theta")
  profile_loglik <- proflik_theta_vec(theta = theta_grid, obj = obj, 
                                      family = family, int_terms = int_terms, DImodel = DImodel, 
                                      nSpecies = nSpecies, FGnames = FGnames)
  profile_loglik = data.frame(prof = c(profile_loglik, optimum$objective), 
                              grid = c(theta_grid,theta_hat))
  profile_loglik = profile_loglik[order(profile_loglik$grid),]
  return(list(theta_hat = theta_hat, profile_loglik = profile_loglik))
}

# The actual function that is optimized to get theta value
# Copied as is from DImodels
proflik_theta <- function (theta, obj, family, int_terms, DImodel, nSpecies, FGnames) 
{
  mm <- model.matrix(obj)
  if (DImodel %in% c("E", "AV")) {
    data_theta_E_AV <- obj$data
    data_theta <- data.frame(mm)
    data_theta[1:nSpecies] <- data_theta[1:nSpecies]^theta
    new_E_AV <- DI_data_E_AV_internal(prop = 1:nSpecies, 
                                      data = data_theta)
    data_theta_E_AV$E <- new_E_AV$E 
    data_theta_E_AV$AV <-new_E_AV$AV
    
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
    data_theta_FG <- obj$data
    data_theta <- data.frame(mm)
    data_theta[1:nSpecies] <- data_theta[1:nSpecies]^theta
    new_FG <- DI_data_FG_internal(prop = 1:nSpecies, FG = FGnames, 
                                  data = data_theta)
    FG_ <- new_FG$FG
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
    new_ADD <- DI_data_ADD_theta(prop = 1:nSpecies, data = data_theta, 
                                 theta = theta)
    ADD_cols_in_the_data <- int_terms
    if (length(ADD_cols_in_the_data) > 0) {
      j <- 1
      for (i in ADD_cols_in_the_data) {
        data_theta[, i] <- new_ADD$ADD_theta[, j]
        j <- j + 1
      }
    }
    old_formula <- formula(obj)
    new_formula <- paste0(old_formula[2], " ~ ", old_formula[3])
    data_theta_ADD <- data_theta
    data_theta_ADD$y <- obj$y
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
    mm[, int_terms] <- (mm[, int_terms])^theta
    ndata <- obj$data
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

# The following functions are same as DImodels and are there just to make 
# sure everything runs smoothly
DI_data_E_AV_internal<- function (prop, data) 
{
  if (!is.character(prop)) {
    prop <- get_P_indices(prop = prop, data = data)$prop
  }
  nSpecies <- length(prop)
  if (nSpecies <= 1) 
    stop("must have at least 2 species to fit DI models")
  nComb <- choose(nSpecies, 2)
  pairwise_fmla <- as.formula(paste("~", "0+", 
                                    "(", paste(prop, collapse = "+"), ")^2"))
  normE <- 2 * nSpecies/(nSpecies - 1)
  Ematrix <- model.matrix(pairwise_fmla, data = data)
  if (nSpecies > 2) {
    AV <- rowSums(Ematrix[, (nSpecies + 1):ncol(Ematrix)])
    E <- normE * AV
    even_flag <- FALSE
    return(list(AV = AV, E = E, even_flag = even_flag))
  }
  else {
    even_flag <- TRUE
    return(list(even_flag = even_flag))
  }
}

DI_data_FG_internal <- function (prop, FG, data) 
{
  FG_name_check(FG = FG)
  nfg <- length(unique(FG))
  n_check <- choose(nfg, 2) + nfg
  Pind <- get_P_indices(prop = prop, data = data)$Pind
  if (any(!is.character(FG))) 
    stop("FG argument takes character strings with functional", 
         " group names referring to each species, in order")
  fg_index <- FG
  fg_index_names <- levels(as.factor(fg_index))
  testdata <- data[, Pind]
  prop <- colnames(testdata)
  prop <- paste(prop, fg_index, sep = "_infg_")
  testdata2 <- testdata
  colnames(testdata2) <- prop
  nSpecies <- length(prop)
  pairwise_fmla <- as.formula(paste("~", "0+", 
                                    "(", paste(prop, collapse = "+"), ")^2"))
  Pmatrix <- model.matrix(pairwise_fmla, data = testdata2)
  FGmatrix_raw <- Pmatrix[, (nSpecies + 1):ncol(Pmatrix)]
  FGrawnames <- colnames(FGmatrix_raw)
  FGmatch <- list()
  for (i in 1:nfg) {
    FGmatch[[i]] <- grep(fg_index_names[i], FGrawnames)
  }
  FGnames <- rep(NA, length(FGrawnames))
  for (i in 1:nfg) {
    for (j in 1:nfg) {
      if (i < j) {
        FGnames[FGmatch[[i]][FGmatch[[i]] %in% FGmatch[[j]]]] <- paste("bfg_", 
                                                                       fg_index_names[i], ":", fg_index_names[j], 
                                                                       sep = "")
      }
    }
  }
  FGmatch2 <- list()
  for (i in 1:nfg) {
    FGmatch2[[i]] <- grep(paste(fg_index_names[i], ":", 
                                sep = ""), FGrawnames)
  }
  for (i in 1:nfg) {
    helper_index <- FGmatch[[i]][FGmatch[[i]] %in% FGmatch2[[i]]]
    FGnames[helper_index[is.na(FGnames[helper_index])]] <- paste("wfg_", 
                                                                 fg_index_names[i], sep = "")
  }
  FGeffects <- levels(as.factor(FGnames))
  FG <- matrix(NA, nrow = nrow(FGmatrix_raw), ncol = length(FGeffects))
  for (i in 1:length(FGeffects)) {
    FG[, i] <- rowSums(as.matrix(FGmatrix_raw[, FGnames == 
                                                FGeffects[i]]))
  }
  FGeffects <- gsub(":", "_", FGeffects)
  colnames(FG) <- FGeffects
  if (any(table(fg_index) < 2)) {
    return(list(FG = FG))
  }
  else if (n_check != ncol(FG)) {
    stop("Expected ", n_check, " terms, but have ", 
         ncol(FG), ". Please give your functional groups a different name.", 
         " One or more of the options are being used internally.", 
         " Perhaps use upper-case names?")
  }
  return(list(FG = FG))
}

FG_name_check <- function (FG) 
{
  cond1 <- length(grep(":", FG)) > 0
  cond2 <- any(FG == "_")
  cond3 <- any(FG == "i")
  cond4 <- any(FG == "n")
  cond5 <- any(FG == "f")
  cond6 <- any(FG == "g")
  cond7 <- any(FG == "_i")
  cond8 <- any(FG == "in")
  cond9 <- any(FG == "nf")
  cond10 <- any(FG == "fg")
  cond11 <- any(FG == "g_")
  cond12 <- any(FG == "_in")
  cond13 <- any(FG == "inf")
  cond14 <- any(FG == "nfg")
  cond15 <- any(FG == "fg_")
  cond16 <- any(FG == "_inf")
  cond17 <- any(FG == "infg")
  cond18 <- any(FG == "nfg_")
  cond19 <- any(FG == "_infg")
  cond20 <- any(FG == "infg_")
  cond21 <- any(FG == "_infg_")
  if (cond1 | cond2 | cond3 | cond4 | cond5 | cond6 | cond7 | 
      cond8 | cond9 | cond10 | cond11 | cond12 | cond13 | cond14 | 
      cond15 | cond16 | cond17 | cond18 | cond19 | cond20 | 
      cond21) {
    stop("Please give your functional groups a different name.", 
         " Names should not include colons (':'), or any single or multiple", 
         " character combination of the expression '_infg_'.", 
         " This expression is reserved for computing functional groups internally.")
  }
}

get_int_terms_FULL <- function (mm_names, prop_names) 
{
  mm_names <- gsub("`", "", mm_names)
  all_pairwise_ints <- which(lengths(regmatches(mm_names, gregexpr(":", 
                                                                   mm_names))) == 1)
  prop_match <- lapply(strsplit(mm_names, ":"), function(x) {
    #print(x)
    x %in% prop_names
  })
  prop_pairwise <- which(unlist(lapply(prop_match, sum)) == 
                           2)
  int_terms <- intersect(all_pairwise_ints, prop_pairwise)
  return(int_terms)
}

get_P_indices <- function (prop, data) 
{
  if (!is.character(prop)) {
    Pind <- prop
    prop <- names(data[, prop])
  }
  else {
    vec_grep <- Vectorize(grep, "pattern")
    Pind <- as.numeric(vec_grep(paste("\\<", prop, 
                                      "\\>", sep = ""), names(data)))
  }
  return(list(Pind = Pind, prop = prop))
}

DI_data_ADD_theta <- function (prop, data, theta) 
{
  Pind <- get_P_indices(prop = prop, data = data)$Pind
  nSpecies <- length(prop)
  P_matrix <- data[, Pind]
  Pi_theta <- P_matrix^theta
  ADD_vars_theta <- matrix(NA, ncol = nSpecies, nrow = nrow(data))
  for (i in 1:nSpecies) {
    sum_Pj_theta <- apply(Pi_theta[, -i], 1, sum)
    ADD_vars_theta[, i] <- Pi_theta[, i] * sum_Pj_theta
  }
  ADD_vars_theta <- as.data.frame(ADD_vars_theta)
  prop_names <- names(data)[Pind]
  names(ADD_vars_theta) <- paste(prop_names, "_add", 
                                 sep = "")
  return(list(ADD_theta = ADD_vars_theta))
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


theta_CI <- function(obj, conf = .95, n = 100) {
  threshold <- max(obj$profile_loglik$prof) - qchisq(conf, 1)/2
  if(threshold < min(obj$profile_loglik$prof) | threshold > max(obj$profile_loglik$prof)) {
    stop("CI cannot be computed. This is because the profile log-likelihood function is flat or displays unusual behaviour at the interval theta = (0.01, 2.5).") 
  }
  CI_finder <- approxfun(x = obj$profile_loglik$grid,
                         y = obj$profile_loglik$prof - threshold)
  CI <- rootSolve::uniroot.all(CI_finder, interval = range(obj$profile_loglik$grid), n = n)
  alpha <- 1 - conf
  names(CI) <- c("lower","upper")
  return(CI)
}

get_CI <- function(mod, n = 100){
  CI <- c(0,0)
  names(CI) <- c('lower','upper')
  tryCatch(
    expr = {
      CI <- theta_CI(mod, n = n)
    },
    error = function(e){ 
      CI <- c(0,0)
    },
    warning = function(w){
    },
    finally = {
    }
  )
  
  if (length(CI)!=2){
    CI <- c(0,0)
  }
  
  return (CI)
}

get_WaldCI <- function(mod, theta_est, FG, DImodel, nSpecies, species, conf = 0.95){
  CI <- c(0,0)
  names(CI) <- c('lower','upper')
  tryCatch(
    expr = {
      mm <- model.matrix(mod)
      int_terms <- switch(EXPR = DImodel, AV = grep("AV", 
                                                    colnames(mm)), E = grep("E", colnames(mm)), FG = grep("FG", 
                                                                                                          colnames(mm)), ADD = grep("_add", colnames(mm)), 
                          FULL = get_int_terms_FULL(mm_names = colnames(mm), prop_names = species))
      if (missing(FG)) {
        FG <- NULL
      }
      
      hess <- numDeriv::hessian(proflik_theta, x = theta_est, 
                                obj = mod, family = 'gaussian', 
                                int_terms = int_terms, 
                                DImodel = DImodel,
                                nSpecies = nSpecies, 
                                FGnames = FG)
      if((hess > 0)){
        SE <- sqrt(diag(solve(hess)))
      } else if ((hess) < 0){
        SE <- sqrt(diag(solve(-hess)))
      } else {
        stop('Error hessian is 0')
      }
      
      CI <- theta_est + c(-1, 1)*SE*1.96
      
    },
    error = function(e){ 
      #print('e')
      CI <- c(0,0)
    },
    warning = function(w){
    },
    finally = {
    }
  )
  
  return (CI)
}

