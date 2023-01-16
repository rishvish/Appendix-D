required_packages <- c('tidyverse', 'DImodels', 'ggpubr', 'readr')

for (package in required_packages){
  if (!(package %in% installed.packages())){
    install.packages(package)
  }
}

library(tidyverse)
library(DImodels)
library(ggpubr)

dir <- dirname(rstudioapi::getSourceEditorContext()$path)
#source(paste0(dir,'/Pipeline_Helper_Code.R'))

#' Function to create the experiment design, which will constitute the communities to use in the simulation
#'
#' @param no_species Number of species in the experiment  
#' @param reps Number of repetitions for each community in the design
#' @param treatments Vector of treatments, if any
#' @param random Whether to perform a random selection of communities at each richness level. Useful for situations when there are more than four species in the experiment
#' @param ncomms Number of random communities to select at each level of richness
#'
#' @return A dataframe consisting of communities in the experiment
create_communities <- function(no_species=4, reps=1, treatments= NA, 
                               random=FALSE, ncomms=NA){
  if(missing(no_species) | !is.numeric(no_species)){
    stop('Number of species is mandatory and should be numeric')
  }
  
  if(!is.numeric(reps)){
    stop('Number of repititions should be numeric')
  }
  
  if(!is.na(ncomms) & !is.numeric(ncomms)){
    stop('Number of random communities should be numeric')
  }
  
  if(!any(is.na(treatments)) & !is.vector(treatments)){
    stop('Treatments should be a vector of different treatments')
  }
  
  if(!is.logical(random)){
    stop('Random should be a boolean flag indicating whether or not to perform random selection')
  }
  
  if(random==TRUE & is.na(ncomms)){
    stop('Specify the number of random communities to select at each richness level')
  }
  
  if(!is.na(ncomms) & ncomms < (no_species/2)){
    stop(glue('Number of random communities should be atleast { round(no_species/2) + 1 }'))
  }
    
  species <- paste0('p',1:no_species)
  grid <- expand.grid(rep(list(0:1), no_species))
  colnames(grid) <- species
  grid <- grid/rowSums(grid)
  props <- as_tibble(grid[-1,]) %>% mutate(Richness=rowSums(.[,]!=0))
  
  if (random==TRUE & !is.na(ncomms)){
    props_random <- props %>% filter(Richness==1)
    for (i in 2:length(species)){
      samp <- props %>% filter(Richness==i) %>% slice_sample(n=ncomms) 
      
      which_miss <- which(colSums(samp[, species])==0)
      
      while(length(which_miss)!=0){
        srow <- props %>% filter(Richness==i & get(names(which_miss[1]))!=0) %>% 
          slice_sample(n=1)
        samp[sample(nrow(samp), size=1),] <- srow
        which_miss <- which(colSums(samp[, species])==0)
      }
      props_random <- bind_rows(props_random, samp )
    }
    props <- props_random
  }
  
  props <- props %>% slice(rep(1:n(), times = reps*length(treatments))) %>% 
    mutate(Community=rownames(.), rep=rep(1:reps, each=nrow(props)*length(treatments)),
           Treatment=rep(treatments, each=nrow(props), times=reps)) %>% 
    dplyr::select(Community, Richness, Treatment, everything()) 
  
  if (any(is.na(treatments))){
    props <- props %>% dplyr::select(-Treatment)  
  }
  
  return(props)
}


#' Simulate Response Function from dataset generate by create_communities
#'
#' @param data The dataframe consisting of species communities in the design
#' @param theta Vector of theta values
#' @param sds Vector of sigma values for the error term sd
#' @param seed Initial seed for the random normal error term
#' @param num_datasets Number of datasets for each unique theta, sd and seed combination
#' @param model The true underlying interaction structure used to simulate the response. Should be one of AV, FG, ADD or FULL.
#' @param FG A vector with length equal to the number of species in the experiment giving their functional groupings
#' @param coeff Vector of coefficients for the DImodel used to simulate the response
#' @param fileName Name of the csv file in which the simulated datasets will be saved.
#'
#' @return String giving name of csv file contatining the simulated datasets
simulate_response <- function(data, theta, sds, seed=100, 
                              num_datasets=1, model, FG=NA, 
                              coeff, fileName = NA){
  
  if(missing(data)){
    stop('Please specify a dataframe with species proportions')
  }
  
  if(missing(theta) | ! is.numeric(theta)){
    stop('Thetas is mandatory and should be numeric or a numeric vector with values in between 0 and 1.5')
  }
  
  if(missing(sds) | ! is.numeric(sds)){
    stop('Sds is mandatory and should be numeric or a numeric vector')
  }
  
  if(missing(model) | ! is.character(model)){
    stop("Model is mandatory and should be one from c('ADD','AV','FG','FULL')")
  }
  
  species <- grep('^p',colnames(data), value=TRUE)
  
  if (!all(near(rowSums(data[,species]),1))){
    stop('Species proportions don\'t sum to one')
  }
  
  if(model=='FG' & (missing(FG) | ! is.vector(FG, mode='character') | length(FG)!=length(species))){
    stop('FG argument is mandatory for functional group model and takes character strings with functional group names referring to each species, in order')
  }
  
  #Checking if files already exists and we don't overwrite an existing one
  if (is.na(fileName)){
    fileName <- paste0("SimulatedData",model,".csv")
  } else {
    fileName <- paste0(fileName,".csv")
  }
    
  i=1
  while(file.exists(fileName)){
    file <- strsplit(fileName, '.csv')[[1]][1]
    fileName <- paste0(file,' (',i,')',".csv")
    i=i+1
  }
  
  # if (is.na(data$Treatment)){
  #   data <- data %>% select(-Treatment)
  # }
  oldw <- getOption('warn')
  for (i in 1:length(theta)){
    options(warn=-1)  #Supress warnings
    
    
    if (model=='FG'){
      sim_data_iter <- cbind(data,DI_data(prop = species, theta = theta[i],FG=FG,data = data, what='FG'))
      sim_data_iter$theta <- theta[i] 
      sim_data_iter$model <- model
      
      X <- sim_data_iter %>% select(starts_with('p', ignore.case = F) | contains('fg') | contains('Treat')) 
      X <- model.matrix(~.-1, data=X)
      
      if(length(coeff)!=ncol(X)){
        print(colnames(X))
        stop(paste0('Number of coefficents should be equal to ', ncol(X),' for ',model,' model'))
      }
      
      sim_data_iter$Res <- as.numeric(X %*% coeff)
      
    } else if (model=='FULL'){
      
      sim_data_iter <- cbind(data,DI_data(prop = species, theta = theta[i],data = data, what='FULL'))
      sim_data_iter$theta <- theta[i] 
      sim_data_iter$model <- model
      
      sim_data_iter <- sim_data_iter %>% rename_at(vars(contains(":")), funs(toupper(.)))
      
      X <- sim_data_iter %>% select(starts_with('p') | contains('Treat')) 
      X <- model.matrix(~.-1, data=X)
      
      if(length(coeff)!=ncol(X)){
        print(colnames(X))
        stop(paste0('Number of coefficents should be equal to ', ncol(X),' for ',model,' model'))
      }
      
      sim_data_iter$Res <- as.numeric(X %*% coeff)
      
    } else if (model=='AV'){
      AV <- DI_data(prop = species, theta = theta[i],data = data, what='AV')
      sim_data_iter <- cbind(data,AV)
      sim_data_iter$theta <- theta[i] 
      sim_data_iter$model <- model
      
      X <- sim_data_iter %>% select(starts_with('p') | contains('Treat'), AV) 
      X <- model.matrix(~.-1, data=X)
      
      if(length(coeff)!=ncol(X)){
        print(colnames(X))
        stop(paste0('Number of coefficents should be equal to ', ncol(X),' for ',model,' model'))
      }
      
      sim_data_iter$Res <- as.numeric(X %*% coeff)
      
    } else if (model=='E'){
      E <- DI_data(prop = species, theta = theta[i],data = data, what='E')
      sim_data_iter <- cbind(data,E)
      sim_data_iter$theta <- theta[i] 
      sim_data_iter$model <- model
      
      X <- sim_data_iter %>% select(starts_with('p') | contains('Treat'), E) 
      X <- model.matrix(~.-1, data=X)
      
      if(length(coeff)!=ncol(X)){
        print(colnames(X))
        stop(paste0('Number of coefficents should be equal to ', ncol(X),' for ',model,' model'))
      }
      
      sim_data_iter$Res <- as.numeric(X %*% coeff)
      
    } else if (model=='ADD'){
      sim_data_iter <- cbind(data,DI_data(prop = species, theta = theta[i],data = data, what='ADD'))
      sim_data_iter$theta <- theta[i] 
      sim_data_iter$model <- model
      
      X <- sim_data_iter %>% select(starts_with('p') | contains('add') | contains('Treat')) 
      X <- model.matrix(~.-1, data=X)
      
      if(length(coeff)!=ncol(X)){
        print(colnames(X))
        stop(paste0('Number of coefficents should be equal to ', ncol(X),' for ',model,' model'))
      }
      
      sim_data_iter$Res <- as.numeric(X %*% coeff)
      
    }
    
    for (j in 1:length(sds)){
      sim_data_iter$sd <- sds[j]
      for (k in 1:num_datasets){
        set.seed(seed+k)
        r <- rnorm(n=nrow(sim_data_iter), mean=0, sd=sds[j]) 
        sim_data_iter$seed <- seed+k
        sim_data_iter$Response <- round(sim_data_iter$Res + r, digits = 3)
        
        
        if(!file.exists(fileName)){
          write_csv(sim_data_iter %>% select(-Res),  fileName, col_names = TRUE) #Highly efficient and safeguards against system failure
        } else{
          write_csv(sim_data_iter %>% select(-Res),  fileName, append=TRUE) #Highly efficient and safeguards against system failure
        }
        
      }
      
      print(paste0('SD:',sds[j],', Theta:',theta[i],'. ', ((length(sds)*(i-1))+j)*100/(length(theta)*length(sds)),'% Work Done.'))
    }
  }
  
  options(warn=oldw)  ## Restore warnings
  print(paste0('Data created and written in ', fileName))
  return (fileName)
}



#' Estimate Theta from output of simulate_response
#'
#' @param data Dataframe created by simulate_response function
#' @param FG Vector of functional groupings to use when fitting the FG interaction structure
#'
#' @return A dataframe with the estimated theta values and other statistics for each unique combination of theta, sd and seed values in the data
estimate_theta <- function(data, FG){
  if(missing(data)){
    stop('Please specify a dataframe with species proportions')
  }
  
  species <- data %>% select(starts_with('p', ignore.case = F) & !contains('_add')) %>% colnames()
  #print(species)
  if (sum(grep('_add',colnames(data))) != 0){
    data <- data %>% select(-contains('_add'))
    
  }
  
  if('Treatment' %in% colnames(data)){
    print('Adding treatment to model')
    treat = 'Treatment'
  } else {
    treat = NA
  }
  
  if(missing(FG)){
    stop('FG argument is mandatory for functional group model and takes character strings with functional group names referring to each species, in order of apperance in data')
  }
  if(!is.vector(FG, mode='character') | length(FG)!=length(species)){
    stop('FG argument is a vector with same length as number of species and takes character strings with functional group names referring to each species, in order of appearance in data')
  }
  
  thetas <-  unique(data$theta)    #Uncomment this to check for all thetas, but it will take some time to run
  sds <- unique(data$sd)
  seeds <- unique(data$seed)
  models <- c('Average','Functional','Additive','Full Pairwise')  ## Checking for only these models as 
  ## they're the only ones having a contribution of theta
  
  ### Creating with null values dataframe first as it would be faster to just modify the rows
  ### rather than adding rows to an empty dataframe
  
  len <- length(seeds)*length(thetas)*length(sds)*length(models)   
  estimates <- tibble(true_theta = vector(mode = 'numeric', length=len),
                      sd = vector(mode='numeric', length=len),
                      seed = vector(mode = 'numeric', length=len),
                      model = vector(mode='character', length = len),
                      theta_est = vector(mode='numeric', length=len),
                      theta_lower= vector(mode='numeric', length=len),
                      theta_upper= vector(mode='numeric', length=len),
                      waldCI_lower = vector(mode='numeric', length=len),
                      waldCI_uppper = vector(mode='numeric', length=len),
                      aic= vector(mode='numeric', length=len),
                      bic= vector(mode='numeric', length=len))
  
  ## Main execution loop
  for (i in 1:length(thetas)){
    upper_boundary <- ifelse(thetas[i]>1.1, 2.5, 1.5)
    for (j in 1:length(sds)){
      for (k in 1:length(seeds)){
        data_iter <- data %>% filter(theta==thetas[i] & sd==sds[j] & seed==seeds[k])  ## Taking a specific subset of data
        
        modAV <- DI(y='Response', prop=species, treat = treat, DImodel = 'AV', data=data_iter)
        cust_modAV <- DI_theta(obj=modAV, DImodel = 'AV', prop=species, nSpecies= length(species), family='gaussian', lower_boundary = 0.01, upper_boundary = upper_boundary)
        AV_CI <- get_CI(cust_modAV, n =10000)
        thetaAV <- cust_modAV$coefficients['theta']
        AV_WaldCI <- get_WaldCI(cust_modAV, thetaAV, DImodel = 'AV', nSpecies = length(species), species = species)
        AV_aic <- AIC(cust_modAV)
        AV_bic <- BIC(cust_modAV)
        
        modFG <- DI(y='Response', prop=species, FG=FG,  treat = treat, DImodel = 'FG', data=data_iter)
        cust_modFG <- DI_theta(obj=modFG, DImodel = 'FG',FGnames=FG ,prop=species, nSpecies= length(species), family='gaussian', lower_boundary = 0.01, upper_boundary = upper_boundary)
        FG_CI <- get_CI(cust_modFG, n =10000)
        thetaFG <- cust_modFG$coefficients['theta']
        FG_WaldCI <- get_WaldCI(cust_modFG, thetaFG, FG = FG, DImodel = 'FG', nSpecies = length(species), species = species)
        FG_aic <- AIC(cust_modFG)
        FG_bic <- BIC(cust_modFG)
        
        modADD <- DI(y='Response', prop=species,  treat = treat, DImodel = 'ADD', data=data_iter)
        cust_modADD <- DI_theta(obj=modADD, DImodel = 'ADD', prop=species, nSpecies= length(species), family='gaussian', lower_boundary = 0.01, upper_boundary = upper_boundary)
        ADD_CI <- get_CI(cust_modADD, n =10000)
        thetaADD <- cust_modADD$coefficients['theta']
        ADD_WaldCI <- get_WaldCI(cust_modADD, thetaADD, DImodel = 'ADD', nSpecies = length(species), species = species)
        ADD_aic <- AIC(cust_modADD)
        ADD_bic <- BIC(cust_modADD)
        
        modFULL <- DI(y='Response', prop=species,  treat = treat, DImodel = 'FULL', data=data_iter)
        cust_modFULL <- DI_theta(obj=modFULL, DImodel = 'FULL', prop=species, nSpecies= length(species), family='gaussian', lower_boundary = 0.01, upper_boundary = upper_boundary)
        FULL_CI <- get_CI(cust_modFULL, n =10000)
        thetaFULL <- cust_modFULL$coefficients['theta']
        FULL_WaldCI <- get_WaldCI(cust_modFULL, thetaFULL, DImodel = 'FULL', nSpecies = length(species), species = species)
        FULL_aic <- AIC(cust_modFULL)
        FULL_bic <- BIC(cust_modFULL)
        
        index <- (length(models)*length(seeds)*length(sds)*(i-1) + length(models)*length(seeds)*(j-1) + length(models)*(k-1))
        
        estimates[index+1,] <- tibble_row(thetas[i], sds[j], seeds[k], models[1],
                                 thetaAV, AV_CI[1], AV_CI[2], 
                                 AV_WaldCI[1], AV_WaldCI[2],
                                 AV_aic, AV_bic)
        estimates[index+2,] <- tibble_row(thetas[i], sds[j], seeds[k], models[2],
                                 thetaFG, FG_CI[1], FG_CI[2],
                                 FG_WaldCI[1], FG_WaldCI[2],
                                 FG_aic, FG_bic)
        estimates[index+3,] <- tibble_row(thetas[i], sds[j], seeds[k], models[3],
                                 thetaADD, ADD_CI[1], ADD_CI[2],
                                 ADD_WaldCI[1], ADD_WaldCI[2],
                                 ADD_aic, ADD_bic)
        estimates[index+4,] <- tibble_row(thetas[i], sds[j], seeds[k], models[4],
                                 thetaFULL, FULL_CI[1], FULL_CI[2],
                                 FULL_WaldCI[1], FULL_WaldCI[2],
                                 FULL_aic, FULL_bic)
        
        print(paste0(round((index + 4)*100/nrow(estimates),3), '% Done.'))
          
      }
    }
  }
  return (estimates)
}


#' Efficient version of estimate theta function for examples with high number of species (will not fit the FULL pairwise interaction structure)
#'
#' @param data Dataframe created by simulate_response function
#' @param FG Vector of functional groupings to use when fitting the FG interaction structure
#'
#' @return A dataframe with the estimated theta values and other statistics for each unique combination of theta, sd and seed values in the data
estimate_theta_high <- function(data, FG){
  if(missing(data)){
    stop('Please specify a dataframe with species proportions')
  }
  
  species <- data %>% select(starts_with('p', ignore.case = F) & !contains('_add')) %>% colnames()
  #print(species)
  if (sum(grep('_add',colnames(data))) != 0){
    data <- data %>% select(-contains('_add'))
    
  }
  
  if(missing(FG)){
    stop('FG argument is mandatory for functional group model and takes character strings with functional group names referring to each species, in order of apperance in data')
  }
  if(!is.vector(FG, mode='character') | length(FG)!=length(species)){
    stop('FG argument is a vector with same length as number of species and takes character strings with functional group names referring to each species, in order of appearance in data')
  }
  
  thetas <-  unique(data$theta)    #Uncomment this to check for all thetas, but it will take some time to run
  sds <- unique(data$sd)
  seeds <- unique(data$seed)
  models <- c('Average','Functional','Additive')  ## Checking for only these models as 
  ## they're the only ones having a contribution of theta
  
  ### Creating with null values dataframe first as it would be faster to just modify the rows
  ### rather than adding rows to an empty dataframe
  
  len <- length(seeds)*length(thetas)*length(sds)*length(models)   
  estimates <- tibble(true_theta = vector(mode = 'numeric', length=len),
                      sd = vector(mode='numeric', length=len),
                      seed = vector(mode = 'numeric', length=len),
                      model = vector(mode='character', length = len),
                      theta_est = vector(mode='numeric', length=len),
                      theta_lower= vector(mode='numeric', length=len),
                      theta_upper= vector(mode='numeric', length=len),
                      aic= vector(mode='numeric', length=len),
                      bic= vector(mode='numeric', length=len))
  
  ## Main execution loop
  for (i in 1:length(thetas)){
    upper_boundary <- ifelse(thetas[i]>1.1, 2.5, 1.5)
    data_kter <- data %>% filter(theta==thetas[i])
    for (j in 1:length(sds)){
      data_jter <- data_kter %>% filter(sd==sds[j])
      for (k in 1:length(seeds)){
        data_iter <- data_jter %>% filter(seed==seeds[k])  ## Taking a specific subset of data
        
        modAV <- DI(y='Response', prop=species,  DImodel = 'AV', data=data_iter)
        cust_modAV <- DI_theta(obj=modAV, DImodel = 'AV', prop=species, nSpecies= length(species), family='gaussian', lower_boundary = 0.01, upper_boundary = upper_boundary)
        AV_CI <- get_CI(cust_modAV, n =10000)
        thetaAV <- cust_modAV$coefficients['theta']
        AV_aic <- AIC(cust_modAV)
        AV_bic <- BIC(cust_modAV)
        
        modFG <- DI(y='Response', prop=species, FG=FG, DImodel = 'FG', data=data_iter)
        cust_modFG <- DI_theta(obj=modFG, DImodel = 'FG',FGnames=FG ,prop=species, nSpecies= length(species), family='gaussian', lower_boundary = 0.01, upper_boundary = upper_boundary)
        FG_CI <- get_CI(cust_modFG, n =10000)
        thetaFG <- cust_modFG$coefficients['theta']
        FG_aic <- AIC(cust_modFG)
        FG_bic <- BIC(cust_modFG)
        
        modADD <- DI(y='Response', prop=species, DImodel = 'ADD', data=data_iter)
        cust_modADD <- DI_theta(obj=modADD, DImodel = 'ADD', prop=species, nSpecies= length(species), family='gaussian', lower_boundary = 0.01, upper_boundary = upper_boundary)
        ADD_CI <- get_CI(cust_modADD, n =10000)
        thetaADD <- cust_modADD$coefficients['theta']
        ADD_aic <- AIC(cust_modADD)
        ADD_bic <- BIC(cust_modADD)
        
        index <- (length(models)*length(seeds)*length(sds)*(i-1) + length(models)*length(seeds)*(j-1) + length(models)*(k-1))
        
        estimates[index+1,] <- tibble_row(thetas[i], sds[j], seeds[k], models[1],
                                          thetaAV, AV_CI[1], AV_CI[2],
                                          AV_aic, AV_bic)
        estimates[index+2,] <- tibble_row(thetas[i], sds[j], seeds[k], models[2],
                                          thetaFG, FG_CI[1], FG_CI[2],
                                          FG_aic, FG_bic)
        estimates[index+3,] <- tibble_row(thetas[i], sds[j], seeds[k], models[3],
                                          thetaADD, ADD_CI[1], ADD_CI[2],
                                          ADD_aic, ADD_bic)
        
        print(paste0(round((index + 3)*100/nrow(estimates),3), '% Done.'))
        
      }
    }
  }
  return (estimates)
}

#' Function to test robustness of theta across different functional groupings for the same data
#'
#' @param data Dataframe created by simulate_response function
#' @param FGs named vector or list giving the functional groupings to test
#'
#' @return A dataframe with the estimated theta values and other statistics for each unique combination of theta, sd and seed values in the data
test_FG <- function(data, FGs){
  if(missing(data)){
    stop('Please specify a dataframe with species proportions')
  }
  
  species <- data %>% select(starts_with('p', ignore.case = F) & !contains('_add')) %>% colnames()
  #print(species)
  if (sum(grep('_add',colnames(data))) != 0){
    data <- data %>% select(-contains('_add'))
    
  }
  
  if(missing(FGs)){
    stop('FG argument is mandatory for functional group model and takes character strings with functional group names referring to each species, in order of apperance in data')
  }
  # if(!is.vector(FG, mode='character') | length(FG)!=length(species)){
  #   stop('FG argument is a vector with same length as number of species and takes character strings with functional group names referring to each species, in order of appearance in data')
  # }
  
  thetas <-  unique(data$theta)    #Uncomment this to check for all thetas, but it will take some time to run
  sds <- unique(data$sd)
  seeds <- unique(data$seed)
  
  
  ### Creating with null values dataframe first as it would be faster to just modify the rows
  ### rather than adding rows to an empty dataframe
  
  len <- length(seeds)*length(thetas)*length(sds)*length(FGs)   
  estimates <- tibble(true_theta = vector(mode = 'numeric', length=len),
                      sd = vector(mode='numeric', length=len),
                      seed = vector(mode = 'numeric', length=len),
                      FG = vector(mode='character', length = len),
                      theta_est = vector(mode='numeric', length=len),
                      theta_lower= vector(mode='numeric', length=len),
                      theta_upper= vector(mode='numeric', length=len),
                      AIC = vector(mode='numeric', length=len),
                      BIC = vector(mode='numeric', length=len))
  ## Main execution loop
  for (i in 1:length(thetas)){
    upper_boundary <- ifelse(thetas[i]>1.1, 2.5, 1.5)
    for (j in 1:length(sds)){
      for (k in 1:length(seeds)){
        data_iter <- data %>% filter(theta==thetas[i] & sd==sds[j] & seed==seeds[k])  ## Taking a specific subset of data
        
        for (FG in 1:length(FGs)){
          group = FGs[[FG]]
          modFG <- DI(y='Response', prop=species, FG=group, DImodel = 'FG', data=data_iter)
          cust_modFG <- DI_theta(obj=modFG, DImodel = 'FG',FGnames=group , prop=species, nSpecies= length(species), family='gaussian', lower_boundary = 0.01, upper_boundary = upper_boundary)
          FG_CI <- get_CI(cust_modFG, n =10000)
          thetaFG <- cust_modFG$coefficients['theta']
          aic <- cust_modFG$aic
          bic <- BIC(cust_modFG)
          
          index <- (length(FGs)*length(seeds)*length(sds)*(i-1) + length(FGs)*length(seeds)*(j-1) + length(FGs)*(k-1)) + FG
          
          estimates[index,] <- tibble_row(thetas[i], sds[j], seeds[k], names(FGs)[FG],
                                          thetaFG, FG_CI[1], FG_CI[2],
                                          aic, bic)
          
          print(paste0(round((index)*100/nrow(estimates),3), '% Done.'))
        }
      }
    }
  }
  return (estimates)
}


#' Model selection for the different FGs using proposed model selection procedure (b): Estimate theta first and then substitute that value for all other models
#'
#' @param data Dataframe created by simulate_response function
#' @param FGs named vector or list giving the functional groupings to test
#'
#' @return A dataframe with the estimated theta values and other statistics for each unique combination of theta, sd and seed values in the data
select_FG <- function(data, FGs){
  if(missing(data)){
    stop('Please specify a dataframe with species proportions')
  }
  
  species <- data %>% select(starts_with('p', ignore.case = F) & !contains('_add')) %>% colnames()
  #print(species)
  if (sum(grep('_add',colnames(data))) != 0){
    data <- data %>% select(-contains('_add'))
    
  }
  
  if(missing(FGs)){
    stop('FG argument is mandatory for functional group model and takes character strings with functional group names referring to each species, in order of apperance in data')
  }
  # if(!is.vector(FG, mode='character') | length(FG)!=length(species)){
  #   stop('FG argument is a vector with same length as number of species and takes character strings with functional group names referring to each species, in order of appearance in data')
  # }
  
  thetas <-  unique(data$theta)    #Uncomment this to check for all thetas, but it will take some time to run
  sds <- unique(data$sd)
  seeds <- unique(data$seed)
  
  
  ### Creating with null values dataframe first as it would be faster to just modify the rows
  ### rather than adding rows to an empty dataframe
  
  len <- length(seeds)*length(thetas)*length(sds)*length(FGs)   
  estimates <- tibble(true_theta = vector(mode = 'numeric', length=len),
                      sd = vector(mode='numeric', length=len),
                      seed = vector(mode = 'numeric', length=len),
                      FG = vector(mode='character', length = len),
                      theta_est = vector(mode='numeric', length=len),
                      AIC = vector(mode='numeric', length=len),
                      BIC = vector(mode='numeric', length=len))
  ## Main execution loop
  for (i in 1:length(thetas)){
    upper_boundary <- ifelse(thetas[i]>1.1, 2.5, 1.5)
    for (j in 1:length(sds)){
      for (k in 1:length(seeds)){
        data_iter <- data %>% filter(theta==thetas[i] & sd==sds[j] & seed==seeds[k])  ## Taking a specific subset of data
        
        group = FGs[[1]]
        modFG <- DI(y='Response', prop=species, FG=group, DImodel = 'FG', data=data_iter)
        cust_modFG <- DI_theta(obj=modFG, DImodel = 'FG',FGnames=group , prop=species, nSpecies= length(species), family='gaussian', lower_boundary = 0.01, upper_boundary = upper_boundary)
        # Test if theta differs from 1 significantly
        test <- anova(modFG, cust_modFG, test = 'F')
        include_theta <- (test$`Pr(>F)`[2] <= 0.05)
        thetaAV <- ifelse(include_theta, cust_modFG$coefficients['theta'], 1)
        
        aicAV <- ifelse(include_theta, cust_modFG$aic, modFG$aic)
        bicAV <- ifelse(include_theta, BIC(cust_modFG), BIC(modFG))
        
        thetaFG <- cust_modFG$coefficients['theta']
        aic <- cust_modFG$aic
        bic <- BIC(cust_modFG)
        
        index <- (length(FGs)*length(seeds)*length(sds)*(i-1) + length(FGs)*length(seeds)*(j-1) + length(FGs)*(k-1)) + 1
        
        estimates[index,] <- tibble_row(thetas[i], sds[j], seeds[k], names(FGs)[1],
                                        thetaFG,
                                        aic, bic)
        
        for (FG in 2:length(FGs)){
          group = FGs[[FG]]
          modFG <- DI(y='Response', prop=species, FG=group, DImodel = 'FG', data=data_iter, theta = thetaFG)
          #thetaFG <- modFG$coefficients['theta']
          aic <- modFG$aic
          bic <- BIC(modFG)
          index <- (length(FGs)*length(seeds)*length(sds)*(i-1) + length(FGs)*length(seeds)*(j-1) + length(FGs)*(k-1)) + FG
          estimates[index,] <- tibble_row(thetas[i], sds[j], seeds[k], names(FGs)[FG],
                                          thetaFG,
                                          aic, bic)
          
          print(paste0(round((index)*100/nrow(estimates),3), '% Done.'))
        }
        
      }
    }
  }
  return (estimates)
}

#' Model selection for the selecting the best interaction structure using proposed model selection procedure (b): Estimate theta first and then substitute that value for all other models
#'
#' @param data Dataframe created by simulate_response function
#' @param FG Vector giving the functional groupings of species
#'
#' @return A dataframe with the estimated theta values and other statistics for each unique combination of theta, sd and seed values in the data
model_selection_b <- function(data, FG){
  if(missing(data)){
    stop('Please specify a dataframe with species proportions')
  }
  
  species <- data %>% select(starts_with('p', ignore.case = F) & !contains('_add')) %>% colnames()
  #print(species)
  if (sum(grep('_add',colnames(data))) != 0){
    data <- data %>% select(-contains('_add'))
    
  }
  
  if(missing(FG)){
    stop('FG argument is mandatory for functional group model and takes character strings with functional group names referring to each species, in order of apperance in data')
  }
  if(!is.vector(FG, mode='character') | length(FG)!=length(species)){
    stop('FG argument is a vector with same length as number of species and takes character strings with functional group names referring to each species, in order of appearance in data')
  }
  
  if('Treatment' %in% colnames(data)){
    print('Adding treatment to model')
    treat = 'Treatment'
  } else {
    treat = NA
  }
  
  thetas <-  unique(data$theta)    #Uncomment this to check for all thetas, but it will take some time to run
  sds <- unique(data$sd)
  seeds <- unique(data$seed)
  models <- c('Average','Functional','Additive','Full Pairwise')  ## Checking for only these models as 
  ## they're the only ones having a contribution of theta
  
  ### Creating with null values dataframe first as it would be faster to just modify the rows
  ### rather than adding rows to an empty dataframe
  
  len <- length(seeds)*length(thetas)*length(sds)*length(models)   
  estimates <- tibble(true_theta = vector(mode = 'numeric', length=len),
                      sd = vector(mode='numeric', length=len),
                      seed = vector(mode = 'numeric', length=len),
                      model = vector(mode='character', length = len),
                      theta_est = vector(mode='numeric', length=len),
                      aic= vector(mode='numeric', length=len),
                      bic= vector(mode='numeric', length=len))
  
  ## Main execution loop
  for (i in 1:length(thetas)){
    upper_boundary <- ifelse(thetas[i]>1.1, 2.5, 1.5)
    data_kter <- data %>% filter(theta==thetas[i])
    for (j in 1:length(sds)){
      data_jter <- data_kter %>% filter(sd==sds[j])
      for (k in 1:length(seeds)){
        data_iter <- data_jter %>% filter(seed==seeds[k])  ## Taking a specific subset of data
        
        modAV <- DI(y='Response', prop=species, treat = treat, DImodel = 'AV', data=data_iter)
        cust_modAV <- DI_theta(obj=modAV, DImodel = 'AV', prop=species, nSpecies= length(species), family='gaussian', lower_boundary = 0.01, upper_boundary = upper_boundary)
        
        # Test if theta differs from 1 significantly
        test <- anova(modAV, cust_modAV, test = 'F')
        include_theta <- (test$`Pr(>F)`[2] <= 0.05)
        thetaAV <- ifelse(include_theta, cust_modAV$coefficients['theta'], 1)
        
        aicAV <- ifelse(include_theta, cust_modAV$aic, modAV$aic)
        bicAV <- ifelse(include_theta, BIC(cust_modAV), BIC(modAV))
        
        modFG <- DI(y='Response', prop=species, treat = treat, FG=FG, DImodel = 'FG', data=data_iter, theta = thetaAV)
        thetaFG <- modFG$coefficients['theta']
        aicFG <- modFG$aic
        bicFG <- BIC(modFG)
        
        modADD <- DI(y='Response', prop=species, treat = treat, DImodel = 'ADD', data=data_iter, theta = thetaAV)
        thetaADD <- modADD$coefficients['theta']
        aicADD <- modADD$aic
        bicADD <- BIC(modADD)
        
        modFULL <- DI(y='Response', prop=species, treat = treat, DImodel = 'FULL', data=data_iter, theta = thetaAV)
        thetaFULL <- modFULL$coefficients['theta']
        aicFULL <- modFULL$aic
        bicFULL <- BIC(modFULL)
        
        index <- (length(models)*length(seeds)*length(sds)*(i-1) + length(models)*length(seeds)*(j-1) + length(models)*(k-1))
        
        estimates[index+1,] <- tibble_row(thetas[i], sds[j], seeds[k], models[1],
                                          thetaAV, aicAV, bicAV)
        estimates[index+2,] <- tibble_row(thetas[i], sds[j], seeds[k], models[2],
                                          thetaAV, aicFG, bicFG)
        estimates[index+3,] <- tibble_row(thetas[i], sds[j], seeds[k], models[3],
                                          thetaAV, aicADD, bicADD)
        estimates[index+4,] <- tibble_row(thetas[i], sds[j], seeds[k], models[4],
                                          thetaAV, aicFULL, bicFULL)
        
        print(paste0(round((index + 4)*100/nrow(estimates),3), '% Done.'))
        
      }
    }
  }
  return (estimates)
}

#' Model selection for the selecting the best interaction structure using model selection procedure (a): Select interaciton structure first and then estimate theta
#'
#' @param data Dataframe created by simulate_response function
#' @param FG Vector giving the functional groupings of species
#'
#' @return A dataframe with the estimated theta values and other statistics for each unique combination of theta, sd and seed values in the data
model_selection_a <- function(data, FG){
  if(missing(data)){
    stop('Please specify a dataframe with species proportions')
  }
  
  species <- data %>% select(starts_with('p', ignore.case = F) & !contains('_add')) %>% colnames()
  #print(species)
  if (sum(grep('_add',colnames(data))) != 0){
    data <- data %>% select(-contains('_add'))
    
  }
  
  if(missing(FG)){
    stop('FG argument is mandatory for functional group model and takes character strings with functional group names referring to each species, in order of apperance in data')
  }
  if(!is.vector(FG, mode='character') | length(FG)!=length(species)){
    stop('FG argument is a vector with same length as number of species and takes character strings with functional group names referring to each species, in order of appearance in data')
  }
  
  if('Treatment' %in% colnames(data)){
    print('Adding treatment to model')
    treat = Treatment
  } else {
    treat = NA
  }
  
  thetas <-  unique(data$theta)    #Uncomment this to check for all thetas, but it will take some time to run
  sds <- unique(data$sd)
  seeds <- unique(data$seed)
  models <- c('Average','Functional','Additive','Full Pairwise')  ## Checking for only these models as 
  ## they're the only ones having a contribution of theta
  
  ### Creating with null values dataframe first as it would be faster to just modify the rows
  ### rather than adding rows to an empty dataframe
  
  len <- length(seeds)*length(thetas)*length(sds)*length(models)   
  estimates <- tibble(true_theta = vector(mode = 'numeric', length=len),
                      sd = vector(mode='numeric', length=len),
                      seed = vector(mode = 'numeric', length=len),
                      model = vector(mode='character', length = len),
                      aic= vector(mode='numeric', length=len),
                      bic= vector(mode='numeric', length=len))
  
  ## Main execution loop
  for (i in 1:length(thetas)){
    upper_boundary <- ifelse(thetas[i]>1.1, 2.5, 1.5)
    data_kter <- data %>% filter(theta==thetas[i])
    for (j in 1:length(sds)){
      data_jter <- data_kter %>% filter(sd==sds[j])
      for (k in 1:length(seeds)){
        data_iter <- data_jter %>% filter(seed==seeds[k])  ## Taking a specific subset of data
        
        modAV <- DI(y='Response', prop=species,  DImodel = 'AV', data=data_iter)
        aicAV <- modAV$aic
        bicAV <- BIC(modAV)
        
        modFG <- DI(y='Response', prop=species, FG=FG, DImodel = 'FG', data=data_iter)
        aicFG <- modFG$aic
        bicFG <- BIC(modFG)
        
        modADD <- DI(y='Response', prop=species, DImodel = 'ADD', data=data_iter)
        aicADD <- modADD$aic
        bicADD <- BIC(modADD)
        
        modFULL <- DI(y='Response', prop=species,DImodel = 'FULL', data=data_iter)
        aicFULL <- modFULL$aic
        bicFULL <- BIC(modFULL)
        
        index <- (length(models)*length(seeds)*length(sds)*(i-1) + length(models)*length(seeds)*(j-1) + length(models)*(k-1))
        
        estimates[index+1,] <- tibble_row(thetas[i], sds[j], seeds[k], models[1],
                                          aicAV, bicAV)
        estimates[index+2,] <- tibble_row(thetas[i], sds[j], seeds[k], models[2],
                                          aicFG, bicFG)
        estimates[index+3,] <- tibble_row(thetas[i], sds[j], seeds[k], models[3],
                                          aicADD, bicADD)
        estimates[index+4,] <- tibble_row(thetas[i], sds[j], seeds[k], models[4],
                                          aicFULL, bicFULL)
        print(paste0(round((index + 4)*100/nrow(estimates),3), '% Done.'))
        
      }
    }
  }
  return (estimates)
}

#' Model selection using model selection procedure 3
#'
#' @param theta Dataframe consisting of interaction structures fitting with theta estimated
#' @param notheta Dataframe consisting of interaction structures fitting with theta not estimated
#'
#' @return A dataframe with statistics for each unique combination of theta, sd and seed values in the data
model_selection_c <- function(theta, notheta, metric = 'aic'){
  theta %>% select(true_theta, sd, seed, model, aic, bic) %>% 
    mutate('Version' = 'Theta') %>% 
    bind_rows(notheta %>% mutate('Version' = 'No theta')) %>% 
    arrange(true_theta, sd, seed, model) %>% 
    group_by(true_theta, sd, seed, model) %>% 
    slice_min(get(metric))
}




