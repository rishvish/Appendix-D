dir <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(dir,'/Helpers/Final_Pipeline.R'))
source(paste0(dir,'/Helpers/Pipeline_Helper_Code.R'))

# Set your desired working directory
setwd(paste0(dir, '/Datasets/9 species'))


########## This code chunk from line 14 to 167 can take upto a day to run
########## If it is not desirable to run this code for that long
########## Skip to line 169 to load the pre-computed objects and verify results
################################################################################

library(DImodels)
data(sim3)

# Get species communities
props9 <- sim3[, paste0('p',1:9)] %>% 
            distinct() %>% 
            slice(rep(1:100, times = 3))

#Define coefficients

# Identity effects
set.seed(388)
coeff_idenFULL <- rnorm(9, mean=7, sd =2)

# FULL interaction effects coefficients
set.seed(389)
coeff_interactionsFULL <- round(rnorm(36, mean = 9, sd = 6), 2)

# Identity effects
set.seed(793)
coeff_idenADD <- rnorm(9, mean=7, sd =2)

# ADD interaction effects coefficients
set.seed(827)
coeff_interactionsADD <- round(rnorm(9, mean = 9, sd = 3), 2)

# Identity effects
set.seed(722)
coeff_idenFG <- rnorm(9, mean=7, sd =2)

# FG interaction effects coefficients
set.seed(146)
coeff_interactionsFG <- round(c(rnorm(3, mean = 10, sd = 2), 
                                rnorm(3, mean = 5, sd = 3)), 2)

# Identity effects
set.seed(637)
coeff_idenAV <- rnorm(9, mean=7, sd =2)

# AV interaction effects coefficients
set.seed(990)
coeff_interactionsAV <- round(rnorm(1, mean = 10, sd = 3), 2)


# Simulating responses for each interaction structure 
FULL_file <- simulate_response(data=props9,         #Dataset containing proportion of species
                               theta=  c(0.05, 0.19, 0.35, 0.48, 0.63, 0.77, 0.91, 1.00, 1.17, 1.33),  #List of thetas between 0 and 1.5
                               sds=c(0.8, 0.9, 1, 1.1, 1.2),     #List of standard deviations
                               seed=1000,         #Specify initial seed
                               num_datasets = 200, #Specify number of different datasets to create per combination of theta and sd
                               model='FULL',       #Choose model to estimate
                               coeff=c(coeff_idenFULL, coeff_interactionsFULL), # Specify coefficients for the model
                               fileName = '9_species_full_pairwise')        # Filename for data file 

ADD_file <- simulate_response(data=props9,         #Dataset containing proportion of species
                              theta=  c(0.05, 0.19, 0.35, 0.48, 0.63, 0.77, 0.91, 1.00, 1.17, 1.33), #List of thetas between 0 and 1.5
                              sds=c(0.8, 0.9, 1, 1.1, 1.2),     #List of standard deviations
                              seed=1200,         #Specify initial seed
                              num_datasets = 200, #Specify number of different datasets to create per combination of theta and sd
                              model='ADD',       #Choose model to estimate
                              coeff=c(coeff_idenADD, coeff_interactionsADD), # Coefficiens for model
                              fileName = '9_species_add_species') # Filename for data file

FG_file <- simulate_response(data=props9,         #Dataset containing proportion of species
                             theta=  c(0.05, 0.19, 0.35, 0.48, 0.63, 0.77, 0.91, 1.00, 1.17, 1.33),  #List of thetas between 0 and 1.5
                             sds=c(0.8, 0.9, 1, 1.1, 1.2),     #List of standard deviations
                             seed=1400,            #Specify initial seed
                             num_datasets = 200,  #Specify number of different datasets to create per combination of theta and sd
                             model='FG',          #Choose model to estimate
                             FG=c('G','G','G','G','G','H','H','L','L'), #Specify Functional groups. (Needed only if model is FG)
                             coeff=c(coeff_idenFG, coeff_interactionsFG), #Specify coefficients for the model
                             fileName = '9_species_fg_group')      #Filename for data file

AV_file <- simulate_response(data=props9,         #Dataset containing proportion of species
                             theta=  c(0.05, 0.19, 0.35, 0.48, 0.63, 0.77, 0.91, 1.00, 1.17, 1.33),    #List of thetas between 0 and 1.5
                             sds=c(0.8, 0.9, 1, 1.1, 1.2),     #List of standard deviations
                             seed=1600,           #Specify initial seed
                             num_datasets = 200, #Specify number of different datasets to create per combination of theta and sd
                             model='AV',         #Choose model to estimate
                             coeff=c(coeff_idenAV, coeff_interactionsAV), # Specify coefficients for the model
                             fileName = '9_species_av_pairwise')      # Filename for data file

#file <- 'SimulatedDataFULL (16 species high reps).csv'
#Reading the csv file created
responsesFULL <- read_csv(FULL_file)
responsesAV <- read_csv(AV_file)
responsesFG <- read_csv(FG_file)
responsesADD <- read_csv(ADD_file)

# Estimating theta for each true interaction structure

FG <- c('G','G','G','G','G','H','H','L','L')
#Estimating theta. Function returns summary table containing estimates of theta
estimatesFULL <- estimate_theta(data=responsesFULL,     #File containing the simulated response
                                FG=FG) #Same FG as specified in simulate_response
estimatesAV <- estimate_theta(data=responsesAV,     #File containing the simulated response
                              FG=FG)
estimatesFG <- estimate_theta(data=responsesFG,     #File containing the simulated response
                              FG=FG)
estimatesADD <- estimate_theta(data=responsesADD,     #File containing the simulated response
                               FG=FG)
View(estimatesFULL)

#Incase you want to write estimates to csv
#write_csv(estimatesFULL, 'Estimates (9_species_FULL).csv')
#write_csv(estimatesFULL, 'Estimates (9_species_FG).csv')
#write_csv(estimatesFULL, 'Estimates (9_species_ADD).csv')
#write_csv(estimatesFULL, 'Estimates (9_species_AV).csv')


# Model selection results
# Method a: Select interaction structure first and then estimate theta
FULL_model_selection_a <- model_selection_a(data=responsesFULL,     #File containing the simulated response
                                            FG=FG) #Same FG as specified in simulate_response
AV_model_selection_a <- model_selection_a(data=responsesAV,     #File containing the simulated response
                                          FG=FG)
FG_model_selection_a <- model_selection_a(data=responsesFG,     #File containing the simulated response
                                          FG=FG)
ADD_model_selection_a <- model_selection_a(data=responsesADD,     #File containing the simulated response
                                           FG=FG)
# Method b: Estimate theta first and then select best interaction structure
FULL_model_selection_b <- model_selection_b(data=responsesFULL,     #File containing the simulated response
                                            FG=FG) #Same FG as specified in simulate_response
AV_model_selection_b <- model_selection_b(data=responsesAV,     #File containing the simulated response
                                          FG=FG)
FG_model_selection_b <- model_selection_b(data=responsesFG,     #File containing the simulated response
                                          FG=FG)
ADD_model_selection_b <- model_selection_b(data=responsesADD,     #File containing the simulated response
                                           FG=FG)
# Method c: Estimate theta for all interaction structures and then select best structure
# Can use same data used to assessing robustness of theta as we estimate theta for each interaction structure there
FULL_model_selection_c <- model_selection_c(estimatesFULL, FULL_model_selection_a, metric =  'aic')
AV_model_selection_c <- model_selection_c(estimatesAV, AV_model_selection_a, metric =  'aic')
FG_model_selection_c <- model_selection_c(estimatesFG, FG_model_selection_a, metric =  'aic')
ADD_model_selection_c <- model_selection_c(estimatesADD, ADD_model_selection_a, metric =  'aic')


# write data if you want to save
#write_csv(FULL_model_selection_a, 'Model Selection (9_species_FULL_a).csv')
#write_csv(AV_model_selection_a, 'Model Selection (9_species_FG_a).csv')
#write_csv(FG_model_selection_a, 'Model Selection (9_species_ADD_a).csv')
#write_csv(ADD_model_selection_a, 'Model Selection (9_species_AV_a).csv')

#write_csv(FULL_model_selection_b, 'Model Selection (9_species_FULL_b).csv')
#write_csv(AV_model_selection_b, 'Model Selection (9_species_FG_b).csv')
#write_csv(FG_model_selection_b, 'Model Selection (9_species_ADD_b).csv')
#write_csv(ADD_model_selection_b, 'Model Selection (9_species_AV_b).csv')

#write_csv(FULL_model_selection_c, 'Model Selection (9_species_FULL_c).csv')
#write_csv(AV_model_selection_c, 'Model Selection (9_species_FG_c).csv')
#write_csv(FG_model_selection_c, 'Model Selection (9_species_ADD_c).csv')
#write_csv(ADD_model_selection_c, 'Model Selection (9_species_AV_c).csv')

#####################################################################################

#### Loading pre-computed results

# Theta robustness results
load('Robustness.rda')

# Model selection results
load('ModelSelection.RData')

##### Visualising results
# Visualizing results for only one true interaction structure. Same results for all
thetas <- unique(estimatesFULL$true_theta)

#Visualizing theta distribution
for (i in 1:length(thetas)){
  print(
    estimatesFULL %>% filter(sd == 1, true_theta==thetas[i]) %>% 
    ggplot(data=.)+
    geom_histogram(aes(x=theta_est, y=..density..), fill='light green', bins = 30)+
    geom_density(aes(x=theta_est), colour='blue', size=1)+
    geom_vline(aes(xintercept=true_theta), linetype='dashed', colour='red', size=1)+
    facet_wrap(~model)+
    theme_pubr()
  )
}



# Summary stats
summ <- estimatesFULL %>% group_by(true_theta, model) %>% 
  summarize(sd_est=sd(theta_est),
            mean_est=mean(theta_est),
            median_est=median(theta_est),
            Q2.5 = quantile(theta_est, 0.025),
            Q97.5 = quantile(theta_est, 0.975))
#View(summ)


# Analysing theta Coverage
coverage <- estimatesFULL %>% filter(theta_upper!=0) %>%
  mutate('In'=ifelse(true_theta<theta_upper & true_theta>theta_lower, 1, 0))

coverage_per <- coverage %>% filter(sd ==1) %>% group_by(true_theta, model) %>% 
  summarize('Coverage'=mean(In))
View(coverage_per)


## True vs estimated theta
summ <- estimatesFULL %>% group_by(true_theta, model) %>% 
  summarize(sd_est=sd(theta_est),
            mean_est=mean(theta_est),
            median_est=median(theta_est),
            Q2.5 = quantile(theta_est, 0.025),
            Q97.5 = quantile(theta_est, 0.975)) %>% 
  mutate(model = fct_relevel(model, 'Average', 'Functional',
                             'Additive', 'Full Pairwise')) 
#View(summ)

levels(summ$model) <- c('Average Pairwise', 'Functional Groups',
                        'Additive Species', 'Full Pairwise')

summ %>%
  ggplot(.)+
  geom_point(aes(y=mean_est, x=true_theta), size = 1.3)+
  geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5, x=true_theta), size =0.7, alpha=1)+
  geom_abline(aes(slope=1, intercept=0))+
  theme_pubr()+
  #ggtitle('Predicted vs True Theta')+
  labs(x = expression(bold(paste('True ', theta))),
       y = expression(bold(paste('Mean ', theta, ' Estimate'))))+
  scale_x_continuous(breaks=thetas,labels = thetas)+
  labs(color='Model')+
  theme(axis.title = element_text(face='bold', size='12'),
        plot.title = element_text(hjust = 0.5, face='bold', size='14'),
        axis.text = element_text(face='bold', size='10'),
        plot.subtitle = element_text(face='bold', size='11', hjust=0.5 ),
        legend.text = element_text(face='bold', size='10'),
        legend.title = element_text(face='bold', size='12'))+
  facet_wrap(~model)


# Scaled SD
estimatesFULL %>% group_by(true_theta, model) %>%
  summarise('Scaled' = sd(theta_est)/IQR(theta_est)) %>% 
  ggplot(data=.)+
  geom_col(aes(x=true_theta, y=Scaled, fill=model), 
           position = position_dodge2(0.5), alpha = .7) +
  theme_pubr()+
  #ggtitle('Predicted vs True Theta')+
  xlab('Theta')+
  ylab('Scaled Standard Deviation')+
  scale_x_continuous(breaks=thetas,labels = thetas)+
  labs(fill='Model')+
  theme(axis.title = element_text(face='bold', size='11'),
        plot.title = element_text(hjust = 0.5, face='bold', size='14'),
        axis.text = element_text(face='bold', size='10'),
        plot.subtitle = element_text(face='bold', size='11', hjust=0.5 ),
        legend.text = element_text(face='bold', size='10'),
        legend.title = element_text(face='bold', size='11'))+
  scale_fill_manual(labels = c("Average Pairwise", 
                               'Functional Groups',
                               'Additive Species',
                               'Full Pairwise'),
                    values = c('#F15A60','#7AC36A',
                               '#FAA75B','#9E67AB'))
# Effect of sigma on theta
estimatesFULL %>% filter(sd %in% c(.8,1,1.2), true_theta %in% c(0.05, 0.35, 0.77, 1, 1.33), model == 'Average') %>% 
  mutate('true_theta1' = factor(true_theta, labels = paste0('theta==', unique(true_theta))),
         'sd' = factor(sd, labels= paste0('sigma==', unique(sd)))) %>% 
  ggplot()+
  geom_histogram(aes(x=theta_est, y=..density..), fill='light green', bins = 30)+
  geom_density(aes(x=theta_est), colour='blue', size=1)+
  geom_vline(aes(xintercept=true_theta), linetype='dashed', colour='red', size=1)+
  facet_grid(true_theta1 ~ sd, scales = 'free', labeller = label_parsed)+
  ggpubr::theme_pubr()+
  xlab('Theta Estimate')+
  ylab('Density')+
  theme(axis.title = element_text(face='bold', size='12'),
        axis.text = element_text(face='bold', size='11'),
        plot.subtitle = element_text(face='bold', size='12', hjust=0.5 ),
        strip.text = element_text(face='bold', size='12'))

#' Plot model selection results
#'
#' @param data Data from model seleciton functions
#' @param true_model True underlying model
#' @param metric metric to use for selecting best model. One of aic or bic
#'
#' @return A ggplot object
model_selection_analysis <- function(data, true_model, metric = 'aic'){
  data %>% group_by(true_theta, sd, seed) %>% 
    slice_min(get(metric)) %>% 
    ungroup() %>% 
    mutate('Selected Model' = model) %>% 
    group_by(true_theta, sd, `Selected Model`) %>% 
    summarise('Prop' = n()/length(unique(.$seed))) %>% 
    arrange(desc(Prop), .by_group= T) %>% 
    mutate('Position' = rev(cumsum(rev(Prop)) - c(Prop[length(Prop)]/2, diff(cumsum(rev(Prop)))/2)),
           'Selected Model' = factor(`Selected Model`,
                                     levels = c('Average', 'Functional',
                                                'Additive', 'Full Pairwise', 'ID'))) %>%
    mutate('label' = ifelse(Prop %in% c(Prop[1], Prop[2]), paste0(Prop*100, '%'), ''),
           'sd' = paste0('sigma==', sd),
           'true_theta' = paste0('theta==', true_theta)) %>% 
    ggplot()+
    geom_col(aes(x = 1, y = Prop, fill = `Selected Model`), alpha = .7)+
    geom_text(aes(x = 1, y = Prop, label = label,
                  group = `Selected Model`),
              colour = 'black', position = position_stack(vjust = 0.5))+
    coord_flip()+
    facet_grid(true_theta ~ sd, switch = 'y', labeller = label_parsed)+
    labs(subtitle = paste0('True Model = ', true_model))+
    ylab('Proportion')+
    ggpubr::theme_pubr()+
    theme(plot.subtitle = element_text(face = 'bold', size =11, hjust =0.5),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_text(face = 'bold', size =11),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),#element_text(label = 'Proportion', face = 'bold', size = 11),
          axis.text.x = element_blank(),#element_text(face = 'bold', size =9),
          legend.text = element_text(face = 'bold', size =10),
          legend.title = element_text(face = 'bold', size =11))+
    scale_fill_manual(labels = as_labeller(
      c("Average" = "Average Pairwise", 
        "Functional" = 'Functional Groups',
        "Additive" = 'Additive Species',
        "Full Pairwise" = 'Full Pairwise',
        "ID" = 'No interaction')),
      values = c('#F15A60','#7AC36A',
                 '#FAA75B','#9E67AB',
                 'black'), drop = F)
  
}

a <- model_selection_analysis(FULL_model_selection_a, '9 species: True Model = Full Pairwise', metric = 'aic')
b <- model_selection_analysis(FULL_model_selection_b, '9 species: True Model = Full Pairwise', metric = 'aic')
c <- model_selection_analysis(FULL_model_selection_c, '9 species: True Model = Full Pairwise', metric = 'aic')

ggarrange(a,b,c, nrow = 3)
