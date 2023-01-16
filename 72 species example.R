dir <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(dir,'/Helpers/Final_Pipeline.R'))
source(paste0(dir,'/Helpers/Pipeline_Helper_Code.R'))

# Set your desired working directory
setwd(paste0(dir, '/Datasets/72 species'))


########## This code chunk from line 13 to 56 can take upto a day to run
########## If it is not desirable to run this code for that long
########## Skip to line 58 to load the pre-computed objects and verify results
################################################################################
# Get species communities
data('Bell')
species <- paste0('p',1:72)
props72 <- Bell %>% select(id, community, 'Richness' = richness, species)

#Define coefficients

# Identity effects
set.seed(687)
coeff_iden <- round(rnorm(72, mean=8, sd =2),2)

# Interaction Effect
set.seed(259)
coeff_interactionsAV <- round(rnorm(1, mean = 11, sd = 2), 2)

# Generating datasets
AV_file <- simulate_response(data=props72,         #Dataset containing proportion of species
                             theta=  c(0.05, 0.35, 0.77, 1.00, 1.33),    #List of thetas between 0 and 1.5
                             sds=c(0.8, 1, 1.2),     #List of standard deviations
                             seed=7000,           #Specify initial seed
                             num_datasets = 200, #Specify number of different datasets to create per combination of theta and sd
                             model='AV',         #Choose model to estimate
                             coeff=c(coeff_iden, coeff_interactionsAV), # Specify coefficients for the model
                             fileName = '72_species_av_pairwise')      # Filename for data file

#Reading the csv file created
#AV_file <- '72_species_av_pairwise.csv'
responsesAV <- read_csv(AV_file)

#Estimating theta. Function returns summary table containing estimates of theta
estimatesAV <- estimate_theta_high(data=responsesAV,     #File containing the simulated response
                                   FG=c(rep('FG1', 12),rep('FG2', 12),rep('FG3', 12),
                                        rep('FG4', 12),rep('FG5', 12),rep('FG6', 12)))

#Incase you want to write estimates to csv
write_csv(estimatesAV, 'Estimates (72_species_AV).csv')

# Model selection
modelsAV <- model_selection_b(data=responsesAV,     #File containing the simulated response
                               FG=c(rep('FG1', 12),rep('FG2', 12),rep('FG3', 12),
                                    rep('FG4', 12),rep('FG5', 12),rep('FG6', 12)))
write_csv(modelsAV, 'Model Selection (72_species_AV).csv')

#####################################################################################

#### Loading pre-computed results

# Theta robustness results
load('Robustness.rda')

# Model selection results
load('ModelSelection.rda')

thetas <- unique(estimatesAV$true_theta)

#Visualizing theta distribution
for (i in 1:length(thetas)){
  print(
    estimatesAV %>% filter(sd == 1, true_theta==thetas[i]) %>% 
      ggplot(data=.)+
      geom_histogram(aes(x=theta_est, y=..density..), fill='light green', bins = 30)+
      geom_density(aes(x=theta_est), colour='blue', size=1)+
      geom_vline(aes(xintercept=true_theta), linetype='dashed', colour='red', size=1)+
      facet_wrap(~model, nrow = 2)+
      theme_pubr()
  )
}



# Summary stats
summ <- estimatesAV %>% group_by(true_theta, model) %>% 
  summarize(sd_est=sd(theta_est),
            mean_est=mean(theta_est),
            median_est=median(theta_est),
            Q2.5 = quantile(theta_est, 0.025),
            Q97.5 = quantile(theta_est, 0.975))

summ %>% 
  ggplot(.)+
  geom_point(aes(y=mean_est, x=true_theta), size = 1.3)+
  geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5, x=true_theta), size =0.7, alpha=1)+
  geom_abline(aes(slope=1, intercept=0))+
  theme_pubr()+
  ggtitle('Predicted vs True Theta')+
  xlab('True Theta')+
  ylab('Mean Theta Estimate')+
  scale_x_continuous(breaks=thetas,labels = thetas)+
  labs(color='Model')+
  theme(axis.title = element_text(face='bold', size='12'),
        plot.title = element_text(hjust = 0.5, face='bold', size='14'),
        axis.text = element_text(face='bold', size='10'),
        plot.subtitle = element_text(face='bold', size='11', hjust=0.5 ),
        legend.text = element_text(face='bold', size='10'),
        legend.title = element_text(face='bold', size='12'))+
  facet_wrap(~model, nrow = 2)

# Analysing theta Coverage
coverage <- estimatesAV %>% #filter(theta_upper!=0) %>%
  mutate('In'=ifelse(true_theta<theta_upper & true_theta>theta_lower, 1, 0))

coverage_per <- coverage %>% group_by(true_theta, model) %>% 
  summarize('Coverage'=mean(In))
View(coverage_per)


# Scaled SD
estimatesAV %>% group_by(true_theta, model) %>%
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
estimatesAV %>% filter(sd %in% c(.8,1,1.2), true_theta %in% c(0.05, 0.35, 0.77, 1, 1.33), model == 'Average') %>% 
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

# Model selection figure for appendix C
# Model selection
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

model_selection_analysis(modelsAV, 'Average Pairwise', metric = 'aic')
