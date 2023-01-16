dir <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(dir,'/Helpers/Final_Pipeline.R'))
source(paste0(dir,'/Helpers/Pipeline_Helper_Code.R'))


# Set your desired working directory
setwd(paste0(dir, '/Datasets/6 species'))

########## This code chunk from line 13 to 60 can take upto a day to run
########## If it is not desirable to run this code for that long
########## Skip to line 64 to load the pre-computed objects and verify results
################################################################################
# Get species communities
species <- paste0('p',1:6)
design <- read_csv('expdesign.csv')
props6 <- bind_rows(design %>% mutate('Treatment' = 'A'),
                    design %>% mutate('Treatment' = 'B'))
  

#Define coefficients

# Identity effects
set.seed(141)
coeff_iden <- round(rnorm(6, mean=7, sd =2), 2)

# Interaction Effect
set.seed(12)
coeff_interactionsADD <- c(round(rnorm(6, mean = 9, sd = 3), 2))

# treatment effects
treats <- c(5, 8)

# Generating datasets
ADD_file <- simulate_response(data=props6,         #Dataset containing proportion of species
                             theta=  c(0.05, 0.35, 0.77, 1.00, 1.33),    #List of thetas between 0 and 1.5
                             sds=c(0.8, 1, 1.2),     #List of standard deviations
                             seed=11000,           #Specify initial seed
                             num_datasets = 200, #Specify number of different datasets to create per combination of theta and sd
                             model='ADD',         #Choose model to estimate
                             coeff=c(coeff_iden, coeff_interactionsADD, treats), # Specify coefficients for the model
                             fileName = '6_species_add_treatment')      # Filename for data file

#Reading the csv file created
responsesADD <- read_csv(ADD_file)

FG = c('G','G','H','H','L','L')
#Estimating theta. Function returns summary table containing estimates of theta
estimatesADD <- estimate_theta(data=responsesADD,     #File containing the simulated response
                               FG=FG)
#Incase you want to write estimates to csv
write_csv(estimatesADD, 'Estimates (6_species_add_treatment).csv')
View(estimatesADD)
save(estimatesADD, file = 'Robustness.rda')
# Model selection
modelsADD <- model_selection_b(data=responsesADD,     #File containing the simulated response
                               FG=FG)
write_csv(modelsADD, 'Model Selection (6_species_ADD_treat).csv')

#####################################################################################

#### Loading pre-computed results

# Theta robustness results
load('Robustness.rda')

# Model selection results
load('ModelSelection.rda')

thetas <- unique(estimatesADD$true_theta)

for (i in 1:length(thetas)){
  print(
    estimatesADD %>% filter(sd == 1, true_theta==thetas[i]) %>% 
    ggplot(data=.)+
    geom_histogram(aes(x=theta_est, y=..density..), fill='light green', bins = 30)+
    geom_density(aes(x=theta_est), colour='blue', size=1)+
    geom_vline(aes(xintercept=true_theta), linetype='dashed', colour='red', size=1)+
    facet_wrap(~model)+
    theme_pubr()
  )
}

# Summary stats
summ <- estimatesADD %>% group_by(true_theta, model) %>% 
  summarize(sd_est=sd(theta_est),
            mean_est=mean(theta_est),
            median_est=median(theta_est))
View(summ)



summ %>% mutate(max=mean_est+1.96*sd_est, min=mean_est-1.96*sd_est) %>%
  filter(true_theta<1.5) %>% 
  ggplot(.)+
  geom_point(aes(y=mean_est, x=true_theta), size = 1.5, position = position_dodge(width=0.1))+
  geom_errorbar(aes(ymin=min, ymax=max, x=true_theta), size =0.75, alpha=1, width=0.1, position = 'dodge')+
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
  facet_wrap(~model)

# Analysing theta Coverage
coverage <- estimatesADD %>% filter(theta_upper!=0) %>%
  mutate('In'=ifelse(true_theta<theta_upper & true_theta>theta_lower, 1, 0))

coverage_per <- coverage %>% group_by(true_theta, model) %>% 
  summarize('Coverage'=mean(In))
View(coverage_per)

estimatesADD %>% group_by(true_theta, model) %>%
  mutate(sd=sd(theta_est), lower=theta_est-1.96*sd, upper=theta_est+1.96*sd) %>% 
  mutate('In'=ifelse(true_theta<upper & true_theta>lower, 1, 0)) %>% 
  group_by(true_theta, model) %>% 
  summarize('Coverage'=mean(In)) 


# Predicted vs True Theta
estimatesADD %>% group_by(true_theta, model) %>%
  summarise('Scaled' = sd(theta_est)/IQR(theta_est)) %>% 
  ggplot(data=.)+
  geom_col(aes(x=true_theta, y=Scaled, fill=model), 
           position = position_dodge2(0.5)) +
  theme_pubr()+
  #ggtitle('Predicted vs True Theta')+
  xlab('Theta')+
  ylab('Scaled Standard Deviation')+
  scale_x_continuous(breaks=thetas,labels = thetas)+
  labs(fill='Model')+
  theme(axis.title = element_text(face='bold', size='12'),
        plot.title = element_text(hjust = 0.5, face='bold', size='14'),
        axis.text = element_text(face='bold', size='10'),
        plot.subtitle = element_text(face='bold', size='11', hjust=0.5 ),
        legend.text = element_text(face='bold', size='10'),
        legend.title = element_text(face='bold', size='12'))

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

model_selection_analysis(modelsADD, 'ADD', metric = 'aic')

# Compare with method c
model_selection_analysis(estimatesADD, 'ADD', metric = 'aic')
