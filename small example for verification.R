# This is an example which runs on a small subset of the simulation study
# If you don't want to run the entire simulation study, running this
# example should give a flavour of the entire study and would help in 
# verifying the results as they were similar across the board

# The entire code should take anywhere between 20 and 45 minutes to execute
# depending on the specifications of the machine

install.packages('rstudioapi')

dir <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(dir,'/Helpers/Final_Pipeline.R'))
source(paste0(dir,'/Helpers/Pipeline_Helper_Code.R'))

# Set your desired working directory
setwd(paste0(dir, '/Datasets/4 species'))

# Get species communities
props4 <- read_csv('props4.csv')

#Define coefficients

# Identity effects
coeff_iden <- c(5,7,6,3) #rnorm(9, mean=8, sd =2)

# FULL interaction effects coefficients
set.seed(652)
coeff_interactionsFULL <- round(rnorm(6, mean = 8, sd = 4), 2)

# Simulating responses for the full pairwise interaction structure 
FULL_file <- simulate_response(data=props4,         #Dataset containing proportion of species
                               theta=  c(0.05, 0.35, 0.48, 0.77, 1.00),  #List of thetas between 0 and 1.5
                               sds=c(1),     #List of standard deviations
                               seed=100,         #Specify initial seed
                               num_datasets = 200, #Specify number of different datasets to create per combination of theta and sd
                               model='FULL',       #Choose model to estimate
                               coeff=c(coeff_iden, coeff_interactionsFULL), # Specify coefficients for the model
                               fileName = '4_species_tiny_example')        # Filename for data file 

#Reading the csv file created
responsesFULL <- read_csv(FULL_file)

#Estimating theta. Function returns summary table containing estimates of theta
estimatesFULL <- estimate_theta(data=responsesFULL,     #File containing the simulated response
                                FG=c('G','G','H','H')) #Same FG as specified in simulate_response

# Model selection results
# Method a: Select interaction structure first and then estimate theta
FULL_model_selection_a <- model_selection_a(data=responsesFULL,     #File containing the simulated response
                                            FG=c('G','G','H','H')) #Same FG as specified in simulate_response

# Method b: Estimate theta first and then select best interaction structure
FULL_model_selection_b <- model_selection_b(data=responsesFULL,     #File containing the simulated response
                                            FG=c('G','G','H','H')) #Same FG as specified in simulate_response

# Method c: Estimate theta for all interaction structures and then select best structure
# Can use same data used to assessing robustness of theta as we estimate theta for each interaction structure there
FULL_model_selection_c <- model_selection_c(estimatesFULL, FULL_model_selection_a, metric =  'aic')


##### Visualising results
# Visualizing results for only one true interaction structure. Same results for all
thetas <- unique(estimatesFULL$true_theta)

#Visualizing theta distribution
for (i in 1:length(thetas)){
  print(estimatesFULL %>% filter(true_theta==thetas[i]) %>% 
          ggplot(data=.)+
          geom_histogram(aes(x=theta_est, y=after_stat(density)), fill='light green', bins = 30)+
          geom_density(aes(x=theta_est), colour='blue', size=1)+
          geom_vline(aes(xintercept=true_theta), linetype='dashed', colour='red', size=1)+
          facet_wrap(~model)+
          theme_pubr())
}

# Summary stats
summ <- estimatesFULL %>% group_by(true_theta, model) %>% 
  summarize(mean_est=mean(theta_est),
            median_est=median(theta_est),
            sd_est=sd(theta_est),
            Q2.5 = quantile(theta_est, 0.025),
            Q97.5 = quantile(theta_est, 0.975))
View(summ)


# Analysing theta Coverage
coverage <- estimatesFULL %>% filter(theta_upper!=0) %>%
  mutate('In'=ifelse(true_theta<theta_upper & true_theta>theta_lower, 1, 0))

coverage_per <- coverage %>% filter(sd ==1) %>% group_by(true_theta, model) %>% 
  summarize('Coverage'=mean(In))

coverage_per




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
estimatesFULL %>% filter( model == 'Average') %>% 
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

a <- model_selection_analysis(FULL_model_selection_a, 'FULL', metric = 'aic')
b <- model_selection_analysis(FULL_model_selection_b, 'FULL', metric = 'aic')
c <- model_selection_analysis(FULL_model_selection_c, 'FULL', metric = 'aic')

ggarrange(a,b,c, nrow = 3)
