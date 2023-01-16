dir <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(dir,'/Helpers/Final_Pipeline.R'))
source(paste0(dir,'/Helpers/Pipeline_Helper_Code.R'))

# Set your desired working directory
setwd(paste0(dir, '/Datasets/16 species'))


########## This code chunk from line 14 to 125 can take upto a day to run
########## If it is not desirable to run this code for that long
########## Skip to line 132 to load the pre-computed objects and verify results
################################################################################

# Get species communities
species <- paste0('p',1:16)

# Create experiment design for 16 species example
dat <- create_communities(no_species=16, reps=1, treatments= NA, 
                               random=FALSE, ncomms=NA) %>% 
  filter(Richness %in% c(1,2,4,8,16)) 
  
# Sample communities specific number of communities from each level of richness
props16 <- dat %>%
  group_by(Richness) %>% 
  nest() %>%            
  ungroup() %>% 
  mutate(n = c(16, 80, 40, 20, 10)) %>% 
  mutate(samp = map2(data, n, sample_n, replace = T)) %>% 
  select(-data) %>%
  unnest(samp)

# 2 repititions for each community
props16 <- rbind(props16, props16)

#Define coefficients

# Identity effects
set.seed(573)
coeff_iden <- round(rnorm(16, mean=9, sd =4), 2)

# Interaction Effect
set.seed(764)
coeff_interactionsFG <- c(round(rnorm(6, mean = 11, sd = 2), 2),
                          round(rnorm(4, mean = 5, sd = 2), 2))

# Generating datasets
FG_file <- simulate_response(data=props16,         #Dataset containing proportion of species
                             theta=  c(0.05, 0.35, 0.77, 1.00, 1.33),    #List of thetas between 0 and 1.5
                             sds=c(0.8, 1, 1.2),     #List of standard deviations
                             seed=9000,           #Specify initial seed
                             num_datasets = 1, #Specify number of different datasets to create per combination of theta and sd
                             model='FG',         #Choose model to estimate
                             FG=c(rep('FG1', 4),rep('FG2', 4),
                                  rep('FG3', 4),rep('FG4', 4)),
                             coeff=c(coeff_iden, coeff_interactionsFG), # Specify coefficients for the model
                             fileName = '16_species_fg_groups')      # Filename for data file


#Reading the csv file created
responsesFG <- read_csv(FG_file)

#Estimating theta. Function returns summary table containing estimates of theta
estimatesFG <- test_FG(data=responsesFG,     #File containing the simulated response
                       FGs=list('FG4' = c(rep('FG1', 4),rep('FG2', 4), rep('FG3', 4),rep('FG4', 4)),
                                'FG2' = c(rep('FG1', 8),rep('FG2', 8)),
                                'FG1' = c(rep('FG1', 16)),
                                'FG8' =  c(rep('FG1', 2),rep('FG2', 2), rep('FG3', 2),rep('FG4', 2),
                                           rep('FG5', 2),rep('FG6', 2), rep('FG7', 2),rep('FG8', 2))
                                ))
write_csv(estimatesFG, 'Estimates (16_species_FG).csv')
#View(estimatesFG)

# Test robustness of theta for different functional groupings for this data
modelsFG <- select_FG(data=responsesFG,     #File containing the simulated response
                       FGs=list('FG4' = c(rep('FG1', 4),rep('FG2', 4), rep('FG3', 4),rep('FG4', 4)),
                                'FG2' = c(rep('FG1', 8),rep('FG2', 8)),
                                'FG1' = c(rep('FG1', 16)),
                                'FG8' =  c(rep('FG1', 2),rep('FG2', 2), rep('FG3', 2),rep('FG4', 2),
                                           rep('FG5', 2),rep('FG6', 2), rep('FG7', 2),rep('FG8', 2))
                       ))

write_csv(modelsFG, 'Model Selection (16_species_FG).csv')


#####################################################################################

#### Loading pre-computed results

# Theta robustness results
load('Robustness.rda')

# Model selection results
load('ModelSelection.rda')

thetas <- unique(estimatesFG$true_theta)
#Visualizing theta distribution
for (i in 1:length(thetas)){
  plot(
    estimatesFG %>% filter(sd == 1, true_theta==thetas[i]) %>% 
    ggplot(data=.)+
    geom_histogram(aes(x=theta_est, y=..density..), fill='light green')+
    geom_density(aes(x=theta_est), colour='blue', size=1)+
    geom_vline(aes(xintercept=true_theta), linetype='dashed', colour='red', size=1)+
    facet_wrap(~FG)+
    theme_pubr()
  )
}

# Summary stats
summ <- estimatesFG %>% group_by(true_theta, FG) %>% 
  summarize(sd_est=sd(theta_est),
            mean_est=mean(theta_est),
            median_est=median(theta_est))
View(summ)

# Estimated vs True Theta
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
  facet_wrap(~FG)

# Analysing theta Coverage
coverage <- estimatesFG %>% filter(theta_upper!=0) %>%
  mutate('In'=ifelse(true_theta<theta_upper & true_theta>theta_lower, 1, 0))

coverage_per <- coverage %>% group_by(true_theta, FG) %>% 
  summarize('Coverage'=mean(In))
View(coverage_per)

# Scaled Sd figure
estimatesFG %>% group_by(true_theta, FG) %>%
  summarise('Scaled' = sd(theta_est)/IQR(theta_est)) %>% 
  ggplot(data=.)+
  geom_col(aes(x=true_theta, y=Scaled, fill=FG), 
           position = position_dodge2(0.5), alpha = .5) +
  theme_pubr()+
  #ggtitle('Predicted vs True Theta')+
  xlab('Theta')+
  ylab('Scaled Standard Deviation')+
  scale_fill_manual(labels = c("FG1", 
                               'FG2',
                               'FG4',
                               'FG8'),
                    values = c('#2D16BC','#BC1652',
                               '#139164','#A5BC16'))+
  scale_x_continuous(breaks=thetas,labels = thetas)+
  labs(fill='Model')+
  theme(axis.title = element_text(face='bold', size='12'),
        plot.title = element_text(hjust = 0.5, face='bold', size='14'),
        axis.text = element_text(face='bold', size='10'),
        plot.subtitle = element_text(face='bold', size='11', hjust=0.5 ),
        legend.text = element_text(face='bold', size='10'),
        legend.title = element_text(face='bold', size='12'))

# Model selection figures
true_model <- 'FG4'
c <- modelsFG %>% group_by(true_theta, sd, seed) %>% slice_min(aic) %>% 
  group_by(true_theta, sd, FG) %>% 
  summarise('Prop' = n()/length(unique(.$seed))) %>% 
  arrange(desc(Prop), .by_group= T) %>% 
  mutate('Position' = rev(cumsum(rev(Prop)) - c(Prop[length(Prop)]/2, diff(cumsum(rev(Prop)))/2)),
         'Selected Model' = factor(FG,
                                   levels = c("FG1", 
                                              'FG2',
                                              'FG4',
                                              'FG8'))) %>%
  mutate('label' = ifelse(Prop %in% c(Prop[1],Prop[2]), paste0(Prop*100, '%'), ''),
         'sd' = paste0('sigma==', sd),
         'true_theta' = paste0('theta==', true_theta)) %>% 
  ggplot()+
  geom_col(aes(x = 1, y = Prop, fill = `Selected Model`), alpha = .5)+
  geom_text(aes(x = 1, y = Prop, label = label,
                group = `Selected Model`),
            colour = 'black', position = position_stack(vjust = 0.5))+
  coord_flip()+
  facet_grid(true_theta ~ sd, switch = 'y', labeller = label_parsed)+
  labs(subtitle = paste0('True Model = ', true_model))+
  ylab('Proportion')+
  ggpubr::theme_pubr()+
  theme(plot.subtitle = element_text(face = 'bold', size =16, hjust =0.5),
        axis.title.y = element_blank(),
        strip.text = element_text(face = 'bold', size =12),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),#element_text(label = 'Proportion', face = 'bold', size = 11),
        axis.text.x = element_blank(),#element_text(face = 'bold', size =9),
        legend.text = element_text(face = 'bold', size =12),
        legend.title = element_text(face = 'bold', size =14))+
  # scale_fill_manual(values = c('#F2AFAD','#D9E4AA',
  #                              '#D5B2D4','#F3D1B0'))+
  scale_fill_manual(labels = c("FG1", 
                               'FG2',
                               'FG4',
                               'FG8'),
                    values = c('#2D16BC','#BC1652',
                               '#139164','#A5BC16'))
c
