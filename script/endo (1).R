
# Assessment of relationship between AMF diversity and tree divers --------

library(tidyverse)  # version 2.0.0
library(vegan)      # version 2.6.4
library(agricolae)  # version 1.3.7
library(corrplot)   # version 0.92
library(grDevices)  # version 4.3.2


# Loading Data ------------------------------------------------------------

endo <- read_csv('data/endo.csv')

# Checking first few rows of the data -------------------------------------

head(endo)

# Dropping the plot collumn -----------------------------------------------

endo_div <- endo %>% select (-plot_no)


# Running Shannon Diversity Index -----------------------------------------

endo_shannon <- diversity(endo_div, 'shannon')#simpson and inv also available


# Adding endo_shannon to endo datasheet -----------------------------------

mutated_endo <- endo %>% mutate(diversity = endo_shannon)

# Visualizing the data ----------------------------------------------------

png(filename = 'output/fragments_endo.png', width = 1800, height = 2200, res = 300)
mutated_endo %>% select(plot_no, diversity) %>% 
  mutate(fragment = case_when(startsWith(as.character(plot_no),'NG') ~ 'Ngangao',
                              startsWith(as.character(plot_no),'CH') ~ 'Chawia',
                              startsWith(as.character(plot_no),'FUR') ~ 'Fururu')) %>% 
  ggplot(aes(reorder(fragment, -diversity), diversity)) + 
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_classic() +
  labs(x = 'Forest fragments',
       y = 'Shannon diversity index') +
  theme(axis.text.x = element_text(colour = 'black', face = 'bold')) +
  theme(axis.text.y = element_text(colour = 'black', face = 'bold')) +
  theme(text = element_text(size = 12))
dev.off()  


# Testing significance difference (ANOVA) ---------------------------------

endo_aov <- mutated_endo %>% select(plot_no, diversity) %>% 
  mutate(fragment = case_when(startsWith(as.character(plot_no),'NG') ~ 'Ngangao',
                              startsWith(as.character(plot_no),'CH') ~ 'Chawia',
                              startsWith(as.character(plot_no),'FUR') ~ 'Fururu'))

model_aov <- aov(formula = diversity ~ fragment, data = endo_aov)
summary(aov(formula = diversity ~ fragment, data = endo_aov))

# Mean comparison of significance difference ------------------------------

TukeyHSD(model_aov)
b <- LSD.test(model_aov, trt = 'fragment')
b$groups 
  
#######################################################
  # TREE DIVERSITY
#######################################################

# Loading tree data -------------------------------------------------------

tree <- read_csv('data/trees.csv')
head(tree)

# Dropping columns --------------------------------------------------------

tree_div <- tree %>% select(-`Plot no`, -`AMF comp`, -`Altitudes(m)`, -Disturbances)

# Running Shannon diversity Index -----------------------------------------

tree_shannon <- diversity(tree_div, 'shannon')
tree_shannon


# Mutating Shannon to tree data -------------------------------------------

mutated_tree <- tree %>% mutate(tree_diversity = tree_shannon)


# Visualization of data ---------------------------------------------------

png(filename = 'output/fragments_tree.png', width = 1800, height = 2200, res = 300)
mutated_tree %>% select(`Plot no`, tree_diversity) %>% 
  mutate(fragment = case_when(startsWith(as.character(`Plot no`),'NG') ~ 'Ngangao',
                              startsWith(as.character(`Plot no`),'CH') ~ 'Chawia',
                              startsWith(as.character(`Plot no`),'FUR') ~ 'Fururu')) %>% 
  ggplot(aes(reorder(fragment, -tree_diversity), tree_diversity)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  theme_classic() +
  labs(x = 'Forest fragments',
       y = 'Shannon diversity index') +
  theme(axis.text.x = element_text(colour = 'black', face = 'bold')) +
  theme(axis.text.y = element_text(colour = 'black', face = 'bold')) +
  theme(text = element_text(size = 12))
dev.off()

# Testing for significance difference -------------------------------------

tree_aov <- mutated_tree %>% select(`Plot no`, tree_diversity) %>% 
  mutate(fragment = case_when(startsWith(as.character(`Plot no`),'NG') ~ 'Ngangao',
                              startsWith(as.character(`Plot no`),'CH') ~ 'Chawia',
                              startsWith(as.character(`Plot no`),'FUR') ~ 'Fururu'))

model_aov <- aov(formula = tree_diversity ~ fragment, data = tree_aov)
summary(aov(formula = tree_diversity ~ fragment, data = tree_aov))

# Mean comparison of significance difference ------------------------------

TukeyHSD(model_aov)
b <- LSD.test(model_aov, trt = 'fragment')
b$groups 

tree_aov %>% 
  mutate(tree_aov = (formula = tree_diversity)) %>% 
  #group_by(fragments) %>% 
  summarise(sd = sd(tree_aov))


# Correlation  ------------------------------------------------------------

cor_data <- data.frame(endo_diversity = mutated_endo$diversity,
                       tree_diversity = mutated_tree$tree_diversity,
                       altitude = tree$`Altitudes(m)`,
                       disturbance = tree$Disturbances,
                       fragment = tree_aov$fragment)

cor_data_Ngangao <- cor_data %>% filter(fragment == 'Ngangao')
cor_data_Chawia <- cor_data %>% filter(fragment == 'Chawia')
cor_data_Fururu <- cor_data %>% filter(fragment == 'Fururu')
cor_data_Ngangao <- cor_data_Ngangao %>% 
  rename(AMF_diversity = endo_diversity)
cor_Ngangao <- cor(cor_data_Ngangao[1:3 ])
corrplot::corrplot(cor_Ngangao, method = 'ellipse', diag = FALSE, 
         addCoef.col = 'black', number.cex = 0.8, type = 'lower',
         tl.cex = 1.3, tl.col = 'black')

cor_Chawia <- cor(cor_data_Chawia[1:3 ])
corrplot(cor_Chawia, method = 'ellipse', diag = FALSE, 
         addCoef.col = 'black', number.cex = 0.9, type = 'lower',
         tl.cex = 1.3, tl.col = 'black')

cor_Fururu <- cor(cor_data_Fururu[1:3 ])
corrplot(cor_Fururu, method = 'ellipse', diag = FALSE, 
         addCoef.col = 'black', number.cex = 0.9, type = 'lower', 
         tl.col = 'black', tl.cex = 1.3) 

# Working with disturbance column -----------------------------------------

cor_mat_Ngangao <- cor_data_Ngangao %>% 
  mutate(disturbance_level = 
  case_when(disturbance == 'none' ~ 0,
  disturbance == 'low' ~ 1,
  disturbance == 'moderate' ~ 2)) %>% 
  select(-disturbance, -fragment) %>% cor()

png(filename = 'output/Ngangao_cor.png', width = 1600, height = 2200, res = 300) 
corrplot(cor_mat_Ngangao, method = 'ellipse', diag = FALSE, 
         addCoef.col = 'black', number.cex = 0.9, type = 'lower',
         tl.cex = 1.3, tl.col = 'black',sig.level = 'label_sig')
dev.off()

cor_data_Fururu <- cor_data_Fururu %>% 
  rename(AMF_diversity = endo_diversity)
cor_mat_Fururu <- cor_data_Fururu %>% mutate(disturbance_level = 
                                      case_when(disturbance == 'low' ~ 0,
                                      disturbance == 'moderate' ~ 1)) %>% 
  select(-disturbance, -fragment) %>% cor()
png(filename = 'output/Fururu_cor.png', width = 1600, height = 2200, res = 300)
corrplot(cor_mat_Fururu, method = 'ellipse', diag = FALSE, 
         addCoef.col = 'black', number.cex = 0.9, type = 'lower',
         tl.cex = 1.3, tl.col = 'black',sig.level = 'label_sig')
dev.off()

cor_data_Chawia <- cor_data_Chawia %>% 
  rename(AMF_diversity = endo_diversity)
cor_mat_Chawia <- cor_data_Chawia %>% mutate(disturbance_level = 
                                               case_when(disturbance == 'none' ~ 0,
                                                         disturbance == 'low' ~ 1,
                                                         disturbance == 'moderate' ~ 2,
                                                         disturbance == 'intense' ~ 3)) %>% 
  select(-disturbance, -fragment) %>% cor()


png(filename = 'output/Chawia_cor.png', width = 1600, height = 2200, res = 300)
corrplot(cor_mat_Chawia, method = 'ellipse', diag = FALSE, 
         addCoef.col = 'black', number.cex = 0.9, type = 'lower',
         tl.cex = 1.3, tl.col = 'black',sig.level = 'label_sig')



















