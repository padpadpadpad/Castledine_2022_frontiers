# Analysis of Castledine et al. 2022: Greater phage genotypic diversity constrains arms-race coevolution

# load packages ####
library(lme4)
library(emmeans)
library(tidyverse)
library(flextable)

# load in data time shift assay
d <- read.csv("data/timeseries_assay.csv", header = TRUE)

#-----------------------------#
# wrangle and process data ####
#-----------------------------#

# calculate mean infectivity for each population and phage and bacterial time point combination
d_means <- group_by(d, treatment, phage_tp, bacteria_tp, phage_tp_cat, bacteria_tp_cat, replicate_culture) %>%
  summarise(mean_infectivity = mean(infective),
            mean_resistance = 1 - mean_infectivity,
            .groups = 'drop') %>%
  mutate(., time_shift = bacteria_tp - phage_tp)

# get data in correct format for binomial glm
# calculate number of infective and total number of clones tested per combination (always 12)
d_sum <- group_by(d, treatment, bacteria_tp, phage_tp, phage_tp_cat, bacteria_tp_cat, replicate_culture, time_shift) %>%
  summarise(., infective = sum(infective),
            not_infective = 12 - infective,
            .groups = 'drop') %>%
  mutate(., obs = 1:n(),
         tot_clones = 12)

#-----------------------------------------------------#
# analyses on contemporary resistance through time ####
#-----------------------------------------------------#

#subset data for contemporary populations
d_contemporary <- filter(d_sum, bacteria_tp == phage_tp)

# run binomial model - look for differences in infectivity through time between treatments
cont_mod <- glmer(cbind(infective, not_infective) ~ bacteria_tp_cat * treatment + (1|replicate_culture), 
                  d_contemporary, 
                  family = binomial('logit'), 
                  control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)), 
                  na.action = na.fail)

blmeco::dispersion_glmer(cont_mod) #not overdispersed (should not be over 1.4)

# do model simplification
cont_mod2 <- glmer(cbind(infective, not_infective) ~ bacteria_tp_cat + treatment + (1|replicate_culture), 
                   d_contemporary, 
                   family = binomial('logit'), 
                   control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)), 
                   na.action = na.fail)

anova(cont_mod, cont_mod2) #significant interaction

# look at pairwise differences
emmeans::emmeans(cont_mod, pairwise ~ bacteria_tp_cat|treatment, type = "response")

#-----------------------------------------------#
# analyses on bacteria and phage coevolution ####
#-----------------------------------------------#

# analysing whether bacteria and phage coevolve in each treatment

# model 1
mod1_glm_overdisp <- glmer(cbind(infective, not_infective) ~ bacteria_tp_cat * treatment * phage_tp_cat + (1|replicate_culture), 
                           d_sum, 
                           family = binomial)

# model 2 with observation level random effect to account for overdispersion
mod1_glm <- glmer(cbind(infective, not_infective) ~ bacteria_tp_cat * treatment * phage_tp_cat + (1 |replicate_culture) + (1|obs), 
                  d_sum, family = binomial, 
                  control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))

# check for overdispersion (it should not be over 1.4)
blmeco::dispersion_glmer(mod1_glm_overdisp) # 1.57
blmeco::dispersion_glmer(mod1_glm) # 0.68 - much better

# do model selection
# remove 3-way interaction
mod2_glm <- glmer(cbind(infective, not_infective) ~ bacteria_tp_cat + treatment + phage_tp_cat + phage_tp_cat:treatment + bacteria_tp_cat:treatment + bacteria_tp_cat:phage_tp_cat + (1|replicate_culture) + (1|obs), 
                  d_sum, family = binomial, 
                  control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
anova(mod2_glm) # remove bacteria_tp:phage_tp as has lowest F value

# remove bacteria_tp:phage_tp
mod3_glm <- glmer(cbind(infective, not_infective) ~ bacteria_tp_cat + treatment + phage_tp_cat + phage_tp_cat:treatment + bacteria_tp_cat:treatment + (1|replicate_culture) + (1|obs), 
                  d_sum, 
                  family = binomial, 
                  control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
anova(mod3_glm) # remove treatment:phage_tp as has lowest F value

mod4_glm <- glmer(cbind(infective, not_infective) ~ bacteria_tp_cat + treatment + phage_tp_cat + bacteria_tp_cat:treatment + (1|replicate_culture) + (1|obs), 
                  d_sum,
                  family = binomial,
                  control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
anova(mod4_glm) # remove bacteria_tp:treatment as last interaction

mod5_glm <- glmer(cbind(infective, not_infective) ~ bacteria_tp_cat + treatment + phage_tp_cat + (1|replicate_culture) + (1|obs), 
                  d_sum,
                  family = binomial,
                  control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
anova(mod5_glm)

mod6_glm <- glmer(cbind(infective, not_infective) ~ bacteria_tp_cat + treatment + bacteria_tp_cat:treatment + (1|replicate_culture) + (1|obs), 
                  d_sum, 
                  family = binomial, 
                  control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))

#model simplification
anova(mod1_glm, mod2_glm)
anova(mod2_glm, mod3_glm)
anova(mod3_glm, mod4_glm)
anova(mod4_glm, mod5_glm) #sig bact:treat interaction
anova(mod4_glm, mod6_glm) #sig phage time

# model mod4_glm is the best

#post-hoc Tukey-HSDs
emmeans::emmeans(mod4_glm, pairwise ~ phage_tp_cat, type = "response") #phage increase infectivity through time
emmeans::emmeans(mod4_glm, pairwise ~ bacteria_tp_cat|treatment, type = "response") # bact generally become more resistant through time with phage_1 and ish for phage_2, not both_phage

#-----------------#
# Make Figures ####
#-----------------#

# create correct panels for treatment
d_means <- mutate(d_means, treatment2 = case_when(treatment == 'both_phage' ~ '(c) Phage 1 and 2',
                                                  treatment == 'phage_1' ~ '(a) Phage 1',
                                                  treatment == 'phage_2' ~ '(b) Phage 2'))

# labels for x axis
labs <- c("Day 4", "Day 8", "Day 12")

# calculate averages for all data
d_averages <- group_by(d_means, treatment, treatment2, phage_tp, bacteria_tp, bacteria_tp_cat, phage_tp_cat) %>%
  summarise(infect = mean(mean_infectivity),
            resist = mean(mean_resistance),
            se = sd(mean_infectivity)/sqrt(length(mean_infectivity)),
            .groups = 'drop')

# Figure 1 
d_contemporary_averages <- filter(d_averages, phage_tp == bacteria_tp)

ggplot(d_contemporary_means, aes(x = bacteria_tp_cat, y = infect)) + #bacteria tp and phage tp now the same
  theme_bw(base_size = 14) +
  geom_point(data = filter(d_means, bacteria_tp == phage_tp), aes(x = bacteria_tp, y = mean_infectivity), size = 2, position = position_jitter(0.1, height = 0), alpha = 0.5) +
  facet_wrap(~treatment2) +
  geom_errorbar(aes(x = bacteria_tp, ymin = infect-se, ymax = infect+se), width = 0.2) +
  geom_point(size = 3.5, shape = 21, col = "black", fill = "white") +
  xlab("Time point") +
  ylab("Phage infectivity")  +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 14, hjust = 0),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "black")) +
  scale_x_discrete(labels = labs)

ggsave('plots/Figure_1.pdf', last_plot(), width = 10, height = 4)
ggsave('plots/Figure_1.png', last_plot(), width = 10, height = 4)

# Figure 2
ggplot(d_averages, aes(x = bacteria_tp_cat, y = infect, group = phage_tp_cat)) +
  theme_bw() +
  geom_point(data = d_means, aes(x = bacteria_tp_cat, y = mean_infectivity, group = phage_tp_cat, col = phage_tp_cat), size = 1, position = position_jitterdodge(0.8), alpha = 0.5) +
  facet_wrap(~treatment2) +
  geom_errorbar(aes(x = bacteria_tp_cat, ymin = infect-se, ymax = infect+se), width = 0.2, position = position_dodge(0.8)) +
  geom_point(aes(fill = phage_tp_cat), size = 2, position = position_dodge(0.8), shape = 21, col = "black") +
  xlab("Bacterial time point") +
  ylab("Phage infectivity")  +
  palettetown::scale_color_poke(pokemon = 'charizard', spread = 3, labels = labs) +
  palettetown::scale_fill_poke(pokemon = 'charizard', spread = 3, labels = labs) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14, colour = "black"), strip.background = element_blank(), strip.text = element_text(size = 14, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  guides(fill = guide_legend("Phage time point"), col = guide_legend("Phage time point")) +
  scale_x_discrete(labels = labs) 

ggsave('plots/Figure_2.pdf', last_plot(), width = 10, height = 6)
ggsave('plots/Figure_2.png', last_plot(), width = 10, height = 6)

#----------------#
# Make Tables ####
#----------------#

# Table S1

# grab contrasts
contrast_contemporary <- emmeans::emmeans(cont_mod, pairwise ~ bacteria_tp_cat|treatment, type = "response")$contrasts %>%
  data.frame() %>%
  mutate(contrast = gsub('b_tp_', '', contrast),
         across(where(is.factor), as.character),
         treatment2 = case_when(treatment == 'both_phage' ~ 'Phage 1 and 2',
                                                treatment == 'phage_1' ~ 'Phage 1',
                                                treatment == 'phage_2' ~ 'Phage 2'),
         order = parse_number(treatment)) %>%
  arrange(order) %>%
  select(treatment2, contrast, odds.ratio, SE, z.ratio, p.value) %>%
  mutate(across(odds.ratio:p.value, signif, digits = 3),
         across(odds.ratio:p.value, as.character))
  
names(contrast_contemporary) <- c("Treatment", "Time contrast", "Odds-ratio", "SE", "z-ratio", "p-value")

contrast_table <- flextable(contrast_contemporary) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit() %>%
  hline(i = c(3, 6, 9), border = officer::fp_border(color="black")) %>%
  align(j = c(1:6), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')  %>%
  bold(i = c(1, 7:9))

contrast_table

# save
save_as_image(contrast_table, 'tables/Table_S1.png', webshot = 'webshot2')

# Table S2 - infectivity estimates for all bacteria-phage time point combinations across all days

d_averages_table <- mutate(d_averages, phage_tp2 = case_when(phage_tp == 1 ~ 'Day 4',
                                                       phage_tp == 2 ~ 'Day 8',
                                                       phage_tp == 3 ~ 'Day 12'),
                     bacteria_tp2 = case_when(bacteria_tp == 1 ~ 'Day 4',
                                              bacteria_tp == 2 ~ 'Day 8',
                                              bacteria_tp == 3 ~ 'Day 12'),
                     treatment2 = case_when(treatment == 'both_phage' ~ 'Phage 1 and 2',
                                            treatment == 'phage_1' ~ 'Phage 1',
                                            treatment == 'phage_2' ~ 'Phage 2'),
                     order = parse_number(treatment)) %>%
  arrange(order) %>%
  select(treatment2, bacteria_tp2, phage_tp2, infect, resist, se) %>%
  mutate(across(infect:se, signif, digits = 3),
         across(infect:se, as.character))

names(d_averages_table) <- c("Treatment", "Bacteria Tp", "Phage Tp", "Susceptible", "Resistant", "SE")

means_table <- flextable(d_averages_table) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit() %>%
  hline(i = c(9, 18, 27), border = officer::fp_border(color="black")) %>%
  align(j = c(1:6), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')

means_table

# save
save_as_image(contrast_table, 'tables/Table_S2.png', webshot = 'webshot2')

# Table S3 - bacterial susceptibility over time across all phage time points
resist_tab <- emmeans::emmeans(mod4_glm, pairwise ~ bacteria_tp_cat|treatment, type = "response")$contrasts %>%
  data.frame() %>%
  mutate(contrast = gsub('b_tp_', '', contrast),
         across(where(is.factor), as.character),
         treatment2 = case_when(treatment == 'both_phage' ~ 'Phage 1 and 2',
                                treatment == 'phage_1' ~ 'Phage 1',
                                treatment == 'phage_2' ~ 'Phage 2'),
         order = parse_number(treatment)) %>%
  arrange(order) %>%
  select(treatment2, contrast, odds.ratio, SE, z.ratio, p.value) %>%
  mutate(across(odds.ratio:p.value, signif, digits = 3),
         across(odds.ratio:p.value, as.character))

names(contrast_contemporary) <- c("Treatment", "Time contrast", "Odds-ratio", "SE", "z-ratio", "p-value")

resist_flex <- flextable(resist_tab) %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit() %>%
  hline(i = c(3, 6, 9), border = officer::fp_border(color="black")) %>%
  align(j = c(1:6), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header') %>%
  bold(i = c(1, 2, 5))

resist_flex

# save
save_as_image(contrast_table, 'tables/Table_S3.png', webshot = 'webshot2')
