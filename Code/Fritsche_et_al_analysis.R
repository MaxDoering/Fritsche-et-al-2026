
# Fritsche et al. 2026 #
# Microplastic ingestion induces changes in coelomocyte composition of Eisenia fetida #



# 1 Load packages ----

library("tidyr")
library("ggplot2")
library("ggpattern")
library("dplyr")
library("ggpubr")
library("lme4")
library("DHARMa")
library("car") 
library("multcomp") 
library("RColorBrewer")
library("dunn.test")



# 2 Read in data ---- 


in_vivo_data = read.csv2("Data/Fritsche_et_al_in_vivo_data.csv")
sub_g1_data = read.csv2("Data/Fritsche_et_al_subG1_data.csv")
ex_vivo_data = read.csv2("Data/Fritsche_et_al_ex_vivo_data.csv")



# 3 Adjust data type ----


in_vivo_data$day = as.factor(in_vivo_data$day)
in_vivo_data$treatment = as.factor(in_vivo_data$treatment)
in_vivo_data$starting_weight = as.numeric(in_vivo_data$starting_weight)
in_vivo_data$ending_weight = as.numeric(in_vivo_data$ending_weight)
in_vivo_data$number_peaks = as.factor(in_vivo_data$number_peaks)
in_vivo_data$total_number_cells = as.numeric(in_vivo_data$total_number_cells)
in_vivo_data$total_number_cells_per_mg = as.numeric(in_vivo_data$total_number_cells_per_mg)
in_vivo_data$cell_viability = as.numeric(in_vivo_data$cell_viability)
in_vivo_data$total_number_living_cells = as.numeric(in_vivo_data$total_number_living_cells)
in_vivo_data$total_number_living_cells_per_mg = as.numeric(in_vivo_data$total_number_living_cells_per_mg)
in_vivo_data$eleozytes_percent = as.numeric(in_vivo_data$eleozytes_percent)
in_vivo_data$eleozytes_count = as.numeric(in_vivo_data$eleozytes_count)
in_vivo_data$eleozytes_count_per_mg = as.numeric(in_vivo_data$eleozytes_count_per_mg)
in_vivo_data$amoebozytes_percent = as.numeric(in_vivo_data$amoebozytes_percent)
in_vivo_data$amoebozytes_count = as.numeric(in_vivo_data$amoebozytes_count)
in_vivo_data$amoebozytes_count_per_mg = as.numeric(in_vivo_data$amoebozytes_count_per_mg)
in_vivo_data$intermediate_population_count = as.numeric(in_vivo_data$intermediate_population_count)
in_vivo_data$intermediate_population_count_per_mg = as.numeric(in_vivo_data$intermediate_population_count_per_mg)
in_vivo_data$intermediate_population_percent = as.numeric(in_vivo_data$intermediate_population_percent)
in_vivo_data$treatment = factor(in_vivo_data$treatment, levels = c("control", "mock", "ps", "additive"))

# Remove one dead worm (no worm died in the uptake experiment)
in_vivo_data = subset(in_vivo_data, !is.na(number_peaks))


sub_g1_data$treatment = as.factor(sub_g1_data$treatment)
sub_g1_data$cell_population = as.factor(sub_g1_data$cell_population)
sub_g1_data$percentage = as.numeric(sub_g1_data$percentage)


ex_vivo_data$timepoint = as.factor(ex_vivo_data$timepoint)
ex_vivo_data$particle_concentration = as.factor(ex_vivo_data$particle_concentration)
ex_vivo_data$treatment = as.factor(ex_vivo_data$treatment)
ex_vivo_data$absorption = as.numeric(ex_vivo_data$absorption)
ex_vivo_data$relative_metabolic_activity = as.numeric(ex_vivo_data$relative_metabolic_activity)



# 4 Plot data ----

## 4.1 Generate single plots ----


# Total number of isolated cells per mg freshweight (Fig. 3A)
plot_number_cells_per_mg_fw = 
  ggplot(in_vivo_data, aes(x = treatment, y = total_number_cells_per_mg, fill = treatment)) +  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0), size = 1.5, color = 1) +
  guides(fill = "none", color = "none") +
  theme_classic() +
  scale_fill_brewer() +
  labs(x = "", y = "total cells per mg fw") +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 4, color = "gold2") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  scale_x_discrete(labels = c("control", "mock", "PS", expression(PS[irgafos])))

plot_number_cells_per_mg_fw


# Cell viability (Fig. 3B)
plot_cell_viability =
  ggplot(in_vivo_data, aes(x = treatment, y = cell_viability, fill = treatment)) +  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0), size = 1.5, color = 1) +
  guides(fill = "none", color = "none") +
  theme_classic() +
  scale_fill_brewer() +
  labs(x = "", y = "cell viability [%]") +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 4, color = "gold2") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  scale_x_discrete(labels = c("control", "mock", "PS", expression(PS[irgafos]))) +
  scale_y_continuous(breaks = c(50,60,70,80,90,100), labels = c(50,60,70,80,90,100), expand = c(0,0), limits = c(48,102)) +
  geom_text(
    data = in_vivo_data %>% group_by(treatment) %>% summarise(y = max(cell_viability, na.rm = TRUE)) %>%
      cbind(label = c("","a","b","b")),
    aes(x = treatment, y = max(y), label = label),
    position = position_dodge(width = 1),
    vjust = -1,
    size = 7,
    stat = "unique",
    color = "black")

plot_cell_viability


# Total number living cells per mg freshweight (Fig. 3C)
plot_living_cells_per_mg_fw =
  ggplot(in_vivo_data, aes(x = treatment, y = total_number_living_cells_per_mg, fill = treatment)) +  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0), size = 1.5, color = 1) +
  guides(fill = "none", color = "none") +
  theme_classic() +
  scale_fill_brewer() +
  labs(x = "", y = "viable cells per mg fw") +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 4, color = "gold2") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  scale_x_discrete(labels = c("control", "mock", "PS", expression(PS[irgafos]))) +
  scale_y_continuous(breaks = c(0,250,500,750,1000,1250), labels = c(0,250,500,750,1000,1250), expand = c(0,0), limits = c(-50,1450)) + 
  geom_text(
    data = in_vivo_data %>% group_by(treatment) %>% summarise(y = max(total_number_living_cells_per_mg, na.rm = TRUE)) %>%
      cbind(label = c("","a","ab","b")),
    aes(x = treatment, y = max(y), label = label),
    position = position_dodge(width = 1),
    vjust = -1,
    size = 7,
    stat = "unique",
    color = "black")

plot_living_cells_per_mg_fw


# Number of peaks (Fig. 4)
number_peaks = as.data.frame(table(in_vivo_data$number_peaks, in_vivo_data$treatment))
colnames(number_peaks) = c("peak_number", "treatment", "frequency")
number_peaks$peak_number = as.numeric(number_peaks$peak_number)
number_peaks$peak_number = replace(number_peaks$peak_number, number_peaks$peak_number == 1, "unimodal")
number_peaks$peak_number = replace(number_peaks$peak_number, number_peaks$peak_number == 2, "bimodal")
number_peaks$peak_number = replace(number_peaks$peak_number, number_peaks$peak_number == 3, "trimodal")
number_peaks$peak_number = factor(number_peaks$peak_number, levels = c("unimodal", "bimodal", "trimodal")) 

number_peaks = number_peaks %>% 
  group_by(treatment) %>%
  mutate(rel_frequency = (frequency/sum(frequency))*100)

number_peaks$rel_frequency = round(number_peaks$rel_frequency, 1)

plot_number_peaks =
  ggplot(number_peaks, aes(x = treatment, y = rel_frequency, fill = peak_number, label = rel_frequency)) + 
  geom_bar(position = "fill", stat = "identity") +
  theme_classic() + 
  geom_text(aes(label = ifelse(rel_frequency == 0, "", rel_frequency)),
            size = 4, position = position_fill(vjust = 0.5), fontface = "bold")+ 
  scale_fill_manual(values = c("grey90","grey70","grey50")) +
  labs(x = "", y = "cell subpopulations [%]", fill = "") +
  scale_x_discrete(labels = c("control","mock", "PS", expression(PS[irgafos]))) +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14),
        legend.title = element_text(size = 12), legend.text = element_text(size = 12), legend.key.size = unit(0.7, 'cm')) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 25, 50, 75, 100)) 

plot_number_peaks


# Percentage of amoebocytes (Fig. 5B)
plot_amoebocytes_percent =
  ggplot(in_vivo_data, aes(x = treatment, y = amoebozytes_percent, fill = treatment)) +  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0), size = 1.5, color = 1) +
  guides(fill = "none", color = "none") +
  theme_classic() +
  scale_fill_brewer() +
  labs(x = "", y = expression(amoebocytes*" "*i.e.*" "*Riboflavin^Dim*" [%]")) +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 4, color = "gold2") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  scale_x_discrete(labels = c("control", "mock", "PS", expression(PS[irgafos]))) 

plot_amoebocytes_percent


# Percentage of intermediate population (Fig. 5C)
plot_intermediate_population_percentage =
  ggplot(in_vivo_data, aes(x = treatment, y = intermediate_population_percent, fill = treatment)) +  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0), size = 1.5, color = 1) +
  guides(fill = "none", color = "none") +
  theme_classic() +
  scale_fill_brewer() +
  labs(x = "", y = expression(Riboflavin^Inter*" [%]")) +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 4, color = "gold2") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  scale_x_discrete(labels = c("control", "mock", "PS", expression(PS[irgafos]))) +
  scale_y_continuous(breaks = c(seq(0, 50, 10)), labels = c(seq(0, 50, 10)), expand = c(0,0), limits = c(-2,48)) 

plot_intermediate_population_percentage


# Percentage of eleocytes (Fig. 5D)
plot_eleocytes_percent = 
  ggplot(in_vivo_data, aes(x = treatment, y = eleozytes_percent, fill = treatment)) +  
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0), size = 1.5, color = 1) +
  guides(fill = "none", color = "none") +
  theme_classic() +
  scale_fill_brewer() +
  labs(x = "", y = expression(eleocytes*" "*i.e.*" "*Riboflavin^Bright*" [%]")) +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 4, color = "gold2") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  scale_x_discrete(labels = c("control", "mock", "PS", expression(PS[irgafos])))+
  scale_y_continuous(breaks = c(seq(40, 90, 10)), labels = c(seq(40, 90, 10)), expand = c(0,0), limits = c(36,100))

plot_eleocytes_percent


# PS ex vivo  
# This plot is similar to Fig. 6A, which was ultimately done in Origin in the final manuscript
std.error = function(x) sd(na.omit(x))/sqrt(length(na.omit(x)))
ex_vivo_data_mean = data.frame(ex_vivo_data %>%
                                 group_by(treatment,timepoint, particle_concentration) %>%
                                 summarize(mean_value = mean(relative_metabolic_activity, na.rm = TRUE), 
                                           median_value = median(relative_metabolic_activity, na.rm = TRUE), 
                                           sd_value = sd(relative_metabolic_activity, na.rm = TRUE),
                                           se_value = std.error(relative_metabolic_activity)))

ex_vivo_data_mean_PS = subset(ex_vivo_data_mean, treatment == "PS" & particle_concentration != "500")
ex_vivo_data_mean_PS$timepoint_treatment = paste(ex_vivo_data_mean_PS$timepoint, ex_vivo_data_mean_PS$treatment, sep = "_")
ex_vivo_data_mean_PS$timepoint_treatment = factor(ex_vivo_data_mean_PS$timepoint_treatment,
                                                  levels = c("24_PS", "72_PS", "144_PS"))

ex_vivo_data_mean_PS$particle_concentration = as.numeric(as.character(ex_vivo_data_mean_PS$particle_concentration))

ex_vivo_data_mean_PS_plot = 
  ggplot(ex_vivo_data_mean_PS, aes(x = particle_concentration, y = mean_value, colour = timepoint_treatment, shape = timepoint_treatment)) +  
  geom_point(size = 3) +
  theme_classic() + 
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.1, linewidth = 1) +
  scale_y_continuous(breaks = c(seq(0,150,10)), labels = c(seq(0,150,10)), expand = c(0,0), limits = c(0, 151)) + 
  scale_x_log10(breaks = ex_vivo_data_mean_PS$particle_concentration,
                labels = ex_vivo_data_mean_PS$particle_concentration) +
  labs(x = expression("particle concentration ["*mu*"g/ml]"), y = "metabolic activity [%]") +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme(legend.position = c(0.9, 0.2),
        legend.key.size = unit(0.9, units = "cm"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  scale_colour_manual(values = c("steelblue1","steelblue3","steelblue4"),
                      labels = c("24 h", "72 h", "144 h"),
                      name = "") +
  scale_shape_manual(values = c(15,19,17),
                     labels = c("24 h", "72 h", "144 h"),
                     name = "")

ex_vivo_data_mean_PS_plot


# PS irgafos ex vivo
# This plot is similar to Fig. 6B, which was ultimately done in Origin in the final manuscript
ex_vivo_data_mean_PS_Irgafos = subset(ex_vivo_data_mean, treatment == "PS_Irgafos" & particle_concentration != "500")
ex_vivo_data_mean_PS_Irgafos$timepoint_treatment = paste(ex_vivo_data_mean_PS_Irgafos$timepoint, ex_vivo_data_mean_PS_Irgafos$treatment, sep = "_")
ex_vivo_data_mean_PS_Irgafos$timepoint_treatment = factor(ex_vivo_data_mean_PS_Irgafos$timepoint_treatment,
                                                  levels = c("24_PS_Irgafos", "72_PS_Irgafos", "144_PS_Irgafos"))

ex_vivo_data_mean_PS_Irgafos$particle_concentration = as.numeric(as.character(ex_vivo_data_mean_PS_Irgafos$particle_concentration))

ex_vivo_data_mean_PS_Irgafos_plot = 
  ggplot(ex_vivo_data_mean_PS_Irgafos, aes(x = particle_concentration, y = mean_value, colour = timepoint_treatment, shape = timepoint_treatment)) +  
  geom_point(size = 3) +
  theme_classic() + 
  theme(legend.position = "bottom") + 
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), width = 0.03, linewidth = 0.1) +
  scale_y_continuous(breaks = c(seq(50,150,10)), labels = c(seq(50,150,10)), expand = c(0,0), limits = c(50, 151)) + 
  scale_x_log10(breaks = ex_vivo_data_mean_PS$particle_concentration,
                labels = ex_vivo_data_mean_PS$particle_concentration) +
  labs(x = expression("particle concentration ["*mu*"g/ml]"), y = "metabolic activity [%]") +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme(legend.position = c(0.9, 0.95),
        legend.key.size = unit(0.9, units = "cm"),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  scale_colour_manual(values = c("steelblue1","steelblue3","steelblue4"),
                      labels = c("24 h", "72 h", "144 h"),
                      name = "") +
  scale_shape_manual(values = c(15,19,17),
                     labels = c("24 h", "72 h", "144 h"),
                     name = "")

ex_vivo_data_mean_PS_Irgafos_plot


## 4.2 Combine plots ----


# Fig 3
ggarrange(plot_number_cells_per_mg_fw, plot_cell_viability, plot_living_cells_per_mg_fw, 
          nrow = 1, ncol = 3, 
          labels = c("A", "B", "C"), font.label = list(size = 22, color = "black", face = "bold", family = NULL),
          hjust = -0.3)

# Fig 5
ggarrange(plot_amoebocytes_percent, plot_intermediate_population_percentage, plot_eleocytes_percent, 
          nrow = 1, ncol = 3, 
          labels = c("B", "C", "D"), font.label = list(size = 22, color = "black", face = "bold", family = NULL),
          hjust = -0.3)

# Fig 6
ggarrange(ex_vivo_data_mean_PS_plot, ex_vivo_data_mean_PS_Irgafos_plot,
          nrow = 1, ncol = 2, 
          labels = c("A", "B"), font.label = list(size = 22, color = "black", face = "bold", family = NULL),
          hjust = -0.2, common.legend	= F)



# 5 Statistical analysis ----

## 5.1 Effects of MPP on E. fetida and coelomocyte quantity and viability ----


data_without_control = subset(in_vivo_data, !treatment == "control")


# Starting weight between treatments
model = lmer(starting_weight~treatment + (1|day), data = data_without_control)
plot(simulateResiduals(model))
Anova(model)  


# Ending weight between treatments
model = lmer(ending_weight~treatment + (1|day), data = data_without_control)
plot(simulateResiduals(model))
Anova(model) 


# Total number of isolated cells per mg freshweight
model = lmer(total_number_cells_per_mg~treatment + (1|day), data = data_without_control)
plot(simulateResiduals(model))
hist(in_vivo_data$total_number_cells_per_mg)
hist(resid(model))
Anova(model)

data.frame(in_vivo_data %>%
             group_by(treatment) %>%
             summarize(mean_value = mean(total_number_cells_per_mg, na.rm = TRUE), 
                       median_value = median(total_number_cells_per_mg, na.rm = TRUE), 
                       sd_value = sd(total_number_cells_per_mg, na.rm = TRUE)))


# Cell viability
model = lmer(cell_viability~treatment + (1|day), data = data_without_control)
plot(simulateResiduals(model))
hist(in_vivo_data$cell_viability)
hist(resid(model)) 
Anova(model)

PH = glht(model, mcp(treatment = "Tukey"))
summary(PH, test = adjusted("BH"))

data.frame(in_vivo_data %>%
             group_by(treatment) %>%
             summarize(mean_value = mean(cell_viability, na.rm = TRUE), 
                       median_value = median(cell_viability, na.rm = TRUE), 
                       sd_value = sd(cell_viability, na.rm = TRUE)))


# Total number living cells per mg freshweight
model = lmer(total_number_living_cells_per_mg~treatment + (1|day), data = data_without_control)
plot(simulateResiduals(model))
hist(in_vivo_data$total_number_living_cells_per_mg)
hist(resid(model)) 
Anova(model)

PH = glht(model, mcp(treatment = "Tukey"))
summary(PH, test = adjusted("BH"))

data.frame(in_vivo_data %>%
             group_by(treatment) %>%
             summarize(mean_value = mean(total_number_living_cells_per_mg, na.rm = TRUE), 
                       median_value = median(total_number_living_cells_per_mg, na.rm = TRUE), 
                       sd_value = sd(total_number_living_cells_per_mg, na.rm = TRUE)))



## 5.2 Effects of MPP on coelomocyte subpopulations ----

# Number of peaks between treatments
number_peaks = data.frame(table(in_vivo_data$number_peaks, in_vivo_data$treatment))
colnames(number_peaks) = c("peak_number", "treatment", "frequency")
number_peaks = subset(number_peaks, treatment != "control")
number_peaks = data.frame(mock = subset(number_peaks, treatment == "mock")$frequency,
                          ps = subset(number_peaks, treatment == "ps")$frequency,
                          additive = subset(number_peaks, treatment == "additive")$frequency)

chisq.test(number_peaks) # overall test
chisq.test(number_peaks[, c(1,2)], simulate.p.value = TRUE, B = 1e6) # mock vs ps
chisq.test(number_peaks[, c(1,3)], simulate.p.value = TRUE, B = 1e6) # mock vs additive
chisq.test(number_peaks[, c(2,3)], simulate.p.value = TRUE, B = 1e6) # ps vs additive


# Percentage of eleocytes
model = lmer(eleozytes_percent~treatment + (1|day), data = data_without_control)
plot(simulateResiduals(model))
hist(in_vivo_data$eleozytes_percent)
hist(resid(model)) 
Anova(model)

data.frame(in_vivo_data %>%
             group_by(treatment) %>%
             summarize(mean_value = mean(eleozytes_percent, na.rm = TRUE), 
                       median_value = median(eleozytes_percent, na.rm = TRUE), 
                       sd_value = sd(eleozytes_percent, na.rm = TRUE)))


# Percentage of amoebocytes
model = lmer(amoebozytes_percent~treatment + (1|day), data = data_without_control)
plot(simulateResiduals(model))
hist(in_vivo_data$amoebozytes_percent)
hist(resid(model)) 
Anova(model)

data.frame(in_vivo_data %>%
             group_by(treatment) %>%
             summarize(mean_value = mean(amoebozytes_percent, na.rm = TRUE), 
                       median_value = median(amoebozytes_percent, na.rm = TRUE), 
                       sd_value = sd(amoebozytes_percent, na.rm = TRUE)))


# Percentage of intermediate population
model = lmer(intermediate_population_percent~treatment + (1|day), data = data_without_control)
plot(simulateResiduals(model))
hist(in_vivo_data$intermediate_population_percent)
hist(resid(model)) 
Anova(model)

data.frame(in_vivo_data %>%
             group_by(treatment) %>%
             summarize(mean_value = mean(intermediate_population_percent, na.rm = TRUE), 
                       median_value = median(intermediate_population_percent, na.rm = TRUE), 
                       sd_value = sd(intermediate_population_percent, na.rm = TRUE)))



## 5.3 Effects of ex vivo MPP exposure ----


# PS
ex_vivo_ps = subset(ex_vivo_data, treatment == "control_PS" | treatment == "PS")
ex_vivo_ps = subset(ex_vivo_ps, !is.na(absorption))
ex_vivo_ps$time_concentration = paste(ex_vivo_ps$timepoint, ex_vivo_ps$particle_concentration, sep = "_")
ex_vivo_ps$time_concentration = as.factor(ex_vivo_ps$time_concentration)

data.frame(ex_vivo_ps %>%
             group_by(particle_concentration,timepoint) %>%
             summarize(mean_value = mean(relative_metabolic_activity, na.rm = TRUE), 
                       sd_value = sd(relative_metabolic_activity, na.rm = TRUE),
                       percent = (sd_value/mean_value)*100))


kruskal.test(relative_metabolic_activity~time_concentration, data = ex_vivo_ps)
dunn.test(ex_vivo_ps$relative_metabolic_activity, ex_vivo_ps$time_concentration, 
          method = "bh", kw = F, list = T, table = F, altp = T)


# PS with Irgafos
ex_vivo_ps_irgafos = subset(ex_vivo_data, treatment == "control_PS_Irgafos" | treatment == "PS_Irgafos" & particle_concentration != 500)
ex_vivo_ps_irgafos = subset(ex_vivo_ps_irgafos, !is.na(absorption))
ex_vivo_ps_irgafos$time_concentration = paste(ex_vivo_ps_irgafos$timepoint, ex_vivo_ps_irgafos$particle_concentration, sep = "_")
ex_vivo_ps_irgafos$time_concentration = as.factor(ex_vivo_ps_irgafos$time_concentration)

data.frame(ex_vivo_ps_irgafos %>%
             group_by(particle_concentration,timepoint) %>%
             summarize(mean_value = mean(relative_metabolic_activity, na.rm = TRUE), 
                       sd_value = sd(relative_metabolic_activity, na.rm = TRUE),
                       percent = (sd_value/mean_value)*100))


kruskal.test(relative_metabolic_activity~time_concentration, data = ex_vivo_ps_irgafos)
dunn.test(ex_vivo_ps_irgafos$relative_metabolic_activity, ex_vivo_ps_irgafos$time_concentration, 
              method = "bh", kw = F, list = T, table = F, altp = T)



# 6 Supplement ----

# Worm weight 
weight_data = subset(in_vivo_data, !is.na(ending_weight))
weight_data = pivot_longer(weight_data, names_to = "weight_measurement", values_to = "weight", cols = c(starting_weight, ending_weight))
weight_data$weight_measurement = factor(weight_data$weight_measurement, levels = c("starting_weight", "ending_weight"))

ggplot(weight_data, aes(x = treatment, y = weight, fill = treatment, pattern = weight_measurement)) + 
  geom_boxplot_pattern(
    pattern_density = 0.05,
    pattern_angle = 45,
    pattern_spacing = 0.02,
    position = position_dodge(width = 0.8)) +
  scale_pattern_manual(values = c("none", "stripe")) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0),
              size = 1.5, color = 1) +
  guides(fill = "none", color = "none") +
  theme_classic() + 
  scale_fill_manual(values = c("#BDD7E7", "#6BAED6", "#2171B5")) +
  labs(x = "", y = "Earthworm weight [g]") +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 4, color = "gold2", 
               position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)) +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14),  legend.position = "none") +
  scale_x_discrete(labels = c("mock", "PS", expression(PS[irgafos]))) +
  scale_y_continuous(breaks = c(seq(0.2, 0.5, 0.1)), labels = c(seq(0.2, 0.5, 0.1)), expand = c(0,0), limits = c(0.14,0.54))


# SubG1 peak
ggplot(sub_g1_data, aes(x = treatment, y = percentage, fill = treatment, pattern = cell_population)) +  
  geom_boxplot_pattern(
    pattern_density = 0.05,
    pattern_angle = 45,
    pattern_spacing = 0.02,
    position = position_dodge(width = 0.8)) +
  scale_pattern_manual(values = c("none", "stripe")) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0),
              size = 1.5, color = 1) +
  labs(x = "", y = "Sub G1 [%]") +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 4, color = "gold2",
               position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)) +
  scale_x_discrete(labels = c("mock", "PS", expression(PS[irgafos]))) +
  scale_fill_manual(values = c("#BDD7E7", "#6BAED6", "#2171B5")) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.position = "none") 

