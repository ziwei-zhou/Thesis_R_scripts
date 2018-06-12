install.packages("devtools") # install if needed
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")
install.deps(c("tidyverse", "broom", "RColorBrewer", "tidyr", "ggrepel"))
install.deps("NormqPCR", repo = "bioc")
devtools::install_github("slowkow/ggrepel@0.6.12", force=TRUE)


# Define plotting theme
plot_theme <-  theme_grey(base_size=22) +
  theme(axis.title.y=element_text(face="bold", vjust = 1.5, size=rel(0.8)),
        axis.title.x=element_text(face="bold", vjust = 0.1, size=rel(0.8)),
        legend.title=element_text(size=rel(0.8), hjust=0.5),
        legend.text=element_text(size = rel(0.7),lineheight = 1.5),
        #panel.grid.minor=element_blank(),
        strip.text=element_text(size=rel(0.8)),
        strip.background=element_rect(fill = "lightskyblue"))#, co

##### Primer Efficiency ######
qpcr_plate <- readxl::read_excel("RGA3_14_211117/RGA3_14_Quantification_Cq_Results.xlsx")

sum_qpcr <- qpcr_plate %>% mutate(LogQty=log10(`Starting Quantity (SQ)`)) %>% 
  filter(is.na(Ignore)) %>% # 
  group_by(Target, Content, LogQty) %>% 
  summarise(Cq_Mean=mean(Cq, na.rm=TRUE), Cq_Sd=sd(Cq, na.rm=TRUE)) %>%
  filter(!is.na(Cq_Mean), Content!="NTC") # (Cq_Sd<1 | is.na(Cq_Sd)), 

gene_model <- sum_qpcr %>% group_by(Target) %>% do(model = lm(Cq_Mean ~ LogQty, data = .))

tidyModel <- tidy(gene_model, model) %>% ungroup() %>% filter(term=="LogQty") %>% 
  mutate(Intercept=sapply(gene_model$model, function(m) m$coefficients[1]),  
         Efficiency=10^(-1/estimate)-1, R2=sapply(gene_model$model, function(m) summary(m)$r.squared))


plot_genes <- unique(sum_qpcr$Target)[1:6]
plot_data <- sum_qpcr %>% filter(Target %in% plot_genes)

primer_annot <- tidyModel %>% filter(Target %in% plot_genes) 
primer_annot$eqn <- sprintf("italic(y) == \"%.2f\" + \"%.2f\" * italic(x) * \",\" ~ ~italic(r)^2 ~ \"=\" ~ \"%.2f\" * \",\" ~ ~italic(Eff) ~ \"=\" ~ \"%.2f\"", primer_annot$Intercept, primer_annot$estimate, primer_annot$R2, primer_annot$Efficiency) #   
# Produce plots
ggplot(plot_data, aes(x=LogQty, y=Cq_Mean, shape=Target, colour=Target)) +
  geom_point(size=3) +
  # scale_shape_manual(values=c(1,2)) +
  scale_colour_brewer(palette="Dark2") + 
  # stat_smooth(method=lm, se=FALSE , colour="black") +
  geom_errorbar(aes(ymin=Cq_Mean-Cq_Sd, ymax=Cq_Mean+Cq_Sd), width=0.05) +
  geom_text_repel(aes(x = 0.25, y = 25, label=eqn), segment.alpha = 0,
                  data=primer_annot, parse = TRUE, nudge_y = c(0.5)) + # repel , segment.size = 0,  segment.alpha=1, 
  geom_smooth(method=lm, se=FALSE) + labs(y = "Mean Cq") + plot_theme

ggsave(filedate("RGA3-8_calibration_curves", ".pdf" ,"plots" ), width=11, height=8)
# RGA9-13_calibration_curves # RGA3-8_calibration_curves

#### qPCR Plate Quantification ####
qpcr_plate <- readxl::read_excel("RGA123_2017-10-20 13-41-41_CT026919 -  Quantification Cq Results.xlsx", "0")
Sample <- strsplit(qpcr_plate$Sample, split="")
gene="RGA3"
qpcr_data <- qpcr_plate %>% filter(Target==gene, Cq>0) %>% select(Target, Sample, Cq) %>%
  separate(Sample, into=c("Genotype", "Treatment", "Time"), sep=c(1,2), extra="merge") %>% mutate(Time=as.numeric(Time))




# Produce plots
ggplot(qpcr_data, aes(x=Time, y=Cq, lty=Treatment, colour=Genotype, shape=Genotype)) +
  geom_point(size=3.5) + scale_colour_brewer(palette="Set1") + geom_line() + 
  scale_x_continuous(breaks = unique(qpcr_data$Time)) +
  labs(x = "Time (hpi)") + plot_theme
  # scale_shape_manual(values=c(1,2)) 

ggsave(filedate("RGA3_test_run", ".pdf" ,"plots" ), width=8, height=6)



##### Primer Efficiency ######
sum_qpcr <- qpcr_plate2 %>% rename(LogQty=`Log Starting Quantity`) %>% 
  group_by(Target, Content, LogQty) %>% 
  summarise(Cq_Mean=mean(Cq, na.rm=TRUE), Cq_Sd=sd(Cq, na.rm=TRUE)) %>%
  filter(Cq_Sd<1, !is.na(Cq_Mean), Content!="NTC")

gene_model <- sum_qpcr %>% group_by(Target) %>% do(model = lm(Cq_Mean ~ LogQty, data = .))
# summary(gene_model$model[1][[1]])
# 
# model <- gene_model$model[1][[1]]
# model$coefficients
# 
# sapply(gene_model$model, function(m) format(summary(m)$r.squared, digits=2))

tidyModel <- tidy(gene_model, model) %>% ungroup() %>% filter(term=="LogQty") %>% 
  mutate(Intercept=sapply(gene_model$model, function(m) m$coefficients[1]),  Efficiency=exp(-1/estimate), R2=sapply(gene_model$model, function(m) summary(m)$r.squared))

good_genes <- c("PR2", "RING", "RBP")
plot_data <- sum_qpcr %>% filter(Target %in% good_genes)



primer_annot <- tidyModel %>% filter(Target %in% good_genes) 
primer_annot$eqn <- sprintf("italic(y) == \"%.2f\" + \"%.2f\" * italic(x) * \",\" ~ ~italic(r)^2 ~ \"=\" ~ \"%.2f\" * \",\" ~ ~italic(Eff) ~ \"=\" ~ \"%.2f\"", primer_annot$Intercept, primer_annot$estimate, primer_annot$R2, primer_annot$Efficiency) #   
# Produce plots
ggplot(plot_data, aes(x=LogQty, y=Cq_Mean, shape=Target, colour=Target)) +
  geom_point(size=3) +
  # scale_shape_manual(values=c(1,2)) +
  scale_colour_brewer(palette="Set1") + 
 # stat_smooth(method=lm, se=FALSE , colour="black") +
  geom_errorbar(aes(ymin=Cq_Mean-Cq_Sd, ymax=Cq_Mean+Cq_Sd), width=0.1) +
  geom_text_repel(aes(x = 2, y = 31, label=eqn), data=primer_annot, parse = TRUE, segment.size = 0,
                  segment.alpha=1, nudge_x = c(0,0.5,0.5)) +
  geom_smooth(method=lm, se=FALSE)

ggsave(filedate("Plate2_good_genes", ".pdf" ,"plots" ), width=8, height=6)

