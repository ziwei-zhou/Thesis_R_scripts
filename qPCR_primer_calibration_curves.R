install.deps(c("tidyverse", "broom", "RColorBrewer", "tidyr", "ggrepel"))
# devtools::install_github("slowkow/ggrepel@0.6.12", force=TRUE)

# Define plotting theme
eff_plot_theme <-  theme_grey(base_size=22) +
  theme(axis.title.y=element_text(face="bold", vjust = 1.5, size=rel(0.8)),
        axis.title.x=element_text(face="bold", vjust = 0.1, size=rel(0.8)),
        legend.title=element_text(size=rel(0.8), hjust=0.5),
        legend.text=element_text(size = rel(0.7),lineheight = 1.5),
        #panel.grid.minor=element_blank(),
        strip.text=element_text(size=rel(0.8)),
        strip.background=element_rect(fill = "lightskyblue"))#, co

pal <- brewer.pal(7, "Set1")[-6]
##### Primer Efficiency ######

plot_primer_eff <- function(plate_file, gene_range, plot_name, eq_x=0.25, eq_y=25,
                            plot_width=11, plot_height=8, plot_theme=eff_plot_theme, plot_pal=NULL){
  qpcr_plate <- readxl::read_excel(plate_file)
  
  sum_qpcr <- qpcr_plate %>% mutate(LogQty=log10(`Starting Quantity (SQ)`)) %>% 
    filter(is.na(Ignore)) %>% # 
    group_by(Target, Content, LogQty) %>% 
    summarise(Cq_Mean=mean(Cq, na.rm=TRUE), Cq_Sd=sd(Cq, na.rm=TRUE)) %>%
    filter(!is.na(Cq_Mean), Content!="NTC") # (Cq_Sd<1 | is.na(Cq_Sd)), 
  # Check number of points for each gene
  for (g in unique(sum_qpcr$Target)){
    target_sum <- sum_qpcr %>% filter(Target==g)
    if (nrow(target_sum) < 3) {
      cat(sprintf("Warning: Gene %s has less than 3 points, do you want to remove it from the plot?\n", 
              g), file=stderr())
      response <- readline("Press (Y/N) then Enter: ")
      if (toupper(response)=="Y") {
        sum_qpcr <- sum_qpcr %>% filter(Target!=g)
      }
    }
  }
  
  gene_model <- sum_qpcr %>% group_by(Target) %>% do(model = lm(Cq_Mean ~ LogQty, data = .))
  
  tidyModel <- tidy(gene_model, model) %>% ungroup() %>% filter(term=="LogQty") %>% 
    mutate(Intercept=sapply(gene_model$model, function(m) m$coefficients[1]),  
           Efficiency=10^(-1/estimate)-1, R2=sapply(gene_model$model, function(m) summary(m)$r.squared))
  # Defining plotting palette if it wasn't supplied
  if (is.null(plot_pal)) plot_pal <- brewer.pal(length(gene_range), "Dark2")
  if (gene_range[length(gene_range)]>length(unique(sum_qpcr$Target))) {
    stop(sprintf("Error: gene range exceeds maximum number of available genes after filtration (%d)", length(unique(sum_qpcr$Target))))
  }
  plot_genes <- unique(sum_qpcr$Target)[gene_range]
  plot_data <- sum_qpcr %>% filter(Target %in% plot_genes)
  
  primer_annot <- tidyModel %>% filter(Target %in% plot_genes) 
  primer_annot$eqn <- sprintf("italic(y) == \"%.2f\" - \"%.2f\" * italic(x) * \",\" ~ ~italic(R)^2 ~ \"=\" ~ \"%.2f\" * \",\" ~ ~italic(Eff) ~ \"=\" ~ \"%.2f\"", primer_annot$Intercept, 
                              -1*primer_annot$estimate, primer_annot$R2, primer_annot$Efficiency) #   
  # Produce plots
  ggplot(plot_data, aes(x=LogQty, y=Cq_Mean, shape=Target, colour=Target)) +
    geom_point(size=3) +
    # scale_shape_manual(values=c(1,2)) +
    scale_colour_manual(values=plot_pal) + 
    # stat_smooth(method=lm, se=FALSE , colour="black") +
    geom_errorbar(aes(ymin=Cq_Mean-Cq_Sd, ymax=Cq_Mean+Cq_Sd), width=0.05) +
    geom_text_repel(aes(x = eq_x, y = eq_y, label=eqn), segment.alpha = 0,
                    data=primer_annot, parse = TRUE, nudge_y = c(0.5)) + # repel , segment.size = 0,  segment.alpha=1, 
    geom_smooth(method=lm, se=FALSE) + labs(y = "Mean Cq") + plot_theme
  
  ggsave(plot_name, width=plot_width, height=plot_height)
}

# Run everything above this point

###  Plot calibration curves  ####
# For RGAs 03-08
plot_primer_eff("RGA3_14_211117/RGA3_14_Quantification_Cq_Results.xlsx", gene_range = 1:6, plot_name = filedate("RGA3-8_calibration_curves", ".pdf" ,"plots" ), plot_pal = pal)
# For RGAs 09-13
plot_primer_eff("RGA3_14_211117/RGA3_14_Quantification_Cq_Results.xlsx", gene_range = 7:11, plot_name = filedate("RGA9-13_calibration_curves", ".pdf" ,"plots" ), plot_pal = pal)
# For RGAs 15-24
plot_primer_eff(plate_file = "RGA15_25 efficiency test_221117/RGA15_25_test_221117_Cq_Results.xlsx", gene_range = 1:8, plot_name = filedate("RGA15-24_calibration_curves_Dark", ".pdf" ,"plots" ), plot_pal = brewer.pal(8, "Dark2"))
# For RGAs 23-25
plot_primer_eff(plate_file = "RGA15-25 test_221117 -  Quantification Cq Results.xlsx", gene_range = 7:9, plot_name = filedate("RGA23-25_calibration_curves", ".pdf" ,"plots" ), plot_pal = pal)
# For Ref 1-3
plot_primer_eff(plate_file = "Ref123 efficiency test_231117/Ref123_efficiency_test_231117Cq_Results.xlsx", gene_range = 1:3, eq_x = 0.5, eq_y = 22, plot_name = filedate("Ref1-3_calibration_curves", ".pdf" ,"plots" ), plot_pal = pal)
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




