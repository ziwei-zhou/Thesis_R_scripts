devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")
install.deps(c("tidyverse", "xlsx", "ggplot2", "RColorBrewer", "car"))

# install.deps(c("EasyqpcR", "NormqPCR"), repo = "bioc")

# plot_theme <-  function(def_theme="grey", baseSize=22){
#   require("RColorBrewer")
#   my_theme <- eval(parse(text = sprintf("theme_%s(base_size=%d)", def_theme,baseSize)))  +
#     theme(axis.title.y=element_text(face="bold", vjust = 1.5, size=rel(0.8)),
#           axis.title.x=element_text(face="bold", vjust = 0.1, size=rel(0.8)),
#           legend.title=element_text(size=rel(0.8), hjust=0.5),
#           legend.text=element_text(size = rel(0.7),lineheight = 1.5),
#           #panel.grid.minor=element_blank(),
#           strip.text=element_text(size=rel(0.8)),
#           strip.background=element_rect(fill = "lightskyblue"))
#   return(my_theme)
# }

options(stringsAsFactors = FALSE)

#### Process LinRegPCR results #####
# read.xlsx(combined_file, sheetName = "Cor_data")
linreg_data <- xlsx::read.xlsx("RGA_qPCR_raw_data_18_01_2018.xlsx", sheetName = "combined_data_compact",
                            colIndex = 1:10, startRow = 4)  %>%
  mutate(Genotype= substring(tissue, 1, 1), Treatment=substring(tissue, 2, 2),
         hpi=as.numeric(substring(tissue, 3))) %>% arrange(Amplicon, Genotype, Treatment, hpi) %>%
  mutate(Bio_rep=rep(rep(c(rep(1:3, 6), rep(4:6, each=2)), length(unique(.$Genotype))), length(unique(.$Amplicon))), 
         Tech_rep=rep(rep(c(rep(rep(c("i", "ii"), each=3), length(unique(.$hpi))), rep(c("i", "ii"),  length(unique(.$hpi)))),  length(unique(.$Genotype))) , length(unique(.$Amplicon)))  , ID=paste(tissue, Bio_rep, Tech_rep, sep="_"))


# Extract Reference genes data and put side-by-side
# Ref1 data
Ref1_data <- linreg_data %>% filter(grepl("^Ref1", name), Cq<40, Cq>5)  %>% 
  filter(!str_detect(Quality_checks, "[1-4]")) %>% 
  dplyr::select(ID, Ref1_N0=N0, Ref1_Cq=Cq) 
# Combined Ref data
Ref_data <- linreg_data %>% filter(grepl("^Ref3", name), Cq<40, Cq>5)  %>% 
  filter(!str_detect(Quality_checks, "[1-4]")) %>% 
  dplyr::select(ID, Ref3_N0=N0, Ref3_Cq=Cq) %>%
  full_join(Ref1_data) %>% 
  mutate(mean_ref_Cq=rowMeans(cbind(Ref1_Cq, Ref3_Cq), na.rm = TRUE), 
         mean_ref_N0=rowMeans(cbind(Ref1_N0, Ref3_N0), na.rm = TRUE))

# Compare reference genes
Cq_cols <- colnames(Ref_data)[grepl("Cq", colnames(Ref_data))]
# se <- function(x) {
#   sqrt(var(x, na.rm=TRUE)/length(x[!is.na(x)]))
# }

Ref_stats <- sapply(Cq_cols, function(col) c(mean(Ref_data[,col], na.rm=TRUE), se(Ref_data[,col])))
rownames(Ref_stats) <- c("Mean", "SE")

# Prepare Target genes ddCq data

RGA_data <- linreg_data %>% filter(grepl("^RGA", name), Cq<40, Cq>5) %>% 
  filter(!str_detect(Quality_checks, "[1-4]")) %>% 
  full_join(Ref_data[c("ID", "mean_ref_Cq", "mean_ref_N0")]) %>%
  mutate(dCq=mean_ref_Cq-Cq, dN0=mean_ref_N0-N0)
# Summarise PCR efficiencies
PCR_eff <- linreg_data %>% group_by(Amplicon) %>% summarise(PCR_eff=mean(mean_PCR_eff)) %>%
  ungroup()
PCR_eff_dict <- setNames(PCR_eff$PCR_eff, PCR_eff$Amplicon)

#### Calculate significance between T and U for each gene ####
RGAs <- PCR_eff %>% filter(!grepl("Ref", Amplicon)) %>% .$Amplicon
signif_df <- NULL
for (target in RGAs){
  # target="RGA12"
  target_data <- RGA_data %>% filter(Amplicon==target) %>% select(tissue, dCq)
  untreated_data <- target_data %>% filter(grepl("U", tissue)) %>% group_by(tissue) %>%
    summarise(mean_dCq=mean(dCq, na.rm=TRUE)) 
  U_dict <- setNames(untreated_data$mean_dCq, sub("U", "", untreated_data$tissue))
  
  target_ddcq <- target_data %>% filter(!grepl("U", tissue), !is.na(dCq)) %>% 
    mutate(U=U_dict[sub("T", "", tissue)], ddCq=dCq-U, tissue=as.factor(tissue)) 
  # Compute the analysis of variance
  res.aov <- aov(ddCq ~ tissue, data = target_ddcq)
  # Summary of the analysis
  aov_sum <- summary(res.aov)
  pvals <- aov_sum[[1]][[1, "Pr(>F)"]]
  
  # Check Homogeneity of variances (P>0.05) Levene
  L_test <- leveneTest(ddCq ~ tissue, data = target_ddcq)
  L_pval <- L_test[1, "Pr(>F)"]
  # Extract the residuals
  aov_residuals <- residuals(object = res.aov )
  # Run Shapiro-Wilk test (if p-value>0.05, then passes Normality test)
  S_test <- shapiro.test(x = aov_residuals )
  if (S_test$p.value<0.05) {
    K_test <- kruskal.test(ddCq ~ tissue, data = target_ddcq)
    pvals <- K_test$p.value
  }
  
  if (pvals<0.05 & L_pval>0.05) {
    # Tukey multiple pairwise-comparisons
    tukey_test <- TukeyHSD(res.aov, conf.level = 0.05)$tissue %>% as.data.frame(.) %>% 
      rownames_to_column("Comparison") %>% 
      mutate(signif=symnum(`p adj`, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                           symbols = c("***", "**", "*", ".", "")), target=target)
    
  } else {
    tukey_test <- c(rep(NA, 5), "Gene data does not fulfill assumptions for ANOVA test", target)
  }
  signif_df <- rbind(signif_df, tukey_test)
}


contrasts <- c('HT24-HT2','HT6-HT2','IT2-HT2','ST2-HT2','HT6-HT24','IT24-HT24','KT24-HT24','ST24-HT24','IT6-HT6','KT6-HT6','ST6-HT6','IT6-IT2','KT2-IT2','ST2-IT2','IT6-IT24','KT24-IT24','ST24-IT24','KT6-IT6','ST6-IT6','KT24-KT2','KT6-KT2','ST2-KT2','KT6-KT24','ST24-KT24','ST6-KT6','ST24-ST2','ST6-ST2','ST6-ST24')
gene_signif <- signif_df %>% filter(Comparison %in% contrasts) %>% 
  .[c("Comparison", "signif", "target")] %>%
  spread(key="target", value="signif")  
xlsx::write.xlsx(gene_signif, filedate("statistical_analysis_ddCq_values", ".xlsx"), row.names = FALSE,
                                                        sheetName = "signif")

signif_df %>% filter(grepl("ANOVA", signif))

# , hpi==time_point
xlsx::write.xlsx(signif_df, filedate("statistical_analysis_ddCq_values", ".xlsx"), row.names = FALSE,
                  sheetName = "stats", append = TRUE)


