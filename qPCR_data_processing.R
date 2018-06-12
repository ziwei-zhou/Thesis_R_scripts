devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")
install.deps(c("tidyverse", "xlsx", "ggplot2", "RColorBrewer"))

install.deps(c("EasyqpcR", "NormqPCR"), repo = "bioc")

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
#### Combine Cq results from multiple files (plates) #####
qPCR_files <- list.files(path = "./", pattern = "*Quantification Summary.xlsx", full.names = TRUE, recursive = TRUE) %>% .[!grepl("~$", ., fixed = TRUE)] 

combined_file <- filedate("Shane_RGA_qPCR_Cq_data", ".xlsx")

if (file.exists(combined_file)){
  LogMsg(sprintf("file '%s' already exists!", combined_file))
  x <- NULL
  x <- readline("Do you wish to overwrite it? (Yes/No) ")
  if (tolower(substr(x, 1, 1))=="n")  {
    combined_file <- readline("Please enter a new file name: ")
  } else unlink(combined_file)
} 
LogMsg(sprintf("Writing qPCR data from plates as sheets in file '%s'", combined_file))
Cq_data <- NULL
for (f in qPCR_files){
  
  plate_name <- gsub(".", "_", sub(" -  Quantification Summary.xlsx", "", basename(f)), fixed = TRUE)
  s <- xlsx::read.xlsx(f, sheetName = "0", colIndex = 1:7) %>% rename_("Plate"="NA.") %>% 
    mutate(Plate=plate_name, Sample=gsub(".","_", Sample, fixed = TRUE))
  append_to_file <- ifelse(file.exists(combined_file), TRUE, FALSE)
  write.xlsx(s, combined_file, sheetName = plate_name, showNA = FALSE,
               col.names = TRUE, row.names = FALSE, append = append_to_file)
  Cq_data <- bind_rows(Cq_data, s)
  LogMsg(sprintf("Finished processing data from plate '%s'.", plate_name))
}

# Append the combined data to one sheet
write.xlsx(Cq_data, combined_file, sheetName = "combined_data", showNA = FALSE,
           col.names = TRUE, row.names = FALSE, append = TRUE)
LogMsg(sprintf("Combined all data from plates into file '%s'.", combined_file))

# Summarise factor-corrected data (use either Factor_qPCR or NormqPCR) and plot reference genes
ref_genes <- c("Ref1", "Ref3")
# read.xlsx(combined_file, sheetName = "Cor_data")
##### Read Corrected Cq data #####
cor_data <- xlsx::read.xlsx("RGA_dCq_data.xlsx", sheetName = "Cor_data", colIndex = 1:10)  
plot_data <- cor_data %>% 
  filter(r_Cq>10, !grepl("Pool", Sample), !is.na(Sample)) %>% group_by(Target, Sample, Bio_rep) %>% summarise(mean_Cq=mean(r_Cq), Cq_se=se(r_Cq)) %>% ungroup() %>% 
  mutate(Genotype= substring(Sample, 1, 1), Treatment=substring(Sample, 2, 2),
                                             hpi=as.numeric(substring(Sample, 3)))
##### Plot RGA genes ddCq #####
plot_genes <- unique(plot_data$Target)#paste0("RGA", seq(4,12,2))
for (p in plot_genes){
   plot_se <- plot_data %>% dplyr::select(-mean_Cq,-Sample) %>%  filter(Target %in% p) %>% 
      spread(Treatment, Cq_se) %>% mutate(U=if_else(is.na(U), 0, U), Diff_se=U+T)
   
   plot_refs <- plot_data %>% dplyr::select(-Cq_se,-Sample) %>%  filter(Target %in% p) %>% 
      spread(Treatment, mean_Cq) %>% mutate(Diff=U-T,Diff_se=plot_se$Diff_se )
   
   # group_by(Target, Genotype, hpi) %>% summarise(Cq_se=se(mean_Cq), mean_Cq=mean(mean_Cq)) %>% 
   #   ungroup() 
   
   ggplot(plot_refs, aes(x=hpi, y=Diff, color=Genotype)) +  # , line=Treatment
      geom_point(size=3) +  scale_color_brewer(palette = "Set1") +  geom_line() +
      scale_x_continuous("Time post inoculation", breaks=c(2, 6, 24)) + 
      labs(y="dCq", title=sprintf("%s differential Cq plot", p)) +
      geom_errorbar(data = plot_refs, aes(x=hpi, ymin=Diff-Diff_se, ymax=Diff+Diff_se), width=0.25, alpha=0.7) +
      plot_theme(baseSize = 24) 
   ggsave(filedate(sprintf("%s_dCq_plot", p), ".pdf", "plots"))
}


##### Prepare factor-corrected data for loading to REST   #####
cor_data %>% filter(!grepl("NTC", Content), Sample!="Pool") %>% 
  mutate(Genotype= substring(Sample, 1, 1), Treatment=substring(Sample, 2, 2),
         hpi=as.numeric(substring(Sample, 3)))  %>%
  group_by( Target, Genotype, Treatment, hpi , Bio_rep)   %>% 
  summarise(mean_Cq=mean(r_Cq,na.rm = TRUE)) %>% 
  tidyr::spread(Target, mean_Cq ) %>% arrange(Genotype, hpi, Bio_rep, Treatment) %>%
  write_tsv(., file.path(getwd(),"RGA_corrected_qPCR_data_for_REST.txt"))

#### Prepare factor-corrected data for loading to ReadqPCR  #####
cor_data %>% filter(!grepl("NTC", Content)) %>% 
  mutate(Sample=sub("^[^ ]+", "", comb_cond))  %>% select(Well, Plate, Sample, Detector=Target, Cq=r_Cq) %>%
  write_tsv(., "corrected_qPCR_data_matrix.txt")
qPCRBatch.qPCR <- read.qPCR("corrected_qPCR_data_matrix.txt")

#### Prepare data for LinRegPCR  ####
qPCR_files <- list.files(path = "./", pattern = "*Quantification Summary.xlsx", full.names = TRUE, recursive = TRUE) %>% .[!grepl("~$", ., fixed = TRUE)] 
combined_file <- filedate("RGA_qPCR_raw_data", ".xlsx")

if (file.exists(combined_file)){
  LogMsg(sprintf("file '%s' already exists!", combined_file))
  x <- NULL
  x <- readline("Do you wish to overwrite it? (Yes/No) ")
  if (tolower(substr(x, 1, 1))=="n")  {
    combined_file <- readline("Please enter a new file name: ")
  } else unlink(combined_file)
} 
LogMsg(sprintf("Writing qPCR data from plates as sheets in file '%s'", combined_file))
raw_data <- NULL

# f=qPCR_files[4]
for (f in qPCR_files){
  gene_raw_data <- NULL
  # f <- qPCR_files[3]
  # gene <- substr(basename(f), 1, str_locate(basename(f), "\\(")[1]-1) #-2
  plate_name <- gsub(".", "_", sub(" -  .+", "", basename(f)), fixed = TRUE)
  s <- xlsx::read.xlsx(f, sheetName = "0", colIndex = 1:7) %>% rename_("Plate"="NA.") %>% 
    #filter(!grepl("NTC", Content), Sample!="Pool") %>%
    mutate(Plate=plate_name, Sample=paste(Target, Sample, sep="-"))
             
             #ifelse((grepl("^Ref",Target) && is.na(Sample)), paste(Target, "Pool", sep = " "), 
                #         paste(Target, Sample, sep="-")))
  samples <- s$Sample
  qPCR_raw_files <- list.files(path = "./", pattern = "*Quantification Amplification Results.xlsx", full.names = TRUE, recursive = TRUE)
  qPCR_raw_file <- qPCR_raw_files[grepl(plate_name, qPCR_raw_files, fixed = TRUE)]
  gene_raw_data <- rbind(c(NA, "Cycle", samples), xlsx::read.xlsx(qPCR_raw_file, sheetName = "SYBR"))
  if (is.null(raw_data)) raw_data <- gene_raw_data[1:2]
  
  # colnames(gene_raw_data)[3:ncol(gene_raw_data)] <- samples
  export_data <- gene_raw_data[, grepl("R.+-.+\\d$",gene_raw_data[1,] )]
  append_to_file <- ifelse(file.exists(combined_file), TRUE, FALSE)
  
  write.xlsx(gene_raw_data, combined_file, sheetName = plate_name, showNA = FALSE,
             col.names = TRUE, row.names = FALSE, append = append_to_file)
  raw_data <- bind_cols(raw_data, export_data)
  LogMsg(sprintf("Finished processing data from plate '%s'.", plate_name))
}

# Append the combined data to one sheet
pheno_data <- data.frame(Samples=as.character(raw_data[1,3:ncol(raw_data)])) %>% separate(., Samples, c("Amplicon", "Group"),sep = "-", remove = FALSE) %>% 
   mutate(Genotype= substring(Group, 1, 1), Treatment=substring(Group, 2, 2), hpi=as.numeric(substring(Group, 3)))

write.xlsx(pheno_data, combined_file, sheetName = "pheno_data", showNA = FALSE,
           col.names = TRUE, row.names = FALSE, append = TRUE)
write.xlsx(raw_data, combined_file, sheetName = "combined_data", showNA = FALSE,
           col.names = FALSE, row.names = FALSE, append = TRUE)
LogMsg(sprintf("Combined all data from plates into file '%s'.", combined_file))

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
# Prepare side-by-side treatment for ddCq
ddCq_data <- RGA_data %>% 
  group_by(Amplicon, Genotype, Treatment, hpi) %>%
  summarise(mean_dCq=mean(dCq, na.rm=TRUE)) %>%
  spread(Treatment, mean_dCq) %>% ungroup() %>%
  mutate(se_dCq=(RGA_data %>% 
           group_by(Amplicon, Genotype, hpi) %>%
           summarise(SE=se(dCq)) %>% .$SE ), ddCq=T-U  )
# Write to excel
write.xlsx(as.data.frame(ddCq_data), filedate("ddCq_from_LinRegQPCR_data", ".xlsx"), 
           sheetName = "ddCq_data", showNA = FALSE,
           col.names = TRUE, row.names = FALSE, append = FALSE)
#### Plot ddCq data  #####
plot_genes <- unique(ddCq_data$Amplicon)#paste0("RGA", seq(4,12,2))
for (p in plot_genes){
  # p=plot_genes[1]
  plot_ddCq_data <- ddCq_data %>% filter(Amplicon %in% p)
  ggplot(plot_ddCq_data, aes(x=hpi, y=ddCq, color=Genotype)) +  # , line=Treatment
  geom_point(size=3) +  scale_color_brewer(palette = "Set1") +  geom_line() +
  scale_x_continuous("Time post inoculation [hpi]", breaks=c(2, 6, 24)) + 
  labs(y="Relative expression [ddCq]", title=sprintf("%s differential expression plot",p),
       subtitle="(Treated vs. Un-treated)") +
  geom_errorbar(aes(x=hpi, ymin=ddCq-se_dCq, ymax=ddCq+se_dCq), 
                width=0.25, alpha=0.7) + # data = plot_ddCq_data, 
  plot_theme(baseSize = 24) 
  ggsave(filedate(sprintf("%s_ddCq_plot", p), ".pdf", "plots"))
}
colnames(RGA_data)


