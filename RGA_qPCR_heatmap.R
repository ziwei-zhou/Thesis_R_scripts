devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")
install.deps(c("tidyverse", "xlsx", "ggplot2", "RColorBrewer", "gplots"))

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

#### load ddCq results from file #####
ddCq_data <- readxl::read_excel("ddCq_from_LinRegQPCR_data_07_02_2018.xlsx", 
           sheet = "ddCq_data")
ddcq_mat <- as.data.frame(ddCq_data) %>% unite("Condition", Genotype, hpi) %>% 
  select(Amplicon, Condition,  ddCq) %>% 
  spread(key = Condition, value = ddCq) %>% write_tsv(., path = "ddCq_data_matrix.txt") %>% 
  remove_rownames() %>% column_to_rownames("Amplicon") %>%
  as.matrix()

# define a few color palettes  
my_cols <- colorspace::diverge_hcl(75, h = c(255, 330), l=c(40,90)) # c('#d01c8b','#f1b6da','#f7f7f7','#b8e186','#4dac26')
my_cols <- brewer.pal(11, "RdYlGn")
my_cols <- redgreen(75)
pdf(filedate("ddCq_heatmap_blue_pink", ".pdf", "plots"),width=6, height=8)
heatmap.2(t(ddcq_mat), col=my_cols, trace = "none", density.info="none")
dev.off()
  

