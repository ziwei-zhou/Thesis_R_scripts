#library(dplyr)
#library(sqldf)
install.deps(c("dplyr", "tidyverse", "sqldf", "devtools"))
options(stringsAsFactors = FALSE)
# Update reutils to the recent version (>0.2.3), to handle the new https addresses at NCBI
install.deps("gschofl/reutils", "git")

# Install and load bioconductor packages
bioc_packages <- c("GenomicRanges", "Biostrings")
install.deps(bioc_packages, repo="bioc")
# Reload dplyr (must be last package loaded)
detach("package:dplyr", unload=TRUE)
library(dplyr)

# Read primer sequences
new_primers <- readxl::read_excel("RGA information.xlsx", sheet = "RGA11-28")
new_primers_fasta <- new_primers[1:2] %>% mutate(RGA_family=paste(RGA_family, "fw", sep="_")) %>% rename(Primer_name=RGA_family, Sequence=Fw_primer) %>% bind_rows(., data.frame(Primer_name=paste(new_primers$RGA_family, "rv", sep="_"), Sequence=new_primers$Rv_primer)) %>% arrange(Primer_name)
# Convert to fasta format
write.fasta(new_primers_fasta,filename = "RGA11-28_primers.fasta")

# Read input files
primers_blast <- read.delim("RGA_family_primers_blastn_results.txt", comment.char = "#", header = FALSE)
colnames(primers_blast) <- c("query_id", "subject_ids", "query_acc.ver", "subject_acc.ver", "identity", "alignment_length", "mismatches", "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit_score")

new_primers_blast <- read.delim("RGA_11-28_primers_blastn_results.txt", header = FALSE)
colnames(new_primers_blast) <- c("query_id", "subject_ids", "query_acc.ver", "subject_acc.ver", "identity", "alignment_length", "mismatches", "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit_score")
new_primers_blast <- rbind(new_primers_blast, primers_blast)

filtered_blast <- new_primers_blast %>%  group_by(query_id) %>% filter(bit_score==max(bit_score)) %>%  ungroup(.) %>% mutate(RGA_family=gsub("_[A-z][a-z]", "", .$query_id)) %>% arrange(RGA_family, subject_acc.ver)

#write.table(filtered_blast, "filtered_RGA_family_primers_blast_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
accesions <- NULL
for (i in 2:nrow(filtered_blast)) {
  if (filtered_blast$subject_acc.ver[i]==filtered_blast$subject_acc.ver[i-1]){
    accesions <- bind_rows(accesions, filtered_blast[i,c("RGA_family", "subject_ids", "subject_acc.ver")])
  }
}

accesions <- unique(accesions)
# acc.id=c("AC145027.17", "CT010459.8")
get_descr <- function(acc.id) {
  p <- efetch(acc.id, "nuccore", rettype="fasta")
  desc <- p$xmlValue("//TSeq_defline")
  return(desc)
}


accesions <- accesions %>% mutate(description=get_descr(.$subject_acc.ver))
filtered_primers <- filtered_blast %>% filter(.$subject_acc.ver %in% accesions$subject_acc.ver) %>% mutate(description=accesions$description[match(.$subject_acc.ver, accesions$subject_acc.ver)])

write.table(filtered_primers, filedate("matching_all_RGA_family_primer_pairs_blast_results", ext = ".txt"), sep = "\t", row.names = FALSE, quote = FALSE,append = FALSE)

unique(filtered_primers)

#################### get the amplicon sequence between the primers ############
#(use min and max of the range)

primers <- read.delim("matching_all_RGA_family_primer_pairs_blast_results_27_02_2017.txt", stringsAsFactors = FALSE)
i=1
positions <- c(primers$s.start[i], primers$s.start[i+1], primers$s.end[i], primers$s.end[i+1] )
targets <- data.frame("subject_acc.ver"=primers$subject_acc.ver[i], "description"=primers$description[i], "Fw_primer"=primers$query_id[i], "Rv_primer"=primers$query_id[i+1], "Amplicon_start"=min(positions), "Amplicon_end"=max(positions), "Amplicon_strand"=ifelse(positions[1]>positions[3], "-", "+"))

for (i in seq(3,nrow(primers), by=2)){
  positions <- c(primers$s.start[i], primers$s.start[i+1], primers$s.end[i], primers$s.end[i+1] )
  targets <- bind_rows(targets, data.frame("subject_acc.ver"=primers$subject_acc.ver[i], "description"=primers$description[i], "Fw_primer"=primers$query_id[i], "Rv_primer"=primers$query_id[i+1], "Amplicon_start"=min(positions), "Amplicon_end"=max(positions), "Amplicon_strand"=ifelse(positions[1]>positions[3], "-", "+")))
}

# acc.id=targets$subject_acc.ver[1]
# start_pos=targets$Amplicon_start[1]
# end_pos=targets$Amplicon_end[1]
# strand=targets$Amplicon_strand[1]
get_amplicon <- function(acc.id, start_pos, end_pos, strand) {
  acc_info <- efetch(acc.id, "nuccore", rettype="gb", retmode = "xml")
  acc_features <- setNames(acc_info$xmlValue("//GBFeature_location"),acc_info$xmlValue("//GBFeature_key"))
  #p <- efetch(acc.id, "nuccore", rettype="fasta")
  sequence <- toupper(acc_info$xmlValue("//GBSeq_sequence"))

  amplicon <- substr(sequence,start_pos,end_pos)
  amp_seq <- ifelse(strand=="+", amplicon, chartr("ATGC","TACG", reverse(amplicon)))
  cdsString <- unlist(strsplit(acc_features["CDS"], split = ",", fixed=TRUE))
  cds <- gsub("[^0-9\\.]", "", cdsString)
  in_exon <- "Unknown"
  if (!is.na(cds)) {
    exons <- sapply(cds, function(x) as.numeric(unlist(strsplit(x, split = "..", fixed = TRUE))))
    amp_cross_exons <- apply(exons, 2, function(clm) (start_pos>clm[1] & end_pos<clm[2]) )
    in_exon <- ifelse(any(amp_cross_exons), "Yes", "No")
  }
  return(setNames(amp_seq, in_exon))
}



amplicons <- sapply(1:nrow(targets), function(j) get_amplicon(targets$subject_acc.ver[j], targets$Amplicon_start[j], targets$Amplicon_end[j], strand=targets$Amplicon_strand[j]), USE.NAMES = FALSE)

target_amplicons <- targets %>%
  mutate(Amplicon_seq=amplicons, in_exon=names(amplicons), amp_length=nchar(amplicons)) %>% distinct(Fw_primer, Rv_primer, Amplicon_seq, .keep_all = TRUE)

write.table(target_amplicons, filedate("RGA_all_family_blast_results_unique_amplicons", ".txt"), sep = "\t", row.names = FALSE, quote = FALSE,append = FALSE)

# Read updated blast results, sort by primer name and export to fasta
amplicon_table <- readxl::read_excel("RGA_all_family_blast_results_unique_amplicons_27_02_2017.xlsx", sheet = "RGA_all_family_blast_results_un") %>% arrange(Fw_primer) %>% select(Primer=Fw_primer, Amplicon_seq)
write.fasta(amplicon_table, filename = filedate("amplicon_sequences", ".fasta"))
