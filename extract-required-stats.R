library(data.table)
library(dplyr)
library(tidyr)
library(Rsamtools)
library(GenomicAlignments)
library(stringr)

path_to_dir <- "./bam_BK000964.3.looped_3008"
pattern_to_rm <- "_1_val_1_bismark_bt2_pe.bam"

setwd(path_to_dir)

list_files <- list.files(pattern = ".bam$")

stat_all <- tibble(methylA = as.integer(), unmethylA = as.integer(), methylC = as.integer(), unmethylC = as.integer(), read_pairs = as.integer(), sample = as.character())
for (i in 1:length(list_files)){

###
### Load rDNA reads  
###
  
    pick_file <- list_files[i]
    pick_sample <- gsub(pattern_to_rm, "", pick_file)

    print(paste0("Processing sample ", pick_sample, "."))
    path_file <- paste0(path_to_dir, "/", pick_file)
    
    bam <- BamFile(path_file)
    aln <- readGAlignments(bam, use.names = FALSE, param = ScanBamParam(
        what = c("qname", "mapq", "seq"), tag = c("XM"),
        flag = scanBamFlag(
            isUnmappedQuery = FALSE,
            isSecondaryAlignment = FALSE,
            isDuplicate = FALSE,
            isNotPassingQualityControls = FALSE)))
    aln <- as_tibble(aln)
  
###
### Connect forward/reverse reads under the same qname 
###
  
    aln_mod <- tibble(qname = aln %>% pull(qname) %>% unique()) %>% 
        left_join(aln %>% filter(strand == "+") %>% 
            dplyr::select(qname, start, end, cigar, seq, XM) %>% 
            dplyr::rename(s_plus = start, e_plus = end, cigar_plus = cigar, seq_plus = seq, xm_plus = XM), by = "qname") %>%
        left_join(aln %>% filter(strand == "-") %>% 
            dplyr::select(qname, start, end, cigar, seq, XM) %>% 
            dplyr::rename(s_minus = start, e_minus = end, cigar_minus = cigar, seq_minus = seq, xm_minus = XM), by = "qname") 

###
### Change -133/-104 coordinates to modified coordinates (as rDNA was modified for alignment)
###
  
    ref_pos <- c(-133, -132, -104) # 2876 & 2905
    ref_pos <- ifelse(ref_pos >= 1, ref_pos + 3008, ref_pos + 1 + 3008)   
        min_pos <- min(ref_pos)
        max_pos <- max(ref_pos)   

###
### Pick only reads covering both -133 (methylation) and -104 (SNV) coordinates, modify wrt CIGARs
###
  
  aln_mod <- aln_mod %>% filter(pmin(s_plus, s_minus) <= min_pos & pmax(e_plus, e_minus) >= max_pos)
    
    aln_mod <- aln_mod %>% mutate(
        seq_plus = as.character(sequenceLayer(DNAStringSet(aln_mod$seq_plus), aln_mod$cigar_plus)),
        seq_minus = as.character(sequenceLayer(DNAStringSet(aln_mod$seq_minus), aln_mod$cigar_minus)),
        xm_plus = as.character(sequenceLayer(BStringSet(xm_plus), aln_mod$cigar_plus)),
        xm_minus = as.character(sequenceLayer(BStringSet(xm_minus), aln_mod$cigar_minus))) %>%
        dplyr::select(-c(cigar_plus, cigar_minus))
        
if (nrow(aln_mod) != 0){
aln_merge <- tibble(id = 1:nrow(aln_mod))
for (j in 1:length(ref_pos)){

###
### Check that +/- nucleotides at -133 & -104 are the same
###
  
    aln_pos <- aln_mod %>% 
        mutate(r_plus = ref_pos[j] - s_plus + 1, r_minus = ref_pos[j] - s_minus + 1) %>%
        mutate(r_plus = ifelse(r_plus <= 0, NA, r_plus), r_minus = ifelse(r_minus <= 0, NA, r_minus)) %>%
        dplyr::select(-c(s_plus, s_minus, e_plus, e_minus))

    aln_pos <- aln_pos %>% 
        mutate(b_plus = str_sub(seq_plus, r_plus, r_plus)) %>% 
        mutate(b_minus = str_sub(seq_minus, r_minus, r_minus)) %>% 
        mutate(t_plus = str_sub(xm_plus, r_plus, r_plus)) %>% 
        mutate(t_minus = str_sub(xm_minus, r_minus, r_minus)) %>% 
        dplyr::select(-c(seq_plus, xm_plus, seq_minus, xm_minus, r_plus, r_minus))

    aln_pos <- aln_pos %>% mutate(
        b_plus = paste0(b_plus, t_plus), b_plus = ifelse(b_plus == "NANA", NA, b_plus),
        b_minus = paste0(b_minus, t_minus), b_minus = ifelse(b_minus == "NANA", NA, b_minus)) %>%
        dplyr::select(c(b_plus, b_minus))
        
    aln_pos <- aln_pos %>% mutate(
        b = NA, 
        b = ifelse(!is.na(b_plus) & is.na(b_minus), b_plus, b),
        b = ifelse(is.na(b_plus) & !is.na(b_minus), b_minus, b),
        b = ifelse(!is.na(b_plus) & !is.na(b_minus) & b_plus != b_minus, "-", b),
        b = ifelse(!is.na(b_plus) & !is.na(b_minus) & b_plus == b_minus, b_plus, b)) %>% 
        dplyr::select(-c(b_plus, b_minus))
    aln_merge <- aln_merge %>% cbind(
        aln_pos %>% dplyr::rename(!!paste0("b", ref_pos[j]) := b)) %>%
        as_tibble()}

###
### Overlap with Bismark methylation tag
###
  
# -133:CZ or -132:GZ methyl
# -133:Tz or -133:Az unmethyl
        
stat_sample <- aln_merge %>% mutate(Status = NA, 
    Status = ifelse(b2876 == "CZ" | b2877 == "GZ", "methyl", Status),
    Status = ifelse(b2876 == "Tz" | b2877 == "Az", "unmethyl", Status))
stat_sample <- stat_sample %>% mutate(Variant = NA,
    Variant = ifelse(substr(b2905,1,1) == "A", "A", Variant),    
    Variant = ifelse(substr(b2905,1,1) == "T", "T", Variant),    
    Variant = ifelse(substr(b2905,1,1) == "C", "C", Variant))      
stat_sample <- stat_sample %>% filter(!is.na(Status) & !is.na(Variant)) %>% summarize(
    methylA = sum(Status == "methyl" & Variant == "A"),
    unmethylA = sum(Status == "unmethyl" & Variant == "A"),
    methylC = sum(Status == "methyl" & Variant == "C"),
    unmethylC = sum(Status == "unmethyl" & Variant == "C"),
    methylT = sum(Status == "methyl" & Variant == "T"),
    unmethylT = sum(Status == "unmethyl" & Variant == "T"),
    read_pairs = nrow(stat_sample), sample = pick_sample)
} else {
stat_sample <- tibble(
    methylA = 0, unmethylA = 0, methylC = 0, unmethylC = 0, methylT = 0, unmethylT = 0,
    read_pairs = 0, sample = pick_sample)}

    stat_all <- stat_all %>% rbind(stat_sample) %>% as_tibble()}

###
### Report
###

write.csv(stat_all, "./variant-methylation-per-read-pair-level.csv")

