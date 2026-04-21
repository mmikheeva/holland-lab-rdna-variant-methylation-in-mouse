# holland-lab-rdna-variant-methylation-in-mouse

An R script developed for an ongoing project. This script aims to extract the methylation status of the A/C variant in mouse rDNA.

Input: Reads aligned to the mouse "modified" rDNA sequence in BAM format. All files end with "_1_val_1_bismark_bt2_pe.bam" pattern.

Output: Counts of methylated/unmethylated A/C rDNA variants. Methylated/unmethylated CpG is allocated at -133, A/C SNV is allocated at -104.

Please, check path_to_dir (path to input) before running.
