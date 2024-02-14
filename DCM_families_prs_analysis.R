# Analysis of DCM-PRS in DCM families
# Authors: Ingrid Tarr and Jamie-Lee Thompson

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(data.table)
library(magrittr)
library(ggridges)
library(GMMAT)

# Load score data
score_fl <- list.files(pattern = "dcm_families_prs.csv", full.names = T)
scores <- read.csv(score_fl, header = T)

# Create standardised scores
all_scores <- mutate(scores, stSCORESUM = (SCORESUM - mean(SCORESUM))/sd(SCORESUM))

# Load relatedness matrix, check IDs
if (!file.exists("rel_mat.rel.id_edited")){
  grm_samples <- read.delim("rel_mat.rel.id", header = F)
  grep(":", grm_samples$V1) # find and deal with any duplicate samples
  grm_samples[grep(":", grm_samples$V1), ] <- gsub("2:", "", grm_samples[grep(":", grm_samples$V1), ])
  write.table(grm_samples, "rel_mat.rel.id_edited", quote = F, col.names = F, row.names = F, sep = "\t")
}
#
grm_samples <- read.delim("rel_mat.rel.id_edited", header = F)
grm <- read.delim("rel_mat.rel", header = F)
grm_mat <- as.matrix(grm)
colnames(grm_mat) <- grm_samples$V1
rownames(grm_mat) <- grm_samples$V1

dummy_infile_cols <- data.frame(SNP = "grs")

# Perform modelling per pairwise group comparison: controls vs. affected patients
all_scores_control_vs_affected <- all_scores %>% 
  filter(affected_relative != "n")
# Check numbers
group_by(all_scores_control_affected, affected) %>% tally 
# Create file of scores for modelling
all_scores_control_affected %>% 
  select(IID, stSCORESUM) %>% 
  pivot_wider(names_from = IID, values_from = stSCORESUM) %>% 
  cbind(dummy_infile_cols, .) %>% 
  write.table(., file = "control_vs_affected_prs_file_infile_for_wald.txt", row.names = F, sep = "\t", quote = F)

# Run model
wald_control_vs_affected <- glmm.wald(as.factor(affected) ~ 1, data = all_scores_control_affected, id = "NCI_ID",
                                      kins = grm_mat, family = binomial(link = "logit"),
                                      infile = "control_vs_affected_prs_file_infile_for_wald.txt", snps = "grs",
                                      infile.nrow = 2, infile.nrow.skip = 1, snp.col = 1,infile.ncol.skip = 1)
wald_control_vs_affected
exp(wald_control_vs_affected$BETA)

# Perform modelling per pairwise group comparison: controls vs. unaffected patients
all_scores_control_vs_unaffected <- all_scores %>% 
  filter(affected_relative != "y")
# Check numbers
group_by(all_scores_control_unaffected, affected) %>% tally 
# Create file of scores for modelling
all_scores_control_unaffected %>% 
  select(IID, stSCORESUM) %>% 
  pivot_wider(names_from = IID, values_from = stSCORESUM) %>% 
  cbind(dummy_infile_cols, .) %>% 
  write.table(., file = "control_vs_unaffected_prs_file_infile_for_wald.txt", row.names = F, sep = "\t", quote = F)

# Run model
wald_control_vs_unaffected <- glmm.wald(as.factor(affected) ~ 1, data = all_scores_control_unaffected, id = "NCI_ID",
                                      kins = grm_mat, family = binomial(link = "logit"),
                                      infile = "control_vs_unaffected_prs_file_infile_for_wald.txt", snps = "grs",
                                      infile.nrow = 2, infile.nrow.skip = 1, snp.col = 1,infile.ncol.skip = 1)
wald_control_vs_unaffected
exp(wald_control_vs_unaffected$BETA)

# Perform modelling per pairwise group comparison: unaffected vs. affected patients
all_scores_unaffected_vs_affected <- all_scores %>% 
  filter(control != "y")
# Check numbers
group_by(all_scores_unaffected_affected, affected) %>% tally 
# Create file of scores for modelling
all_scores_unaffected_affected %>% 
  select(IID, stSCORESUM) %>% 
  pivot_wider(names_from = IID, values_from = stSCORESUM) %>% 
  cbind(dummy_infile_cols, .) %>% 
  write.table(., file = "unaffected_vs_affected_prs_file_infile_for_wald.txt", row.names = F, sep = "\t", quote = F)

# Run model
wald_unaffected_vs_affected <- glmm.wald(as.factor(affected) ~ 1, data = all_scores_unaffected_affected, id = "NCI_ID",
                                      kins = grm_mat, family = binomial(link = "logit"),
                                      infile = "unaffected_vs_affected_prs_file_infile_for_wald.txt", snps = "grs",
                                      infile.nrow = 2, infile.nrow.skip = 1, snp.col = 1,infile.ncol.skip = 1)
wald_unaffected_vs_affected
exp(wald_unaffected_vs_affected$BETA)


# Benjamini and Hochberg correction for multiple testing 
p.adjust(c(wald_control_vs_affected$PVAL, wald_control_vs_unaffected$PVAL, wald_unaffected_vs_affected$PVAL), method = "BH")

# Plot results
all_scores_plot <- all_scores %>% 
  mutate(affected = ifelse(affected_relative == "c", "Controls\nn=1127", 
                           ifelse(affected_relative == "y", "Affected\nn=176", "Unaffected\nn=153")
  ))

p <- ggplot(all_scores_sub_plot, aes(x = stSCORESUM, fill = affected, y = affected)) + 
  geom_density_ridges(alpha = 1, bandwidth = 0.18) +
  guides(fill = "none") + 
  scale_fill_manual(values = col_bli[c(1:3)]) +
  theme_minimal() +
  labs(y = "", x = "Standardised DCM PRS (units of SD)")
p + theme(text = element_text(size = 19))

ggsave("DCM_families.pdf", width=6, height=5)
