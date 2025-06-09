library(TwoSampleMR)
bmi2_file <- system.file("extdata/bmi.csv",package="TwoSampleMR")
bmi_exp_dat <- read_exposure_data(
  filename = bmi2_file,
  sep = ",",
  snp_col = "rsid",
  beta_col = "effect",
  se_col = "SE",
  effect_allele_col = "a1",
  other_allele_col = "a2",
  eaf_col = "a1_freq",
  pval_col = "p-value",
  units_col = "Units",
  gene_col = "Gene",
  samplesize_col = "n"
)
bmi_exp_dat <- clump_data(bmi_exp_dat,clump_r2=0.01,pop = "EUR")
chd_out_dat <- extract_outcome_data(
  snps = bmi_exp_dat$SNP,
  outcomes = 'ieu-a-7'
)
dat <- harmonise_data(
  exposure_dat = bmi_exp_dat, 
  outcome_dat = chd_out_dat
)
res <- mr(dat)
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
