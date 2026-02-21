
############################################################
# Reproducibility Script for LassoHiDFastGibbs
#
# IMPORTANT NOTE:
#
# Running all R Markdown files in this script with the full
# MCMC settings (e.g., 11000 iterations per chain) may take
# several days on a standard desktop machine.
#
# The long runtime is due to:
#   - High-dimensional MCMC sampling
#   - Large benchmark datasets
#   - Multiple competing methods fitted for comparison
#
############################################################


results_dir <- "results"
stopifnot(dir.exists(results_dir))

datasets_ngtp <- c("BostonHousing2", "diabetes2", "Hitters2", "Crime", "energy2", "Kakadu2")

#####################
# ----- Running models with horseshoe prior and datasets n>p -------
#####################

for (d in datasets_ngtp) {
  print(d)
  rmarkdown::render(
    "benchmark_horseshoe_ngtp.Rmd",
    params = list(dataset_name = d, results_dir = "results"),
    envir = new.env(parent = globalenv())
  )
}

fun = median
load(file.path(results_dir,"BostonHousing2_horseshoe_results_ngtp.Rdata"))

res_brms = res_brms[-4]

stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_brms      <- res_brms
stat_horseshoe <- apply(res_horseshoe$mStat,2,fun)


table_bh2 = rbind(
  stat_4BG,
  stat_3BG,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayeslm,
  stat_bayesreg,
  stat_brms,
  stat_horseshoe)



load(file.path(results_dir,"Crime_horseshoe_results_ngtp.Rdata"))

res_brms = res_brms[-4]

stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_brms      <- res_brms
stat_horseshoe <- apply(res_horseshoe$mStat,2,fun)


table_crime = rbind(
  stat_4BG,
  stat_3BG,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayeslm,
  stat_bayesreg,
  stat_brms,
  stat_horseshoe)


load(file.path(results_dir,"diabetes2_horseshoe_results_ngtp.Rdata"))

res_brms = res_brms[-4]

stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_brms      <- res_brms
stat_horseshoe <- apply(res_horseshoe$mStat,2,fun)


table_diabetes2 = rbind(
  stat_4BG,
  stat_3BG,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayeslm,
  stat_bayesreg,
  stat_brms,
  stat_horseshoe)


load(file.path(results_dir,"energy2_horseshoe_results_ngtp.Rdata"))
res_brms = res_brms[-4]

stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_brms      <- res_brms
stat_horseshoe <- apply(res_horseshoe$mStat,2,fun)


table_energy2 = rbind(
  stat_4BG,
  stat_3BG,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayeslm,
  stat_bayesreg,
  stat_brms,
  stat_horseshoe)

load(file.path(results_dir,"Hitters2_horseshoe_results_ngtp.Rdata"))

res_brms = res_brms[-4]

stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_brms      <- res_brms
stat_horseshoe <- apply(res_horseshoe$mStat,2,fun)


table_Hitters2 = rbind(
  stat_4BG,
  stat_3BG,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayeslm,
  stat_bayesreg,
  stat_brms,
  stat_horseshoe)



load(file.path(results_dir,"Kakadu2_horseshoe_results_ngtp.Rdata"))

res_brms = res_brms[-4]

stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_brms      <- res_brms
stat_horseshoe <- apply(res_horseshoe$mStat,2,fun)


table_Kakadu2 = rbind(
  stat_4BG,
  stat_3BG,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayeslm,
  stat_bayesreg,
  stat_brms,
  stat_horseshoe)


rel_val = 3

reff_crime = table_crime[,6]/table_crime[rel_val,6]
rank_crime = rank(-reff_crime)
reff_crime = round(reff_crime,2)

reff_diabetes2 = table_diabetes2[,6]/table_diabetes2[rel_val,6]
rank_diabetes2 = rank(-reff_diabetes2)
reff_diabetes2 = round(reff_diabetes2,2)

reff_hitters2 = table_Hitters2[,6]/table_Hitters2[rel_val,6]
rank_hitters2 = rank(-reff_hitters2)
reff_hitters2 = round(reff_hitters2,2)

reff_bh2 = table_bh2[,6]/table_bh2[rel_val,6]
rank_bh2 = rank(-reff_bh2)
reff_bh2 = round(reff_bh2,2)

reff_Kakadu2 = table_Kakadu2[,6]/table_Kakadu2[rel_val,6]
rank_Kakadu2 = rank(-reff_Kakadu2)
reff_Kakadu2 = round(reff_Kakadu2,2)

reff_energy2 = table_energy2[,6]/table_energy2[rel_val,6]
rank_energy2 = rank(-reff_energy2)
reff_energy2 = round(reff_energy2,2)

ave_rank = round( (rank_crime + rank_diabetes2 + rank_hitters2 + rank_bh2 + rank_Kakadu2 + rank_energy2)/6, 1)

table3  = cbind(reff_bh2, reff_crime, reff_diabetes2, reff_energy2, reff_hitters2, reff_Kakadu2,
                rank_bh2, rank_crime, rank_diabetes2, rank_energy2, rank_hitters2, rank_Kakadu2,  ave_rank)

table3

table4 = rbind(
  table_bh2,
  table_crime,
  table_diabetes2,
  table_energy2,
  table_Hitters2,
  table_Kakadu2
)

table4[table4[,1]>100,1] = 100
table4[table4[,3]>100,3] = 100
table4[table4[,5]>100,5] = 100

round(table4,1)


#####################
# ---- Running models with Lasso prior and datasets n>p -----
#####################

for (d in datasets_ngtp) {
  print(d)
  rmarkdown::render(
    "benchmark_lasso_ngtp.Rmd",
    params = list(dataset_name = d, results_dir = "results"),
    envir = new.env(parent = globalenv())
  )
}


load(file.path(results_dir,"Hitters2_results_ngtp.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_bl    <- apply(res_2BG_bl$mStat,2,fun)
stat_2BG_bs    <- apply(res_2BG_bs$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_la    <- apply(res_pcg_la$mStat,2,fun)
stat_pcg_sa    <- apply(res_pcg_sa$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_monomvn   <- apply(res_monomvn$mStat,2,fun)
stat_rstan     <- apply(res_rstan$mStat,2,fun)

table_Hitters2 = rbind(
  stat_hans,
  stat_4BG,
  stat_3BG,
  stat_2BG_bl,
  stat_2BG_bs,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_la,
  stat_pcg_sa,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayeslm,
  stat_bayesreg,
  stat_monomvn,
  stat_rstan
)


load(file.path(results_dir,"diabetes2_results_ngtp.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_bl    <- apply(res_2BG_bl$mStat,2,fun)
stat_2BG_bs    <- apply(res_2BG_bs$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_la    <- apply(res_pcg_la$mStat,2,fun)
stat_pcg_sa    <- apply(res_pcg_sa$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_monomvn   <- apply(res_monomvn$mStat,2,fun)
stat_rstan     <- apply(res_rstan$mStat,2,fun)

table_diabetes2 = rbind(
  stat_hans,
  stat_4BG,
  stat_3BG,
  stat_2BG_bl,
  stat_2BG_bs,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_la,
  stat_pcg_sa,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayeslm,
  stat_bayesreg,
  stat_monomvn,
  stat_rstan
)


load(file.path(results_dir,"Crime_results_ngtp.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_bl    <- apply(res_2BG_bl$mStat,2,fun)
stat_2BG_bs    <- apply(res_2BG_bs$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_la    <- apply(res_pcg_la$mStat,2,fun)
stat_pcg_sa    <- apply(res_pcg_sa$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_monomvn   <- apply(res_monomvn$mStat,2,fun)
stat_rstan     <- apply(res_rstan$mStat,2,fun)

table_crime = rbind(
  stat_hans,
  stat_4BG,
  stat_3BG,
  stat_2BG_bl,
  stat_2BG_bs,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_la,
  stat_pcg_sa,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayeslm,
  stat_bayesreg,
  stat_monomvn,
  stat_rstan
)

load(file.path(results_dir,"BostonHousing2_results_ngtp.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_bl    <- apply(res_2BG_bl$mStat,2,fun)
stat_2BG_bs    <- apply(res_2BG_bs$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_la    <- apply(res_pcg_la$mStat,2,fun)
stat_pcg_sa    <- apply(res_pcg_sa$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_monomvn   <- apply(res_monomvn$mStat,2,fun)
stat_rstan     <- apply(res_rstan$mStat,2,fun)

table_bh2 = rbind(
  stat_hans,
  stat_4BG,
  stat_3BG,
  stat_2BG_bl,
  stat_2BG_bs,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_la,
  stat_pcg_sa,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayeslm,
  stat_bayesreg,
  stat_monomvn,
  stat_rstan
)


load(file.path(results_dir,"energy2_results_ngtp.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_bl    <- apply(res_2BG_bl$mStat,2,fun)
stat_2BG_bs    <- apply(res_2BG_bs$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_la    <- apply(res_pcg_la$mStat,2,fun)
stat_pcg_sa    <- apply(res_pcg_sa$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_monomvn   <- apply(res_monomvn$mStat,2,fun)
stat_rstan     <- apply(res_rstan$mStat,2,fun)

table_energy2 = rbind(
  stat_hans,
  stat_4BG,
  stat_3BG,
  stat_2BG_bl,
  stat_2BG_bs,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_la,
  stat_pcg_sa,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayeslm,
  stat_bayesreg,
  stat_monomvn,
  stat_rstan
)


load(file.path(results_dir,"Kakadu2_results_ngtp.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_bl    <- apply(res_2BG_bl$mStat,2,fun)
stat_2BG_bs    <- apply(res_2BG_bs$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_la    <- apply(res_pcg_la$mStat,2,fun)
stat_pcg_sa    <- apply(res_pcg_sa$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayeslm   <- apply(res_bayeslm$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)
stat_monomvn   <- apply(res_monomvn$mStat,2,fun)
stat_rstan     <- apply(res_rstan$mStat,2,fun)

table_Kakadu2 = rbind(
  stat_hans,
  stat_4BG,
  stat_3BG,
  stat_2BG_bl,
  stat_2BG_bs,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_la,
  stat_pcg_sa,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayeslm,
  stat_bayesreg,
  stat_monomvn,
  stat_rstan
)


reff_crime = table_crime[,6]/table_crime[6,6]
rank_crime = rank(-reff_crime)
reff_crime = round(reff_crime,2)

reff_diabetes2 = table_diabetes2[,6]/table_diabetes2[6,6]
rank_diabetes2 = rank(-reff_diabetes2)
reff_diabetes2 = round(reff_diabetes2,2)

reff_hitters2 = table_Hitters2[,6]/table_Hitters2[6,6]
rank_hitters2 = rank(-reff_hitters2)
reff_hitters2 = round(reff_hitters2,2)

reff_bh2 = table_bh2[,6]/table_bh2[6,6]
rank_bh2 = rank(-reff_bh2)
reff_bh2 = round(reff_bh2,2)

reff_Kakadu2 = table_Kakadu2[,6]/table_Kakadu2[6,6]
rank_Kakadu2 = rank(-reff_Kakadu2)
reff_Kakadu2 = round(reff_Kakadu2,2)

reff_energy2 = table_energy2[,6]/table_energy2[6,6]
rank_energy2 = rank(-reff_energy2)
reff_energy2 = round(reff_energy2,2)

ave_rank = round( (rank_crime + rank_diabetes2 + rank_hitters2 + rank_bh2 + rank_Kakadu2 + rank_energy2)/6, 1)

table1  = cbind(reff_bh2, reff_crime, reff_diabetes2, reff_energy2, reff_hitters2, reff_Kakadu2,
                rank_bh2, rank_crime, rank_diabetes2, rank_energy2, rank_hitters2, rank_Kakadu2,  ave_rank)

table1


table2 = rbind(
  table_bh2,
  table_crime,
  table_diabetes2,
  table_energy2,
  table_Hitters2,
  table_Kakadu2
)


table2[table2[,1]>100,1] = 100
table2[table2[,3]>100,3] = 100
table2[table2[,5]>100,5] = 100

round(table2,1)



################################################################################
################################################################################
#  Running the Rmd files for datasets with p>n
################################################################################
################################################################################



datasets_pgtn <- c("eyedata", "cookie", "riboflavin", "covid", "qtl")

#####################
# ---- Running models with horseshoe prior and datasets p>n -----
#####################

for (d in datasets_pgtn) {
  print(d)
  rmarkdown::render(
    "benchmark_horseshoe_pgtn.Rmd",
    params = list(dataset_name = d, results_dir = "results"),
    envir = new.env(parent = globalenv())
  )
}


load(file.path(results_dir,"eyedata_horseshoe_results_pgtn.Rdata"))

stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng2       <- apply(res_ng2$mStat,2,fun)
stat_ng3       <- apply(res_ng3$mStat,2,fun)
stat_ng4       <- apply(res_ng4$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_EHS       <- apply(res_EHS$mStat,2,fun)
stat_AHS       <- apply(res_AHS$mStat,2,fun)

table_eye <- rbind(
  stat_4BG,
  stat_3BG,
  stat_2BG_31,
  stat_ng2,
  stat_ng3,
  stat_ng4,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_EHS,
  stat_AHS
)




load(file.path(results_dir,"cookie_horseshoe_results_pgtn.Rdata"))

stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng2       <- apply(res_ng2$mStat,2,fun)
stat_ng3       <- apply(res_ng3$mStat,2,fun)
stat_ng4       <- apply(res_ng4$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_EHS       <- apply(res_EHS$mStat,2,fun)
stat_AHS       <- apply(res_AHS$mStat,2,fun)

table_cookie <- rbind(
  stat_4BG,
  stat_3BG,
  stat_2BG_31,
  stat_ng2,
  stat_ng3,
  stat_ng4,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_EHS,
  stat_AHS
)



load(file.path(results_dir,"riboflavin_horseshoe_results_pgtn.Rdata"))

stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng2       <- apply(res_ng2$mStat,2,fun)
stat_ng3       <- apply(res_ng3$mStat,2,fun)
stat_ng4       <- apply(res_ng4$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_EHS       <- apply(res_EHS$mStat,2,fun)
stat_AHS       <- apply(res_AHS$mStat,2,fun)

table_ribo <- rbind(
  stat_4BG,
  stat_3BG,
  stat_2BG_31,
  stat_ng2,
  stat_ng3,
  stat_ng4,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_EHS,
  stat_AHS
)


load(file.path(results_dir,"covid_horseshoe_results_pgtn.Rdata"))

stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng2       <- apply(res_ng2$mStat,2,fun)
stat_ng3       <- apply(res_ng3$mStat,2,fun)
stat_ng4       <- apply(res_ng4$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_EHS       <- apply(res_EHS$mStat,2,fun)
stat_AHS       <- apply(res_AHS$mStat,2,fun)

table_covid <- rbind(
  stat_4BG,
  stat_3BG,
  stat_2BG_31,
  stat_ng2,
  stat_ng3,
  stat_ng4,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_EHS,
  stat_AHS
)


load(file.path(results_dir,"qtl_horseshoe_results_pgtn.Rdata"))

stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng2       <- apply(res_ng2$mStat,2,fun)
stat_ng3       <- apply(res_ng3$mStat,2,fun)
stat_ng4       <- apply(res_ng4$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_EHS       <- apply(res_EHS$mStat,2,fun)
stat_AHS       <- apply(res_AHS$mStat,2,fun)

table_qtl <- rbind(
  stat_4BG,
  stat_3BG,
  stat_2BG_31,
  stat_ng2,
  stat_ng3,
  stat_ng4,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_EHS,
  stat_AHS
)


val = 3

reff_eye = table_eye[,6]/table_eye[val,6]
rank_eye = rank(-reff_eye)
reff_eye = round(reff_eye,2)

reff_cookie = table_cookie[,6]/table_cookie[val,6]
rank_cookie = rank(-reff_cookie)
reff_cookie = round(reff_cookie,2)

reff_ribo = table_ribo[,6]/table_ribo[val,6]
rank_ribo = rank(-reff_ribo)
reff_ribo = round(reff_ribo,2)

reff_covid = table_covid[,6]/table_covid[val,6]
rank_covid = rank(-reff_covid)
reff_covid = round(reff_covid,2)

reff_qtl = table_qtl[,6]/table_qtl[val,6]
rank_qtl = rank(-reff_qtl)
reff_qtl = round(reff_qtl,2)

ave_rank = round( (rank_eye + rank_cookie + rank_ribo + rank_covid + rank_qtl)/5, 1)

table7  = cbind(reff_eye, reff_cookie, reff_ribo, reff_covid, reff_qtl,
                rank_eye, rank_cookie, rank_ribo, rank_covid, rank_qtl,  ave_rank)

table7









#####################
# ---- Running models with Lasso prior and datasets p>n -----
#####################
for (d in datasets_pgtn) {
  print(d)
  rmarkdown::render(
    "benchmark_lasso_pgtn.Rmd",
    params = list(dataset_name = d, results_dir = "results"),
    envir = new.env(parent = globalenv())
  )
}



load(file.path(results_dir,"eyedata_results_pgtn.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_bl    <- apply(res_2BG_bl$mStat,2,fun)
stat_2BG_bs    <- apply(res_2BG_bs$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_la    <- apply(res_pcg_la$mStat,2,fun)
stat_pcg_sa    <- apply(res_pcg_sa$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)


table_eye = rbind(
  stat_hans,
  stat_4BG,
  stat_3BG,
  stat_2BG_bl,
  stat_2BG_bs,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_la,
  stat_pcg_sa,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayesreg

)



load(file.path(results_dir,"cookie_results_pgtn.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_bl    <- apply(res_2BG_bl$mStat,2,fun)
stat_2BG_bs    <- apply(res_2BG_bs$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_la    <- apply(res_pcg_la$mStat,2,fun)
stat_pcg_sa    <- apply(res_pcg_sa$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)


table_cookie = rbind(
  stat_hans,
  stat_4BG,
  stat_3BG,
  stat_2BG_bl,
  stat_2BG_bs,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_la,
  stat_pcg_sa,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayesreg
)


load(file.path(results_dir,"riboflavin_results_pgtn.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_bl    <- apply(res_2BG_bl$mStat,2,fun)
stat_2BG_bs    <- apply(res_2BG_bs$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_la    <- apply(res_pcg_la$mStat,2,fun)
stat_pcg_sa    <- apply(res_pcg_sa$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)


table_ribo = rbind(
  stat_hans,
  stat_4BG,
  stat_3BG,
  stat_2BG_bl,
  stat_2BG_bs,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_la,
  stat_pcg_sa,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayesreg
)


load(file.path(results_dir,"covid_results_pgtn.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_bl    <- apply(res_2BG_bl$mStat,2,fun)
stat_2BG_bs    <- apply(res_2BG_bs$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_la    <- apply(res_pcg_la$mStat,2,fun)
stat_pcg_sa    <- apply(res_pcg_sa$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)


table_covid = rbind(
  stat_hans,
  stat_4BG,
  stat_3BG,
  stat_2BG_bl,
  stat_2BG_bs,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_la,
  stat_pcg_sa,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayesreg

)


load(file.path(results_dir,"qtl_results_pgtn.Rdata"))

stat_hans      <- apply(res_hans$mStat,2,fun)
stat_4BG       <- apply(res_4BG$mStat,2,fun)
stat_3BG       <- apply(res_3BG$mStat,2,fun)
stat_2BG_bl    <- apply(res_2BG_bl$mStat,2,fun)
stat_2BG_bs    <- apply(res_2BG_bs$mStat,2,fun)
stat_2BG_31    <- apply(res_ng1$mStat,2,fun)
stat_ng5       <- apply(res_ng5$mStat,2,fun)
stat_ng10      <- apply(res_ng10$mStat,2,fun)
stat_pcg_ls    <- apply(res_pcg_ls$mStat,2,fun)
stat_pcg_sl    <- apply(res_pcg_sl$mStat,2,fun)
stat_pcg_la    <- apply(res_pcg_la$mStat,2,fun)
stat_pcg_sa    <- apply(res_pcg_sa$mStat,2,fun)
stat_pcg_bs    <- apply(res_pcg_bs$mStat,2,fun)
stat_pcg_sb    <- apply(res_pcg_sb$mStat,2,fun)
stat_bayesreg  <- apply(res_bayesreg$mStat,2,fun)


table_qtl = rbind(
  stat_hans,
  stat_4BG,
  stat_3BG,
  stat_2BG_bl,
  stat_2BG_bs,
  stat_2BG_31,
  stat_ng5,
  stat_ng10,
  stat_pcg_ls,
  stat_pcg_sl,
  stat_pcg_la,
  stat_pcg_sa,
  stat_pcg_bs,
  stat_pcg_sb,
  stat_bayesreg

)



reff_eye = table_eye[,6]/table_eye[6,6]
rank_eye = rank(-reff_eye)
reff_eye = round(reff_eye,2)

reff_cookie = table_cookie[,6]/table_cookie[6,6]
rank_cookie = rank(-reff_cookie)
reff_cookie = round(reff_cookie,2)

reff_ribo = table_ribo[,6]/table_ribo[6,6]
rank_ribo = rank(-reff_ribo)
reff_ribo = round(reff_ribo,2)

reff_covid = table_covid[,6]/table_covid[6,6]
rank_covid = rank(-reff_covid)
reff_covid = round(reff_covid,2)

reff_qtl = table_qtl[,6]/table_qtl[6,6]
rank_qtl = rank(-reff_qtl)
reff_qtl = round(reff_qtl,2)

ave_rank = round( (rank_eye + rank_cookie + rank_ribo + rank_covid + rank_qtl)/5, 1)

table5  = cbind(reff_eye, reff_cookie, reff_ribo, reff_covid, reff_qtl,
                rank_eye, rank_cookie, rank_ribo, rank_covid, rank_qtl,  ave_rank)

table5


