## scripts to assess relationship of CyTOF phenotype in T1D to C-peptide levels and age

##### set up environment: load/save data objects #####

# rm(list=ls())
# opar <- par()


##### set up environment: load packages #####

## load general packages
library(tidyverse)
theme_set(
  theme_bw(20) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, size=1),
          axis.text=element_text(colour="black"),
          axis.ticks=element_line(colour="black")))
update_geom_defaults("point", list(shape=16))

## load analysis-specific packages
library(lmerTest)

# load Matt Dufort's useful packages(mjdufort/miscHelpers on gitHub)
library(miscHelpers)


##### load data from other sources (usually only do this when initially running analyses) #####

# load Tmr+ CyTOF data
data.CyTOF <-
  readxl::read_xlsx(
    "allSubjectsSummary_Compiled_10K_Samples_3_Metas_with_Age_Rapid_vs_Slow_Plus_Total_CD8_180918.xlsx",
    sheet=1, skip=1) %>%
  standardize_dimnames() %>%
  dplyr::rename(subject=id) %>%
  as.data.frame() %>%
  remove_all_NA_rowcols()

## alternatively, load total CD8 CyTOF data
# data.CyTOF <-
#   readxl::read_xlsx(
#     "allSubjectsSummary_Compiled_10K_Samples_3_Metas_with_Age_Rapid_vs_Slow_Plus_Total_CD8_180918.xlsx",
#     sheet=2, skip=1) %>%
#   standardize_dimnames() %>%
#   dplyr::rename(subject=id) %>%
#   as.data.frame() %>%
#   remove_all_NA_rowcols()
 
# extract meta_1_perc, meta_11_perc, and cpeptide_group from the columns on the right
pre_data_columns <- 2
colnames(data.CyTOF) <-
  colnames(data.CyTOF) %>%
  str_replace_all("slow_1", "slow_too_short") %>%
  str_replace_all("slow_2", "rapid_too_long")
data.CyTOF$cpeptide_group <-
  factor(NA, levels=c("slow", "rapid", "slow_too_short", "rapid_too_long"))
data.CyTOF[,c("meta_1_perc", "meta_11_perc", "meta_12_perc")] <- as.numeric(NA)
for (i in 1:nrow(data.CyTOF)) {
  # skip ones that have no value in any of the % columns
  if (all(is.na(data.CyTOF[i, str_detect(colnames(data.CyTOF), "meta")])))
    next
  data.CyTOF$cpeptide_group[i] <- # pull the cpeptide group by the non-NA meta_1 column
    str_replace(
      colnames(data.CyTOF)[
        which(!is.na(data.CyTOF[i, pre_data_columns+(1:12)]))[1]+pre_data_columns],
      "^meta_1_", "")
  
  data.CyTOF$meta_1_perc[i] <- # pull the % for meta_1
    data.CyTOF[i, paste0("meta_1_", data.CyTOF$cpeptide_group[i])]
  data.CyTOF$meta_11_perc[i] <- # pull the % for meta_11
    data.CyTOF[i, paste0("meta_11_", data.CyTOF$cpeptide_group[i])]
  data.CyTOF$meta_12_perc[i] <- # pull the % for meta_12
    data.CyTOF[i, paste0("meta_12_", data.CyTOF$cpeptide_group[i])]
}

data.CyTOF <-
  data.CyTOF %>%
  dplyr::select(
    -contains("rapid"), -contains("slow"))

glimpse(data.CyTOF)


##### plot meta_1 and meta_11 by age and cpeptide_group #####

pdf(paste0("plots/meta_1_percent_vs_age_at_draw.by_cpeptide_group.", data_date, ".pdf"),
    w=9, h=6)
data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  ggplot(
    mapping=aes(
      x=age_years,
      y=meta_1_perc,
      color=cpeptide_group)) +
  geom_point(size=3, shape=16) +
  geom_smooth(method="lm", se=FALSE) +
  scale_color_manual(
    "C-peptide\ngroup",
    values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 1 %")
dev.off()

pdf(paste0("plots/meta_11_percent_vs_age_at_draw.by_cpeptide_group.", data_date, ".pdf"),
    w=9, h=6)
data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  ggplot(
    mapping=aes(
      x=age_years,
      y=meta_11_perc,
      color=cpeptide_group)) +
  geom_point(size=3, shape=16) +
  geom_smooth(method="lm", se=FALSE) +
  scale_color_manual(
    "C-peptide\ngroup",
    values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 11 %")
dev.off()

pdf(paste0("plots/meta_12_percent_vs_age_at_draw.by_cpeptide_group.", data_date, ".pdf"),
    w=9, h=6)
data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  ggplot(
    mapping=aes(
      x=age_years,
      y=meta_12_perc,
      color=cpeptide_group)) +
  geom_point(size=3, shape=16) +
  geom_smooth(method="lm", se=FALSE) +
  scale_color_manual(
    "C-peptide\ngroup",
    values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 12 %")
dev.off()


# with log-transformed age at draw
pdf(paste0("plots/meta_1_percent_vs_log_age_at_draw.by_cpeptide_group.", data_date, ".pdf"),
    w=9, h=6)
data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  ggplot(
    mapping=aes(
      x=age_years,
      y=meta_1_perc,
      color=cpeptide_group)) +
  geom_point(size=3, shape=16) +
  geom_smooth(method="lm", se=FALSE) +
  scale_x_log10(breaks=c(10,20,30)) +
  scale_color_manual(
    "C-peptide\ngroup",
    values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 1 %")
dev.off()

pdf(paste0("plots/meta_11_percent_vs_log_age_at_draw.by_cpeptide_group.", data_date, ".pdf"),
    w=9, h=6)
data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  ggplot(
    mapping=aes(
      x=age_years,
      y=meta_11_perc,
      color=cpeptide_group)) +
  geom_point(size=3, shape=16) +
  geom_smooth(method="lm", se=FALSE) +
  scale_x_log10(breaks=c(10,20,30)) +
  scale_color_manual(
    "C-peptide\ngroup",
    values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 11 %")
dev.off()

pdf(paste0("plots/meta_12_percent_vs_log_age_at_draw.by_cpeptide_group.", data_date, ".pdf"),
    w=9, h=6)
data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  ggplot(
    mapping=aes(
      x=age_years,
      y=meta_12_perc,
      color=cpeptide_group)) +
  geom_point(size=3, shape=16) +
  geom_smooth(method="lm", se=FALSE) +
  scale_x_log10(breaks=c(10,20,30)) +
  scale_color_manual(
    "C-peptide\ngroup",
    values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 12 %")
dev.off()


data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  ggplot(
    mapping=aes(
      x=age_years,
      y=meta_1_perc,
      color=cpeptide_group)) +
  geom_point(size=3, shape=16) +
  geom_smooth(method="lm", se=FALSE) +
  ggthemes::scale_color_colorblind("C-peptide\ngroup") +
  # scale_color_manual(
  #   "C-peptide\ngroup",
  #   values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 1 %")

data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  ggplot(
    mapping=aes(
      x=age_years,
      y=meta_11_perc,
      color=cpeptide_group)) +
  geom_point(size=3, shape=16) +
  geom_smooth(method="lm", se=FALSE) +
  ggthemes::scale_color_colorblind("C-peptide\ngroup") +
  # scale_color_manual(
  #   "C-peptide\ngroup",
  #   values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 11 %")

data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  ggplot(
    mapping=aes(
      x=age_years,
      y=meta_12_perc,
      color=cpeptide_group)) +
  geom_point(size=3, shape=16) +
  geom_smooth(method="lm", se=FALSE) +
  ggthemes::scale_color_colorblind("C-peptide\ngroup") +
  # scale_color_manual(
  #   "C-peptide\ngroup",
  #   values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 12 %")


##### fit models of metacluster to C-peptide group, adjusting for age_years #####

## models with interaction
lm.meta_1_perc_vs_cpeptide_group.age_years.with_interaction.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  lm(meta_1_perc ~ cpeptide_group * age_years,
     data=.)
summary(lm.meta_1_perc_vs_cpeptide_group.age_years.with_interaction.slow_rapid_only)
# plot(lm.meta_1_perc_vs_cpeptide_group.age_years.with_interaction.slow_rapid_only)
# interaction term is not significant; proceed without interaction terms

lm.meta_11_perc_vs_cpeptide_group.age_years.with_interaction.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  lm(meta_11_perc ~ cpeptide_group * age_years,
     data=.)
summary(lm.meta_11_perc_vs_cpeptide_group.age_years.with_interaction.slow_rapid_only)
# plot(lm.meta_11_perc_vs_cpeptide_group.age_years.with_interaction.slow_rapid_only)
# interaction term is marginally significant; what to do?

lm.meta_12_perc_vs_cpeptide_group.age_years.with_interaction.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  lm(meta_12_perc ~ cpeptide_group * age_years,
     data=.)
summary(lm.meta_12_perc_vs_cpeptide_group.age_years.with_interaction.slow_rapid_only)
# plot(lm.meta_12_perc_vs_cpeptide_group.age_years.with_interaction.slow_rapid_only)
# interaction term is not significant; proceed without interaction terms


## models without interaction term
lm.meta_1_perc_vs_cpeptide_group.age_years.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  lm(meta_1_perc ~ cpeptide_group + age_years,
     data=.)
summary(lm.meta_1_perc_vs_cpeptide_group.age_years.slow_rapid_only)
# plot(lm.meta_1_perc_vs_cpeptide_group.age_years.slow_rapid_only)

lm.meta_11_perc_vs_cpeptide_group.age_years.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  lm(meta_11_perc ~ cpeptide_group + age_years,
     data=.)
summary(lm.meta_11_perc_vs_cpeptide_group.age_years.slow_rapid_only)
# plot(lm.meta_11_perc_vs_cpeptide_group.age_years.slow_rapid_only)

lm.meta_12_perc_vs_cpeptide_group.age_years.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  lm(meta_12_perc ~ cpeptide_group + age_years,
     data=.)
summary(lm.meta_12_perc_vs_cpeptide_group.age_years.slow_rapid_only)
# plot(lm.meta_12_perc_vs_cpeptide_group.age_years.slow_rapid_only)


##### fit models of metacluster to C-peptide group, adjusting for log_age_years #####

## models with interaction
lm.meta_1_perc_vs_cpeptide_group.log_age_years.with_interaction.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  lm(meta_1_perc ~ cpeptide_group * log(age_years),
     data=.)
summary(lm.meta_1_perc_vs_cpeptide_group.log_age_years.with_interaction.slow_rapid_only)
# plot(lm.meta_1_perc_vs_cpeptide_group.log_age_years.with_interaction.slow_rapid_only)
# interaction term is not significant; proceed without interaction terms

lm.meta_11_perc_vs_cpeptide_group.log_age_years.with_interaction.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  lm(meta_11_perc ~ cpeptide_group * log(age_years),
     data=.)
summary(lm.meta_11_perc_vs_cpeptide_group.log_age_years.with_interaction.slow_rapid_only)
# plot(lm.meta_11_perc_vs_cpeptide_group.log_age_years.with_interaction.slow_rapid_only)
# interaction term is marginally significant; what to do?

lm.meta_12_perc_vs_cpeptide_group.log_age_years.with_interaction.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  lm(meta_12_perc ~ cpeptide_group * log(age_years),
     data=.)
summary(lm.meta_12_perc_vs_cpeptide_group.log_age_years.with_interaction.slow_rapid_only)
# plot(lm.meta_12_perc_vs_cpeptide_group.log_age_years.with_interaction.slow_rapid_only)
# interaction term is not significant; proceed without interaction terms


## models without interaction term
lm.meta_1_perc_vs_cpeptide_group.log_age_years.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  lm(meta_1_perc ~ cpeptide_group + log(age_years),
     data=.)
summary(lm.meta_1_perc_vs_cpeptide_group.log_age_years.slow_rapid_only)
# plot(lm.meta_1_perc_vs_cpeptide_group.log_age_years.slow_rapid_only)

lm.meta_11_perc_vs_cpeptide_group.log_age_years.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  lm(meta_11_perc ~ cpeptide_group + log(age_years),
     data=.)
summary(lm.meta_11_perc_vs_cpeptide_group.log_age_years.slow_rapid_only)
# plot(lm.meta_11_perc_vs_cpeptide_group.log_age_years.slow_rapid_only)

lm.meta_12_perc_vs_cpeptide_group.log_age_years.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  lm(meta_12_perc ~ cpeptide_group + log(age_years),
     data=.)
summary(lm.meta_12_perc_vs_cpeptide_group.log_age_years.slow_rapid_only)
# plot(lm.meta_12_perc_vs_cpeptide_group.log_age_years.slow_rapid_only)


##### make plots with common slope #####

data.tmp <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  cbind(
    meta_1_perc_predicted =
      predict(lm.meta_1_perc_vs_cpeptide_group.age_years.slow_rapid_only),
    meta_11_perc_predicted =
      predict(lm.meta_11_perc_vs_cpeptide_group.age_years.slow_rapid_only),
    meta_12_perc_predicted =
      predict(lm.meta_12_perc_vs_cpeptide_group.age_years.slow_rapid_only))

pdf(paste0(
  "plots/meta_1_percent_vs_age_at_draw.by_cpeptide_group.common_slope.", data_date, ".pdf"),
  w=9,h=6)
ggplot(
  data=data.tmp,
  mapping=aes(
    x=age_years, y=meta_1_perc,
    colour=cpeptide_group, group=cpeptide_group)) +
  geom_point(size=4, shape=16) +
  geom_smooth(
    method="lm",
    mapping=aes(y=meta_1_perc_predicted),
    size=3.0, se=FALSE, color="black", linetype="42") + # control dash length to match
  geom_smooth(
    method="lm",
    mapping=aes(y=meta_1_perc_predicted),
    size=2.0, se=FALSE, linetype="63") + # control dash length to match
  scale_color_manual(
    "C-peptide\ngroup",
    values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 1 %")
dev.off()

pdf(paste0(
  "plots/meta_11_percent_vs_age_at_draw.by_cpeptide_group.common_slope.", data_date, ".pdf"),
  w=9,h=6)
ggplot(
  data=data.tmp,
  mapping=aes(
    x=age_years, y=meta_11_perc,
    colour=cpeptide_group, group=cpeptide_group)) +
  geom_point(size=4, shape=16) +
  geom_smooth(
    method="lm",
    mapping=aes(y=meta_11_perc_predicted),
    size=3.0, se=FALSE, color="black", linetype="42") + # control dash length to match
  geom_smooth(
    method="lm",
    mapping=aes(y=meta_11_perc_predicted),
    size=2.0, se=FALSE, linetype="63") + # control dash length to match
  scale_color_manual(
    "C-peptide\ngroup",
    values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 11 %")
dev.off()

pdf(paste0(
  "plots/meta_12_percent_vs_age_at_draw.by_cpeptide_group.common_slope.", data_date, ".pdf"),
  w=9,h=6)
ggplot(
  data=data.tmp,
  mapping=aes(
    x=age_years, y=meta_12_perc,
    colour=cpeptide_group, group=cpeptide_group)) +
  geom_point(size=4, shape=16) +
  geom_smooth(
    method="lm",
    mapping=aes(y=meta_12_perc_predicted),
    size=3.0, se=FALSE, color="black", linetype="42") + # control dash length to match
  geom_smooth(
    method="lm",
    mapping=aes(y=meta_12_perc_predicted),
    size=2.0, se=FALSE, linetype="63") + # control dash length to match
  scale_color_manual(
    "C-peptide\ngroup",
    values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 12 %")
dev.off()


##### fit logistic regression model of metacluster 11 to C-peptide group #####

## binarize the metacluster 11 levels
data.CyTOF$meta_11_binary <-
  data.CyTOF$meta_11_perc > 0

## with an interaction term
logit_lm.meta_11_perc_vs_cpeptide_group.age_years.with_interaction.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  glm(meta_11_binary ~ cpeptide_group * age_years, family="binomial",
      data=.)
summary(logit_lm.meta_11_perc_vs_cpeptide_group.age_years.with_interaction.slow_rapid_only)

## plot the logistic regression (from https://stackoverflow.com/a/36943518)

# Predict probability
pred_data.with_interaction.tmp <-
  expand.grid(
    age_years=
      seq(min(data.CyTOF$age_years[data.CyTOF$cpeptide_group %in% c("slow", "rapid")]),
          max(data.CyTOF$age_years[data.CyTOF$cpeptide_group %in% c("slow", "rapid")]),
          length.out=1000),
    cpeptide_group=c("slow", "rapid"))
pred_data.with_interaction.tmp$prob <-
  predict(logit_lm.meta_11_perc_vs_cpeptide_group.age_years.with_interaction.slow_rapid_only,
          pred_data.with_interaction.tmp, type="response")

pdf(paste0(
  "plots/meta_11_binary_logistic_vs_age_at_draw.by_cpeptide_group.", data_date, ".pdf"),
  w=9,h=6)
ggplot(pred_data.with_interaction.tmp,
       mapping=aes(x=age_years, y=prob, color=cpeptide_group, group=cpeptide_group)) +
  geom_line(size=1.5) +
  geom_point(
    data.CyTOF %>%
      dplyr::filter(cpeptide_group %in% c("slow","rapid")),
    mapping=aes(y=as.numeric(meta_11_binary)),
    size=3) +
  scale_color_manual(
    "C-peptide\ngroup",
    values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 11 cells > 0")
dev.off()

## without an interaction term
logit_lm.meta_11_perc_vs_cpeptide_group.age_years.slow_rapid_only <-
  data.CyTOF %>%
  dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  glm(meta_11_binary ~ cpeptide_group + age_years, family="binomial",
      data=.)
summary(logit_lm.meta_11_perc_vs_cpeptide_group.age_years.slow_rapid_only)
# plot(logit_lm.meta_11_perc_vs_cpeptide_group.age_years.slow_rapid_only)


## plot the logistic regression (from https://stackoverflow.com/a/36943518)

# Predict probability
pred_data.no_interaction.tmp <-
  expand.grid(
    age_years=
      seq(min(data.CyTOF$age_years[data.CyTOF$cpeptide_group %in% c("slow", "rapid")]),
          max(data.CyTOF$age_years[data.CyTOF$cpeptide_group %in% c("slow", "rapid")]),
          length.out=1000),
    cpeptide_group=c("slow", "rapid"))
pred_data.no_interaction.tmp$prob <-
  predict(logit_lm.meta_11_perc_vs_cpeptide_group.age_years.slow_rapid_only,
          pred_data.no_interaction.tmp, type="response")

pdf(paste0(
  "plots/meta_11_binary_logistic_vs_age_at_draw.by_cpeptide_group.", data_date, ".common_age_relationship.pdf"),
  w=9,h=6)
ggplot(pred_data.no_interaction.tmp,
       mapping=aes(x=age_years, y=prob, color=cpeptide_group, group=cpeptide_group)) +
  geom_line(size=1.5) +
  geom_point(
    data.CyTOF %>%
      dplyr::filter(cpeptide_group %in% c("slow","rapid")),
    mapping=aes(y=as.numeric(meta_11_binary)),
    size=3) +
  scale_color_manual(
    "C-peptide\ngroup",
    values=c("rapid"="red", "slow"="blue")) +
  labs(x="Age at draw", y="Metacluster 11 cells > 0")
dev.off()


rm_tmp(ask=FALSE)


##### fit logistic regression model of metacluster 11 to C-peptide group, with other groups #####

## with an interaction term
logit_lm.meta_11_perc_vs_cpeptide_group.age_years.with_interaction.all_groups <-
  data.CyTOF %>%
  # dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  glm(meta_11_binary ~ cpeptide_group * age_years, family="binomial",
      data=.)
summary(logit_lm.meta_11_perc_vs_cpeptide_group.age_years.with_interaction.all_groups)

## plot the logistic regression (from https://stackoverflow.com/a/36943518)

# Predict probability
pred_data.with_interaction.tmp <-
  expand.grid(
    age_years=
      seq(min(data.CyTOF$age_years[!is.na(data.CyTOF$cpeptide_group)]),
          max(data.CyTOF$age_years[!is.na(data.CyTOF$cpeptide_group)]),
          length.out=1000),
    cpeptide_group=levels(data.CyTOF$cpeptide_group))
pred_data.with_interaction.tmp$prob <-
  predict(logit_lm.meta_11_perc_vs_cpeptide_group.age_years.with_interaction.all_groups,
          pred_data.with_interaction.tmp, type="response")

pdf(paste0(
  "plots/meta_11_binary_logistic_vs_age_at_draw.by_cpeptide_group.all_groups.", data_date, ".pdf"),
  w=9,h=6)
ggplot(pred_data.with_interaction.tmp,
       mapping=aes(x=age_years, y=prob, color=cpeptide_group, group=cpeptide_group)) +
  geom_line(size=1.5) +
  geom_point(
    data.CyTOF %>%
      dplyr::filter(cpeptide_group %in% c("slow","rapid", "slow_too_short", "rapid_too_long")),
    mapping=aes(y=as.numeric(meta_11_binary)),
    size=3) +
  scale_color_manual(
    "C-peptide\ngroup",
    values=
      c("rapid"="red", "slow"="blue",
        "slow_too_short"="deepskyblue", "rapid_too_long"="firebrick4")) +
  labs(x="Age at draw", y="Metacluster 11 cells > 0")
dev.off()

## without an interaction term
logit_lm.meta_11_perc_vs_cpeptide_group.age_years.all_groups <-
  data.CyTOF %>%
  # dplyr::filter(cpeptide_group %in% c("slow","rapid")) %>%
  glm(meta_11_binary ~ cpeptide_group + age_years, family="binomial",
      data=.)
summary(logit_lm.meta_11_perc_vs_cpeptide_group.age_years.all_groups)
# plot(logit_lm.meta_11_perc_vs_cpeptide_group.age_years.all_groups)


## plot the logistic regression (from https://stackoverflow.com/a/36943518)

# Predict probability
pred_data.no_interaction.tmp <-
  expand.grid(
    age_years=
      seq(min(data.CyTOF$age_years[!is.na(data.CyTOF$cpeptide_group)]),
          max(data.CyTOF$age_years[!is.na(data.CyTOF$cpeptide_group)]),
          length.out=1000),
    cpeptide_group=levels(data.CyTOF$cpeptide_group))
pred_data.no_interaction.tmp$prob <-
  predict(logit_lm.meta_11_perc_vs_cpeptide_group.age_years.all_groups,
          pred_data.no_interaction.tmp, type="response")

pdf(paste0(
  "plots/meta_11_binary_logistic_vs_age_at_draw.by_cpeptide_group.all_groups.", data_date, ".common_age_relationship.pdf"),
  w=9,h=6)
ggplot(pred_data.no_interaction.tmp,
       mapping=aes(x=age_years, y=prob, color=cpeptide_group, group=cpeptide_group)) +
  geom_line(size=1.5) +
  geom_point(
    data.CyTOF %>%
      dplyr::filter(cpeptide_group %in% c("slow","rapid", "slow_too_short", "rapid_too_long")),
    mapping=aes(y=as.numeric(meta_11_binary)),
    size=3) +
  scale_color_manual(
    "C-peptide\ngroup",
    values=
      c("rapid"="red", "slow"="blue",
        "slow_too_short"="deepskyblue", "rapid_too_long"="firebrick4")) +
  labs(x="Age at draw", y="Metacluster 11 cells > 0")
dev.off()


rm_tmp(ask=FALSE)
