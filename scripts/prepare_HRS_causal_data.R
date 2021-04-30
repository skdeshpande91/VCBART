##################
# Causal VCBART
##################
library(tidyverse)
load("data/raw_joined_data.RData")
############################

# treatments: physical inactivity & depression
# depression is operationalized as "a score of 4 or more on the modified 8-item CES-D scale administed by the HRS"
# see the report: https://hrs.isr.umich.edu/sites/default/files/biblio/dr-005.pdf (Section VI.B)

data_tbl <-
  raw_data %>%
  filter(inw5 == 1 & inw6 & cogfunction1998 == 1 & cogfunction2000 == 1 & r6agem_b >= 720 & r6agem_b <= 1020) %>%
  rename(Y = r6cogtot,
         AGE = r6agem_b,
         SMOKE = r5smokev,
         cSEP = cses_index,
         EDUC = raedyrs) %>%
  mutate(PHYS_INACT = 1 - r5vigact,
         DEPRESS = ifelse(r5cesd >= 4, 1, 0),
         CHLD_HLTH = ifelse(6-f992 > 0, 6-f992, NA), # reverse code childhood health
         GENDER = if_else(ragender == 1, 0, 1), # male = 0, female = 1
         RACE = case_when(raracem == 1 & rahispan == 0 ~ 0, # non-hispanic white,
                          raracem == 2 & rahispan == 0 ~ 1, # non-hispanic black,
                          raracem == 1 ~ 2, # hispanic
                          raracem == 3 & rahispan == 0 ~ 3), # non-hispanic other
         SOUTHERN = if_else(rabplace %in% c(5,6,7), 1, 0),
         FOREIGN = if_else(rabplace == 11, 1, 0)) %>%
  filter(!is.na(SMOKE) & !is.na(cSEP) & !is.na(EDUC) & !is.na(PHYS_INACT) & !is.na(SMOKE) &  
         !is.na(DEPRESS) & !is.na(CHLD_HLTH) & !is.na(GENDER) & !is.na(RACE) & !is.na(SOUTHERN) & !is.na(FOREIGN)) %>%
  select(hhidpn, Y, PHYS_INACT, DEPRESS, AGE, cSEP, EDUC, CHLD_HLTH, GENDER, RACE, SOUTHERN, FOREIGN, SMOKE)

data_df <- as.data.frame(data_tbl)
rownames(data_df) <- paste0("id.", data_df[,"hhidpn"])
data_df <- data_df[,colnames(data_df) != "hhidpn"]

data_df[,"CHLD_HLTH"] <- as.factor(data_df[,"CHLD_HLTH"])
data_df[,"GENDER"] <- as.factor(data_df[,"GENDER"])
data_df[,"RACE"] <- as.factor(data_df[,"RACE"])
data_df[,"SOUTHERN"] <- as.factor(data_df[,"SOUTHERN"])
data_df[,"FOREIGN"] <- as.factor(data_df[,"FOREIGN"])
data_df[,"SMOKE"] <- as.factor(data_df[,"SMOKE"])


attrition_df <- data_df
attrition_df[,"Y"] <- 1*is.na(data_df[,"Y"])

attrition_fit <- glm(Y ~ ., data = attrition_df, family = binomial)
attrition_coef <- summary(attrition_fit)[["coefficients"]][c("PHYS_INACT", "DEPRESS"),]
if(any(attrition_coef[,"Pr(>|z|)"] < 0.05)){
  print(round(attrition_coef, digits = 3))
  stop("At least one treatment indicator is significantly associated with missing outcomes")
}

data_df <- data_df[-which(is.na(data_df[,"Y"])),]
save(attrition_df, data_df, file = "~/Dropbox/vcbart_revision_sims/data/HRS_causal_df.RData")
