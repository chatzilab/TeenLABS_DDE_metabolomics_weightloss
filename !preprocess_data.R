# Prepare data for analyses
## Subset TeenLABS data by AT POPs availability
chemical_names <- c("HCB", "24_DDE", "44_DDE", "PCB_118", "2_4_DDT", "PCB_153", "44_DDT", "TCDD", "PBDE_47", "PBDE_85")
chemical_variables <- grep(paste(tolower(chemical_names),collapse="|"), colnames(data_wide), value=TRUE)
chemical_variables <- grep("\\_0$", chemical_variables, value = T)

subject_w_POPs <- data_wide %>%
  select(key, any_of(chemical_variables)) %>%
  filter(rowSums(is.na(.)) < length(chemical_variables)) %>%
  select(key) %>% unlist

subject_w_plasma <- unique(data_plasma_metabo[[1]]$key)

subject_include <- intersect(subject_w_POPs, subject_w_plasma) 

data_wide <- data_wide %>%
  filter(key %in% subject_include)
data_long <- data_long %>%
  filter(key %in% subject_include)
table(data_long$visit_new, data_long$visit)
sapply(data_wide[,c("visit_year_0", "visit_year_6", "visit_year_12", "visit_year_36", "visit_year_60")], table)

## Recoding
# Create BMI change variable: bmi_c
data_wide <- data_wide %>% 
  mutate(
    bmi0_0 = 0,
    bmi0_6 = (bmi_6 - bmi_0),
    bmi6_12 =(bmi_12 - bmi_6),
    bmi12_36 =(bmi_36 - bmi_12),
    bmi36_60 =(bmi_60 - bmi_36))

temp1 <- data_wide %>% pivot_longer(
  cols = bmi0_0:bmi36_60,
  values_to = "bmi_c",
  names_to = "name"
) %>% 
  separate(name, c(NA, "visit_new"), 
           sep = "_",
           remove = TRUE) %>% 
  select(key, visit_new, bmi_c) 

data_l2 <- data_long%>%
  tidylog::left_join(temp1, by = c("key","visit_new")) %>%
  dplyr::select(key:agemos, bmi_c,everything())

data_long <- data_l2 %>% janitor::clean_names() %>% 
  arrange(key, visit_new) %>%
  mutate(visit = factor(visit_new, levels = c("0", "6", "12", "36", "60"), ordered = TRUE))

data_wide <- data_wide %>% 
  mutate(
    bmi0_6 = (bmi_6 - bmi_0)/bmi_0,
    bmi0_12 =(bmi_12 - bmi_0)/bmi_0,
    bmi0_36 =(bmi_36 - bmi_0)/bmi_0,
    bmi0_60 =(bmi_60 - bmi_0)/bmi_0)

temp1 <- data_wide %>% pivot_longer(
  cols = c("bmi0_0", "bmi0_6", "bmi0_12", "bmi0_36", "bmi0_60"),
  values_to = "bmi_percent_c",
  names_to = "name"
) %>% 
  separate(name, c(NA, "visit_new"), 
           sep = "_",
           remove = TRUE) %>% 
  select(key, visit_new, bmi_percent_c) 

data_l2 <- data_long%>%
  tidylog::left_join(temp1, by = c("key","visit_new")) %>%
  dplyr::select(key:agemos, bmi_percent_c,everything())

data_long <- data_l2 %>% janitor::clean_names() %>% 
  arrange(key, visit_new) %>%
  mutate(visit = factor(visit_new, levels = c("0", "6", "12", "36", "60"), ordered = TRUE))

# Re-order visit ###############################################################
levels(data_long$visit_new) <- as.numeric(c('0','0.5','1','3','5'))

# Create parents_income variable with "unknown" as level #######################
data_long <- data_long %>%
  group_by(key) %>%
  fill(parents_income)
data_long <- data_long %>% mutate(parents_income_new = ifelse(is.na(parents_income), "unknown", parents_income))

# IMPORTANT: Exclude participant with >300 for HOMA-IR, as advised by Todd Jenkins #####

data_long <- data_long %>% 
  replace_with_na_at(.vars = "homa_0",
                     condition = ~ (.x) >300)

data_wide <- data_wide %>% 
  replace_with_na_at(.vars = "homa_0",
                     condition = ~ (.x) >300)

# Impute missing lip POPs levels ###############################################
# lower than LOD: pbde47: 1; pcb118: 1; pcb153: 2; 
lower_lod = 0.1/sqrt(2)
data_long = data_long %>% mutate(across(c(hcb_adipose_0,x44_dde_adipose_0
                                          ,x44_ddt_adipose_0,pcb_118_adipose_0,pcb_153_adipose_0,
                                          pbde_47_adipose_0),
                                        function(x) ifelse(x == 0, lower_lod,x)))

data_wide = data_wide %>% mutate(across(c(hcb_adipose_0,x44_dde_adipose_0
                                          ,x44_ddt_adipose_0,pcb_118_adipose_0,pcb_153_adipose_0,
                                          pbde_47_adipose_0),
                                        function(x) ifelse(x == 0, lower_lod,x)))

# Repeat values of AT POPs for all visits ######################################

data_long <- data_long %>%
  group_by(key) %>%
  fill(hcb_adipose_0)

data_long <- data_long %>%
  group_by(key) %>%
  fill(x44_dde_adipose_0)

data_long <- data_long %>%
  group_by(key) %>%
  fill(pcb_118_adipose_0)

data_long <- data_long %>%
  group_by(key) %>%
  fill(pcb_153_adipose_0)

data_long <- data_long %>%
  group_by(key) %>%
  fill(pbde_47_adipose_0)

data_long <- data_long %>%
  group_by(key) %>%
  fill(x44_ddt_adipose_0)


# Rename datasets to restrict to those that have AT POPs #######################

data <- data_long

data_w <- data_wide

data <- data %>%
  mutate(bmi_baseline = case_when(visit_new == 0 ~  bmi))

data <- data %>%
  group_by(key) %>%
  fill(bmi_baseline)

data <- data %>%
  mutate(age_baseline = case_when(visit_new == 0 ~  agemos))

data <- data %>%
  group_by(key) %>%
  fill(age_baseline)

data$age_baseline <- data$age_baseline/12

# Restrict to those we have AT POPs only #######################################
data <- data %>% drop_na(hcb_adipose_0)
data <- data %>% drop_na(key)

data_w <- data_w %>% drop_na(hcb_adipose_0)
data_w <- data_w %>% drop_na(key)

# Make sure variables are correctly classified correctly ######################

data$key <- as.factor(data$key)
data$sex <- as.factor(data$sex)
data$parents_income_new <- as.factor(data$parents_income_new)
data$race_binary <- as.factor(data$race_binary)
data$site <- as.factor(data$site)

data_w$homa_0 <- as.numeric(data_w$homa_0)
data_w$homa_6 <- as.numeric(data_w$homa_6)
data_w$homa_12 <- as.numeric(data_w$homa_12)
data_w$homa_36 <- as.numeric(data_w$homa_36)
data_w$homa_60 <- as.numeric(data_w$homa_60)

# Revalue visit variable #######################################################

data$visit <- as.numeric(revalue(as.character(data$visit_new), c("0" = "0", "6" = "0.5", "12" = "1", "36" = "3", "60" = "5")))
data$visit <- as.numeric(data$visit)

# Create total DDT
data$total_ddt_adipose_0 <- data$x44_dde_adipose_0 + data$x44_ddt_adipose_0

# Create groups for different levels of DDE exposure, defined by tertiles #######

## Find tertiles for  DDE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(x44_dde_tertile <- quantile(data_w$x44_dde_adipose_0, c(0:3/3)))
(x44_ddt_tertile <- quantile(data_w$x44_ddt_adipose_0, c(0:3/3)))
(total_ddt_tertile <- quantile(data %>% dplyr::filter(visit_new == 0) %>% ungroup %>% select(total_ddt_adipose_0) %>% unlist, c(0:3/3)))

## Categorise DDE/DDT variable, create new variable "dde cat" ~~~~~~~~~~~~~~~~~~~~~~~~~
data$dde_cat <- as.factor(ifelse(data$x44_dde_adipose_0 <= x44_dde_tertile[2], 'Low (<=14.0ng/g)',
                                 ifelse(data$x44_dde_adipose_0 > x44_dde_tertile[2] & data$x44_dde_adipose_0 <= x44_dde_tertile[3], 'Moderate (14.0-18.7ng/g)', 
                                        ifelse(data$x44_dde_adipose_0 > x44_dde_tertile[3], 'High (>=18.7ng/g)',0))))

data <- data %>% 
  mutate(dde_cat = fct_relevel(dde_cat, 
                               "Low (<=14.0ng/g)", "Moderate (14.0-18.7ng/g)", "High (>=18.7ng/g)"))
table(data$visit, data$dde_cat)

data$ddt_cat <- as.factor(ifelse(data$x44_ddt_adipose_0 <= x44_ddt_tertile[2], 'Low (<=0.6ng/g)',
                                 ifelse(data$x44_ddt_adipose_0 > x44_ddt_tertile[2] & data$x44_ddt_adipose_0 <= x44_ddt_tertile[3], 'Moderate (0.6-1.0ng/g)', 
                                        ifelse(data$x44_ddt_adipose_0 > x44_ddt_tertile[3], 'High (>=1.0ng/g)',0))))

data <- data %>% 
  mutate(ddt_cat = fct_relevel(ddt_cat, 
                               "Low (<=0.6ng/g)", "Moderate (0.6-1.0ng/g)", "High (>=1.0ng/g)"))
table(data$visit, data$ddt_cat)

data$total_ddt_cat <- as.factor(ifelse(data$total_ddt_adipose_0 <= total_ddt_tertile[2], 'Low (<=14.3ng/g)',
                                       ifelse(data$total_ddt_adipose_0 > total_ddt_tertile[2] & data$total_ddt_adipose_0 <= total_ddt_tertile[3], 'Moderate (14.3-19.4ng/g)', 
                                              ifelse(data$total_ddt_adipose_0 > total_ddt_tertile[3], 'High (>=19.4ng/g)',0))))

data <- data %>% 
  mutate(total_ddt_cat = fct_relevel(total_ddt_cat, 
                                     "Low (<=14.3ng/g)", "Moderate (14.3-19.4ng/g)", "High (>=19.4ng/g)"))
table(data$visit, data$total_ddt_cat)

## Relevel variables ###############################################################

data$parents_income_new <- relevel(data$parents_income_new, ref = "less than 25000")
data$dde_cat <- relevel(data$dde_cat, ref = "Low (<=14.0ng/g)")
data$ddt_cat <- relevel(data$ddt_cat, ref = "Low (<=0.6ng/g)")
data$site <- relevel(data$site, ref = "CIN")

data$bmi_baseline_cat<-ifelse(data$bmi_baseline>=50,"Above 50","Below 50")

# Label categories
data$race_binary <- factor(data$race_binary,
                           levels = c(0, 1),
                           labels = c("Others", "White or Caucasian"))

data$parents_income_new <- factor(data$parents_income_new,
                                  levels = c("less than 25000", "25000-74999", "75000 or more", "unknown"),
                                  labels = c("Less than $25000", "$25000 to $74000", "$75000 or more", "Unknown"))

data$site <- factor(data$site,
                    levels = c("CIN", "BCM"),
                    labels = c("A", "B"))

# Label variables
label(data$race_binary) <- "Race"
label(data$sex) <-"Sex"
label(data$age_baseline) <- "Age at baseline (in years)"
label(data$parents_income_new) <- "Parents income category"
label(data$site) <- "Site"
label(data$dde_cat) <- "DDE level category"
label(data$ddt_cat) <- "DDT level category"

# Create new variable for DDE: high vs low. Low = 1st and 2nd tertiles, High = 3rd tertile
data$dde_hl <- as.factor(ifelse(data$x44_dde_adipose_0 <= x44_dde_tertile[3], 'Low',
                                ifelse(data$x44_dde_adipose_0 > x44_dde_tertile[3], 'High',0)))
table(data$dde_hl)
data$dde_hl <- relevel(data$dde_hl, ref = "Low")

# Create a numeric dde_cat
data$dde_cat_num <- as.numeric(data$dde_cat)

# Plasma Metabolome
## Save feature mz and retention time separately
data_plasma_metabo_mzrt <- lapply(data_plasma_metabo, function(x) x[,3:dim(x)[2]])
data_plasma_metabo_mzrt <- lapply(data_plasma_metabo_mzrt, function(x) as.data.frame(colnames(x)))
for(i in 1:length(data_plasma_metabo_mzrt))
{
  data_plasma_metabo_mzrt[[i]]$feat_label <- paste0(names(data_plasma_metabo_mzrt)[i], "_feat_", 1:nrow(data_plasma_metabo_mzrt[[i]]))
  data_plasma_metabo_mzrt[[i]]$mz <- as.numeric(sub('_.*', '', data_plasma_metabo_mzrt[[i]]$`colnames(x)`))
  data_plasma_metabo_mzrt[[i]]$rt <- as.numeric(sub('.*_', '', data_plasma_metabo_mzrt[[i]]$`colnames(x)`))
}

## Create visit_new consistent with the one in survey dataset
for(i in 1:length(data_plasma_metabo))
{
  data_plasma_metabo[[i]]$visit_new <- as.integer(ifelse(data_plasma_metabo[[i]]$visit == 1, 0, data_plasma_metabo[[i]]$visit)) 
  data_plasma_metabo[[i]] <- data_plasma_metabo[[i]] %>%
    relocate(visit_new, .after=key)
}

## Remove the visit record at 24m
for(i in 1:length(data_plasma_metabo))
{
  data_plasma_metabo[[i]] <- data_plasma_metabo[[i]][data_plasma_metabo[[i]]$visit_new != 24, ]
}

## Convert the visit index to match the cohort data
for(i in 1:length(data_plasma_metabo))
{
  data_plasma_metabo[[i]] <- data_plasma_metabo[[i]] %>%
    arrange(key, visit_new) %>%
    mutate(visit = factor(visit_new, levels = c("0", "6", "12", "36"), ordered = TRUE))
  
  levels(data_plasma_metabo[[i]]$visit_new) <- as.numeric(c('0','0.5','1','3'))
  
  data_plasma_metabo[[i]]$visit <- as.numeric(revalue(as.character(data_plasma_metabo[[i]]$visit_new), c("0" = "0", "6" = "0.5", "12" = "1", "36" = "3")))
  data_plasma_metabo[[i]]$visit <- as.numeric(data_plasma_metabo[[i]]$visit)
}

## Remove the subjects without vAT POPs data
for(i in 1:length(data_plasma_metabo))
{
  data_plasma_metabo[[i]] <- data_plasma_metabo[[i]] %>%
    filter(key %in% subject_include) 
}

## Match the key and visit between cohort data and plasma data
data <- data %>%
  filter(visit != 5) %>%
  arrange(key, visit)

for(i in 1:length(data_plasma_metabo))
{
  data_plasma_metabo[[i]] <- merge(data[,c("key", "visit", "visit_new")], data_plasma_metabo[[i]],
                                   by = c("key", "visit", "visit_new"), all.x = T) 
} #automatically generate NAs for plasma datasets where there was no record for that visit

## Generate an intensity only dataset
data_abundance <- lapply(data_plasma_metabo, as.data.frame)
for(i in 1:length(data_plasma_metabo))
{
  rownames(data_abundance[[i]]) <- paste(data_plasma_metabo[[i]]$key, "_", data_plasma_metabo[[i]]$visit_new)
  data_abundance[[i]] <- data_abundance[[i]][,-c(1:3)]
}

## Generate plasma L1 data
feat_annot_l_plasma_endo_l1 <- modify(
  feat_annot_l,
  ~ .x %>% 
    filter(Standard_grp == "Endogenous") %>% 
    select(`colnames(x)`,Metabolite_Name)
)
feat_annot_l_plasma_endo_l1 <- modify2(
  feat_annot_l_plasma_endo_l1,
  names(feat_annot_l_plasma_endo_l1),
  ~ cbind(.x, mode = .y)
) %>% bind_rows() %>% as_tibble() %>% 
  group_by(Metabolite_Name) %>% 
  slice(-1) %>% ungroup
data_plasma_metabo_endo_l1 <- modify2(
  data_plasma_metabo,
  names(data_plasma_metabo),
  ~ .x %>% 
    select(key, visit, visit_new, all_of(feat_annot_l_plasma_endo_l1$`colnames(x)`[which(feat_annot_l_plasma_endo_l1$mode == .y)]))
)
data_plasma_metabo_endo_l1 <- Reduce(
  function(x,y) left_join(x, y, by = c("key", "visit", "visit_new")),
  data_plasma_metabo_endo_l1)
