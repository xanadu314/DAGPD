setwd("~/DAGPD/AMG/")
AMG_metadata <- read.table("./AMG_metadata.txt",header = T)
result <- AMG_metadata %>%
  group_by(scaffold, ko) %>%
  summarise(count = n(), .groups = "drop")

phage_host <- read.table("./phage_host.txt",header = T)

test <- merge(result,phage_host,by = "scaffold")
test$host <- gsub("\\.2$", "", test$host)
phage_lenth <- read.table("./phage_lenth.txt",header = T)
phage_AMG_count_length <- merge(test,phage_lenth,by="scaffold")
phage_AMG_count_length$phage_AMH_den <- phage_AMG_count_length$count/phage_AMG_count_length$phage_length*1000000

total_host_AMG <- read.table("./total_host_AMG.90_50.txt",header = F)
colnames(total_host_AMG) <- c("qseqid", "AMG","qlen", "slen", "length", "evalue", "pident", "qcovs")

total_host_AMG$host <- gsub("GEM@","",total_host_AMG$qseqid,)
total_host_AMG$host <- gsub("@.*","",total_host_AMG$host,)
total_host_AMG$host <- gsub("MAGs_co_837_.*", "MAGs_co_837",total_host_AMG$host)
total_host_AMG$host <- gsub("MAGs_co_97_.*", "MAGs_co_97",total_host_AMG$host)

total_host_AMG_ko <- merge(total_host_AMG,AMG_metadata,by="AMG")

total_host_AMG_ko_uni <- total_host_AMG_ko %>%
  distinct(qseqid, .keep_all = TRUE)
host_ko <- total_host_AMG_ko_uni[,c(9,11)]

host_ko_count <- host_ko %>%
  group_by(host, ko) %>%
  summarise(count = n(), .groups = "drop")

host_ko_count$host <- gsub("\\.1$", "", host_ko_count$host)
host_ko_count$host <- gsub("\\.2$", "", host_ko_count$host)

host_length <- read.table("./host_length.txt",header = T)
host_ko_count_length <- merge(host_ko_count,host_length,by="host")
host_ko_count_length$host_AMG_den <- host_ko_count_length$count/host_ko_count_length$length*1000000

host_ko_count_length$host_ko <- paste(host_ko_count_length$host,host_ko_count_length$ko,sep = "@")
host_ko_den <- host_ko_count_length[,c(6,5)]

phage_AMG_count_length$host_ko <- paste(phage_AMG_count_length$host,phage_AMG_count_length$ko,sep = "@" )

phage_ko_den <- phage_AMG_count_length[,c(1,7,6)]
phage_host_ko_den <- left_join(phage_ko_den,host_ko_den,by="host_ko")

phage_host_ko_den$enrich <- phage_host_ko_den$phage_AMH_den/phage_host_ko_den$host_AMG_den

phage_host_ko_den_uni <- phage_host_ko_den %>% distinct()


host_ko <- host_ko_count_length[,c(1,2)]
phage_ko_host <- phage_AMG_count_length[,c(1,2,4)]

host_ko_wide <- host_ko %>%
  pivot_wider(names_from = ko, values_from = ko, values_fn = list(ko = ~1), values_fill = list(ko = 0))
colnames(host_ko_wide)=c("host",paste0("host_",colnames(host_ko_wide)[-1],""))
phage_ko_host_wide <- phage_ko_host %>% 
  pivot_wider(names_from = ko, values_from = ko, values_fn = list(ko = ~1), values_fill = list(ko = 0))

phage_source <- read.table("./phage_source.txt",header = T)
phage_source_ko_wide <- merge(phage_source,phage_ko_host_wide,by="scaffold")
colnames(phage_source_ko_wide) <- c("scaffold","source","host",paste0("AMG_",colnames(phage_source_ko_wide)[c(-1,-2,-3)]))

test2 <- merge(phage_source_ko_wide,host_ko_wide,by="host")

find_mismatched_amgs <- function(data, amg_prefix = "AMG_", host_prefix = "host_") {
  amg_columns <- grep(paste0("^", amg_prefix), names(data), value = TRUE)
  
  host_columns <- grep(paste0("^", host_prefix), names(data), value = TRUE)
  
  mismatched_amgs <- amg_columns[!sapply(amg_columns, function(amg) {
    host_col <- gsub(amg_prefix, host_prefix, amg)
    host_col %in% host_columns
  })]
  
  return(mismatched_amgs)
}

find_mismatched_amgs(test2)

host_AMG_df <- data.frame(matrix(0, nrow = dim(test2)[1], ncol = length(gsub("AMG_","host_",find_mismatched_amgs(test2)))))
colnames(host_AMG_df) <- gsub("AMG_","host_",find_mismatched_amgs(test2))

AMG_phage_host_df <- cbind(test2,host_AMG_df)
  
  
library(tidyverse)
library(broom)

batch_fisher_test <- function(data, amg_prefix = "AMG_", host_prefix = "host_") {
  
  amg_names <- grep(paste0("^", amg_prefix), names(data), value = TRUE)
  
  results <- data.frame(
    AMG = character(),
    p_value = numeric(),
    odds_ratio = numeric(),
    conf_int_low = numeric(),
    conf_int_high = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (amg in amg_names) {
    host_col <- gsub(amg_prefix, host_prefix, amg)
    
    if (!host_col %in% colnames(data)) {
      warning(paste("Missing host column for", amg))
      next
    }
    
    xtab <- table(
      AMG = data[[amg]],
      Host = data[[host_col]]
    )
    
      warning(paste("Insufficient data for", amg))
      next
    }
    
    tryCatch({
      test <- fisher.test(xtab, conf.int = TRUE)
      
      results <- rbind(results, data.frame(
        AMG = amg,
        p_value = test$p.value,
        odds_ratio = test$estimate,
        conf_int_low = test$conf.int[1],
        conf_int_high = test$conf.int[2]
      ))
    }, error = function(e) {
      warning(paste("Error in", amg, ":", e$message))
    })
  }
  
  results$adj_p <- p.adjust(results$p_value, method = "BH")
  
  return(results)
}


results_AMG_host <- batch_fisher_test(AMG_phage_host_df)

results_AMG_host_90 <- batch_fisher_test(AMG_phage_host_df_90)
write.csv(results_AMG_host_90,"results_AMG_host_90_related_bac_host.csv")

results_AMG_host_95 <- batch_fisher_test(AMG_phage_host_df_95)
write.csv(results_AMG_host_95,"results_AMG_host_95_related_bac_host.csv")
batch_fisher_test2 <- function(data, source_col = "source", amg_prefix = "AMG_") {
  library(tidyverse)
  
  amg_columns <- grep(paste0("^", amg_prefix), names(data), value = TRUE)
  
  results <- tibble(
    AMG = character(),
    p_value = numeric(),
    odds_ratio = numeric(),
    conf_int_low = numeric(),
    conf_int_high = numeric(),
    method = character(),
    stringsAsFactors = FALSE
  )
  
  for (amg in amg_columns) {
    xtab <- table(data[[amg]], data[[source_col]])
    
    if (any(dim(xtab) < 2)) {
      warning(paste("Skipping", amg, ": insufficient dimensions (", paste(dim(xtab), collapse = "x"), ")"))
      next
    }
    
    tryCatch({
      test <- fisher.test(xtab, simulate.p.value = TRUE, B = 1e5, conf.int = TRUE)
      
      or <- ifelse(is.null(test$estimate), NA, test$estimate)
      ci <- ifelse(is.null(test$conf.int), c(NA, NA), test$conf.int)
      
      results <- results %>% add_row(
        AMG = amg,
        p_value = test$p.value,
        odds_ratio = or,
        conf_int_low = ci[1],
        conf_int_high = ci[2],
        method = "Fisher's Exact Test"
      )
    }, error = function(e) {
      warning(paste("Error in", amg, ":", e$message))
    })
  }
  
  results$adj_p <- p.adjust(as.numeric(results$p_value), method = "BH")
  
  results <- results[, c(
    
    "AMG",
    
    "p_value",
    
    "odds_ratio",
    
    "conf_int_low",
    
    "conf_int_high",
    
    "adj_p",
    
    "method"
    
  )]
  
  return(results)
}

results2_host_homolog <- batch_fisher_test2(AMG_phage_host_df)
results2_host_homolog_90 <- batch_fisher_test2(AMG_phage_host_df_90)
write.csv(results2_host_homolog_90,"results2_host_homolog_90_relate_animal_host.csv")

results2_host_homolog_95 <- batch_fisher_test2(AMG_phage_host_df_95)
write.csv(results2_host_homolog_95,"results2_host_homolog_95_relate_animal_host.csv")

library(broom)
set.seed(123)
test_data <- data.frame(
  AMG_1 = sample(c(0,1), 100, replace=TRUE),
  Host_Homolog_AMG_1 = sample(c(0,1), 100, replace=TRUE),
  Source = factor(rep(c("Chicken", "Pig", "Ruminant"), length.out = 100))
)

colnames(test2)[1] <- "Host"
batch_glm_analysis <- function(data, 
                               amg_prefix = "AMG_", 
                               host_prefix = "host_",
                               source_var = "source") {
  
  library(tidyverse)
  library(broom)
  
  amg_columns <- grep(paste0("^", amg_prefix), names(data), value = TRUE)
  
  results <- tibble(
    AMG = character(),
    Term = character(),
    Estimate = numeric(),
    Std_Error = numeric(),
    z_value = numeric(),
    p_value = numeric(),
    Odds_Ratio = numeric(),
    CI_low = numeric(),
    CI_high = numeric()
  )
  
  for (amg in amg_columns) {
    host_col <- gsub(amg_prefix, host_prefix, amg)
    
    if (!host_col %in% colnames(data)) {
      warning(paste("Missing host column for", amg))
      next
    }
    
    formula <- as.formula(paste(amg, "~", host_col, "+", source_var))
    
    tryCatch({
      model <- glm(formula, 
                   family = binomial(link = "logit"), 
                   data = data)
      
      model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
        rename(
          estimate = estimate,
          std.error = std.error,
          statistic = statistic
        ) %>%
        select(term, estimate, std.error, statistic, p_value, conf.low, conf.high)
      
      model_summary <- model_summary %>%
        filter(term != "(Intercept)") %>%
        mutate(
          AMG = amg,
          Odds_Ratio = estimate,
          CI_low = conf.low,
          CI_high = conf.high
        ) %>%
        select(AMG, term, estimate, std.error, statistic, p_value, Odds_Ratio, CI_low, CI_high)
      
      results <- bind_rows(results, model_summary)
    }, error = function(e) {
      warning(paste("Error in", amg, ":", e$message))
    })
  }
  
  results <- results %>%
    group_by(term) %>%
    mutate(adj_p = p.adjust(p_value, method = "BH")) %>%
    ungroup()
  
  return(results)
}

glm_results_chicken <- batch_glm_analysis(AMG_phage_host_df)

test3 <- AMG_phage_host_df
test3$source=as.factor(gsub("pig","a_pig",AMG_phage_host_df$source))

glm_results_pig <- batch_glm_analysis(test3)

test4 <- AMG_phage_host_df
test4$source=as.factor(gsub("ruminant","a_ruminant",test2$source))
glm_results_ruminant <- batch_glm_analysis(test4)

chick_test <- AMG_phage_host_df
chick_test$source <- as.factor(gsub("pig","a_other",AMG_phage_host_df$source))
chick_test$source <- as.factor(gsub("ruminant","a_other",chick_test$source))
glm_results_chick <- batch_glm_analysis(chick_test)

pig_test <- AMG_phage_host_df
pig_test$source <- as.factor(gsub("chicken","a_other",AMG_phage_host_df$source))
pig_test$source <- as.factor(gsub("ruminant","a_other",pig_test$source))
glm_results_pig <- batch_glm_analysis(pig_test)

ruminant_test <- AMG_phage_host_df
ruminant_test$source <- as.factor(gsub("chicken","a_other",AMG_phage_host_df$source))
ruminant_test$source <- as.factor(gsub("pig","a_other",ruminant_test$source))
glm_results_ruminant <- batch_glm_analysis(ruminant_test)


batch_glm_analysis_inter <- function(data, 
                               amg_prefix = "AMG_", 
                               host_prefix = "host_",
                               source_var = "source") {
  
  library(tidyverse)
  library(broom)
  
  amg_columns <- grep(paste0("^", amg_prefix), names(data), value = TRUE)
  
  results <- tibble(
    AMG = character(),
    Term = character(),
    Estimate = numeric(),
    Std_Error = numeric(),
    z_value = numeric(),
    p_value = numeric(),
    Odds_Ratio = numeric(),
    CI_low = numeric(),
    CI_high = numeric()
  )
  
  for (amg in amg_columns) {
    host_col <- gsub(amg_prefix, host_prefix, amg)
    
    if (!host_col %in% colnames(data)) {
      warning(paste("Missing host column for", amg))
      next
    }
    
    formula <- as.formula(paste(amg, "~", host_col, "*", source_var ))
    
    tryCatch({
      model <- glm(formula, 
                   family = binomial(link = "logit"), 
                   data = data)
      
      model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
        rename(
          estimate = estimate,
          std.error = std.error,
          statistic = statistic
        ) %>%
        select(term, estimate, std.error, statistic, p_value, conf.low, conf.high)
      
      model_summary <- model_summary %>%
        filter(term != "(Intercept)") %>%
        mutate(
          AMG = amg,
          Odds_Ratio = estimate,
          CI_low = conf.low,
          CI_high = conf.high
        ) %>%
        select(AMG, term, estimate, std.error, statistic, p_value, Odds_Ratio, CI_low, CI_high)
      
      results <- bind_rows(results, model_summary)
    }, error = function(e) {
      warning(paste("Error in", amg, ":", e$message))
    })
  }
  
  results <- results %>%
    group_by(term) %>%
    mutate(adj_p = p.adjust(p_value, method = "BH")) %>%
    ungroup()
  
  return(results)
}

glm_results_inter_chick <- batch_glm_analysis_inter(chick_test)

glm_results_inter_pig <- batch_glm_analysis_inter(pig_test)

glm_results_inter_ruminant <- batch_glm_analysis_inter(ruminant_test)


library(lme4)
library(broom.mixed)

batch_glmm_analysis <- function(data, 
                                amg_prefix = "AMG_", 
                                host_prefix = "host_",
                                random_effect = "(1 | Host)") {
  
  amg_columns <- grep(paste0("^", amg_prefix), names(data), value = TRUE)
  
  results <- tibble(
    AMG = character(),
    Term = character(),
    Estimate = numeric(),
    Std_Error = numeric(),
    z_value = numeric(),
    p_value = numeric(),
    Odds_Ratio = numeric(),
    CI_low = numeric(),
    CI_high = numeric()
  )
  
  for (amg in amg_columns) {
    host_col <- gsub(amg_prefix, host_prefix, amg)
    
    if (!host_col %in% colnames(data)) {
      warning(paste("Missing host column for", amg))
      next
    }
    
    formula <- reformulate(
      termlabels = c(host_col, fixed_effects, random_effect),
      response = amg
    )
    
    tryCatch({
      model <- glmer(
        formula,
        family = binomial(link = "logit"),
        data = data,
      )
      
      model_summary <- tidy(model, effects = "fixed", conf.int = TRUE) %>%
        mutate(
          AMG = amg,
          Odds_Ratio = exp(estimate),
          CI_low = exp(conf.low),
          CI_high = exp(conf.high)
        ) %>%
        select(AMG, term, estimate, std.error, statistic, p.value, Odds_Ratio, CI_low, CI_high)
      
      results <- bind_rows(results, model_summary)
    }, error = function(e) {
      warning(paste("Error in", amg, ":", e$message))
    })
  }
  
   
   return(results)
}



amg_data <- AMG_phage_host_df %>% select(contains(head(amg_name,3)))
amg_data <- cbind(AMG_phage_host_df[,c(1,2,3)],amg_data)

colnames(amg_data)[1] <- "Host"
glmm_results <- batch_glmm_analysis(amg_data)


