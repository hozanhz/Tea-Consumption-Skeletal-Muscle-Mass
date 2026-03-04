library(biomaRt)
library(UniProt.ws)
library(AnnotationDbi)
library(haven)
library(dplyr)
library(readxl)
library(stringr)
library(ggplot2)
library(tidyr)
library(writexl)
library(tidyverse)
library(ggplot2)
library(cluster)
library(factoextra)
library(randomForest)
library(glmnet)
library(pracma)
library(lmerTest)
library(lme4)
library(lcmm)
library(figpatch)
library(gridExtra)
library(gridGraphics)
library(grid)
library(png)
library(maaslin3)
library(pacman)
library(microeco)
library(ggpubr)
library(scales)
library(viridis)
library(patchwork)
library(cowplot)
library(ggsankey)
library(ggrepel)
library(lubridate)
pacman::p_load(tidyverse,microeco,magrittr)
library(randomForest)
library(forcats)
library(compareGroups)
library(ggalluvial)
library(Hmisc)
library(circlize)
library(colorspace)
library(purrr)
library(stringr)
library(tidytext)
library(ggsci)
library(openxlsx)
library(ggtext)

rm(list=ls())
####***共用的函数***####
####数据前处理####
#******标准化
scale_columns <- function(data, cols) {
  for (col in cols) {
    new_col <- paste0(col, "_z")
    data[[new_col]] <- as.numeric(scale(data[[col]], center = TRUE, scale = TRUE))
  }
  return(data)
}
#******分组标准化
scale_columns_group <- function(data, vars, group_var) {
  data %>%
    group_by(.data[[group_var]]) %>%
    mutate(
      across(
        all_of(vars),
        ~ as.numeric(scale(.x), center = TRUE, scale = TRUE),
        .names = "{.col}_z"
      )
    ) %>%
    ungroup()
}
#******代谢物数据转换
prepare_metabolite_data2 <- function(data1, Metabolites_filtered, exclude_col) {
  
  Met <-setdiff(colnames(Metabolites_filtered), exclude_col)
  
  ## log2 转换
  data1_log <- data1 %>%
    mutate(across(
      all_of(Met),
      ~ log2(. + 1e-6)
    ))
  
  ## Z 分数标准化
  data1_log_z <- scale_columns(
    data1_log,
    Met
  )
  
  return(data1_log_z)
}

#******菌群过滤
perform_micro_filter <- function(
    data_analyse1,              # 含菌群+metadata的数据框
    Bacteria_species_F2,        # 仅用于筛选判断的数据框
    exclude_cols = c("ID", "id", "followup"),
    metadata_cols = c("ID","followup"),
    min_abundance = 0.0001,     # 平均丰度阈值
    min_sample_detection = 0.1  # 最低非零比例
) { 
  
  # 1. 提取菌群列
  species_cols <- setdiff(colnames(Bacteria_species_F2), exclude_cols)
  
  # 转换为数值型
  data_analyse1[, species_cols] <- lapply(
    data_analyse1[, species_cols],
    as.numeric
  )
  
  Bacteria_species_F2[, species_cols] <- lapply(
    Bacteria_species_F2[, species_cols],
    as.numeric
  )
  
  # 2. 低丰度筛选（按平均丰度）
  mean_abundance_all <- colMeans(
    Bacteria_species_F2[, species_cols],
    na.rm = TRUE
  )
  
  valid_species_abundance <- species_cols[
    mean_abundance_all > min_abundance
  ]
  
  # 3. 存在率筛选（非0比例）
  detection_rate <- colMeans(
    Bacteria_species_F2[, valid_species_abundance] > 0,
    na.rm = TRUE
  )
  
  valid_species_detection <- valid_species_abundance[
    detection_rate >= min_sample_detection
  ]
  
  # 4. 保留列
  all_cols_to_keep <- union(metadata_cols, valid_species_detection)
  
  data_qc <- data_analyse1[, all_cols_to_keep, drop = FALSE]
  
  # 5. 统计每个菌的 0 值情况

  micro_data_final <- data_qc[, valid_species_detection, drop = FALSE]
  
  zero_summary <- data.frame(
    Species = valid_species_detection,
    N_total = nrow(micro_data_final),
    N_zero = colSums(micro_data_final == 0, na.rm = TRUE),
    Zero_rate = colMeans(micro_data_final == 0, na.rm = TRUE),
    Detection_rate = colMeans(micro_data_final > 0, na.rm = TRUE),
    Mean_abundance = colMeans(micro_data_final, na.rm = TRUE)
  )
  
  # 按 Zero_rate 排序（从高到低）
  zero_summary <- zero_summary[order(-zero_summary$Zero_rate), ]
  
  # 6. 打印筛选信息

  cat("原始菌种数:", length(species_cols), "\n")
  cat("平均丰度筛选后:", length(valid_species_abundance), "\n")
  cat("存在率筛选后:", length(valid_species_detection), "\n")
  
  # 7. 返回结果

  return(list(
    filtered_data = data_qc,
    zero_summary = zero_summary
  ))
}

#******菌群clr转换
perform_micro_clr <- function(
    species_data,
    exclude_cols = c("ID", "Tea_freq", "Times")  # Times 放进默认排除列
) {
  library(dplyr)
  library(microbiome)
  
  ##  判断是否有随访次数变量
  has_times <- "Times" %in% colnames(species_data)
  
  ##  提取需要保留的元信息列
  meta_cols <- intersect(c("ID", "Times"), colnames(species_data))
  
  ##  提取物种丰度矩阵并加 pseudo-count
  species_only <- species_data %>%
    dplyr::select(-all_of(exclude_cols)) %>%
    as.matrix()
  
  species_only <- species_only + 1e-6
  
  ##  CLR 转换
  species_clr <- microbiome::transform(species_only, "clr") %>%
    as.data.frame()
  
  ## 合并回 ID / Times
  species_clr_out <- bind_cols(
    species_data[, meta_cols, drop = FALSE],
    species_clr
  )
  
  return(species_clr_out)
}
####计算新变量####
#******计算AUC
calculate_auc_slope <- function(prefixes = c("F0", "F2", "F3", "F4"), 
                                data, 
                                base_names) {
  # 创建结果数据框
  Metabolites_target_AUC_slope <- data.frame(ID = data$ID)
  
  # 遍历每个代谢物
  for (metab in base_names) {
    # 提取该代谢物在不同时间点的变量名
    metab_vars <- paste0(prefixes, metab)
    
    # 确保变量名都存在
    if (all(metab_vars %in% names(data))) {
      auc_values <- numeric(nrow(data))
      slope_values <- numeric(nrow(data))
      
      # 对每个受试者计算 AUC 和 slope
      for (i in 1:nrow(data)) {
        # 根据 prefixes 动态生成对应的 age 列名
        age_vars <- paste0(prefixes, "age")  # 去掉 "_" 后再加 age
        ages <- unlist(data[i, age_vars])
        
        # 取代谢物值
        values <- unlist(data[i, metab_vars])
        
        # 检查是否有缺失值
        if (any(is.na(ages)) || any(is.na(values))) {
          auc_values[i] <- NA
          slope_values[i] <- NA
        } else {
          # AUC：用梯形法积分
          auc_values[i] <- trapz(ages, values)
          # slope：线性回归的斜率
          slope_values[i] <- coef(lm(values ~ ages))[2]
        }
      }
      
      # 将结果添加到输出数据框
      Metabolites_target_AUC_slope[[paste0(metab, "_AUC")]] <- auc_values
      Metabolites_target_AUC_slope[[paste0(metab, "_slope")]] <- slope_values
    }
  }
  
  return(Metabolites_target_AUC_slope)
}

#******计算菌群等指数
# 如果只有有害菌或有益菌也可以计算
calculate_GMSI <- function(data,
                           outcomes_pos,
                           outcomes_neg,
                           GMSI_name,
                           time_var = "Times") {
  
  taxa_all <- unique(c(outcomes_pos, outcomes_neg))
  
  ## 判断是否为长数据
  has_time <- time_var %in% colnames(data)
  
  keep_cols <- c("ID", taxa_all)
  if (has_time) keep_cols <- c("ID", time_var, taxa_all)
  
  data_GMSI <- data[, keep_cols, drop = FALSE]
  
  ## 构建 abundance matrix（行 = 样本，而不是 ID）
  abund_mat <- as.matrix(
    data_GMSI[, taxa_all, drop = FALSE]
  )
  
  ## 行名唯一化（防止 ID 重复）
  rownames(abund_mat) <- if (has_time) {
    paste(data_GMSI$ID, data_GMSI[[time_var]], sep = "_")
  } else {
    data_GMSI$ID
  }
  
  MG <- unique(outcomes_pos)
  MP <- unique(outcomes_neg)
  
  calc_weighted_abundance <- function(x, taxa_set) {
    taxa_set <- intersect(taxa_set, names(x))
    if (length(taxa_set) == 0) return(NA_real_)
    
    abund_vec <- x[taxa_set]
    mean_abund <- mean(abund_vec, na.rm = TRUE)
    presence_ratio <- mean(abund_vec > 0, na.rm = TRUE)
    
    mean_abund * presence_ratio
  }
  
  epsilon <- 1e-8
  
  GMSI <- apply(abund_mat, 1, function(x) {
    
    MG_w <- calc_weighted_abundance(x, MG)
    MP_w <- calc_weighted_abundance(x, MP)
    
    has_MG <- !is.na(MG_w)
    has_MP <- !is.na(MP_w)
    
    if (has_MG & has_MP) {
      log10((MG_w + epsilon) / (MP_w + epsilon))
    } else if (has_MG & !has_MP) {
      log10((MG_w + epsilon) / epsilon)
    } else if (!has_MG & has_MP) {
      -log10((MP_w + epsilon) / epsilon)
    } else {
      NA_real_
    }
  })
  
  ## ===== 构造返回数据框
  if (has_time) {
    GMSI_df <- data.frame(
      ID = data_GMSI$ID,
      Times = data_GMSI[[time_var]],
      stringsAsFactors = FALSE
    )
  } else {
    GMSI_df <- data.frame(
      ID = data_GMSI$ID,
      stringsAsFactors = FALSE
    )
  }
  
  GMSI_df[[GMSI_name]] <- as.numeric(GMSI)
  
  return(GMSI_df)
}
####数据分析####
process_lmer <- function(X, Y, data, covariates, gender_label) {
  
  library(lme4)
  library(lmerTest)
  library(dplyr)
  
  # 用反引号保护变量名
  safe_var <- function(v) paste0("`", v, "`")
  
  results_list <- list()
  counter <- 1
  
  for (x in X) {
    for (y in Y) {
      
      # ========= 构建公式
      x_safe  <- safe_var(x)
      y_safe  <- safe_var(y)
      cov_safe <- sapply(covariates, safe_var)
      
      formula_text <- paste(
        y_safe, "~",
        x_safe, "+",
        paste(cov_safe, collapse = " + "),
        "+ (1|ID)"
      )
      
      formula_full <- as.formula(formula_text)
      
      # ========= 运行模型 =========
      lme_model <- tryCatch(
        suppressMessages(
          suppressWarnings(
            lmer(
              formula_full,
              data = data,
              control = lmerControl(
                optimizer = "bobyqa",
                optCtrl = list(maxfun = 1e5)
              )
            )
          )
        ),
        error = function(e) NULL
      )
      
      if (is.null(lme_model)) next
      
      is_singular <- isSingular(lme_model)
      
      # ========= 提取系数
      sum_model <- summary(lme_model)
      coef_df <- as.data.frame(sum_model$coefficients)
      coef_df$Term <- rownames(coef_df)
      
      # 去掉反引号
      coef_df$Term_clean <- gsub("`", "", coef_df$Term)
      
      # 匹配该变量所有行（连续或分类）
      coef_x <- coef_df[startsWith(coef_df$Term_clean, x), ]
      if (nrow(coef_x) == 0) next
      
      # ========= 提取 Level（不用正则）
      coef_x$Level <- substring(
        coef_x$Term_clean,
        nchar(x) + 1
      )
      
      coef_x$Level <- ifelse(
        coef_x$Level == "",
        "Continuous",
        coef_x$Level
      )
      
      # ========= 判断变量类型 
      variable_type <- if (is.factor(data[[x]])) {
        "Categorical"
      } else {
        "Continuous"
      }
      
      # ========= 计算CI和显著性
      coef_x <- coef_x %>%
        mutate(
          Outcome = y,
          Predictor = x,
          Variable_Type = variable_type,
          
          CI_low  = Estimate - 1.96 * `Std. Error`,
          CI_high = Estimate + 1.96 * `Std. Error`,
          P_value = `Pr(>|t|)`,
          
          Significance = ifelse(
            P_value < 0.05,
            "Significant",
            "Not Significant"
          ),
          
          Singular_Fit = is_singular,
          Warning_Message = ifelse(
            is_singular,
            "boundary (singular) fit",
            NA_character_
          )
        )
      
      # ========= 方差信息
      vc_df <- as.data.frame(VarCorr(lme_model))
      
      coef_x$Random_Intercept_Var <-
        vc_df$vcov[vc_df$grp == "ID"]
      
      coef_x$Residual_Var <-
        vc_df$vcov[vc_df$grp == "Residual"]
      
      # ========= 性别标签
      coef_x$Gender <- switch(
        as.character(gender_label),
        "0" = "Female",
        "1" = "Male",
        "2" = "All",
        NA
      )
      
      # ========= 存储
      results_list[[counter]] <-
        coef_x[, c("Outcome", "Predictor", "Variable_Type",
                   "Level", "Estimate", "Std. Error",
                   "CI_low", "CI_high",
                   "P_value", "Significance",
                   "Gender", "Warning_Message",
                   "Singular_Fit",
                   "Random_Intercept_Var",
                   "Residual_Var")]
      
      counter <- counter + 1
    }
  }
  
  bind_rows(results_list)
}

process_glm1 <- function(X, Y, data, covariates, gender_label) {
  
  gender_forest_data <- data.frame()
  
  bt <- function(x) paste0("`", x, "`")
  
  regex_escape <- function(x) {
    gsub("([][{}()+*^$.|\\\\?\\-])", "\\\\\\1", x)
  }
  
  for (x in X) {
    for (y in Y) {
      
      formula <- as.formula(
        paste(
          bt(y), "~",
          bt(x), "+",
          paste(bt(covariates), collapse = "+")
        )
      )
      
      model <- glm(formula, data = data, family = gaussian())
      
      coef_df <- as.data.frame(summary(model)$coefficients)
      coef_df$Term <- rownames(coef_df)
      
      # ⭐ 关键：用“裸 x”，而不是 `x`
      x_esc <- regex_escape(x)
      
      coef_x <- coef_df[
        grepl(paste0("^", x_esc), coef_df$Term),
        ,
        drop = FALSE
      ]
      
      if (nrow(coef_x) == 0) next
      
      coef_x$Outcome   <- y
      coef_x$Predictor <- x
      
      coef_x$Level <- gsub(paste0("^", x_esc), "", coef_x$Term)
      coef_x$Level <- ifelse(coef_x$Level == "", "Continuous", coef_x$Level)
      
      coef_x$CI_low  <- coef_x$Estimate - 1.96 * coef_x$`Std. Error`
      coef_x$CI_high <- coef_x$Estimate + 1.96 * coef_x$`Std. Error`
      
      coef_x$Significance <- ifelse(
        coef_x$`Pr(>|t|)` < 0.05,
        "Significant",
        "Not Significant"
      )
      
      coef_x$Gender <- switch(
        as.character(gender_label),
        "0" = "Female",
        "1" = "Male",
        "2" = "All",
        NA
      )
      
      gender_forest_data <- rbind(
        gender_forest_data,
        coef_x[, c(
          "Outcome", "Predictor", "Level",
          "Estimate", "Std. Error",
          "CI_low", "CI_high",
          "Pr(>|t|)", "Significance", "Gender"
        )]
      )
    }
  }
  
  colnames(gender_forest_data)[
    colnames(gender_forest_data) == "Pr(>|t|)"
  ] <- "P_value"
  
  gender_forest_data
}

process_logistic <- function(X, Y, data, covariates, gender_label) {
  
  library(dplyr)
  library(nnet)
  
  results <- data.frame()
  
  ## ---------- 统一变量名（关键）
  name_map <- data.frame(
    raw   = colnames(data),
    clean = make.names(colnames(data), unique = TRUE),
    stringsAsFactors = FALSE
  )
  
  colnames(data) <- name_map$clean
  
  X_clean  <- name_map$clean[match(X, name_map$raw)]
  Y_clean  <- name_map$clean[match(Y, name_map$raw)]
  cov_clean <- name_map$clean[match(covariates, name_map$raw)]
  
  ## ---------- 工具函数
  bt <- function(x) paste0("`", x, "`")
  
  regex_escape <- function(x) {
    gsub("([][{}()+*^$.|\\\\?\\-])", "\\\\\\1", x)
  }
  
  for (i in seq_along(X_clean)) {
    x <- X_clean[i]
    x_raw <- X[i]
    
    for (j in seq_along(Y_clean)) {
      y <- Y_clean[j]
      y_raw <- Y[j]
      
      y_var <- data[[y]]
      y_nlevel <- nlevels(factor(y_var))
      
      ## ---------- 公式（现在已经是 syntactic name，其实不用反引号也行）
      fml <- as.formula(
        paste(y, "~", x, "+", paste(cov_clean, collapse = "+"))
      )
      
      x_esc <- regex_escape(x)
      
      ## ===== 1. 二分类
      if (y_nlevel == 2) {
        
        model <- glm(fml, data = data, family = binomial())
        coef_summary <- coef(summary(model))
        rn <- rownames(coef_summary)
        
        x_rows <- rn[grepl(paste0("^", x_esc), rn)]
        if (length(x_rows) == 0) next
        
        for (coef_name in x_rows) {
          
          beta <- coef_summary[coef_name, "Estimate"]
          se   <- coef_summary[coef_name, "Std. Error"]
          pval <- coef_summary[coef_name, "Pr(>|z|)"]
          
          results <- rbind(
            results,
            data.frame(
              Outcome = y_raw,
              Outcome_level = levels(factor(y_var))[2],
              Predictor = x_raw,
              Term = coef_name,
              Beta = beta,
              OR = exp(beta),
              OR_95CI_low = exp(beta - 1.96 * se),
              OR_95CI_high = exp(beta + 1.96 * se),
              P_value = pval,
              P_signif = ifelse(pval < 0.05, "*", ""),
              Gender = c("Female", "Male", "All")[gender_label + 1],
              stringsAsFactors = FALSE
            )
          )
        }
      }
      
      ## ===== 2. 多分类
      if (y_nlevel >= 3) {
        
        model <- nnet::multinom(fml, data = data, trace = FALSE)
        summ <- summary(model)
        
        coef_mat <- summ$coefficients
        se_mat   <- summ$standard.errors
        
        out_levels <- rownames(coef_mat)
        x_cols <- colnames(coef_mat)[grepl(paste0("^", x_esc), colnames(coef_mat))]
        if (length(x_cols) == 0) next
        
        for (ol in out_levels) {
          for (cn in x_cols) {
            
            beta <- coef_mat[ol, cn]
            se   <- se_mat[ol, cn]
            zval <- beta / se
            pval <- 2 * (1 - pnorm(abs(zval)))
            
            results <- rbind(
              results,
              data.frame(
                Outcome = y_raw,
                Outcome_level = paste(ol, "vs ref"),
                Predictor = x_raw,
                Term = cn,
                Beta = beta,
                OR = exp(beta),
                OR_95CI_low = exp(beta - 1.96 * se),
                OR_95CI_high = exp(beta + 1.96 * se),
                P_value = pval,
                P_signif = ifelse(pval < 0.05, "*", ""),
                Gender = c("Female", "Male", "All")[gender_label + 1],
                stringsAsFactors = FALSE
              )
            )
          }
        }
      }
    }
  }
  
  results
}

#******重复测量中介效应分析
run_mediation_glmer_batch <- function(
    data,
    id_var = "ID",
    exposures,      # list of exposures
    mediators,      # list of mediators
    outcomes,       # list of outcomes
    covariates,
    exposure_high = 3,
    exposure_low  = 0,
    sims = 1000,
    seed = 908
) {
  
  library(dplyr)
  library(lme4)
  library(mediation)
  
  # 给变量名加反引号
  safe_var <- function(v) paste0("`", v, "`")
  
  results <- list()
  
  for(x in exposures){
    for(m in mediators){
      for(y in outcomes){
        
        # 1. 构造 3 vs 0 的暴露
        df_med <- data %>%
          filter(.data[[x]] %in% c(exposure_low, exposure_high)) %>%
          mutate(
            Treat_bin = as.numeric(.data[[x]] == exposure_high)
          )
        
        # 跳过样本量太少的组合
        if(nrow(df_med) < 10 || length(unique(df_med$Treat_bin)) < 2) next
        
        # 2. 构建协变量字符串（加反引号）
        cov_string <- paste(safe_var(covariates), collapse = " + ")
        
        # 3. 构建中介模型公式
        formula_m <- as.formula(
          paste0(
            safe_var(m), " ~ Treat_bin",
            if(length(covariates) > 0) paste0(" + ", cov_string),
            " + (1 | ", safe_var(id_var), ")"
          )
        )
        
        # 4. 构建结局模型公式
        formula_y <- as.formula(
          paste0(
            safe_var(y), " ~ Treat_bin + ", safe_var(m),
            if(length(covariates) > 0) paste0(" + ", cov_string),
            " + (1 | ", safe_var(id_var), ")"
          )
        )
        
        # 5. 运行模型
        model.m <- tryCatch(
          suppressWarnings(
            glmer(formula_m, data = df_med, family = gaussian())
          ),
          error = function(e) NULL
        )
        if(is.null(model.m)) next
        
        model.y <- tryCatch(
          suppressWarnings(
            glmer(formula_y, data = df_med, family = gaussian())
          ),
          error = function(e) NULL
        )
        if(is.null(model.y)) next
        
        # 6. 中介分析
        set.seed(seed)
        med_res <- tryCatch(
          mediate(model.m, model.y, treat = "Treat_bin", mediator = m, sims = sims),
          error = function(e) NULL
        )
        if(is.null(med_res)) next
        
        # 7. 提取结果
        res <- data.frame(
          Exposure   = x,
          Mediator   = m,
          Outcome    = y,
          ACME       = med_res$d0,
          ACME_p     = med_res$d0.p,
          ADE        = med_res$z0,
          ADE_p      = med_res$z0.p,
          Total      = med_res$tau.coef,
          Total_p    = med_res$tau.p,
          PropMediated = med_res$n0,
          PropMediated_p = med_res$n0.p
        )
        
        results[[length(results)+1]] <- res
      }
    }
  }
  
  # 合并成总表
  final_res <- do.call(rbind, results)
  return(final_res)
}
#******FC的计算（基于均数）
run_FC_pairwise_analysis <- function(
    data,
    group_var = "Tea_freq",
    col,
    fc_cutoff = 0.6,
    p_cutoff = 0.05,
    fdr_cutoff = 0.25,
    pseudocount = 1e-8,
    use_fdr = TRUE,
    ref_level = 0
) {
  
  ## 转长
  species_long <- data %>%
    pivot_longer(
      cols = col,
      names_to = "Species",
      values_to = "Abundance"
    )
  
  ## 自动比较组
  group_levels <- sort(unique(data[[group_var]]))
  group_levels <- group_levels[!is.na(group_levels)]
  test_levels  <- setdiff(group_levels, ref_level)
  comparisons  <- lapply(test_levels, function(g) c(g, ref_level))
  
  ## 统计检验
  pairwise_test_auto <- function(df, g1, g2) {
    
    x1 <- df %>% filter(.data[[group_var]] == g1) %>% pull(Abundance)
    x2 <- df %>% filter(.data[[group_var]] == g2) %>% pull(Abundance)
    
    ## 去掉 NA
    x1 <- x1[!is.na(x1)]
    x2 <- x2[!is.na(x2)]
    
    ## 如果任一组样本量太小，直接 Wilcoxon
    if (length(x1) < 3 || length(x2) < 3) {
      test <- wilcox.test(x1, x2)
      return(tibble(P_value = test$p.value, Method = "wilcox"))
    }
    
    ## 如果没有变异（所有值一样），不能做 Shapiro
    var1 <- stats::var(x1)
    var2 <- stats::var(x2)
    
    if (is.na(var1) || is.na(var2) || var1 == 0 || var2 == 0) {
      test <- wilcox.test(x1, x2)
      return(tibble(P_value = test$p.value, Method = "wilcox"))
    }
    
    ## 正态性检验
    normal1 <- shapiro.test(x1)$p.value > 0.05
    normal2 <- shapiro.test(x2)$p.value > 0.05
    
    if (normal1 && normal2) {
      test <- t.test(x1, x2)
      method <- "t.test"
    } else {
      test <- wilcox.test(x1, x2)
      method <- "wilcox"
    }
    
    tibble(
      P_value = test$p.value,
      Method  = method
    )
  }
  
  pairwise_summary_stats <- function(df, g1, g2) {
    tibble(
      mean_g1   = mean(df$Abundance[df[[group_var]] == g1], na.rm = TRUE),
      mean_g2   = mean(df$Abundance[df[[group_var]] == g2], na.rm = TRUE),
      median_g1 = median(df$Abundance[df[[group_var]] == g1], na.rm = TRUE),
      median_g2 = median(df$Abundance[df[[group_var]] == g2], na.rm = TRUE)
    )
  }
  
  ## 统计部分
  stat_results <- species_long %>%
    group_by(Species) %>%
    group_modify(~ {
      map_dfr(comparisons, function(comp) {
        g1 <- comp[1]; g2 <- comp[2]
        bind_cols(
          pairwise_test_auto(.x, g1, g2),
          pairwise_summary_stats(.x, g1, g2)
        ) %>%
          mutate(Comparison = paste0(g1, "_vs_", g2))
      })
    }) %>%
    ungroup()
  
  ## FC
  fc_long <- species_long %>%
    group_by(Species, .data[[group_var]]) %>%
    summarise(
      mean_abundance = mean(Abundance, na.rm = TRUE) + pseudocount,
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = .data[[group_var]],
      values_from = mean_abundance,
      names_prefix = "Tea_"
    ) %>%
    {
      fc_df <- .
      map_dfr(test_levels, function(g) {
        tibble(
          Species = fc_df$Species,
          Comparison = paste0(g, "_vs_", ref_level),
          FC = fc_df[[paste0("Tea_", g)]] / fc_df[[paste0("Tea_", ref_level)]]
        )
      })
    } %>%
    mutate(log2FC = log2(FC))
  
  ## 合并 + FDR
  final_results <- stat_results %>%
    left_join(fc_long, by = c("Species", "Comparison")) %>%
    group_by(Comparison) %>%
    mutate(FDR = p.adjust(P_value, method = "BH")) %>%
    ungroup()
  
  final_sig <- final_results %>%
    filter(
      abs(log2FC) >= fc_cutoff,
      P_value < p_cutoff,
      if (use_fdr) FDR < fdr_cutoff else TRUE
    )
  
  list(
    all_results = final_results,
    sig_results = final_sig
  )
}
####可视化####
#******火山图
plot_volcano_species <- function(
    df,
    title = NULL,
    fc_cutoff = 0.6,
    p_cutoff = 0.02,
    fdr_cutoff = 0.25,
    label = TRUE,
    comparison_labels = c(
      "1_vs_0" = "Occasionally vs. Never",
      "2_vs_0" = "Regularly vs. Never",
      "2_vs_1" = "Regularly vs. Occasionally"
    )
) {
  
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(ggrepel)
  
  plot_df <- df %>%
    mutate(
      # sig = case_when(
      #   abs(log2FC) >= fc_cutoff & P_value < p_cutoff & FDR < fdr_cutoff & log2FC > 0 ~ "Up",
      #   abs(log2FC) >= fc_cutoff & P_value < p_cutoff & FDR < fdr_cutoff & log2FC < 0 ~ "Down",
      #   TRUE ~ "Not significant"
      # ),
      sig = case_when(
        abs(log2FC) >= fc_cutoff & P_value < p_cutoff & log2FC > 0 ~ "Up",
        abs(log2FC) >= fc_cutoff & P_value < p_cutoff & log2FC < 0 ~ "Down",
        TRUE ~ "Not significant"
      ),
      Species = Species %>%
        str_remove("^s__") %>%
        str_replace_all("_", " "),
      Comparison = factor(
        recode(Comparison, !!!comparison_labels),
        levels = unname(comparison_labels)
      )
    )
  
  p <- ggplot(plot_df, aes(x = log2FC, y = -log10(P_value))) +
    geom_point(aes(color = sig), size = 2, alpha = 0.8) +
    
    scale_color_manual(
      values = c(
        "Up" = "#D73027",
        "Down" = "#4575B4",
        "Not significant" = "grey70"
      )
    ) +
    
    geom_vline(
      xintercept = c(-fc_cutoff, fc_cutoff),
      linetype = "dashed",color = "grey80"
    ) +
    geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed",color = "grey80"
    ) +
    
    facet_wrap(~ Comparison, nrow = 1) +
    
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 20),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 18),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 20)
    ) +
    labs(
      title = title,
      x = "log2(Fold Change)",
      y = "-log10(P value)",
      color = NULL
    )
  
  if (label) {
    p <- p +
      geom_text_repel(
        data = plot_df %>%
          filter(
            P_value < p_cutoff,
            abs(log2FC) >= fc_cutoff#,
            #FDR < fdr_cutoff
          ),
        aes(label = Species),
        size = 3,
        force = 15,
        fontface = "italic",
        max.overlaps = 1000
      )
  }
  
  return(p)
}
#******组间比较
plot_wilcox_violin3 <- function(
    data,
    group_var,
    y_var,
    group_levels = NULL,
    group_labels = NULL,
    ylab = NULL,
    y_unit = NULL,
    xlab = NULL,
    title = NULL,
    base_size = 14,
    violin_fill = "#E6E6E6"
) {
  
  library(ggplot2)
  library(dplyr)
  library(rlang)
  
  group_sym <- ensym(group_var)
  group_chr <- as_string(group_sym)
  
  if (is.character(y_var)) {
    y_chr <- y_var
  } else {
    y_chr <- as_string(ensym(y_var))
  }
  
  ## 分组 factor
  if (!is.null(group_levels)) {
    data <- data %>%
      mutate(
        !!group_sym := factor(
          !!group_sym,
          levels = group_levels,
          labels = group_labels %||% group_levels
        )
      )
  } else {
    data <- data %>%
      mutate(!!group_sym := factor(!!group_sym))
  }
  
  ## 组数
  n_group <- nlevels(data[[group_chr]])
  
  ## 总体检验
  p_val <- if (n_group == 2) {
    wilcox.test(
      data[[y_chr]] ~ data[[group_chr]],
      exact = FALSE
    )$p.value
  } else {
    kruskal.test(
      data[[y_chr]] ~ data[[group_chr]]
    )$p.value
  }
  
  ## P 值文本（plotmath）
  p_text <- if (p_val < 0.001) {
    "italic(P)<0.001"
  } else {
    paste0("italic(P)==\"", sprintf("%.3f", p_val), "\"")
  }
  
  ## y 轴范围
  ymax   <- max(data[[y_chr]], na.rm = TRUE)
  yrange <- diff(range(data[[y_chr]], na.rm = TRUE))
  
  ## 标签位置
  x_p <- mean(seq_len(n_group))
  y_p <- ymax + 0.15 * yrange
  
  ## y 轴标题
  y_label_final <- if (!is.null(ylab)) {
    if (!is.null(y_unit)) {
      paste0(ylab, " (", y_unit, ")")
    } else {
      ylab
    }
  } else {
    y_chr
  }
  
  ggplot(data, aes(x = !!group_sym, y = .data[[y_chr]])) +
    
    ## violin
    geom_violin(
      fill = violin_fill,
      color = NA,
      alpha = 0.6,
      width = 0.85
    ) +
    
    ## box
    geom_boxplot(
      width = 0.28,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.6,
      median.linewidth = 0.9
    ) +
    
    ## median 点
    stat_summary(
      fun = median,
      geom = "point",
      size = 1.6,
      shape = 16,
      color = "black"
    ) +
    
    ## 灰色底框
    annotate(
      "rect",
      xmin = x_p - 0.6,
      xmax = x_p + 0.6,
      ymin = y_p - 0.05 * yrange,
      ymax = y_p + 0.05 * yrange,
      fill = "grey95",
      color = "black"
    ) +
    
    ## 斜体 P 值（关键）
    annotate(
      "text",
      x = x_p,
      y = y_p,
      label = p_text,
      parse = TRUE,
      size = 4.2,
      color = "black"
    ) +
    
    labs(
      title = title,
      x = xlab %||% group_chr,
      y = y_label_final
    ) +
    
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      axis.line  = element_line(linewidth = 0.6),
      plot.margin = ggplot2::margin(
        t = 10, r = 10, b = 10, l = 10, unit = "pt"
      )
    )
}
####***theme主题***####
theme_myself <- theme_classic(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #strip.background = element_rect(fill = "#333333", color = NA),
    strip.background = element_rect(fill = "grey85",color = NA),
    strip.text.x = element_text(face = "bold", size = 16),
    strip.text.y = element_text(face = "bold",size = 16,angle = 0),
    axis.text.y = element_markdown(face =  "italic",size = 12),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )
####*********混杂*********####
####只取基线的混杂####
#*****************************************4000之前基线混杂
Cov_F0 <- read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_基本信息_F01234_20230214_VF.xlsx",sheet="F0")
# 计算准确的F0年龄
F0_surveyTimes <- read_sav("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/5010调查时间及随访间隔（重新编制）__6-1-2017_VF.sav")
Cov_F0 <- merge(Cov_F0, F0_surveyTimes[,c("编号","调查时间_F0")], by.x = "ID", by.y = "编号")
Cov_F0$`调查时间_F0` <- as.Date(Cov_F0$`调查时间_F0`)
head(Cov_F0$`出生日期_F0`)
Cov_F0$`出生日期_F0` <- as.Date(
  as.numeric(Cov_F0$`出生日期_F0`),
  origin = "1899-12-30"
)
Cov_F0$Age_F0 <- time_length(
  difftime(Cov_F0$调查时间_F0, Cov_F0$出生日期_F0), 
  unit = "years"
)
# 只保留想要的变量
Cov_F0 <- Cov_F0[,c("ID","Age_F0","出生日期_F0","性别_F0","本人教育程度分类_3_F0","家庭人均收入_3_F0","是否吸烟_F0","是否饮酒_F0","是否喝茶_F0","过去一年是否钙片_F0","过去一年是否复合维生素_F0",
                    "脑卒中_F0","心脏病_F0","癫痫_F0","帕金森_F0","老年痴呆_F0","糖尿病_F0","癌症_F0","骨折_F0","雌激素_F0","绝经年限_F0"
)]
colnames(Cov_F0) <- dplyr::recode(colnames(Cov_F0),
                                  "出生日期_F0" = "Birthday_F0",
                                  "性别_F0" = "Sex_F0",
                                  "本人教育程度分类_3_F0" = "Education_F0",
                                  "家庭人均收入_3_F0" = "Income_F0",
                                  "是否吸烟_F0" = "Smoke_F0",
                                  "是否饮酒_F0" = "Alcohol_F0",
                                  "是否喝茶_F0" = "Tea_F0",
                                  "过去一年是否钙片_F0" = "Calcium_F0",
                                  "过去一年是否复合维生素_F0" = "Vitamin_F0",
                                  "骨折_F0" = "Fracture_F0")
#绝经方法：女(0)：取值为0或NA的话即还没绝经赋值为0，如果有年龄数字的话则赋值为1；男(0)：都赋值为0
Cov_F0$Menopause_F0 <- ifelse(Cov_F0$Sex_F0 == 1,0, 
                              ifelse(Cov_F0$`绝经年限_F0`=="NA" | Cov_F0$`绝经年限_F0` == 0, 0, 1))
#只保留4000以内的观测
Cov_F0 <- Cov_F0 %>%
  filter(!grepl("^NL4", ID))
# a <- Cov_F0
# a$Age_F0 <- round(a$Age_F0,0)
# writexl::write_xlsx(a[,c("ID","Age_F0","Sex_F0","Education_F0")], "Document.xlsx")
#*****************************************4000之后基线混杂
Cov_F1 <- read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_基本信息_F01234_20230214_VF.xlsx",sheet="F1")
Cov_F1$ID <- sub("^F1", "", Cov_F1$CODE_F1)
Cov_F1 <- Cov_F1 %>%
  filter(grepl("^NL4", ID))
Cov_F1 <- Cov_F1[,c("ID","性别_F1","家庭人均收入_3_F1","是否吸烟_F1","是否饮酒_F1","是否喝茶_F1","过去一年是否钙片_F1","过去一年是否复合维生素_F1",
                    "脑卒中_F1","心脏病_F1","癫痫_F1","帕金森病_F1","老年痴呆_F1","糖尿病_F1","癌症_F1","骨折_F1","绝经年限_F1", "雌激素_F1"
)]
Birthday_NL4 <- read_sav("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/其它/F1NL教育20210806_VF.sav") 
#openxlsx::write.xlsx(Birthday_NL4, "F1NL教育20210806_VF.xlsx")
Cov_F1 <- merge(Cov_F1, Birthday_NL4[,c("CODE", "出生日期","本人教育程度_3")], by.x = "ID", by.y = "CODE")
Cov_F1$`出生日期` <- as.Date(Cov_F1$`出生日期`)
#绝经方法：女(0)：取值为0或NA或-1的话即还没绝经赋值为0，如果有年龄数字的话则赋值为1；男(0)：都赋值为0
Cov_F1$Menopause_F0 <- ifelse(Cov_F1$`性别_F1` == 1,0, 
                              ifelse(Cov_F1$`绝经年限_F1`=="NA" | Cov_F1$`绝经年限_F1` == 0 | Cov_F1$`绝经年限_F1` == -1, 0, 1))
colnames(Cov_F1) <- dplyr::recode(colnames(Cov_F1),
                                  "出生日期" = "Birthday_F0",
                                  "性别_F1" = "Sex_F0",
                                  "本人教育程度_3" = "Education_F0",
                                  "家庭人均收入_3_F1" = "Income_F0",
                                  "是否吸烟_F1" = "Smoke_F0",
                                  "是否饮酒_F1" = "Alcohol_F0",
                                  "是否喝茶_F1" = "Tea_F0",
                                  "过去一年是否钙片_F1" = "Calcium_F0",
                                  "过去一年是否复合维生素_F1" = "Vitamin_F0",
                                  "骨折_F1" = "Fracture_F0",
                                  "脑卒中_F1" = "脑卒中_F0",
                                  "心脏病_F1" = "心脏病_F0",
                                  "癫痫_F1" = "癫痫_F0",
                                  "帕金森病_F1" = "帕金森_F0",
                                  "老年痴呆_F1" = "老年痴呆_F0",
                                  "糖尿病_F1" = "糖尿病_F0",
                                  "癌症_F1" = "癌症_F0",
                                  "绝经年限_F1" = "绝经年限_F0",
                                  "雌激素_F1" = "雌激素_F0")
# 计算F1的年龄
Cov_F1$CODE <- paste0("F1",Cov_F1$ID)
Follow_time <- read.csv("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/FollowData_VF.csv")
Follow_time$ID <- substr(Follow_time$CODE, 3, nchar(Follow_time$CODE))
Follow_time_BD  <- merge(Follow_time, Cov_F1[,c("CODE","Birthday_F0")], by = "CODE")
Follow_time_BD <- Follow_time_BD[complete.cases(Follow_time_BD),]

Follow_time_BD$FollowDate <- as.Date(Follow_time_BD$FollowDate)
Follow_time_BD$Birthday_F0 <- as.Date(Follow_time_BD$Birthday_F0)
Follow_time_BD$Age_F0 <- time_length(
  difftime(Follow_time_BD$FollowDate, Follow_time_BD$Birthday_F0), 
  unit = "years"
)
# 合并年龄
Cov_F1 <- merge(Cov_F1, Follow_time_BD[,c("ID","Age_F0")], by = "ID")
Cov_F1 <- dplyr::select(Cov_F1,-c("CODE"))
#*****************************************合并基线4000之前之后的基线混杂
Cov_F0_all <- rbind(Cov_F0, Cov_F1)
#雌激素为NA的是男性，所以除了女性有1外，其它直接赋值为0即可
Cov_F0_all$Estrogen_F0 <- ifelse(Cov_F0_all$`雌激素_F0` == 1,1,0)
table(Cov_F0_all$`雌激素_F0`)
table(Cov_F0_all$Estrogen)
#疾病只要之一为1，则新变量Disease就为1
Cov_F0_all$Disease_F0 <- ifelse(Cov_F0_all$`脑卒中_F0` == 1 | Cov_F0_all$`心脏病_F0` == 1 | Cov_F0_all$`癫痫_F0` == 1 | Cov_F0_all$`帕金森_F0` == 1 | Cov_F0_all$`老年痴呆_F0` == 1 | Cov_F0_all$`糖尿病_F0` == 1 | Cov_F0_all$`癌症_F0` == 1, 1, 0)
colnames(Cov_F0_all)
Cov_F0_all <- Cov_F0_all[,c("ID","Age_F0","Birthday_F0","Sex_F0","Education_F0","Smoke_F0", "Alcohol_F0","Tea_F0","Calcium_F0",
                            "Vitamin_F0","Fracture_F0","Menopause_F0","Estrogen_F0","Disease_F0")]
#取值为9也定义为无骨折
Cov_F0_all$Fracture_F0[is.na(Cov_F0_all$Fracture_F0) | Cov_F0_all$Fracture_F0 == "NA" | Cov_F0_all$Fracture_F0 == 9] <- 0
#剔除变量值取值缺失或者为NA的
Cov_F0_all <- Cov_F0_all[!apply(Cov_F0_all, 1, function(row) any(is.na(row) | row == "NA")), ]
####查看每个变量的取值情况####
Cov_F0_all_subset <- Cov_F0_all[, !colnames(Cov_F0_all) %in% c("ID", "Birthday_F0")]
# 创建一个空的数据框，用于存储所有变量的频率表
table_df <- data.frame(Variable = character(), Level = character(), Count = numeric(), stringsAsFactors = FALSE)
# 循环处理每一列，生成频率表
for (col in colnames(Cov_F0_all_subset)) {
  freq_table <- as.data.frame(table(Cov_F0_all_subset[[col]]))
  colnames(freq_table) <- c("Level", "Count")
  freq_table$Variable <- col
  table_df <- rbind(table_df, freq_table)
}
####多次随访的混杂####
#*************************************Met-宽数据
Cov2 <- read_excel('/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_身高体重MET_VF.xlsx', 
                   sheet=1)
Cov2_wide <- Cov2 %>%
  mutate(across(-ID, as.numeric))
# 4000之前的取F0、F1、F3的均值；4000之后的取F1、F2、F3的均值;允许有缺失值。
Cov2_wide$Met_mean <- ifelse(
  grepl("^NL4", Cov2_wide$ID), 
  rowMeans(Cov2_wide[, c("一天总MET_F1", "一天总MET_F2", "一天总MET_F3")], na.rm = TRUE),  
  rowMeans(Cov2_wide[, c("一天总MET_F0", "一天总MET_F1", "一天总MET_F3")], na.rm = TRUE) 
)
#*************************************营养素-宽数据
Nutrition <- read_excel('/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_F0123膳食营养素数据_20230206_VF.xlsx',
                        sheet="GNHS_膳食数据_20230206")
Nutrition$ID <- substr(Nutrition$`编号`, 3, nchar(Nutrition$`编号`))
colnames(Nutrition)
Nutrition <- Nutrition[,c("ID","Followup","能量摄入")]
####宽数据
Nutrition_wide <- Nutrition %>%
  # 使用 pivot_wider 转换为宽格式
  pivot_wider(
    names_from = Followup, 
    values_from = c(`能量摄入`),
    names_glue = "F{Followup}{.value}"
  ) %>%
  # 保留 ID 列
  dplyr::select(ID, everything()) %>%
  as.data.frame()
####计算均值
colnames(Nutrition_wide)
Nutrition_wide <- Nutrition_wide %>%
  mutate(across(where(is.list), ~ as.character(.))) %>%  # 先转换 list 为字符
  mutate(across(-ID, as.numeric))  # 再转换数值

Nutrition_wide$Energy_mean <- ifelse(
  grepl("^NL4", Nutrition_wide$ID), 
  rowMeans(Nutrition_wide[, c("F1能量摄入","F2能量摄入","F3能量摄入")], na.rm = TRUE),  
  rowMeans(Nutrition_wide[, c("F0能量摄入", "F1能量摄入", "F3能量摄入")], na.rm = TRUE) 
)
####鱼油补充剂####
Oil1 <- read_excel('/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_基本信息_F01234_20230214_VF.xlsx',
                   sheet=1)
Oil1 <- Oil1[,c("ID","是否使用深海鱼油胶囊_F2","过去一年是否服用深海鱼油_F3","否服用过深海鱼油胶囊_F1")]

Oil2 <- read_excel('/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_基本信息_F01234_20230214_VF.xlsx',
                   sheet="F4")
Oil2 <- Oil2[,c("CODE1","深海鱼油")]
colnames(Oil2) <- recode(colnames(Oil2),
                         "深海鱼油" = "深海鱼油_F4")

Oil2$CODE1 <- substr(Oil2$CODE1, 3, nchar(Oil2$CODE1))

Oil <- dplyr::left_join(Oil1, Oil2,  by = c("ID" = "CODE1"))
Oil_complete <- Oil %>%
  filter(!if_all(
    c("是否使用深海鱼油胶囊_F2", 
      "过去一年是否服用深海鱼油_F3", 
      "否服用过深海鱼油胶囊_F1", 
      "深海鱼油_F4"),
    ~ is.na(.) | . == "NA"
  ))
colnames(Oil_complete)
# 宽数据
fish_cols <- grep("鱼油", colnames(Oil), value = TRUE)
Oil_complete_wide <- Oil_complete
Oil_complete_wide$Oil_sup_follow <- ifelse(
  rowSums(Oil_complete_wide[, fish_cols] == 1, na.rm = TRUE) > 0, 
  1, 0
)

# 长数据
factor_name<-c("否服用过深海鱼油胶囊_F1", 
               "是否使用深海鱼油胶囊_F2", 
               "过去一年是否服用深海鱼油_F3", 
               "深海鱼油_F4")
idx <- which(names(Oil_complete)%in% factor_name)
for(i in idx ){
  Oil_complete[[i]] <-  as.factor(Oil_complete[[i]])
}

Oil_long <- Oil_complete %>%
  pivot_longer(
    cols = c("否服用过深海鱼油胶囊_F1", 
             "是否使用深海鱼油胶囊_F2", 
             "过去一年是否服用深海鱼油_F3", 
             "深海鱼油_F4"),
    names_to = "Times",
    values_to = "Oil_sup_follow_lmer"
  ) %>%
  mutate(
    Times = case_when(
      grepl("_F1$", Times) ~ "F1",
      grepl("_F2$", Times) ~ "F2",
      grepl("_F3$", Times) ~ "F3",
      grepl("_F4$", Times) ~ "F4",
      TRUE ~ NA_character_
    )
  )

Oil_long <- Oil_long[complete.cases(Oil_long) & Oil_long$Oil_sup_follow_lmer != "NA",]
####合并所有混杂，并定义定性定量变量####
data_list <- list(Cov_F0_all,Cov2_wide[,c("ID", "Met_mean")],Nutrition_wide[,c("ID","Energy_mean")], Oil_complete_wide[,c("ID","Oil_sup_follow")])
Cov_F0 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID")
Cov_F0 <- Cov_F0[complete.cases(Cov_F0),]
#*************************定性定量定义好
numeric_name <- c("Met_mean","Energy_mean","Age_F0")
for(n in numeric_name){
  Cov_F0[[n]] <- as.numeric(Cov_F0[[n]])
}

factor_name<-c("Sex_F0","Smoke_F0","Alcohol_F0","Tea_F0","Calcium_F0","Vitamin_F0","Disease_F0",
               "Fracture_follow","Fracture_Medicine","Estrogen_F0","Menopause_F0","Fracture_F0","Oil_sup_follow")
idx <- which(names(Cov_F0)%in% factor_name)
for(i in idx ){
  Cov_F0[[i]] <-  as.factor(Cov_F0[[i]])
}
####身高数据####
Height_F123 <- Cov2_wide[, c("ID", "身高_F1","身高_F2","身高_F3")]
Height_F4 <- read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_基本信息_F01234_20230214_VF.xlsx",sheet="F4")
colnames(Height_F4) <- recode(colnames(Height_F4),
                              "身高" = "身高_F4")
Height_F4$ID <- substr(Height_F4$CODE1, 3, nchar(Height_F4$CODE1))
Height <- left_join(Height_F123, Height_F4[,c("ID","身高_F4")])
Height_long <- Height %>%
  pivot_longer(cols = starts_with("身高_F"),
               names_to = "Followup",
               values_to = "Height") %>%
  mutate(Followup = sub("身高_", "", Followup))  %>%
  drop_na()
Height_long <- Height_long[Height_long$Height != 0,]
####年龄数据####
#****************长数据
Age_long <- left_join(Follow_time, Cov_F0_all[,c("ID","Birthday_F0")], by = "ID")
Age_long$FollowDate <- as.Date(Age_long$FollowDate)
Age_long$Birthday_F0 <- as.Date(Age_long$Birthday_F0)
Age_long$Age <- time_length(
  difftime(Age_long$FollowDate, Age_long$Birthday_F0), 
  unit = "years" 
)

Age_long$Times <- substr(Age_long$CODE, 1, 2)

Age_long1 <- Age_long[,c("Age","Times","ID")]
Age_long2 <- Cov_F0_all[,c("ID","Age_F0")]
colnames(Age_long2) <- recode(colnames(Age_long2),
                              "Age_F0" = "Age")
Age_long2$Times <- "F0"
Age_long_all <- rbind(Age_long1, Age_long2)
#****************转换为宽数据
# Age_wide <- Age_long %>%
#   select(ID, Times, Age) %>%   # 只保留需要的列
#   pivot_wider(
#     names_from  = Times,
#     values_from = Age,
#     names_prefix = "Age_"
#   )
Age_wide <- Age_long %>%
  dplyr::select(ID, Times, Age) %>%
  pivot_wider(
    names_from  = Times,
    values_from = Age,
    names_glue  = "{Times}age"
  )
####*********暴露*********####
####饮茶频率-综合指标####
# 因为变量"每周冲茶次数_F2"前面很多缺失值，所以会被识别为分类变量，所以要加上col_types = c(每周冲茶次数_F2 = "numeric")
Tea <- as.data.frame(read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_饮茶相关变量索取_更新_VF.xlsx",sheet=1, guess_max = 10000))
# F2的两个变量都不为缺失值的：没有
Tea_both <- Tea[!is.na(Tea$每周冲茶次数_F2) & !is.na(Tea$平均每周泡茶次数_F2), ]
# 取两个变量全集
Tea$每周泡茶次数_F2 <- dplyr::coalesce(
  Tea$每周冲茶次数_F2,
  Tea$平均每周泡茶次数_F2
)


Tea_Follow <- Tea[,c("ID","是否喝茶_F0","每周冲茶次数_F0","是否喝茶_F1","每周冲茶次数_F1","过去一年是否经常喝茶_F2","每周冲茶次数_F2","平均每周泡茶次数_F2","每周泡茶次数_F2","是否喝茶_F3","每周冲茶次数_F3")]
# Tea_Follow <- Tea_Follow[!is.na(Tea_Follow$是否喝茶_F0),]
Tea_Follow <- Tea[,c("ID","是否喝茶_F0","每周冲茶次数_F0","是否喝茶_F1","每周冲茶次数_F1","过去一年是否经常喝茶_F2","每周泡茶次数_F2","是否喝茶_F3","每周冲茶次数_F3")]
Tea_Follow <- Tea_Follow %>%
  mutate(
    每周冲茶次数_F0 = if_else(是否喝茶_F0 == 0, 0, 每周冲茶次数_F0),
    每周冲茶次数_F1 = if_else(是否喝茶_F1 == 0, 0, 每周冲茶次数_F1),
    每周泡茶次数_F2 = if_else(过去一年是否经常喝茶_F2 == 0, 0, 每周泡茶次数_F2),
    每周冲茶次数_F3 = if_else(是否喝茶_F3 == 0, 0, 每周冲茶次数_F3)
  )

a <- Tea_Follow %>%
  summarise(
    across(
      c(每周冲茶次数_F0,
        每周冲茶次数_F1,
        每周泡茶次数_F2,
        每周冲茶次数_F3),
      ~ sum(!is.na(.x))
    )
  )
# tea_cut <- function(x) {
#   dplyr::case_when(
#     x >= 1 & x <= 3 ~ 1,
#     x >= 4 & x <= 7 ~ 2,
#     x > 7 ~ 3,
#     TRUE            ~ 0
#   )
# }

# tea_cut <- function(x) {
#   dplyr::case_when(
#     x >= 1 & x <= 3 ~ 1,
#     x >= 4 & x < 7 ~ 2,
#     x == 7 ~ 3,
#     x > 7 ~ 4,
#     TRUE ~ 0
#   )
# }

tea_cut <- function(x) {
  dplyr::case_when(
    x >= 1 & x < 7 ~ 1,
    x == 7 ~ 2,
    x > 7 ~ 3,
    TRUE ~ 0
  )
}

calc_tea_freq <- function(df, vars) {
  df %>%
    filter(complete.cases(across(all_of(vars)))) %>%
    mutate(across(all_of(vars), tea_cut)) %>%
    rowwise() %>%
    mutate(
      Tea_freq = case_when(
        #max(c_across(all_of(vars))) == 4 ~ 4,
        max(c_across(all_of(vars))) == 3 ~ 3,
        max(c_across(all_of(vars))) == 2 ~ 2,
        max(c_across(all_of(vars))) == 1 ~ 1,
        min(c_across(all_of(vars))) == 0 ~ 0,
        TRUE ~ NA_real_
      )
    ) %>%
    ungroup() #%>%
    #dplyr::select(ID, Tea_freq)
}

Tea_Follow_all <- bind_rows(
  calc_tea_freq(
    filter(Tea_Follow, !grepl("^NL4", ID)),
    c("每周冲茶次数_F0","每周冲茶次数_F1","每周泡茶次数_F2","每周冲茶次数_F3")
  ),
  calc_tea_freq(
    filter(Tea_Follow, grepl("^NL4", ID)),
    c("每周冲茶次数_F1","每周泡茶次数_F2","每周冲茶次数_F3")
  )
)


Tea_Follow <- Tea_Follow_all[,c("ID","Tea_freq")]
Tea_Follow <- Tea_Follow[complete.cases(Tea_Follow),]
table(Tea_Follow$Tea_freq)

# 转换为分类变量
factor_name<-c("Tea_freq")
Tea_Follow <- Tea_Follow %>%
  mutate(across(all_of(factor_name), as.factor))
####验证茶暴露赋值####
# Tea_Follow_noNL4 <- Tea_Follow[!grepl("^NL4", Tea_Follow$ID), ]
# Tea_Follow_noNL4 <- Tea_Follow_noNL4[complete.cases(Tea_Follow_noNL4),]
# 
# Tea_Follow_NL4 <- Tea_Follow[grepl("^NL4", Tea_Follow$ID), c("ID","是否喝茶_F1","每周冲茶次数_F1","过去一年是否经常喝茶_F2","每周泡茶次数_F2","是否喝茶_F3","每周冲茶次数_F3")]
# Tea_Follow_NL4 <- Tea_Follow_NL4[complete.cases(Tea_Follow_NL4),]
# 
# 
# tea_cut <- function(x) {
#   dplyr::case_when(
#     x >= 1 & x <= 3 ~ 1,
#     x >= 4 & x <= 7 ~ 2,
#     x > 7 ~ 3,
#     TRUE            ~ 0
#   )
# }
# 
# # tea_cut <- function(x) {
# #   dplyr::case_when(
# #     x >= 1 & x < 7 ~ 1,
# #     x >= 7           ~ 2,
# #     TRUE            ~ 0
# #   )
# # }
# 
# Tea_Follow_noNL4 <- Tea_Follow_noNL4 %>%
#   mutate(
#     Tea_freq_F0 = tea_cut(每周冲茶次数_F0),
#     Tea_freq_F1 = tea_cut(每周冲茶次数_F1),
#     Tea_freq_F2 = tea_cut(每周泡茶次数_F2),
#     Tea_freq_F3 = tea_cut(每周冲茶次数_F3)
#   )
# 
# # 根据随访情况赋值
# Tea_Follow_noNL4 <- Tea_Follow_noNL4 %>%
#   mutate(Tea_freq = case_when(
#     pmax(Tea_freq_F0, Tea_freq_F1, Tea_freq_F2, Tea_freq_F3, na.rm = TRUE) == 3 ~ 3,
#     pmax(Tea_freq_F0, Tea_freq_F1, Tea_freq_F2, Tea_freq_F3, na.rm = TRUE) == 2 ~ 2,
#     pmax(Tea_freq_F0, Tea_freq_F1, Tea_freq_F2, Tea_freq_F3, na.rm = TRUE) == 1 ~ 1,
#     pmin(Tea_freq_F0, Tea_freq_F1, Tea_freq_F2, Tea_freq_F3, na.rm = TRUE) == 0 ~ 0,
#     TRUE ~ NA_real_
#   ))
# 
# 
# Tea_Follow_NL4 <- Tea_Follow_NL4 %>%
#   mutate(
#     Tea_freq_F1 = tea_cut(每周冲茶次数_F1),
#     Tea_freq_F2 = tea_cut(每周泡茶次数_F2),
#     Tea_freq_F3 = tea_cut(每周冲茶次数_F3)
#   )
# 
# # 根据随访情况赋值
# Tea_Follow_NL4 <- Tea_Follow_NL4 %>%
#   mutate(Tea_freq = case_when(
#     pmax(Tea_freq_F1, Tea_freq_F2, Tea_freq_F3, na.rm = TRUE) == 3 ~ 3,
#     pmax(Tea_freq_F1, Tea_freq_F2, Tea_freq_F3, na.rm = TRUE) == 2 ~ 2,
#     pmax(Tea_freq_F1, Tea_freq_F2, Tea_freq_F3, na.rm = TRUE) == 1 ~ 1,
#     pmin(Tea_freq_F1, Tea_freq_F2, Tea_freq_F3, na.rm = TRUE) == 0 ~ 0,
#     TRUE ~ NA_real_
#   ))
# 
# 
# Tea_Follow <- rbind(Tea_Follow_noNL4[,c("ID","Tea_freq")], Tea_Follow_NL4[,c("ID","Tea_freq")])
####F0血清儿茶素####
Catechin <- read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/血清氧化应激炎症因子血黄酮尿黄酮/20250527_GNHS_F0-3_VF.xlsx",sheet=1)
Catechin <- Catechin[,c("ID","serum_I_catechin_F0", "serum_I_epicatechin_F0", "serum_I_EGC_F0", "serum_I_EGCG_F0", "serum_I_ECG_F0")]
Catechin <- Catechin[complete.cases(Catechin),]
Catechin_sex <- merge(Catechin, Cov_F0[,c("ID","Sex_F0")], by = "ID")

zero_tertile4 <- function(data, vars, suffix = "_T") {
  
  data %>%
    dplyr::mutate(
      dplyr::across(
        all_of(vars),
        ~ {
          x <- .
          out <- rep(NA_integer_, length(x))
          
          ## 0 单独一类
          out[x == 0] <- 0
          
          ## 非 0 且非缺失（逻辑索引）
          nz <- x != 0 & !is.na(x)
          
          if (sum(nz) > 0) {
            q <- stats::quantile(x[nz], probs = c(1/3, 2/3), na.rm = TRUE)
            
            out[nz & x <= q[1]] <- 1
            out[nz & x >  q[1] & x <= q[2]] <- 2
            out[nz & x >  q[2]] <- 3
          }
          
          out
        },
        .names = paste0("{.col}", suffix)
      )
    )
}

quartile_group <- function(x) {
  
  out <- rep(NA_integer_, length(x))
  
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  out[x <= q[1]]                 <- 0  # Q1
  out[x > q[1] & x <= q[2]]      <- 1  # Q2
  out[x > q[2] & x <= q[3]]      <- 2  # Q3
  out[x > q[3]]                  <- 3  # Q4
  
  return(out)
}
####*********DXA-骨骼肌*********####
DXA_all <- read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_DXA_follow_10807条记录_20241014_2024.10.18recheck_VF.xlsx")
DXA_all <- DXA_all[,-1]

DXA_all$ASM <- rowSums(DXA_all[, c("LARM_LEAN", "RARM_LEAN", "L_LEG_LEAN", "R_LEG_LEAN")] - 
                         DXA_all[, c("LARM_BMC", "RARM_BMC", "LLEG_BMC", "RLEG_BMC")])
####DXA长数据####
ASM <- DXA_all[,c("CODE","ASM","WBTOT_FAT")]
ASM <- ASM %>%
  mutate(Followup = substr(CODE, 1, 2),
         CODE = substr(CODE, 3, nchar(CODE))) %>%
  filter(!Followup %in% c("F5")) %>%
  filter(grepl("^NL", CODE))

# 计算ASMI kg/m2
ASM_Height_long <- merge(ASM, Height_long, by.x = c("CODE","Followup"), by.y = c("ID","Followup"))
ASM_Height_long$ASMI <- (ASM_Height_long$ASM/1000)/(ASM_Height_long$Height/100)/(ASM_Height_long$Height/100)
colnames(ASM_Height_long) <- dplyr::recode(colnames(ASM_Height_long),
                                    "CODE" = "ID")
# 合并年龄
ASM_Height_long <- merge(ASM_Height_long, Age_long, by.x = c("ID","Followup"), by.y = c("ID","Times"))
ASM_Height_long$WBTOT_FAT <- ASM_Height_long$WBTOT_FAT*0.001
ASM_Height_long <- ASM_Height_long[complete.cases(ASM_Height_long),]
table(ASM_Height_long$Followup)
####骨骼肌轨迹####
data_trac <- merge(ASM_Height_long,Cov_F0[,c("ID","Sex_F0")], by = "ID")
data_trac <- data_trac[complete.cases(data_trac),]
data_trac$ID_num <- as.numeric(sub("^NL", "", data_trac$ID))
data_trac$Times <- as.numeric(sub("^F", "", data_trac$Followup))
data_trac_female <- data_trac[data_trac$Sex_F0 == 0,]
data_trac_male <- data_trac[data_trac$Sex_F0 == 1,]
#********************************函数
run_hlme_traj_plot <- function(
    data,
    max_ng = 5,
    ng_plot = 3,
    sex_label = "Female",
    legend.x = legend.x,
    legend.y = legend.y,
    pdf_file = NULL,
    width = 6,
    height = 6,
    colors = c("#1f77b4","#ff7f0e","#2ca02c")
){
  
  library(lcmm)
  library(dplyr)
  
  
  ## 固定变量名（不作为参数）
  outcome  <- "ASMI"
  time_var <- "Times"
  id_var   <- "ID_num"
  
  ## 1. ng = 1 基线模型
  model_list <- list()
  
  model_list[[1]] <- hlme(
    ASMI ~ Times,
    subject = id_var,
    ng = 1,
    data = data
  )
  
  ## 2. ng = 2 ~ max_ng
  if (max_ng >= 2) {
    for (k in 2:max_ng) {
      model_list[[k]] <- hlme(
        ASMI ~ Times,
        mixture = ~ Times,
        subject = id_var,
        ng = k,
        data = data,
        B = model_list[[1]]
      )
    }
  }
  
  ## 3. BIC + %class 汇总
  summary_tab <- do.call(
    summarytable,
    c(model_list, list(which = c("BIC", "%class")))
  )
  
  ## 4. 选定画图模型
  model_plot <- model_list[[ng_plot]]
  
  ## 后验概率
  pprob <- model_plot$pprob
  prob_cols <- paste0("prob", 1:ng_plot)
  
  # 获得每一个ID的分类
  Class <- as.data.frame(pprob)
  ID_unique <- data[,c("ID","ID_num")] %>%
    group_by(ID) %>%
    dplyr::slice(1) %>%
    ungroup()
  Class <- left_join(Class, ID_unique, by = "ID_num")
  
  ## class 占比
  class_percent <- round(colMeans(pprob[, prob_cols]) * 100, 1)
  class_labels <- paste0(
    "Class ", 1:ng_plot, " (", class_percent, "%)"
  )
  
  ## 5. 构造预测数据
  newdata <- data.frame(
    Times = seq(
      min(data$Times, na.rm = TRUE),
      max(data$Times, na.rm = TRUE),
      length = 15
    ),
    ASMI = mean(data$ASMI, na.rm = TRUE)
  )
  
  ## 6. 画轨迹图
  if (!is.null(pdf_file)) {
    pdf(pdf_file, width = width, height = height)
  }
  
  # 自定义横坐标刻度：1,2,3,4
  x_ticks <- 1:4
  
  plot_obj <- plot(
    predictY(model_plot, newdata, var.time = "Times", draws = TRUE),
    ylab = "ASMI (kg/m²)",
    xlab = "Times",
    main = "",
    cex = 0.8,
    lwd = 2,
    shade = TRUE,
    legend.loc = FALSE,
    col = colors[1:ng_plot],
    font.lab = 2,
    cex.lab = 1.2,
    xaxt = "n"  # 不显示默认横坐标
  )
  
  mtext(
    sex_label,
    side = 3,     # 上边
    adj  = 0,     # 左对齐（关键）
    line = 0.5,   # 距离图框的高度，可微调
    font = 2,      # 加粗（可选）
    cex = 1.4
  )
  
  # 添加自定义横坐标
  axis(
    side = 1,
    at = x_ticks,
    labels = x_ticks
  )
  
  legend(
    x = legend.x, 
    y = legend.y, 
    legend = class_labels,
    col = colors[1:ng_plot],
    lty = 1:ng_plot,
    lwd = 2,
    bty = "n"
  )
  
  if (!is.null(pdf_file)) {
    dev.off()
  }
  
  
  ## 7. 返回结果
  return(list(
    model_list    = model_list,
    summary_table = summary_tab,
    model_plot    = model_plot,
    class_percent = class_percent,
    plot          = plot_obj,
    Class         = Class
  ))
}

# pdf("Figure_ASMI_Trajectory_Female_Male.pdf", width = 14, height = 6)
# par(
#   mfrow = c(1, 2),          # 一行两列
#   mar = c(5, 5, 3, 2)       # 下左上右边距，防止挤
# )
# 
# data_trac_female <- data_trac[data_trac$Sex_F0 == 0,]
# run_hlme_traj_plot(
#   data = data_trac_female,
#   max_ng = 5,
#   ng_plot = 3,
#   sex_label = "Female",
#   legend.x = 2.8,
#   legend.y = 7.0,
#   pdf_file = NULL,
#   colors = c("#1f77b4","#ff7f0e","#2ca02c")
# )
# 
# data_trac_male <- data_trac[data_trac$Sex_F0 == 1,]
# run_hlme_traj_plot(
#   data = data_trac_male,
#   max_ng = 5,
#   ng_plot = 3,
#   sex_label = "Male",
#   legend.x = 2.8,
#   legend.y = 8.2,
#   pdf_file = NULL,
#   colors = c( "#8c564b","#9467bd","#d62728")
# )
# 
# dev.off()

Class_female <- run_hlme_traj_plot(
  data = data_trac_female,
  max_ng = 5,
  ng_plot = 3,
  sex_label = "Female",
  legend.x = 2.8,
  legend.y = 7.0,
  pdf_file = "Female_tra.pdf",
  colors = c("#1f77b4","#ff7f0e","#2ca02c")
)
Class_female$Class$class <- as.factor(Class_female$Class$class)

Class_male <- run_hlme_traj_plot(
  data = data_trac_male,
  max_ng = 5,
  ng_plot = 3,
  sex_label = "Male",
  legend.x = 2.8,
  legend.y = 8.2,
  pdf_file = "Male_tra.pdf",
  colors = c( "#8c564b","#9467bd","#d62728")
)
Class_male$Class$class <- as.factor(Class_male$Class$class)

####验证轨迹####
# data_trac_female <- data_trac[data_trac$Sex_F0 == 0,]
# class(data_trac_female$ID_num)
# model1 <- hlme(ASMI~Times, subject="ID_num",ng=1,data=data_trac_female)
# model2 <- hlme(ASMI~Times, mixture= ~Times, subject="ID_num",ng=2, data=data_trac_female, B=model1)
# model3 <- hlme(ASMI~Times, mixture= ~Times, subject="ID_num",ng=3, data=data_trac_female, B=model1)
# model4 <- hlme(ASMI~Times, mixture= ~Times, subject="ID_num",ng=4, data=data_trac_female, B=model1)
# model5 <- hlme(ASMI~Times, mixture= ~Times, subject="ID_num",ng=5, data=data_trac_female, B=model1)
# summarytable(model1,model2, model3, model4, model5, which=c("BIC","%class"))
# postprob(model3)
# 
# Class_female <- as.data.frame(model3$pprob)
# ID_unique <- data_trac_female[,c("ID","ID_num")] %>%
#   group_by(ID) %>%
#   slice(1) %>%
#   ungroup()
# Class_female <- left_join(Class_female, ID_unique, by = "ID_num")
# # 画图
# newdata <- data.frame(Times=seq(0, 4,length=15), ASMI=rep(10, 15))
# pprob <- model3$pprob
# class_percent <- round(colMeans(pprob[,c("prob1","prob2","prob3")]) * 100, 1)
# class_labels <- paste0("Class ", 1:ncol(pprob[,c("prob1","prob2","prob3")]), " (", class_percent, "%)")
# 
# pdf("Female_trac.pdf",width = 8, height = 6)
# plot(
#   predictY(model3, newdata, var.time="Times", draws=TRUE),
#   ylab = "ASMI (kg/m²)",
#   xlab = "Times",
#   main = "Female",
#   cex = 0.8,          # 点大小
#   lwd = 2,            # 线宽
#   shade = TRUE,       # 阴影
#   legend.loc = FALSE,
#   col = c("#1f77b4","#ff7f0e","#2ca02c"), # 自定义颜色
#   font.lab = 2,       # 坐标轴标题加粗
#   cex.lab = 1.2       # 坐标轴标题大小
# )
# 
# # 添加自定义 legend（带占比）
# legend(
#   #"topright",
#   x = 2.8, y = 7.0,
#   legend = class_labels,
#   col = c("#1f77b4","#ff7f0e","#2ca02c"),
#   lty = 1:3,
#   lwd = 2,
#   bty = "n"  # 去掉边框
# )
# dev.off()
####DXA宽数据####
ASMI_wide <- ASM_Height_long %>%
  pivot_wider(
    id_cols = ID,  
    names_from = Followup,   
    values_from = c(ASM, ASMI, WBTOT_FAT),  
    names_glue = "{Followup}{.value}", 
    values_fill = NA 
  )
ASMI_wide$WBTOT_FAT_mean_F1234 <- rowMeans(ASMI_wide[, c("F1WBTOT_FAT","F2WBTOT_FAT","F3WBTOT_FAT","F4WBTOT_FAT")], na.rm = TRUE)


ASM_wide_Age <- merge(Age_wide[, c("ID","F1age","F2age","F3age","F4age")], ASMI_wide, by = "ID")
# 只保留F1-F4都不缺失ASM的观测
ASM_wide_Age_filter_F14 <- drop_na(ASM_wide_Age)
ASM_wide_Age_filter_F14$WBTOT_FAT_mean <- rowMeans(ASM_wide_Age_filter_F14[, c("F1WBTOT_FAT","F2WBTOT_FAT","F3WBTOT_FAT","F4WBTOT_FAT")], na.rm = TRUE)
# 只保留F2-F4都不缺失ASM的观测
ASM_wide_Age_F24 <- ASM_wide_Age[,c("ID","F2age","F3age","F4age","F2ASM","F3ASM","F4ASM","F2ASMI","F3ASMI","F4ASMI","F2WBTOT_FAT","F3WBTOT_FAT","F4WBTOT_FAT")]
ASM_wide_Age_filter_F24 <- drop_na(ASM_wide_Age_F24)
ASM_wide_Age_filter_F24$WBTOT_FAT_mean <- rowMeans(ASM_wide_Age_filter_F24[, c("F2WBTOT_FAT","F3WBTOT_FAT","F4WBTOT_FAT")], na.rm = TRUE)
####骨骼肌AUC和slope####
vars <- c("F1ASMI","F2ASMI","F3ASMI","F4ASMI","F1ASM","F2ASM","F3ASM","F4ASM","F1WBTOT_FAT","F2WBTOT_FAT","F3WBTOT_FAT","F4WBTOT_FAT")
base_names <- unique(sub("^F[0-4]", "", vars))

Muscle_AUC_slope_F14 <- calculate_auc_slope(
  prefixes = c("F1", "F2", "F3", "F4"),
  data = ASM_wide_Age_filter_F14,
  base_names = base_names
)
Muscle_AUC_slope_F14 <- Muscle_AUC_slope_F14 %>%
  rename_with(~ ifelse(.x == "ID", .x, paste0(.x, "_F14")))

Muscle_AUC_slope_F24 <- calculate_auc_slope(
  prefixes = c("F2", "F3", "F4"),
  data = ASM_wide_Age_filter_F24,
  base_names = base_names
)
Muscle_AUC_slope_F24 <- Muscle_AUC_slope_F24 %>%
  rename_with(~ ifelse(.x == "ID", .x, paste0(.x, "_F24")))
####*****************肌少症评定****************####
#（1）ASM（DXA）:男<7.0 kg/m2;女<5.4 kg/m2。（2）握力：男<28 kg；女<18 kg。（3）5次坐立试验：≥12 s
SPPB1 <- read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/SPPB/20251206_GNHS_F0-3.xlsx")
colnames(SPPB1)
SPPB2 <- read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/SPPB/20251206_GNHS_F4.xlsx")
SPPB2$ID <- substr(SPPB2$CODE, 3, nchar(SPPB2$CODE))
# 判定函数
compute_sarcopenia_F <- function(SPPB_data, 
                                 grip1_col, grip2_col, chair_col,
                                 followup_time, output_col) {
  
  # 1️⃣ 数据合并
  data_list <- list(
    Cov_F0[, c("ID", "Sex_F0")],
    SPPB_data[, c("ID", grip1_col, grip2_col, chair_col)],
    ASM_Height_long[ASM_Height_long$Followup == followup_time, c("ID", "ASMI")]
  )
  
  SPPB_merged <- purrr::reduce(data_list, dplyr::inner_join, by = "ID")
  
  # 2️⃣ 肌少症判定
  SPPB_processed <- SPPB_merged %>%
    mutate(
      # 平均握力
      grip_mean = rowMeans(dplyr::select(., all_of(c(grip1_col, grip2_col))), na.rm = TRUE),
      
      # 低肌肉量
      low_ASMI = case_when(
        Sex_F0 == 1 & ASMI < 7.0 ~ TRUE,
        Sex_F0 == 0 & ASMI < 5.4 ~ TRUE,
        TRUE ~ FALSE
      ),
      
      # 握力低
      low_grip = case_when(
        Sex_F0 == 1 & grip_mean < 28 ~ TRUE,
        Sex_F0 == 0 & grip_mean < 18 ~ TRUE,
        TRUE ~ FALSE
      ),
      
      # 坐凳试验慢
      slow_chair = ifelse(!is.na(.data[[chair_col]]) & .data[[chair_col]] >= 12, TRUE, FALSE),
      
      # 肌少症判定
      !!output_col := if_else((low_ASMI & low_grip) | (low_ASMI & low_grip & slow_chair), 1, 0)
    ) %>%
    dplyr::select(-grip_mean, -low_ASMI, -low_grip, -slow_chair)  # 删除中间变量
  
  return(SPPB_processed)
}

SPPB_F2 <- compute_sarcopenia_F(
  SPPB_data = SPPB1,
  grip1_col = "握力1_F2",
  grip2_col = "握力2_F2",
  chair_col = "重复坐凳试验时间_F2",
  followup_time = "F2",
  output_col = "sarcopenia_F2"
)

SPPB_F3 <- compute_sarcopenia_F(
  SPPB_data = SPPB1,
  grip1_col = "握力1_F3",
  grip2_col = "握力2_F3",
  chair_col = "重复坐凳试验时间_F3",
  followup_time = "F3",
  output_col = "sarcopenia_F3"
)
SPPB_F3 <- SPPB_F3[complete.cases(SPPB_F3),]
table(SPPB_F3$sarcopenia_F3)
SPPB_F3$sarcopenia_F3 <- as.factor(SPPB_F3$sarcopenia_F3)

SPPB_F4 <- compute_sarcopenia_F(
  SPPB_data = SPPB2,
  grip1_col = "握力1",
  grip2_col = "握力2",
  chair_col = "重复坐凳试验时间",
  followup_time = "F4",
  output_col = "sarcopenia_F4"
)
SPPB_F4 <- SPPB_F4[complete.cases(SPPB_F4),]
table(SPPB_F4$sarcopenia_F4)
SPPB_F4$sarcopenia_F4 <- as.factor(SPPB_F4$sarcopenia_F4)

data_list <- list(SPPB_F2[,c("ID","sarcopenia_F2")], 
                  SPPB_F3[,c("ID","sarcopenia_F3")], 
                  SPPB_F4[,c("ID","sarcopenia_F4")])
SPPB <- data_list %>% purrr::reduce(dplyr::left_join, by = "ID")
SPPB <- SPPB %>%
  mutate(
    sarcopenia = if_else(
      pmax(sarcopenia_F2, sarcopenia_F3, sarcopenia_F4, na.rm = TRUE) == 1,
      1, 0
    ),
    sarcopenia_F234 = if_else(
      sarcopenia_F2 == 1 & sarcopenia_F3 == 1 & sarcopenia_F4 == 1,
      1, 0
    ),
    sarcopenia_F34 = if_else(
      sarcopenia_F3 == 1 & sarcopenia_F4 == 1,
      1, 0
    )
  )

SPPB_complete <- SPPB[complete.cases(SPPB),]
table(SPPB_complete$sarcopenia)
table(SPPB_complete$sarcopenia_F234)
table(SPPB_complete$sarcopenia_F34)
####验证肌少症评定函数####
# data_list <- list(Cov_F0[,c("ID","Sex_F0")], 
#                   SPPB1[,c("ID","握力1_F2","握力2_F2","重复坐凳试验时间_F2")], 
#                   ASM_Height_long[ASM_Height_long$Followup == "F2",c("ID","ASMI")])
# SPPB_F2 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID")
# 
# SPPB_F2_processed <- SPPB_F2 %>%
#   mutate(
#     # 计算两次握力平均值
#     grip_mean = rowMeans(dplyr::select(., 握力1_F2, 握力2_F2), na.rm = TRUE),
#     
#     # 按性别判断 ASMI 是否低
#     low_ASMI = case_when(
#       Sex_F0 == 1 & ASMI < 7.0 ~ TRUE,   # 男
#       Sex_F0 == 0 & ASMI < 5.4 ~ TRUE,   # 女
#       TRUE ~ FALSE
#     ),
#     
#     # 按性别判断握力是否低
#     low_grip = case_when(
#       Sex_F0 == 1 & grip_mean < 28 ~ TRUE,  # 男
#       Sex_F0 == 0 & grip_mean < 18 ~ TRUE,  # 女
#       TRUE ~ FALSE
#     ),
#     
#     # 重复坐凳试验是否达标
#     slow_chair = ifelse(!is.na(重复坐凳试验时间_F2) & 重复坐凳试验时间_F2 >= 12, TRUE, FALSE),
#     
#     # 肌少症判定：同时满足(1)+(2) 或 (1)+(3)
#     sarcopenia_F2 = if_else(
#       (low_ASMI & low_grip) | (low_ASMI & slow_chair),
#       1, 0
#     )
#   )

####*********菌群*********####
####各水平####
Micro <- as.data.frame(read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_metaphlan4_taxa_3678_IDcorrect_251125_VF.xlsx",sheet=1))
Micro_all <- Micro %>%
  mutate(
    Times = substr(CODE, 1, 2),
    ID    = substr(CODE, 3, nchar(CODE))
  ) %>%
  relocate(ID, Times, .after = 1)
Micro_all_F3 <- Micro_all[Micro_all$Times == "F3",]
Micro_all_F3 <- dplyr::select(Micro_all_F3,-c("CODE","Times"))
# 筛选出菌种
Micro_selected <- Micro[, colnames(Micro)[colnames(Micro) == "CODE" | grepl("^s_", colnames(Micro))]]
Micro_selected <- Micro_selected %>%
  mutate(
    Times = substr(CODE, 1, 2),
    ID    = substr(CODE, 3, nchar(CODE))
  ) %>%
  relocate(ID, Times, .after = 1)
table(Micro_selected$Times)
Micro_selected <- dplyr::select(Micro_selected,-CODE)
#****************菌群过滤
Micro_selected_filter <- perform_micro_filter(Micro_selected,
                                              Micro_selected,
                                              exclude_cols = c("ID", "Times"),
                                              metadata_cols = c("ID", "Times"),
                                              min_abundance = 0.0001,  
                                              min_sample_detection = 0.1 
)
a <- Micro_selected_filter$filtered_data
#****************clr转换，因为clr是按行转换，所以跟整体分布无关
Micro_selected_filter_clr <- perform_micro_clr(Micro_selected_filter$filtered_data, exclude_cols = c("ID", "Times"))

####多样性####
Diversity <- as.data.frame(read.csv("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_metaphlan3_diversity_2664_221028_VF.csv",header = TRUE))
Diversity <- Diversity %>%
  mutate(
    Times = substr(id, 1, 2),
    ID    = substr(id, 3, nchar(id))
  ) %>%
  relocate(ID, Times, .after = 1)
####微生物功能####
Micro_func <- readr::read_csv(
  "/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_humann3_pathabundance_2659_unstra_230606_VF.tsv",
  locale = locale(encoding = "UTF-8")
)
Micro_func <- Micro_func %>%
  mutate(
    Times = substr(id, 1, 2),
    ID    = substr(id, 3, nchar(id))
  ) %>%
  relocate(ID, Times, .after = 1)
# 进行clr转换或过滤时要把一些变量去掉
Micro_func <- dplyr::select(Micro_func,- c("...1","id","batch","tag","UNMAPPED","UNINTEGRATED"))
#****************菌群过滤
Micro_func_filter <- perform_micro_filter(Micro_func,
                                          Micro_func,
                                          exclude_cols = c("ID", "Times"),
                                          metadata_cols = c("ID", "Times"),
                                          min_abundance = 0.0001,
                                          min_sample_detection = 0.1 
)
#****************clr转换，因为clr是按行转换，所以跟整体分布无关
Micro_func_filter_clr <- perform_micro_clr(Micro_func_filter$filtered_data, exclude_cols = c("ID", "Times"))
####*********血清蛋白质*********####
F0protein <- as.data.frame(read.csv("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_血清蛋白_基线3415人_VF.csv",header = TRUE))
F2protein <- as.data.frame(read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_血清蛋白F2_2567人_VF.xlsx",sheet = 1))
warnings()
F3protein1 <- as.data.frame(read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_F3蛋白组数据索取_更新_VF.xlsx",sheet = 1))
F3protein2 <- as.data.frame(read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_F3蛋白组数据索取_更新_VF.xlsx",sheet = 2))
colnames(F3protein2) <- sub("_.*", "", colnames(F3protein2))
#*************************F3蛋白质合并
phase1 <- as.data.frame(read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_F3蛋白组数据索取_更新_VF.xlsx",sheet = 3))
phase1 <- phase1[,c("patient_ID", "sample_collection_batch","phase")]
phase1 <- phase1[phase1$sample_collection_batch == "F3",]
phase1_unique <- phase1[!duplicated(phase1$patient_ID), ]
F3protein1_phase <- merge(phase1_unique, F3protein1, by.x = c("patient_ID", "sample_collection_batch"), by.y = c("sampleid","followup"))
colnames(F3protein1_phase) <- recode(colnames(F3protein1_phase),
                                     "patient_ID" = "ID",
                                     "sample_collection_batch" = "followup",
                                     "phase" = "Phase")
F3protein1_phase <- dplyr::select(F3protein1_phase,-c("time"))

phase2 <- as.data.frame(read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_F3蛋白组数据索取_更新_VF.xlsx",sheet = 4))
phase2$ID <- substr(phase2$ID, 1, 6)
phase2 <- phase2[,c("ID","time","phase")]
phase2 <- phase2[complete.cases(phase2),]
phase2 <- phase2[phase2$time == "F3",]
phase2_unique <- phase2[!duplicated(phase2$ID), ]
F3protein2_phase <- merge(phase2_unique, F3protein2, by.x = c("ID", "time"), by.y = c("sampleid","followup"))
colnames(F3protein2_phase) <- recode(colnames(F3protein2_phase),
                                     "time"  = "followup",
                                     "phase" = "Phase")
F3protein2_phase <- dplyr::select(F3protein2_phase,-c("id", "time.y"))

Val_intersect <- intersect(colnames(F3protein1_phase), colnames(F3protein2_phase))
F3protein_all <- rbind(F3protein1_phase[,Val_intersect], F3protein2_phase[,Val_intersect])
a <- intersect(intersect(colnames(F0protein),colnames(F2protein)),colnames(F3protein_all))
#*************************看看各个蛋白质的缺失值个数
missing_rate_F3 <- F3protein_all %>%
  dplyr::select(-ID, -followup, -Phase) %>%   # 排除非蛋白质列
  summarise(across(everything(), ~ mean(is.na(.)))) %>%  # 计算每列缺失率
  pivot_longer(cols = everything(),
               names_to = "Protein",
               values_to = "Missing_rate")%>%
  arrange(desc(Missing_rate))

missing_rate_F2 <- F2protein %>%
  dplyr::select(-ID, -followup, -Phase, -CODE) %>%   # 排除非蛋白质列
  summarise(across(everything(), ~ mean(is.na(.)))) %>%  # 计算每列缺失率
  pivot_longer(cols = everything(),
               names_to = "Protein",
               values_to = "Missing_rate")%>%
  arrange(desc(Missing_rate))

missing_rate_F0 <- F0protein %>%
  dplyr::select(-ID, -followup, -Phase, -CODE) %>%   # 排除非蛋白质列
  summarise(across(everything(), ~ mean(is.na(.)))) %>%  # 计算每列缺失率
  pivot_longer(cols = everything(),
               names_to = "Protein",
               values_to = "Missing_rate")%>%
  arrange(desc(Missing_rate))
#*************************蛋白质剔除和填补
process_protein_data <- function(F3protein_all, exclude_cols, detection_cutoff = 0.80) {
  
  # 获取蛋白质列（排除 exclude_cols 中的列）
  protein_cols <- setdiff(colnames(F3protein_all), exclude_cols)
  
  # 删除检测率低于 80% 样本数的蛋白质
  F3protein_all_filtered <- F3protein_all %>%
    mutate(across(all_of(protein_cols), ~ !is.na(.))) %>%  # 标记缺失值为FALSE
    summarise(across(all_of(protein_cols), ~ mean(.))) %>%  # 计算每个蛋白质的缺失率
    pivot_longer(cols = everything(), names_to = "protein", values_to = "detection_rate") %>%  # 转换为长格式
    filter(detection_rate >= detection_cutoff) %>%  # 筛选检测率 >= 80% 的蛋白质
    pull(protein)  # 获取满足条件的蛋白质名
  
  # 筛选掉检测率低于 80% 的蛋白质列
  F3protein_all_filtered <- F3protein_all %>%
    dplyr::select(all_of(c(exclude_cols, F3protein_all_filtered)))
  
  # 在进行任何操作前，先将蛋白质列转换为数值型
  F3protein_all_filtered_fill <- F3protein_all_filtered %>%
    mutate(across(all_of(setdiff(colnames(F3protein_all_filtered), exclude_cols)), 
                  as.numeric)) %>%
    # 填补缺失值，使用最小浓度的一半填补
    mutate(across(all_of(setdiff(colnames(F3protein_all_filtered), exclude_cols)), 
                  ~ {
                    # 计算非缺失值的最小值
                    min_val <- min(., na.rm = TRUE)
                    # 如果所有值都是NA，min_val会是Inf，需要处理这种情况
                    if (is.finite(min_val)) {
                      ifelse(is.na(.), min_val / 2, .)
                    } else {
                      .  # 如果所有值都是NA，保持原样
                    }
                  }))
  
  return(F3protein_all_filtered_fill)
}

F3protein_fill <- process_protein_data(F3protein_all, c("ID", "followup", "Phase"), 0.5)

F2protein[F2protein == "NA" & !is.na(F2protein)] <- NA
F2protein_fill <- process_protein_data(F2protein, c("CODE", "ID", "followup", "Phase"), 0.5)

F0protein_fill <- process_protein_data(F0protein, c("CODE", "ID", "followup", "Phase"), 0.5)
####蛋白质长数据####
Val_intersect <- intersect(intersect(colnames(F0protein_fill), colnames(F2protein_fill)),colnames(F3protein_fill))

Protein_long_all <- rbind(F3protein_fill[, Val_intersect], F2protein_fill[, Val_intersect], F0protein_fill[, Val_intersect])
Protein_long_all$followup <- recode(Protein_long_all$followup,
                                       "basline" = "F0")
Protein_long_all$followup <- recode(Protein_long_all$followup,
                                       "baseline" = "F0")
table(Protein_long_all$followup)
####蛋白质名称####
# Protein_names <- toupper(setdiff(colnames(Protein_long_all),c("ID", "followup", "Phase")))
# mapping_all <- getBM(
#   attributes = c(
#     "uniprotswissprot",
#     "hgnc_symbol",
#     "description"
#   ),
#   filters = "uniprotswissprot",
#   values = Protein_names,
#   mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# )
# colnames(mapping_all)
# Protein_names_output <- left_join(as.data.frame(Protein_names),mapping_all, by = c("Protein_names" = "uniprotswissprot"))
# # 找出有重复的Protein_names
# dup_proteins <- Protein_names_output %>%
#   group_by(Protein_names) %>%
#   filter(n() >= 2) %>%
#   ungroup()
# 
# a <- Protein_names_output[is.na(Protein_names_output$hgnc_symbol),]
# a$Protein_names
# # 只保留重复观测的第一条
# Protein_names_unique <- Protein_names_output %>%
#   distinct(Protein_names, .keep_all = TRUE)
# # 缺失的蛋白质hgnc手动补充
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "P69905"] <-  "HBA1"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "Q8IVJ8"] <-  "APRG1"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "P02746"] <-  "C1QB"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "P00736"] <-  "C1R"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "P0CG22"] <-  "DHRS4L1"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "Q8NFI4"] <-  "ST13P5"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "P26927"] <-  "MST1"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "P0DOX2"] <-  "IGHA2"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "P0DOX3"] <-  "IGHA1"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "P0DOX4"] <-  "IGHE"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "P0DOX5"] <-  "IGHG1"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "P01860"] <-  "IGHG3"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "P0DOX7"] <-  "IGKC"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "P01834"] <-  "IGKC"
# Protein_names_unique$hgnc_symbol[Protein_names_unique$Protein_names == "A6NGW2"] <-  "STRCP1"
# 
# writexl::write_xlsx(Protein_names_unique, "/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/9 饮茶频率&骨骼肌/Database/Protein_names_unique.xlsx")

# 导入蛋白质名称
Protein_names <- read_excel("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/9 饮茶频率&骨骼肌/Database/Protein_names_unique.xlsx",sheet=1)

####验证蛋白质剔除和填补####
# protein_cols <- setdiff(colnames(F3protein_all), c("ID", "followup", "Phase"))
# # 删除检测率低于 95% 样本数的蛋白质
# F3protein_all_filtered <- F3protein_all %>%
#   mutate(across(protein_cols, ~ !is.na(.))) %>%  # 标记缺失值为FALSE
#   summarise(across(protein_cols, ~ mean(.))) %>%  # 计算每个蛋白质的缺失率
#   gather(key = "protein", value = "detection_rate") %>%  # 转换为长格式
#   filter(detection_rate >= 0.80) %>%  # 筛选检测率 >= 80% 的蛋白质
#   pull(protein)  # 获取满足条件的蛋白质名
# # 筛选掉检测率低于 95% 的蛋白质列
# F3protein_all_filtered <- F3protein_all %>%
#   dplyr::select(c("ID", "Times", "phase", all_of(F3protein_all_filtered)))
# # min(F3protein_all_filtered$a0a075b6i0, na.rm = TRUE)
# # 填补缺失值，使用最小浓度的一半填补
# F3protein_all_filtered_fill <- F3protein_all_filtered %>%
#   mutate(across(setdiff(colnames(F3protein_all_filtered), c("ID", "Times", "phase")), 
#                 ~ ifelse(is.na(.), min(., na.rm = TRUE) / 2, .)))
# min(F3protein_all_filtered_fill$a0a075b6i0, na.rm = TRUE)

#*************************蛋白质ID找蛋白质名称
Dm_Acc <- read.delim("/Users/hozan/Library/CloudStorage/OneDrive-Personal/Papers/8 多组学&骨骼肌/重要文件/uniprotkb_organism_id_7227_AND_reviewed_2025_10_06.tsv", header = TRUE, sep = "\t")
####*********混杂定义*********####
Covariates_all = c("Age_F0","Sex_F0","WBTOT_FAT_mean_F1234","Met_mean","Energy_mean","Alcohol_F0","Smoke_F0","Calcium_F0","Vitamin_F0","Disease_F0","Estrogen_F0","Menopause_F0","Oil_sup_follow")
Covariates_female = c("Age_F0","Met_mean","WBTOT_FAT_mean_F1234","Energy_mean","Alcohol_F0","Calcium_F0","Vitamin_F0","Disease_F0","Estrogen_F0","Menopause_F0","Oil_sup_follow")
Covariates_male = c("Age_F0","Met_mean","WBTOT_FAT_mean_F1234","Energy_mean","Alcohol_F0","Smoke_F0","Calcium_F0","Vitamin_F0","Disease_F0","Oil_sup_follow")
# 线形混合效应混杂
Covariates_all_lmer = c("Age","Sex_F0","WBTOT_FAT","Met_mean","Energy_mean","Alcohol_F0","Smoke_F0","Calcium_F0","Vitamin_F0","Disease_F0","Estrogen_F0","Menopause_F0","Oil_sup_follow")
Covariates_female_lmer = c("Age","WBTOT_FAT","Met_mean","Energy_mean","Alcohol_F0","Calcium_F0","Vitamin_F0","Disease_F0","Estrogen_F0","Menopause_F0","Oil_sup_follow")
Covariates_male_lmer = c("Age","WBTOT_FAT","Met_mean","Energy_mean","Alcohol_F0","Smoke_F0","Calcium_F0","Vitamin_F0","Disease_F0","Oil_sup_follow")
# 不加能量摄入
Covariates_female_med = c("Age_F0","Met_mean","WBTOT_FAT","Alcohol_F0","Calcium_F0","Vitamin_F0","Disease_F0","Estrogen_F0","Menopause_F0","Oil_sup_follow")
Covariates_male_med = c("Age_F0","Met_mean","WBTOT_FAT","Alcohol_F0","Smoke_F0","Calcium_F0","Vitamin_F0","Disease_F0","Oil_sup_follow")
####儿茶素变量####
Catechin_vars <- c(
  "serum_I_catechin_F0",
  "serum_I_epicatechin_F0",
  "serum_I_EGC_F0",
  "serum_I_EGCG_F0",
  "serum_I_ECG_F0"
)

Catechin_vars3 <- c(
  "serum_I_catechin_F0",
  "serum_I_epicatechin_F0",
  "serum_I_EGC_F0"
)

Catechin_vars2 <- c(
  "serum_I_EGCG_F0",
  "serum_I_ECG_F0"
)

Catechin_vars_T <- c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
Catechin_vars_all <- c("Tea_freq","serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
####*****************************饮茶&骨骼肌-data1*****************************####
####数据集处理合并####
data_list <- list(Cov_F0, Tea_Follow, Catechin, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], ASM_Height_long)
data1_all <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))

data1_all_z <- scale_columns_group(data1_all, c("ASM","ASMI"), "Sex_F0")

data1_CatT4_all <- zero_tertile4(
  data = data1_all_z,
  vars = Catechin_vars3,
  suffix = "_T"
)

data1_CatT4_all <- data1_CatT4_all %>%
  mutate(
    across(
      all_of(Catechin_vars2),
      ~ quartile_group(.x),
      .names = "{.col}_T"
    )
  )

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(data1_CatT4_all)%in% factor_name)
for(i in idx ){
  data1_CatT4_all[[i]] <-  as.factor(data1_CatT4_all[[i]])
}

a1 <- unique(data1_CatT4_all$ID)
table(data1_CatT4_all$Followup)
####lmer####
Lmer_results1_all <- process_lmer(c("Tea_freq","serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T"),
                                     c("ASMI_z","ASM_z"), 
                                     data1_CatT4_all, 
                                     Covariates_all_lmer, 
                                     2)

####验证lmm####
# y_safe <- "ASMI_z"
# x_safe <- "serum_I_ECG_F0_T"
# data <- data1_CatT4_all
# cov_safe <- Covariates_all_lmer
# 
# formula_text <- paste(y_safe, "~",x_safe, "+", paste(cov_safe, collapse = " + "),"+ (1|ID)")
# 
# formula_full <- as.formula(formula_text)
# model <- lmer(formula_full, data = data)
# summary(model)
####正文写作####
a <- Lmer_results1_all[Lmer_results1_all$P_value < 0.05,]
a1 <- a[a$Predictor == "Tea_freq",]
round(max(a1$Estimate),3)
round(min(a1$Estimate),3)
round(max(a1$P_value),3)
round(min(a1$P_value),3)

a2 <- a[!a$Predictor %in% c("serum_I_catechin_F0_T", "Tea_freq"),]
round(max(a2$Estimate),3)
round(min(a2$Estimate),3)
round(max(a2$P_value),3)
round(min(a2$P_value),3)
####Lmer可视化####
level_order_tea <- c(
  "Non",
  "Less-than-daily",
  "Daily",
  "More-than-daily"
)

level_order_catechin <- c(
  "Undetectable",
  "Low",
  "Medium",
  "High"
)

level_order_quartile <- c(
  "Quartile 1",
  "Quartile 2",
  "Quartile 3",
  "Quartile 4"
)

df_plot <- Lmer_results1_all %>%
  mutate(
    Level_label = case_when(
      ## Tea frequency
      Predictor == "Tea_freq" & Level == 0 ~ "Non",
      Predictor == "Tea_freq" & Level == 1 ~ "Less-than-daily",
      Predictor == "Tea_freq" & Level == 2 ~ "Daily",
      Predictor == "Tea_freq" & Level == 3 ~ "More-than-daily",
      
      ## Catechins（0 = Undetectable）
      Predictor %in% c(
        "serum_I_catechin_F0_T",
        "serum_I_epicatechin_F0_T",
        "serum_I_EGC_F0_T"
      ) & Level == 0 ~ "Undetectable",
      Predictor %in% c(
        "serum_I_catechin_F0_T",
        "serum_I_epicatechin_F0_T",
        "serum_I_EGC_F0_T"
      ) & Level == 1 ~ "Low",
      Predictor %in% c(
        "serum_I_catechin_F0_T",
        "serum_I_epicatechin_F0_T",
        "serum_I_EGC_F0_T"
      ) & Level == 2 ~ "Medium",
      Predictor %in% c(
        "serum_I_catechin_F0_T",
        "serum_I_epicatechin_F0_T",
        "serum_I_EGC_F0_T"
      ) & Level == 3 ~ "High",
      
      ## EGCG / ECG（四分位）
      Predictor %in% c(
        "serum_I_EGCG_F0_T",
        "serum_I_ECG_F0_T"
      ) & Level == 0 ~ "Quartile 1",
      Predictor %in% c(
        "serum_I_EGCG_F0_T",
        "serum_I_ECG_F0_T"
      ) & Level == 1 ~ "Quartile 2",
      Predictor %in% c(
        "serum_I_EGCG_F0_T",
        "serum_I_ECG_F0_T"
      ) & Level == 2 ~ "Quartile 3",
      Predictor %in% c(
        "serum_I_EGCG_F0_T",
        "serum_I_ECG_F0_T"
      ) & Level == 3 ~ "Quartile 4",
      
      TRUE ~ NA_character_
    )
  ) %>%
  ## 按 Predictor 分别设 factor 顺序
  mutate(
    Level_label = case_when(
      Predictor == "Tea_freq" ~
        factor(Level_label, levels = level_order_tea),
      
      Predictor %in% c(
        "serum_I_catechin_F0_T",
        "serum_I_epicatechin_F0_T",
        "serum_I_EGC_F0_T"
      ) ~
        factor(Level_label, levels = level_order_catechin),
      
      Predictor %in% c(
        "serum_I_EGCG_F0_T",
        "serum_I_ECG_F0_T"
      ) ~
        factor(Level_label, levels = level_order_quartile),
      
      TRUE ~ factor(Level_label)
    )
  )


df_plot$Predictor <- recode(
  df_plot$Predictor,
  "Tea_freq" = "Tea consumption frequency",
  "serum_I_catechin_F0_T"    = "Catechin",
  "serum_I_epicatechin_F0_T" = "Epicatechin",
  "serum_I_EGC_F0_T"         = "Epigallocatechin",
  "serum_I_ECG_F0_T"         = "Epicatechin gallate",
  "serum_I_EGCG_F0_T"        = "Epigallocatechin gallate"
)

df_plot$Predictor <- factor(
  df_plot$Predictor,
  levels = c(
    "Tea consumption frequency",
    "Catechin",
    "Epicatechin",
    "Epigallocatechin",
    "Epicatechin gallate",
    "Epigallocatechin gallate"
  )
)

range_ASM  <- range(
  df_plot$Estimate[df_plot$Outcome == "ASM_z"],
  na.rm = TRUE
)

range_ASMI <- range(
  df_plot$Estimate[df_plot$Outcome == "ASMI_z"],
  na.rm = TRUE
)

scale_factor <- diff(range_ASM) / diff(range_ASMI)
scale_shift  <- range_ASM[1] - range_ASMI[1] * scale_factor

df_plot <- df_plot %>%
  mutate(
    ## 显著性
    Sig = ifelse(P_value < 0.05, "P < 0.05", "P ≥ 0.05"),
    
    ## ASMI 映射到 ASM 量级
    Estimate_plot = ifelse(
      Outcome == "ASMI_z",
      Estimate * scale_factor + scale_shift,
      Estimate
    ),
    CI_low_plot = ifelse(
      Outcome == "ASMI_z",
      CI_low * scale_factor + scale_shift,
      CI_low
    ),
    CI_high_plot = ifelse(
      Outcome == "ASMI_z",
      CI_high * scale_factor + scale_shift,
      CI_high
    ),
    
    ## 点位左右偏移
    x_nudge = ifelse(Outcome == "ASM_z", -0.15, 0.15)
  )

plot <- ggplot(
  df_plot,
  aes(
    x = Level_label,
    y = Estimate_plot,
    color = Outcome,
    shape = Sig,
    alpha = Sig
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  
  geom_errorbar(
    aes(
      ymin = CI_low_plot,
      ymax = CI_high_plot
    ),
    width = 0.15,
    linewidth = 0.8,
    position = position_dodge(width = 0.4)
  ) +
  
  geom_point(
    size = 3,
    position = position_dodge(width = 0.4)
  ) +
  
  facet_wrap(
    ~ Predictor,
    nrow = 2,
    scales = "free_x"
  ) +
  
  ## 双 y 轴
  scale_y_continuous(
    name = "ASM_z Estimate (95% CI)",
    sec.axis = sec_axis(
      ~ (. - scale_shift) / scale_factor,
      name = "ASMI_z Estimate (95% CI)"
    )
  ) +
  
  scale_color_manual(
    values = c("ASM_z" = "#3C5488", "ASMI_z" = "#E64B35")
  ) +
  scale_shape_manual(values = c("P < 0.05" = 16, "P ≥ 0.05" = 1)) +
  scale_alpha_manual(values = c("P < 0.05" = 1, "P ≥ 0.05" = 0.6)) +
  
  labs(
    x = "",
    color = "Outcome",
    shape = "Significance",
    alpha = "Significance"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    strip.text.x = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    axis.title = element_text(face = "bold", size = 14)
  )

pdf("Tea&ASM.pdf",height = 6, width = 12)
print(plot)
dev.off()
####Table 1####
data_first <- data1_CatT4_all %>%
  arrange(ID, Followup) %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup()

Table1a <- descrTable(Tea_freq ~ Age_F0+Sex_F0+
                        WBTOT_FAT_mean_F1234+
                        Met_mean+Energy_mean+
                        Smoke_F0+Alcohol_F0+Calcium_F0+Vitamin_F0+Oil_sup_follow+Fracture_F0+Disease_F0+Estrogen_F0+Menopause_F0,
                      data = data_first ,method =NA,show.all = FALSE)#show.all = TRUE :显示all
export2word(Table1a, file='Table1.docx')
####血清儿茶素类统计描述####
data_first_flo <- data_first %>%
  mutate(
    Tea_freq = factor(
      recode(
        Tea_freq,
        "0" = "Non",
        "1" = "Less-than-daily",
        "2" = "Daily",
        "3" = "More-than-daily"
      ),
      levels = c(
        "Non",
        "Less-than-daily",
        "Daily",
        "More-than-daily"
      )
    )
  ) %>%
  mutate(
    across(
      c(
        serum_I_catechin_F0_T,
        serum_I_epicatechin_F0_T,
        serum_I_EGC_F0_T
      ),
      ~ factor(
        recode(
          .x,
          "0" = "Undetectable",
          "1" = "Low",
          "2" = "Medium",
          "3" = "High"
        ),
        levels = c(
          "Undetectable",
          "Low",
          "Medium",
          "High"
        )
      )
    )
  ) %>%
  mutate(
    across(
      c(
        serum_I_EGCG_F0_T,
        serum_I_ECG_F0_T
      ),
      ~ factor(
        recode(
          .x,
          "0" = "Quartile 1",
          "1" = "Quartile 2",
          "2" = "Quartile 3",
          "3" = "Quartile 4"
        ),
        levels = c(
          "Quartile 1",
          "Quartile 2",
          "Quartile 3",
          "Quartile 4"
        )
      )
    )
  ) %>%
  rename(
    `Tea consumption frequency` = Tea_F0,
    Catechin = serum_I_catechin_F0_T,
    Epicatechin = serum_I_epicatechin_F0_T,
    Epigallocatechin = serum_I_EGC_F0_T,
    `Epigallocatechin gallate` = serum_I_EGCG_F0_T,
    `Epicatechin gallate` = serum_I_ECG_F0_T
  )

#**********************可视化
# 1. 提取并统计数据
target_vars <- c("Catechin", "Epicatechin", "Epigallocatechin", 
                 "Epigallocatechin gallate", "Epicatechin gallate")

summary_table <- data_first_flo %>%
  dplyr::select(all_of(target_vars)) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Category") %>%
  filter(!is.na(Category)) %>%
  group_by(Variable, Category) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Variable) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# 2. 【关键步骤】自定义变量排列顺序
# 在这里按你想要的顺序排列这五个变量
custom_order <- c(
  "Catechin", 
  "Epicatechin", 
  "Epigallocatechin", 
  "Epigallocatechin gallate",
  "Epicatechin gallate"
)

summary_table$Variable <- factor(summary_table$Variable, levels = custom_order)

# 3. 绘图
p <- ggplot(summary_table, aes(x = Category, y = Count, fill = Category)) +
  # 绘制带填充的条形图
  geom_bar(stat = "identity", width = 0.75, color = "white", linewidth = 0.3) +
  # 添加数值标签 (n 和 %)
  geom_text(
    aes(label = paste0(Count, "\n(", sprintf("%.1f", Percentage), "%)")),
    vjust = -0.3, 
    size = 3.2,
    lineheight = 0.9,
    fontface = "bold",
    color = "grey10"
  ) +
  # 分面排列：通过 ncol=3 或其他参数控制布局
  facet_wrap(~ Variable, scales = "free_x", ncol = 3) +
  # 使用学术感强的 Lancet 配色或 D3 配色
  scale_fill_lancet() +
  # 扩展 Y 轴顶部空间以容纳标签
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = NULL,
    x = NULL,
    y = "Number of Participants"
  ) +
  # 深度定制主题
  theme_minimal(base_size = 14) +
  theme(
    # 标题样式
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = ggplot2::margin(b=15)),
    # 分面标签美化
    strip.background = element_rect(fill = "#34495E", color = NA),
    strip.text = element_text(face = "bold", size = 15, color = "white"),
    # 坐标轴美化
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", color = "grey20"),
    axis.text.y = element_text(color = "grey30"),
    axis.title.y = element_text(face = "bold", size = 15, margin = ggplot2::margin(r=10)),
    # 网格线淡化
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linetype = "dotted"),
    # 图例
    legend.position = "none",
    panel.spacing = unit(1.5, "lines"),
    plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )

# 4. 输出预览与保存
print(p)
ggsave("Custom_Catechin_Distribution.pdf", p, width = 14, height = 9)
####基线饮茶习惯&儿茶素浓度####
Table1A_all <- descrTable(Tea_freq ~ Catechin + Epicatechin + Epigallocatechin + `Epigallocatechin gallate` + `Epicatechin gallate`,
                          data = data_first_flo, method =NA,show.all = FALSE) # show.all = TRUE :显示all
export2word(Table1A_all, file='Table2.docx')
####*****************************饮茶&蛋白质&骨骼肌-data2*****************************####
####数据集处理合并####
Protein_Age_all <- merge(Protein_long_all, Age_long_all, by.x = c("ID","followup"), by.y = c("ID","Times"))
data_list <- list(Cov_F0, Tea_Follow, Catechin, Protein_Age_all,ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")]) 
data2_all <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))

data2_CatT4_all <- zero_tertile4(
  data = data2_all,
  vars = Catechin_vars3,
  suffix = "_T"
)

data2_CatT4_all <- data2_CatT4_all %>%
  mutate(
    across(
      all_of(Catechin_vars2),
      ~ quartile_group(.x),
      .names = "{.col}_T"
    )
  )

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(data2_CatT4_all)%in% factor_name)
for(i in idx ){
  data2_CatT4_all[[i]] <-  as.factor(data2_CatT4_all[[i]])
}

table(data2_CatT4_all$Tea_freq)
####蛋白质数据转换####
data2_CatT4_all_z <- prepare_metabolite_data2(data2_CatT4_all, Protein_long_all, c("ID","followup","Phase"))
table(data2_CatT4_all_z$Phase)
a <- unique(data2_CatT4_all_z$ID)
####Lmer####
Lmer_results2_all <- process_lmer(c("Tea_freq","serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T"),
                                     paste0(setdiff(colnames(Protein_long_all),c("ID","followup","Phase")),"_z"), 
                                     data2_CatT4_all_z, 
                                     c("Phase",setdiff(Covariates_all_lmer,"WBTOT_FAT")),#"WBTOT_FAT_mean_F1234",
                                     2)
a <- unique(Lmer_results2_all$Outcome)

Lmer_results2_all_level <- Lmer_results2_all[Lmer_results2_all$Level == 3,]
Lmer_results2_all_level <- Lmer_results2_all_level %>%
  group_by(Predictor) %>%
  mutate(
    P_FDR = p.adjust(P_value, method = "fdr")
  ) %>%
  ungroup()

Lmer_results2_all_level_sig <- Lmer_results2_all_level[Lmer_results2_all_level$P_value < 0.05 & Lmer_results2_all_level$P_FDR < 0.05,]
Lmer_results2_all_level_sig$Outcome_name <- sub("_z", "", Lmer_results2_all_level_sig$Outcome)
Lmer_results2_all_level_sig$Outcome_name <- toupper(Lmer_results2_all_level_sig$Outcome_name)
Lmer_results2_all_level_sig <- merge(Lmer_results2_all_level_sig, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Outcome_name", by.y = "Protein_names")

unique(Lmer_results2_all_level_sig$Predictor)
Lmer_results2_all_level_sig_Tea <- Lmer_results2_all_level_sig[Lmer_results2_all_level_sig$Predictor == "Tea_freq",]
Lmer_results2_all_level_sig_epicatechin <- Lmer_results2_all_level_sig[Lmer_results2_all_level_sig$Predictor == "serum_I_epicatechin_F0_T",]
Lmer_results2_all_level_sig_ECG <- Lmer_results2_all_level_sig[Lmer_results2_all_level_sig$Predictor == "serum_I_ECG_F0_T",]
Lmer_results2_all_level_sig_EGCG <- Lmer_results2_all_level_sig[Lmer_results2_all_level_sig$Predictor == "serum_I_EGCG_F0_T",]
#Lmer_results2_all_level_sig_catechin <- Lmer_results2_all_level_sig2[Lmer_results2_all_level_sig2$Predictor == "serum_I_catechin_F0_T",]
#Lmer_results2_all_level_sig_EGC <- Lmer_results2_all_level_sig[Lmer_results2_all_level_sig$Predictor == "serum_I_EGC_F0_T",]
####饮茶&蛋白质可视化####
plot_circos_level <- function(
    df_input,      
    gender = "All",  
    level = "3",        
    file = "circos_plot.pdf",  
    width = 6, 
    height = 6,
    total_degree = 90,   
    gap_degree = 0.2,
    show_labels = TRUE    
) {
  
  library(dplyr)
  library(circlize)
  library(stringr)
  
  # 1. 数据整理
  df <- df_input %>%
    filter(Gender == gender, Level == level) %>%
    transmute(
      protein_raw = hgnc_symbol, 
      beta        = Estimate,
      p_fdr       = P_FDR,
      sig = case_when(
        p_fdr < 0.001 ~ "***",
        p_fdr < 0.01  ~ "**",
        p_fdr < 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  
  n_protein <- nrow(df)
  if(n_protein == 0) return(message("No data found for the specified filters."))
  
  # 处理重复名称
  df$protein_unique <- make.unique(as.character(df$protein_raw))
  df$protein_unique <- factor(df$protein_unique, levels = df$protein_unique)
  
  # 2. PDF 输出
  pdf(file, width = width, height = height)
  circos.clear()
  
  # ===== 【核心修改：颜色映射逻辑】
  # 自动获取最大绝对值，确保 0 点处于颜色条的正中心，颜色对比更均匀
  max_val <- max(abs(df$beta), na.rm = TRUE)
  
  col_fun <- colorRamp2(
    c(-max_val, 0, max_val),           # 定义三个关键点：最小值、中值(0)、最大值
    c("#4575B4", "white", "#D73027")  # 分别对应：蓝色(负)、白色(0)、红色(正)
  )
  
  # 计算角度逻辑
  big_gap <- 360 - total_degree
  
  circos.par(
    start.degree = 90,           
    gap.after = c(rep(gap_degree, n_protein - 1), big_gap), 
    track.margin = c(0.002, 0.002),
    cell.padding = c(0, 0, 0, 0),
    clock.wise = TRUE,
    canvas.xlim = c(-1.3, 1.3),  
    canvas.ylim = c(-1.3, 1.3)
  )
  
  # 3. 初始化
  circos.initialize(
    factors = df$protein_unique,
    xlim = cbind(rep(0, n_protein), rep(1, n_protein))
  )
  
  # 4. 绘制轨道
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.20,
    bg.border = NA,
    panel.fun = function(x, y) {
      p_unique <- CELL_META$sector.index
      i <- find_row_index <- which(df$protein_unique == p_unique)
      
      # 绘制热图块
      circos.rect(
        0, 0, 1, 1,
        col = col_fun(df$beta[i]), # 这里会根据 beta 正负自动选择红/蓝
        border = "white"
      )
      
      # 绘制显著性星号
      if (df$sig[i] != "") {
        circos.text(0.5, 0.5, df$sig[i], cex = 1.3, font = 2)
      }
      
      # 控制标签显示
      if (show_labels) {
        circos.text(
          CELL_META$xcenter,
          CELL_META$ylim[2] + mm_y(2.5),
          df$protein_raw[i], 
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0, 0.5),
          cex = 1.1
        )
      }
    }
  )
  
  circos.clear()
  dev.off()
  
  message("Circos plot saved. Red: Estimate > 0, Blue: Estimate < 0.")
}

plot_circos_level(
  df_input = Lmer_results2_all_level_sig_Tea,
  level = "3",
  file = "Tea_circos.pdf",
  total_degree = 360,   
  gap_degree = 0.2,       # 变量之间的微小间距
  show_labels = TRUE
)

plot_circos_level(
  df_input = Lmer_results2_all_level_sig_ECG,
  level = "3",
  file = "ECG_circos.pdf",
  total_degree = 360,   
  gap_degree = 0.2      
)

plot_circos_level(
  df_input = rbind(Lmer_results2_all_level_sig_epicatechin,  Lmer_results2_all_level_sig_EGCG),
  level = "3",
  file = "epicatechin_circos.pdf",
  total_degree = 360,   
  gap_degree = 0.2      
)
####正文写作####
a1 <- Lmer_results2_all_level_sig_Tea
a11 <- a1[a1$Estimate >0,]
a11 <- a11[order(a11$hgnc_symbol),]
a11$hgnc_symbol
round(max(a11$P_FDR),3)
round(min(a11$P_FDR),3)
round(max(a11$Estimate),3)
round(min(a11$Estimate),3)

a12 <- a1[a1$Estimate <0,]
a12 <- a12[order(a12$hgnc_symbol),]
a12$hgnc_symbol
round(max(a12$P_FDR),3)
round(min(a12$P_FDR),3)
round(max(a12$Estimate),3)
round(min(a12$Estimate),3)

a2 <- Lmer_results2_all_level_sig_epicatechin
a21 <- a2[a2$Estimate >0,]
a21 <- a21[order(a21$hgnc_symbol),]
a21$hgnc_symbol
round(max(a21$P_FDR),3)
round(min(a21$P_FDR),3)
round(max(a21$Estimate),3)
round(min(a21$Estimate),3)

a22 <- a2[a2$Estimate <0,]
a22 <- a22[order(a22$hgnc_symbol),]
a22$hgnc_symbol
round(max(a22$P_FDR),3)
round(min(a22$P_FDR),3)
round(max(a22$Estimate),3)
round(min(a22$Estimate),3)

a3 <- Lmer_results2_all_level_sig_EGCG
a31 <- a3[a3$Estimate >0,]
a31 <- a31[order(a31$hgnc_symbol),]
a31$hgnc_symbol
round(max(a31$P_FDR),3)
round(min(a31$P_FDR),3)
round(max(a31$Estimate),3)
round(min(a31$Estimate),3)

a32 <- a3[a3$Estimate <0,]
a32 <- a32[order(a32$hgnc_symbol),]
a32$hgnc_symbol
round(max(a32$P_FDR),3)
round(min(a32$P_FDR),3)
round(max(a32$Estimate),3)
round(min(a32$Estimate),3)

a4 <- Lmer_results2_all_level_sig_ECG
a41 <- a4[a4$Estimate >0,]
a41 <- a41[order(a41$hgnc_symbol),]
a41$hgnc_symbol
round(max(a41$P_FDR),3)
round(min(a41$P_FDR),3)
round(max(a41$Estimate),3)
round(min(a41$Estimate),3)

a42 <- a4[a4$Estimate <0,]
a42 <- a42[order(a42$hgnc_symbol),]
a42$hgnc_symbol
round(max(a42$P_FDR),3)
round(min(a42$P_FDR),3)
round(max(a42$Estimate),3)
round(min(a42$Estimate),3)
####导出附表####
Lmer_protein_sm <- Lmer_results2_all_level
Lmer_protein_sm$Outcome_name <- sub("_z", "", Lmer_protein_sm$Outcome)
Lmer_protein_sm$Outcome_name <- toupper(Lmer_protein_sm$Outcome_name)

Lmer_protein_sm2 <- merge(Lmer_protein_sm, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Outcome_name", by.y = "Protein_names")


Lmer_protein_sm2_output <- Lmer_protein_sm2[,c("Outcome_name","hgnc_symbol", "Predictor", "Level", "Estimate", "Std. Error", "P_value", "P_FDR")]
colnames(Lmer_protein_sm2_output) <- recode(colnames(Lmer_protein_sm2_output),
                                            "Outcome_name" = "Outcome",
                                            "hgnc_symbol" = "HGNC")
unique(Lmer_protein_sm2_output$Predictor)
Lmer_protein_sm2_output$Predictor <- recode(Lmer_protein_sm2_output$Predictor,
                                            "Tea_freq" = "Tea consumption frequency",
                                            "serum_I_catechin_F0_T" = "Catechin",
                                            "serum_I_epicatechin_F0_T" = "Epicatechin",
                                            "serum_I_EGC_F0_T" = "Epigallocatechin",
                                            "serum_I_EGCG_F0_T" = "Epigallocatechin gallate",
                                            "serum_I_ECG_F0_T" = "Epicatechin gallate")

Lmer_protein_sm2_output$Level[Lmer_protein_sm2_output$Predictor == "Tea consumption frequency"] <-  "More-than-daily vs. Non"
Lmer_protein_sm2_output$Level[Lmer_protein_sm2_output$Predictor %in% c("Catechin","Epicatechin","Epigallocatechin")] <-  "High vs. Undetectable"
Lmer_protein_sm2_output$Level[Lmer_protein_sm2_output$Predictor %in% c("Epigallocatechin gallate","Epicatechin gallate")] <-  "Quartile 4 vs. Quartile 1"


# 创建一个工作簿
wb <- createWorkbook()

# 添加工作表
addWorksheet(wb, "Significant_Results")

# 写入数据
writeData(wb, sheet = 1, Lmer_protein_sm2_output, rowNames = FALSE)

# 找出P_FDR < 0.05的行
significant_rows <- which(Lmer_protein_sm2_output$P_FDR < 0.05) + 1  # +1 因为有标题行
p_fdr_col <- which(names(Lmer_protein_sm2_output) == "P_FDR")

# 创建红色填充样式
red_style <- createStyle(fgFill = "#FFC7CE", fontColour = "#9C0006")  # 浅红背景，深红文字

# 对显著的P_FDR单元格应用红色样式
if(length(significant_rows) > 0) {
  addStyle(wb, 
           sheet = 1, 
           style = red_style, 
           rows = significant_rows, 
           cols = p_fdr_col,
           gridExpand = FALSE)
}

# 可选：添加条件格式规则（更动态）
conditionalFormatting(wb, 
                      sheet = 1,
                      cols = p_fdr_col,
                      rows = 2:(nrow(Lmer_protein_sm2_output)+1),  # 从第2行开始
                      rule = "<0.05",
                      style = red_style)

# 保存Excel文件
saveWorkbook(wb, "Supplemental Table2.xlsx", overwrite = TRUE)

####验证lmer####
# y_safe <- paste0(setdiff(colnames(Protein_long_female),c("ID","followup","Phase")),"_z")[100]
# x_safe <- "serum_I_ECG_F0_T"
# cov_safe <- c("Phase","WBTOT_FAT_mean_F1234",setdiff(Covariates_all_lmer,"WBTOT_FAT"))
# data <- data2_CatT4_all_z
# formula <- as.formula(
#   paste(
#     y_safe, "~",
#     x_safe, "+",
#     paste(cov_safe, collapse = "+"),
#     "+ (1|ID)"
#   )
# )
# model <- lmer(formula, data = data,
#               control = lmerControl(
#                 optimizer = "bobyqa",
#                 optCtrl = list(maxfun = 1e5)
#               ))

####计算GMSI####
Pro_GMSI_all_tea <- calculate_GMSI(
  data = data2_CatT4_all,
  outcomes_pos = sub("_z","", Lmer_results2_all_level_sig_Tea[Lmer_results2_all_level_sig_Tea$Estimate > 0,]$Outcome),
  outcomes_neg = sub("_z","", Lmer_results2_all_level_sig_Tea[Lmer_results2_all_level_sig_Tea$Estimate < 0,]$Outcome),
  "PI_Tea",
  time_var = "followup"
)

# Pro_GMSI_all_catechin <- calculate_GMSI(
#   data = data2_CatT4_all,
#   outcomes_pos = sub("_z","", Lmer_results2_all_level_sig_catechin[Lmer_results2_all_level_sig_catechin$Estimate > 0,]$Outcome),
#   outcomes_neg = sub("_z","", Lmer_results2_all_level_sig_catechin[Lmer_results2_all_level_sig_catechin$Estimate < 0,]$Outcome),
#   "PI_catechin",
#   time_var = "followup"
# )

Pro_GMSI_all_epicatechin <- calculate_GMSI(
  data = data2_CatT4_all,
  outcomes_pos = sub("_z","", Lmer_results2_all_level_sig_epicatechin[Lmer_results2_all_level_sig_epicatechin$Estimate > 0,]$Outcome),
  outcomes_neg = sub("_z","", Lmer_results2_all_level_sig_epicatechin[Lmer_results2_all_level_sig_epicatechin$Estimate < 0,]$Outcome),
  "PI_Epicatechin",
  time_var = "followup"
)

# Pro_GMSI_all_EGC <- calculate_GMSI(
#   data = data2_CatT4_all,
#   outcomes_pos = sub("_z","", Lmer_results2_all_level_sig_EGC[Lmer_results2_all_level_sig_EGC$Estimate > 0,]$Outcome),
#   outcomes_neg = sub("_z","", Lmer_results2_all_level_sig_EGC[Lmer_results2_all_level_sig_EGC$Estimate < 0,]$Outcome),
#   "PI_EGC",
#   time_var = "followup"
# )

Pro_GMSI_all_EGCG <- calculate_GMSI(
  data = data2_CatT4_all,
  outcomes_pos = sub("_z","", Lmer_results2_all_level_sig_EGCG[Lmer_results2_all_level_sig_EGCG$Estimate > 0,]$Outcome),
  outcomes_neg = sub("_z","", Lmer_results2_all_level_sig_EGCG[Lmer_results2_all_level_sig_EGCG$Estimate < 0,]$Outcome),
  "PI_EGCG",
  time_var = "followup"
)

Pro_GMSI_all_ECG <- calculate_GMSI(
  data = data2_CatT4_all,
  outcomes_pos = sub("_z","", Lmer_results2_all_level_sig_ECG[Lmer_results2_all_level_sig_ECG$Estimate > 0,]$Outcome),
  outcomes_neg = sub("_z","", Lmer_results2_all_level_sig_ECG[Lmer_results2_all_level_sig_ECG$Estimate < 0,]$Outcome),
  "PI_ECG",
  time_var = "followup"
)

data_list <- list(Pro_GMSI_all_tea,Pro_GMSI_all_epicatechin,Pro_GMSI_all_EGCG, Pro_GMSI_all_ECG) #Pro_GMSI_all_catechin, Pro_GMSI_all_epicatechin,Pro_GMSI_all_EGC,
Pro_GMSI_all <- data_list %>%
  reduce(dplyr::inner_join, by = c("ID", "Times"))
####验证指数计算####
# a <- data2_CatT4_all[,c("ID","followup",sub("_z","", Lmer_results2_all_level_sig_epicatechin[Lmer_results2_all_level_sig_epicatechin$Estimate > 0,]$Outcome),sub("_z","", Lmer_results2_all_level_sig_epicatechin[Lmer_results2_all_level_sig_epicatechin$Estimate < 0,]$Outcome))]
# log10(((476.6410 + 43.76440) / 2 + 1e-8) / (71.52110 + 1e-8))
####GMSI&Tea####
GMSI_Tea <- merge(Pro_GMSI_all, data2_CatT4_all_z[,c("ID","followup",Catechin_vars_all)], by.x = c("ID","Times"), by.y = c("ID","followup"))  %>%
  dplyr::filter(stats::complete.cases(.))


plot_trend_2group <- function(
    data,
    id_var,
    time_var,
    group_var,
    y_var,
    group_levels = c(0, 3),
    group_labels = c("Non", "More-than-daily"),
    time_levels = c("F0","F2","F3"),
    xlab = "Follow-up time",
    ylab = NULL,
    legend_title = NULL,
    base_size = 14,
    y_digits = 2  # 新增参数：统一 y 轴显示几位小数
) {
  
  library(dplyr)
  library(ggplot2)
  library(rlang)
  library(scales) # 必须加载 scales 包来处理数字格式
  
  id_sym    <- ensym(id_var)
  time_sym  <- ensym(time_var)
  group_sym <- ensym(group_var)
  
  if (is.character(y_var)) {
    y_chr <- y_var
  } else {
    y_chr <- as_string(ensym(y_var))
  }
  
  # 数据整理
  data_plot <- data %>%
    filter(!!group_sym %in% group_levels) %>%
    mutate(
      !!group_sym := factor(
        !!group_sym,
        levels = group_levels,
        labels = group_labels
      ),
      !!time_sym := factor(
        !!time_sym,
        levels = time_levels
      )
    )
  
  # y轴名称
  y_label_final <- ifelse(is.null(ylab), y_chr, ylab)
  
  # 画图
  p <- ggplot(data_plot,
              aes(x = !!time_sym,
                  y = .data[[y_chr]],
                  color = !!group_sym,
                  group = !!group_sym)) +
    
    stat_summary(fun = mean,
                 geom = "line",
                 linewidth = 1.2) +
    
    stat_summary(fun = mean,
                 geom = "point",
                 size = 3) +
    
    stat_summary(fun.data = mean_se,
                 geom = "errorbar",
                 width = 0.12) +
    
    # 核心修改：强制 y 轴刻度格式统一
    scale_y_continuous(labels = label_number(accuracy = 10^(-y_digits))) +
    
    labs(
      x = xlab,
      y = y_label_final,
      color = legend_title %||% as_string(group_sym)
    ) +
    
    theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold", size = 18),
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 14),
      legend.position = "bottom",
      axis.text = element_text(size = 16)
    )
  
  return(p)
}

p_trend1 <- plot_trend_2group(
  data = GMSI_Tea,
  id_var = ID,
  time_var = Times,
  group_var = Tea_freq,
  y_var = "PI_Tea",
  group_levels = c(0,3),
  group_labels = c("Non","More-than-daily"),
  ylab = "PI_Tea",
  legend_title = "Tea consumption frequency"
)

p_trend2 <- plot_trend_2group(
  data = GMSI_Tea,
  id_var = ID,
  time_var = Times,
  group_var = serum_I_epicatechin_F0_T,
  y_var = "PI_Epicatechin",
  group_levels = c(0,3),
  group_labels = c("Undetectable","High"),
  ylab = "PI_Epicatechin",
  legend_title = "Epicatechin"
)

# p_trend3 <- plot_trend_2group(
#   data = GMSI_Tea,
#   id_var = ID,
#   time_var = Times,
#   group_var = serum_I_EGC_F0_T,
#   y_var = "PI_EGC",
#   group_levels = c(0,3),
#   group_labels = c("Undetectable","High"),
#   ylab = "PI_EGC",
#   legend_title = "Epigallocatechin"
# )

p_trend4 <- plot_trend_2group(
  data = GMSI_Tea,
  id_var = ID,
  time_var = Times,
  group_var = serum_I_EGCG_F0_T,
  y_var = "PI_EGCG",
  group_levels = c(0,3),
  group_labels = c("Quartile 1","Quartile 4"),
  ylab = "PI_EGCG",
  legend_title = "Epigallocatechin gallate"
)

p_trend5 <- plot_trend_2group(
  data = GMSI_Tea,
  id_var = ID,
  time_var = Times,
  group_var = serum_I_ECG_F0_T,
  y_var = "PI_ECG",
  group_levels = c(0,3),
  group_labels = c("Quartile 1","Quartile 4"),
  ylab = "PI_ECG",
  legend_title = "Epicatechin gallate"
)

pdf("trend.pdf", width = 16, height = 10)
p_trend1 + p_trend2 + p_trend4 + p_trend5 + #p_trend3 + 
  plot_layout(
    nrow = 2
  )
dev.off()
####GMSI&ASMI-Lmer####
Pro_GMSI_all_ASMI <- merge(Pro_GMSI_all, ASM_Height_long[,c("ID","Followup","ASM","ASMI","WBTOT_FAT")], by.x = c("ID","Times"), by.y = c("ID","Followup"))
Pro_GMSI_all_ASMI <- merge(Pro_GMSI_all_ASMI, data2_CatT4_all_z, by.x = c("ID","Times"), by.y = c("ID","followup"))  %>%
  dplyr::filter(stats::complete.cases(.))
Pro_GMSI_all_ASMI <- scale_columns_group(Pro_GMSI_all_ASMI, c("ASM","ASMI"),"Sex_F0")
Pro_GMSI_all_ASMI <- scale_columns(Pro_GMSI_all_ASMI,c("PI_Tea", "PI_Epicatechin","PI_EGCG","PI_ECG"))

Lmer_results3_all <- process_lmer(c("PI_Tea_z", "PI_Epicatechin_z","PI_EGCG_z","PI_ECG_z"), #"PI_catechin", "PI_epicatechin", "PI_EGC",
                                     c("ASMI_z","ASM_z"), 
                                     Pro_GMSI_all_ASMI, 
                                     c(Covariates_all_lmer,"Tea_freq"), 
                                     2)
####相关系数可视化####
plot_corr_lower_triangle <- function(
    data,
    vars,
    method = "spearman",
    digits = 3,
    title = NULL,
    sig_level = c(0.05, 0.01, 0.001)
) {
  
  library(dplyr)
  library(ggplot2)
  library(Hmisc)
  
  ## 数据准备
  df_use <- data %>%
    dplyr::select(all_of(vars)) %>%
    na.omit()
  
  ## 相关计算
  res <- rcorr(as.matrix(df_use), type = method)
  
  ## 整理成长格式数据框
  df_plot <- expand.grid(
    Var1 = vars,
    Var2 = vars
  ) %>%
    mutate(
      rho = res$r[cbind(match(Var1, vars), match(Var2, vars))],
      p   = res$P[cbind(match(Var1, vars), match(Var2, vars))]
    ) %>%
    filter(match(Var1, vars) > match(Var2, vars)) %>%
    mutate(
      label = sprintf(paste0("%.", digits, "f"), rho),
      sig = case_when(
        p < sig_level[3] ~ "***",
        p < sig_level[2] ~ "**",
        p < sig_level[1] ~ "*",
        TRUE ~ ""
      )
    )
  
  ## 作图
  p <- ggplot(df_plot, aes(Var1, Var2)) +
    
    geom_tile(
      aes(fill = rho),
      color = "white",
      linewidth = 0.8
    ) +
    
    geom_text(
      aes(label = label),
      size = 4.8,
      color = "black"
    ) +
    
    geom_text(
      aes(label = sig),
      vjust = 1.8,
      size = 4,
      color = "black"
    ) +
    
    scale_fill_gradient2(
      low = "#E8EEF6",
      mid = "white",
      high = "#D6E3F3",
      midpoint = 0,
      limits = c(-1, 1),
      guide = "none"
    ) +
    
    coord_fixed() +
    
    labs(
      x = NULL,
      y = NULL,
      title = title
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        face = "bold",
        size = 12,
        color = "black"
      ),
      axis.text.y = element_text(
        face = "bold",
        size = 12,
        color = "black"
      ),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  ## 返回列表：图 + 结果数据框
  return(list(
    plot = p,
    results = df_plot
  ))
}

vars <- c(
  "PI_Tea", "PI_Epicatechin","PI_EGCG","PI_ECG", "ASM_z", "ASMI_z" #"PI_EGC",
)

plot1 <- plot_corr_lower_triangle(
  data   = Pro_GMSI_all_ASMI %>% filter(Times == "F2"),
  vars   = vars,
  method = "spearman",
  title  = "F2"
)

plot2 <- plot_corr_lower_triangle(
  data   = Pro_GMSI_all_ASMI %>% filter(Times == "F3"),
  vars   = vars,
  method = "spearman",
  title  = "F3"
)

pdf("Spearman.pdf", width = 12, height = 5)
plot1$plot + plot2$plot +
  plot_layout(
    nrow = 1
  )
dev.off()
#******************正文写作
a <- rbind(plot1$results,  plot2$results)
round(max(a$rho),3)
round(min(a$rho),3)
round(max(a$p),3)
round(min(a$p),3)

####蛋白质&ASM####
Pro_ASMI_all <- merge(ASM_Height_long[,c("ID","Followup","ASM","ASMI","WBTOT_FAT")], data2_CatT4_all_z, by.x = c("ID","Followup"), by.y = c("ID","followup")) %>%
  dplyr::filter(stats::complete.cases(.))
Pro_ASMI_all <- scale_columns_group(Pro_ASMI_all, c("ASM","ASMI"),"Sex_F0")

a <- unique(Pro_ASMI_all$ID)
Lmer_results7_all <- process_lmer(unique(Lmer_results2_all_level_sig$Outcome),
                                     c("ASMI_z","ASM_z"), 
                                     Pro_ASMI_all, 
                                     c(Covariates_all_lmer,"Tea_freq"), 
                                     0)
Lmer_results7_all$P_FDR <-  p.adjust(Lmer_results7_all$P_value, method = "fdr")
Lmer_results7_all_sig <- Lmer_results7_all[Lmer_results7_all$P_value < 0.05,]
Lmer_results7_all_sigFDR <- Lmer_results7_all[Lmer_results7_all$P_FDR < 0.05,]
head(Lmer_results7_all_sigFDR)
#***************可视化
Lmer_results7_all$Protein_name <- sub("_z", "", Lmer_results7_all$Predictor)
Lmer_results7_all$Protein_name <- toupper(Lmer_results7_all$Protein_name)
Lmer_results7_all <- merge(Lmer_results7_all, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Protein_name", by.y = "Protein_names")


df <- Lmer_results7_all %>%
  mutate(
    Significant_FDR = ifelse(P_FDR < 0.05, "Significant", "Not Significant"),
    hgnc_symbol = factor(hgnc_symbol, levels = unique(hgnc_symbol))
  )

plot <- ggplot(df, aes(x = hgnc_symbol, y = Estimate, color = Significant_FDR)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high),
                  position = position_dodge(width = 0.5),
                  size = 0.5, fatten = 2) +
  facet_wrap(~ Outcome, scales = "free_y", ncol = 1, strip.position = "top") +
  scale_color_manual(values = c("Significant" = "#E64B35", "Not Significant" = "#4DBBD5")) +
  theme_light(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", color = "#D0D0D0"),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  labs(
    x = "Predictor",
    y = "Estimate (95%CI)",
  )

pdf("plot.pdf", width = 10, height = 6)
print(plot)
dev.off()
####正文写作####
a1 <- df[df$P_FDR < 0.05,]
a12 <- a1[a1$Estimate >0,]
a12 <- a12[order(a12$hgnc_symbol),]
unique(a12$hgnc_symbol)
round(max(a12$P_FDR),3)
round(min(a12$P_FDR),3)
round(max(a12$Estimate),3)
round(min(a12$Estimate),3)

a11 <- a1[a1$Estimate <0,]
a11 <- a11[order(a11$hgnc_symbol),]
unique(a11$hgnc_symbol)
round(max(a11$P_FDR),3)
round(min(a11$P_FDR),3)
round(max(a11$Estimate),3)
round(min(a11$Estimate),3)
####中介效应分析####
exposures <- c("Tea_freq")
mediators <- c("PI_Tea")
outcomes  <- c("ASMI_z","ASM_z")

med_results1_pro <- run_mediation_glmer_batch(
  data        = Pro_GMSI_all_ASMI,
  id_var      = "ID",
  exposures   = exposures,
  mediators   = mediators,
  outcomes    = outcomes,
  covariates  = Covariates_all_lmer,
  exposure_high = 3,
  exposure_low  = 0,
  sims = 1000
)

exposures <- c("serum_I_epicatechin_F0_T")
mediators <- c("PI_Epicatechin")
outcomes  <- c("ASMI_z","ASM_z")

med_results2_pro <- run_mediation_glmer_batch(
  data        = Pro_GMSI_all_ASMI,
  id_var      = "ID",
  exposures   = exposures,
  mediators   = mediators,
  outcomes    = outcomes,
  covariates  = Covariates_all_lmer,
  exposure_high = 3,
  exposure_low  = 0,
  sims = 1000
)

# exposures <- c("serum_I_EGC_F0_T")
# mediators <- c("PI_EGC")
# outcomes  <- c("ASMI_z","ASM_z")
# 
# med_results3_pro <- run_mediation_glmer_batch(
#   data        = Pro_GMSI_all_ASMI,
#   id_var      = "ID",
#   exposures   = exposures,
#   mediators   = mediators,
#   outcomes    = outcomes,
#   covariates  = Covariates_all_lmer,
#   exposure_high = 3,
#   exposure_low  = 0,
#   sims = 1000
# )


# exposures <- c("serum_I_EGCG_F0_T")
# mediators <- c("PI_EGCG")
# outcomes  <- c("ASMI_z","ASM_z")
# 
# med_results4_pro <- run_mediation_glmer_batch(
#   data        = Pro_GMSI_all_ASMI,
#   id_var      = "ID",
#   exposures   = exposures,
#   mediators   = mediators,
#   outcomes    = outcomes,
#   covariates  = Covariates_all_lmer,
#   exposure_high = 3,
#   exposure_low  = 0,
#   sims = 1000
# )


exposures <- c("serum_I_ECG_F0_T")
mediators <- c("PI_ECG")
outcomes  <- c("ASMI_z","ASM_z")

med_results5_pro <- run_mediation_glmer_batch(
  data        = Pro_GMSI_all_ASMI,
  id_var      = "ID",
  exposures   = exposures,
  mediators   = mediators,
  outcomes    = outcomes,
  covariates  = Covariates_all_lmer,
  exposure_high = 3,
  exposure_low  = 0,
  sims = 1000
)

med_results_all <- rbind(med_results1_pro, med_results2_pro,med_results5_pro) #med_results4_pro,
####验证中介效应分析####
# a <- Pro_GMSI_all_ASMI[Pro_GMSI_all_ASMI$serum_I_ECG_F0_T == 3 | Pro_GMSI_all_ASMI$serum_I_ECG_F0_T == 0,]
# data <- Pro_GMSI_all_ASMI
# covariates <- Covariates_all_lmer
# cov_string <- paste(Covariates_all_lmer, collapse = " + ")
# safe_var <- function(v) paste0("`", v, "`")
# m <- "PI_ECG"
# x <- "serum_I_ECG_F0_T"
# y <- "ASMI_z"
# exposure_high <- 3
# exposure_low <- 0
# df_med <- data %>%
#   filter(.data[[x]] %in% c(exposure_low, exposure_high)) %>%
#   mutate(
#     Treat_bin = as.numeric(.data[[x]] == exposure_high)
#   )
# 
# id_var <- "ID"
# formula_m <- as.formula(
#   paste0(
#     safe_var(m), " ~ Treat_bin",
#     if(length(covariates) > 0) paste0(" + ", cov_string),
#     " + (1 | ", safe_var(id_var), ")"
#   )
# )
# model.m <- glmer(formula_m, data = df_med, family = gaussian())
# 
# formula_y <- as.formula(
#   paste0(
#     safe_var(y), " ~ Treat_bin + ", safe_var(m),
#     if(length(covariates) > 0) paste0(" + ", cov_string),
#     " + (1 | ", safe_var(id_var), ")"
#   )
# )
# model.y <- glmer(formula_y, data = df_med, family = gaussian())
# 
# set.seed(908)
# med_res <- mediation::mediate(model.m, model.y, treat = "Treat_bin", mediator = m, sims = 1000)
# summary(med_res)
####**********************其他**********************####
####年龄分布图和性别饼图####
# a <- data2_CatT4_all %>%
#   group_by(ID) %>%
#   slice(1) %>%
#   ungroup()
# table(a$Tea_freq)
# a1 <- unique(data1_CatT4_all$ID)


data_ga_all <- data1_CatT4_all %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup()
table(data_ga_all$Tea_freq)
data_ga_male <- data1_CatT4_male %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup()
data_ga <- rbind(data_ga_all, data_ga_male)
#********************年龄分布图
p_age <- ggplot(data_ga, aes(x = Age_F0)) +
  geom_density(
    fill = "#9ECAE1",
    color = "#3182BD",
    alpha = 0.7,
    linewidth = 1
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Age (years)",
    y = "Density",
    title = "Age distribution"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black")
  )

pdf("Age.pdf", height = 4, width = 5)
print(p_age)
dev.off()
#*********************性别分布图
data_ga$Sex_F0 <- recode(data_ga$Sex_F0,
                         "0" = "Female",
                         "1" = "Male")
df_sex <- data_ga %>%
  count(Sex_F0) %>%
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(Sex_F0, "\n", sprintf("%.1f%%", percent))
  )

p_sex <- ggplot(df_sex, aes(x = "", y = n, fill = Sex_F0)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4.5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#F4A6B8",
      "Male"   = "#9EC9E2"
    )
  ) +
  theme_void(base_size = 14) +
  labs(
    title = "Sex distribution",
    fill = "Sex"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )
pdf("Sex.pdf", height = 4, width = 5)
print(p_sex)
dev.off()
####随访间隔####
Follow_up_time <- Age_wide[Age_wide$ID %in% data_ga$ID,]
Follow_up_time <- merge(Follow_up_time, Cov_F0_all[,c("ID","Age_F0")], by = "ID")

Follow_up_time_filtered <- Follow_up_time %>%
  filter(!grepl("^NL4", ID))
Follow_up_time_filtered$F01 <- Follow_up_time_filtered$F1age-Follow_up_time_filtered$Age_F0
Follow_up_time$F12 <- Follow_up_time$F2age-Follow_up_time$F1age
Follow_up_time$F23 <- Follow_up_time$F3age-Follow_up_time$F2age
Follow_up_time$F34 <- Follow_up_time$F4age-Follow_up_time$F3age

round(mean(Follow_up_time_filtered$F01, na.rm = TRUE),1)
round(mean(Follow_up_time$F12, na.rm = TRUE),1)
round(mean(Follow_up_time$F23, na.rm = TRUE),1)
round(mean(Follow_up_time$F34, na.rm = TRUE),1)
####*****************************饮茶&菌群&骨骼肌-data3*****************************####
####儿茶素数据集处理合并####
Micro_selected2 <- Micro_selected_filter_clr[Micro_selected_filter_clr$Times %in% c("F2","F3"),]
Micro_ASMI <- merge(Micro_selected2, ASM_Height_long, by.x = c("ID","Times"), by.y = c("ID","Followup"))
data_list <- list(Cov_F0, Tea_Follow, Catechin, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], Micro_ASMI) #
data3_all_B <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
species_names <- grep("^s__", colnames(data3_all_B), value = TRUE)

data3_CatT4_all <- zero_tertile4(
  data = data3_all_B,
  vars = Catechin_vars3,
  suffix = "_T"
)

data3_CatT4_all <- data3_CatT4_all %>%
  mutate(
    across(
      all_of(Catechin_vars2),
      ~ quartile_group(.x),
      .names = "{.col}_T"
    )
  )

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(data3_CatT4_all)%in% factor_name)
for(i in idx ){
  data3_CatT4_all[[i]] <-  as.factor(data3_CatT4_all[[i]])
}

#***************ASMI根据性别Z分数，菌群Z分数
data3_all_B_final <- scale_columns_group(data3_CatT4_all, c("ASM", "ASMI"), c("Sex_F0"))
data3_all_B_final[, species_names] <- scale(data3_all_B_final[, species_names], center = TRUE, scale = TRUE)
a <- unique(data3_all_B_final$ID)
####FC计算数据集####
data3_all_B_final_FC <- merge(data3_all_B_final[,c("ID","Times", Catechin_vars_all)], Micro_selected_filter$filtered_data, by = c("ID","Times"))

param_all <- tibble::tribble(
  ~Times, ~group_var, ~y_var,  ~xlab,  ~ylab, ~labels,
  "F2",   "Tea_freq", "s__Erysipelatoclostridium_ramosum", "Tea consumption frequency",  "Erysipelatoclostridium ramosum", c("Never","Low","Moderate","High"),
  "F3",   "Tea_freq", "s__Erysipelatoclostridium_ramosum", "Tea consumption frequency",  "Erysipelatoclostridium ramosum", c("Never","Low","Moderate","High"),
  
  "F2",   "serum_I_ECG_F0_T", "s__Erysipelatoclostridium_ramosum", "Epicatechin gallate", "Erysipelatoclostridium ramosum",  c("Zero","Low","Moderate","High"),
  "F3",   "serum_I_ECG_F0_T", "s__Erysipelatoclostridium_ramosum", "Epicatechin gallate", "Erysipelatoclostridium ramosum", c("Zero","Low","Moderate","High"),
)

plots_all <- param_all %>%
  mutate(
    plot = pmap(
      list(Times, group_var, y_var, xlab, ylab, labels),
      ~ plot_wilcox_violin3(
        data = data3_all_B_final[data3_all_B_final$Times == ..1, ],
        group_var = ..2,          # ← 关键修改点
        y_var     = ..3,
        group_levels = c(0, 1, 2, 3),
        group_labels = ..6,
        xlab = ..4,
        ylab = ..5,
        title = ..1
      )
    )
  )
png("Erysipelatoclostridium ramosum.png", width = 12, height = 8, units = "in", res = 500)
wrap_plots(plots_all$plot, ncol = 3)
dev.off()
####FC计算####
exposure_vars1 <- c(
  Tea = "Tea_freq"
)

exposure_vars2 <- c(
  catechin     = "serum_I_catechin_F0_T",
  epicatechin  = "serum_I_epicatechin_F0_T",
  EGC          = "serum_I_EGC_F0_T",
  EGCG         = "serum_I_EGCG_F0_T",
  ECG          = "serum_I_ECG_F0_T"
)

FC_results1_F2_tea <- lapply(exposure_vars1, function(x) {
  run_FC_pairwise_analysis(
    data3_all_B_final_FC[data3_all_B_final_FC$Times == "F2",],
    x,
    col = species_names,
    p_cutoff = 0.05,
    fc_cutoff = 1,
    use_fdr = FALSE
  )
})

FC_results1_F3_tea <- lapply(exposure_vars1, function(x) {
  run_FC_pairwise_analysis(
    data3_all_B_final_FC[data3_all_B_final_FC$Times == "F3",],
    x,
    col = species_names,
    p_cutoff = 0.05,
    fc_cutoff = 1,
    use_fdr = FALSE
  )
})

FC_results1_F2_flo <- lapply(exposure_vars2, function(x) {
  run_FC_pairwise_analysis(
    data3_all_B_final_FC[data3_all_B_final_FC$Times == "F2",],
    x,
    col = species_names,
    p_cutoff = 0.05,
    fc_cutoff = 1,
    use_fdr = FALSE
  )
})

FC_results1_F3_flo <- lapply(exposure_vars2, function(x) {
  run_FC_pairwise_analysis(
    data3_all_B_final_FC[data3_all_B_final_FC$Times == "F3",],
    x,
    col = species_names,
    p_cutoff = 0.05,
    fc_cutoff = 1,
    use_fdr = FALSE
  )
})


Tea_intersect <- c(FC_results1_F2_tea$Tea$sig_results[FC_results1_F2_tea$Tea$sig_results$Comparison == "3_vs_0",]$Species,
                   FC_results1_F3_tea$Tea$sig_results[FC_results1_F3_tea$Tea$sig_results$Comparison == "3_vs_0",]$Species)

catechin_intersect <- c(FC_results1_F2_flo$catechin$sig_results[FC_results1_F2_flo$catechin$sig_results$Comparison == "3_vs_0",]$Species,
                        FC_results1_F3_flo$catechin$sig_results[FC_results1_F3_flo$catechin$sig_results$Comparison == "3_vs_0",]$Species)

epicatechin_intersect <- c(FC_results1_F2_flo$epicatechin$sig_results[FC_results1_F2_flo$epicatechin$sig_results$Comparison == "3_vs_0",]$Species,
                           FC_results1_F3_flo$epicatechin$sig_results[FC_results1_F3_flo$epicatechin$sig_results$Comparison == "3_vs_0",]$Species)

EGC_intersect <- c(FC_results1_F2_flo$EGC$sig_results[FC_results1_F2_flo$EGC$sig_results$Comparison == "3_vs_0",]$Species,
                   FC_results1_F3_flo$EGC$sig_results[FC_results1_F3_flo$EGC$sig_results$Comparison == "3_vs_0",]$Species)

EGCG_intersect <- c(FC_results1_F2_flo$EGCG$sig_results[FC_results1_F2_flo$EGCG$sig_results$Comparison == "3_vs_0",]$Species,
                    FC_results1_F3_flo$EGCG$sig_results[FC_results1_F3_flo$EGCG$sig_results$Comparison == "3_vs_0",]$Species)

ECG_intersect <- c(FC_results1_F2_flo$ECG$sig_results[FC_results1_F2_flo$ECG$sig_results$Comparison == "3_vs_0",]$Species,
                   FC_results1_F3_flo$ECG$sig_results[FC_results1_F3_flo$ECG$sig_results$Comparison == "3_vs_0",]$Species)

####FC可视化####
plot_list <- list(
  list(data = FC_results1_F2_tea$Tea$all_results, 
       title = "Tea frequency (F2)", 
       comparison_labels = c("3_vs_0" = "More-than-daily vs. Non")),
  list(data = FC_results1_F3_tea$Tea$all_results, 
       title = "Tea frequency (F3)", 
       comparison_labels = c("3_vs_0" = "More-than-daily vs. Non")),
  list(data = FC_results1_F2_flo$catechin$all_results, 
       title = "Catechin (F2)", 
       comparison_labels = c("3_vs_0" = "High vs. Undetectable")),
  list(data = FC_results1_F3_flo$catechin$all_results, 
       title = "Catechin (F3)", 
       comparison_labels = c("3_vs_0" = "High vs. Undetectable")),
  list(data = FC_results1_F2_flo$epicatechin$all_results, 
       title = "Epicatechin (F2)", 
       comparison_labels = c("3_vs_0" = "High vs. Undetectable")),
  list(data = FC_results1_F3_flo$epicatechin$all_results, 
       title = "Epicatechin (F3)", 
       comparison_labels = c("3_vs_0" = "High vs. Undetectable")),
  list(data = FC_results1_F2_flo$EGC$all_results, 
       title = "Epigallocatechin (F2)", 
       comparison_labels = c("3_vs_0" = "High vs. Undetectable")),
  list(data = FC_results1_F3_flo$EGC$all_results, 
       title = "Epigallocatechin (F3)", 
       comparison_labels = c("3_vs_0" = "High vs. Undetectable")),
  list(data = FC_results1_F2_flo$EGCG$all_results, 
       title = "Epigallocatechin gallate (F2)", 
       comparison_labels = c("3_vs_0" = "Quartile 1 vs. Quartile 4")),
  list(data = FC_results1_F3_flo$EGCG$all_results, 
       title = "Epigallocatechin gallate (F3)", 
       comparison_labels = c("3_vs_0" = "Quartile 1 vs. Quartile 4")),
  list(data = FC_results1_F2_flo$ECG$all_results, 
       title = "Epicatechin gallate (F2)", 
       comparison_labels = c("3_vs_0" = "Quartile 1 vs. Quartile 4")),
  list(data = FC_results1_F3_flo$ECG$all_results, 
       title = "Epicatechin gallate (F3)", 
       comparison_labels = c("3_vs_0" = "Quartile 1 vs. Quartile 4"))
)

# 循环绘图
plot_objs <- list()

for (i in seq_along(plot_list)) {
  p <- plot_list[[i]]
  plot_objs[[i]] <- plot_volcano_species(
    p$data[p$data$Comparison == "3_vs_0", ],
    title = p$title,
    label = TRUE,
    comparison_labels = p$comparison_labels,
    fc_cutoff = 1,
    p_cutoff = 0.05
  )
}

combined_plot <- wrap_plots(plot_objs, ncol = 4, guides = "collect") & theme(legend.position = "bottom")
pdf("volcano.pdf", width = 20, height = 20)
print(combined_plot)
dev.off()


# combined_plot <- wrap_plots(plot_objs[1:4], ncol = 2, guides = "collect") & theme(legend.position = "bottom")
# png("volcano1.png", width = 16, height = 16, units = "in", res = 500)
# print(combined_plot)
# dev.off()
# 
# combined_plot <- wrap_plots(plot_objs[5:8], ncol = 2, guides = "collect") & theme(legend.position = "bottom")
# png("volcano2.png", width = 14, height = 14, units = "in", res = 500)
# print(combined_plot)
# dev.off()
# 
# combined_plot <- wrap_plots(plot_objs[9:12], ncol = 2, guides = "collect") & theme(legend.position = "bottom")
# png("volcano3.png", width = 14, height = 14, units = "in", res = 500)
# print(combined_plot)
# dev.off()
####Lmer####
run_lmer_level3_fdr <- function(
    exposure_vars,
    intersect_var,
    data = data3_CatT4_all,
    covariates = Covariates_all_lmer,
    gender = 2,
    level_keep = 3,
    fdr_cutoff = 0.05
) {
  
  ## 1. 运行 lmer
  lmer_res <- process_lmer(
    exposure_vars,
    intersect_var,
    data,
    covariates,
    gender
  )
  
  ## 2. 只保留指定 Level
  lmer_level <- lmer_res %>%
    dplyr::filter(Level == level_keep)
  
  ## 3. FDR 校正
  lmer_level <- lmer_level %>%
    dplyr::mutate(
      P_FDR = p.adjust(P_value, method = "fdr")
    )
  
  ## 4a. nominal 显著（P < 0.05）
  res_p05 <- lmer_level %>%
    dplyr::filter(P_value < 0.05)
  
  ## 4b. FDR 显著
  res_fdr <- lmer_level %>%
    dplyr::filter(P_FDR < fdr_cutoff)
  
  ## 5. 返回两个结果
  return(
    list(
      all = lmer_level,
      P_value_0.05 = res_p05,
      FDR_0.05 = res_fdr
    )
  )
}


Tea_data3_sigFDR <- run_lmer_level3_fdr(
  exposure_vars = c("Tea_freq"),
  intersect_var = unique(Tea_intersect),
  data = data3_all_B_final
)

catechin_data3_sigFDR <- run_lmer_level3_fdr(
  exposure_vars = c("serum_I_catechin_F0_T"),
  intersect_var = unique(catechin_intersect),
  data = data3_all_B_final
)

epicatechin_data3_sigFDR <- run_lmer_level3_fdr(
  exposure_vars = c("serum_I_epicatechin_F0_T"),
  intersect_var = unique(epicatechin_intersect),
  data = data3_all_B_final
)

EGC_data3_sigFDR <- run_lmer_level3_fdr(
  exposure_vars = c("serum_I_EGC_F0_T"),
  intersect_var = unique(EGC_intersect),
  data = data3_all_B_final
)

EGCG_data3_sigFDR <- run_lmer_level3_fdr(
  exposure_vars = c("serum_I_EGCG_F0_T"),
  intersect_var = unique(EGCG_intersect),
  data = data3_all_B_final
)

ECG_data3_sigFDR <- run_lmer_level3_fdr(
  exposure_vars = c("serum_I_ECG_F0_T"),
  intersect_var = unique(ECG_intersect),
  data = data3_all_B_final
)

Lmer_results3_sigFDR <- rbind(Tea_data3_sigFDR$FDR_0.05, catechin_data3_sigFDR$FDR_0.05, epicatechin_data3_sigFDR$FDR_0.05, EGC_data3_sigFDR$FDR_0.05, EGCG_data3_sigFDR$FDR_0.05, ECG_data3_sigFDR$FDR_0.05)
Lmer_results3_sig <- rbind(Tea_data3_sigFDR$P_value_0.05, catechin_data3_sigFDR$P_value_0.05, epicatechin_data3_sigFDR$P_value_0.05, EGC_data3_sigFDR$P_value_0.05, EGCG_data3_sigFDR$P_value_0.05, ECG_data3_sigFDR$P_value_0.05)
Lmer_results3_all_plot <- rbind(Tea_data3_sigFDR$all, catechin_data3_sigFDR$all, epicatechin_data3_sigFDR$all, EGC_data3_sigFDR$all, EGCG_data3_sigFDR$all, ECG_data3_sigFDR$all)

####正文写作####
a1 <- Tea_data3_sigFDR$FDR_0.05
a12 <- a1[a1$Estimate >0, ]
a12$Outcome
round(max(a12$P_FDR),3)
round(min(a12$P_FDR),3)
round(max(a12$Estimate),3)
round(min(a12$Estimate),3)

a11 <- a1[a1$Estimate <0, ]
a11$Outcome
round(max(a11$P_FDR),3)
round(min(a11$P_FDR),3)
round(max(a11$Estimate),3)
round(min(a11$Estimate),3)

a2 <- Lmer_results3_sigFDR
a21 <- a2[a2$Outcome== "s__Faecalimonas_umbilicata",]
round(max(a21$Estimate),3)
round(min(a21$Estimate),3)
round(max(a21$P_FDR),3)
round(min(a21$P_FDR),3)

a22 <- a2[a2$Outcome== "s__Erysipelatoclostridium_ramosum",]
round(max(a22$Estimate),3)
round(min(a22$Estimate),3)
round(max(a22$P_FDR),3)
round(min(a22$P_FDR),3)

a23 <- a2[a2$Outcome== "s__Adlercreutzia_equolifaciens",]
round(max(a23$Estimate),3)
round(min(a23$Estimate),3)
round(max(a23$P_FDR),3)
round(min(a23$P_FDR),3)
####Lmer可视化####
# 1. 数据预处理
df_plot <- Lmer_results3_all_plot %>%
  mutate(
    # 区分正负向
    Direction = ifelse(Estimate > 0, "Positive", "Negative"),
    # 关键：根据 P_FDR 区分显著性（通常阈值为 0.05）
    Is_Sig = ifelse(P_FDR < 0.05, "Significant", "Non-significant"),
    # 清理菌群名称
    Species = str_replace(Outcome, "^s__", ""),
    Species = str_replace_all(Species, "_", " ")
  )

# 2. 统一 Y 轴排序
species_order <- df_plot %>%
  group_by(Species) %>%
  summarise(order_val = median(Estimate, na.rm = TRUE)) %>%
  arrange(order_val) %>%
  pull(Species)

df_plot$Species <- factor(df_plot$Species, levels = species_order)

# 3. 重新编码 Predictor
df_plot$Predictor <- recode(
  df_plot$Predictor,
  "Tea_freq" = "Tea consumption",
  "serum_I_catechin_F0_T"    = "Catechin",
  "serum_I_epicatechin_F0_T" = "Epicatechin",
  "serum_I_EGC_F0_T"         = "Epigallocatechin",
  "serum_I_ECG_F0_T"         = "Epicatechin gallate",
  "serum_I_EGCG_F0_T"        = "Epigallocatechin gallate"
)

df_plot$Predictor <- factor(
  df_plot$Predictor,
  levels = c("Tea consumption", "Catechin", "Epicatechin", "Epigallocatechin", "Epigallocatechin gallate", "Epicatechin gallate")
)

# 4. 绘图
p_forest <- ggplot(
  df_plot,
  aes(
    x = Estimate,
    y = Species,
    xmin = CI_low,
    xmax = CI_high,
    color = Direction,
    shape = Is_Sig  # 映射形状到显著性
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.6) +
  
  # 误差线：如果不显著，可以让线细一点或颜色淡一点
  geom_errorbarh(
    aes(alpha = Is_Sig), # 显著性映射到透明度，让非显著项不抢眼
    height = 0.2, 
    linewidth = 0.8
  ) +
  
  # 数据点：实心代表显著，空心代表不显著
  geom_point(
    aes(fill = Direction, alpha = Is_Sig), 
    size = 2.8, 
    stroke = 1  # 边框粗细
  ) +
  
  facet_wrap(~ Predictor, nrow = 1) +
  
  # 形状设置：19 是实心圆，21 是可以填充颜色的空心圆
  scale_shape_manual(values = c("Significant" = 19, "Non-significant" = 21)) +
  
  # 透明度设置：非显著项设为 0.4，显著项 1
  scale_alpha_manual(values = c("Significant" = 1, "Non-significant" = 0.4), guide = "none") +
  
  scale_color_manual(
    values = c("Positive" = "#C04A3A", "Negative" = "#3B7EA1")
  ) +
  
  # 保证空心圆的填充色逻辑
  scale_fill_manual(
    values = c("Positive" = "#C04A3A", "Negative" = "#3B7EA1"),
    guide = "none"
  ) +
  
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  
  labs(
    x = "Estimate (95% CI)",
    y = NULL,
    color = "Effect Direction",
    shape = "FDR Significance"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "italic", size = 12.5),
    axis.title.x = element_text(face = "bold", margin = ggplot2::margin(t = 15)),
    strip.text = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "gray95"),
    panel.spacing.x = unit(1.5, "lines"),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(t = 15, r = 15, b = 0, l = 15, unit = "pt")
    #legend.box = "vertical" # 标签多时垂直排列更整齐
  )

# 5. 保存
pdf("Tea_Microbiome.pdf", height = 20, width = 19)
print(p_forest)
dev.off()
####计算GMSI####
unique(Lmer_results3_sigFDR$Predictor)
Micro_GMSI_all_tea <- calculate_GMSI(
  data = data3_all_B_final_FC,
  outcomes_pos = Tea_data3_sigFDR$FDR_0.05[Tea_data3_sigFDR$FDR_0.05$Estimate > 0,]$Outcome,
  outcomes_neg = Tea_data3_sigFDR$FDR_0.05[Tea_data3_sigFDR$FDR_0.05$Estimate < 0,]$Outcome,
  "GMI_Tea",
  time_var = "Times"
)

# Micro_GMSI_all_catechin <- calculate_GMSI(
#   data = data3_all_B_final_FC,
#   outcomes_pos = catechin_data3_sigFDR$FDR_0.05[catechin_data3_sigFDR$FDR_0.05$Estimate > 0,]$Outcome,
#   outcomes_neg = catechin_data3_sigFDR$FDR_0.05[catechin_data3_sigFDR$FDR_0.05$Estimate < 0,]$Outcome,
#   "GMI_catechin",
#   time_var = "Times"
# )

Micro_GMSI_all_epicatechin <- calculate_GMSI(
  data = data3_all_B_final_FC,
  outcomes_pos = epicatechin_data3_sigFDR$FDR_0.05[epicatechin_data3_sigFDR$FDR_0.05$Estimate > 0,]$Outcome,
  outcomes_neg = epicatechin_data3_sigFDR$FDR_0.05[epicatechin_data3_sigFDR$FDR_0.05$Estimate < 0,]$Outcome,
  "GMI_Epicatechin",
  time_var = "Times"
)

Micro_GMSI_all_ECG <- calculate_GMSI(
  data = data3_all_B_final_FC,
  outcomes_pos = ECG_data3_sigFDR$FDR_0.05[ECG_data3_sigFDR$FDR_0.05$Estimate > 0,]$Outcome,
  outcomes_neg = ECG_data3_sigFDR$FDR_0.05[ECG_data3_sigFDR$FDR_0.05$Estimate < 0,]$Outcome,
  "GMI_ECG",
  time_var = "Times"
)

Micro_GMSI_all_EGC <- calculate_GMSI(
  data = data3_all_B_final_FC,
  outcomes_pos = EGC_data3_sigFDR$FDR_0.05[EGC_data3_sigFDR$FDR_0.05$Estimate > 0,]$Outcome,
  outcomes_neg = EGC_data3_sigFDR$FDR_0.05[EGC_data3_sigFDR$FDR_0.05$Estimate < 0,]$Outcome,
  "GMI_EGC",
  time_var = "Times"
)

Micro_GMSI_all_EGCG <- calculate_GMSI(
  data = data3_all_B_final_FC,
  outcomes_pos = EGCG_data3_sigFDR$FDR_0.05[EGCG_data3_sigFDR$FDR_0.05$Estimate > 0,]$Outcome,
  outcomes_neg = EGCG_data3_sigFDR$FDR_0.05[EGCG_data3_sigFDR$FDR_0.05$Estimate < 0,]$Outcome,
  "GMI_EGCG",
  time_var = "Times"
)

#data_list <- list(Micro_GMSI_all_tea, Micro_GMSI_all_catechin, Micro_GMSI_all_epicatechin, Micro_GMSI_all_ECG, Micro_GMSI_all_EGC, Micro_GMSI_all_EGCG)
data_list <- list(Micro_GMSI_all_tea, Micro_GMSI_all_epicatechin, Micro_GMSI_all_ECG, Micro_GMSI_all_EGC, Micro_GMSI_all_EGCG)
Micro_GMSI_all <- data_list %>%
  reduce(dplyr::inner_join, by = c("ID", "Times"))

Micro_GMSI_all_ASMI <- merge(Micro_GMSI_all, ASM_Height_long[,c("ID","Followup","ASM","ASMI","WBTOT_FAT")], by.x = c("ID","Times"), by.y = c("ID","Followup"))
Micro_GMSI_all_ASMI <- merge(Micro_GMSI_all_ASMI, dplyr::select(data3_CatT4_all,-c("ASM","ASMI","WBTOT_FAT")), by = c("ID","Times")) %>%
  dplyr::filter(stats::complete.cases(.))
table(Micro_GMSI_all_ASMI$Times)

Micro_GMSI_all_ASMI <- scale_columns_group(Micro_GMSI_all_ASMI, c("ASM", "ASMI"), c("Sex_F0"))
Micro_GMSI_all_ASMI <- scale_columns(Micro_GMSI_all_ASMI,c("GMI_Tea","GMI_Epicatechin", "GMI_EGC", "GMI_ECG","GMI_EGCG"))

####GMSI&Tea####
plot_wilcox_violin3 <- function(
    data,
    group_var,
    y_var,
    group_levels = NULL,
    group_labels = NULL,
    compare_levels = NULL,   # ⭐ 新增参数
    ylab = NULL,
    y_unit = NULL,
    xlab = NULL,
    title = NULL,
    base_size = 14,
    violin_fill = "#E6E6E6"
) {
  
  library(ggplot2)
  library(dplyr)
  library(rlang)
  
  group_sym <- ensym(group_var)
  group_chr <- as_string(group_sym)
  
  if (is.character(y_var)) {
    y_chr <- y_var
  } else {
    y_chr <- as_string(ensym(y_var))
  }
  
  ## 如果指定 compare_levels → 先筛选
  if (!is.null(compare_levels)) {
    data <- data %>%
      filter(!!group_sym %in% compare_levels)
  }
  
  ## 分组 factor
  if (!is.null(group_levels)) {
    
    ## 如果只选部分组，要同步截取 levels & labels
    if (!is.null(compare_levels)) {
      keep_index <- which(group_levels %in% compare_levels)
      group_levels <- group_levels[keep_index]
      if (!is.null(group_labels)) {
        group_labels <- group_labels[keep_index]
      }
    }
    
    data <- data %>%
      mutate(
        !!group_sym := factor(
          !!group_sym,
          levels = group_levels,
          labels = group_labels %||% group_levels
        )
      )
  } else {
    data <- data %>%
      mutate(!!group_sym := factor(!!group_sym))
  }
  
  ## 组数
  n_group <- nlevels(data[[group_chr]])
  
  ## 统计检验
  p_val <- if (n_group == 2) {
    wilcox.test(
      data[[y_chr]] ~ data[[group_chr]],
      exact = FALSE
    )$p.value
  } else {
    kruskal.test(
      data[[y_chr]] ~ data[[group_chr]]
    )$p.value
  }
  
  ## P 值文本
  p_text <- if (p_val < 0.001) {
    "italic(P)<0.001"
  } else {
    paste0("italic(P)==\"", sprintf("%.3f", p_val), "\"")
  }
  
  ## y 轴范围
  ymax   <- max(data[[y_chr]], na.rm = TRUE)
  yrange <- diff(range(data[[y_chr]], na.rm = TRUE))
  
  x_p <- mean(seq_len(n_group))
  y_p <- ymax + 0.15 * yrange
  
  ## y 轴标题
  y_label_final <- if (!is.null(ylab)) {
    if (!is.null(y_unit)) {
      paste0(ylab, " (", y_unit, ")")
    } else {
      ylab
    }
  } else {
    y_chr
  }
  
  ggplot(data, aes(x = !!group_sym, y = .data[[y_chr]])) +
    
    geom_violin(
      fill = violin_fill,
      color = NA,
      alpha = 0.6,
      width = 0.85
    ) +
    
    geom_boxplot(
      width = 0.28,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.6,
      median.linewidth = 0.9
    ) +
    
    stat_summary(
      fun = median,
      geom = "point",
      size = 1.6,
      shape = 16,
      color = "black"
    ) +
    
    annotate(
      "rect",
      xmin = x_p - 0.45,
      xmax = x_p + 0.45,
      ymin = y_p - 0.05 * yrange,
      ymax = y_p + 0.05 * yrange,
      fill = "grey95",
      color = "black"
    ) +
    
    annotate(
      "text",
      x = x_p,
      y = y_p,
      label = p_text,
      parse = TRUE,
      size = 4.2,
      color = "black"
    ) +
    
    labs(
      title = title,
      x = xlab %||% group_chr,
      y = y_label_final
    ) +
    
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      axis.line  = element_line(linewidth = 0.6),
      plot.margin = ggplot2::margin(
        t = 10, r = 10, b = 10, l = 10, unit = "pt"
      )
    )
}


param_all <- tibble::tribble(
  ~Times, ~group_var, ~y_var,  ~xlab,  ~ylab, ~labels,
  "F2",   "Tea_freq", "GMI_Tea", "Tea consumption frequency",  "GMI_Tea", c("Non","Less-than-daily","Daily","More-than-daily"),
  "F3",   "Tea_freq", "GMI_Tea", "Tea consumption frequency",  "GMI_Tea", c("Non","Less-than-daily","Daily","More-than-daily"),
  
  "F2",   "serum_I_epicatechin_F0_T", "GMI_Epicatechin", "Epicatechin",  "GMI_Epicatechin", c("Undetectable","Low","Medium","High"),
  "F3",   "serum_I_epicatechin_F0_T", "GMI_Epicatechin", "Epicatechin",  "GMI_Epicatechin", c("Undetectable","Low","Medium","High"),
  
  "F2",   "serum_I_EGC_F0_T", "GMI_EGC", "Epigallocatechin",  "GMI_EGC", c("Undetectable","Low","Medium","High"),
  "F3",   "serum_I_EGC_F0_T", "GMI_EGC", "Epigallocatechin",  "GMI_EGC", c("Undetectable","Low","Medium","High"),
  
  "F2",   "serum_I_EGCG_F0_T", "GMI_EGCG", "Epigallocatechin gallate", "GMI_EGCG",  c("Quartile 1","Quartile 2","Quartile 3","Quartile 4"),
  "F3",   "serum_I_EGCG_F0_T", "GMI_EGCG", "Epigallocatechin gallate", "GMI_EGCG", c("Quartile 1","Quartile 2","Quartile 3","Quartile 4"),
  
  "F2",   "serum_I_ECG_F0_T", "GMI_ECG", "Epicatechin gallate", "GMI_ECG",  c("Quartile 1","Quartile 2","Quartile 3","Quartile 4"),
  "F3",   "serum_I_ECG_F0_T", "GMI_ECG", "Epicatechin gallate", "GMI_ECG", c("Quartile 1","Quartile 2","Quartile 3","Quartile 4")
)

plots_all <- param_all %>%
  mutate(
    plot = pmap(
      list(Times, group_var, y_var, xlab, ylab, labels),
      ~ plot_wilcox_violin3(
        data = Micro_GMSI_all_ASMI[Micro_GMSI_all_ASMI$Times == ..1, ],
        group_var = ..2,
        y_var     = ..3,
        group_levels = c(0, 1, 2, 3),
        group_labels = ..6,
        compare_levels = c(0, 3),   # ⭐ 关键：只比较0和3
        xlab = ..4,
        ylab = ..5,
        title = ..1
      )
    )
  )

pdf("GMI_tea.pdf", width = 14, height = 10)
wrap_plots(plots_all$plot, ncol = 4)
dev.off()
####指数之间的相关性####
vars <- c(
  "GMI_Tea", "GMI_Epicatechin","GMI_EGC", "GMI_EGCG",
  "GMI_ECG", "ASM_z", "ASMI_z"
)

plot1 <- plot_corr_lower_triangle(
  data   = Micro_GMSI_all_ASMI %>% filter(Times == "F2"),
  vars   = vars,
  method = "spearman",
  title  = "F2"
)

plot2 <- plot_corr_lower_triangle(
  data   = Micro_GMSI_all_ASMI %>% filter(Times == "F3"),
  vars   = vars,
  method = "spearman",
  title  = "F3"
)

pdf("Spearman2.pdf", height = 5.5, width = 11)
plot1$plot + plot2$plot +
  plot_layout(
    nrow = 1
  )
dev.off()

#******************正文写作
a <- rbind(plot1$results,  plot2$results)
a1 <- a[a$p < 0.05,]
round(max(a1$rho),3)
round(min(a1$rho),3)
round(max(a1$p),3)
round(min(a1$p),3)
####GMSI&ASMI-Lmer####
Lmer_results5_all <- process_lmer(c("GMI_Tea_z","GMI_Epicatechin_z", "GMI_EGC_z", "GMI_ECG_z","GMI_EGCG_z"),  
                                  c("ASMI_z","ASM_z"), 
                                  Micro_GMSI_all_ASMI, 
                                  c(Covariates_all_lmer,"Tea_freq"), 
                                  2)
####GMI/PI&ASMI一起可视化####
predictor_order <- c(
  "PI_Tea", "PI_Epicatechin", "PI_EGCG", "PI_ECG",
  "GMI_Tea", "GMI_Epicatechin", "GMI_EGC", "GMI_EGCG", "GMI_ECG"
)
outcome_order <- c("ASMI_z", "ASM_z")
data_scaled <- rbind(Lmer_results3_all, Lmer_results5_all) %>%
  mutate(
    Predictor_type = case_when(
      grepl("^PI_", Predictor)  ~ "PI",
      grepl("^GMI_", Predictor) ~ "GMI"
    ),
   Predictor = sub("_z","",Predictor),
   Predictor = factor(Predictor, levels = predictor_order),
   Outcome   = factor(Outcome, levels = outcome_order)
  )%>%
  arrange(Predictor, Outcome)

# 计算范围
gmi_range <- range(data_scaled$Estimate[data_scaled$Predictor_type == "GMI"])
pi_range  <- range(data_scaled$Estimate[data_scaled$Predictor_type == "PI"])

g_min <- gmi_range[1]
g_max <- gmi_range[2]
p_min <- pi_range[1]
p_max <- pi_range[2]

# 对 PI 做线性缩放
data_scaled <- data_scaled %>%
  mutate(
    Estimate_plot = ifelse(
      Predictor_type == "PI",
      (Estimate - p_min) / (p_max - p_min) * (g_max - g_min) + g_min,
      Estimate
    ),
    CI_low_plot = ifelse(
      Predictor_type == "PI",
      (CI_low - p_min) / (p_max - p_min) * (g_max - g_min) + g_min,
      CI_low
    ),
    CI_high_plot = ifelse(
      Predictor_type == "PI",
      (CI_high - p_min) / (p_max - p_min) * (g_max - g_min) + g_min,
      CI_high
    )
  )


plot_lmer_forest_dualx_optimized <- function(
    data, 
    facet_var = "Outcome", 
    base_size = 14
) {
  
  data <- data %>%
    mutate(
      Predictor = factor(Predictor, levels = rev(unique(Predictor))),
      Sig_flag = ifelse(Significance == "Significant", "Yes", "No")
    )
  
  ggplot(data, aes(x = Estimate_plot, y = Predictor)) +
    
    # 2. 误差线：稍微加宽，增加专业感
    geom_errorbarh(
      aes(xmin = CI_low_plot, xmax = CI_high_plot),
      height = 0.2,
      linewidth = 0.7,
      color = "#2C3E50" # 深灰蓝，比纯黑更有质感
    ) +
    
    # 3. 散点：使用更现代的配色
    geom_point(
      aes(fill = Sig_flag),
      size = 4,
      shape = 21,
      color = "white", # 白色描边，产生立体感
      stroke = 0.5
    ) +
    
    facet_wrap(
      as.formula(paste("~", facet_var))#,
      #scales = "free_x",
      
    ) +
    
    # 4. 配色优化
    scale_fill_manual(
      values = c("Yes" = "#E74C3C", "No" = "#BDC3C7"), # 显著为红，不显著为浅灰
      name = "Significant"
    ) +
    
    scale_x_continuous(
      name = "Estimate (95%CI)-GMI",
      expand = expansion(mult = 0.1), # 给两侧留点呼吸空间
      sec.axis = sec_axis(
        ~ (. - g_min) / (g_max - g_min) * (p_max - p_min) + p_min,
        name = "Estimate (95%CI)-PI"
      )
    ) +
    
    labs(y = NULL) +
    # 5. 精细化主题定制
    theme_minimal(base_size = base_size) + # 换成 minimal 基础，去掉多余边框
    theme(
      panel.grid.major = element_blank(), # 移除横向主网格线
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80", linewidth = 0.5), # 加上细边框包围分面
      
      strip.background = element_rect(fill = "gray95", color = NA), # 分面标题背景
      strip.text = element_text(face = "bold", color = "#2C3E50"),
      
      axis.title.x.bottom = element_text(face = "bold",margin = ggplot2::margin(t = 15)),
      axis.title.x.top    = element_text(face = "bold",margin = ggplot2::margin(b = 15)),
      
      axis.text.y = element_text(face = "bold", color = "#2C3E50"),
      panel.spacing = unit(2, "lines"),
      legend.position = "right",
      legend.background = element_blank()
    )
}

lmer_plot <- plot_lmer_forest_dualx_optimized(data_scaled)
pdf("plot2.pdf", height = 5, width = 10)
print(lmer_plot)
dev.off()
####正文写作####
a1 <- Lmer_results3_all
a11 <- a1[a1$P_value <0.05,]
round(max(a11$Estimate),3)
round(min(a11$Estimate),3)
round(max(a11$P_value),3)
round(min(a11$P_value),3)

a2 <- Lmer_results5_all
a21 <- a2[a2$P_value <0.05,]
round(max(a21$Estimate),3)
round(min(a21$Estimate),3)
round(max(a21$P_value),3)
round(min(a21$P_value),3)
####菌群&ASM--lmer可视化####
Lmer_results8_all <- process_lmer(unique(Lmer_results3_sigFDR$Outcome), # Lmer_results3_sigFDR$Outcome c("s__Bifidobacterium_pseudocatenulatum")
                                 c("ASMI_z","ASM_z"), 
                                 data3_all_B_final, 
                                 c(Covariates_all_lmer, "Tea_freq"), 
                                 2)

Lmer_results8_all$P_FDR <-  p.adjust(Lmer_results8_all$P_value, method = "fdr")
Lmer_results8_all_sig <- Lmer_results8_all[Lmer_results8_all$P_value < 0.05,]
Lmer_results8_all_sigFDR <- Lmer_results8_all[Lmer_results8_all$P_FDR < 0.05,]

#*******************可视化
plot_df <- Lmer_results8_all %>%
  mutate(
    Direction = ifelse(Estimate > 0, "Positive", "Negative"),
    Predictor = gsub("^s__", "", Predictor),
    Signif_star = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

plot_df2 <- plot_df %>%
  mutate(
    Predictor_reorder = reorder_within(Predictor, Estimate, Outcome)
  )

plot <- ggplot(plot_df2,
               aes(x = Estimate,
                   y = Predictor_reorder,
                   color = Direction)) +
  
  ## 棒棒糖的“杆”（从 0 到 Estimate）
  geom_segment(
    aes(x = 0,
        xend = Estimate,
        y = Predictor_reorder,
        yend = Predictor_reorder),
    size = 1
  ) +
  
  ## 棒棒糖的“头”
  geom_point(size = 3) +
  
  ## 著性星号（贴在点旁边）
  geom_text(
    aes(
      label = Signif_star,
      x = Estimate + sign(Estimate) * 0.0006
    ),
    color = "red",
    size = 7,
    hjust = ifelse(plot_df2$Estimate > 0, 0, 1)
  ) +
  
  ## 0 参考线
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  ## 分面
  facet_wrap(~ Outcome, scales = "free_y") +
  scale_y_reordered() +
  
  ## 配色
  scale_color_manual(
    values = c(
      "Positive" = "#D55E00",
      "Negative" = "#0072B2"
    )
  ) +
  
  ## 标签
  labs(
    x = "Estimate",
    y = NULL,
    color = "Effect direction"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.15, 0.25))
  ) +
  ## 主题
  theme_myself

pdf("Micro&ASM.pdf", width = 14, height = 10)
print(plot)
dev.off()
####正文写作####
a3 <- Lmer_results8_all_sig
a3$Estimate <- round(a3$Estimate,3)
a3$P_value <- round(a3$P_value,3)

a31 <- a3[a3$Estimate >0,]
a32 <- a3[a3$Estimate <0,]

####中介效应分析####
exposures <- c("Tea_freq")
mediators <- c("GMI_Tea")
outcomes  <- c("ASMI_z","ASM_z")

med_results1_micro <- run_mediation_glmer_batch(
  data        = Micro_GMSI_all_ASMI,
  id_var      = "ID",
  exposures   = exposures,
  mediators   = mediators,
  outcomes    = outcomes,
  covariates  = Covariates_all_lmer,
  exposure_high = 3,
  exposure_low  = 0,
  sims = 1000
)

exposures <- c("serum_I_EGCG_F0_T")
mediators <- c("GMI_EGCG")
outcomes  <- c("ASMI_z","ASM_z")

med_results3_micro <- run_mediation_glmer_batch(
  data        = Micro_GMSI_all_ASMI,
  id_var      = "ID",
  exposures   = exposures,
  mediators   = mediators,
  outcomes    = outcomes,
  covariates  = Covariates_all_lmer,
  exposure_high = 3,
  exposure_low  = 0,
  sims = 1000
)

exposures <- c("serum_I_ECG_F0_T")
mediators <- c("GMI_ECG")
outcomes  <- c("ASMI_z","ASM_z")

med_results2_micro <- run_mediation_glmer_batch(
  data        = Micro_GMSI_all_ASMI,
  id_var      = "ID",
  exposures   = exposures,
  mediators   = mediators,
  outcomes    = outcomes,
  covariates  = Covariates_all_lmer,
  exposure_high = 3,
  exposure_low  = 0,
  sims = 1000
)


med_results_all_micro <- rbind(med_results1_micro, med_results2_micro, med_results3_micro)
####中介效应表格####
med_results_all_Pro_Micro <- rbind(med_results_all, med_results_all_micro)

med_results_all_Pro_Micro$Exposure <- recode(med_results_all_Pro_Micro$Exposure,
                                             "Tea_freq" = "Tea consumption frequency",
                                             "serum_I_epicatechin_F0_T" = "Epicatechin",
                                             "serum_I_EGC_F0_T" = "Epigallocatechin",
                                             "serum_I_EGCG_F0_T" = "Epigallocatechin gallate",
                                             "serum_I_ECG_F0_T" = "Epicatechin gallate")

med_results_all_Pro_Micro$Level[med_results_all_Pro_Micro$Exposure == "Tea consumption frequency"] <-  "More-than-daily vs. Non"
med_results_all_Pro_Micro$Level[med_results_all_Pro_Micro$Exposure %in% c("Catechin","Epicatechin","Epigallocatechin")] <-  "High vs. Undetectable"
med_results_all_Pro_Micro$Level[med_results_all_Pro_Micro$Exposure %in% c("Epigallocatechin gallate","Epicatechin gallate")] <-  "Quartile 4 vs. Quartile 1"

med_results_all_Pro_Micro$ACME <- round(med_results_all_Pro_Micro$ACME,3)
med_results_all_Pro_Micro$ACME_p <- round(med_results_all_Pro_Micro$ACME_p,3)
med_results_all_Pro_Micro$PropMediated <- round(med_results_all_Pro_Micro$PropMediated*100,1)
colnames(med_results_all_Pro_Micro)
med_results_all_Pro_Micro_output <- med_results_all_Pro_Micro[,c("Exposure","Level", "Mediator", "Outcome", "ACME", "ACME_p", "PropMediated")]
#writexl::write_xlsx(med_results_all_Pro_Micro_output, "med_results_all_Pro_Micro_output.xlsx")
####*****************************蛋白质分数&菌群分数*****************************####
PI_GMI <- merge(Pro_GMSI_all[Pro_GMSI_all$Times == "F2",], Micro_GMSI_all[Micro_GMSI_all$Times == "F2",], by = c("ID"))
PI_GMI <- merge(merge(Pro_GMSI_all, Micro_GMSI_all, by = c("ID","Times")), Age_long, by = c("ID","Times"))
PI_GMI_Cov <- merge(PI_GMI, Cov_F0, by = "ID")
colnames(PI_GMI_Cov)
a <- process_lmer(c("GMI_Tea","GMI_epicatechin","GMI_ECG","GMI_EGC","GMI_EGCG"), # Lmer_results3_sigFDR$Outcome c("s__Bifidobacterium_pseudocatenulatum")
                  c("PI_Tea", "PI_catechin", "PI_ECG", "PI_epicatechin"), 
                  PI_GMI_Cov, 
                  c("Age","Sex_F0"), 
                  2)

plot1 <- plot_corr_lower_triangle(
  data   = PI_GMI,
  vars   = setdiff(colnames(PI_GMI_F2),c("ID","Times.x", "Times.y","Times")),
  method = "spearman",
  title  = "F2"
)

####*****************************多样性计算*****************************####
####组间alpha多样性比较####
data_list <- list(Cov_F0, Tea_Follow, Catechin, Diversity[,c("ID","Times","shannon","simpson","pielou_evenness")])
data_alpha_all <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID")
data_alpha_all <- data_alpha_all[complete.cases(data_alpha_all),]
data_alpha_all <- data_alpha_all[order(data_alpha_all$ID, data_alpha_all$Times),]

# data_alpha_CatT4_all <- zero_tertile4(
#   data = data_alpha_all,
#   vars = Catechin_vars,
#   suffix = "_T"
# )
data_alpha_CatT4_all <- zero_tertile4(
  data = data_alpha_all,
  vars = Catechin_vars3,
  suffix = "_T"
)

data_alpha_CatT4_all <- data_alpha_CatT4_all %>%
  mutate(
    across(
      all_of(Catechin_vars2),
      ~ quartile_group(.x),
      .names = "{.col}_T"
    )
  )

tab_results_all <- setNames(
  lapply(Catechin_vars_T, function(v) table(data_alpha_CatT4_all[data_alpha_CatT4_all$Times == "F3",][[v]])),
  Catechin_vars_T
)



param_all <- tibble::tribble(
  ~Times, ~group_var, ~y_var,  ~xlab,  ~ylab, ~labels,
  "F2",   "Tea_freq", "shannon", "Tea consumption frequency",  "Shannon", c("Non","Less-than-daily","Daily","More-than-daily"),
  "F3",   "Tea_freq", "shannon", "Tea consumption frequency",  "Shannon", c("Non","Less-than-daily","Daily","More-than-daily"),
  
  "F2",   "serum_I_catechin_F0_T", "shannon", "Catechin", "Shannon", c("Undetectable","Low", "Medium","High"),
  "F3",   "serum_I_catechin_F0_T", "shannon", "Catechin", "Shannon", c("Undetectable","Low", "Medium","High"),
  
  "F2",   "serum_I_epicatechin_F0_T", "shannon", "Epicatechin", "Shannon", c("Undetectable","Low", "Medium","High"),
  "F3",   "serum_I_epicatechin_F0_T", "shannon", "Epicatechin", "Shannon", c("Undetectable","Low", "Medium","High"),
  
  "F2",   "serum_I_EGC_F0_T", "shannon", "Epigallocatechin", "Shannon",  c("Undetectable","Low", "Medium","High"),
  "F3",   "serum_I_EGC_F0_T", "shannon", "Epigallocatechin", "Shannon", c("Undetectable","Low", "Medium","High"),

  "F2",   "serum_I_EGCG_F0_T", "shannon", "Epigallocatechin gallate", "Shannon",  c("Quartile 1","Quartile 2","Quartile 3","Quartile 4"),
  "F3",   "serum_I_EGCG_F0_T", "shannon", "Epigallocatechin gallate", "Shannon", c("Quartile 1","Quartile 2","Quartile 3","Quartile 4"),

  "F2",   "serum_I_ECG_F0_T", "shannon", "Epicatechin gallate", "Shannon",  c("Quartile 1","Quartile 2","Quartile 3","Quartile 4"),
  "F3",   "serum_I_ECG_F0_T", "shannon", "Epicatechin gallate", "Shannon", c("Quartile 1","Quartile 2","Quartile 3","Quartile 4"),
)

plots_all <- param_all %>%
  mutate(
    plot = pmap(
      list(Times, group_var, y_var, xlab, ylab, labels),
      ~ {
        
        df_sub <- data_alpha_CatT4_all[
          data_alpha_CatT4_all$Times == ..1 &
            data_alpha_CatT4_all[[..2]] %in% c(0, 3),
        ]
        
        plot_wilcox_violin3(
          data = df_sub,
          group_var    = ..2,
          y_var        = ..3,
          group_levels = c(0, 3),
          group_labels = ..6[c(1, 4)],
          xlab  = ..4,
          ylab  = ..5,
          title = ..1
        )
      }
    )
  )

pdf("alpha_shannon.pdf", height = 14, width = 16)
wrap_plots(plots_all$plot, ncol = 4)
dev.off()


param_all <- tibble::tribble(
  ~Times, ~group_var, ~y_var,  ~xlab,  ~ylab, ~labels,
  "F2",   "Tea_freq", "simpson", "Tea consumption frequency",  "Simpson", c("Non","Less-than-daily","Daily","More-than-daily"),
  "F3",   "Tea_freq", "simpson", "Tea consumption frequency",  "Simpson", c("Non","Less-than-daily","Daily","More-than-daily"),
  
  "F2",   "serum_I_catechin_F0_T", "simpson", "Catechin", "Simpson", c("Undetectable","Low", "Medium","High"),
  "F3",   "serum_I_catechin_F0_T", "simpson", "Catechin", "Simpson", c("Undetectable","Low", "Medium","High"),
  
  "F2",   "serum_I_epicatechin_F0_T", "simpson", "Epicatechin", "Simpson", c("Undetectable","Low", "Medium","High"),
  "F3",   "serum_I_epicatechin_F0_T", "simpson", "Epicatechin", "Simpson", c("Undetectable","Low", "Medium","High"),
  
  "F2",   "serum_I_EGC_F0_T", "simpson", "Epigallocatechin", "Simpson",  c("Undetectable","Low", "Medium","High"),
  "F3",   "serum_I_EGC_F0_T", "simpson", "Epigallocatechin", "Simpson", c("Undetectable","Low", "Medium","High"),
  
  "F2",   "serum_I_EGCG_F0_T", "simpson", "Epigallocatechin gallate", "Simpson",  c("Quartile 1","Quartile 2","Quartile 3","Quartile 4"),
  "F3",   "serum_I_EGCG_F0_T", "simpson", "Epigallocatechin gallate", "Simpson", c("Quartile 1","Quartile 2", "Quartile 3", "Quartile 4"),
  
  "F2",   "serum_I_ECG_F0_T", "simpson", "Epicatechin gallate", "Simpson",  c("Quartile 1","Quartile 2", "Quartile 3", "Quartile 4"),
  "F3",   "serum_I_ECG_F0_T", "simpson", "Epicatechin gallate", "Simpson", c("Quartile 1","Quartile 2", "Quartile 3", "Quartile 4"),
)

plots_all <- param_all %>%
  mutate(
    plot = pmap(
      list(Times, group_var, y_var, xlab, ylab, labels),
      ~ {
        
        df_sub <- data_alpha_CatT4_all[
          data_alpha_CatT4_all$Times == ..1 &
            data_alpha_CatT4_all[[..2]] %in% c(0, 3),
        ]
        
        plot_wilcox_violin3(
          data = df_sub,
          group_var    = ..2,
          y_var        = ..3,
          group_levels = c(0, 3),
          group_labels = ..6[c(1, 4)],
          xlab  = ..4,
          ylab  = ..5,
          title = ..1
        )
      }
    )
  )

pdf("alpha_Simpson.pdf", height = 14, width = 16)
wrap_plots(plots_all$plot, ncol = 4)
dev.off()
####组间beta多样性比较####
# data_list <- list(Cov_F0[Cov_F0$Sex_F0 == 0,], Tea_Follow[,c("ID","Tea_freq")], Micro_all_F3)
# data_beta_all <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID")
# data_beta_all <- data_beta_all[complete.cases(data_beta_all),]
# data_beta_all <- data_beta_all[, unique(c(colnames(Tea_Follow[,c("ID","Tea_freq")]), colnames(Micro_all_F3)))]
plot_pcoa_permanova <- function(
    data,
    id_var,
    group_var,
    otu_vars = NULL, 
    group_levels = NULL,
    group_labels = NULL,     
    dist_method = "bray",
    n_perm = 999,
    title = NULL,
    base_size = 14
) {
  
  ## tidy evaluation
  id_var    <- rlang::ensym(id_var)
  group_var <- rlang::ensym(group_var)
  
  ## meta & otu
  meta <- data %>%
    dplyr::select(!!id_var, !!group_var)
  
  if (is.null(otu_vars)) {
    otu <- data %>%
      dplyr::select(-!!id_var, -!!group_var)
  } else {
    otu <- data %>%
      dplyr::select(all_of(otu_vars))
  }
  
  otu <- as.data.frame(lapply(otu, as.numeric))
  rownames(otu) <- meta[[rlang::as_string(id_var)]]
  
  ## 分组 factor（levels + labels）
  if (!is.null(group_levels)) {
    
    if (!is.null(group_labels)) {
      stopifnot(length(group_levels) == length(group_labels))
      
      meta[[rlang::as_string(group_var)]] <-
        factor(
          meta[[rlang::as_string(group_var)]],
          levels = group_levels,
          labels = group_labels
        )
    } else {
      meta[[rlang::as_string(group_var)]] <-
        factor(
          meta[[rlang::as_string(group_var)]],
          levels = group_levels
        )
    }
    
  } else {
    meta[[rlang::as_string(group_var)]] <-
      factor(meta[[rlang::as_string(group_var)]])
  }
  
  rownames(meta) <- meta[[rlang::as_string(id_var)]]
  
  ## 距离矩阵
  library(vegan)
  dist_mat <- vegdist(otu, method = dist_method)
  
  ## PERMANOVA
  set.seed(2028)
  adonis_res <- adonis2(
    dist_mat ~ meta[[rlang::as_string(group_var)]],
    permutations = n_perm
  )
  
  p_adonis  <- adonis_res$`Pr(>F)`[1]
  r2_adonis <- adonis_res$R2[1]
  
  lab_adonis <- if (p_adonis < 0.001) {
    "italic(P)<0.001"
  } else {
    paste0(
      "italic(P)==\"",
      sprintf("%.3f", p_adonis),
      "\""
    )
  }
  
  # lab_adonis <- paste0(
  #   "italic(P)==\"",
  #   sprintf("%.3f", p_adonis),
  #   "\""
  # )
  ## PCoA
  pcoa_res <- cmdscale(dist_mat, k = 2, eig = TRUE)
  
  pcoa_df <- data.frame(
    ID  = rownames(pcoa_res$points),
    PC1 = pcoa_res$points[, 1],
    PC2 = pcoa_res$points[, 2]
  ) %>%
    dplyr::left_join(meta, by = c("ID" = rlang::as_string(id_var)))
  
  var_exp <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 1)
  
  ## 作图
  p <- ggplot(
    pcoa_df,
    aes(
      x = PC1,
      y = PC2,
      color = .data[[rlang::as_string(group_var)]]
    )
  ) +
    geom_point(size = 3, alpha = 0.85) +
    stat_ellipse(type = "norm", level = 0.68, linewidth = 1) +
    annotate(
      "label",
      x = Inf, y = Inf,
      label = lab_adonis,
      parse = TRUE,
      size = 4.2,
      hjust = 1.9,
      vjust = 1.05,
      label.size = 0,      # 去掉边框
      fill = "grey95",     # 底色
      color = "black"
    ) +
    labs(
      title = title,
      x = paste0("PCoA1 (", var_exp[1], "%)"),
      y = paste0("PCoA2 (", var_exp[2], "%)"),
      color = NULL
    ) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0,
        size = base_size + 2
      ),
      axis.title   = element_text(face = "bold"),
      axis.text    = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.direction = "horizontal",
      legend.position = "bottom"
    )
  
  return(p)
}

species_data_all <- merge(data3_CatT4_all[,c("ID","Times",Catechin_vars_all)], data3_all_filter, by = c("ID","Times"))
species_names <- grep("^s__", colnames(species_data_all), value = TRUE)

species_data_all_F3 <- species_data_all[species_data_all$Times == "F3",]
species_data_all_F2 <- species_data_all[species_data_all$Times == "F2",]


plot_configs <- list(
  list(df = species_data_all_F2, var = "Tea_freq",                 labs = c("Non", "More-than-daily"), tit = "Tea frequency (F2)"),
  list(df = species_data_all_F3, var = "Tea_freq",                 labs = c("Non", "More-than-daily"), tit = "Tea frequency (F3)"),
  
  list(df = species_data_all_F2, var = "serum_I_catechin_F0_T",    labs = c("Undetectable", "High"),   tit = "Catechin (F2)"),
  list(df = species_data_all_F3, var = "serum_I_catechin_F0_T",    labs = c("Undetectable", "High"),   tit = "Catechin (F3)"),
  
  list(df = species_data_all_F2, var = "serum_I_epicatechin_F0_T", labs = c("Undetectable", "High"),   tit = "Epicatechin (F2)"),
  list(df = species_data_all_F3, var = "serum_I_epicatechin_F0_T", labs = c("Undetectable", "High"),   tit = "Epicatechin (F3)"),
  
  list(df = species_data_all_F2, var = "serum_I_EGC_F0_T", labs = c("Undetectable", "High"),   tit = "Epigallocatechin (F2)"),
  list(df = species_data_all_F3, var = "serum_I_EGC_F0_T", labs = c("Undetectable", "High"),   tit = "Epigallocatechin (F3)"),
  
  list(df = species_data_all_F2, var = "serum_I_EGCG_F0_T", labs = c("Undetectable", "High"),   tit = "Epigallocatechin gallate (F2)"),
  list(df = species_data_all_F3, var = "serum_I_EGCG_F0_T", labs = c("Undetectable", "High"),   tit = "Epigallocatechin gallate (F3)"),
  
  list(df = species_data_all_F2, var = "serum_I_ECG_F0_T", labs = c("Undetectable", "High"),   tit = "Epicatechin gallate (F2)"),
  list(df = species_data_all_F3, var = "serum_I_ECG_F0_T", labs = c("Undetectable", "High"),   tit = "Epicatechin gallate (F3)")
)


all_plots <- map(plot_configs, function(conf) {
  
  # 统一过滤 0 和 3 的逻辑
  plot_data <- conf$df %>% filter(.data[[conf$var]] %in% c(0, 3))
  
  plot_pcoa_permanova(
    data         = plot_data,
    id_var       = ID,                
    group_var    = !!sym(conf$var),   
    otu_vars     = species_names,
    group_levels = c(0, 3),
    group_labels = conf$labs,
    title        = conf$tit
  )
})

combined_plot <- wrap_plots(all_plots, ncol = 4) + 
  plot_layout(guides = "collect") &        # 合并图例
  theme(legend.position = "bottom")

pdf("Beta.pdf", width = 16, height = 12)
print(combined_plot)
dev.off()

####验证beta多样性####
# meta <- species_data_all %>%
#   dplyr::select(ID, Tea_freq)
# 
# otu <- species_data_all %>%
#   dplyr::select(-ID, -Tea_freq)
# 
# # 确保是数值型
# otu <- as.data.frame(lapply(otu, as.numeric))
# rownames(otu) <- species_data$ID
# 
# meta$Tea_freq <- factor(meta$Tea_freq, levels = c(0, 1, 2))
# rownames(meta) <- meta$ID
# 
# library(vegan)
# 
# bray_dist <- vegdist(otu, method = "bray")
# 
# adonis_res <- adonis2(
#   bray_dist ~ Tea_freq,
#   data = meta,
#   permutations = 999
# )
# 
# p_adonis  <- adonis_res$`Pr(>F)`[1]
# r2_adonis <- adonis_res$R2[1]
# 
# lab_adonis <- paste0(
#   "PERMANOVA (Bray–Curtis)\n",
#   "R² = ", round(r2_adonis, 3),
#   ", P = ", signif(p_adonis, 3)
# )
# 
# 
# disp <- betadisper(bray_dist, meta$Tea_freq)
# anova(disp)
# 
# 
# pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)
# 
# pcoa_df <- data.frame(
#   ID = rownames(pcoa_res$points),
#   PC1 = pcoa_res$points[,1],
#   PC2 = pcoa_res$points[,2]
# ) %>%
#   left_join(meta, by = "ID")
# 
# var_exp <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 1)
# 
# 
# library(ggplot2)
# 
# p_pcoa <- ggplot(pcoa_df, aes(PC1, PC2, color = Tea_freq)) +
#   geom_point(size = 3, alpha = 0.85) +
#   stat_ellipse(type = "norm", level = 0.68, linewidth = 1) +
#   annotate(
#     "text",
#     x = Inf, y = Inf,
#     label = lab_adonis,
#     hjust = 1.05, vjust = 1.3,
#     size = 4.2
#   ) +
#   labs(
#     x = paste0("PCoA1 (", var_exp[1], "%)"),
#     y = paste0("PCoA2 (", var_exp[2], "%)"),
#     color = "Tea frequency"
#   ) +
#   theme_classic(base_size = 14)
# 
# p_pcoa