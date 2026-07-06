.libPaths("E:/R_library")
# options(pkgType = "source")
# BiocManager::install("ropls")
# BiocManager::install("UniProt.ws")
# BiocManager::install("biomaRt",force = TRUE)
# remotes::install_github("hrbrmstr/ggchicklet")
# BiocManager::install("microbiome")
# BiocManager::install("factoextra")
# BiocManager::install("pacman")
# BiocManager::install("microeco")
# BiocManager::install("ggalluvial")
# BiocManager::install("sva")
library(ggrepel)
library(stringr)
library(dplyr)
library(purrr)
library(patchwork)
library(microbiome)
library(biomaRt)
#library(UniProt.ws)
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
#library(maaslin3)
library(pacman)
library(microeco)
library(ggpubr)
library(scales)
library(viridis)
library(patchwork)
library(cowplot)
#library(ggsankey)
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
library(MESS)
library(ggalluvial)
library(stringr)
library(forcats)
library(WeightIt)
library(cobalt)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(performance)
library(purrr)
library(stringr)
library(tibble)
library(stringr)
library(forcats)
library(patchwork)

rm(list=ls())
setwd("D:/OneDrive/Papers/9 饮茶频率&骨骼肌/Output")
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
    exclude_cols = c("ID", "Tea_freq", "Times")
) {
  
  library(dplyr)
  library(microbiome)
  
  ## 需要保留的元信息列
  meta_cols <- intersect(c("ID", "Times"), colnames(species_data))
  
  ## 原始物种数据
  species_raw <- species_data %>%
    dplyr::select(-all_of(exclude_cols))
  
  ## CLR 转换矩阵
  species_clr <- species_raw %>%
    as.matrix()
  
  species_clr <- species_clr + 1e-6
  
  species_clr <- microbiome::transform(species_clr, "clr") %>%
    as.data.frame()
  
  ## CLR变量添加后缀
  colnames(species_clr) <- paste0(colnames(species_clr), "_clr")
  
  ## 合并：
  ## ID/Times + 原始变量 + CLR变量
  species_clr_out <- bind_cols(
    species_data[, meta_cols, drop = FALSE],
    species_raw,
    species_clr
  )
  
  return(species_clr_out)
}
#******函数分位
group_variables4 <- function(data,
                            zero_vars = NULL,
                            quartile_vars = NULL,
                            suffix = "_G") {
  
  data %>%
    
    ## 0 单独一组
    dplyr::mutate(
      dplyr::across(
        all_of(zero_vars),
        ~ {
          x <- .
          
          out <- rep(NA_integer_, length(x))
          
          ## 0 单独一组
          out[x == 0] <- 0
          
          ## 非0部分
          nz <- x != 0 & !is.na(x)
          
          if(sum(nz) > 0) {
            
            q <- quantile(
              x[nz],
              probs = c(1/3, 2/3),
              na.rm = TRUE
            )
            
            out[nz & x <= q[1]] <- 1
            out[nz & x > q[1] & x <= q[2]] <- 2
            out[nz & x > q[2]] <- 3
          }
          
          out
        },
        .names = paste0("{.col}", suffix)
      )
    ) %>%
    
    ## 普通四分位
    dplyr::mutate(
      dplyr::across(
        all_of(quartile_vars),
        ~ {
          x <- .
          
          out <- rep(NA_integer_, length(x))
          
          q <- quantile(
            x,
            probs = c(0.25, 0.5, 0.75),
            na.rm = TRUE
          )
          
          out[x <= q[1]] <- 0
          out[x > q[1] & x <= q[2]] <- 1
          out[x > q[2] & x <= q[3]] <- 2
          out[x > q[3]] <- 3
          
          out
        },
        .names = paste0("{.col}", suffix)
      )
    )
}

group_variables3 <- function(data,
                                 zero_vars = NULL,
                                 tertile_vars = NULL,
                                 suffix = "_T") {
  
  data %>%
    
    ## 0 单独一类 + 非0二分位
    dplyr::mutate(
      dplyr::across(
        all_of(zero_vars),
        ~ {
          x <- .
          
          out <- rep(NA_integer_, length(x))
          
          ## 0 单独一组
          out[x == 0] <- 0
          
          ## 非0部分
          nz <- x != 0 & !is.na(x)
          
          if(sum(nz) > 0) {
            
            q <- quantile(
              x[nz],
              probs = 0.5,
              na.rm = TRUE
            )
            
            out[nz & x <= q] <- 1
            out[nz & x > q] <- 2
          }
          
          out
        },
        .names = paste0("{.col}", suffix)
      )
    ) %>%
    
    ## 普通三分位
    dplyr::mutate(
      dplyr::across(
        all_of(tertile_vars),
        ~ {
          x <- .
          
          out <- rep(NA_integer_, length(x))
          
          q <- quantile(
            x,
            probs = c(1/3, 2/3),
            na.rm = TRUE
          )
          
          out[x <= q[1]] <- 0
          out[x > q[1] & x <= q[2]] <- 1
          out[x > q[2]] <- 2
          
          out
        },
        .names = paste0("{.col}", suffix)
      )
    )
}
#******查看数据框中各变量情况
check_df_summary <- function(data) {
  
  library(dplyr)
  
  # 判断变量类型
  get_var_type <- function(x) {
    if (is.numeric(x)) {
      return("定量变量")
    } else {
      return("定性变量")
    }
  }
  
  # 检测特殊字符取值（非纯数字的字符值，以及 99、999 等特殊编码）
  detect_special_chars <- function(x) {
    
    # 定义需要标记为特殊的数值编码
    special_codes <- c("99", "999")
    
    if (!is.numeric(x)) {
      # 对于已经是字符/因子的变量
      unique_vals <- unique(as.character(x))
      unique_vals <- unique_vals[!is.na(unique_vals)]
      
      # 找出不能转为数值的值，或者属于特殊编码的值
      special_vals <- unique_vals[
        is.na(suppressWarnings(as.numeric(unique_vals))) | 
          unique_vals %in% special_codes
      ]
      
      if (length(special_vals) > 0) {
        return(paste0(special_vals, collapse = "; "))
      }
    } else {
      # 对于数值型变量
      x_char <- as.character(x)
      unique_vals <- unique(x_char)
      unique_vals <- unique_vals[!is.na(unique_vals)]
      
      # 找出不能转为数值的值，或者属于特殊编码的值
      special_vals <- unique_vals[
        is.na(suppressWarnings(as.numeric(unique_vals))) | 
          unique_vals %in% special_codes
      ]
      
      if (length(special_vals) > 0) {
        return(paste0(special_vals, collapse = "; "))
      }
    }
    return(NA)
  }
  
  # 主循环
  result_list <- lapply(names(data), function(var) {
    
    x <- data[[var]]
    
    # 基本信息
    n_total   <- length(x)
    n_missing <- sum(is.na(x))
    n_zero    <- if (is.numeric(x)) sum(x == 0, na.rm = TRUE) else NA
    
    var_type <- get_var_type(x)
    
    # 检测特殊字符取值
    special_values <- detect_special_chars(x)
    
    # ===== 定量变量
    if (is.numeric(x)) {
      
      res <- data.frame(
        Variable       = var,
        Type           = var_type,
        N              = n_total,
        Missing_N      = n_missing,
        Missing_P      = round(n_missing / n_total * 100, 2),
        Zero_N         = n_zero,
        Mean           = round(mean(x, na.rm = TRUE), 3),
        SD             = round(sd(x, na.rm = TRUE), 3),
        Median         = round(median(x, na.rm = TRUE), 3),
        Q1             = round(quantile(x, 0.25, na.rm = TRUE), 3),
        Q3             = round(quantile(x, 0.75, na.rm = TRUE), 3),
        Min            = round(min(x, na.rm = TRUE), 3),
        Max            = round(max(x, na.rm = TRUE), 3),
        Special_Values = special_values,
        Levels         = NA
      )
      
    } else {
      
      # ===== 定性变量
      tb <- table(x, useNA = "ifany")
      
      level_info <- paste0(
        names(tb),
        " (n=",
        as.vector(tb),
        ")",
        collapse = "; "
      )
      
      res <- data.frame(
        Variable       = var,
        Type           = var_type,
        N              = n_total,
        Missing_N      = n_missing,
        Missing_P      = round(n_missing / n_total * 100, 2),
        Zero_N         = NA,
        Mean           = NA,
        SD             = NA,
        Median         = NA,
        Q1             = NA,
        Q3             = NA,
        Min            = NA,
        Max            = NA,
        Special_Values = special_values,
        Levels         = level_info
      )
    }
    
    return(res)
  })
  
  # 合并结果
  result <- bind_rows(result_list)
  
  return(result)
}

#******对列的一些特殊取值的处理
clean_var <- function(data, col) {
  
  data[[col]] <- data[[col]] %>%
    as.character() %>%
    trimws() %>%
    
    # 删除所有空格
    gsub("\\s+", "", .) %>%
    
    # 中文句号 → .
    gsub("。", ".", .) %>%
    # 中文逗号 → .
    gsub("，", ".", .) %>%
    
    # ;. → .
    gsub(";\\.", ".", .) %>%
    
    # 分号 → .
    gsub(";", ".", .) %>%
    
    # 单引号删除
    gsub("'", "", .) %>%
    
    # 去掉结尾多余 .
    gsub("\\.$", "", .)
  
  return(data)
}
#******中介效应数据导出的处理
format_mediation_table <- function(
    data,
    digits_acme = 3,
    digits_p = 3,
    digits_prop = 1,
    prop_as_percent = TRUE,
    include_fdr = NULL,          # 是否显示FDR；NULL表示自动判断
    fdr_col = "ACME_FDR"         # FDR列名
){
  
  library(dplyr)
  
  ## 自动判断是否包含FDR列
  if(is.null(include_fdr)){
    include_fdr <- fdr_col %in% colnames(data)
  }
  
  ## 如果要求显示FDR，但数据中没有对应列，则报错
  if(include_fdr && !(fdr_col %in% colnames(data))){
    stop(
      paste0(
        "include_fdr = TRUE，但数据中没有列：", fdr_col,
        "。请设置 include_fdr = FALSE，或检查 fdr_col 参数。"
      )
    )
  }
  
  ## P值格式化函数
  format_p <- function(p){
    ifelse(
      is.na(p),
      NA_character_,
      ifelse(
        p < 0.001,
        "<0.001",
        sprintf(
          paste0("%.", digits_p, "f"),
          p
        )
      )
    )
  }
  
  ## 基础数据整理：这里不要引用 ACME_FDR
  data2 <- data %>%
    dplyr::mutate(
      ACME = as.numeric(ACME),
      ACME_low = as.numeric(ACME_low),
      ACME_high = as.numeric(ACME_high),
      ACME_p = as.numeric(ACME_p),
      PropMediated = as.numeric(PropMediated)
    )
  
  ## 如果需要FDR，再单独添加 ACME_FDR_tmp
  if(include_fdr){
    data2 <- data2 %>%
      dplyr::mutate(
        ACME_FDR_tmp = as.numeric(.data[[fdr_col]])
      )
  }
  
  ## 输出基础表
  out <- data2 %>%
    dplyr::transmute(
      
      Exposure,
      Mediator,
      Outcome,
      
      `ACME (95% CI)` =
        sprintf(
          paste0(
            "%.", digits_acme, "f (%.", digits_acme,
            "f, %.", digits_acme, "f)"
          ),
          ACME,
          ACME_low,
          ACME_high
        ),
      
      `p` = format_p(ACME_p),
      
      `Proportion mediated (%)` =
        if(prop_as_percent){
          sprintf(
            paste0("%.", digits_prop, "f"),
            PropMediated * 100
          )
        } else {
          sprintf(
            paste0("%.", digits_prop, "f"),
            PropMediated
          )
        }
    )
  
  ## 如果需要FDR，则插入 p_FDR 列
  if(include_fdr){
    out <- out %>%
      dplyr::mutate(
        `p_FDR` = format_p(data2$ACME_FDR_tmp),
        .after = "p"
      )
  }
  
  return(out)
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
calculate_GMSI <- function(
    data,
    outcomes_pos,
    outcomes_neg,
    GMSI_name,
    time_var = "Times",
    remove_suffix = "_clr_z$",
    epsilon = 1e-8
) {
  
  ## 1. 提取菌名

  MG <- unique(
    gsub(remove_suffix, "", outcomes_pos)
  )
  
  MP <- unique(
    gsub(remove_suffix, "", outcomes_neg)
  )
  
  taxa_all <- unique(c(MG, MP))
  
  ## 2. 判断是否为长数据

  has_time <- time_var %in% colnames(data)
  
  ## 3. 保留变量

  keep_cols <- c("ID", taxa_all)
  
  if (has_time) {
    
    keep_cols <- c(
      "ID",
      time_var,
      taxa_all
    )
  }
  
  keep_cols <- unique(keep_cols)
  
  ## 检查变量是否存在
  missing_taxa <- setdiff(
    taxa_all,
    colnames(data)
  )
  
  if (length(missing_taxa) > 0) {
    
    warning(
      paste(
        "The following taxa are missing:",
        paste(missing_taxa,
              collapse = ", ")
      )
    )
  }
  
  keep_cols <- intersect(
    keep_cols,
    colnames(data)
  )
  
  data_GMSI <- data[
    ,
    keep_cols,
    drop = FALSE
  ]
  
  ## 4. abundance matrix

  abund_mat <- as.matrix(
    data_GMSI[
      ,
      intersect(
        taxa_all,
        colnames(data_GMSI)
      ),
      drop = FALSE
    ]
  )
  
  ## 转数值
  abund_mat <- apply(
    abund_mat,
    2,
    as.numeric
  )
  
  abund_mat <- as.matrix(abund_mat)
  
  ## 5. 行名

  rownames(abund_mat) <- if (has_time) {
    
    paste(
      data_GMSI$ID,
      data_GMSI[[time_var]],
      sep = "_"
    )
    
  } else {
    
    data_GMSI$ID
  }
  
  ## 6. 加权丰度函数

  calc_weighted_abundance <- function(
    x,
    taxa_set
  ) {
    
    taxa_use <- intersect(
      taxa_set,
      names(x)
    )
    
    ## 无有效菌
    if (length(taxa_use) == 0) {
      
      return(list(
        mean_abund = NA_real_,
        presence_ratio = NA_real_,
        weighted_abund = NA_real_,
        n_taxa = 0
      ))
    }
    
    ## 丰度
    abund_vec <- as.numeric(
      x[taxa_use]
    )
    
    ## 平均丰度
    mean_abund <- mean(
      abund_vec,
      na.rm = TRUE
    )
    
    ## presence ratio
    ## 原始数据中 >0 代表存在
    presence_ratio <- mean(
      abund_vec > 0,
      na.rm = TRUE
    )
    
    ## weighted abundance
    weighted_abund <- (
      mean_abund *
        presence_ratio
    )
    
    return(list(
      mean_abund = mean_abund,
      presence_ratio = presence_ratio,
      weighted_abund = weighted_abund,
      n_taxa = length(taxa_use)
    ))
  }
  
  ## 7. 逐样本计算

  results_list <- lapply(
    
    seq_len(nrow(abund_mat)),
    
    function(i) {
      
      x <- abund_mat[i, ]
      
      ## positive taxa
      MG_res <- calc_weighted_abundance(
        x = x,
        taxa_set = MG
      )
      
      ## negative taxa
      MP_res <- calc_weighted_abundance(
        x = x,
        taxa_set = MP
      )
      
      ## 计算 GMSI

      ## 两边都有
      if (
        length(MG) > 0 &
        length(MP) > 0
      ) {
        
        GMSI_value <- log10(
          (
            MG_res$weighted_abund +
              epsilon
          ) /
            (
              MP_res$weighted_abund +
                epsilon
            )
        )
        
        ## 只有 positive
      } else if (
        length(MG) > 0 &
        length(MP) == 0
      ) {
        
        GMSI_value <- log10(
          MG_res$weighted_abund +
            epsilon
        )
        
        ## 只有 negative
      } else if (
        length(MG) == 0 &
        length(MP) > 0
      ) {
        
        GMSI_value <- -log10(
          MP_res$weighted_abund +
            epsilon
        )
        
      } else {
        
        GMSI_value <- NA_real_
      }
      
      ## 输出

      data.frame(
        
        MG_n_taxa =
          MG_res$n_taxa,
        
        MG_mean_abund =
          MG_res$mean_abund,
        
        MG_presence_ratio =
          MG_res$presence_ratio,
        
        MG_weighted_abund =
          MG_res$weighted_abund,
        
        MP_n_taxa =
          MP_res$n_taxa,
        
        MP_mean_abund =
          MP_res$mean_abund,
        
        MP_presence_ratio =
          MP_res$presence_ratio,
        
        MP_weighted_abund =
          MP_res$weighted_abund,
        
        GMSI = GMSI_value
      )
    }
  )
  
  ## 8. 合并结果

  results_df <- do.call(
    rbind,
    results_list
  )
  
  ## 9. 返回数据框

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
  
  ## 合并结果
  GMSI_df <- cbind(
    GMSI_df,
    results_df
  )
  
  ## 修改列名
  colnames(GMSI_df)[
    colnames(GMSI_df) == "GMSI"
  ] <- GMSI_name
  
  return(GMSI_df)
}


calculate_signature_index <- function(
    data,
    outcomes_pos,
    outcomes_neg,
    index_name,
    time_var = "Times",
    remove_suffix = NULL,
    method = c(
      "difference",
      "log_ratio"
    ),
    epsilon = 1e-8
) {
  
  library(dplyr)
  
  ## method

  method <- match.arg(method)
  
  ## 去后缀

  if (!is.null(remove_suffix)) {
    
    pos_vars <- unique(
      gsub(remove_suffix, "", outcomes_pos)
    )
    
    neg_vars <- unique(
      gsub(remove_suffix, "", outcomes_neg)
    )
    
  } else {
    
    pos_vars <- unique(outcomes_pos)
    neg_vars <- unique(outcomes_neg)
  }
  
  ## 全部变量

  feature_all <- unique(
    c(pos_vars, neg_vars)
  )
  
  ## 是否long format

  has_time <- time_var %in% colnames(data)
  
  ## 保留列

  keep_cols <- c("ID", feature_all)
  
  if (has_time) {
    
    keep_cols <- c(
      "ID",
      time_var,
      feature_all
    )
  }
  
  keep_cols <- unique(keep_cols)
  
  df <- data[, keep_cols, drop = FALSE]
  
  ## abundance matrix

  mat <- as.matrix(
    df[, feature_all, drop = FALSE]
  )
  
  ## 行名

  rownames(mat) <- if (has_time) {
    
    paste(
      df$ID,
      df[[time_var]],
      sep = "_"
    )
    
  } else {
    
    df$ID
  }
  
  ## 单样本计算函数

  calc_group_score <- function(
    x,
    taxa_set
  ) {
    
    taxa_use <- intersect(
      taxa_set,
      names(x)
    )
    
    ## 没有变量
    if (length(taxa_use) == 0) {
      
      return(list(
        mean_value = NA_real_,
        presence_ratio = NA_real_
      ))
    }
    
    values <- x[taxa_use]
    
    ## 平均值
    mean_value <- mean(
      values,
      na.rm = TRUE
    )
    
    ## 存在率
    presence_ratio <- mean(
      values > 0,
      na.rm = TRUE
    )
    
    return(list(
      mean_value = mean_value,
      presence_ratio = presence_ratio
    ))
  }
  
  ## 逐行计算

  results_list <- lapply(
    
    seq_len(nrow(mat)),
    
    function(i) {
      
      x <- mat[i, ]
      
      pos_res <- calc_group_score(
        x,
        pos_vars
      )
      
      neg_res <- calc_group_score(
        x,
        neg_vars
      )
      
      ## 指数计算

      ## 只有正向signature
      if (length(neg_vars) == 0) {
        
        index_value <- pos_res$mean_value
        
        ## 只有负向signature
      } else if (length(pos_vars) == 0) {
        
        index_value <- -neg_res$mean_value
        
      } else {
        
        ## Difference
        ## 适合：
        ## CLR / z-score /
        ## proteomics /
        ## metabolomics

        if (method == "difference") {
          
          index_value <-
            pos_res$mean_value -
            neg_res$mean_value
        }
        
        ## log-ratio
        ## 适合：
        ## 原始relative abundance

        if (method == "log_ratio") {
          
          index_value <- log10(
            (pos_res$mean_value + epsilon) /
              (neg_res$mean_value + epsilon)
          )
        }
      }
      
      data.frame(
        
        POS_mean =
          pos_res$mean_value,
        
        POS_presence =
          pos_res$presence_ratio,
        
        NEG_mean =
          neg_res$mean_value,
        
        NEG_presence =
          neg_res$presence_ratio,
        
        INDEX =
          index_value
      )
    }
  )
  
  ## 合并结果

  results_df <- do.call(
    rbind,
    results_list
  )
  
  ## 输出

  if (has_time) {
    
    out_df <- data.frame(
      ID = df$ID,
      Times = df[[time_var]],
      stringsAsFactors = FALSE
    )
    
  } else {
    
    out_df <- data.frame(
      ID = df$ID,
      stringsAsFactors = FALSE
    )
  }
  
  out_df <- cbind(
    out_df,
    results_df
  )
  
  ## 修改指数列名
  colnames(out_df)[
    colnames(out_df) == "INDEX"
  ] <- index_name
  
  return(out_df)
}


#******计算轨迹
run_hlme_traj_plot <- function(
    data,
    outcome,
    time_var,
    id_var,
    max_ng = 5,
    ng_plot = 3,
    sex_label = NULL,
    xlab_label = "Time",
    ylab_label = "Outcome",
    legend.x,
    legend.y,
    pdf_file = NULL,
    width = 6,
    height = 6,
    class_names = NULL,
    legend_order = NULL,
    colors = c(
      Class1 = "#1f77b4", Class2 = "#ff7f0e", Class3 = "#2ca02c",
      Class4 = "#d62728", Class5 = "#9467bd"
    ),
    linetypes = c(
      Class1 = 1, Class2 = 2, Class3 = 3, Class4 = 4, Class5 = 5
    )
){
  
  library(lcmm)
  library(dplyr)
  
  data <- as.data.frame(data)
  
  fixed_f   <- reformulate(time_var, response = paste0("`", outcome, "`"))
  mixture_f <- reformulate(time_var)
  
  ## 1 模型拟合
  
  model_list <- list()
  
  model_list[[1]] <- hlme(
    fixed = fixed_f,
    subject = id_var,
    ng = 1,
    data = data
  )
  
  model_list[[1]]$call$fixed <- fixed_f
  
  if (max_ng >= 2) {
    
    for (k in 2:max_ng) {
      
      model_list[[k]] <- hlme(
        fixed = fixed_f,
        mixture = mixture_f,
        subject = id_var,
        ng = k,
        data = data,
        B = model_list[[1]]
      )
      
      model_list[[k]]$call$fixed <- fixed_f
      model_list[[k]]$call$mixture <- mixture_f
      
    }
    
  }
  
  ## 2 模型评价
  
  summary_tab <- do.call(summarytable, c(model_list))
  
  model_plot <- model_list[[ng_plot]]
  
  pprob_data <- as.data.frame(model_plot$pprob)
  
  ## 平均后验概率
  
  avg_pp <- sapply(1:ng_plot, function(i){
    
    sub_probs <- pprob_data[pprob_data$class == i, paste0("prob", i)]
    
    if(length(sub_probs) > 0){
      round(mean(sub_probs), 4)
    }else{
      0
    }
    
  })
  
  names(avg_pp) <- paste0("Class", 1:ng_plot)
  
  ## 类别比例
  
  class_counts <- table(factor(pprob_data$class, levels = 1:ng_plot))
  
  class_percent <- round(as.numeric(class_counts) / sum(class_counts) * 100, 1)
  
  ## 类别名称
  
  default_names <- paste0("Class ", 1:ng_plot)
  
  if(!is.null(class_names)){
    
    for(i in seq_along(class_names)){
      
      class_id <- as.numeric(gsub("Class", "", names(class_names)[i]))
      
      if(class_id <= ng_plot){
        default_names[class_id] <- class_names[i]
      }
      
    }
    
  }
  
  class_labels <- paste0(default_names, " (", class_percent, "%)")
  
  ## 颜色与线型
  
  class_cols <- colors[paste0("Class", 1:ng_plot)]
  class_lty  <- linetypes[paste0("Class", 1:ng_plot)]
  
  ## legend顺序
  
  legend_labels <- class_labels
  legend_cols   <- class_cols
  legend_lty    <- class_lty
  
  if(!is.null(legend_order)){
    
    legend_labels <- class_labels[legend_order]
    legend_cols   <- class_cols[legend_order]
    legend_lty    <- class_lty[legend_order]
    
  }
  
  ## 3 绘图
  
  if (!is.null(pdf_file)) {
    pdf(pdf_file, width = width, height = height)
  }
  
  time_range <- range(data[[time_var]], na.rm = TRUE)
  
  newdata <- data.frame(seq(time_range[1], time_range[2], length.out = 100))
  colnames(newdata) <- time_var
  
  x_ticks <- sort(unique(data[[time_var]]))
  
  predY <- predictY(model_plot, newdata, var.time = time_var, draws = TRUE)
  
  op <- par(mar = c(5,4,3,2)+0.1)
  
  ## 基础参数
  
  plot_args <- list(
    
    x = predY,
    ylab = ylab_label,
    xlab = xlab_label,
    main = "",
    lwd = 2,
    shade = TRUE,
    col = class_cols,
    lty = class_lty,
    legend = NULL, # 控制原有的legend显不显示
    xaxt = "n",
    font.lab = 2,
    cex.lab = 1.2
  )
  
  
  do.call(plot, plot_args)
  
  ## X轴
  
  axis(
    side = 1,
    at = x_ticks,
    labels = x_ticks,
    cex.axis = 0.9
  )
  
  ## 顶部标签
  
  if(!is.null(sex_label)){
    
    mtext(
      sex_label,
      side = 3,
      adj = 0,
      line = 0.5,
      font = 2,
      cex = 1.2
    )
    
  }
  
  ## 自定义legend
  
  legend(
    x = legend.x,
    y = legend.y,
    legend = legend_labels,
    col = legend_cols,
    lty = legend_lty,
    lwd = 3,
    seg.len = 2,
    bty = "n",
    cex = 1.0
  )
  
  if (!is.null(pdf_file)) {
    dev.off()
  }
  
  par(op)
  
  ## 4 输出
  
  pprob_data$class_label <- default_names[pprob_data$class]
  
  return(
    list(
      summary_table = as.data.frame(summary_tab),
      selected_avg_pp = avg_pp,
      model_plot = model_plot,
      Class = pprob_data
    )
  )
  
}
####数据分析####
process_lmer <- function(X, Y, data,
                         covariates = NULL,
                         gender_label) {
  
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
      x_safe <- safe_var(x)
      y_safe <- safe_var(y)
      
      # 协变量可有可无
      if (!is.null(covariates) && length(covariates) > 0) {
        
        cov_safe <- sapply(covariates, safe_var)
        
        formula_text <- paste(
          y_safe, "~",
          x_safe, "+",
          paste(cov_safe, collapse = " + "),
          "+ (1|ID)"
        )
        
      } else {
        
        formula_text <- paste(
          y_safe, "~",
          x_safe,
          "+ (1|ID)"
        )
      }
      
      formula_full <- as.formula(formula_text)
      
      # ========= 运行模型
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
      
      # 匹配该变量所有行
      coef_x <- coef_df[startsWith(coef_df$Term_clean, x), ]
      
      if (nrow(coef_x) == 0) next
      
      # ========= 提取 Level
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
        coef_x[, c(
          "Outcome", "Predictor",
          "Variable_Type",
          "Level",
          "Estimate",
          "Std. Error",
          "CI_low",
          "CI_high",
          "P_value",
          "Significance",
          "Gender",
          "Warning_Message",
          "Singular_Fit",
          "Random_Intercept_Var",
          "Residual_Var"
        )]
      
      counter <- counter + 1
    }
  }
  
  bind_rows(results_list)
}

#************性别交互
process_lmer_interaction <- function(
    X,
    Y,
    data,
    covariates = NULL,
    interaction_var = "Sex_F0"
) {
  
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
      int_safe <- safe_var(interaction_var)
      
      # 去掉 interaction_var，避免重复进入模型
      cov_use <- setdiff(covariates, interaction_var)
      
      # 协变量
      if (!is.null(cov_use) && length(cov_use) > 0) {
        
        cov_safe <- sapply(cov_use, safe_var)
        
        formula_text <- paste(
          y_safe, "~",
          
          # 交互项
          paste0(x_safe, " * ", int_safe),
          
          "+",
          paste(cov_safe, collapse = " + "),
          
          "+ (1|ID)"
        )
        
      } else {
        
        formula_text <- paste(
          y_safe, "~",
          paste0(x_safe, " * ", int_safe),
          "+ (1|ID)"
        )
      }
      
      formula_full <- as.formula(formula_text)
      
      # ========= 运行模型
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
      
      # 只提取：
      # x主效应
      # x:Sex_F0交互项
      
      coef_x <- coef_df[
        grepl(
          paste0(
            "^", x, "$|",
            "^", x, ".*:", interaction_var, "|",
            "^", interaction_var, ".*:", x
          ),
          coef_df$Term_clean
        ),
      ]
      
      if (nrow(coef_x) == 0) next
      
      # ========= 判断交互 or 主效应
      coef_x$Effect_Type <- ifelse(
        grepl(":", coef_x$Term_clean),
        "Interaction",
        "Main Effect"
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
      
      # ========= 存储
      results_list[[counter]] <-
        coef_x[, c(
          "Outcome",
          "Predictor",
          "Variable_Type",
          "Term_clean",
          "Effect_Type",
          "Estimate",
          "Std. Error",
          "CI_low",
          "CI_high",
          "P_value",
          "Significance",
          "Warning_Message",
          "Singular_Fit",
          "Random_Intercept_Var",
          "Residual_Var"
        )]
      
      counter <- counter + 1
    }
  }
  
  bind_rows(results_list)
}

#************glm分析
process_glm1 <- function(X, Y, data, covariates = NULL, gender_label = 2,
                         min_n = 10) {
  
  library(dplyr)
  
  # 加反引号，保护特殊变量名
  bt <- function(v) paste0("`", v, "`")
  
  # 去掉反引号，用于匹配
  remove_bt <- function(v) gsub("`", "", v, fixed = TRUE)
  
  # 正则转义
  regex_escape <- function(x) {
    gsub("([][{}()+*^$.|\\\\?\\-])", "\\\\\\1", x)
  }
  
  # 去重
  X <- unique(X)
  Y <- unique(Y)
  covariates <- unique(covariates)
  
  # 检查变量是否存在
  all_vars <- unique(c(X, Y, covariates))
  missing_vars <- setdiff(all_vars, colnames(data))
  
  if (length(missing_vars) > 0) {
    stop(
      paste0(
        "These variables are not in data:\n",
        paste(missing_vars, collapse = "\n")
      )
    )
  }
  
  gender_text <- switch(
    as.character(gender_label),
    "0" = "Female",
    "1" = "Male",
    "2" = "All",
    NA
  )
  
  results_list <- list()
  counter <- 1
  
  for (x in X) {
    for (y in Y) {
      
      # 当前模型需要的变量
      rhs_vars <- unique(c(x, covariates))
      vars_needed <- unique(c(y, rhs_vars))
      
      dat_sub <- data[, vars_needed, drop = FALSE]
      dat_sub <- na.omit(dat_sub)
      
      if (nrow(dat_sub) < min_n) next
      
      # 构建公式
      formula_text <- paste(
        bt(y), "~",
        paste(bt(rhs_vars), collapse = " + ")
      )
      
      fml <- as.formula(formula_text)
      
      # 建模
      model <- try(
        glm(fml, data = dat_sub, family = gaussian()),
        silent = TRUE
      )
      
      if (inherits(model, "try-error")) next
      
      # 提取模型系数
      coef_df <- as.data.frame(summary(model)$coefficients)
      coef_df$Term <- rownames(coef_df)
      coef_df$Term_clean <- remove_bt(coef_df$Term)
      
      # 用 model.matrix 自动识别 x 对应的模型项
      mm <- model.matrix(fml, data = dat_sub)
      assign_vec <- attr(mm, "assign")
      term_labels <- attr(terms(fml), "term.labels")
      term_labels_clean <- remove_bt(term_labels)
      
      x_index <- which(term_labels_clean == x)
      
      if (length(x_index) == 0) next
      
      x_terms <- colnames(mm)[assign_vec %in% x_index]
      x_terms <- setdiff(x_terms, "(Intercept)")
      
      # 先用原始 term 精确匹配
      coef_x <- coef_df %>%
        filter(Term %in% x_terms)
      
      # 如果因为反引号问题没有匹配到，再用 clean term 匹配
      if (nrow(coef_x) == 0) {
        x_terms_clean <- remove_bt(x_terms)
        coef_x <- coef_df %>%
          filter(Term_clean %in% x_terms_clean)
      }
      
      # 兜底：用前缀匹配
      if (nrow(coef_x) == 0) {
        x_esc <- regex_escape(x)
        coef_x <- coef_df %>%
          filter(grepl(paste0("^", x_esc), Term_clean))
      }
      
      if (nrow(coef_x) == 0) next
      
      # 判断 P 值列名
      p_col <- grep("^Pr\\(", colnames(coef_x), value = TRUE)[1]
      
      # 生成 Level
      x_esc <- regex_escape(x)
      
      coef_x <- coef_x %>%
        mutate(
          Outcome = y,
          Predictor = x,
          Level = gsub(paste0("^", x_esc), "", Term_clean),
          Level = ifelse(Level == "" | Term_clean == x, "Continuous", Level),
          CI_low = Estimate - 1.96 * `Std. Error`,
          CI_high = Estimate + 1.96 * `Std. Error`,
          P_value = .data[[p_col]],
          Significance = ifelse(P_value < 0.05, "Significant", "Not Significant"),
          Gender = gender_text,
          N = nrow(dat_sub)
        ) %>%
        dplyr::select(
          Outcome,
          Predictor,
          Level,
          Estimate,
          `Std. Error`,
          CI_low,
          CI_high,
          P_value,
          Significance,
          Gender,
          N
        )
      
      results_list[[counter]] <- coef_x
      counter <- counter + 1
    }
  }
  
  if (length(results_list) == 0) {
    return(data.frame())
  }
  
  results <- bind_rows(results_list)
  
  return(results)
}

process_logistic <- function(X, Y, data, covariates = NULL, gender_label = 2,
                             min_n = 10) {
  
  library(nnet)
  library(dplyr)
  library(broom)
  
  # 1. 统一变量名
  original_names <- colnames(data)
  clean_names <- make.names(original_names, unique = TRUE)
  colnames(data) <- clean_names
  
  name_map <- setNames(clean_names, original_names)
  reverse_map <- setNames(original_names, clean_names)
  
  # 检查原始变量是否存在
  input_vars <- c(X, Y, covariates)
  not_found <- setdiff(input_vars, names(name_map))
  
  if (length(not_found) > 0) {
    stop(
      paste(
        "These variables are not in data:",
        paste(not_found, collapse = ", ")
      )
    )
  }
  
  # 同步更新变量名
  X_clean <- unname(name_map[X])
  Y_clean <- unname(name_map[Y])
  
  if (!is.null(covariates) && length(covariates) > 0) {
    cov_clean <- unname(name_map[covariates])
  } else {
    cov_clean <- character(0)
  }
  
  gender_text <- c("Female", "Male", "All")[gender_label + 1]
  
  # 2. 辅助函数：提取某个 predictor 对应的模型项
  get_predictor_terms <- function(fml, dat_sub, x) {
    
    mm <- model.matrix(fml, data = dat_sub)
    assign_vec <- attr(mm, "assign")
    term_labels <- attr(terms(fml), "term.labels")
    
    x_index <- which(term_labels == x)
    
    if (length(x_index) == 0) {
      return(character(0))
    }
    
    x_terms <- colnames(mm)[assign_vec == x_index]
    x_terms <- setdiff(x_terms, "(Intercept)")
    
    return(x_terms)
  }
  
  # 3. 辅助函数：生成可读的 predictor level
  get_predictor_level <- function(term, x, dat_sub) {
    
    if (is.factor(dat_sub[[x]]) || is.character(dat_sub[[x]])) {
      level <- sub(paste0("^", x), "", term)
      if (level == "") level <- NA_character_
      return(level)
    } else {
      return("per unit increase")
    }
  }
  
  # 4. 初始化结果
  results_list <- list()
  counter <- 1
  
  # 5. 主循环
  for (x in X_clean) {
    for (y in Y_clean) {
      
      vars_needed <- c(x, y, cov_clean)
      dat_sub <- data[, vars_needed, drop = FALSE]
      dat_sub <- na.omit(dat_sub)
      
      if (nrow(dat_sub) < min_n) next
      
      # 结局变量转为 factor
      dat_sub[[y]] <- factor(dat_sub[[y]])
      y_nlevel <- nlevels(dat_sub[[y]])
      
      # 如果结局只有一个水平，跳过
      if (y_nlevel < 2) next
      
      # 公式
      rhs_vars <- c(x, cov_clean)
      fml <- as.formula(
        paste(y, "~", paste(rhs_vars, collapse = " + "))
      )
      
      # 获取该 X 在模型中对应的所有系数项
      x_terms <- get_predictor_terms(fml, dat_sub, x)
      
      if (length(x_terms) == 0) next
      
      # 5.1 二分类 logistic
      if (y_nlevel == 2) {
        
        model <- try(
          glm(fml, data = dat_sub, family = binomial()),
          silent = TRUE
        )
        
        if (inherits(model, "try-error")) next
        
        coef_tab <- as.data.frame(summary(model)$coefficients)
        coef_tab$Term <- rownames(coef_tab)
        
        colnames(coef_tab)[1:4] <- c(
          "Beta", "SE", "Z_value", "P_value"
        )
        
        x_terms_exist <- intersect(x_terms, coef_tab$Term)
        
        if (length(x_terms_exist) == 0) next
        
        for (term in x_terms_exist) {
          
          tmp <- coef_tab %>%
            filter(Term == term) %>%
            mutate(
              Outcome = reverse_map[y],
              Outcome_ref = levels(dat_sub[[y]])[1],
              Outcome_level = levels(dat_sub[[y]])[2],
              Predictor = reverse_map[x],
              Predictor_term = term,
              Predictor_level = get_predictor_level(term, x, dat_sub),
              Beta = Beta,
              OR = exp(Beta),
              OR_95CI_low = exp(Beta - 1.96 * SE),
              OR_95CI_high = exp(Beta + 1.96 * SE),
              P_signif = ifelse(P_value < 0.05, "*", ""),
              Gender = gender_text,
              N = nrow(dat_sub),
              Model_type = "Binary logistic"
            ) %>%
            dplyr::select(
              Outcome, Outcome_ref, Outcome_level,
              Predictor, Predictor_term, Predictor_level,
              Beta, OR, OR_95CI_low, OR_95CI_high,
              P_value, P_signif, Gender, N, Model_type
            )
          
          results_list[[counter]] <- tmp
          counter <- counter + 1
        }
      }
      
      # 5.2 多分类 logistic
      if (y_nlevel >= 3) {
        
        model <- try(
          nnet::multinom(fml, data = dat_sub, trace = FALSE),
          silent = TRUE
        )
        
        if (inherits(model, "try-error")) next
        
        summ <- summary(model)
        coef_mat <- summ$coefficients
        se_mat <- summ$standard.errors
        
        # 防止只有一行时被转成 vector
        if (is.vector(coef_mat)) {
          coef_mat <- matrix(
            coef_mat,
            nrow = 1,
            dimnames = list(levels(dat_sub[[y]])[-1], names(coef_mat))
          )
        }
        
        if (is.vector(se_mat)) {
          se_mat <- matrix(
            se_mat,
            nrow = 1,
            dimnames = list(levels(dat_sub[[y]])[-1], names(se_mat))
          )
        }
        
        x_terms_exist <- intersect(x_terms, colnames(coef_mat))
        
        if (length(x_terms_exist) == 0) next
        
        for (ol in rownames(coef_mat)) {
          for (term in x_terms_exist) {
            
            beta <- coef_mat[ol, term]
            se <- se_mat[ol, term]
            zval <- beta / se
            pval <- 2 * (1 - pnorm(abs(zval)))
            
            tmp <- data.frame(
              Outcome = reverse_map[y],
              Outcome_ref = levels(dat_sub[[y]])[1],
              Outcome_level = ol,
              Predictor = reverse_map[x],
              Predictor_term = term,
              Predictor_level = get_predictor_level(term, x, dat_sub),
              Beta = beta,
              OR = exp(beta),
              OR_95CI_low = exp(beta - 1.96 * se),
              OR_95CI_high = exp(beta + 1.96 * se),
              P_value = pval,
              P_signif = ifelse(pval < 0.05, "*", ""),
              Gender = gender_text,
              N = nrow(dat_sub),
              Model_type = "Multinomial logistic",
              stringsAsFactors = FALSE
            )
            
            results_list[[counter]] <- tmp
            counter <- counter + 1
          }
        }
      }
    }
  }
  
  # 6. 合并结果
  if (length(results_list) == 0) {
    return(data.frame())
  }
  
  results <- bind_rows(results_list)
  
  return(results)
}
#******重复测量中介效应分析，个别菌群/蛋白质
run_microbiome_mediation_batch <- function(
    data,
    micro_sig,
    exposure_list,
    covariates,
    id_var = "ID",
    predictor_col = "Predictor",
    mediator_col = "Outcome",
    outcome_col = "Outcome_muscle",
    exposure_high = 2,
    exposure_low = 0,
    sims = 1000,
    output_file = NULL
){
  
  med_results_all <- list()
  
  for(exp_var in exposure_list){
    
    cat(
      "\n=============================\n",
      "Running:", exp_var,
      "\n=============================\n"
    )
    
    ## 提取 exposure 对应显著 mediator-outcome pair
    
    pair_use <- micro_sig %>%
      
      filter(
        .data[[predictor_col]] == exp_var
      ) %>%
      
      dplyr::select(
        mediator = all_of(mediator_col),
        outcome  = all_of(outcome_col)
      ) %>%
      
      distinct()
    
    ## 没有结果
    
    if(nrow(pair_use) == 0){
      
      cat(
        "No mediator-outcome pairs found for:",
        exp_var,
        "\n"
      )
      
      next
    }
    
    ## 循环 pair
    
    res_exp <- list()
    
    for(i in 1:nrow(pair_use)){
      
      mediator_i <- pair_use$mediator[i]
      outcome_i  <- pair_use$outcome[i]
      
      cat(
        "Mediator:",
        mediator_i,
        "| Outcome:",
        outcome_i,
        "\n"
      )
      
      tmp_res <- run_mediation_glmer_batch(
        
        data = data,
        
        id_var = id_var,
        
        exposures = exp_var,
        
        mediators = mediator_i,
        
        outcomes = outcome_i,
        
        covariates = covariates,
        
        exposure_high = exposure_high,
        
        exposure_low = exposure_low,
        
        sims = sims
      )
      
      tmp_res$Exposure <- exp_var
      
      res_exp[[i]] <- tmp_res
    }
    
    med_results_all[[exp_var]] <- bind_rows(res_exp)
  }
  
  ## 合并
  
  med_results_all_df <- bind_rows(
    med_results_all
  )
  
  ## 导出
  
  if(!is.null(output_file)){
    
    write.csv(
      med_results_all_df,
      output_file,
      row.names = FALSE
    )
  }
  
  return(
    med_results_all_df
  )
}

#******重复测量中介效应分析，指数
run_mediation_glmer_batch <- function(
    data,
    id_var = "ID",
    exposures,
    mediators,
    outcomes,
    covariates,
    exposure_high = 3,
    exposure_low  = 0,
    sims = 1000,
    seed = 908
) {
  
  library(dplyr)
  library(lme4)
  library(mediation)
  
  safe_var <- function(v) paste0("`", v, "`")
  
  results <- list()
  
  for(x in exposures){
    for(m in mediators){
      for(y in outcomes){
        
        df_med <- data %>%
          filter(.data[[x]] %in% c(exposure_low, exposure_high)) %>%
          mutate(Treat_bin = as.numeric(.data[[x]] == exposure_high))
        
        if(nrow(df_med) < 10 || length(unique(df_med$Treat_bin)) < 2) next
        
        cov_string <- paste(safe_var(covariates), collapse = " + ")
        
        formula_m <- as.formula(
          paste0(
            safe_var(m), " ~ Treat_bin",
            if(length(covariates) > 0) paste0(" + ", cov_string),
            " + (1 | ", safe_var(id_var), ")"
          )
        )
        
        formula_y <- as.formula(
          paste0(
            safe_var(y), " ~ Treat_bin + ", safe_var(m),
            if(length(covariates) > 0) paste0(" + ", cov_string),
            " + (1 | ", safe_var(id_var), ")"
          )
        )
        
        model.m <- tryCatch(
          suppressWarnings(glmer(formula_m, data = df_med, family = gaussian())),
          error = function(e) NULL
        )
        if(is.null(model.m)) next
        
        model.y <- tryCatch(
          suppressWarnings(glmer(formula_y, data = df_med, family = gaussian())),
          error = function(e) NULL
        )
        if(is.null(model.y)) next
        
        set.seed(seed)
        med_res <- tryCatch(
          mediate(model.m, model.y, treat = "Treat_bin", mediator = m, sims = sims),
          error = function(e) NULL
        )
        if(is.null(med_res)) next
        
        res <- data.frame(
          Exposure   = x,
          Mediator   = m,
          Outcome    = y,
          
          ACME       = med_res$d0,
          ACME_low   = med_res$d0.ci[1],
          ACME_high  = med_res$d0.ci[2],
          ACME_p     = med_res$d0.p,

          ADE        = med_res$z0,
          ADE_low    = med_res$z0.ci[1],
          ADE_high   = med_res$z0.ci[2],
          ADE_p      = med_res$z0.p,
          
          Total      = med_res$tau.coef,
          Total_low  = med_res$tau.ci[1],
          Total_high = med_res$tau.ci[2],
          Total_p    = med_res$tau.p,
          
          PropMediated     = med_res$n0,
          PropMed_low      = med_res$n0.ci[1],
          PropMed_high     = med_res$n0.ci[2],
          PropMediated_p   = med_res$n0.p
        )
        
        results[[length(results)+1]] <- res
      }
    }
  }
  
  final_res <- do.call(rbind, results)
  return(final_res)
}
#******中介效应分析-非重复测量
run_mediation_lm_once <- function(
    data,
    exposure,
    mediator,
    outcome,
    covariates = NULL,
    exposure_high = 2,
    exposure_low = 0,
    sims = 1000
){
  
  library(dplyr)
  library(mediation)
  
  ## ========= 统一的空结果函数 
  empty_result <- function(N = NA, Error = NA){
    
    data.frame(
      Exposure = exposure,
      Mediator = mediator,
      Outcome = outcome,
      
      ACME = NA,
      ACME_low = NA,
      ACME_high = NA,
      ACME_p = NA,
      
      ADE = NA,
      ADE_low = NA,
      ADE_high = NA,
      ADE_p = NA,
      
      Total = NA,
      Total_low = NA,
      Total_high = NA,
      Total_p = NA,
      
      PropMediated = NA,
      PropMediated_low = NA,
      PropMediated_high = NA,
      PropMediated_p = NA,
      
      N = N,
      Error = Error
    )
  }
  
  ## ========= 需要变量 
  vars_needed <- c(
    exposure,
    mediator,
    outcome,
    covariates
  )
  
  vars_needed <- unique(vars_needed)
  
  ## ========= complete-case + 只保留 exposure_low / exposure_high
  data_use <- data %>%
    filter(
      if_all(all_of(vars_needed), ~ !is.na(.))
    ) %>%
    filter(
      .data[[exposure]] %in% c(exposure_low, exposure_high)
    )
  
  ## ========= 样本量太小 
  if(nrow(data_use) < 10){
    
    return(
      empty_result(
        N = nrow(data_use),
        Error = "Sample size < 10"
      )
    )
  }
  
  ## ========= exposure 没有两个水平
  if(length(unique(data_use[[exposure]])) < 2){
    
    return(
      empty_result(
        N = nrow(data_use),
        Error = "Exposure has less than two levels"
      )
    )
  }
  
  ## ========= 重命名变量，避免特殊字符问题 
  model_data <- data_use %>%
    mutate(
      X = .data[[exposure]],
      M = .data[[mediator]],
      Y = .data[[outcome]]
    )
  
  if(!is.null(covariates) && length(covariates) > 0){
    
    for(j in seq_along(covariates)){
      model_data[[paste0("C", j)]] <- model_data[[covariates[j]]]
    }
    
    cov_terms <- paste0("C", seq_along(covariates))
    
  } else {
    
    cov_terms <- NULL
  }
  
  ## ========= 构建公式
  if(!is.null(cov_terms)){
    
    formula_m <- as.formula(
      paste(
        "M ~ X +",
        paste(cov_terms, collapse = " + ")
      )
    )
    
    formula_y <- as.formula(
      paste(
        "Y ~ X + M +",
        paste(cov_terms, collapse = " + ")
      )
    )
    
  } else {
    
    formula_m <- as.formula("M ~ X")
    formula_y <- as.formula("Y ~ X + M")
  }
  
  ## ========= 拟合模型和中介分析
  set.seed(920)
  res <- tryCatch({
    
    med_fit <- lm(
      formula_m,
      data = model_data
    )
    
    out_fit <- lm(
      formula_y,
      data = model_data
    )
    
    med_out <- mediate(
      model.m = med_fit,
      model.y = out_fit,
      treat = "X",
      mediator = "M",
      control.value = exposure_low,
      treat.value = exposure_high,
      sims = sims
    )
    
    s <- summary(med_out)
    
    ## ========= 结果整合为一行
    data.frame(
      Exposure = exposure,
      Mediator = mediator,
      Outcome = outcome,
      
      ACME = s$d.avg,
      ACME_low = s$d.avg.ci[1],
      ACME_high = s$d.avg.ci[2],
      ACME_p = s$d.avg.p,
      
      ADE = s$z.avg,
      ADE_low = s$z.avg.ci[1],
      ADE_high = s$z.avg.ci[2],
      ADE_p = s$z.avg.p,
      
      Total = s$tau.coef,
      Total_low = s$tau.ci[1],
      Total_high = s$tau.ci[2],
      Total_p = s$tau.p,
      
      PropMediated = s$n.avg,
      PropMediated_low = s$n.avg.ci[1],
      PropMediated_high = s$n.avg.ci[2],
      PropMediated_p = s$n.avg.p,
      
      N = nrow(model_data),
      Error = NA
    )
    
  }, error = function(e){
    
    empty_result(
      N = nrow(model_data),
      Error = e$message
    )
  })
  
  return(res)
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
#*******批量进行中介效应分析
run_microbiome_mediation_batch <- function(
    data,
    micro_sig,
    covariates,
    id_var = "ID",
    predictor_col = "Predictor",
    mediator_col = "Mediator",
    outcome_col = "Outcome_muscle",
    exposure_high = 2,
    exposure_low = 0,
    sims = 1000,
    output_file = NULL
){
  
  library(dplyr)
  
  med_results_all <- list()
  
  for(i in 1:nrow(micro_sig)){
    
    cat(
      "\n=============================\n",
      "Running row:", i,
      "\n=============================\n"
    )
    
    ## 提取当前组合
    
    exp_var <- micro_sig[[predictor_col]][i]
    
    mediator_i <- micro_sig[[mediator_col]][i]
    
    outcome_i <- micro_sig[[outcome_col]][i]
    
    cat(
      "Exposure :", exp_var, "\n",
      "Mediator :", mediator_i, "\n",
      "Outcome  :", outcome_i, "\n"
    )
    
    ## complete-case
    
    vars_needed <- c(
      exp_var,
      mediator_i,
      outcome_i,
      covariates
    )
    
    vars_needed <- unique(vars_needed)
    
    data_use <- data %>%
      filter(
        if_all(all_of(vars_needed), ~ !is.na(.))
      ) %>%
      filter(
        !!sym(exp_var) %in% c(exposure_low, exposure_high)
      )
    # data_use <- data %>%
    #   
    #   filter(
    #     if_all(
    #       all_of(vars_needed),
    #       ~ !is.na(.)
    #     )
    #   )
    cat(
      "Sample size:",
      nrow(data_use),
      "\n"
    )
    
    ## mediation
    
    tmp_res <- run_mediation_glmer_batch(
      
      data = data_use,
      
      id_var = id_var,
      
      exposures = exp_var,
      
      mediators = mediator_i,
      
      outcomes = outcome_i,
      
      covariates = covariates,
      
      exposure_high = exposure_high,
      
      exposure_low = exposure_low,
      
      sims = sims
    )
    
    ## 添加信息
    
    tmp_res$Exposure <- exp_var
    
    tmp_res$Mediator <- mediator_i
    
    tmp_res$Outcome <- outcome_i
    
    ## 保存
    
    med_results_all[[i]] <- tmp_res
  }
  
  ## 合并
  
  med_results_all_df <- bind_rows(
    med_results_all
  )
  
  ## 导出
  
  if(!is.null(output_file)){
    
    write.csv(
      med_results_all_df,
      output_file,
      row.names = FALSE
    )
  }
  
  return(
    med_results_all_df
  )
}

run_microbiome_mediation_batch_glm <- function(
    data,
    micro_sig,
    covariates = NULL,
    id_var = NULL,   # 保留这个参数只是为了兼容旧代码，非重复测量模型中不再使用
    predictor_col = "Predictor",
    mediator_col = "Mediator",
    outcome_col = "Outcome_muscle",
    exposure_high = 2,
    exposure_low = 0,
    sims = 1000,
    output_file = NULL
){
  
  library(dplyr)
  
  ## ========= 检查 micro_sig 中的列是否存在
  cols_needed_micro_sig <- c(
    predictor_col,
    mediator_col,
    outcome_col
  )
  
  missing_cols_micro_sig <- setdiff(
    cols_needed_micro_sig,
    colnames(micro_sig)
  )
  
  if(length(missing_cols_micro_sig) > 0){
    
    stop(
      "These columns are missing in micro_sig: ",
      paste(missing_cols_micro_sig, collapse = ", ")
    )
  }
  
  ## ========= 空结果函数：保证出错时也返回一行
  empty_result_batch <- function(
    exp_var,
    mediator_i,
    outcome_i,
    N = NA,
    Error = NA
  ){
    
    data.frame(
      Exposure = exp_var,
      Mediator = mediator_i,
      Outcome = outcome_i,
      
      ACME = NA,
      ACME_low = NA,
      ACME_high = NA,
      ACME_p = NA,
      
      ADE = NA,
      ADE_low = NA,
      ADE_high = NA,
      ADE_p = NA,
      
      Total = NA,
      Total_low = NA,
      Total_high = NA,
      Total_p = NA,
      
      PropMediated = NA,
      PropMediated_low = NA,
      PropMediated_high = NA,
      PropMediated_p = NA,
      
      N = N,
      Error = Error
    )
  }
  
  ## ========= 保存结果
  med_results_all <- list()
  
  for(i in 1:nrow(micro_sig)){
    
    cat(
      "\n=============================\n",
      "Running row:", i,
      "\n=============================\n"
    )
    
    ## ========= 提取当前组合
    exp_var <- micro_sig[[predictor_col]][i]
    mediator_i <- micro_sig[[mediator_col]][i]
    outcome_i <- micro_sig[[outcome_col]][i]
    
    cat(
      "Exposure :", exp_var, "\n",
      "Mediator :", mediator_i, "\n",
      "Outcome  :", outcome_i, "\n"
    )
    
    ## ========= 如果当前变量名缺失
    if(
      is.na(exp_var) ||
      is.na(mediator_i) ||
      is.na(outcome_i)
    ){
      
      tmp_res <- empty_result_batch(
        exp_var = exp_var,
        mediator_i = mediator_i,
        outcome_i = outcome_i,
        N = NA,
        Error = "Exposure, mediator, or outcome is NA"
      )
      
      tmp_res$Row <- i
      med_results_all[[i]] <- tmp_res
      
      next
    }
    
    ## ========= 检查变量是否在 data 中存在
    vars_needed <- c(
      exp_var,
      mediator_i,
      outcome_i,
      covariates
    )
    
    vars_needed <- unique(vars_needed)
    
    missing_vars_data <- setdiff(
      vars_needed,
      colnames(data)
    )
    
    if(length(missing_vars_data) > 0){
      
      cat(
        "Missing variables:",
        paste(missing_vars_data, collapse = ", "),
        "\n"
      )
      
      tmp_res <- empty_result_batch(
        exp_var = exp_var,
        mediator_i = mediator_i,
        outcome_i = outcome_i,
        N = NA,
        Error = paste(
          "Missing variables in data:",
          paste(missing_vars_data, collapse = ", ")
        )
      )
      
      tmp_res$Row <- i
      med_results_all[[i]] <- tmp_res
      
      next
    }
    
    ## ========= complete-case + 只保留 exposure_low / exposure_high
    data_use <- data %>%
      filter(
        if_all(all_of(vars_needed), ~ !is.na(.))
      ) %>%
      filter(
        .data[[exp_var]] %in% c(exposure_low, exposure_high)
      )
    
    cat(
      "Sample size:",
      nrow(data_use),
      "\n"
    )
    
    ## ========= 调用新的非重复测量中介分析函数
    tmp_res <- tryCatch({
      
      run_mediation_lm_once(
        data = data_use,
        exposure = exp_var,
        mediator = mediator_i,
        outcome = outcome_i,
        covariates = covariates,
        exposure_high = exposure_high,
        exposure_low = exposure_low,
        sims = sims
      )
      
    }, error = function(e){
      
      empty_result_batch(
        exp_var = exp_var,
        mediator_i = mediator_i,
        outcome_i = outcome_i,
        N = nrow(data_use),
        Error = e$message
      )
    })
    
    ## ========= 添加行号
    tmp_res$Row <- i
    
    ## ========= 保存
    med_results_all[[i]] <- tmp_res
  }
  
  ## ========= 合并结果
  med_results_all_df <- bind_rows(
    med_results_all
  ) %>%
    dplyr::select(
      Row,
      dplyr::everything()
    )
  
  ## ========= 导出
  if(!is.null(output_file)){
    
    write.csv(
      med_results_all_df,
      output_file,
      row.names = FALSE
    )
  }
  
  return(
    med_results_all_df
  )
}

#*******十折交叉
run_cv_stability <- function(data,
                             outcomes,
                             predictor,
                             covariates,
                             gender_label = 2,
                             k_fold = 10,
                             fdr_cutoff = 0.05,
                             freq_cutoff = 60,
                             seed = 920) {
  
  library(dplyr)
  library(caret)
  
  set.seed(seed)
  
  ## 创建 participant-level folds
  
  ids_all <- unique(data$ID)
  
  folds <- createFolds(
    ids_all,
    k = k_fold,
    list = TRUE,
    returnTrain = FALSE
  )
  
  ## 保存每折结果
  cv_results <- list()
  
  ## 循环 CV
  
  for (i in seq_len(k_fold)) {
    
    cat("Processing fold:", i, "/", k_fold, "\n")
    
    ## test/train IDs
    test_ids <- ids_all[folds[[i]]]
    
    train_ids <- setdiff(
      ids_all,
      test_ids
    )
    
    ## train data
    train_data <- data[
      data$ID %in% train_ids,
    ]
    
    ## LMER
    res_train <- process_lmer(
      X = predictor,
      Y = outcomes,
      data = train_data,
      covariates = covariates,
      gender_label = gender_label
    )
    
    ## 保留 Level==2
    res_train <- res_train[
      res_train$Level == 2,
    ]
    
    ## FDR
    res_train$P_FDR <- p.adjust(
      res_train$P_value,
      method = "fdr"
    )
    
    ## 显著结果
    sig_train <- res_train[
      res_train$P_FDR < fdr_cutoff,
    ]
    
    ## 保存
    if (nrow(sig_train) > 0) {
      
      cv_results[[i]] <- data.frame(
        fold = i,
        feature = sig_train$Outcome,
        direction = sign(
          sig_train$Estimate
        )
      )
    }
  }
  
  ## 合并结果
  
  cv_df <- bind_rows(cv_results)
  
  ## stability
  stability_df <- cv_df %>%
    group_by(feature, direction) %>%
    summarise(
      frequency = n() / k_fold * 100,
      .groups = "drop"
    ) %>%
    arrange(desc(frequency))
  
  ## 高频 feature
  stability_sig <- stability_df[
    stability_df$frequency >= freq_cutoff,
  ]
  
  ## 返回
  return(list(
    folds = folds,
    raw = cv_df,
    stability = stability_df,
    stability_sig = stability_sig
  ))
}

#********spls
run_lmer_pipeline <- function(
    
  predictor,
  
  outcomes,
  
  data,
  
  cov_simple = c("Age","Sex_F0","Phase_group"),
  
  cov_full = c(
    setdiff(Covariates_all_lmer, c("Protein_mean")),
    "Phase_group"
  ),
  
  level_keep = 2
){
  
  cat(
    "\n========================\n",
    "Running:",
    predictor,
    "\n========================\n"
  )
  
  #### Step1 简单调整
  
  res_simple <- process_lmer(
    
    c(predictor),
    
    outcomes,
    
    data,
    
    cov_simple,
    
    level_keep
  )
  
  ## 保留 level
  res_simple <- res_simple %>%
    
    filter(
      Level == level_keep
    )
  
  ## FDR
  res_simple <- res_simple %>%
    
    mutate(
      P_FDR = p.adjust(
        P_value,
        method = "fdr"
      )
    )
  
  ## 显著蛋白
  sig_outcomes <- res_simple %>%
    
    filter(
      P_FDR < 0.05
    ) %>%
    
    pull(Outcome) %>%
    
    unique()
  
  cat(
    "Significant proteins:",
    length(sig_outcomes),
    "\n"
  )
  
  #### 没有显著结果
  
  if(length(sig_outcomes) == 0){
    
    return(list(
      simple = res_simple,
      full = NULL,
      full_sig = NULL
    ))
  }
  
  #### Step2 全调整
  
  res_full <- process_lmer(
    
    c(predictor),
    
    sig_outcomes,
    
    data,
    
    cov_full,
    
    level_keep
  )
  
  ## 保留 level
  res_full <- res_full %>%
    
    filter(
      Level == level_keep
    )
  
  ## FDR
  res_full <- res_full %>%
    
    group_by(Predictor) %>%
    
    mutate(
      P_FDR = p.adjust(
        P_value,
        method = "fdr"
      )
    ) %>%
    
    ungroup()
  
  ## 最终显著
  res_full_sig <- res_full %>%
    
    filter(
      P_FDR < 0.05
    )
  
  cat(
    "Final significant proteins:",
    nrow(res_full_sig),
    "\n"
  )
  
  #### 输出
  
  return(list(
    
    simple = res_simple,
    
    full = res_full,
    
    full_sig = res_full_sig
  ))
}

#****************模型效能提高
safe_var <- function(v){
  paste0("`", v, "`")
}

make_lmer_formula <- function(y, fixed_terms, id_var = "ID"){
  
  fixed_terms <- fixed_terms[!is.na(fixed_terms)]
  fixed_terms <- fixed_terms[fixed_terms != ""]
  
  fixed_part <- paste(
    safe_var(fixed_terms),
    collapse = " + "
  )
  
  formula_text <- paste0(
    safe_var(y),
    " ~ ",
    fixed_part,
    " + (1|", safe_var(id_var), ")"
  )
  
  as.formula(formula_text)
}

calc_lmer_metrics <- function(model){
  
  y_obs <- model.response(model.frame(model))
  
  pred_fixed <- predict(
    model,
    re.form = NA
  )
  
  r2_res <- performance::r2_nakagawa(model)
  
  data.frame(
    Marginal_R2 = as.numeric(r2_res$R2_marginal),
    RMSE = sqrt(mean((y_obs - pred_fixed)^2, na.rm = TRUE))
  )
}

compare_lmer_predictors_by_pair <- function(
    pair_df,                         # 结果数据框，例如 Lmer_results7_all_ASM_sig
    data,
    covariates,
    predictor_col = "Predictor",
    outcome_col = "Outcome",
    predictor_type_col = "Variable_Type",
    id_var = "ID",
    continuous_vars = NULL,
    exposure_low = 0,
    exposure_high = 2,
    predictor_rename = NULL,
    outcome_rename = NULL,
    remove_duplicate_pairs = TRUE,
    output_file_raw = "Lmer_model_performance_raw.csv",
    output_file_format = "Lmer_model_performance_percent_format.csv"
){
  
  library(dplyr)
  library(lme4)
  library(performance)
  
  ## ========= 检查 pair_df 中是否有 Predictor 和 Outcome 列
  cols_needed_pair <- c(predictor_col, outcome_col)
  
  missing_pair_cols <- setdiff(
    cols_needed_pair,
    colnames(pair_df)
  )
  
  if(length(missing_pair_cols) > 0){
    stop(
      "These columns are missing in pair_df: ",
      paste(missing_pair_cols, collapse = ", ")
    )
  }
  
  ## ========= 提取 Predictor-Outcome 对应关系
  pair_use <- pair_df %>%
    dplyr::mutate(
      Pair_Row = dplyr::row_number(),
      Predictor_pair = .data[[predictor_col]],
      Outcome_pair = .data[[outcome_col]]
    )
  
  ## 如果有 Variable_Type，则保留
  if(predictor_type_col %in% colnames(pair_use)){
    pair_use <- pair_use %>%
      dplyr::mutate(
        Variable_Type_pair = .data[[predictor_type_col]]
      )
  } else {
    pair_use <- pair_use %>%
      dplyr::mutate(
        Variable_Type_pair = NA_character_
      )
  }
  
  ## 去掉重复的 Predictor-Outcome 组合
  if(remove_duplicate_pairs){
    pair_use <- pair_use %>%
      dplyr::distinct(
        Predictor_pair,
        Outcome_pair,
        .keep_all = TRUE
      )
  }
  
  results_all <- list()
  counter <- 1
  
  for(i in 1:nrow(pair_use)){
    
    x <- pair_use$Predictor_pair[i]
    y <- pair_use$Outcome_pair[i]
    variable_type_i <- pair_use$Variable_Type_pair[i]
    
    cat("\n============================\n")
    cat("Running pair:", i, "\n")
    cat("Outcome  :", y, "\n")
    cat("Predictor:", x, "\n")
    cat("============================\n")
    
    ## ========= 判断当前 Predictor 是连续变量还是分类变量
    if(!is.null(continuous_vars)){
      
      is_continuous_x <- x %in% continuous_vars
      
    } else if(!is.na(variable_type_i)){
      
      is_continuous_x <- variable_type_i %in% c(
        "Continuous",
        "continuous",
        "CONTINUOUS"
      )
      
    } else {
      
      ## 如果没有 continuous_vars，也没有 Variable_Type，就根据 data 中变量类型自动判断
      is_continuous_x <- is.numeric(data[[x]])
    }
    
    cat(
      "Predictor type:",
      ifelse(is_continuous_x, "Continuous", "Categorical"),
      "\n"
    )
    
    ## ========= 检查变量是否存在
    vars_needed <- unique(c(
      y,
      x,
      covariates,
      id_var
    ))
    
    missing_vars <- setdiff(
      vars_needed,
      colnames(data)
    )
    
    if(length(missing_vars) > 0){
      
      results_all[[counter]] <- data.frame(
        Pair_Row = pair_use$Pair_Row[i],
        Outcome = y,
        Predictor = x,
        Predictor_type = ifelse(is_continuous_x, "Continuous", "Categorical"),
        Contrast = ifelse(
          is_continuous_x,
          "Per unit increase",
          paste0(exposure_high, " vs ", exposure_low)
        ),
        N = NA,
        N_ID = NA,
        N_low = NA,
        N_high = NA,
        Base_Marginal_R2 = NA,
        Extended_Marginal_R2 = NA,
        Delta_Marginal_R2_percent = NA,
        Base_RMSE = NA,
        Extended_RMSE = NA,
        RMSE_reduction_percent = NA,
        LRT_P_value = NA,
        Error = paste(
          "Missing variables:",
          paste(missing_vars, collapse = ", ")
        )
      )
      
      counter <- counter + 1
      next
    }
    
    ## ========= 数据准备
    if(is_continuous_x){
      
      data_yx <- data %>%
        dplyr::filter(
          dplyr::if_all(
            dplyr::all_of(vars_needed),
            ~ !is.na(.)
          )
        ) %>%
        dplyr::mutate(
          X_compare = as.numeric(.data[[x]])
        )
      
      n_low <- NA
      n_high <- NA
      contrast_label <- "Per unit increase"
      predictor_type <- "Continuous"
      
    } else {
      
      data_yx <- data %>%
        dplyr::filter(
          dplyr::if_all(
            dplyr::all_of(vars_needed),
            ~ !is.na(.)
          )
        ) %>%
        dplyr::filter(
          as.character(.data[[x]]) %in%
            as.character(c(exposure_low, exposure_high))
        ) %>%
        dplyr::mutate(
          X_compare = factor(
            as.character(.data[[x]]),
            levels = as.character(c(exposure_low, exposure_high))
          )
        )
      
      n_low <- sum(data_yx$X_compare == as.character(exposure_low), na.rm = TRUE)
      n_high <- sum(data_yx$X_compare == as.character(exposure_high), na.rm = TRUE)
      contrast_label <- paste0(exposure_high, " vs ", exposure_low)
      predictor_type <- "Categorical"
    }
    
    cat("Sample size:", nrow(data_yx), "\n")
    cat("Number of IDs:", length(unique(data_yx[[id_var]])), "\n")
    
    if(!is_continuous_x){
      cat("N low:", n_low, "; N high:", n_high, "\n")
    }
    
    ## ========= 样本量或变量变异检查
    if(
      nrow(data_yx) < 20 ||
      dplyr::n_distinct(data_yx$X_compare) < 2 ||
      (!is_continuous_x && (n_low == 0 || n_high == 0))
    ){
      
      results_all[[counter]] <- data.frame(
        Pair_Row = pair_use$Pair_Row[i],
        Outcome = y,
        Predictor = x,
        Predictor_type = predictor_type,
        Contrast = contrast_label,
        N = nrow(data_yx),
        N_ID = length(unique(data_yx[[id_var]])),
        N_low = n_low,
        N_high = n_high,
        Base_Marginal_R2 = NA,
        Extended_Marginal_R2 = NA,
        Delta_Marginal_R2_percent = NA,
        Base_RMSE = NA,
        Extended_RMSE = NA,
        RMSE_reduction_percent = NA,
        LRT_P_value = NA,
        Error = "Sample size too small or predictor has insufficient variation"
      )
      
      counter <- counter + 1
      next
    }
    
    ## ========= 构建模型公式
    formula_base <- make_lmer_formula(
      y = y,
      fixed_terms = covariates,
      id_var = id_var
    )
    
    formula_extended <- make_lmer_formula(
      y = y,
      fixed_terms = c(covariates, "X_compare"),
      id_var = id_var
    )
    
    ## ========= 拟合模型并比较
    tmp_res <- tryCatch({
      
      model_base <- lme4::lmer(
        formula_base,
        data = data_yx,
        REML = FALSE,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 2e5)
        )
      )
      
      model_extended <- lme4::lmer(
        formula_extended,
        data = data_yx,
        REML = FALSE,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 2e5)
        )
      )
      
      perf_base <- calc_lmer_metrics(model_base)
      perf_extended <- calc_lmer_metrics(model_extended)
      
      lrt <- anova(
        model_base,
        model_extended
      )
      
      delta_r2_percent <- (
        perf_extended$Marginal_R2 -
          perf_base$Marginal_R2
      ) * 100
      
      rmse_reduction_percent <- (
        perf_base$RMSE -
          perf_extended$RMSE
      ) / perf_base$RMSE * 100
      
      data.frame(
        Pair_Row = pair_use$Pair_Row[i],
        Outcome = y,
        Predictor = x,
        Predictor_type = predictor_type,
        Contrast = contrast_label,
        N = nrow(data_yx),
        N_ID = length(unique(data_yx[[id_var]])),
        N_low = n_low,
        N_high = n_high,
        Base_Marginal_R2 = perf_base$Marginal_R2,
        Extended_Marginal_R2 = perf_extended$Marginal_R2,
        Delta_Marginal_R2_percent = delta_r2_percent,
        Base_RMSE = perf_base$RMSE,
        Extended_RMSE = perf_extended$RMSE,
        RMSE_reduction_percent = rmse_reduction_percent,
        LRT_P_value = lrt[2, "Pr(>Chisq)"],
        Error = NA
      )
      
    }, error = function(e){
      
      data.frame(
        Pair_Row = pair_use$Pair_Row[i],
        Outcome = y,
        Predictor = x,
        Predictor_type = predictor_type,
        Contrast = contrast_label,
        N = nrow(data_yx),
        N_ID = length(unique(data_yx[[id_var]])),
        N_low = n_low,
        N_high = n_high,
        Base_Marginal_R2 = NA,
        Extended_Marginal_R2 = NA,
        Delta_Marginal_R2_percent = NA,
        Base_RMSE = NA,
        Extended_RMSE = NA,
        RMSE_reduction_percent = NA,
        LRT_P_value = NA,
        Error = e$message
      )
    })
    
    results_all[[counter]] <- tmp_res
    counter <- counter + 1
  }
  
  ## ========= 合并结果：数值版
  results_raw <- dplyr::bind_rows(results_all) %>%
    dplyr::group_by(Outcome) %>%
    dplyr::mutate(
      Rank_by_Delta_R2 = rank(
        -Delta_Marginal_R2_percent,
        ties.method = "min",
        na.last = "keep"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Base_Marginal_R2 = round(Base_Marginal_R2, 5),
      Extended_Marginal_R2 = round(Extended_Marginal_R2, 5),
      Delta_Marginal_R2_percent = round(Delta_Marginal_R2_percent, 3),
      Base_RMSE = round(Base_RMSE, 5),
      Extended_RMSE = round(Extended_RMSE, 5),
      RMSE_reduction_percent = round(RMSE_reduction_percent, 3),
      LRT_P_value = signif(LRT_P_value, 3)
    )
  
  ## ========= 格式化版
  results_format <- results_raw
  
  if(!is.null(predictor_rename)){
    results_format <- results_format %>%
      dplyr::mutate(
        Predictor_raw = Predictor,
        Predictor = dplyr::recode(
          Predictor,
          !!!predictor_rename
        )
      )
  }
  
  if(!is.null(outcome_rename)){
    results_format <- results_format %>%
      dplyr::mutate(
        Outcome_raw = Outcome,
        Outcome = dplyr::recode(
          Outcome,
          !!!outcome_rename
        )
      )
  }
  
  results_format <- results_format %>%
    dplyr::mutate(
      Delta_Marginal_R2_percent = ifelse(
        is.na(Delta_Marginal_R2_percent),
        NA,
        paste0(sprintf("%.3f", Delta_Marginal_R2_percent), "%")
      ),
      RMSE_reduction_percent = ifelse(
        is.na(RMSE_reduction_percent),
        NA,
        paste0(sprintf("%.3f", RMSE_reduction_percent), "%")
      ),
      LRT_P_value = dplyr::case_when(
        is.na(LRT_P_value) ~ NA_character_,
        LRT_P_value < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", LRT_P_value)
      )
    )
  
  ## ========= 导出
  if(!is.null(output_file_raw)){
    write.csv(
      results_raw,
      output_file_raw,
      row.names = FALSE
    )
  }
  
  if(!is.null(output_file_format)){
    write.csv(
      results_format,
      output_file_format,
      row.names = FALSE
    )
  }
  
  return(
    list(
      raw = results_raw,
      format = results_format
    )
  )
}

#*********lasso
lasso_auto <- function(data, outcome, predictors, nfolds = 10, seed = 920) {
  
  library(glmnet)
  
  # 提取 x 和 y
  x <- as.matrix(data[, predictors])
  y <- data[[outcome]]
  
  #  自动识别结局类型
  
  if (is.factor(y)) {
    
    n_class <- nlevels(y)
    
    if (n_class == 2) {
      family_type <- "binomial"
      cat("结局类型：二分类（binomial）\n")
      
    } else if (n_class > 2) {
      family_type <- "multinomial"
      cat("结局类型：多分类（multinomial）\n")
    }
    
  } else if (length(unique(y)) == 2) {
    
    family_type <- "binomial"
    y <- as.factor(y)
    cat("结局类型：二分类（binomial）\n")
    
  } else {
    
    family_type <- "gaussian"
    y <- as.numeric(y)
    cat("结局类型：连续型（gaussian）\n")
  }
  
  #  LASSO 路径图
  
  lasso_model <- glmnet(x, y, family = family_type, alpha = 1)
  plot(lasso_model, xvar = "lambda", label = FALSE)
  
  #  交叉验证
  
  set.seed(seed)
  cv_model <- cv.glmnet(
    x, y,
    family = family_type,
    alpha = 1,
    nfolds = nfolds
  )
  
  lambda_min <- cv_model$lambda.min
  lambda_1se <- cv_model$lambda.1se
  
  #  提取系数函数（兼容 multinomial）
  
  extract_coef <- function(model) {
    
    coef_raw <- coef(model)
    
    # multinomial 情况
    if (family_type == "multinomial") {
      
      results <- data.frame()
      
      for (cls in names(coef_raw)) {
        
        mat <- as.matrix(coef_raw[[cls]])
        df <- data.frame(
          Var = rownames(mat),
          Class = cls,
          Coef = mat[,1]
        )
        
        df <- df[df$Var != "(Intercept)" & df$Coef != 0, ]
        results <- rbind(results, df)
      }
      
      return(results)
      
    } else {
      
      mat <- as.matrix(coef_raw)
      df <- data.frame(
        Var = rownames(mat),
        Coef = mat[,1]
      )
      
      df <- df[df$Var != "(Intercept)" & df$Coef != 0, ]
      return(df)
    }
  }
  
  #  最优模型
  
  model_min <- glmnet(x, y, family = family_type,
                      alpha = 1, lambda = lambda_min)
  
  model_1se <- glmnet(x, y, family = family_type,
                      alpha = 1, lambda = lambda_1se)
  
  coef_min <- extract_coef(model_min)
  coef_1se <- extract_coef(model_1se)
  
  #  输出
  
  return(list(
    family = family_type,
    lambda_min = lambda_min,
    lambda_1se = lambda_1se,
    coef_min = coef_min,
    coef_1se = coef_1se,
    cv_model = cv_model,
    model_min = model_min,
    model_1se = model_1se
  ))
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
#******组间指数比较
plot_index_bar <- function(data,
                           lmer_result,
                           index_var = "PI_Tea_z",
                           group_var = "tea_freq_group",
                           time_var = "Times",
                           y_title = "Tea-Related Proteomic Index (Z-score) ",
                           group_levels = c(0, 2),
                           group_labels = c(
                             "Non-drinker",
                             "≥ 7 times/week"
                           ),
                           palette = c(
                             "#D8DEE9",
                             "#C06C84"
                           )) {
  
  library(dplyr)
  library(ggplot2)
  
  ## 数据整理
  
  plot_df <- data %>%
    
    filter(
      .data[[group_var]] %in% group_levels
    ) %>%
    
    mutate(
      Group = factor(
        .data[[group_var]],
        levels = group_levels,
        labels = group_labels
      )
    ) %>%
    
    group_by(
      .data[[time_var]],
      Group
    ) %>%
    
    summarise(
      mean_value = mean(
        .data[[index_var]],
        na.rm = TRUE
      ),
      
      se = sd(
        .data[[index_var]],
        na.rm = TRUE
      ) /
        sqrt(n()),
      
      .groups = "drop"
    )
  
  ## 提取模型结果（2 vs 0）
  
  model_res <- lmer_result %>%
    
    filter(Level == 2)
  
  beta_text <- paste0(
    "β = ",
    sprintf("%.3f", model_res$Estimate),
    " (",
    sprintf("%.3f", model_res$CI_low),
    " to ",
    sprintf("%.3f", model_res$CI_high),
    ")"
  )
  
  ## Y轴上限
  
  ymax <- max(
    plot_df$mean_value +
      plot_df$se,
    na.rm = TRUE
  )
  
  ## 作图
  
  p <- ggplot(
    plot_df,
    aes(
      x = .data[[time_var]],
      y = mean_value,
      fill = Group
    )
  ) +
    
    ## 柱状图
    geom_col(
      position = position_dodge(
        width = 0.72
      ),
      width = 0.62,
      color = "white",
      linewidth = 0.7
    ) +
    
    ## 误差线
    geom_errorbar(
      aes(
        ymin = mean_value - se,
        ymax = mean_value + se
      ),
      position = position_dodge(
        width = 0.72
      ),
      width = 0.16,
      linewidth = 0.65,
      color = "#4B5563"
    ) +
    
    ## 顶部显著性文字
    annotate(
      "label",
      x = 1.8,
      y = ymax + 0.3,
      label = beta_text,
      size = 4.2,
      label.size = 0,
      fill = "#F8FAFC",
      color = "#111827",
      fontface = "bold",
      lineheight = 1.15
    ) +
    
    ## 配色
    scale_fill_manual(
      values = palette
    ) +
    
    ## 标签
    labs(
      x = NULL,
      y = y_title,
      fill = NULL
    ) +
    
    ## Y轴留白
    expand_limits(
      y = ymax + 0.32
    ) +
    
    ## 高级主题
    theme_classic(
      base_size = 15
    ) +
    
    theme(
      
      axis.title.y = element_text(
        face = "bold",
        size = 15,
        color = "#111827"
      ),
      
      axis.text = element_text(
        size = 13,
        color = "#374151"
      ),
      
      legend.position = "top",
      
      legend.text = element_text(
        size = 12,
        face = "bold"
      ),
      
      legend.key.width = unit(
        1.3,
        "cm"
      ),
      
      panel.border = element_blank(),
      
      axis.line = element_line(
        linewidth = 0.7,
        color = "#374151"
      ),
      
      plot.margin = ggplot2::margin(
        18,
        18,
        18,
        18
      )
    )
  
  return(p)
}
#####菌/功能/蛋白质&骨骼肌
plot_forest_facet <- function(
    data,
    estimate_col = "Estimate",
    ci_low_col = "CI_low",
    ci_high_col = "CI_high",
    pathway_col = "Pathway",
    outcome_col = "Outcome",
    
    ## P值和FDR列名
    p_col = "P_value",
    fdr_col = "P_FDR",
    
    ## 显著性规则
    sig_rule = c("p_fdr", "p_only", "none"),
    
    ## 阈值参数
    p_level = 0.05,
    fdr_level = 0.10,
    
    ncol = 3,
    point_size = 5.8,
    ci_width = 1.4,
    base_size = 15,
    low_color = "#4C78A8",
    mid_color = "white",
    high_color = "#D65F5F",
    y_text_face = "italic"
) {
  
  library(ggplot2)
  library(dplyr)
  library(rlang)
  library(grid)
  
  sig_rule <- match.arg(sig_rule)
  
  ## 1. 非标准求值
  estimate_sym <- rlang::sym(estimate_col)
  ci_low_sym   <- rlang::sym(ci_low_col)
  ci_high_sym  <- rlang::sym(ci_high_col)
  pathway_sym  <- rlang::sym(pathway_col)
  outcome_sym  <- rlang::sym(outcome_col)
  p_sym        <- rlang::sym(p_col)
  
  ## 2. 根据显著性规则生成星号
  if (sig_rule == "p_fdr") {
    
    if (!(p_col %in% colnames(data))) {
      stop("p_col 不在 data 中，请检查 p_col 参数。")
    }
    
    if (!(fdr_col %in% colnames(data))) {
      stop("fdr_col 不在 data 中，请检查 fdr_col 参数。")
    }
    
    fdr_sym <- rlang::sym(fdr_col)
    
    plot_df <- data %>%
      dplyr::mutate(
        sig = dplyr::case_when(
          !!p_sym < p_level & !!fdr_sym < fdr_level ~ "*",
          TRUE ~ ""
        )
      )
    
  } else if (sig_rule == "p_only") {
    
    if (!(p_col %in% colnames(data))) {
      stop("p_col 不在 data 中，请检查 p_col 参数。")
    }
    
    plot_df <- data %>%
      dplyr::mutate(
        sig = dplyr::case_when(
          !!p_sym < p_level ~ "*",
          TRUE ~ ""
        )
      )
    
  } else {
    
    plot_df <- data %>%
      dplyr::mutate(
        sig = ""
      )
  }
  
  ## 3. 绘图
  p <- ggplot(
    plot_df,
    aes(
      x = !!estimate_sym,
      y = !!pathway_sym
    )
  ) +
    
    geom_segment(
      aes(
        x = !!ci_low_sym,
        xend = !!ci_high_sym,
        yend = !!pathway_sym,
        color = !!estimate_sym
      ),
      linewidth = ci_width,
      alpha = 0.9
    ) +
    
    geom_point(
      aes(
        fill = !!estimate_sym
      ),
      shape = 21,
      size = point_size,
      color = "black",
      stroke = 0.25
    ) +
    
    geom_text(
      aes(label = sig),
      color = "black",
      size = 4.5,
      fontface = "bold"
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey70",
      linewidth = 0.5
    ) +
    
    facet_wrap(
      vars(!!outcome_sym),
      ncol = ncol
    ) +
    
    scale_fill_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = 0,
      name = "Estimate"
    ) +
    
    scale_color_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = 0
    ) +
    
    guides(
      color = "none"
    ) +
    
    labs(
      x = "β (95% CI)",
      y = NULL
    ) +
    
    theme_minimal(
      base_size = base_size
    ) +
    
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.text.y = element_text(
        size = 12,
        color = "black",
        face = y_text_face
      ),
      
      axis.text.x = element_text(
        size = 13,
        color = "black"
      ),
      
      axis.title = element_text(
        face = "bold",
        size = 16
      ),
      
      strip.text = element_text(
        face = "bold",
        size = 15
      ),
      
      strip.background = element_rect(
        fill = "#F3F3F3",
        color = NA
      ),
      
      legend.position = "right",
      
      panel.spacing = unit(1.5, "lines")
    )
  
  return(p)
}
####中介效应可视化
plot_mediation_alluvial <- function(
    data,
    mediator_col = "Mediator",
    title = "Mediation Effects",
    fill_low = "#CFE8F3",
    fill_high = "#0072B2",
    stratum_fill = "grey92",
    stratum_color = "grey50"
){
  
  p <- ggplot(
    data,
    aes(
      axis1 = Exposure,
      axis2 = .data[[mediator_col]],
      axis3 = Outcome,
      y = ACME
    )
  ) +
    
    geom_alluvium(
      aes(fill = PropMediated),
      width = 1/12,
      knot.pos = 0.4
    ) +
    
    geom_stratum(
      width = 0.35,
      fill = stratum_fill,
      color = stratum_color
    ) +
    
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 4,
      family = "serif"
    ) +
    
    scale_fill_gradient(
      low = fill_low,
      high = fill_high,
      name = "Prop. mediated"
    ) +
    
    labs(
      title = title,
      y = NULL,
      caption = NULL
    ) +
    
    theme_minimal(base_size = 12) +
    
    theme(
      
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      
      panel.grid = element_blank(),
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14,
        margin = ggplot2::margin(b = 0)
      ),
      
      legend.position = "right",
      
      legend.title = element_text(
        size = 11,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 10
      ),
      
      plot.margin = ggplot2::margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    )
  
  return(p)
}

# "Outlined flows indicate ACME FDR < 0.05"
####森林图可视化
plot_forest_simple <- function(data,
                               x,
                               y,
                               p_col,
                               facet = "Outcome",
                               ci_low = "CI_low",
                               ci_high = "CI_high",
                               x_lab = "β (95% CI)",
                               y_lab = NULL,
                               output_file = NULL,
                               width = 8,
                               height = 5) {
  
  library(ggplot2)
  library(dplyr)
  
  data_plot <- data %>%
    mutate(
      .x = .data[[x]],
      .y = .data[[y]],
      .ci_low = .data[[ci_low]],
      .ci_high = .data[[ci_high]],
      .facet = .data[[facet]],
      .p = .data[[p_col]],
      .sig = if_else(.p < 0.05, "Significant", "Not Significant"),
      .sig = factor(.sig, levels = c("Not Significant", "Significant"))
    )
  
  p <- ggplot(
    data_plot,
    aes(
      x = .x,
      y = .y,
      color = .sig,
      shape = .sig
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "#888888",
      linewidth = 0.6
    ) +
    geom_errorbarh(
      aes(
        xmin = .ci_low,
        xmax = .ci_high
      ),
      height = 0.2,
      linewidth = 0.8
    ) +
    geom_point(size = 3.5) +
    facet_wrap(
      ~ .facet,
      scales = "free_x"
    ) +
    scale_color_manual(
      values = c(
        "Significant" = "#D95F02",
        "Not Significant" = "#7570B3"
      ),
      name = "Significance"
    ) +
    scale_shape_manual(
      values = c(
        "Significant" = 15,
        "Not Significant" = 16
      ),
      name = "Significance"
    ) +
    labs(
      x = x_lab,
      y = y_lab
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.background = element_rect(
        fill = "#F2F4F4",
        color = NA
      ),
      strip.text = element_text(
        face = "bold",
        size = 12,
        color = "#2C3E50"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.spacing = unit(2, "lines"),
      axis.title = element_text(
        face = "bold",
        size = 12
      ),
      axis.text.y = element_text(
        face = "bold",
        size = 11,
        color = "#333333"
      ),
      axis.text.x = element_text(size = 10),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      legend.background = element_blank()
    )
  
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = p,
      width = width,
      height = height,
      dpi = 600
    )
  }
  
  return(p)
}
####变量名转换####
predictor_rename <- c(
  "tea_freq_group" = "Tea consumption",
  "serum_I_catechin_F0_T" = "Catechin",
  "serum_I_EGC_F0_T" = "Epigallocatechin",
  "serum_I_EGCG_F0_T" = "Epigallocatechin gallate",
  "serum_I_epicatechin_F0_T" = "Epicatechin",
  "serum_I_ECG_F0_T" = "Epicatechin gallate"
)

outcome_rename <- c(
  "ASM_z"        = "ASM",
  "ASMI_z"       = "ASMI",
  "grip_max_z"   = "Handgrip strength",
  "gait_speed_z" = "Walking speed"
)


get_level_label <- function(predictor, level_value) {
  case_when(
    predictor == "tea_freq_group" ~ "≥ 7 times/week vs. Non-drinker",
    predictor %in% c("serum_I_epicatechin_F0_T", "serum_I_EGC_F0_T") ~ "High vs. Undetectable",
    predictor %in% c("serum_I_catechin_F0_T", "serum_I_EGCG_F0_T", "serum_I_ECG_F0_T") ~ "Tertile 3 vs. Tertile 1",
    TRUE ~ as.character(level_value)
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
Cov_F0 <- read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_基本信息_F01234_20230214_VF.xlsx",sheet="F0")
# 计算准确的F0年龄
F0_surveyTimes <- read_sav("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/5010调查时间及随访间隔（重新编制）__6-1-2017_VF.sav")
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
colnames(Cov_F0)
# 只保留想要的变量
Cov_F0 <- Cov_F0[,c("ID","Age_F0","出生日期_F0","性别_F0","本人教育程度分类_3_F0","家庭人均月入_4_F0","是否吸烟_F0","是否饮酒_F0","是否喝茶_F0","过去一年是否钙片_F0","过去一年是否复合维生素_F0",
                    "脑卒中_F0","心脏病_F0","癫痫_F0","帕金森_F0","老年痴呆_F0","糖尿病_F0","癌症_F0","骨折_F0","雌激素_F0","绝经年限_F0"
)]
colnames(Cov_F0) <- dplyr::recode(colnames(Cov_F0),
                                  "出生日期_F0" = "Birthday_F0",
                                  "性别_F0" = "Sex_F0",
                                  "本人教育程度分类_3_F0" = "Education_F0",
                                  "家庭人均月入_4_F0" = "Income_F0",
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
Cov_F1 <- read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_基本信息_F01234_20230214_VF.xlsx",sheet="F1")
Cov_F1$ID <- sub("^F1", "", Cov_F1$CODE_F1)
Cov_F1 <- Cov_F1 %>%
  filter(grepl("^NL4", ID))

Cov_F1 <- Cov_F1[,c("ID","性别_F1","家庭人均收入_4_F1","是否吸烟_F1","是否饮酒_F1","是否喝茶_F1","过去一年是否钙片_F1","过去一年是否复合维生素_F1",
                    "脑卒中_F1","心脏病_F1","癫痫_F1","帕金森病_F1","老年痴呆_F1","糖尿病_F1","癌症_F1","骨折_F1","绝经年限_F1", "雌激素_F1"
)]
Birthday_NL4 <- read_sav("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/其它/F1NL教育20210806_VF.sav") 
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
                                  "家庭人均收入_4_F1" = "Income_F0",
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
Follow_time <- read.csv("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/FollowData_VF.csv")
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
colnames(Cov_F0)
Cov_F0_all <- rbind(Cov_F0, Cov_F1)
#雌激素为NA的是男性，所以除了女性有1外，其它直接赋值为0即可
Cov_F0_all$Estrogen_F0 <- ifelse(Cov_F0_all$`雌激素_F0` == 1,1,0)
table(Cov_F0_all$`雌激素_F0`)
table(Cov_F0_all$Estrogen)
#疾病只要之一为1，则新变量Disease就为1
Cov_F0_all$Disease_F0 <- ifelse(Cov_F0_all$`脑卒中_F0` == 1 | Cov_F0_all$`心脏病_F0` == 1 | Cov_F0_all$`癫痫_F0` == 1 | Cov_F0_all$`帕金森_F0` == 1 | Cov_F0_all$`老年痴呆_F0` == 1 | Cov_F0_all$`糖尿病_F0` == 1 | Cov_F0_all$`癌症_F0` == 1, 1, 0)
colnames(Cov_F0_all)
Cov_F0_all <- Cov_F0_all[,c("ID","Age_F0","Birthday_F0","Sex_F0","Education_F0","Smoke_F0", "Alcohol_F0","Tea_F0","Calcium_F0","Income_F0",
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
Cov2 <- read_excel('D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_身高体重MET_VF.xlsx', 
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
Nutrition <- read_excel('D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_F0123膳食营养素数据_20230206_VF.xlsx',
                        sheet="GNHS_膳食数据_20230206")
colnames(Nutrition)
a <- Nutrition[Nutrition$编号 == "F3NL1221",]
Nutrition$ID <- substr(Nutrition$`编号`, 3, nchar(Nutrition$`编号`))
Nutrition <- Nutrition[,c("ID","Followup","能量摄入","蛋白质","可溶性纤维")]
colnames(Nutrition)
# 找出重复观测
dup_data <- Nutrition %>%
  group_by(ID, Followup) %>%
  filter(n() > 1) %>%
  arrange(ID, Followup)
# 保留重复观测的第一条
Nutrition <- Nutrition %>%
  distinct(ID, Followup, .keep_all = TRUE)
####宽数据
Nutrition_wide <- Nutrition %>%
  pivot_wider(
    names_from = Followup,
    values_from = c(能量摄入, 蛋白质, 可溶性纤维),
    names_glue = "F{Followup}_{.value}"
  )
####计算均值
Nutrition_wide <- Nutrition_wide %>%
  mutate(across(where(is.list), ~ as.character(.))) %>%  # 先转换 list 为字符
  mutate(across(-ID, as.numeric))  # 再转换数值

Nutrition_wide$Energy_mean <- rowMeans(Nutrition_wide[, c("F0_能量摄入","F1_能量摄入","F2_能量摄入","F3_能量摄入")], na.rm = TRUE)
Nutrition_wide$Protein_mean <- rowMeans(Nutrition_wide[, c("F0_蛋白质","F1_蛋白质","F2_蛋白质","F3_蛋白质")], na.rm = TRUE)
Nutrition_wide$Fiber_soluble_mean <- rowMeans(Nutrition_wide[, c("F0_可溶性纤维","F1_可溶性纤维","F2_可溶性纤维","F3_可溶性纤维")], na.rm = TRUE)
####鱼油补充剂####
Oil1 <- read_excel('D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_基本信息_F01234_20230214_VF.xlsx',
                   sheet=1)
Oil1 <- Oil1[,c("ID","是否使用深海鱼油胶囊_F2","过去一年是否服用深海鱼油_F3","否服用过深海鱼油胶囊_F1")]

Oil2 <- read_excel('D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_基本信息_F01234_20230214_VF.xlsx',
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
####合并所有混杂####
data_list <- list(Cov_F0_all,Cov2_wide[,c("ID", "Met_mean")],Nutrition_wide[,c("ID","Energy_mean","Protein_mean","Fiber_soluble_mean")], Oil_complete_wide[,c("ID","Oil_sup_follow")])
Cov_F0_final <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID")
Cov_F0_final <- Cov_F0_final[complete.cases(Cov_F0_final),]
####******膳食指数******####
####食物组
Food <- read_excel('D:/OneDrive/Papers/2.3 荷兰合作/database/20240327_GNHS膳食摄入数据.xlsx', 
                   sheet=1)
#酒精很多为NA值
Food[is.na(Food)] <- 0
Food <- as.data.frame(lapply(Food, function(x) gsub("NA", "0", x)))
colnames(Food)

exclude_vars <- c("ID", "编号", "Followup")
Food <- Food %>%
  mutate(
    across(
      -all_of(exclude_vars),
      ~ suppressWarnings(as.numeric(as.character(.)))
    )
  )
# 查看变量情况
Food_summary <- check_df_summary(Food)

####营养素
Nutrients <- read_excel('D:/OneDrive/Papers/2.3 荷兰合作/database/20240530_GNHS_膳食数据.xlsx', 
                        sheet=1)
colnames(Nutrients)

exclude_vars <- c("ID", "编号", "Followup")
Nutrients <- Nutrients %>%
  mutate(
    across(
      -all_of(exclude_vars),
      ~ suppressWarnings(as.numeric(as.character(.)))
    )
  )
# 查看变量情况
Nutrients_summary <- check_df_summary(Nutrients)
colnames(Nutrients)
####合并食物和营养素
Food_Nutrients <- merge(dplyr::select(Food,-c("ID","Followup")),dplyr::select(Nutrients,-c("ID","Followup")),by="编号")
# 保留F0NL开头以及F1NL4开头的观测
Food_Nutrients <- Food_Nutrients %>%
  filter(
    grepl("^F0NL|^F1NL4", 编号)
  ) %>%
  mutate(
    ID = substring(编号, 3)
  )
# 合并性别变量
Food_Nutrients <- merge(Cov_F0_final[,c("ID","Sex_F0")], Food_Nutrients, by = "ID")
####验证食物组的计算总量####
# total_grains <- Food_Nutrients[,c("编号", "能量摄入", "米饭生重","稀饭生重","面条生重","全麦面包生重","馒头生重",
#                           "白面面包生重","肉包生重","油条生重","萝卜糕生重","蛋糕生重","水饺生重","饼干生重","粮谷类生重")]
# 
# total_grains <- total_grains %>%
#   mutate(
#     验证粮谷类生重 = rowSums(
#       select(., -编号, -粮谷类生重),
#       na.rm = TRUE
#     )
#   )

# Meat_egg <- Food_Nutrients[,c("编号", "半肥瘦猪肉生重", "瘦猪肉生重", "猪手生重" , "牛羊肉生重" ,"猪肚生重","肝肾脑生重" ,"腊肉等生重","鸡鸭带皮生重","鸡鸭去皮生重","鸡爪生重","红肉类生重" ,"禽肉类生重","蛋生重")]
# Meat_egg$Meat_egg <- Meat_egg$红肉类生重 + Meat_egg$禽肉类生重 + Meat_egg$蛋生重
# 
# Meat_egg_check <- Meat_egg %>%
#   mutate(
#     
#     红肉类_计算值 =
#       coalesce(半肥瘦猪肉生重, 0) +
#       coalesce(瘦猪肉生重, 0) +
#       coalesce(猪手生重, 0) +
#       coalesce(牛羊肉生重, 0) +
#       coalesce(猪肚生重, 0) +
#       coalesce(肝肾脑生重, 0) +
#       coalesce(腊肉等生重, 0),
#     
#     禽肉类_计算值 =
#       coalesce(鸡鸭带皮生重, 0) +
#       coalesce(鸡鸭去皮生重, 0) +
#       coalesce(鸡爪生重, 0)
#   )


# Sugar <- Food_Nutrients[,c("编号", "Total_sugar", "Sucrose","Dextrose", "Fructose", "Lactose", "Maltose" , "Galactose","Starch")]
# Sugar <- Sugar %>%
#   mutate(
#     验证Sugar = rowSums(
#       select(., -编号, -Total_sugar),
#       na.rm = TRUE
#     )
#   )
####进行赋分####
Food_Nutrients <- Food_Nutrients %>%
  mutate(
    # 总粮谷类
    Total_grains = 粮谷类生重 / 能量摄入 * 1000,
    Total_grains_score = pmin(Total_grains / (2.5 * (250 / 5)) * 5, 5),
    
    # 全谷/杂豆
    Whole_grains_bean = (全麦面包生重 + 绿豆等生重 + 新鲜粟米生重)/ 能量摄入 * 1000,
    Whole_grains_bean_score = pmin(Whole_grains_bean / (0.6 * (100 / 2)) * 5, 5),
    
    # 薯类
    Tuber = 淀粉类生重/ 能量摄入 * 1000,
    Tuber_score =pmin((Tuber / 能量摄入 * 1000) / (0.3 * (75 / 0.8)) * 5,5),
    
    # 总蔬菜
    Total_vegetables = 蔬菜类生重 / 能量摄入 * 1000,
    Total_vegetables_score = pmin(Total_vegetables / (1.9 * (450 / 4.5)) * 5, 5),
    
    # 深色蔬菜
    Dark_vegetables = (豆角等生重 + 菜心等生重+  菠菜等生重 + 其它深绿色叶菜生重 + 红萝卜生重 + 西红柿生重 + 青椒生重)/ 能量摄入 * 1000,
    Dark_vegetables_score = pmin(Dark_vegetables / (0.9 * (225 / 2.3)) * 5, 5),
    
    # 水果,满分为10分
    Fruits = 水果类生重 / 能量摄入 * 1000,
    Fruits_score = pmin(Fruits / (1.1 * (300 / 3)) * 10, 10),
    
    # 奶类
    Dairy = 奶制品类生重 / 能量摄入 * 1000,
    Dairy_score = pmin(Dairy / (0.5 * (300 / 1.2)) * 5, 5),
    
    # 大豆
    Soybean = ((coalesce(硬豆腐生重, 0) +coalesce(软豆腐生重, 0) +coalesce(豆腐干生重, 0) +coalesce(豆腐花生重, 0)) / 5 + coalesce(豆浆生重, 0) / 20 +
                 coalesce(鲜黄豆生重, 0)) / 能量摄入 * 1000,
    Soybean_score = pmin(Soybean / (0.4 * (15 / 0.8)) * 5, 5),
    
    # 海产品
    Aquatic_products = 鱼类生重 / 能量摄入 * 1000,
    Aquatic_products_score = pmin(Aquatic_products / (0.6 * (50 / 1.1)) * 5, 5),
    
    # 禽肉类
    Poultry = 禽肉类生重 / 能量摄入 * 1000,
    Poultry_score = pmin(Poultry / (0.3 * (50 / 1.1)) * 5, 5),
    
    # 蛋类
    Egg = 蛋生重 / 能量摄入 * 1000,
    Egg_score = pmin(Egg / (0.5 * (50 / 1.1)) * 5, 5),
    
    # 坚果类
    Nuts = 坚果生重 / 能量摄入 * 1000,
    Nuts_score = pmin(Nuts / (0.4 * (10 / 1.0)) * 5, 5),
    
    #***********************反向赋分
    
    # 红肉类
    Red_meat = 红肉类生重 / 能量摄入 * 1000,
    Red_meat_score = pmax(pmin((( 3.5 * (50 / 1.1)) - Red_meat ) /
                                 ((3.5 * (50 / 1.1)) - (0.4 * (50 / 1.1))) * 5, 5),0),
    
    # 烹调油，10分
    Cooking_oils = 烹调油生重 / 能量摄入 * 1000,
    Cooking_oils_score = pmax(pmin((32.6 - Cooking_oils) / (32.6 - 15.6) * 10, 10),0),
    
    # 钠盐， mg/1000 kcal, 10分
    Sodium = Na / 能量摄入 * 1000,
    Sodium_score = pmax(pmin( (3608 - Sodium) / (3608 - 1000) * 10, 10),0),
    
    # 添加糖,数据库中没有
    
    # 酒
    Alcohol_score =
      case_when(
        # 男性
        Sex_F0 == 1 ~
          pmax(
            pmin(
              (60 - 酒精g) / (60 - 25) * 5,
              5
            ),
            0
          ),
        
        # 女性
        Sex_F0 == 0 ~
          pmax(
            pmin(
              (40 - 酒精g) / (40 - 15) * 5,
              5
            ),
            0
          ),
        
        TRUE ~ NA_real_
      )
  )
# 提取ID和以_score结尾的变量,并计算总分
Food_score <- Food_Nutrients %>%
  dplyr::select(
    ID,
    ends_with("_score")
  ) %>%
  mutate(
    Total_score =
      rowSums(
        dplyr::select(., ends_with("_score")),
        na.rm = TRUE
      )
  )
####汇总膳食要用的变量####
data_food <- merge(Food_Nutrients[,c("ID","咖啡生重")],Food_score[,c("ID","Total_score")], by = "ID")
data_food$coffee <- ifelse(data_food$咖啡生重 == 0, 0, 1)
####*****汇总所有混杂，定性定量定义好*****####
Cov_F0_final <- merge(Cov_F0_final, data_food, by = "ID")

convert_types <- function(data, numeric = NULL, factor = NULL) {
  for(v in numeric) if(v %in% names(data)) data[[v]] <- as.numeric(as.character(data[[v]]))
  for(v in factor) if(v %in% names(data)) data[[v]] <- as.factor(data[[v]])
  return(data)
}

# 使用
Cov_F0_final <- convert_types(Cov_F0_final, 
                              numeric = c("Met_mean","Energy_mean","Age_F0","Protein_mean","Total_score","Fiber_soluble_mean"), 
                              factor = c("Sex_F0","Smoke_F0","Alcohol_F0","Tea_F0","Calcium_F0","Vitamin_F0","Disease_F0","coffee","Income_F0","Education_F0",
                                         "Fracture_follow","Fracture_Medicine","Estrogen_F0","Menopause_F0","Fracture_F0","Oil_sup_follow")
                              )

# 剔除收入取值为99的观测
Cov_F0_final <- Cov_F0_final[Cov_F0_final$Income_F0 != 99,]
####身高数据####
Height_F123 <- Cov2_wide[, c("ID", "身高_F1","身高_F2","身高_F3")]
Height_F4 <- read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_基本信息_F01234_20230214_VF.xlsx",sheet="F4")
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
####*********饮茶频率*********####
# 因为变量"每周冲茶次数_F2"前面很多缺失值，所以会被识别为分类变量，所以要加上col_types = c(每周冲茶次数_F2 = "numeric")
Tea <- as.data.frame(read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_饮茶相关变量索取_更新_VF.xlsx", sheet=1, guess_max = 10000))
# F2的两个变量都不为缺失值的：没有
Tea_both <- Tea[!is.na(Tea$每周冲茶次数_F2) & !is.na(Tea$平均每周泡茶次数_F2), ]
# 取两个变量全集
Tea$每周泡茶次数_F2 <- dplyr::coalesce(
  Tea$每周冲茶次数_F2,
  Tea$平均每周泡茶次数_F2
)

Tea$过去一年茶叶量斤_F2 <- Tea$过去一年茶叶量斤_F2 + Tea$过去一年茶叶量两_F2/10

Tea$过去一年消费茶叶量斤_F2 <- dplyr::coalesce(
  Tea$过去一年茶叶量斤_F2,
  Tea$过去一年消费茶叶数量_F2
)

Tea$喝茶浓淡程度_F2 <- dplyr::coalesce(
  Tea$一般喝茶浓淡程度_F2,
  Tea$喝茶浓淡_F2
)

####F0####
Tea_F0_filter <- Tea %>%
  mutate(
    
    累计喝茶总年限_F0 =
      if_else(
        是否喝茶_F0 == 0,
        0,
        coalesce(累计喝茶年限_F0, 0) +
          coalesce(累计喝茶月份_F0, 0) / 12
      ),
    
    过去一年茶叶总量_斤_F0 =
      if_else(
        过去一年是否喝茶_F0 == 0,
        0,
        coalesce(过去一年茶叶量斤_F0, 0) +
          coalesce(过去一年茶叶量两_F0, 0) / 10
      )
  ) %>%
  
  rowwise() %>%
  mutate(
    常饮用茶种类_F0 =
      str_c(
        na.omit(c(
          ifelse(绿茶_F0 == 1, "绿茶", NA),
          ifelse(红茶_F0 == 1, "红茶", NA),
          ifelse(乌龙_F0 == 1, "乌龙茶", NA),
          ifelse(其他茶_F0 == 1, 其他茶叶名称_F0, NA)
        )),
        collapse = "；"
      )
  ) %>%
  ungroup() %>%
  
  dplyr::select(
    CODE_F0,
    是否喝茶_F0,
    每周冲茶次数_F0,
    喝茶浓淡_F0,
    累计喝茶总年限_F0,
    过去一年茶叶总量_斤_F0,
    常饮用茶种类_F0
  )
min(Tea_F0_filter$累计喝茶总年限_F0[Tea_F0_filter$累计喝茶总年限_F0 != 0], na.rm = TRUE)
colnames(Tea_F0_filter)
####F1####
Tea_F1_filter <- Tea %>%
  
  mutate(
    过去一年茶叶总量_斤_F1 =
      coalesce(过去一年茶叶量斤_F1, 0)
  ) %>%
  
  rowwise() %>%
  mutate(
    常饮用茶种类_F1 =
      str_c(
        na.omit(c(
          ifelse(绿茶_F1 == 1, "绿茶", NA),
          ifelse(红茶_F1 == 1, "红茶", NA),
          ifelse(乌龙_F1 == 1, "乌龙茶", NA),
          ifelse(其他茶_F1 == 1, 其他茶叶名称_F1, NA)
        )),
        collapse = "；"
      )
  ) %>%
  ungroup() %>%
  
  dplyr::select(
    CODE_F1,
    是否喝茶_F1,
    每周冲茶次数_F1,
    喝茶浓淡_F1,
    过去一年茶叶总量_斤_F1,
    常饮用茶种类_F1
  )
colnames(Tea_F1_filter)
####F2####
library(dplyr)
library(stringr)

Tea_F2_filter <- Tea %>%
  mutate(
    
    每周泡茶次数_F2 =
      coalesce(
        每周冲茶次数_F2,
        平均每周泡茶次数_F2
      ),
    
    过去一年茶叶总量_斤_F2 =
      if_else(
        过去一年是否经常喝茶_F2 == 0,
        0,
        过去一年消费茶叶量斤_F2
      ),
    
    # 清理“其他茶叶名称_F2”
    其他茶叶名称_F2_clean = str_trim(as.character(其他茶叶名称_F2)),
    其他茶叶名称_F2_clean = na_if(其他茶叶名称_F2_clean, ""),
    其他茶叶名称_F2_clean = na_if(其他茶叶名称_F2_clean, "NA"),
    其他茶叶名称_F2_clean = na_if(其他茶叶名称_F2_clean, "无"),
    其他茶叶名称_F2_clean = na_if(其他茶叶名称_F2_clean, "没有")
  ) %>%
  rowwise() %>%
  mutate(
    常饮用茶种类_F2 =
      str_c(
        na.omit(c(
          ifelse(绿茶_F2 == 1, "绿茶", NA),
          ifelse(红茶_F2 == 1, "红茶", NA),
          ifelse(乌龙_F2 == 1, "乌龙茶", NA),
          
          # 关键修改：
          # 只要“其他茶叶名称_F2”里有茶名，就加入
          其他茶叶名称_F2_clean
        )),
        collapse = "；"
      ),
    
    # 如果最终是空字符串，改成 NA
    常饮用茶种类_F2 = ifelse(
      常饮用茶种类_F2 == "",
      NA,
      常饮用茶种类_F2
    )
  ) %>%
  ungroup() %>%
  dplyr::select(
    CODE_F2,
    过去一年是否经常喝茶_F2,
    每周泡茶次数_F2,
    喝茶浓淡程度_F2,
    过去一年茶叶总量_斤_F2,
    常饮用茶种类_F2
  )
####F3####
Tea_F3_filter <- Tea %>%
  
  mutate(
    过去一年茶叶总量_斤_F3 =
      coalesce(过去一年茶叶量斤_F3, 0)
  ) %>%
  
  rowwise() %>%
  mutate(
    常饮用茶种类_F3 =
      str_c(
        na.omit(c(
          ifelse(绿茶_F3 == 1, "绿茶", NA),
          ifelse(红茶_F3 == 1, "红茶", NA),
          ifelse(乌龙_F3 == 1, "乌龙茶", NA),
          ifelse(普洱_F3 == 1, "普洱茶", NA),
          ifelse(花茶_F3 == 1, "花茶", NA),
          ifelse(其他茶_F3 == 1, 其他茶叶名称_F3, NA)
        )),
        collapse = "；"
      )
  ) %>%
  ungroup() %>%
  
  dplyr::select(
    CODE_F3,
    是否喝茶_F3,
    每周冲茶次数_F3,
    喝茶浓淡_F3,
    过去一年茶叶总量_斤_F3,
    常饮用茶种类_F3
  )


colnames(Tea_F0_filter)
colnames(Tea_F1_filter)
colnames(Tea_F2_filter)
colnames(Tea_F3_filter)
head(Tea_F0_filter)
####饮茶频率长数据####
process_tea <- function(data,
                        code_col,
                        drink_col,
                        freq_col,
                        visit,
                        remove_prefix = NULL) {
  
  data %>%
    
    transmute(
      
      ID =
        if(is.null(remove_prefix)) {
          .data[[code_col]]
        } else {
          sub(remove_prefix, "", .data[[code_col]])
        },
      
      visit = visit,
      
      tea_drink = .data[[drink_col]],
      
      tea_freq = case_when(
        
        .data[[drink_col]] == 0 ~ 0,
        
        .data[[drink_col]] == 1 ~
          .data[[freq_col]],
        
        TRUE ~ NA_real_
      )
    )
}

Tea_F0_long <- process_tea(
  Tea_F0_filter,
  code_col = "CODE_F0",
  drink_col = "是否喝茶_F0",
  freq_col = "每周冲茶次数_F0",
  visit = "F0"
)

Tea_F1_long <- process_tea(
  Tea_F1_filter,
  code_col = "CODE_F1",
  drink_col = "是否喝茶_F1",
  freq_col = "每周冲茶次数_F1",
  visit = "F1",
  remove_prefix = "^F1"
)

Tea_F2_long <- process_tea(
  Tea_F2_filter,
  code_col = "CODE_F2",
  drink_col = "过去一年是否经常喝茶_F2",
  freq_col = "每周泡茶次数_F2",
  visit = "F2",
  remove_prefix = "^F2"
)

Tea_F3_long <- process_tea(
  Tea_F3_filter,
  code_col = "CODE_F3",
  drink_col = "是否喝茶_F3",
  freq_col = "每周冲茶次数_F3",
  visit = "F3",
  remove_prefix = "^F3"
)

Tea_long_prime <- bind_rows(
  Tea_F0_long,
  Tea_F1_long,
  Tea_F2_long,
  Tea_F3_long
) %>%
  
  dplyr::filter(stats::complete.cases(.)) %>%
  
  mutate(
    time = c(
      F0 = 0,
      F1 = 1,
      F2 = 2,
      F3 = 3
    )[visit]
  )

Tea_long_prime <- Tea_long_prime %>%
  mutate(
    tea_freq_group = case_when(
      
      tea_freq == 0 ~ 0,
      
      tea_freq > 0 & tea_freq < 7 ~ 1,
      
      tea_freq >= 7 ~ 2,
      
      TRUE ~ NA_real_
    )
  )
Tea_long_prime$tea_freq_group <- as.factor(Tea_long_prime$tea_freq_group)

table(Tea_long_prime$tea_freq_group)

# table(Tea_long$visit)
# table(Tea_long$visit, Tea_long$tea_freq_group)
# round(
#   prop.table(
#     table(Tea_long$visit, Tea_long$tea_freq_group),
#     margin = 1
#   ) * 100,
#   1
# )
####简单纵向分组####
Tea_long_prime$ID_numeric <- as.numeric(sub("^NL", "", Tea_long_prime$ID))
Tea_long_prime$tea_freq <- as.numeric(Tea_long_prime$tea_freq)
Tea_long_prime_filter <- Tea_long_prime %>%
  group_by(ID) %>%
  filter(sum(!is.na(tea_freq)) >= 3) %>%
  ungroup()


Tea_summary <- Tea_long_prime_filter %>%
  group_by(ID) %>%
  summarise(
    
    # 统计 tea_freq >= 7 的次数
    n_high = sum(tea_freq >= 7, na.rm = TRUE),
    
    tea_freq_group = case_when(
      
      # 所有随访都不喝茶
      all(tea_freq == 0, na.rm = TRUE) ~ 0,
      
      # 至少两次随访 tea_freq >= 7
      n_high >= 2 ~ 2,
      
      # 其他情况
      TRUE ~ 1
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::select(-n_high)

Tea_summary$tea_freq_group <- as.factor(Tea_summary$tea_freq_group)
colnames(Tea_summary) <- recode(colnames(Tea_summary),
                                "tea_freq_group" = "tea_freq_group_summary")

table(Tea_summary$tea_freq_group_summary)
round(prop.table(table(Tea_summary$tea_freq_group_summary)) * 100,1)
####计算累计平均饮茶频率####
Tea_mean_freq <- Tea_long_prime_filter %>%
  group_by(ID) %>%
  summarise(
    
    # 平均饮茶频率
    tea_freq_mean = mean(tea_freq, na.rm = TRUE),
    
    # 实际参与平均计算的随访次数
    n_visit = sum(!is.na(tea_freq)),
    
    .groups = "drop"
  ) %>%
  mutate(
    tea_freq_group = case_when(
      
      tea_freq_mean == 0 ~ 0,
      
      tea_freq_mean > 0 & tea_freq_mean < 7 ~ 1,
      
      tea_freq_mean >= 7 ~ 2,
      
      TRUE ~ NA_real_
    )
  )

table(Tea_mean_freq$tea_freq_group)
round(prop.table(table(Tea_mean_freq$tea_freq_group)) * 100,1)

Tea_mean_freq$tea_freq_group <- as.factor(Tea_mean_freq$tea_freq_group)
####轨迹分析####
Tea_zero_all <- Tea_long_prime_filter %>%
  group_by(ID_numeric) %>%             # 1. 按照 ID 分组
  filter(all(tea_freq == 0)) %>%       # 2. 筛选：该组内所有的 tea_freq 是否都等于 0
  ungroup()                            # 3. 取消分组（好习惯）
a <- unique(Tea_zero_all$ID)
table(Tea_zero_all$visit)

Tea_zero_all <- Tea_zero_all[!duplicated(Tea_zero_all$ID), ]
Tea_zero_all$class <- "0"

Tea_long_prime_drink <- Tea_long_prime_filter[!Tea_long_prime_filter$ID %in% Tea_zero_all$ID, ]
table(Tea_long_prime_drink$visit)

Class_all <- run_hlme_traj_plot(
  data = Tea_long_prime_drink ,
  outcome = "tea_freq",
  time_var = "time",
  id_var = "ID_numeric",
  xlab_label = "Follow-up",
  ylab_label = "Tea consumption",
  max_ng = 5,
  ng_plot = 2,
  legend.x = 1.2,
  legend.y = 6,
  height = 6,
  pdf_file = "Tea_tra.pdf"
)

Class_all$Class$class <- as.factor(Class_all$Class$class)

Tea_long_prime_drink_ID <- Tea_long_prime_drink[!duplicated(Tea_long_prime_drink$ID), ]
data_list <- list(Tea_long_prime_drink_ID, Class_all$Class)
Tea_tra <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID_numeric") %>%
  dplyr::filter(stats::complete.cases(.)) %>%
  dplyr::select(
    ID,
    class
  )

Tea_tra_all <- rbind(Tea_zero_all[,c("ID","class")], Tea_tra)

table(Tea_tra_all$class)
round(prop.table(table(Tea_tra_all$class)) * 100,1)
####只取基线饮茶频率####
Tea_F0_final <- Tea_F0_long[complete.cases(Tea_F0_long),]
Tea_F0_final_NL4 <- Tea_F1_long %>%
  filter(grepl("^NL4", ID)) %>%
  dplyr::filter(stats::complete.cases(.)) 

Tea_F0_final <- rbind(Tea_F0_final, Tea_F0_final_NL4)
table(Tea_F0_final$tea_freq)

Tea_F0_final <- Tea_F0_final %>%
  mutate(
    tea_freq_group = case_when(
      
      tea_freq == 0 ~ 0,
      
      tea_freq > 0 & tea_freq < 7 ~ 1,
      
      tea_freq >= 7 ~ 2,
      
      TRUE ~ NA_real_
    )
  )
table(Tea_F0_final$tea_freq_group)
round(prop.table(table(Tea_F0_final$tea_freq_group)) * 1400,1)

Tea_F0_final$tea_freq_group <- as.factor(Tea_F0_final$tea_freq_group)
####F0血清儿茶素####
Catechin <- read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/血清氧化应激炎症因子血黄酮尿黄酮/20250527_GNHS_F0-3_VF.xlsx",sheet=1)
Catechin <- Catechin[,c("ID","serum_I_catechin_F0", "serum_I_epicatechin_F0", "serum_I_EGC_F0", "serum_I_EGCG_F0", "serum_I_ECG_F0")]
Catechin <- Catechin[complete.cases(Catechin),]

Catechin_sex <- merge(Catechin, Cov_F0[,c("ID","Sex_F0")], by = "ID")

#******查看零值
vars <- setdiff(names(Catechin), "ID")

# 统计0值个数和百分比
zero_summary <- data.frame(
  Variable = vars,
  Zero_N = sapply(Catechin[vars], function(x) sum(x == 0, na.rm = TRUE)),
  Total_N = nrow(Catechin)
)

zero_summary$Zero_Percent <- round(
  100 * zero_summary$Zero_N / zero_summary$Total_N,
  2
)

####*********DXA-骨骼肌*********####
DXA_all <- read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_DXA_follow_10807条记录_20241014_2024.10.18recheck_VF.xlsx")
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

# 二分类
ASM_Height_long <- merge(ASM_Height_long, Cov_F0_final[,c("ID","Sex_F0")], by = "ID")
colnames(ASM_Height_long)
ASM_Height_long <- ASM_Height_long %>%
  mutate(
    low_ASMI = case_when(
      Sex_F0 == 1 & !is.na(ASMI) & ASMI < 7.0 ~ 1L,  # male
      Sex_F0 == 0 & !is.na(ASMI) & ASMI < 5.4 ~ 1L,  # female
      Sex_F0 %in% c(0, 1) & !is.na(ASMI) ~ 0L,
      TRUE ~ NA_integer_
    ))

ASM_Height_long <- dplyr::select(ASM_Height_long, -Sex_F0)
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
####随访过程中出现低ASMI的观测####
ASMI_wide_low <- ASMI_wide[,c("ID", "F4ASMI", "F3ASMI", "F2ASMI", "F1ASMI","WBTOT_FAT_mean_F1234")]
ASMI_wide_low <- merge(ASMI_wide_low, Cov_F0_final[,c("ID","Sex_F0")], by = "ID")

asmi_cols <- c("F1ASMI", "F2ASMI", "F3ASMI", "F4ASMI")

ASMI_wide_low <- ASMI_wide_low %>%
  mutate(
    low_ASMI = case_when(
      # 四次随访 ASMI 全部缺失，不能定义为 0
      if_all(all_of(asmi_cols), is.na) ~ NA_integer_,
      
      # 男性：任意一次 ASMI < 7.0
      Sex_F0 == 1 & if_any(all_of(asmi_cols), ~ !is.na(.x) & .x < 7.0) ~ 1L,
      
      # 女性：任意一次 ASMI < 5.4
      Sex_F0 == 0 & if_any(all_of(asmi_cols), ~ !is.na(.x) & .x < 5.4) ~ 1L,
      
      # 有至少一次 ASMI 不缺失，但均未达到低 ASMI 标准
      TRUE ~ 0L
    )
  )
####*****************肌少症评定****************####
#（1）ASM（DXA）:男<7.0 kg/m2;女<5.4 kg/m2。（2）握力：男<28 kg；女<18 kg。（3）5次坐立试验：≥12 s
SPPB1 <- read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/SPPB/20251206_GNHS_F0-3.xlsx")
colnames(SPPB1)
a <- SPPB1[,c("ID", "握力1_F2", "握力2_F2")]
a <- a[complete.cases(a),]

SPPB2 <- read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/SPPB/20251206_GNHS_F4.xlsx")
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
table(SPPB_F2$sarcopenia_F2)

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

####SPPB问卷计分####
#***************************F2-3
SPPB1 <- read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/SPPB/20251206_GNHS_F0-3.xlsx")
SPPB1 <- SPPB1 %>%
  dplyr::select(
    ID,
    matches("握力|日常步速|重复坐凳试验时间|站立平衡")
  )
vars <- c(
  "重复坐凳试验时间_F2",
  "日常步速1_F3",
  "日常步速2_F3",
  "重复坐凳试验时间_F3"
)
for(i in vars){
  SPPB1 <- clean_var(SPPB1, i)
}


vars <- c(
  "日常步速1_F2",
  "日常步速2_F2",
  "日常步速1_F3",
  "日常步速2_F3"
)
SPPB1[vars] <- lapply(
  SPPB1[vars],
  function(x) {
    x <- suppressWarnings(as.numeric(as.character(x)))
    x[x > 10] <- NA
    return(x)
  }
)
#***************************F4
SPPB2 <- read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/SPPB/20251206_GNHS_F4.xlsx")
SPPB2 <- SPPB2 %>%
  filter(startsWith(CODE, "F4"))
SPPB2$ID <- substr(SPPB2$CODE, 3, nchar(SPPB2$CODE))
SPPB2 <- SPPB2 %>%
  dplyr::select(
    ID,
    matches("握力|日常步速|重复坐凳试验时间|站立平衡")
  )
a <- SPPB2[SPPB2$站立平衡 == 999,]
# 有一个观测站立平衡取值为“999”,其实就是0分
SPPB2$站立平衡[SPPB2$站立平衡 == 999] <- 0
colnames(SPPB2)

vars <- c(
  "日常步速1",
  "日常步速2",
  "重复坐凳试验时间"
)
for(i in vars){
  SPPB2 <- clean_var(SPPB2, i)
}

vars <- c(
  "日常步速1",
  "日常步速2"
)
SPPB2[vars] <- lapply(
  SPPB2[vars],
  function(x) {
    x <- suppressWarnings(as.numeric(as.character(x)))
    x[x > 10] <- NA
    return(x)
  }
)
#***************************赋值
# 1. 站立平衡计分
score_balance <- function(balance_val) {
  
  balance_val <- suppressWarnings(as.numeric(as.character(balance_val)))
  
  score <- ifelse(
    is.na(balance_val),
    NA,
    ifelse(balance_val >= 1 & balance_val <= 4,
           balance_val,
           NA)
  )
  
  return(score)
}


# 2. 步行速度计分（4米）,取最快的一次
score_gait_speed <- function(gait1, gait2) {
  
  # 保留原始字符
  gait1_raw <- as.character(gait1)
  gait2_raw <- as.character(gait2)
  
  # 转小写并去空格
  gait1_clean <- trimws(tolower(gait1_raw))
  gait2_clean <- trimws(tolower(gait2_raw))
  
  # 特殊值 → 0分
  max_flag <- gait1_clean %in% c("max", "无", "不便", "mx", "min") |
    gait2_clean %in% c("max", "无", "不便", "mx", "min")
  
  # 数值化
  gait1_num <- suppressWarnings(as.numeric(gait1_raw))
  gait2_num <- suppressWarnings(as.numeric(gait2_raw))
  
  # 速度取最快的一次,记得转换为速度，原始数据是时间
  Speed_m_s <- pmax(4/gait1_num, 4/gait2_num, na.rm = TRUE)
  
  # 两个都缺失时修正
  Speed_m_s[is.infinite(Speed_m_s)] <- NA
  
  # 取最快一次（时间最短）
  gait_fast <- pmin(gait1_num, gait2_num, na.rm = TRUE)
  
  # 两个都缺失时修正
  gait_fast[is.infinite(gait_fast)] <- NA
  
  # SPPB评分
  score <- ifelse(
    is.na(gait_fast),
    NA,
    ifelse(gait_fast < 4.82, 4,
           ifelse(gait_fast <= 6.20, 3,
                  ifelse(gait_fast <= 8.70, 2,
                         ifelse(gait_fast > 8.70, 1, NA))))
  )
  
  # 特殊值 → 0分
  score[max_flag] <- 0
  
  # 特殊值时速度也设为NA
  Speed_m_s[max_flag] <- NA
  
  # 返回结果
  return(
    data.frame(
      gait_score = score,
      gait_speed = Speed_m_s
    )
  )
}

# 3. 重复坐凳试验计分
score_chair_stand <- function(chair_time) {
  
  # 保留原始字符
  chair_raw <- as.character(chair_time)
  
  # 转小写并去空格
  chair_clean <- trimws(tolower(chair_raw))
  
  # 特殊值标记 → 0分
  max_flag <- chair_clean %in% c("max","无","不便", "mx", "min")
  
  # 数值化
  chair_time <- suppressWarnings(as.numeric(chair_raw))
  
  # 正常评分
  score <- ifelse(
    is.na(chair_time),
    NA,
    ifelse(chair_time <= 11.19, 4,
           ifelse(chair_time <= 13.69, 3,
                  ifelse(chair_time <= 16.69, 2,
                         ifelse(chair_time > 16.69, 1, NA))))
  )
  
  # 特殊值 → 0分
  score[max_flag] <- 0
  
  result <- data.frame(
    chair_time = chair_time,
    chair_score = score
  )
  
  return(result)
}

# 4. 握力
get_grip_strength <- function(grip1, grip2) {
  
  # 保留原始值
  grip1_raw <- grip1
  grip2_raw <- grip2
  
  # 转换数值
  grip1 <- suppressWarnings(as.numeric(grip1_raw))
  grip2 <- suppressWarnings(as.numeric(grip2_raw))
  
  # 99视为缺失
  grip1[grip1 == 99] <- NA
  grip2[grip2 == 99] <- NA
  
  # 取最大握力
  grip_max <- pmax(grip1, grip2, na.rm = TRUE)
  
  # 两个都缺失时，pmax会返回 -Inf，需要改回NA
  grip_max[is.infinite(grip_max)] <- NA
  
  return(grip_max)
}

# SPPB1 计分（F2 和 F3 两个时间点）
SPPB1_scored <- SPPB1 %>%
  mutate(
    # --- F2 ---
    balance_score_F2 = score_balance(站立平衡_F2),
    gait_speed_score_F2 = score_gait_speed(日常步速1_F2, 日常步速2_F2)$gait_score,
    chair_stand_score_F2 = score_chair_stand(重复坐凳试验时间_F2)$chair_score,
    SPPB_total_F2 = balance_score_F2 + gait_speed_score_F2 + chair_stand_score_F2,
    SPPB_category_F2 = ifelse(SPPB_total_F2 >= 10, "行动能力正常", "行动能力障碍"),
    grip_max_F2 = get_grip_strength(握力1_F2, 握力2_F2),
    gait_speed_F2 = score_gait_speed(日常步速1_F2, 日常步速2_F2)$gait_speed,
    chair_time_F2 = score_chair_stand(重复坐凳试验时间_F2)$chair_time,
    
    
    # --- F3 ---
    balance_score_F3 = score_balance(站立平衡_F3),
    gait_speed_score_F3 = score_gait_speed(日常步速1_F3, 日常步速2_F3)$gait_score,
    chair_stand_score_F3 = score_chair_stand(重复坐凳试验时间_F3)$chair_score,
    SPPB_total_F3 = balance_score_F3 + gait_speed_score_F3 + chair_stand_score_F3,
    SPPB_category_F3 = ifelse(SPPB_total_F3 >= 10, "行动能力正常", "行动能力障碍"),
    grip_max_F3 = get_grip_strength(握力1_F3, 握力2_F3),
    gait_speed_F3 = score_gait_speed(日常步速1_F3, 日常步速2_F3)$gait_speed,
    chair_time_F3 = score_chair_stand(重复坐凳试验时间_F3)$chair_time
    
  )

a1 <- check_df_summary(SPPB1_scored)
a <- SPPB1_scored[,c("ID","grip_max_F2")]
a <- a[complete.cases(a),]

# SPPB2 计分（F4 时间点）
SPPB2_scored <- SPPB2 %>%
  mutate(
    balance_score_F4 = score_balance(站立平衡),
    gait_speed_score_F4 = score_gait_speed(日常步速1, 日常步速2)$gait_score,
    chair_stand_score_F4 = score_chair_stand(重复坐凳试验时间)$chair_score,
    SPPB_total_F4 = balance_score_F4 + gait_speed_score_F4 + chair_stand_score_F4,
    SPPB_category_F4 = ifelse(SPPB_total_F4 >= 10, "行动能力正常", "行动能力障碍"),
    grip_max_F4 = get_grip_strength(握力1, 握力2),
    gait_speed_F4 = score_gait_speed(日常步速1, 日常步速2)$gait_speed,
    chair_time_F4 = score_chair_stand(重复坐凳试验时间)$chair_time
  )
a2 <- check_df_summary(SPPB2_scored)


SPPB_grip <- left_join(SPPB1_scored[,c("ID","SPPB_total_F2","SPPB_total_F3","grip_max_F2","grip_max_F3", "gait_speed_F2", "gait_speed_F3","chair_time_F2","chair_time_F3")], 
                   SPPB2_scored[,c("ID","SPPB_total_F4", "grip_max_F4","gait_speed_F4","chair_time_F4")], by = "ID")

####长数据####
SPPB_grip_long <- SPPB_grip %>%
  pivot_longer(
    cols = -ID,
    names_to = c(".value", "Followup"),
    names_pattern = "(.*)_(F\\d+)"
  ) %>%
  arrange(ID, Followup) 

SPPB_grip_long$SPPB_category <- ifelse(SPPB_grip_long$SPPB_total >= 10, 0, 1)
table(SPPB_grip_long$SPPB_category)

a <- SPPB_grip_long[,c("ID","grip_max","Followup")]
a <- a[complete.cases(a),]
table(a$Followup)

a <- SPPB_grip_long[,c("ID","gait_speed","Followup")]
a <- a[complete.cases(a),]
table(a$Followup)
table(a$Followup)

#******转换为二分类
colnames(SPPB_grip_long)
SPPB_grip_long <- merge(SPPB_grip_long, Cov_F0_final[,c("ID","Sex_F0")], by = "ID")

SPPB_grip_long <- SPPB_grip_long %>%
  mutate(
    low_grip = case_when(
      Sex_F0 == 1 & !is.na(grip_max) & grip_max < 28 ~ 1L,
      Sex_F0 == 0 & !is.na(grip_max) & grip_max < 18 ~ 1L,
      Sex_F0 %in% c(0, 1) & !is.na(grip_max) ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    # AWGS: low walking speed
    # walking speed < 1.0 m/s
    low_gait_speed = case_when(
      !is.na(gait_speed) & gait_speed < 1.0 ~ 1L,
      !is.na(gait_speed) ~ 0L,
      TRUE ~ NA_integer_
    )
  )

SPPB_grip_long <- dplyr::select(SPPB_grip_long, -Sex_F0)
####随访过程中出现低握力和低步速的观测####
SPPB_grip_wide <- SPPB_grip_long %>%
  group_by(ID) %>%
  summarise(
    low_grip = case_when(
      any(low_grip == 1, na.rm = TRUE) ~ 1L,
      all(is.na(low_grip)) ~ NA_integer_,
      TRUE ~ 0L
    ),
    
    low_gait_speed = case_when(
      any(low_gait_speed == 1, na.rm = TRUE) ~ 1L,
      all(is.na(low_gait_speed)) ~ NA_integer_,
      TRUE ~ 0L
    ),
    .groups = "drop"
  )
####*********菌群*********####
####各水平####
Micro <- as.data.frame(read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_metaphlan4_taxa_3678_IDcorrect_251125_VF.xlsx",sheet=1))
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
colnames(Micro_selected)
#****************菌群过滤
Micro_selected_filter <- perform_micro_filter(Micro_selected,
                                              Micro_selected,
                                              exclude_cols = c("ID", "Times"),
                                              metadata_cols = c("ID", "Times"),
                                              min_abundance = 0.0001,  
                                              min_sample_detection = 0.1 
)
#****************clr转换，因为clr是按行转换，所以跟整体分布无关
Micro_selected_filter_clr <- perform_micro_clr(Micro_selected_filter$filtered_data, exclude_cols = c("ID", "Times"))
Micro_selected_filter_clr <-Micro_selected_filter_clr[grepl("^NL", Micro_selected_filter_clr$ID),]

####多样性####
Diversity <- as.data.frame(read.csv("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_metaphlan3_diversity_2664_221028_VF.csv",header = TRUE))
Diversity <- Diversity %>%
  mutate(
    Times = substr(id, 1, 2),
    ID    = substr(id, 3, nchar(id))
  ) %>%
  relocate(ID, Times, .after = 1)
####微生物功能-unstra####
Micro_func <- readr::read_csv(
  "D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_humann3_pathabundance_2659_unstra_230606_VF.tsv",
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
head(Micro_func)
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
Micro_func_filter_clr <- Micro_func_filter_clr[grepl("^NL", Micro_func_filter_clr$ID),]
table(Micro_func_filter_clr$Times)
colnames(Micro_func_filter_clr)
####微生物功能-stra####
Micro_func_stra <- readr::read_csv(
  "D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_humann3_pathabundance_2659_stra_230606_VF.tsv",
  locale = locale(encoding = "UTF-8")
)
colnames(Micro_func_stra)
####*********血清蛋白质*********####
F0protein <- as.data.frame(read.csv("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_血清蛋白_基线3415人_VF.csv",header = TRUE))

F2protein <- as.data.frame(read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_血清蛋白F2_2567人_VF.xlsx",sheet = 1))
F2protein[F2protein == "NA" & !is.na(F2protein)] <- NA

F3protein1 <- as.data.frame(read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_F3蛋白组数据索取_更新_VF.xlsx",sheet = 1))
F3protein2 <- as.data.frame(read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_F3蛋白组数据索取_更新_VF.xlsx",sheet = 2))
colnames(F3protein2) <- sub("_.*", "", colnames(F3protein2))
#*************************F3蛋白质合并
phase1 <- as.data.frame(read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_F3蛋白组数据索取_更新_VF.xlsx",sheet = 3))
phase1 <- phase1[,c("patient_ID", "sample_collection_batch","phase")]
phase1 <- phase1[phase1$sample_collection_batch == "F3",]
phase1_unique <- phase1[!duplicated(phase1$patient_ID), ]
F3protein1_phase <- merge(phase1_unique, F3protein1, by.x = c("patient_ID", "sample_collection_batch"), by.y = c("sampleid","followup"))
colnames(F3protein1_phase) <- recode(colnames(F3protein1_phase),
                                     "patient_ID" = "ID",
                                     "sample_collection_batch" = "followup",
                                     "phase" = "Phase")
F3protein1_phase <- dplyr::select(F3protein1_phase,-c("time"))

phase2 <- as.data.frame(read_excel("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_F3蛋白组数据索取_更新_VF.xlsx",sheet = 4))
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
####看看各个蛋白质的缺失值个数####
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
####蛋白质剔除和填补宽数据和长数据####
process_protein_data <- function(F3protein_all,
                                 exclude_cols,
                                 detection_cutoff = 0.80) {
  
  library(dplyr)
  library(tidyr)
  
  ## 获取蛋白列
  protein_cols <- setdiff(
    colnames(F3protein_all),
    exclude_cols
  )
  
  ## 1. 筛选检测率 >= cutoff 的蛋白

  keep_proteins <- F3protein_all %>%
    
    summarise(
      across(
        all_of(protein_cols),
        ~ mean(!is.na(.))
      )
    ) %>%
    
    pivot_longer(
      cols = everything(),
      names_to = "protein",
      values_to = "detection_rate"
    ) %>%
    
    filter(
      detection_rate >= detection_cutoff
    ) %>%
    
    pull(protein)
  
  ## 2. 保留原始变量

  F3protein_all_filtered <- F3protein_all %>%
    
    dplyr::select(
      all_of(c(exclude_cols, keep_proteins))
    )
  
  ## 3. 转数值

  F3protein_all_filtered <- F3protein_all_filtered %>%
    
    mutate(
      across(
        all_of(keep_proteins),
        as.numeric
      )
    )
  
  ## 4. 缺失值填补（生成新变量 _fill）

  fill_df <- F3protein_all_filtered %>%
    
    transmute(
      
      across(
        all_of(keep_proteins),
        
        ~ {
          min_val <- min(., na.rm = TRUE)
          
          if (is.finite(min_val)) {
            ifelse(is.na(.), min_val / 2, .)
          } else {
            .
          }
        },
        
        .names = "{.col}_fill"
      )
    )
  
  ## 5. log2转换（生成新变量 _fill_log2）

  log2_df <- fill_df %>%
    
    mutate(
      across(
        everything(),
        ~ log2(.),
        .names = "{.col}_log2"
      )
    )
  
  ## 6. 合并所有结果

  final_df <- bind_cols(
    F3protein_all_filtered,
    fill_df,
    log2_df
  )
  
  return(final_df)
}
#******************宽数据
F3protein_fill <- process_protein_data(F3protein_all, c("ID", "followup", "Phase"), 0.5)

F2protein_fill <- process_protein_data(F2protein, c("CODE", "ID", "followup", "Phase"), 0.5)

F0protein_fill <- process_protein_data(F0protein, c("CODE", "ID", "followup", "Phase"), 0.5)
#******************长数据
Val_intersect <- intersect(intersect(colnames(F0protein), colnames(F2protein)),colnames(F3protein_all))

Protein_long_all <- rbind(F0protein[, Val_intersect], F2protein[, Val_intersect], F3protein_all[, Val_intersect])
Protein_long_all$followup <- recode(Protein_long_all$followup,
                                    "basline" = "F0")
Protein_long_all$followup <- recode(Protein_long_all$followup,
                                    "baseline" = "F0")
table(Protein_long_all$followup)
Protein_long_all_fill <- process_protein_data(Protein_long_all, c("ID", "followup", "Phase"), 0.5)

a <- Protein_long_all_fill[Protein_long_all_fill$followup %in% c("F2","F3"),]
####蛋白质名称####
# Protein_names <- toupper(setdiff(colnames(Protein_long_adjusted_data),c("ID", "followup", "Phase")))
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

# writexl::write_xlsx(Protein_names_unique, "D:/OneDrive/Papers/9 饮茶频率&骨骼肌/Database/Protein_names_unique.xlsx")

# 导入蛋白质名称
Protein_names <- read_excel("D:/OneDrive/Papers/9 饮茶频率&骨骼肌/Database/Protein_names_unique.xlsx",sheet=1)

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
Dm_Acc <- read.delim("D:/OneDrive/Papers/8 多组学&骨骼肌/重要文件/uniprotkb_organism_id_7227_AND_reviewed_2025_10_06.tsv", header = TRUE, sep = "\t")
####*****************粪便代谢物*****************####
Fecal_met <- as.data.frame(read.csv("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/Fecal_metabolites.csv",header = TRUE))
Met_names <- as.data.frame(t(Fecal_met[1:2,]))
Met_names <- setNames(Met_names[-1, ], as.character(Met_names[1, ]))

Fecal_met1 <- Fecal_met[,-c(201:ncol(Fecal_met))]
Fecal_met2 <- Fecal_met1[-2,]
Fecal_met2 <- setNames(Fecal_met2[-1, ], as.character(Fecal_met2[1, ]))
colnames(Fecal_met2) <- recode(colnames(Fecal_met2),"sample.x" = "ID")
Fecal_met2$ID <- toupper(Fecal_met2$ID)

library(dplyr)

## 1. 提取 Times 和新的 ID
Fecal_met_final <- Fecal_met2 %>%
  mutate(
    Times = substr(ID, 1, 2),
    ID = substr(ID, 3, nchar(ID))
  ) %>%
  relocate(ID, Times, .before = everything())

## 2. 定义代谢物列：排除 ID 和 Times
met_cols <- setdiff(
  names(Fecal_met_final),
  c("ID", "Times")
)

## 3. 新建填补列：xxx_fill
Fecal_met_final <- Fecal_met_final %>%
  mutate(
    across(
      all_of(met_cols),
      ~ {
        x <- as.numeric(as.character(.x))
        
        ## 取最小正值，避免 0 或 NA 导致 log2 出问题
        min_val <- min(x[x > 0], na.rm = TRUE)
        
        if(is.finite(min_val)){
          ifelse(is.na(x) | x <= 0, min_val / 2, x)
        } else {
          x
        }
      },
      .names = "{.col}_fill"
    )
  )

## 4. 对填补后的列新建 log2 列：xxx_fill_log2
fill_cols <- paste0(met_cols, "_fill")

Fecal_met_final <- Fecal_met_final %>%
  mutate(
    across(
      all_of(fill_cols),
      ~ log2(.x),
      .names = "{.col}_log2"
    )
  )

####*********混杂定义*********####
Covariates_all = c("Age_F0","Sex_F0","WBTOT_FAT_mean_F1234","Total_score","Income_F0","Education_F0","coffee","Met_mean","Protein_mean","Alcohol_F0","Smoke_F0","Calcium_F0","Vitamin_F0","Disease_F0","Estrogen_F0","Menopause_F0","Oil_sup_follow")
Covariates_female = c("Age_F0","WBTOT_FAT_mean_F1234","Total_score","Income_F0","Education_F0","coffee","Met_mean","Protein_mean","Alcohol_F0","Calcium_F0","Vitamin_F0","Disease_F0","Estrogen_F0","Menopause_F0","Oil_sup_follow")
Covariates_male = c("Age_F0","WBTOT_FAT_mean_F1234","Total_score","Income_F0","Education_F0","coffee","Met_mean","Protein_mean","Alcohol_F0","Smoke_F0","Calcium_F0","Vitamin_F0","Disease_F0","Oil_sup_follow")

# 线形混合效应混杂
Covariates_all_lmer = c("Age","Sex_F0","WBTOT_FAT","Total_score","Income_F0","Education_F0","coffee","Met_mean","Protein_mean","Alcohol_F0","Smoke_F0","Calcium_F0","Vitamin_F0","Disease_F0","Estrogen_F0","Menopause_F0","Oil_sup_follow")

#********茶多酚
tea_microbiome_pathways <- c(
  "FERMENTATION-PWY: mixed acid fermentation",
  "P41-PWY: pyruvate fermentation to acetate and lactate I",
  "PWY-5100: pyruvate fermentation to acetate and lactate II",
  "P461-PWY: hexitol fermentation to lactate, formate, ethanol and acetate",
  "ANAEROFRUCAT-PWY: homolactic fermentation"
)

tea_microbiome_pathways_clr <- paste0(tea_microbiome_pathways, "_clr_z")
####儿茶素变量####
Catechin_vars <- c(
  "serum_I_catechin_F0",
  "serum_I_epicatechin_F0",
  "serum_I_EGC_F0",
  "serum_I_EGCG_F0",
  "serum_I_ECG_F0"
)

Catechin_vars1 <- c(
  "serum_I_epicatechin_F0",
  "serum_I_EGC_F0"
)

Catechin_vars2 <- c(
  "serum_I_catechin_F0",
  "serum_I_EGCG_F0",
  "serum_I_ECG_F0"
)

Catechin_vars_T <- c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
Catechin_vars_all <- c("tea_freq_group","serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
####*****************************饮茶&骨骼肌*****************************####
####data1-饮茶纵向&骨骼肌/握力####
data1_long <- merge(ASM_Height_long, Tea_long_prime, by.x = c("ID", "Followup") , by.y = c("ID", "visit"))

data_list <- list(Cov_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], data1_long)
data1_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
# 合并握力等，要用left_join
data1_v2 <- left_join(data1_v1, SPPB_grip_long, by = c("ID","Followup"))
data1_v2 <- left_join(data1_v2, Catechin, by = "ID")

table(SPPB_grip_long$Followup)

data1_v3 <- scale_columns_group(data1_v2, c("ASM","ASMI","grip_max"), "Sex_F0")
data1_v4 <- scale_columns(data1_v3, c("gait_speed","chair_time"))

data1_final <- group_variables3(
  data = data1_v4,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(data1_final)%in% factor_name)
for(i in idx ){
  data1_final[[i]] <-  as.factor(data1_final[[i]])
}


# 统计
table(data1_final$Followup)
a <- unique(data1_final$ID)
a <- data1_final %>%
  dplyr::select(ID, Followup, grip_max) %>%
  filter(complete.cases(.))
table(a$Followup)

a <- data1_final %>%
  dplyr::select(ID, Followup, gait_speed) %>%
  filter(complete.cases(.))
table(a$Followup)
####lmer####
Lmer_results1_all <- process_lmer(c("tea_freq_group","serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T"),
                                  c("ASM_z","ASMI_z", "grip_max_z","gait_speed_z"), #,"SPPB_total"
                                  data1_final, 
                                  Covariates_all_lmer, 
                                  2)
Lmer_results1_all_sig <- Lmer_results1_all[Lmer_results1_all$P_value < 0.05,]
table(Lmer_results1_all$Level)
table(data1_final$tea_freq_group)
#***************正文写作
a <- Lmer_results1_all_sig
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low,3)
a$CI_high <- round(a$CI_high,3)

a1 <- a[a$Predictor == "tea_freq_group",]
min(a1$Estimate)
max(a1$Estimate)
min(a1$CI_low)
max(a1$CI_high)

a2 <- a[a$Predictor != "tea_freq_group",]
min(a2$Estimate)
max(a2$Estimate)
min(a2$CI_low)
max(a2$CI_high)
####主分析可视化####
library(dplyr)
library(ggplot2)
library(stringr)
library(grid)
head(Lmer_results1_all)
plot_data_clean <- Lmer_results1_all %>%
  
  mutate(
    ## 1. 统一简化 Level 命名
    Level = case_when(
      Level == 1 ~ "Group 2",
      Level == 2 ~ "Group 3",
      TRUE ~ as.character(Level)
    )
  ) %>%
  
  mutate(
    ## significance
    sig_label = ifelse(P_value < 0.05, "*", ""),
    
    ## Predictor rename
    Predictor = recode(
      Predictor,
      "tea_freq_group"           = "Tea consumption",
      "serum_I_catechin_F0_T"    = "Catechin",
      "serum_I_epicatechin_F0_T" = "Epicatechin",
      "serum_I_EGC_F0_T"         = "Epigallocatechin",
      "serum_I_ECG_F0_T"         = "Epicatechin gallate",
      "serum_I_EGCG_F0_T"        = "Epigallocatechin gallate"
    ),
    
    ## 自动换行
    Predictor = str_wrap(Predictor, width = 16),
    
    ## Outcome rename
    Outcome = recode(
      Outcome,
      "ASM_z"        = "ASM",
      "ASMI_z"       = "ASMI",
      "grip_max_z"   = "Handgrip strength",
      "gait_speed_z" = "Walking speed"
    ),
    
    ## 2. 明确定义 Level 因子顺序（Group 3 在 Group 2 上方）
    Level = factor(Level, levels = c("Group 2", "Group 3")),
    
    ## Predictor 因子顺序
    Predictor = factor(
      Predictor,
      levels = c(
        "Tea consumption",
        "Catechin",
        "Epicatechin\ngallate",
        "Epigallocatechin\ngallate",
        "Epicatechin",
        "Epigallocatechin"

      )
    ),
    
    ## Outcome 因子顺序
    Outcome = factor(
      Outcome,
      levels = c(
        "ASM",
        "ASMI",
        "Handgrip strength",
        "Walking speed"
      )
    )
  )

p_final <- ggplot(
  plot_data_clean,
  aes(x = Estimate, y = Level)
) +
  
  ## 背景条
  geom_rect(
    aes(
      ymin = as.numeric(Level) - 0.35,
      ymax = as.numeric(Level) + 0.35,
      xmin = -Inf,
      xmax = Inf
    ),
    fill = "#F8F9FA",
    color = NA
  ) +
  
  ## 0线
  geom_vline(
    xintercept = 0,
    color = "grey70",
    linewidth = 0.5
  ) +
  
  ## CI 线性段
  geom_segment(
    aes(
      x = CI_low,
      xend = CI_high,
      yend = Level,
      color = Level
    ),
    linewidth = 1.5,
    alpha = 0.7
  ) +
  
  ## 点
  geom_point(
    aes(color = Level),
    size = 3
  ) +
  
  ## 中间白点
  geom_point(
    color = "white",
    size = 1.8
  ) +
  
  ## 显著性星号
  geom_text(
    aes(
      x = CI_high,
      label = sig_label
    ),
    hjust = -0.4,
    vjust = 0.7,
    size = 5,
    fontface = "bold"
  ) +
  
  ## 分面
  facet_grid(
    Predictor ~ Outcome,
    scales = "free_y"
  ) +
  
  ## 统一两组的配色
  scale_color_manual(
    values = c(
      "Group 2" = "#2A9D8F",  # 之前的青绿色
      "Group 3" = "#E76F51"   # 之前的橙红色
    )
  ) +
  
  ## labels
  labs(
    x = "β (95% CI)",
    y = NULL
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  )+
  
  ## theme
  theme_minimal(base_size = 14) +
  
  theme(
    text = element_text(color = "#212529"),
    
    strip.background = element_rect(
      fill = "#E9ECEF",
      color = "white",
      linewidth = 2
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 11,
      color = "#495057",
      lineheight = 0.95
    ),
    
    panel.border = element_rect(
      color = "#CED4DA",
      fill = NA,
      linewidth = 0.7,
      linetype = "dashed"
    ),
    
    ## 3. 修改 Y 轴标签：改为横向 (angle = 0) 并右对齐 (hjust = 1)
    axis.text.y = element_text(
      face = "bold",
      size = 11,
      angle = 0,
      hjust = 1
    ),
    
    axis.text.x = element_text(size = 10),
    
    axis.title.x = element_text(
      margin = ggplot2::margin(t = 15),
      size = 13,
      face = "bold"
    ),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.grid.major.x = element_line(
      color = "#DEE2E6",
      linetype = "dashed",
      linewidth = 0.3
    ),
    
    legend.position = "none",
    panel.spacing = unit(0.3, "lines"),
    
    ## 适当增加了左侧边距 (40)，确保横向的 "Group 2/3" 标签不会被切掉
    plot.margin = ggplot2::margin(20, 20, 20, 50)
  )

## 保存图片
ggsave(
  "Figure_2.pdf",
  p_final,
  width = 10.0,
  height = 10.0,
  device = cairo_pdf
)
####性别交互效应####
Lmer_results_interaction <- process_lmer_interaction(
  X = c(
    "tea_freq_group",
    "serum_I_catechin_F0_T",
    "serum_I_epicatechin_F0_T",
    "serum_I_EGC_F0_T",
    "serum_I_EGCG_F0_T",
    "serum_I_ECG_F0_T"
  ),
  Y = c(
    "ASM_z",
    "ASMI_z",
    "grip_max_z"
  ),
  data = data1_final,
  covariates = Covariates_all_lmer,
  interaction_var = "Sex_F0"
)
####验证性别交互####
x <- "tea_freq_group"
y  <- "ASMI_z"

formula_text <- paste(
  y, "~",
  x, "+",
  x, "*Sex_F0 +",
  paste(setdiff(Covariates_all_lmer,"Sex_F0"), collapse = " + "),
  "+ (1|ID)"
)
formula_full <- as.formula(formula_text)

model <- lmer(formula_full,data = data1_final)
summary(model)
####可视化#####
plot_data <- Lmer_results_interaction %>%
  
  mutate(
    
    ## 暴露名称
    Exposure = recode(
      Predictor,
      "tea_freq_group"           = "Tea consumption",
      "serum_I_catechin_F0_T"    = "Catechin",
      "serum_I_epicatechin_F0_T" = "Epicatechin",
      "serum_I_EGC_F0_T"         = "Epigallocatechin",
      "serum_I_ECG_F0_T"         = "Epicatechin gallate",
      "serum_I_EGCG_F0_T"        = "Epigallocatechin gallate"
    ),
    
    
    ## 组别比较
    Contrast = case_when(
      
      grepl("T1", Term_clean) ~
        "(Group 2 vs. 1)",
      
      grepl("T2", Term_clean) ~
        "(Group 3 vs. 1)",
      
      grepl("group1", Term_clean) ~
        "(Group 2 vs. 1)",
      
      grepl("group2", Term_clean) ~
        "(Group 3 vs. 1)"
    ),
    
    
    ## outcome名称
    Outcome = recode(
      Outcome,
      "ASM_z"  = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    ),
    
    
    ## 最终y轴标签
    Label = paste0(
      Exposure,
      " ",
      Contrast,
      " × Sex (Male vs. Female)"
    )
    
  )

plot_data <- plot_data %>%
  
  mutate(
    
    Exposure = factor(
      Exposure,
      levels = c(
        "Tea consumption",
        "Catechin",
        "Epicatechin",
        "Epigallocatechin",
        "Epigallocatechin gallate",
        "Epicatechin gallate"
      )
    ),
    
    Contrast = factor(
      Contrast,
      levels = c(
        "(Group 2 vs. 1)",
        "(Group 3 vs. 1)"
      )
    )
    
  ) %>%
  
  arrange(Exposure, Contrast)



# 方法：分别对每个 Outcome 排序，然后合并因子
plot_data <- plot_data %>%
  group_by(Outcome) %>%
  arrange(Outcome, Estimate) %>%        # 按 Outcome 分组，每组内按 Estimate 升序
  mutate(Label_order = paste0(Outcome, "_", row_number())) %>%  # 创建临时排序ID
  ungroup() %>%
  mutate(Label = factor(Label, levels = unique(Label)))  # 按当前顺序锁定因子

# 绘图
p <- ggplot(
  plot_data,
  aes(
    x = Estimate,
    y = Label
  )
) +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.15,
    linewidth = 0.8,
    color = "#444444"
  ) +
  geom_point(size = 3, color = "#222222") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(
    ~Outcome,
    ncol = 3
  ) +
  labs(
    x = "β (95%CI)",
    y = NULL,
    title = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "#F3F3F3"),
    strip.text = element_text(face = "bold", size = 13),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold", size = 13)
  )


ggsave(
  "p_interaction_sex.png",
  plot = p,
  width = 12,
  height = 12,
  dpi = 500
)
####年龄交互效应####
Lmer_results_interaction_age <- process_lmer_interaction(
  X = c(
    "tea_freq_group",
    "serum_I_catechin_F0_T",
    "serum_I_epicatechin_F0_T",
    "serum_I_EGC_F0_T",
    "serum_I_EGCG_F0_T",
    "serum_I_ECG_F0_T"
  ),
  Y = c(
    "ASM_z",
    "ASMI_z",
    "grip_max_z"
  ),
  data = data1_final,
  covariates = setdiff(Covariates_all_lmer,"Age"),
  interaction_var = "Followup"
)

a <- Lmer_results_interaction_age[Lmer_results_interaction_age$P_value < 0.05,]

writexl::write_xlsx(Lmer_results_interaction_age, "Lmer_results_interaction_age.xlsx")

####验证随访交互####
x <- "serum_I_EGC_F0_T"
y  <- "ASMI_z"

formula_text <- paste(
  y, "~",
  x, "+",
  x, "*Followup +",
  paste(setdiff(Covariates_all_lmer,c("Followup","Age")), collapse = " + "),
  "+ (1|ID)"
)
formula_full <- as.formula(formula_text)

model <- lmer(formula_full,data = data1_final)
summary(model)

####可视化#####

## 0. 加载 R 包

library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(patchwork)
library(tibble)
library(scales)


## 1. 基础设置

predictor_list <- c(
  "tea_freq_group",
  "serum_I_catechin_F0_T",
  "serum_I_epicatechin_F0_T",
  "serum_I_EGC_F0_T",
  "serum_I_EGCG_F0_T",
  "serum_I_ECG_F0_T"
)

outcome_model_list <- c(
  "ASM_z",
  "ASMI_z",
  "grip_max_z"
)

outcome_map_df <- tibble::tibble(
  outcome_model = c("ASM_z",    "ASMI_z",       "grip_max_z"),
  y_raw         = c("ASM",      "ASMI",         "grip_max"),
  y_label       = c("ASM (kg)", "ASMI (kg/m\u00B2)", "Handgrip strength (kg)")
)

predictor_rename <- c(
  "tea_freq_group"           = "Tea consumption",
  "serum_I_catechin_F0_T"    = "Catechin",
  "serum_I_epicatechin_F0_T" = "Epicatechin",
  "serum_I_EGC_F0_T"         = "Epigallocatechin",
  "serum_I_EGCG_F0_T"        = "Epigallocatechin gallate",
  "serum_I_ECG_F0_T"         = "Epicatechin gallate"
)


## 2. 数据预处理

data_lmer_plot <- data1_final

## 2.1 严格仅转换 ASM：如果 ASM 单位疑似为 g，则转为 kg
## ASMI 绝对不动
if ("ASM" %in% names(data_lmer_plot)) {
  if (mean(data_lmer_plot[["ASM"]], na.rm = TRUE) > 1000) {
    message("检测到 ASM 原始单位可能为 g，已自动将且仅将 ASM 转换为 kg；ASMI 保持原样。")
    data_lmer_plot$ASM <- data_lmer_plot$ASM / 1000
  }
}

## 2.2 统一 Followup 为 F1/F2/F3
normalize_followup_to_factor <- function(x) {
  
  x_chr <- as.character(x)
  x_non_na <- unique(x_chr[!is.na(x_chr)])
  
  if (all(x_non_na %in% c("0", "1", "2"))) {
    
    factor(
      x_chr,
      levels = c("0", "1", "2"),
      labels = c("F1", "F2", "F3")
    )
    
  } else if (all(x_non_na %in% c("1", "2", "3"))) {
    
    factor(
      x_chr,
      levels = c("1", "2", "3"),
      labels = c("F1", "F2", "F3")
    )
    
  } else if (all(x_non_na %in% c("F1", "F2", "F3"))) {
    
    factor(
      x_chr,
      levels = c("F1", "F2", "F3")
    )
    
  } else if (all(x_non_na %in% c("V1", "V2", "V3"))) {
    
    factor(
      x_chr,
      levels = c("V1", "V2", "V3")
    )
    
  } else {
    
    factor(x_chr)
  }
}

data_lmer_plot$Followup <- normalize_followup_to_factor(data_lmer_plot$Followup)

message("Followup 当前水平为：")
print(levels(data_lmer_plot$Followup))


## 2.3 把暴露变量统一设为分类变量
## 注意：这一步很重要，否则 0/1/2 会被 lmer 当成连续变量
for (v in predictor_list) {
  if (v %in% names(data_lmer_plot)) {
    data_lmer_plot[[v]] <- factor(
      as.character(data_lmer_plot[[v]]),
      levels = c("0", "1", "2")
    )
  }
}


## 3. LMM 交互项函数

process_lmer_interaction2 <- function(
    X,
    Y,
    data,
    covariates = NULL,
    interaction_var = "Sex_F0"
) {
  
  library(lme4)
  library(lmerTest)
  library(dplyr)
  
  safe_var <- function(v) paste0("`", v, "`")
  
  results_list <- list()
  counter <- 1
  
  for (x in X) {
    for (y in Y) {
      
      x_safe   <- safe_var(x)
      y_safe   <- safe_var(y)
      int_safe <- safe_var(interaction_var)
      
      cov_use <- setdiff(covariates, interaction_var)
      
      if (!is.null(cov_use) && length(cov_use) > 0) {
        
        cov_safe <- sapply(cov_use, safe_var)
        
        formula_text <- paste(
          y_safe, "~",
          paste0(x_safe, " * ", int_safe),
          "+",
          paste(cov_safe, collapse = " + "),
          "+ (1|ID)"
        )
        
      } else {
        
        formula_text <- paste(
          y_safe, "~",
          paste0(x_safe, " * ", int_safe),
          "+ (1|ID)"
        )
      }
      
      formula_full <- as.formula(formula_text)
      
      message("正在运行模型：", formula_text)
      
      lme_model <- tryCatch(
        suppressMessages(
          suppressWarnings(
            lmer(
              formula_full,
              data = data,
              REML = FALSE,
              control = lmerControl(
                optimizer = "bobyqa",
                optCtrl = list(maxfun = 1e5)
              )
            )
          )
        ),
        error = function(e) {
          message("模型失败：", formula_text)
          message("错误信息：", e$message)
          return(NULL)
        }
      )
      
      if (is.null(lme_model)) next
      
      is_singular <- isSingular(lme_model)
      
      sum_model <- summary(lme_model)
      
      coef_df <- as.data.frame(sum_model$coefficients)
      coef_df$Term <- rownames(coef_df)
      coef_df$Term_clean <- gsub("`", "", coef_df$Term)
      
      ## 提取：
      ## 1. x 的主效应，包括连续变量 x 或分类变量 x1/x2
      ## 2. x 与 interaction_var 的交互项
      is_main_x <- (
        coef_df$Term_clean == x |
          (
            startsWith(coef_df$Term_clean, x) &
              !grepl(":", coef_df$Term_clean) &
              coef_df$Term_clean != "(Intercept)"
          )
      )
      
      is_interaction_x <- (
        grepl(":", coef_df$Term_clean) &
          grepl(x, coef_df$Term_clean, fixed = TRUE) &
          grepl(interaction_var, coef_df$Term_clean, fixed = TRUE)
      )
      
      coef_x <- coef_df[is_main_x | is_interaction_x, , drop = FALSE]
      
      if (nrow(coef_x) == 0) next
      
      coef_x$Effect_Type <- ifelse(
        grepl(":", coef_x$Term_clean),
        "Interaction",
        "Main Effect"
      )
      
      variable_type <- if (is.factor(data[[x]])) {
        "Categorical"
      } else {
        "Continuous"
      }
      
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
      
      vc_df <- as.data.frame(VarCorr(lme_model))
      
      random_intercept_var <- vc_df$vcov[vc_df$grp == "ID"]
      residual_var <- vc_df$vcov[vc_df$grp == "Residual"]
      
      if (length(random_intercept_var) == 0) random_intercept_var <- NA_real_
      if (length(residual_var) == 0) residual_var <- NA_real_
      
      coef_x$Random_Intercept_Var <- random_intercept_var[1]
      coef_x$Residual_Var <- residual_var[1]
      
      results_list[[counter]] <- coef_x[, c(
        "Outcome",
        "Predictor",
        "Variable_Type",
        "Term_clean",
        "Effect_Type",
        "Estimate",
        "Std. Error",
        "CI_low",
        "CI_high",
        "P_value",
        "Significance",
        "Warning_Message",
        "Singular_Fit",
        "Random_Intercept_Var",
        "Residual_Var"
      )]
      
      counter <- counter + 1
    }
  }
  
  bind_rows(results_list)
}


## 4. 运行 LMM 交互项模型

Lmer_results_interaction_age <- process_lmer_interaction2(
  X = predictor_list,
  Y = outcome_model_list,
  data = data_lmer_plot,
  covariates = setdiff(Covariates_all_lmer, "Age"),
  interaction_var = "Followup"
)

## 检查交互项名称
message("提取到的交互项如下：")
print(
  Lmer_results_interaction_age %>%
    filter(Effect_Type == "Interaction") %>%
    dplyr::select(Predictor, Outcome, Term_clean, Estimate, P_value) %>%
    distinct()
)


## 5. 可视化辅助函数

get_legend_labels <- function(predictor) {
  
  if (predictor == "tea_freq_group") {
    
    c(
      "0" = "Non-drinker",
      "1" = "1-6 times/week",
      "2" = "\u22657 times/week"
    )
    
  } else if (predictor %in% c("serum_I_epicatechin_F0_T", "serum_I_EGC_F0_T")) {
    
    c(
      "0" = "Undetectable",
      "1" = "Low",
      "2" = "High"
    )
    
  } else {
    
    c(
      "0" = "Tertile 1",
      "1" = "Tertile 2",
      "2" = "Tertile 3"
    )
  }
}

get_single_level_name <- function(predictor, num_str) {
  
  num_str <- as.character(num_str)
  
  dplyr::case_when(
    predictor == "tea_freq_group" & num_str == "1" ~ "1-6/wk",
    predictor == "tea_freq_group" & num_str == "2" ~ "\u22657/wk",
    
    predictor %in% c("serum_I_epicatechin_F0_T", "serum_I_EGC_F0_T") & num_str == "1" ~ "Low",
    predictor %in% c("serum_I_epicatechin_F0_T", "serum_I_EGC_F0_T") & num_str == "2" ~ "High",
    
    num_str == "1" ~ "T2",
    num_str == "2" ~ "T3",
    
    TRUE ~ num_str
  )
}

format_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(p < 0.001, "<0.001", sprintf("=%.3f", p))
  )
}

format_beta <- function(beta) {
  ifelse(
    is.na(beta),
    "NA",
    sprintf("%.3f", beta)
  )
}

## 6. 修正版交互项标签函数

extract_predictor_level_from_term <- function(term, predictor) {
  
  parts <- strsplit(term, ":", fixed = TRUE)[[1]]
  pred_part <- parts[startsWith(parts, predictor)]
  
  if (length(pred_part) == 0) {
    return(NA_character_)
  }
  
  pred_part <- pred_part[1]
  pred_level <- substring(pred_part, nchar(predictor) + 1)
  
  if (pred_level == "") {
    return(NA_character_)
  }
  
  pred_level
}

extract_followup_from_term <- function(term, interaction_var = "Followup") {
  
  parts <- strsplit(term, ":", fixed = TRUE)[[1]]
  follow_part <- parts[startsWith(parts, interaction_var)]
  
  if (length(follow_part) == 0) {
    return(NA_character_)
  }
  
  follow_part <- follow_part[1]
  follow_raw <- substring(follow_part, nchar(interaction_var) + 1)
  
  if (follow_raw == "") {
    return(NA_character_)
  }
  
  follow_raw
}

clean_followup_name <- function(follow_raw) {
  
  if (is.na(follow_raw)) {
    return(NA_character_)
  }
  
  ## 如果已经是 F2/F3，直接返回
  if (str_detect(follow_raw, "^F[0-9]+$")) {
    return(follow_raw)
  }
  
  ## 如果是 V2/V3，也直接返回
  if (str_detect(follow_raw, "^V[0-9]+$")) {
    return(follow_raw)
  }
  
  ## 如果是数字，兼容旧写法：
  ## 例如 Followup1 / Followup2 通常来自原始 0/1/2 编码，对应 F2/F3
  if (str_detect(follow_raw, "^[0-9]+$")) {
    return(paste0("F", as.numeric(follow_raw) + 1))
  }
  
  follow_raw
}

make_interaction_label <- function(
    result_df,
    predictor,
    outcome_model,
    interaction_var = "Followup"
) {
  
  df <- result_df %>%
    filter(
      Predictor == predictor,
      Outcome == outcome_model,
      Effect_Type == "Interaction"
    )
  
  if (nrow(df) == 0) {
    return("No interaction term")
  }
  
  df_labels <- df %>%
    rowwise() %>%
    mutate(
      pred_level = extract_predictor_level_from_term(Term_clean, predictor),
      followup_raw = extract_followup_from_term(Term_clean, interaction_var),
      time_name = clean_followup_name(followup_raw),
      group_name = get_single_level_name(predictor, pred_level),
      label_one = paste0(
        time_name,
        " \u00D7 ",
        group_name,
        ": \u03B2=",
        format_beta(Estimate),
        ", P",
        format_p(P_value)
      )
    ) %>%
    ungroup()
  
  paste(df_labels$label_one, collapse = "\n")
}


## 7. 单图绘制函数

plot_raw_longitudinal_interaction <- function(
    data,
    predictor,
    y_raw,
    outcome_model,
    result_df,
    followup_var = "Followup",
    interaction_var = "Followup",
    y_label = NULL,
    title = NULL,
    show_x_text = TRUE
) {
  
  plot_data <- data %>%
    filter(
      !is.na(.data[[predictor]]),
      !is.na(.data[[y_raw]]),
      !is.na(.data[[followup_var]])
    ) %>%
    mutate(
      Followup_plot = factor(
        .data[[followup_var]],
        levels = levels(data[[followup_var]])
      ),
      Group_plot = factor(
        .data[[predictor]],
        levels = c("0", "1", "2")
      )
    )
  
  if (nrow(plot_data) == 0) {
    stop(paste0("没有可用于绘图的数据：", predictor, " - ", y_raw))
  }
  
  label_text <- make_interaction_label(
    result_df = result_df,
    predictor = predictor,
    outcome_model = outcome_model,
    interaction_var = interaction_var
  )
  
  pretty_pred_name <- unname(predictor_rename[predictor])
  if (is.na(pretty_pred_name)) pretty_pred_name <- predictor
  
  current_legend_labels <- get_legend_labels(predictor)
  
  color_values <- c(
    "0" = "#2F4F4F",
    "1" = "#D95F02",
    "2" = "#1B7C3D"
  )
  
  group_levels_present <- levels(droplevels(plot_data$Group_plot))
  
  current_legend_labels <- current_legend_labels[
    names(current_legend_labels) %in% group_levels_present
  ]
  
  current_color_values <- color_values[
    names(color_values) %in% group_levels_present
  ]
  
  x_total_levels <- length(levels(plot_data$Followup_plot))
  
  x_position <- if (y_raw == "grip_max") {
    x_total_levels - 1.85
  } else {
    x_total_levels - 1.85
  }
  
  y_range <- range(plot_data[[y_raw]], na.rm = TRUE)
  y_diff <- diff(y_range)
  
  if (is.na(y_diff) || y_diff == 0) {
    y_position <- y_range[2]
  } else {
    y_position <- y_range[2] - 0.02 * y_diff
  }
  
  p <- ggplot(
    plot_data,
    aes(
      x = Followup_plot,
      y = .data[[y_raw]],
      color = Group_plot,
      group = Group_plot
    )
  ) +
    geom_jitter(
      width = 0.08,
      height = 0,
      alpha = 0.16,
      size = 0.9,
      stroke = 0
    ) +
    stat_summary(
      fun = mean,
      geom = "line",
      linewidth = 1.6
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 4.0,
      shape = 21,
      fill = "white",
      stroke = 1.8
    ) +
    annotate(
      "label",
      x = x_position,
      y = y_position,
      label = label_text,
      hjust = 0,
      vjust = 1,
      size = 3.6,
      fontface = "bold.italic",
      fill = scales::alpha("white", 0.88),
      label.size = NA,
      lineheight = 1.0
    ) +
    labs(
      title = title,
      x = NULL,
      y = y_label,
      color = pretty_pred_name
    ) +
    scale_color_manual(
      values = current_color_values,
      breaks = names(current_legend_labels),
      labels = current_legend_labels,
      drop = TRUE
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 15,
        margin = ggplot2::margin(b = 5)
      ),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(color = "black", size = 12),
      axis.text.x = if (show_x_text) {
        element_text(color = "black", size = 13, face = "bold")
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x_text) {
        element_line()
      } else {
        element_blank()
      },
      axis.line = element_line(color = "gray40", linewidth = 0.6),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 12)
    )
  
  return(p)
}


## 8. 每个暴露变量生成一行图

plot_row_per_predictor <- function(
    predictor,
    data_input,
    result_input,
    outcome_map_input,
    is_last_row = FALSE
) {
  
  row_plots <- purrr::pmap(
    outcome_map_input,
    function(outcome_model, y_raw, y_label) {
      
      plot_raw_longitudinal_interaction(
        data = data_input,
        predictor = predictor,
        y_raw = y_raw,
        outcome_model = outcome_model,
        result_df = result_input,
        y_label = y_label,
        title = NULL,
        show_x_text = is_last_row
      )
    }
  )
  
  row_bundle <- patchwork::wrap_plots(row_plots, ncol = 3) +
    patchwork::plot_layout(guides = "collect")
  
  return(row_bundle)
}


## 9. 生成 6 × 3 总图

message("正在构建 6 × 3 纵向交互项大图...")

all_rows <- list()

for (i in seq_along(predictor_list)) {
  
  pred <- predictor_list[i]
  is_last <- i == length(predictor_list)
  
  all_rows[[pred]] <- plot_row_per_predictor(
    predictor = pred,
    data_input = data_lmer_plot,
    result_input = Lmer_results_interaction_age,
    outcome_map_input = outcome_map_df,
    is_last_row = is_last
  )
}

p_mega_final <- patchwork::wrap_plots(all_rows, ncol = 1) +
  patchwork::plot_annotation(
    title = NULL,
    caption = NULL,
    theme = theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 22,
        margin = ggplot2::margin(b = 15)
      ),
      plot.caption = element_text(
        size = 12,
        face = "italic",
        color = "gray30"
      )
    )
  )


## 10. 保存图片

if (requireNamespace("ragg", quietly = TRUE)) {
  
  ggsave(
    filename = "mega_combined_longitudinal_interaction.png",
    plot = p_mega_final,
    width = 16,
    height = 24,
    dpi = 500,
    device = ragg::agg_png
  )
  
} else {
  ggsave(
    filename = "mega_combined_longitudinal_interaction.png",
    plot = p_mega_final,
    width = 16,
    height = 24,
    dpi = 500
  )
}
####模型性能提高####
Lmer_model_compare_Tea <- compare_lmer_predictors_by_pair(
  pair_df = Lmer_results1_all_sig,
  data = data1_final,
  covariates = c(Covariates_all_lmer),
  predictor_col = "Predictor",
  outcome_col = "Outcome",
  predictor_type_col = "Variable_Type",
  id_var = "ID",
  predictor_rename = predictor_rename ,
  outcome_rename = outcome_rename,
  output_file_raw = "Tea_ASM_model_performance_raw.csv",
  output_file_format = "Tea_ASM_model_performance_format.csv"
)

Lmer_model_compare_Tea_output <- Lmer_model_compare_Tea$format[,c("Outcome", "Predictor","Delta_Marginal_R2_percent","LRT_P_value")]

####data2-简单纵向频率分类&骨骼肌####
data_list <- list(Cov_F0_final, Tea_summary, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], ASM_Height_long)
data2_all <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
# 合并握力等，要用left_join
data2_v2 <- left_join(data2_all, SPPB_grip_long, by = c("ID","Followup"))

data2_final <- scale_columns_group(data2_v2, c("ASM","ASMI","grip_max"), "Sex_F0")
a <- unique(data2_final$ID)
####lmer####
Lmer_results2_all <- process_lmer(c("tea_freq_group_summary"),
                                  c("ASM_z","ASMI_z","grip_max_z"), #,"SPPB_total"
                                  data2_final, 
                                  Covariates_all_lmer, 
                                  2)

####data3-饮茶基线&骨骼肌####
data_list <- list(Cov_F0_final, Tea_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], ASM_Height_long)
data3_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))

# 合并握力等，要用left_join
data3_v2 <- left_join(data3_v1, SPPB_grip_long, by = c("ID","Followup"))
data3_v3 <- left_join(data3_v2, Catechin, by = "ID")

data3_v3 <- group_variables3(
  data = data3_v3,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(data3_v3)%in% factor_name)
for(i in idx ){
  data3_v3[[i]] <-  as.factor(data3_v3[[i]])
}

data3_final <- scale_columns_group(data3_v3, c("ASM","ASMI","grip_max"), "Sex_F0")
data3_final <- scale_columns(data3_final, c("gait_speed"))

a <- unique(data3_final$ID)
table(data3_final$Followup)
####lmer####
Lmer_results3_all <- process_lmer(c("tea_freq_group","serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T"),
                                   c("ASM_z","ASMI_z", "grip_max_z","gait_speed_z"), #,"SPPB_total"
                                   data3_final, 
                                   Covariates_all_lmer, 
                                   2)
Lmer_results3_all_sig <- Lmer_results3_all[Lmer_results3_all$P_value < 0.05,]
####data4-饮茶轨迹&骨骼肌####
data_list <- list(Cov_F0_final, Tea_tra_all, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], ASM_Height_long)
data4_all <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
# 合并握力等，要用left_join
data4_v2 <- left_join(data4_all, SPPB_grip_long, by = c("ID","Followup"))

data4_final <- scale_columns_group(data4_v2, c("ASM","ASMI","grip_max"), "Sex_F0")

a <- unique(data4_final$ID)
####lmer####
Lmer_results4_all <- process_lmer(c("class"),
                                  c("ASM_z","ASMI_z","grip_max_z"), #,"SPPB_total"
                                  data4_final, 
                                  Covariates_all_lmer, 
                                  2)
Lmer_results4_all_sig <- Lmer_results4_all[Lmer_results4_all$P_value < 0.05,]
####data5-平均饮茶频率&骨骼肌####
data_list <- list(Cov_F0_final, Tea_mean_freq, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], ASM_Height_long)
data5_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))

# 合并握力等，要用left_join
data5_v2 <- left_join(data5_v1, SPPB_grip_long, by = c("ID","Followup"))

data5_final <- scale_columns_group(data5_v2, c("ASM","ASMI","grip_max"), "Sex_F0")

a <- unique(data5_final$ID)
####lmer####
Lmer_results5_all <- process_lmer(c("tea_freq_group"),
                                   c("ASM_z","ASMI_z","grip_max_z"), #,"SPPB_total"
                                  data5_final, 
                                  Covariates_all_lmer, 
                                   2)
Lmer_results5_all_sig <- Lmer_results5_all[Lmer_results5_all$P_value < 0.05,]
####查看各个随访的喝茶种类情况####
merge_tea_types <- function(data, tea_col, filter_col = NULL, filter_value = 1) {
  
  # 如果有筛选列，先筛选饮茶者
  if (!is.null(filter_col)) {
    data <- data %>%
      filter(!!sym(filter_col) == filter_value)
  }
  
  # 处理茶种类
  result <- data %>%
    dplyr::select(tea_type = all_of(tea_col)) %>%
    mutate(ID = row_number()) %>%
    separate_rows(tea_type, sep = "[；，、]") %>%
    mutate(
      tea_type = str_trim(tea_type),
      
      # 统一合并规则（所有随访轮次一致）
      tea_type = case_when(
        
        # 1. 普洱类（所有变体）
        str_detect(tea_type, "普|耳|饵|菊普|生普|功夫茶普洱|普洱花茶|普洱菊花|普洱龙井|普洱毛尖|普洱学菊|普洱玫瑰花|普洱绞股蓝|普洱单木丛|普洱，|普洱/|普洱田七|小青柑|小金柑|白茶普洱|花茶普洱|普洱苦丁") ~ "普洱",
        
        # 2. 乌龙茶类（不包括单丛，单丛单独处理）
        str_detect(tea_type, "乌龙茶|乌龙$|大红袍|凤凰茶|清茶|水仙") ~ "乌龙茶",
        
        # 3. 单丛类（单独保留）
        str_detect(tea_type, "单丛|单松|单枞|枞茶|凤凰|单") ~ "单丛",
        
        # 4. 绿茶类
        str_detect(tea_type, "绿茶|龙井") ~ "绿茶",
        
        # 5. 红茶类
        str_detect(tea_type, "^红茶$") ~ "红茶",
        
        # 6. 花茶类
        str_detect(tea_type, "花|花茶|菊花|茉莉|香花茶|花茶芦荟菊|菊花茶寿眉|昆仑雪菊|玫瑰花|枸杞红枣|陈皮茶") ~ "花茶",
        
        # 7. 黑茶类（藏茶归入黑茶）
        str_detect(tea_type, "^黑茶$|藏茶") ~ "黑茶",
        
        # 8. 苦丁茶类
        str_detect(tea_type, "苦丁茶|苦$") ~ "苦丁茶",
        
        # 9. 罗布麻类
        str_detect(tea_type, "罗布麻") ~ "罗布麻茶",
        
        # 10. 绞股蓝类
        str_detect(tea_type, "绞股蓝|绞胶蓝|搞股兰") ~ "绞股蓝茶",
        
        # 11. 丹参类
        str_detect(tea_type, "丹参|丹心|丹参保心|天草丹参") ~ "丹参茶",
        
        # 12. 灵芝类
        str_detect(tea_type, "灵芝") ~ "灵芝茶",
        
        # 13. 谷物茶类（大麦、荞麦、苦荞麦）
        str_detect(tea_type, "大麦茶|荞麦茶|苦荞麦|苦麦茶") ~ "谷物茶",
        
        # 14. 保健茶类
        str_detect(tea_type, "保健茶|养生茶|降压茶") ~ "保健茶",
        
        # 15. 牛蒡茶
        str_detect(tea_type, "牛蒡|牛膀") ~ "牛蒡茶",
        
        # 16. 其他保留原样（清理括号）
        TRUE ~ str_replace_all(tea_type, "（.*", "")
      )
    ) %>%
    filter(
      !is.na(tea_type),
      tea_type != "",
      !str_detect(tea_type, "^[0-9]+$"),      # 剔除纯数字
      !tea_type %in% c("4", "5", "NA", "NULL", "无", "不清楚")  # 剔除无效值
    ) %>%
    count(tea_type) %>%
    mutate(
      Percent = 100 * n / sum(n),
      Percent_label = sprintf("%.2f%%", Percent)
    ) %>%
    arrange(desc(Percent))
  
  return(result)
}

# F0
tea_type_F0 <- merge_tea_types(
  data = Tea_F0_filter,
  tea_col = "常饮用茶种类_F0",
  filter_col = "是否喝茶_F0",
  filter_value = 1
)

# F1
tea_type_F1 <- merge_tea_types(
  data = Tea_F1_filter,
  tea_col = "常饮用茶种类_F1",
  filter_col = "是否喝茶_F1",
  filter_value = 1
)

# F2（注意筛选列名不同）
tea_type_F2 <- merge_tea_types(
  data = Tea_F2_filter,
  tea_col = "常饮用茶种类_F2",
  filter_col = "过去一年是否经常喝茶_F2",
  filter_value = 1
)

# F3
tea_type_F3 <- merge_tea_types(
  data = Tea_F3_filter,
  tea_col = "常饮用茶种类_F3",
  filter_col = "是否喝茶_F3",
  filter_value = 1
)
####对茶的其他种类的统计描述####
# 1. 茶种类合并函数

merge_tea_type_one <- function(tea) {
  
  tea <- stringr::str_trim(as.character(tea))
  
  if (is.na(tea) || tea == "") {
    return(NA_character_)
  }
  
  merged_tea <- dplyr::case_when(
    
    # 普洱类
    stringr::str_detect(
      tea,
      "普|耳|饵|菊普|生普|功夫茶普洱|普洱花茶|普洱菊花|普洱龙井|普洱毛尖|普洱学菊|普洱玫瑰花|普洱绞股蓝|普洱单木丛|普洱，|普洱/|普洱田七|小青柑|小金柑|白茶普洱|花茶普洱|普洱苦丁|善洱"
    ) ~ "普洱",
    
    # 乌龙茶类
    stringr::str_detect(
      tea,
      "乌龙茶|乌龙$|大红袍|凤凰茶|清茶|水仙|功夫茶"
    ) ~ "乌龙茶",
    
    # 单丛类
    stringr::str_detect(
      tea,
      "单丛|单松|单枞|枞茶|凤凰|单"
    ) ~ "单丛",
    
    # 绿茶类
    stringr::str_detect(
      tea,
      "绿茶|龙井|毛尖"
    ) ~ "绿茶",
    
    # 红茶类
    stringr::str_detect(
      tea,
      "^红茶$"
    ) ~ "红茶",
    
    # 花茶类
    stringr::str_detect(
      tea,
      "花|花茶|菊花|茉莉|香花茶|花茶芦荟菊|菊花茶寿眉|昆仑雪菊|玫瑰花|枸杞红枣|陈皮茶|菊|茉|罗汉果茶|猴王牌茶叶|六合茶|山楂叶|洛茶|桑叶|柠檬茶|宝梨"
    ) ~ "花茶",
    
    # 黑茶类
    stringr::str_detect(
      tea,
      "^黑茶$|藏茶"
    ) ~ "黑茶",
    
    # 苦丁茶类
    stringr::str_detect(
      tea,
      "苦丁茶|苦$"
    ) ~ "苦丁茶",
    
    # 罗布麻类
    stringr::str_detect(
      tea,
      "罗布麻"
    ) ~ "罗布麻茶",
    
    # 绞股蓝类
    stringr::str_detect(
      tea,
      "绞股蓝|绞胶蓝|搞股兰"
    ) ~ "绞股蓝茶",
    
    # 丹参类
    stringr::str_detect(
      tea,
      "丹参|丹心|丹参保心|天草丹参"
    ) ~ "丹参茶",
    
    # 灵芝类
    stringr::str_detect(
      tea,
      "灵芝"
    ) ~ "灵芝茶",
    
    # 谷物茶类
    stringr::str_detect(
      tea,
      "大麦茶|荞麦茶|苦荞麦|苦麦茶"
    ) ~ "谷物茶",
    
    # 保健茶类
    stringr::str_detect(
      tea,
      "保健茶|养生茶|降压茶|中药茶|舌斛茶|自制茶叶|祛湿|痛风茶|熊胆|溪黄草|东革阿里|丁香母藤茶"
    ) ~ "保健茶",
    
    # 牛蒡茶
    stringr::str_detect(
      tea,
      "牛蒡|牛膀"
    ) ~ "牛蒡茶",
    
    TRUE ~ stringr::str_replace_all(tea, "（.*", "")
  )
  
  return(merged_tea)
}


# 2. 茶种类英文翻译函数

translate_tea_type <- function(x) {
  
  dplyr::case_when(
    x == "乌龙茶"     ~ "Oolong tea",
    x == "决明子"     ~ "Cassia seed tea",
    x == "减肥茶"     ~ "Slimming tea",
    x == "减肥茶（"   ~ "Slimming tea",
    x == "单丛"       ~ "Dancong tea",
    x == "大麦茶"     ~ "Barley tea",
    x == "谷物茶"     ~ "Grain tea",
    x == "山楂茶"     ~ "Hawthorn tea",
    x == "普洱"       ~ "Pu-erh tea",
    x == "红茶"       ~ "Black tea",
    x == "绿茶"       ~ "Green tea",
    x == "罗布麻茶"   ~ "Luobuma tea",
    x == "花茶"       ~ "Flower tea",
    x == "苦丁茶"     ~ "Kuding tea",
    x == "藏茶"       ~ "Tibetan tea",
    x == "黑茶"       ~ "Dark tea",
    x == "银杏叶"     ~ "Ginkgo leaf tea",
    x == "绞股蓝茶"   ~ "Gynostemma tea",
    x == "丹参茶"     ~ "Danshen tea",
    x == "灵芝茶"     ~ "Ganoderma tea",
    x == "保健茶"     ~ "Health tea",
    x == "牛蒡茶"     ~ "Burdock tea",
    TRUE ~ as.character(x)
  )
}


# 3. 主函数：只返回 list，不自动导出

plot_tea_description_list <- function(
    data,
    
    # 茶相关变量名
    drink_col,
    type_col,
    amount_col,
    strength_col,
    
    # 每周冲茶次数变量，可选
    freq_col = NULL,
    
    # 喝茶年限变量，可选
    duration_col = NULL,
    
    # 参数
    drink_yes_value = 1,
    amount_unit_label = "Tea amount (jin/year)",
    freq_unit_label = "Tea consumption frequency (times/week)",
    duration_unit_label = "Duration (years)",
    strength_levels = c(1, 2, 3, 4, 5),
    strength_labels = c(
      "Strongest",
      "Strong",
      "Moderate",
      "Weak",
      "Weakest"
    ),
    
    # 图标题前缀
    visit_label = "F0"
) {
  
  # 0. 检查变量是否存在
  
  required_cols <- c(
    drink_col,
    type_col,
    amount_col,
    strength_col
  )
  
  if (!is.null(freq_col)) {
    required_cols <- c(required_cols, freq_col)
  }
  
  if (!is.null(duration_col)) {
    required_cols <- c(required_cols, duration_col)
  }
  
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "以下变量在 data 中不存在：",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  # 1. 仅保留饮茶者，并剔除茶种类为纯数字的行
  
  tea_drinkers <- data %>%
    mutate(
      .drink = suppressWarnings(as.numeric(.data[[drink_col]])),
      .tea_type_for_filter = stringr::str_trim(as.character(.data[[type_col]]))
    ) %>%
    filter(
      .drink == drink_yes_value
    ) %>%
    filter(
      !(
        !is.na(.tea_type_for_filter) &
          stringr::str_detect(.tea_type_for_filter, "^[0-9]+$")
      )
    )
  
  # 2. A. 茶种类分布
  
  tea_type_long <- tea_drinkers %>%
    dplyr::select(
      tea_type_raw = all_of(type_col)
    ) %>%
    mutate(
      ID_tmp = row_number(),
      tea_type_raw = as.character(tea_type_raw)
    ) %>%
    tidyr::separate_rows(
      tea_type_raw,
      sep = "[；，、,/]"
    ) %>%
    mutate(
      tea_type_raw = stringr::str_trim(tea_type_raw),
      tea_type_raw = dplyr::na_if(tea_type_raw, ""),
      tea_type_raw = dplyr::na_if(tea_type_raw, "NA"),
      tea_type_raw = dplyr::na_if(tea_type_raw, "无"),
      tea_type_raw = dplyr::na_if(tea_type_raw, "没有")
    ) %>%
    filter(
      !is.na(tea_type_raw),
      tea_type_raw != "",
      !stringr::str_detect(tea_type_raw, "^[0-9]+$"),
      !stringr::str_detect(tea_type_raw, "^[0-9]+\\.[0-9]+$")
    ) %>%
    mutate(
      tea_type_merged = sapply(tea_type_raw, merge_tea_type_one)
    ) %>%
    filter(
      !is.na(tea_type_merged),
      tea_type_merged != "",
      !stringr::str_detect(tea_type_merged, "^[0-9]+$"),
      !stringr::str_detect(tea_type_merged, "^[0-9]+\\.[0-9]+$")
    )
  
  tea_type_summary <- tea_type_long %>%
    count(tea_type_merged, name = "n") %>%
    mutate(
      Percent = 100 * n / sum(n),
      Tea_type_EN = translate_tea_type(tea_type_merged),
      Percent_label = sprintf("%.2f%%", Percent),
      Tea_type_EN = forcats::fct_reorder(Tea_type_EN, Percent)
    ) %>%
    arrange(desc(Percent))
  
  p_tea_type <- ggplot(
    tea_type_summary,
    aes(
      x = Percent,
      y = Tea_type_EN
    )
  ) +
    geom_col(
      fill = "#238B45",
      width = 0.75
    ) +
    geom_text(
      aes(label = Percent_label),
      hjust = -0.15,
      size = 4
    ) +
    expand_limits(
      x = max(tea_type_summary$Percent, na.rm = TRUE) * 1.18
    ) +
    labs(
      title = paste0(visit_label, ". Tea type"),
      x = "Percentage (%)",
      y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold")
    )
  
  # 3. B. 每周冲茶次数柱状图
  #    每周频次 >= 21 的全部合并为一组，横坐标显示为 >=21
  
  freq_data <- NULL
  freq_summary <- NULL
  p_weekly_frequency <- NULL
  
  if (!is.null(freq_col)) {
    
    freq_data <- tea_drinkers %>%
      transmute(
        tea_freq_raw = suppressWarnings(
          as.numeric(.data[[freq_col]])
        )
      ) %>%
      filter(
        !is.na(tea_freq_raw),
        tea_freq_raw > 0
      ) %>%
      mutate(
        tea_freq_group = dplyr::case_when(
          tea_freq_raw >= 21 ~ 21,
          TRUE ~ tea_freq_raw
        ),
        tea_freq_label = dplyr::case_when(
          tea_freq_group == 21 ~ "≥21",
          TRUE ~ as.character(tea_freq_group)
        )
      )
    
    if (nrow(freq_data) > 0) {
      
      freq_levels <- freq_data %>%
        distinct(
          tea_freq_group,
          tea_freq_label
        ) %>%
        arrange(tea_freq_group) %>%
        pull(tea_freq_label)
      
      freq_summary <- freq_data %>%
        count(
          tea_freq_group,
          tea_freq_label,
          name = "n"
        ) %>%
        mutate(
          Percent = 100 * n / sum(n),
          Percent_label = sprintf("%.1f%%", Percent),
          tea_freq_factor = factor(
            tea_freq_label,
            levels = freq_levels
          )
        ) %>%
        arrange(tea_freq_group)
      
      p_weekly_frequency <- ggplot(
        freq_summary,
        aes(
          x = tea_freq_factor,
          y = n
        )
      ) +
        geom_col(
          fill = "#2A9D8F",
          width = 0.75
        ) +
        expand_limits(
          y = max(freq_summary$n, na.rm = TRUE) * 1.15
        ) +
        labs(
          title = paste0(visit_label, ". Tea frequency"),
          x = "Tea consumption frequency (times/week)",
          y = "Count"
        ) +
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          axis.text.x = element_text(
            angle = 0,
            hjust = 0.5
          )
        )
    }
  }
  # 4. C. 茶叶用量分布
  
  amount_data <- tea_drinkers %>%
    transmute(
      tea_amount = suppressWarnings(
        as.numeric(.data[[amount_col]])
      )
    ) %>%
    filter(
      !is.na(tea_amount),
      tea_amount > 0
    )
  
  p_tea_amount <- NULL
  
  if (nrow(amount_data) > 2) {
    
    min_amount <- min(
      amount_data$tea_amount,
      na.rm = TRUE
    )
    
    p_tea_amount <- ggplot(
      amount_data,
      aes(x = tea_amount)
    ) +
      geom_density(
        linewidth = 1.2,
        color = "#E76F51",
        fill = "#E76F51",
        alpha = 0.25,
        trim = TRUE
      ) +
      coord_cartesian(
        xlim = c(min_amount, NA)
      ) +
      scale_x_continuous(
        expand = expansion(mult = c(0, 0.02))
      ) +
      labs(
        title = paste0(visit_label, ". Tea amount"),
        x = amount_unit_label,
        y = "Density"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold")
      )
  }
  
  # 5. D. 喝茶浓淡程度
  
  strength_data <- tea_drinkers %>%
    transmute(
      tea_strength = suppressWarnings(
        as.numeric(.data[[strength_col]])
      )
    ) %>%
    filter(
      !is.na(tea_strength)
    ) %>%
    mutate(
      tea_strength_factor = factor(
        tea_strength,
        levels = strength_levels,
        labels = strength_labels
      )
    )
  
  p_brewing_strength <- ggplot(
    strength_data,
    aes(x = tea_strength_factor)
  ) +
    geom_bar(
      fill = "#457B9D"
    ) +
    labs(
      title = paste0(visit_label, ". Brewing strength"),
      x = NULL,
      y = "Count"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(
        angle = 25,
        hjust = 1
      )
    )
  
  # 6. E. 喝茶年限，可选
  
  p_tea_duration <- NULL
  duration_data <- NULL
  
  if (!is.null(duration_col)) {
    
    duration_data <- tea_drinkers %>%
      transmute(
        tea_duration = suppressWarnings(
          as.numeric(.data[[duration_col]])
        )
      ) %>%
      filter(
        !is.na(tea_duration),
        tea_duration > 0
      )
    
    if (nrow(duration_data) > 2) {
      
      min_duration <- min(
        duration_data$tea_duration,
        na.rm = TRUE
      )
      
      p_tea_duration <- ggplot(
        duration_data,
        aes(x = tea_duration)
      ) +
        geom_density(
          linewidth = 1.2,
          color = "#6A4C93",
          fill = "#6A4C93",
          alpha = 0.25,
          trim = TRUE
        ) +
        coord_cartesian(
          xlim = c(min_duration, NA)
        ) +
        scale_x_continuous(
          expand = expansion(mult = c(0, 0.02))
        ) +
        labs(
          title = paste0(visit_label, ". Tea drinking duration"),
          x = duration_unit_label,
          y = "Density"
        ) +
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold")
        )
    }
  }
  
  # 7. 只返回 list，不自动导出
  
  result <- list(
    
    data = list(
      tea_drinkers = tea_drinkers,
      tea_type_long = tea_type_long,
      tea_type_summary = tea_type_summary,
      freq_data = freq_data,
      freq_summary = freq_summary,
      amount_data = amount_data,
      strength_data = strength_data,
      duration_data = duration_data
    ),
    
    plots = list(
      tea_type = p_tea_type,
      weekly_frequency = p_weekly_frequency,
      tea_amount = p_tea_amount,
      brewing_strength = p_brewing_strength,
      tea_duration = p_tea_duration
    )
  )
  
  return(result)
}

res_F0 <- plot_tea_description_list(
  data = Tea_F0_filter,
  drink_col = "是否喝茶_F0",
  type_col = "常饮用茶种类_F0",
  amount_col = "过去一年茶叶总量_斤_F0",
  strength_col = "喝茶浓淡_F0",
  freq_col = "每周冲茶次数_F0",
  duration_col = "累计喝茶总年限_F0",
  visit_label = "F0"
)

res_F1 <- plot_tea_description_list(
  data = Tea_F1_filter,
  drink_col = "是否喝茶_F1",
  type_col = "常饮用茶种类_F1",
  amount_col = "过去一年茶叶总量_斤_F1",
  strength_col = "喝茶浓淡_F1",
  freq_col = "每周冲茶次数_F1",
  duration_col = NULL,
  visit_label = "F1"
)


res_F2 <- plot_tea_description_list(
  data = Tea_F2_filter,
  drink_col = "过去一年是否经常喝茶_F2",
  type_col = "常饮用茶种类_F2",
  amount_col = "过去一年茶叶总量_斤_F2",
  strength_col = "喝茶浓淡程度_F2",
  freq_col = "每周泡茶次数_F2",
  duration_col = NULL,
  visit_label = "F2"
)

res_F3 <- plot_tea_description_list(
  data = Tea_F3_filter,
  drink_col = "是否喝茶_F3",
  type_col = "常饮用茶种类_F3" ,
  amount_col = "过去一年茶叶总量_斤_F3",
  strength_col = "喝茶浓淡_F3",
  freq_col = "每周冲茶次数_F3",
  duration_col = NULL,
  visit_label = "F3"
)

####全面考虑茶的多方面#####
# 茶种类评分函数（整合合并规则）

get_tea_type_score <- function(x) {
  
  # 不喝茶
  if (is.na(x) || x == "") {
    return(0)
  }
  
  # 拆分多种茶（支持 ； ， 、）
  teas <- unlist(str_split(x, "[；，、]"))
  teas <- str_trim(teas)
  
  # 统一合并规则（与 merge_tea_types 完全一致）
  merged_teas <- character(length(teas))
  
  for (i in seq_along(teas)) {
    tea <- teas[i]
    
    merged_tea <- case_when(
      
      # 普洱类
      stringr::str_detect(
        tea,
        "普|耳|饵|菊普|生普|功夫茶普洱|普洱花茶|普洱菊花|普洱龙井|普洱毛尖|普洱学菊|普洱玫瑰花|普洱绞股蓝|普洱单木丛|普洱，|普洱/|普洱田七|小青柑|小金柑|白茶普洱|花茶普洱|普洱苦丁|善洱"
      ) ~ "普洱",
      
      # 乌龙茶类
      stringr::str_detect(
        tea,
        "乌龙茶|乌龙$|大红袍|凤凰茶|清茶|水仙|功夫茶"
      ) ~ "乌龙茶",
      
      # 单丛类
      stringr::str_detect(
        tea,
        "单丛|单松|单枞|枞茶|凤凰|单"
      ) ~ "单丛",
      
      # 绿茶类
      stringr::str_detect(
        tea,
        "绿茶|龙井|毛尖"
      ) ~ "绿茶",
      
      # 红茶类
      stringr::str_detect(
        tea,
        "^红茶$"
      ) ~ "红茶",
      
      # 花茶类
      stringr::str_detect(
        tea,
        "花|花茶|菊花|茉莉|香花茶|花茶芦荟菊|菊花茶寿眉|昆仑雪菊|玫瑰花|枸杞红枣|陈皮茶|菊|茉|罗汉果茶|猴王牌茶叶|六合茶|山楂叶|洛茶|桑叶|柠檬茶|宝梨"
      ) ~ "花茶",
      
      # 黑茶类
      stringr::str_detect(
        tea,
        "^黑茶$|藏茶"
      ) ~ "黑茶",
      
      # 苦丁茶类
      stringr::str_detect(
        tea,
        "苦丁茶|苦$"
      ) ~ "苦丁茶",
      
      # 罗布麻类
      stringr::str_detect(
        tea,
        "罗布麻"
      ) ~ "罗布麻茶",
      
      # 绞股蓝类
      stringr::str_detect(
        tea,
        "绞股蓝|绞胶蓝|搞股兰"
      ) ~ "绞股蓝茶",
      
      # 丹参类
      stringr::str_detect(
        tea,
        "丹参|丹心|丹参保心|天草丹参"
      ) ~ "丹参茶",
      
      # 灵芝类
      stringr::str_detect(
        tea,
        "灵芝"
      ) ~ "灵芝茶",
      
      # 谷物茶类
      stringr::str_detect(
        tea,
        "大麦茶|荞麦茶|苦荞麦|苦麦茶"
      ) ~ "谷物茶",
      
      # 保健茶类
      stringr::str_detect(
        tea,
        "保健茶|养生茶|降压茶|中药茶|舌斛茶|自制茶叶|祛湿|痛风茶|熊胆|溪黄草|东革阿里|丁香母藤茶"
      ) ~ "保健茶",
      
      # 牛蒡茶
      stringr::str_detect(
        tea,
        "牛蒡|牛膀"
      ) ~ "牛蒡茶",
      
      # 其他
      TRUE ~ str_replace_all(tea, "（.*", "")
    )
    
    merged_teas[i] <- merged_tea
  }
  
  # 权重映射（基于合并后的茶种类）
  tea_type_weights <- c(
    # 高 catechin
    "绿茶"     = 3,
    "白茶"     = 3,
    "黄茶"     = 3,
    "苦丁茶"   = 3,
    
    # 中等
    "乌龙茶"   = 2,
    "单丛"     = 2,    # 单丛属于乌龙茶类，同样中等
    "红茶"     = 2,
    "普洱"     = 2,
    "黑茶"     = 2,
    
    # 低 catechin

    "花茶"     = 1,
    "丹参茶"   = 1,    # 草本茶，归低
    "灵芝茶"   = 1,
    "罗布麻茶" = 1,
    "绞股蓝茶" = 1,
    "谷物茶"   = 1,
    "保健茶"   = 1,
    "牛蒡茶"   = 1,
    "决明子"   = 1,
    "减肥茶"   = 1,
    "山楂茶"   = 1
  )
  
  # 获取权重
  scores <- tea_type_weights[merged_teas]
  scores <- scores[!is.na(scores)]
  
  # 没识别到
  if (length(scores) == 0) {
    return(NA_real_)
  }
  
  # 多种茶取平均
  mean(scores)
}

# 测试单个茶种
get_tea_type_score("乌龙茶")  # 应该返回 2
get_tea_type_score("普洱")    # 应该返回 1
get_tea_type_score("绿茶")    # 应该返回 3
get_tea_type_score("单丛")    # 应该返回 2

# 测试多种茶混合
get_tea_type_score("乌龙茶；绿茶")  # 应该返回 (2+3)/2 = 2.5
get_tea_type_score("普洱，红茶")    # 应该返回 (1+1)/2 = 1

#### 统一处理函数
process_tea_score <- function(
    data,
    code_col,
    drink_col,
    freq_col,
    strength_col,
    amount_col,
    type_col,
    visit,
    duration_col = NULL,
    remove_prefix = NULL
){
  
  data %>%
    mutate(
      .drink_tmp = as.numeric(.data[[drink_col]]),
      .type_tmp  = str_trim(as.character(.data[[type_col]])),
      .type_tmp  = na_if(.type_tmp, ""),
      .type_tmp  = na_if(.type_tmp, "NA")
    ) %>%
    
    # 去掉 ID 缺失的空行
    filter(
      !is.na(.data[[code_col]]),
      as.character(.data[[code_col]]) != ""
    ) %>%
    
    # 只剔除“喝茶者中茶种类为纯数字”的行
    # 不喝茶者 tea_type 可以是 NA，不应该删
    filter(
      !(
        .drink_tmp == 1 &
          !is.na(.type_tmp) &
          str_detect(.type_tmp, "^[0-9]+$")
      )
    ) %>%
    
    # 只剔除“喝茶者但茶种类为空”的行
    # 不喝茶者茶种类为空是正常的
    filter(
      !(
        .drink_tmp == 1 &
          (is.na(.type_tmp) | .type_tmp == "")
      )
    ) %>%
    
    transmute(
      ID = if (is.null(remove_prefix)) {
        .data[[code_col]]
      } else {
        sub(remove_prefix, "", .data[[code_col]])
      },
      
      visit = visit,
      
      tea_drink = .drink_tmp,
      
      tea_freq = case_when(
        tea_drink == 0 ~ 0,
        tea_drink == 1 ~ as.numeric(.data[[freq_col]]),
        TRUE ~ NA_real_
      ),
      
      tea_strength = case_when(
        tea_drink == 0 ~ 0,
        tea_drink == 1 ~ 6 - as.numeric(.data[[strength_col]]),
        TRUE ~ NA_real_
      ),
      
      tea_amount = case_when(
        tea_drink == 0 ~ 0,
        tea_drink == 1 ~ as.numeric(.data[[amount_col]]),
        TRUE ~ NA_real_
      ),
      
      tea_type_score = case_when(
        tea_drink == 0 ~ 0,
        tea_drink == 1 ~ sapply(.type_tmp, get_tea_type_score),
        TRUE ~ NA_real_
      ),
      
      tea_duration = if (is.null(duration_col)) {
        NA_real_
      } else {
        as.numeric(.data[[duration_col]])
      }
    )
}


#### F0
table(Tea_F0_filter$常饮用茶种类_F0)
Tea_F0_long <- process_tea_score(
  
  data = Tea_F0_filter,
  
  code_col = "CODE_F0",
  
  drink_col = "是否喝茶_F0",
  
  freq_col = "每周冲茶次数_F0",
  
  strength_col = "喝茶浓淡_F0",
  
  amount_col = "过去一年茶叶总量_斤_F0",
  
  type_col = "常饮用茶种类_F0",
  
  duration_col = "累计喝茶总年限_F0",
  
  visit = "F0"
)

#### F1

Tea_F1_long <- process_tea_score(
  
  data = Tea_F1_filter,
  
  code_col = "CODE_F1",
  
  drink_col = "是否喝茶_F1",
  
  freq_col = "每周冲茶次数_F1",
  
  strength_col = "喝茶浓淡_F1",
  
  amount_col = "过去一年茶叶总量_斤_F1",
  
  type_col = "常饮用茶种类_F1",
  
  visit = "F1",
  
  remove_prefix = "^F1"
)

#### F2

Tea_F2_long <- process_tea_score(
  
  data = Tea_F2_filter,
  
  code_col = "CODE_F2",
  
  drink_col = "过去一年是否经常喝茶_F2",
  
  freq_col = "每周泡茶次数_F2",
  
  strength_col = "喝茶浓淡程度_F2",
  
  amount_col = "过去一年茶叶总量_斤_F2",
  
  type_col = "常饮用茶种类_F2",
  
  visit = "F2",
  
  remove_prefix = "^F2"
)

#### F3

Tea_F3_long <- process_tea_score(
  
  data = Tea_F3_filter,
  
  code_col = "CODE_F3",
  
  drink_col = "是否喝茶_F3",
  
  freq_col = "每周冲茶次数_F3",
  
  strength_col = "喝茶浓淡_F3",
  
  amount_col = "过去一年茶叶总量_斤_F3",
  
  type_col = "常饮用茶种类_F3",
  
  visit = "F3",
  
  remove_prefix = "^F3"
)

#### 合并长数据

Tea_long <- bind_rows(
  
  Tea_F0_long,
  Tea_F1_long,
  Tea_F2_long,
  Tea_F3_long
)

#### duration 只来自 baseline,意味着没有4000之后的了

duration_base <- Tea_F0_long %>%
  
  dplyr::select(
    ID,
    tea_duration
  )

Tea_long <- Tea_long %>%
  
  dplyr::select(
    -tea_duration
  ) %>%
  
  left_join(
    duration_base,
    by = "ID"
  )

Tea_long <- Tea_long[!is.na (Tea_long$tea_duration),]

Tea_long <- Tea_long[complete.cases(Tea_long),]
#### 标准化

Tea_long <- Tea_long %>%
  
  mutate(
    
    z_freq =
      as.numeric(scale(tea_freq)),
    
    z_amount =
      as.numeric(scale(tea_amount)),
    
    z_strength =
      as.numeric(scale(tea_strength)),
    
    z_duration =
      as.numeric(scale(tea_duration)),
    
    z_type =
      as.numeric(scale(tea_type_score))
  )

#### 原始综合 score

Tea_long <- Tea_long %>%
  
  mutate(
    
    tea_score_raw =
      
      z_freq +
      z_amount +
      z_strength +
      z_duration +
      z_type
  )

#### 仅喝茶者内部标准化

drinkers_index <- 
  
  !is.na(Tea_long$tea_drink) &
  Tea_long$tea_drink == 1

Tea_long$tea_score_z <- NA

Tea_long$tea_score_z[drinkers_index] <-
  
  as.numeric(
    
    scale(
      
      Tea_long$tea_score_raw[
        drinkers_index
      ]
    )
  )

#### 喝茶者中位数

median_score <- median(
  
  Tea_long$tea_score_z[
    drinkers_index
  ],
  
  na.rm = TRUE
)

#### 三分类 tea exposure

Tea_long <- Tea_long %>%
  
  mutate(
    
    tea_exposure_group = case_when(
      
      ## 不喝茶
      
      tea_drink == 0 ~ 0,
      
      ## 低 tea exposure
      
      tea_drink == 1 &
        tea_score_z < median_score ~ 1,
      
      ## 高 tea exposure
      
      tea_drink == 1 &
        tea_score_z >= median_score ~ 2,
      
      TRUE ~ NA_real_
    )
  )

#### factor
Tea_long <- Tea_long %>%
  
  mutate(
    
    tea_exposure_group =
      
      factor(
        
        tea_exposure_group,
        
        levels = c(0,1,2),
        
        labels = c(
          "Non-drinker",
          "Lower tea exposure",
          "Higher tea exposure"
        )
      )
  )
table(Tea_long$visit)
Tea_long_sup <- Tea_long[,c("ID","visit", "tea_exposure_group")]
Tea_long_sup <- Tea_long_sup[complete.cases(Tea_long_sup),]
table(Tea_long_sup$tea_exposure_group)
####data-饮茶纵向&骨骼肌/握力
data_mul_long <- merge(ASM_Height_long, Tea_long_sup, by.x = c("ID", "Followup") , by.y = c("ID", "visit"))

data_list <- list(Cov_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], data_mul_long)
data_mul_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
# 合并握力等，要用left_join
data_mul_v2 <- left_join(data_mul_v1, SPPB_grip_long, by = c("ID","Followup"))
table(SPPB_grip_long$Followup)

data_mul_final <- scale_columns_group(data_mul_v2, c("ASM","ASMI","grip_max"), "Sex_F0")

# 统计
table(data_mul_final$Followup)
a <- unique(data_mul_final$ID)


####lmer####
Lmer_results_mul <- process_lmer(c("tea_exposure_group"),
                                  c("ASM_z","ASMI_z", "grip_max_z"), #,"SPPB_total"
                                  data_mul_final, 
                                  Covariates_all_lmer, 
                                  2)
Lmer_results_mul$Level <- recode(Lmer_results_mul$Level,
                                 "Lower tea exposure" = "1",
                                 "Higher tea exposure" = "2"
                                 )
####茶的多方面与总分的相关性####
# 1. 仅保留饮茶者

Tea_long_sta <- Tea_long %>%
  filter(tea_drink == 1)

# 2. 设置变量

cor_vars <- c(
  "z_freq",
  "z_amount",
  "z_strength",
  "z_duration",
  "z_type",
  "tea_score_raw"
)

cor_labels <- c(
  z_freq        = "Tea frequency\n(z-score)",
  z_amount      = "Tea amount\n(z-score)",
  z_strength    = "Tea strength\n(z-score)",
  z_duration    = "Tea duration\n(z-score)",
  z_type        = "Tea type\n(z-score)",
  tea_score_raw = "Composite tea\nexposure score"
)

cor_method <- "spearman"

# 3. 检查变量是否存在

needed_vars <- c(
  "visit",
  cor_vars
)

missing_vars <- setdiff(
  needed_vars,
  colnames(Tea_long_sta)
)

if (length(missing_vars) > 0) {
  stop(
    paste0(
      "以下变量在 Tea_long_sta 中不存在：",
      paste(missing_vars, collapse = ", ")
    )
  )
}

# 4. 整理 visit 顺序

visit_levels <- Tea_long_sta %>%
  filter(!is.na(visit)) %>%
  distinct(visit) %>%
  mutate(
    visit_chr = as.character(visit),
    visit_num = suppressWarnings(
      as.numeric(stringr::str_extract(visit_chr, "\\d+"))
    )
  ) %>%
  arrange(visit_num, visit_chr) %>%
  pull(visit)

# 5. P值格式函数和显著性星号函数

format_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(
      p < 0.001,
      "<0.001",
      sprintf("%.3f", p)
    )
  )
}

p_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

# 6. 安全相关性函数

cor_test_safe <- function(x, y, method = "spearman") {
  
  df <- tibble::tibble(
    x = suppressWarnings(as.numeric(x)),
    y = suppressWarnings(as.numeric(y))
  ) %>%
    filter(
      !is.na(x),
      !is.na(y)
    )
  
  n <- nrow(df)
  
  if (n < 3) {
    return(
      tibble::tibble(
        n = n,
        r = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  if (
    sd(df$x, na.rm = TRUE) == 0 ||
    sd(df$y, na.rm = TRUE) == 0
  ) {
    return(
      tibble::tibble(
        n = n,
        r = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  test <- suppressWarnings(
    cor.test(
      df$x,
      df$y,
      method = method,
      exact = FALSE
    )
  )
  
  tibble::tibble(
    n = n,
    r = unname(test$estimate),
    p_value = test$p.value
  )
}

# 7. 整理分析数据

cor_data_all <- Tea_long_sta %>%
  dplyr::select(
    visit,
    all_of(cor_vars)
  ) %>%
  mutate(
    visit = factor(
      visit,
      levels = visit_levels
    ),
    across(
      all_of(cor_vars),
      ~ suppressWarnings(as.numeric(.))
    )
  )

# 8. 分 visit 计算 6 × 6 两两相关性

cor_result_list <- list()
counter <- 1

for (v in visit_levels) {
  
  dat_v <- cor_data_all %>%
    filter(visit == v)
  
  pair_grid <- expand.grid(
    var_y = cor_vars,
    var_x = cor_vars,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(pair_grid))) {
    
    vx <- pair_grid$var_x[i]
    vy <- pair_grid$var_y[i]
    
    vx_id <- match(vx, cor_vars)
    vy_id <- match(vy, cor_vars)
    
    if (vx == vy) {
      
      n_i <- sum(!is.na(dat_v[[vx]]))
      
      cor_i <- tibble::tibble(
        n = n_i,
        r = 1,
        p_value = NA_real_
      )
      
    } else {
      
      cor_i <- cor_test_safe(
        x = dat_v[[vx]],
        y = dat_v[[vy]],
        method = cor_method
      )
    }
    
    cor_result_list[[counter]] <- tibble::tibble(
      visit = v,
      var_x = vx,
      var_y = vy,
      var_x_id = vx_id,
      var_y_id = vy_id,
      var_x_label = cor_labels[vx],
      var_y_label = cor_labels[vy],
      pair_key = paste0(
        pmin(vx_id, vy_id),
        "__",
        pmax(vx_id, vy_id)
      ),
      n = cor_i$n,
      r = cor_i$r,
      p_value = cor_i$p_value
    )
    
    counter <- counter + 1
  }
}

cor_results_by_visit <- bind_rows(cor_result_list)

# 9. 每个 visit 内对非对角线相关性做 FDR 校正

p_adjust_df <- cor_results_by_visit %>%
  filter(var_x != var_y) %>%
  distinct(
    visit,
    pair_key,
    .keep_all = TRUE
  ) %>%
  group_by(visit) %>%
  mutate(
    p_FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  dplyr::select(
    visit,
    pair_key,
    p_FDR
  )

cor_results_by_visit <- cor_results_by_visit %>%
  left_join(
    p_adjust_df,
    by = c("visit", "pair_key")
  ) %>%
  mutate(
    p_FDR = ifelse(
      var_x == var_y,
      NA_real_,
      p_FDR
    ),
    sig = p_star(p_FDR),
    r_round = round(r, 2),
    p_label = format_p(p_value),
    p_FDR_label = format_p(p_FDR),
    
    heatmap_label = ifelse(
      var_x == var_y,
      "1.00",
      ifelse(
        is.na(r),
        "NA",
        paste0(
          sprintf("%.2f", r),
          sig
        )
      )
    )
  )

# 10. 只保留下三角矩阵，避免重复

plot_cor_matrix <- cor_results_by_visit %>%
  filter(
    var_y_id >= var_x_id
  ) %>%
  mutate(
    visit = factor(
      visit,
      levels = visit_levels
    ),
    var_x_label = factor(
      var_x_label,
      levels = cor_labels[cor_vars]
    ),
    var_y_label = factor(
      var_y_label,
      levels = rev(cor_labels[cor_vars])
    )
  )

# 查看结果
print(cor_results_by_visit)

write.csv(
  cor_results_by_visit,
  "tea_score_six_variables_pairwise_correlation_by_visit.csv",
  row.names = FALSE
)

# 11. 绘制分 visit 的相关性热图

p_cor_heatmap <- ggplot(
  plot_cor_matrix,
  aes(
    x = var_x_label,
    y = var_y_label,
    fill = r
  )
) +
  geom_tile(
    color = "white",
    linewidth = 1.05
  ) +
  geom_text(
    aes(label = heatmap_label),
    size = 4.2,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman r"
  ) +
  facet_wrap(
    ~ visit,
    nrow = 1
  ) +
  coord_fixed() +
  labs(
    x = NULL,
    y = NULL,
    title = NULL,
    subtitle = NULL,
    caption = "Values are Spearman correlation coefficients. * FDR-adjusted P < 0.05; ** FDR-adjusted P < 0.01; *** FDR-adjusted P < 0.001."
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      color = "black",
      size = 14
    ),
    
    axis.text.x = element_text(
      face = "bold",
      color = "black",
      size = 10,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    axis.text.y = element_text(
      face = "bold",
      color = "black",
      size = 10
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    
    legend.text = element_text(
      color = "black",
      size = 11
    ),
    
    plot.caption = element_text(
      hjust = 0,
      size = 10,
      color = "black"
    ),
    
    plot.margin = ggplot2::margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    )
  )

p_cor_heatmap

# 12. 保存图片

ggsave(
  filename = "tea_score_six_variables_pairwise_correlation_heatmap_by_visit.png",
  plot = p_cor_heatmap,
  width = 18,
  height = 6.5,
  dpi = 500
)

####可视化####
# 1. 空图函数：防止某些图为 NULL 时拼图报错

make_blank_plot <- function(title = "") {
  ggplot() +
    theme_void() +
    labs(title = title) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 13
      )
    )
}

# 2. 统一主题函数

clean_one_plot <- function(p) {
  
  if (is.null(p)) {
    return(make_blank_plot())
  }
  
  p +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 13
      ),
      axis.title = element_text(
        face = "bold",
        size = 15
      ),
      axis.text = element_text(
        color = "black",
        size = 12
      ),
      legend.title = element_text(
        face = "bold"
      )
    )
}

# 3. 加 A、B、C... 图标函数

add_panel_tag <- function(p, tag) {
  
  if (is.null(p)) {
    p <- make_blank_plot()
  }
  
  p +
    labs(tag = tag) +
    theme(
      plot.tag = element_text(
        face = "bold",
        size = 22,
        color = "black"
      ),
      plot.tag.position = c(0.01, 0.98),
      plot.margin = ggplot2::margin(
        t = 8,
        r = 8,
        b = 8,
        l = 8
      )
    )
}

# 4. 第一行：F0-F3 Weekly frequency，A-D

p_A <- add_panel_tag(clean_one_plot(res_F0$plots$weekly_frequency), "A")
p_B <- add_panel_tag(clean_one_plot(res_F1$plots$weekly_frequency), "B")
p_C <- add_panel_tag(clean_one_plot(res_F2$plots$weekly_frequency), "C")
p_D <- add_panel_tag(clean_one_plot(res_F3$plots$weekly_frequency), "D")

row_weekly_frequency <- p_A + p_B + p_C + p_D +
  plot_layout(nrow = 1)

# 5. 第二行：F0-F3 Tea type，E-H

p_E <- add_panel_tag(clean_one_plot(res_F0$plots$tea_type), "E")
p_F <- add_panel_tag(clean_one_plot(res_F1$plots$tea_type), "F")
p_G <- add_panel_tag(clean_one_plot(res_F2$plots$tea_type), "G")
p_H <- add_panel_tag(clean_one_plot(res_F3$plots$tea_type), "H")

row_tea_type <- p_E + p_F + p_G + p_H +
  plot_layout(nrow = 1)

# 6. 第三行：F0-F3 Tea amount，I-L

p_I <- add_panel_tag(clean_one_plot(res_F0$plots$tea_amount), "I")
p_J <- add_panel_tag(clean_one_plot(res_F1$plots$tea_amount), "J")
p_K <- add_panel_tag(clean_one_plot(res_F2$plots$tea_amount), "K")
p_L <- add_panel_tag(clean_one_plot(res_F3$plots$tea_amount), "L")

row_tea_amount <- p_I + p_J + p_K + p_L +
  plot_layout(nrow = 1)

# 7. 第四行：F0-F3 Brewing strength，M-P

p_M <- add_panel_tag(clean_one_plot(res_F0$plots$brewing_strength), "M")
p_N <- add_panel_tag(clean_one_plot(res_F1$plots$brewing_strength), "N")
p_O <- add_panel_tag(clean_one_plot(res_F2$plots$brewing_strength), "O")
p_P <- add_panel_tag(clean_one_plot(res_F3$plots$brewing_strength), "P")

row_brewing_strength <- p_M + p_N + p_O + p_P +
  plot_layout(nrow = 1)

# 8. 第五行：相关性热图单独一行，Q

p_cor_heatmap_clean <- p_cor_heatmap +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 13
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    ),
    axis.title = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      face = "bold",
      color = "black",
      size = 11,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      face = "bold",
      color = "black",
      size = 11
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.text = element_text(
      color = "black"
    ),
    plot.margin = ggplot2::margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )

p_Q <- add_panel_tag(
  p_cor_heatmap_clean,
  "Q"
)

row_cor_heatmap <- p_Q +
  plot_layout(nrow = 1)

# 9. 第六行：F0 Tea duration 单独一行，R
#    右侧加空白，让 duration 图不要被拉得太宽

p_R <- add_panel_tag(
  clean_one_plot(res_F0$plots$tea_duration),
  "R"
)

row_tea_duration <- p_R + plot_spacer() + plot_spacer() + plot_spacer() +
  plot_layout(
    nrow = 1,
    widths = c(1, 1, 1, 1)
  )

# 10. 总拼图

final_tea_plot <- row_weekly_frequency /
  row_tea_type /
  row_tea_amount /
  row_brewing_strength /
  row_cor_heatmap /
  row_tea_duration +
  plot_layout(
    heights = c(
      1,
      1.25,
      1,
      1,
      1.65,
      1
    )
  )

final_tea_plot

# 11. 保存图片

ggsave(
  filename = "final_tea_description_and_correlation.png",
  plot = final_tea_plot,
  width = 20,
  height = 26,
  dpi = 500
)
####排除早期结局有显著变化的观测####
early_change <- ASM_Height_long %>%
  filter(Followup %in% c("F1", "F2")) %>%
  dplyr::select(ID, Followup, ASMI) %>%
  pivot_wider(names_from = Followup, values_from = ASMI, names_prefix = "ASMI_") %>%
  mutate(
    change_F1_to_F2 = ASMI_F2 - ASMI_F1,
    pct_change = (ASMI_F2 - ASMI_F1) / ASMI_F1 * 100,
    significant_decline = pct_change < -10   # 下降超过10%定义为显著下降
  )

# 查看有多少人显著下降
table(early_change$significant_decline)

# 获取需要排除的ID
exclude_ids <- early_change %>%
  filter(significant_decline == TRUE) %>%
  pull(ID)

# 创建排除后的数据集
ASM_Height_excluded1 <- ASM_Height_long %>%
  filter(!ID %in% exclude_ids)

a <- unique(ASM_Height_excluded1$ID)
####lmer####
data_sen1_long <- merge(ASM_Height_excluded1, Tea_long_prime, by.x = c("ID", "Followup") , by.y = c("ID", "visit"))

data_list <- list(Cov_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], data_sen1_long)
data_sen1_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
# 合并握力等，要用left_join
data_sen1_v1 <- left_join(data_sen1_v1, SPPB_grip_long, by = c("ID","Followup"))


data_sen1_final <- scale_columns_group(data_sen1_v1, c("ASM","ASMI","grip_max"), "Sex_F0")

Lmer_results_sen1_all <- process_lmer(c("tea_freq_group"),
                                      c("ASM_z","ASMI_z","grip_max_z"), #,"SPPB_total"
                                      data_sen1_final, 
                                      Covariates_all_lmer, 
                                      2)
####排除F1就有严重ASMI的观测####
ASM_Height_excluded <- merge(ASM_Height_long, Cov_F0_final[,c("ID","Sex_F0")],by = "ID")

low_baseline_IDs <- ASM_Height_excluded %>%
  filter(Followup == "F1") %>%
  mutate(
    is_low_baseline = case_when(
      Sex_F0 == 0 & ASMI < 5.4 ~ TRUE,
      Sex_F0 == 1 & ASMI < 7   ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(is_low_baseline == TRUE) %>%
  pull(ID)

# 排除这些ID（但没有F1数据的ID会被保留）
ASM_Height_excluded2 <- ASM_Height_excluded %>%
  filter(!ID %in% low_baseline_IDs)

a <- unique(ASM_Height_excluded2$ID)

####lmer####
data_sen2_long <- merge(ASM_Height_excluded2, Tea_long_prime, by.x = c("ID", "Followup") , by.y = c("ID", "visit"))

data_list <- list(Cov_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], dplyr::select(data_sen2_long,-Sex_F0))
data_sen2_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
# 合并握力等，要用left_join
data_sen2_v1 <- left_join(data_sen2_v1, SPPB_grip_long, by = c("ID","Followup"))


data_sen2_final <- scale_columns_group(data_sen2_v1, c("ASM","ASMI","grip_max"), "Sex_F0")

Lmer_results_sen2_all <- process_lmer(c("tea_freq_group"),
                                      c("ASM_z","ASMI_z","grip_max_z"), #,"SPPB_total"
                                      data_sen2_final, 
                                      Covariates_all_lmer, 
                                      2)
####*****可视化:饮茶&骨骼肌/握力合并*****####
plot_data1 <- Lmer_results1_all[Lmer_results1_all$Predictor == "tea_freq_group" & Lmer_results1_all$Outcome %in% c("ASM_z","ASMI_z","grip_max_z"),]
plot_data1$Method <- "Method 1"

plot_data2 <- Lmer_results2_all
plot_data2$Method <- "Method 2"

plot_data3 <- Lmer_results3_all[Lmer_results3_all$Predictor == "tea_freq_group" & Lmer_results3_all$Outcome %in% c("ASM_z","ASMI_z","grip_max_z"),]
plot_data3$Method <- "Method 3"

plot_data4 <- Lmer_results4_all
plot_data4$Method <- "Method 4"

plot_data5 <- Lmer_results5_all
plot_data5$Method <- "Method 5"

plot_data8 <- Lmer_results_mul
plot_data8$Method <- "Method 6"

plot_data6 <- Lmer_results_sen2_all
plot_data6$Method <- "Sensitivity Analysis 1"

plot_data7 <- Lmer_results_sen1_all
plot_data7$Method <- "Sensitivity Analysis 2"

plot_data <- rbind(plot_data1, plot_data2, plot_data3, plot_data4, plot_data5, plot_data8, plot_data6, plot_data7)

# 假设数据处理部分保持不变
plot_data2 <- plot_data %>%
  mutate(
    Level = factor(Level, levels = c("1", "2"), labels = c("Group 1", "Group 2")),
    Outcome = factor(Outcome, levels = unique(Outcome)),
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    ),
    Outcome = factor(
      Outcome,
      levels = c(
        "Handgrip strength",
        "ASMI",
        "ASM"
      )),
    sig = ifelse(P_value < 0.05, "*", "")
  )

p_refined <- ggplot(
  plot_data2,
  aes(x = Estimate, y = fct_rev(Outcome), color = Level, group = Level)
) +
  # 增加背景参考线（只留横向，辅助阅读）
  geom_hline(aes(yintercept = Outcome), color = "grey95", linewidth = 4) +
  
  # 零线：用更细的深灰色实线或虚线
  geom_vline(xintercept = 0, linetype = "solid", color = "grey70", size = 0.6) +
  
  # 95% CI：末端加一个小垂直线（可选），加宽线径
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0, # 设置为0更现代，或者保持0.15
    linewidth = 1.2, 
    alpha = 0.7,
    position = position_dodge(width = 0.6)
  ) +
  
  # 点：增加白色描边，更有质感
  geom_point(
    size = 4,
    shape = 21,
    fill = "white", # 内部填充白色，边框显色
    stroke = 1.5,
    position = position_dodge(width = 0.6)
  ) +
  
  # 显著性：位置微调，避免重叠
  geom_text(
    aes(label = sig, x = CI_high), # 放在置信区间上限的右侧
    size = 6,
    hjust = -0.3,
    vjust = 0.8,
    fontface = "bold",
    position = position_dodge(width = 0.6),
    show.legend = FALSE
  ) +
  
  # 分面：修饰 strip 背景
  facet_wrap(~ Method, nrow = 4, scales = "free_y") +
  
  # 配色：使用更清爽的配色
  scale_color_manual(
    values = c("Group 1" = "#21918c", "Group 2" = "pink") # 或者使用 c("#5050FF", "#FF5050")
  ) +
  
  labs(
    x = "β (95% CI)",
    y = NULL,
    color = "Tea Consumption"
  ) +
  
  # 深度定制主题
  theme_minimal(base_size = 14) +
  theme(
    # 图例置于上方并左对齐
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(size = 12, face = "bold"),
    legend.background = element_blank(),
    
    # 坐标轴美化
    axis.line.x = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.text.y = element_text(size = 13, face = "bold", color = "grey20"),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    
    # 分面标签美化
    strip.background = element_rect(fill = "grey20", color = NA),
    strip.text = element_text(color = "white", face = "bold", size = 13),
    
    # 网格线微调
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92"),
    
    # 间距
    panel.spacing = unit(1.5, "lines"),
    plot.margin = ggplot2::margin(20, 30, 20, 20)
  )

# 建议导出为更高质量的 PDF 或带透明度的 PNG
#ggsave("Tea_ASMI_Improved.pdf", p_refined, height =8, width = 12, device = cairo_pdf)
ggsave("Tea_ASMI_Improved.png", p_refined, width = 12, height = 12, dpi = 700)
####最终纳入分析与没纳入分析的受试者比较####
included_ids <- unique(
  data1_final$ID
)
#### 基线队列
all_ids <- unique(
  Cov_F0_final$ID
)
#### 排除受试者
excluded_ids <- setdiff(
  all_ids,
  included_ids
)
length(all_ids)
length(included_ids)
length(excluded_ids)


compare_data <- Cov_F0_final %>%
  
  mutate(
    
    Included =
      ifelse(
        ID %in% included_ids,
        "Included",
        "Excluded"
      )
  )
table(compare_data$Included)

Table_compare <- descrTable(Included ~ Age_F0+Sex_F0+Education_F0+Income_F0+
                        Met_mean+Protein_mean+Fiber_soluble_mean+Total_score+coffee+
                        Smoke_F0+Alcohol_F0+Calcium_F0+Vitamin_F0+Oil_sup_follow+Fracture_F0+Disease_F0+Estrogen_F0+Menopause_F0,
                      data = compare_data,method =NA,show.all = TRUE)#show.all = TRUE :显示all

export2word(Table_compare, file='Table_compare.docx')
####*****************************饮茶&骨骼肌(结局是二分类）*****************************####
####data1-饮茶纵向&骨骼肌/握力####
data_list <- list(Cov_F0_final, Tea_F0_final, ASMI_wide_low[,c("ID","low_ASMI","WBTOT_FAT_mean_F1234")])
data1_OR_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))

data1_OR_v2 <- left_join(data1_OR_v1, SPPB_grip_wide, , by = "ID")
data1_OR_v3 <- left_join(data1_OR_v2, Catechin, by = "ID")

data1_OR_final <- group_variables3(
  data = data1_OR_v3,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)
a <- unique(data1_OR_final$ID)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(data1_OR_final)%in% factor_name)
for(i in idx ){
  data1_OR_final[[i]] <-  as.factor(data1_OR_final[[i]])
}


Binary_results <- process_logistic(
  X = c("tea_freq_group","serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T"),
  Y = c("low_ASMI", "low_grip", "low_gait_speed"), #,"low_ASMI", "low_grip", "low_gait_speed"
  data = data1_OR_final,
  covariates = Covariates_all,
  gender_label = 2
)
#**********正文写作
a <- Binary_results[Binary_results$P_value < 0.05,]
a$OR <- round(a$OR, 3)
a$OR_95CI_low <- round(a$OR_95CI_low, 3)
a$OR_95CI_high <- round(a$OR_95CI_high, 3)

min(a$OR)
max(a$OR)
min(a$OR_95CI_low)
max(a$OR_95CI_high)
####验证logistic回归####
x <- "tea_freq_group"
y <- "low_ASMI"

formula_text <- paste(
  y, "~",
  x, "+",
  paste(Covariates_all, collapse = " + ")
)
my_formula <- as.formula(formula_text)

model <- glm(my_formula, data = data1_OR_v1, family = binomial)
summary(model)
####可视化####
Binary_plot <- Binary_results %>%
  mutate(
    Predictor_level_num = as.numeric(as.character(Predictor_level))
  ) %>%
  group_by(Gender, Outcome, Predictor) %>%
  filter(Predictor_level_num == max(Predictor_level_num, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Predictor_label = recode(
      Predictor,
      !!!predictor_rename,
      .default = Predictor
    ),
    
    Outcome_label = case_when(
      Outcome == "low_ASMI" ~ "Low ASMI",
      Outcome == "low_grip" ~ "Low handgrip strength",
      Outcome == "low_gait_speed" ~ "Low walking speed",
      TRUE ~ Outcome
    ),
    
    Level_label = get_level_label(Predictor, Predictor_level),
    
    Exposure_label = paste0(Predictor_label, "\n", Level_label),
    
    # beta 的 95% CI，仍然用于画图横轴
    Beta_95CI_low = log(OR_95CI_low),
    Beta_95CI_high = log(OR_95CI_high),
    
    # 右侧文本显示 OR 和 95%CI，保留三位小数
    OR_CI_label = sprintf(
      "%.3f (%.3f, %.3f)",
      OR, OR_95CI_low, OR_95CI_high
    ),
    
    Sig_group = case_when(
      P_value < 0.05 ~ "P < 0.05",
      TRUE ~ "P ≥ 0.05"
    )
  )

predictor_order <- unname(predictor_rename)

exposure_order <- Binary_plot %>%
  mutate(
    Predictor_label = factor(Predictor_label, levels = predictor_order)
  ) %>%
  arrange(Predictor_label) %>%
  distinct(Exposure_label) %>%
  pull(Exposure_label)

Binary_plot <- Binary_plot %>%
  mutate(
    Outcome_label = factor(
      Outcome_label,
      levels = c(
        "Low ASMI",
        "Low handgrip strength",
        "Low walking speed"
      )
    ),
    Exposure_label = factor(
      Exposure_label,
      levels = rev(exposure_order)
    )
  )


plot_beta_forest <- function(plot_data,
                             title = NULL,
                             xlab = "Odds ratio (95% CI)",
                             output_file = NULL,
                             width = 10,
                             height = 7) {
  
  x_min <- min(plot_data$OR_95CI_low, na.rm = TRUE)
  x_max <- max(plot_data$OR_95CI_high, na.rm = TRUE)
  x_range <- x_max - x_min
  
  text_x <- x_max + 0.05 * x_range
  x_limit_right <- x_max + 0.75 * x_range
  
  p <- ggplot(
    plot_data,
    aes(
      x = OR,
      y = Exposure_label,
      xmin = OR_95CI_low,
      xmax = OR_95CI_high
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.6,
      color = "grey45"
    ) +
    geom_errorbarh(
      aes(color = Sig_group),
      height = 0.18,
      linewidth = 0.8,
      alpha = 0.95
    ) +
    geom_point(
      aes(fill = Sig_group),
      shape = 21,
      size = 3.6,
      color = "white",
      stroke = 0.45
    ) +
    geom_text(
      aes(label = OR_CI_label),
      x = text_x,
      hjust = 0,
      size = 3.4,
      color = "grey15"
    ) +
    facet_grid(
      .~  Outcome_label,
      scales = "free_y",
      space = "free_y"
    ) +
    scale_color_manual(
      values = c(
        "P < 0.05" = "#B22222",
        "P ≥ 0.05" = "#2C3E50"
      )
    ) +
    scale_fill_manual(
      values = c(
        "P < 0.05" = "#B22222",
        "P ≥ 0.05" = "#2C3E50"
      )
    ) +
    coord_cartesian(
      xlim = c(x_min - 0.08 * x_range, x_limit_right),
      clip = "off"
    ) +
    labs(
      title = title,
      x = xlab,
      y = NULL,
      color = NULL,
      fill = NULL,
      caption = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 10)
      ),
      plot.caption = element_text(
        size = 10,
        color = "grey35",
        hjust = 0,
        margin = ggplot2::margin(t = 8)
      ),
      axis.title.x = element_text(
        size = 12.5,
        face = "bold",
        margin = ggplot2::margin(t = 8)
      ),
      axis.text.x = element_text(
        size = 10.5,
        color = "grey20"
      ),
      axis.text.y = element_text(
        size = 10.5,
        color = "grey10",
        lineheight = 0.95
      ),
      strip.background = element_rect(
        fill = "#34495E",
        color = NA
      ),
      strip.text = element_text(
        size = 13,
        face = "bold",
        color = "white",
        margin = ggplot2::margin(t = 6, b = 6)
      ),
      panel.grid.major.y = element_line(
        color = "grey90",
        linewidth = 0.35
      ),
      panel.grid.major.x = element_line(
        color = "grey88",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.justification = "center",
      plot.margin = ggplot2::margin(10, 130, 10, 10)
    )
  
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = p,
      width = width,
      height = height,
      dpi = 600
    )
  }
  
  return(p)
}

p_beta <- plot_beta_forest(
  Binary_plot,
  title = NULL,
  output_file = "beta_forest_plot.pdf",
  width = 14,
  height = 6
)

ggsave(
  "beta_forest_plot.png",
  p_beta,
  width = 14,
  height = 6,
  dpi = 600
)

####*****************************其他分析*****************************####
####Table 1####
data_table1 <- data1_final %>%
  arrange(ID, Followup) %>%
  group_by(ID) %>%
  dplyr::slice(1) %>%
  ungroup()

Table1a <- descrTable(tea_freq_group ~ Age_F0+Sex_F0+Education_F0+Income_F0+
                        WBTOT_FAT_mean_F1234+
                        Met_mean+Protein_mean+Fiber_soluble_mean+Total_score+coffee+
                        Smoke_F0+Alcohol_F0+Calcium_F0+Vitamin_F0+Oil_sup_follow+Fracture_F0+Disease_F0+Estrogen_F0+Menopause_F0,
                      data = data_table1,method =NA,show.all = TRUE)#show.all = TRUE :显示all

export2word(Table1a, file='Table1.docx')
####血清儿茶素类统计描述####
data_list <- list(Tea_F0_final, Catechin)
data1_flo <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
a <- Catechin[!Catechin$ID %in% data1_flo$ID,]

data1_flo <- group_variables3(
  data = data1_flo,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(data1_flo)%in% factor_name)
for(i in idx ){
  data1_flo[[i]] <-  as.factor(data1_flo[[i]])
}

data1_flo <- data1_flo %>%
  mutate(
    tea_freq_group = factor(
      recode(
        tea_freq_group,
        "0" = "Non-drinker",
        "1" = "0-6 times/week",
        "2" = ">= 7 times/week"
      ),
      levels = c(
        "Non-drinker",
        "0-6 times/week",
        ">= 7 times/week"
      )
    )
  ) %>%
  mutate(
    across(
      c(
        serum_I_epicatechin_F0_T,
        serum_I_EGC_F0_T
      ),
      ~ factor(
        recode(
          .x,
          "0" = "Undetectable",
          "1" = "Low",
          "2" = "High"
        ),
        levels = c(
          "Undetectable",
          "Low",
          "High"
        )
      )
    )
  ) %>%
  mutate(
    across(
      c(
        serum_I_catechin_F0_T,
        serum_I_EGCG_F0_T,
        serum_I_ECG_F0_T
      ),
      ~ factor(
        recode(
          .x,
          "0" = "Tertile 1",
          "1" = "Tertile 2",
          "2" = "Tertile 3"
        ),
        levels = c(
          "Tertile 1",
          "Tertile 2",
          "Tertile 3"
        )
      )
    )
  ) 

colnames(data1_flo) <- recode(colnames(data1_flo),
                              "serum_I_catechin_F0_T" = "Catechin",
                              "serum_I_epicatechin_F0_T"  = "Epicatechin",
                              "serum_I_EGC_F0_T" =  "Epigallocatechin",
                              "serum_I_EGCG_F0_T" = "Epigallocatechin gallate",
                              "serum_I_ECG_F0_T"  = "Epicatechin gallate"
                              )

#**********************可视化
# 1. 提取并统计数据
target_vars <- c("Catechin", "Epicatechin", "Epigallocatechin", 
                 "Epigallocatechin gallate", "Epicatechin gallate")

summary_table <- data1_flo %>%
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
  "Epigallocatechin gallate",
  "Epicatechin gallate",
  "Epicatechin", 
  "Epigallocatechin"
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
    size = 4.5,
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
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", color = "grey20", size = 15),
    axis.text.y = element_text(color = "grey30"),
    axis.title.y = element_text(face = "bold", size = 18, margin = ggplot2::margin(r=10)),
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
# ggsave(
#   "Custom_Catechin_Distribution.png",
#   p,
#   width = 14,
#   height = 9,
#   dpi = 700
# )

# ggsave("Custom_Catechin_Distribution.pdf", p, width = 14, height = 9)
####血清儿茶素类统计描述####
data1_flo_all <- group_variables3(
  data = Catechin,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(data1_flo_all)%in% factor_name)
for(i in idx ){
  data1_flo_all[[i]] <-  as.factor(data1_flo_all[[i]])
}

data1_flo_all <- data1_flo_all %>%
  mutate(
    across(
      c(
        serum_I_epicatechin_F0_T,
        serum_I_EGC_F0_T
      ),
      ~ factor(
        recode(
          .x,
          "0" = "Undetectable",
          "1" = "Low",
          "2" = "High"
        ),
        levels = c(
          "Undetectable",
          "Low",
          "High"
        )
      )
    )
  ) %>%
  mutate(
    across(
      c(
        serum_I_catechin_F0_T,
        serum_I_EGCG_F0_T,
        serum_I_ECG_F0_T
      ),
      ~ factor(
        recode(
          .x,
          "0" = "Tertile 1",
          "1" = "Tertile 2",
          "2" = "Tertile 3"
        ),
        levels = c(
          "Tertile 1",
          "Tertile 2",
          "Tertile 3"
        )
      )
    )
  ) 

colnames(data1_flo_all) <- recode(colnames(data1_flo_all),
                              "serum_I_catechin_F0_T" = "Catechin",
                              "serum_I_epicatechin_F0_T"  = "Epicatechin",
                              "serum_I_EGC_F0_T" =  "Epigallocatechin",
                              "serum_I_EGCG_F0_T" = "Epigallocatechin gallate",
                              "serum_I_ECG_F0_T"  = "Epicatechin gallate"
)

#**********************可视化
# 1. 提取并统计数据
target_vars <- c("Catechin", "Epicatechin", "Epigallocatechin", 
                 "Epigallocatechin gallate", "Epicatechin gallate")

summary_table <- data1_flo_all %>%
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
  "Epigallocatechin gallate",
  "Epicatechin gallate",
  "Epicatechin", 
  "Epigallocatechin"
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
    size = 4.5,
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
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", color = "grey20", size = 15),
    axis.text.y = element_text(color = "grey30"),
    axis.title.y = element_text(face = "bold", size = 18, margin = ggplot2::margin(r=10)),
    # 网格线淡化
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linetype = "dotted"),
    # 图例
    legend.position = "none",
    panel.spacing = unit(1.5, "lines"),
    plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )

ggsave(
  "Custom_Catechin_Distribution.png",
  p,
  width = 14,
  height = 9,
  dpi = 700
)
####基线饮茶习惯&儿茶素浓度####
Table1A_all <- descrTable(tea_freq_group ~ Catechin  + `Epigallocatechin gallate` + `Epicatechin gallate`+ Epicatechin + Epigallocatechin,
                          data = data1_flo, method =NA,show.all = FALSE) # show.all = TRUE :显示all
export2word(Table1A_all, file='Table2.docx')
####基线饮茶分数&儿茶素浓度####
Tea_score_catechin <- merge(Tea_long_sup[Tea_long_sup$visit == "F0",], Catechin, by = "ID")

Tea_score_catechin <- group_variables3(
  data = Tea_score_catechin,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(Tea_score_catechin)%in% factor_name)
for(i in idx ){
  Tea_score_catechin[[i]] <-  as.factor(Tea_score_catechin[[i]])
}

Tea_score_catechin <- Tea_score_catechin %>%
  mutate(
    across(
      c(
        serum_I_epicatechin_F0_T,
        serum_I_EGC_F0_T
      ),
      ~ factor(
        recode(
          .x,
          "0" = "Undetectable",
          "1" = "Low",
          "2" = "High"
        ),
        levels = c(
          "Undetectable",
          "Low",
          "High"
        )
      )
    )
  ) %>%
  mutate(
    across(
      c(
        serum_I_catechin_F0_T,
        serum_I_EGCG_F0_T,
        serum_I_ECG_F0_T
      ),
      ~ factor(
        recode(
          .x,
          "0" = "Tertile 1",
          "1" = "Tertile 2",
          "2" = "Tertile 3"
        ),
        levels = c(
          "Tertile 1",
          "Tertile 2",
          "Tertile 3"
        )
      )
    )
  ) 



colnames(Tea_score_catechin) <- recode(colnames(Tea_score_catechin),
                                  "serum_I_catechin_F0_T" = "Catechin",
                                  "serum_I_epicatechin_F0_T"  = "Epicatechin",
                                  "serum_I_EGC_F0_T" =  "Epigallocatechin",
                                  "serum_I_EGCG_F0_T" = "Epigallocatechin gallate",
                                  "serum_I_ECG_F0_T"  = "Epicatechin gallate"
)
Table1A_all <- descrTable(tea_exposure_group ~ Catechin  + `Epigallocatechin gallate` + `Epicatechin gallate`+ Epicatechin + Epigallocatechin,
                          data = Tea_score_catechin, method =NA,show.all = FALSE) # show.all = TRUE :显示all
export2word(Table1A_all, file='Table3.docx')
####*******茶的种类的详细分析*******####
####茶的种类####
Tea_type_only_F0 <- Tea_F0_filter[,c("CODE_F0", "是否喝茶_F0", "常饮用茶种类_F0" )]
colnames(Tea_type_only_F0) <- recode(colnames(Tea_type_only_F0) ,
                                  "CODE_F0" = "ID",
                                  "是否喝茶_F0" = "Tea_drink",
                                  "常饮用茶种类_F0" = "Tea_type")
Tea_type_only_F0$Times <- "F0"

Tea_type_only_F1 <- Tea_F1_filter[,c("CODE_F1", "是否喝茶_F1", "常饮用茶种类_F1"  )]
colnames(Tea_type_only_F1) <- recode(colnames(Tea_type_only_F1) ,
                                     "CODE_F1" = "ID",
                                     "是否喝茶_F1" = "Tea_drink",
                                     "常饮用茶种类_F1"  = "Tea_type")
Tea_type_only_F1$Times <- "F1"


Tea_type_only_F2 <- Tea_F2_filter[,c("CODE_F2", "过去一年是否经常喝茶_F2", "常饮用茶种类_F2")]
colnames(Tea_type_only_F2) <- recode(colnames(Tea_type_only_F2) ,
                                     "CODE_F2" = "ID",
                                     "过去一年是否经常喝茶_F2" = "Tea_drink",
                                     "常饮用茶种类_F2"  = "Tea_type")
Tea_type_only_F2$Times <- "F2"


Tea_type_only_F3 <- Tea_F3_filter[,c("CODE_F3", "是否喝茶_F3", "常饮用茶种类_F3")]
colnames(Tea_type_only_F3) <- recode(colnames(Tea_type_only_F3) ,
                                     "CODE_F3" = "ID",
                                     "是否喝茶_F3" = "Tea_drink",
                                     "常饮用茶种类_F3"  = "Tea_type")
Tea_type_only_F3$Times <- "F3"
#***********合并
Tea_type_only <- rbind(Tea_type_only_F0, Tea_type_only_F1, Tea_type_only_F2, Tea_type_only_F3)
Tea_type_only <- Tea_type_only %>%
  mutate(
    ID = if_else(
      str_starts(ID, "F"),
      str_sub(ID, 3),
      ID
    )
  )

Tea_type_single <- Tea_type_only %>%
  mutate(
    Tea_type_raw = as.character(Tea_type),
    
    # 去掉首尾空格、合并连续空格
    Tea_type_clean = str_squish(Tea_type_raw),
    Tea_type_clean = na_if(Tea_type_clean, ""),
    
    # 把各种可能的分隔符统一成中文分号
    Tea_type_clean = str_replace_all(
      Tea_type_clean,
      "\\s*(；|;|，|,|、|/|／|\\+|＋|&|和|与)\\s*",
      "；"
    ),
    
    # 去掉开头或结尾多余分号
    Tea_type_clean = str_replace_all(Tea_type_clean, "^；+|；+$", ""),
    Tea_type_clean = str_replace_all(Tea_type_clean, "；+", "；"),
    
    # 计算茶种数量
    n_tea_type = case_when(
      is.na(Tea_type_clean) | Tea_type_clean == "" ~ 0L,
      TRUE ~ str_count(Tea_type_clean, "；") + 1L
    )
  ) %>%
  filter(
    Tea_drink == 0 |
      (
        Tea_drink == 1 &
          !is.na(Tea_type_clean) &
          Tea_type_clean != "" &
          n_tea_type == 1
      )
  ) %>%
  mutate(
    Tea_type = if_else(
      Tea_drink == 0,
      "Non-drinker",
      Tea_type_clean
    )
  )


table(Tea_type_single$Tea_type)
Tea_type_count <- as.data.frame(table(Tea_type_single$Tea_type)) %>%
  dplyr::rename(
    Tea_type = Var1,
    N = Freq
  )



Tea_type_recode <- Tea_type_single %>%
  mutate(
    Tea_type_clean = str_squish(Tea_type),
    
    Tea_type_group = case_when(
      Tea_drink == 0 ~ "Non-drinker",
      
      Tea_drink == 1 & Tea_type_clean %in% c("绿茶") ~ "Green tea",
      
      Tea_drink == 1 & Tea_type_clean %in% c("乌龙茶") ~ "Oolong tea",
      
      Tea_drink == 1 & Tea_type_clean %in% c("红茶") ~ "Black tea",
      
      Tea_drink == 1 & Tea_type_clean %in% c("普洱", "普洱茶", "普洱.","黑茶", "藏茶") ~ "Pu-erh/dark tea",
      
      Tea_drink == 1 ~ "Other tea",
      
      TRUE ~ NA_character_
    ),
    
    Tea_type_group = factor(
      Tea_type_group,
      levels = c(
        "Non-drinker",
        "Green tea",
        "Oolong tea",
        "Black tea",
        "Pu-erh/dark tea",
        "Other tea"
      )
    )
  )
table(Tea_type_recode$Tea_type_group)
####lmer####
# data_type_long <- merge(ASM_Height_long, Tea_type_recode, by.x = c("ID", "Followup") , by.y = c("ID", "Times"))
# 
# data_list <- list(Cov_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], data_type_long)
# data_type_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
#   dplyr::filter(stats::complete.cases(.))

data_list <- list(Cov_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], Tea_type_recode[Tea_type_recode$Times == "F0",], ASM_Height_long)
data_type_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
# 合并握力等，要用left_join
data_type_v2 <- left_join(data_type_v1, SPPB_grip_long, by = c("ID","Followup"))
table(SPPB_grip_long$Followup)

data_type_final <- scale_columns_group(data_type_v2, c("ASM","ASMI","grip_max"), "Sex_F0")
table(data_type_final$Tea_type_group)
a <- unique(data_type_final$ID)

Lmer_results_type <- process_lmer(c("Tea_type_group"),
                                 c("ASM_z","ASMI_z", "grip_max_z"), #,"SPPB_total"
                                 data_type_final, 
                                 Covariates_all_lmer, 
                                 2)
####正文写作
a <- Lmer_results_type[Lmer_results_type$P_value <0.05,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low,3)
a$CI_high <- round(a$CI_high,3)

max(a$Estimate)
min(a$Estimate)
max(a$CI_high)
min(a$CI_low)
####可视化####
# 1. 预处理数据：锁定自变量（Level）的上下显示顺序，避免被默认按字母排序
Lmer_results_type$Level <- factor(
  Lmer_results_type$Level, 
  levels = rev(c("Green tea", "Pu-erh/dark tea", "Oolong tea", "Black tea", "Other tea"))
)
Lmer_results_type$Outcome <- recode(Lmer_results_type$Outcome,
                                    "ASM_z"  = "ASM" ,
                                    "ASMI_z" = "ASMI",
                                    "grip_max_z" = "Handgrip strength"
)

p <- plot_forest_simple(
  data = Lmer_results_type,
  x = "Estimate",
  y = "Level",
  p_col = "P_value",
  facet = "Outcome",
  ci_low = "CI_low",
  ci_high = "CI_high",
  x_lab = "β (95% CI)",
  y_lab = "Tea Type Group",
  output_file = "Tea_type_forest_plot.png",
  width = 8,
  height = 4
)

####*******propensity score weighting*******####
####计算权重####
run_weighted_lmer_multicategory <- function(exposure, outcome, data, covariates) {
  
  # 1. 用每个 ID 的基线信息估计 generalized propensity score
  ps_data <- data %>%
    dplyr::select(ID, all_of(exposure), all_of(covariates)) %>%
    distinct(ID, .keep_all = TRUE) %>%
    drop_na()
  
  ps_data[[exposure]] <- as.factor(ps_data[[exposure]])
  
  # 如果暴露少于两个水平，跳过
  if (nlevels(ps_data[[exposure]]) < 2) {
    return(NULL)
  }
  
  # 2. generalized propensity score-based overlap weighting
  ps_formula <- as.formula(
    paste0(
      exposure,
      " ~ ",
      paste(covariates, collapse = " + ")
    )
  )
  
  w.out <- weightit(
    formula = ps_formula,
    data = ps_data,
    method = "ps",
    estimand = "ATO"
  )
  
  ps_data$overlap_weight <- get.w(w.out)
  
  # 3. 合并权重回纵向数据
  data_w <- data %>%
    inner_join(
      ps_data %>% dplyr::select(ID, overlap_weight),
      by = "ID"
    ) %>%
    filter(!is.na(.data[[outcome]]))
  
  # 4. 加权 LMM
  lmer_formula <- as.formula(
    paste0(
      outcome,
      " ~ ",
      exposure,
      " + ",
      paste(covariates, collapse = " + "),
      " + (1 | ID)"
    )
  )
  
  fit <- lmer(
    lmer_formula,
    data = data_w,
    weights = overlap_weight,
    REML = FALSE
  )
  
  # 5. 提取结果
  result <- broom.mixed::tidy(
    fit,
    effects = "fixed",
    conf.int = TRUE
  ) %>%
    filter(grepl(exposure, term)) %>%
    mutate(
      Outcome = outcome,
      Predictor = exposure,
      Model = "Generalized overlap weighting",
      N = n_distinct(data_w$ID)
    ) %>%
    dplyr::select(
      Outcome,
      Predictor,
      Model,
      N,
      term,
      estimate,
      conf.low,
      conf.high,
      p.value
    )
  
  return(result)
}

Exposure_vars <- c(
  "tea_freq_group",
  "serum_I_catechin_F0_T",
  "serum_I_epicatechin_F0_T",
  "serum_I_EGC_F0_T",
  "serum_I_EGCG_F0_T",
  "serum_I_ECG_F0_T"
)

Outcome_vars <- c(
  "ASM_z",
  "ASMI_z",
  "grip_max_z"
)

Weighted_Lmer_results <- map_dfr(
  Exposure_vars,
  function(x) {
    map_dfr(
      Outcome_vars,
      function(y) {
        run_weighted_lmer_multicategory(
          exposure = x,
          outcome = y,
          data = data1_final,
          covariates = Covariates_all_lmer
        )
      }
    )
  }
)
####可视化####
final_table_data <- Weighted_Lmer_results %>%
  mutate(
    suffix = str_extract(term, "[12]$"),
    
    Predictor_Group = case_when(
      Predictor == "tea_freq_group"             ~ "Tea consumption",
      Predictor == "serum_I_catechin_F0_T"      ~ "Catechin",
      Predictor == "serum_I_epicatechin_F0_T"   ~ "Epicatechin",
      Predictor == "serum_I_EGC_F0_T"           ~ "Epigallocatechin",
      Predictor == "serum_I_ECG_F0_T"           ~ "Epicatechin gallate",
      Predictor == "serum_I_EGCG_F0_T"          ~ "Epigallocatechin gallate",
      TRUE ~ Predictor
    ),
    
    # 关键：指定 Predictor 分面从左到右顺序
    Predictor_Group = factor(
      Predictor_Group,
      levels = c(
        "Tea consumption",
        "Catechin",
        "Epicatechin gallate",
        "Epigallocatechin gallate",
        "Epicatechin",
        "Epigallocatechin"
      )
    ),
    
    Level_Detail = case_when(
      Predictor == "tea_freq_group" & suffix == "1" ~ "  1-6 times/week",
      Predictor == "tea_freq_group" & suffix == "2" ~ "  ≥7 times/week",
      Predictor %in% c("serum_I_catechin_F0_T", "serum_I_EGC_F0_T") & suffix == "1" ~ "  Low",
      Predictor %in% c("serum_I_catechin_F0_T", "serum_I_EGC_F0_T") & suffix == "2" ~ "  High",
      suffix == "1" ~ "  Tertile 2",
      suffix == "2" ~ "  Tertile 3",
      TRUE ~ term
    ),
    
    Outcome_clean = case_when(
      Outcome == "ASM_z"        ~ "ASM",
      Outcome == "ASMI_z"       ~ "ASMI",
      Outcome == "grip_max_z"   ~ "Handgrip strength",
      TRUE ~ Outcome
    ),
    
    Outcome_clean = factor(
      Outcome_clean,
      levels = c(
        "ASM",
        "ASMI",
        "Handgrip strength"
      )
    ),
    
    Significance = if_else(p.value < 0.05, "Significant", "Not Significant")
  ) %>%
  arrange(Predictor_Group, suffix) %>%
  mutate(
    Level_Detail = factor(Level_Detail, levels = rev(unique(Level_Detail)))
  )

# 2. 绘制【表格级·纵向通栏】森林图
p_premium <- ggplot(final_table_data, aes(x = estimate, y = Level_Detail)) +
  
  # 后台现代优雅斑马线（提供绝佳的视线轨道）
  geom_hline(yintercept = seq_along(unique(final_table_data$Level_Detail)), color = "#F8F9FA", size = 7) +
  
  # 统一且神圣的 0 刻度无效线
  geom_vline(xintercept = 0, linetype = "solid", color = "#BDC3C7", size = 0.6) +
  
  # 绘制置信区间（不再使用 position_dodge，每条线拥有独立的清爽轨道）
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = Significance), height = 0.15, size = 0.8) +
  
  # 绘制效应量实心圆点
  geom_point(aes(color = Significance), size = 2.5) +
  
  # 灵魂一步：按 Outcome 进行【纵向分面】(nrow = 4)，并且让 Predictor_Group 作为侧边分类！
  facet_grid(
    Outcome_clean ~ Predictor_Group,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      Predictor_Group = label_wrap_gen(width = 18)
    )
  )+  
  # 配色：学术经典深红（显著）+ 庄重石墨灰（不显著），拒绝廉价花哨感
  scale_color_manual(values = c("Significant" = "#B22222", "Not Significant" = "#5A5A5A")) +
  
  # 顶级期刊主题定制
  theme_minimal(base_size = 13) +
  labs(
    x = "β (95% CI)",
    y = ""
  ) +
  theme(
    # 隐藏一切多余的网格线、背景和框线
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#EAEDED", linetype = "dashed", size = 0.5),
    panel.spacing.y = unit(0.5, "lines"),  # 结局之间的间距
    panel.spacing.x = unit(1, "lines"),    # 自变量大类之间的物理留白
    
    # 彻底美化分面条（让表格的大分类看起来像完美的排版加粗文本）
    strip.background = element_rect(fill = "#F2F4F4", color = NA),
    strip.text.x = element_text(face = "bold", size = 13, color = "#2C3E50"),
    strip.text.y = element_text(axis.text.y = element_text(angle = 0), face = "bold", size = 13),
    
    # 坐标轴与文字
    axis.line.x = element_line(color = "#BDC3C7", size = 0.6),
    axis.text.y = element_text(face = "plain", size = 13, color = "#2C3E50"), 
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(face = "bold", size = 14,  margin = ggplot2::margin(t = 10)),
    
    # 隐藏自带的标签图例（因为颜色已经代表了显著性，无需多余图例占用空间）
    legend.position = "top"
  )

# 3. 保存为极限高分辨率的长图（根据内容给足高度，消除拥挤）

ggsave("ultimate_journal_forest_plot.png", plot = p_premium, width = 15, height = 8, dpi = 600)
####混杂是否平衡####
# 变量名标签
var_labels <- c(
  "Age" = "Age (years)",
  "Sex_F0" = "Sex",
  "Education_F0_1" = "Education: junior high school or below",
  "Education_F0_2" = "Education: senior high school/technical school",
  "Education_F0_3" = "Education: college degree or above",
  "Income_F0_0" = "Monthly income: ≤500",
  "Income_F0_1" = "Monthly income: 500–1500",
  "Income_F0_2" = "Monthly income: 1500–3000",
  "Income_F0_3" = "Monthly income: ≥3000",
  "WBTOT_FAT" = "Whole body fat (kg)",
  "Met_mean" = "Physical activity",
  "Protein_mean" = "Dietary protein intake",
  "Total_score" = "China Healthy Eating Index score",
  "coffee" = "Coffee drinker",
  "Smoke_F0" = "Smoker",
  "Alcohol_F0" = "Alcohol drinker",
  "Calcium_F0" = "Use of calcium supplements",
  "Vitamin_F0" = "Use of vitamin supplements",
  "Oil_sup_follow" = "Use of fish oil supplements",
  "Disease_F0" = "Disease status",
  "Estrogen_F0" = "Estrogen therapy",
  "Menopause_F0" = "Menopause"
)

p_love_tea <- love.plot(
  w_tea,
  stats = "mean.diffs",
  abs = TRUE,
  threshold = 0.1,
  var.order = "unadjusted",
  agg.fun = "max",
  var.names = var_labels,
  sample.names = c("Before weighting", "After weighting"),
  colors = c("grey65", "#1F4E79"),
  shapes = c(16, 17),
  size = 3
) +
  labs(
    title = NULL,
    x = "Maximum absolute standardized mean difference",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title = element_blank()
  )


ggsave(
  "Supplementary_Figure_Balance_Tea_Frequency_Overlap_Weighting.png",
  p_love_tea,
  width = 8,
  height = 6.5,
  dpi = 600
)
####*******未来茶摄入预测之前的骨骼肌量*******####
# 提取未来茶暴露，例如 F3 tea frequency
future_tea <- data1_final %>%
  filter(Followup == "F3") %>%
  dplyr::select(ID, future_tea_freq_group = tea_freq_group) %>%
  distinct(ID, .keep_all = TRUE)

# 合并到早期结局数据，例如 F1-F2
data_negative_control <- data1_final %>%
  filter(Followup %in% c("F1", "F2")) %>%
  left_join(future_tea, by = "ID") %>%
  filter(!is.na(future_tea_freq_group))

data_negative_control$future_tea_freq_group <- as.factor(data_negative_control$future_tea_freq_group)

# 用 future tea exposure 预测 earlier muscle outcomes
Negative_control_results <- process_lmer(
  c("future_tea_freq_group"),
  c("ASM_z", "ASMI_z", "grip_max_z", "gait_speed_z"),
  data_negative_control,
  Covariates_all_lmer,
  2
)
a <- unique(data_negative_control$ID)

Negative_control_results <- Negative_control_results %>%
  mutate(
    Level = recode(
      Level,
      "1" = "1-6 times/week",
      "2" = "≥ 7 times/week"
    ),
    Level = factor(
      Level,
      levels = rev(c("1-6 times/week", "≥ 7 times/week"))
    )
  )

Negative_control_results$Outcome <- recode(Negative_control_results$Outcome,
                                           "ASM_z"  = "ASM" ,
                                           "ASMI_z" = "ASMI"
)

p <- plot_forest_simple(
  data = Negative_control_results,
  x = "Estimate",
  y = "Level",
  p_col = "P_value",
  facet = "Outcome",
  ci_low = "CI_low",
  ci_high = "CI_high",
  x_lab = "β (95% CI)",
  y_lab = "Tea consumption",
  output_file = "future tea.png",
  width = 6,
  height = 3
)

####*******E-value计算*******####
calc_evalue_from_smd <- function(beta, ci_low, ci_high) {
  
  # 定义一个内部函数：SMD -> RR -> E-value
  smd_to_evalue <- function(d) {
    rr <- exp(0.91 * d)
    ifelse(
      d == 0,
      1,
      rr + sqrt(rr * (rr - 1))
    )
  }
  
  # 1. Point estimate 的 E-value
  d_est <- abs(beta)
  evalue_est <- smd_to_evalue(d_est)
  
  # 2. CI limit closest to null
  d_ci_closest <- case_when(
    ci_low <= 0 & ci_high >= 0 ~ 0,
    ci_low > 0 ~ abs(ci_low),
    ci_high < 0 ~ abs(ci_high),
    TRUE ~ 0
  )
  evalue_ci_closest <- smd_to_evalue(d_ci_closest)
  
  # 3. 直接用 CI_high 计算 E-value
  d_ci_high <- abs(ci_high)
  evalue_ci_high <- smd_to_evalue(d_ci_high)
  
  # 4. CI 中离 null 最远的边界
  d_ci_farthest <- max(abs(ci_low), abs(ci_high), na.rm = TRUE)
  evalue_ci_farthest <- smd_to_evalue(d_ci_farthest)
  
  return(
    data.frame(
      E_value_estimate = evalue_est,
      E_value_CI_closest_null = evalue_ci_closest,
      E_value_CI_high = evalue_ci_high,
      E_value_CI_farthest_null = evalue_ci_farthest
    )
  )
}

Evalue_results <- Lmer_results1_all %>%
  #filter(Predictor == "tea_freq_group", P_value < 0.05) %>%
  filter(Predictor == "tea_freq_group") %>%
  rowwise() %>%
  mutate(tmp = list(calc_evalue_from_smd(Estimate, CI_low, CI_high))) %>%
  tidyr::unnest_wider(tmp) %>%
  ungroup() %>%
  mutate(
    across(
      starts_with("E_value"),
      ~ round(.x, 2)
    )
  )

Evalue_results$Outcome <- recode(Evalue_results$Outcome,
                                 "ASM_z"  = "ASM" ,
                                 "ASMI_z" = "ASMI",
                                 "grip_max_z" = "Handgrip strength",
                                 "gait_speed_z" = "Walking speed"
)

Evalue_results$Predictor <- recode(Evalue_results$Predictor,
                                   "tea_freq_group" = "Tea consumption")

Evalue_results$Level <- recode(Evalue_results$Level,
                               "1" = "1-6 times/week",
                               "2" = "≥ 7 times/week"
)

fmt3 <- function(x) {
  sprintf("%.3f", as.numeric(x))
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

Evalue_table <- Evalue_results %>%
  #filter(Significance == "Significant") %>%
  transmute(
    Outcome,
    Predictor,
    Level = Level,
    
    `β (95% CI)` = paste0(
      fmt3(Estimate),
      " (",
      fmt3(CI_low),
      ", ",
      fmt3(CI_high),
      ")"
    ),
    
    `p` = fmt_p(P_value),
    
    `E-value (95% CI)` = paste0(
      fmt3(E_value_estimate),
      " (",
      fmt3(E_value_CI_closest_null),
      ", ",
      fmt3(E_value_CI_high),
      ")"
    )
  )

writexl::write_xlsx(Evalue_table, "Evalue_table.xlsx")
####**********************其他**********************####
####年龄分布图和性别饼图####
data_ga_all <- data1_final %>%
  group_by(ID) %>%
  dplyr::slice(1) %>%
  ungroup()
#********************年龄分布图
p_age <- ggplot(data_ga_all, aes(x = Age_F0)) +
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
data_ga_all$Sex_F0 <- recode(data_ga_all$Sex_F0,
                             "0" = "Female",
                             "1" = "Male")
df_sex <- data_ga_all %>%
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
Follow_up_time <- Age_wide[Age_wide$ID %in% data_ga_all$ID,]
Follow_up_time <- merge(Follow_up_time, Cov_F0_all[,c("ID","Age_F0")], by = "ID")

Follow_up_time_filtered <- Follow_up_time %>%
  filter(!grepl("^NL4", ID))
Follow_up_time_filtered$F01 <- Follow_up_time_filtered$F1age-Follow_up_time_filtered$Age_F0
Follow_up_time$F12 <- Follow_up_time$F2age-Follow_up_time$F1age
Follow_up_time_filtered$F12 <- Follow_up_time_filtered$F2age-Follow_up_time_filtered$F1age

Follow_up_time$F23 <- Follow_up_time$F3age-Follow_up_time$F2age
Follow_up_time$F34 <- Follow_up_time$F4age-Follow_up_time$F3age

round(mean(Follow_up_time_filtered$F01, na.rm = TRUE),1)
round(mean(Follow_up_time$F12, na.rm = TRUE),1)
round(mean(Follow_up_time$F23, na.rm = TRUE),1)
round(mean(Follow_up_time$F34, na.rm = TRUE),1)
round(mean(Follow_up_time_filtered$F12, na.rm = TRUE),1)


Follow_up_time$F02 <- Follow_up_time$F2age-Follow_up_time$Age_F0
Follow_up_time$F03 <- Follow_up_time$F3age-Follow_up_time$Age_F0
Follow_up_time$F04 <- Follow_up_time$F4age-Follow_up_time$Age_F0
round(mean(Follow_up_time_filtered$F01, na.rm = TRUE),1)
round(mean(Follow_up_time$F02, na.rm = TRUE),1)
round(mean(Follow_up_time$F03, na.rm = TRUE),1)
round(mean(Follow_up_time$F04, na.rm = TRUE),1)
####*****************************data7-饮茶&菌群&骨骼肌*****************************####
####数据集处理合并####
Micro_ASMI <- merge(Micro_selected_filter_clr, ASM_Height_long, by.x = c("ID","Times"), by.y = c("ID","Followup"))
Micro_ASMI <- Micro_ASMI[Micro_ASMI$Times %in% c("F2","F3"),]

data_list <- list(Cov_F0_final, Tea_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], Micro_ASMI) #
data7_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
species_names <- grep("^s__.*_clr$", colnames(data7_v1), value = TRUE)
# 合并握力等
data7_v2 <- left_join(data7_v1,SPPB_grip_long,by = c("ID" = "ID","Times" = "Followup"))
data7_v2 <- left_join(data7_v2, Catechin, by = c("ID" = "ID"))

data7_v3 <- group_variables3(
  data = data7_v2,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(data7_v3)%in% factor_name)
for(i in idx ){
  data7_v3[[i]] <-  as.factor(data7_v3[[i]])
}
#***************ASMI根据性别Z分数，菌群Z分数
data7_v4 <- scale_columns_group(data7_v3, c("ASM", "ASMI", "grip_max"), c("Sex_F0"))
data7_v5 <- scale_columns(data7_v4, c("gait_speed","chair_time"))
data7_final <- scale_columns(data7_v5, species_names)
table(data7_final$Times)

a <- unique(data7_final$ID)
a <- unique(data7_final[data7_final$tea_freq_group %in% c(0,2),]$ID)
a <- data7_final[,c("ID","Times","tea_freq_group","ASM","s__Tyzzerella_nexilis_clr_z")]
a <- a[complete.cases(a),]
a1 <- unique(a$ID)
####简单调整####
Lmer_results7_all <- process_lmer(c("tea_freq_group"), 
                                  paste0(species_names,"_z"),
                                  data7_final,
                                  c("Age","Sex_F0"), #
                                  2)
Lmer_results7_all <- Lmer_results7_all[Lmer_results7_all$Level == 2,]
Lmer_results7_all$P_FDR <- p.adjust(Lmer_results7_all$P_value, method = "fdr")
Lmer_results7_all_FDR_sig <- Lmer_results7_all[Lmer_results7_all$P_FDR <0.05,]
a <- unique(data7_final[data7_final$tea_freq_group %in% c(0, 2),]$ID)

# 数据清洗和转换（不使用函数，直接在 mutate 中判断）
Lmer_results7_all_output <- Lmer_results7_all %>%
  # 只保留需要的列
  dplyr::select(Outcome, Predictor, Level, Estimate, CI_low, CI_high, P_value, P_FDR) %>%
  # 清洗 Outcome 列
  mutate(
    Outcome = str_remove_all(Outcome, "s__"),
    Outcome = str_remove_all(Outcome, "_clr_z"),
    Outcome = str_replace_all(Outcome, "_", " "),
    Outcome = str_to_sentence(Outcome),
    
    # 根据原始 Predictor 名称转换 Level 标签
    Level_Label = case_when(
      # tea_freq_group 的处理
      Predictor == "tea_freq_group" ~ "≥ 7 times/week vs. Non-drinker",
      # epicatechin 和 EGC 的处理
      Predictor %in% c("serum_I_epicatechin_F0_T", "serum_I_EGC_F0_T") ~ "High vs. Undetectable",
      # catechin, EGCG, ECG 的处理
      Predictor %in% c("serum_I_catechin_F0_T", "serum_I_EGCG_F0_T", "serum_I_ECG_F0_T") ~ "Tertile 3 vs. Tertile 1",
      TRUE ~ as.character(Level)
    ),
    
    # 重命名 Predictor（在 Level 转换之后进行）
    Predictor = recode(Predictor, !!!predictor_rename),
    
    # 构建 β (95% CI)
    `β (95% CI)` = sprintf("%.3f (%.3f, %.3f)", Estimate, CI_low, CI_high),
    
    # 格式化 P 值
    P_value_formatted = case_when(
      P_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", P_value)
    ),
    P_FDR_formatted = case_when(
      P_FDR < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", P_FDR)
    )
  ) %>%
  # 重新排列列顺序
  dplyr::select(Outcome, Predictor, Level = Level_Label, `β (95% CI)`, 
         p = P_value_formatted, 
         p_FDR = P_FDR_formatted) %>%
  # 按 Predictor 和 Outcome 排序
  arrange(Predictor, Outcome)
####可视化####
library(dplyr)
library(ggplot2)
library(forcats)
plot_df <- Lmer_results7_all_FDR_sig %>%
  mutate(
    Outcome = str_replace_all(Outcome, "^s__|_clr_z$", ""),
    Outcome = str_replace_all(Outcome, "_", " ")
  ) %>%
  arrange(Estimate) %>%
  mutate(
    Outcome = factor(Outcome, levels = Outcome)
  )

p <- ggplot(plot_df, aes(y = Outcome, x = Estimate)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey70") +
  geom_segment(aes(yend = Outcome, x = 0, xend = Estimate, color = Estimate),
               linewidth = 1.5, alpha = 0.9) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high),
                 height = 0.16, linewidth = 0.75, color = "grey40") +
  geom_point(aes(fill = Estimate), shape = 21, size = 4.5, color = "black", stroke = 0.3) +
  scale_fill_gradient2(low = "#3C78D8", mid = "white", high = "#D65F5F", midpoint = 0, guide = "none") +
  scale_color_gradient2(low = "#3C78D8", mid = "white", high = "#D65F5F", midpoint = 0, guide = "none") +
  labs(y = NULL, x = "β (95% CI)") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "italic", size = 10, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = ggplot2::margin(10, 15, 10, 10)
  )

ggsave("Microbiome_Lollipop_Horizontal.png", p, width = 7, height = 14, dpi = 700)
####lmer-全调整####
Lmer_results7_all2 <- process_lmer(c("tea_freq_group", "serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T"),
                                   #unique(Lmer_results7_all_ASM_sig$Predictor),
                                   Lmer_results7_all_FDR_sig$Outcome,
                                   data7_final,
                                   c(setdiff(Covariates_all_lmer,c("Protein_mean")),c("Fiber_soluble_mean")),
                                   2)
Lmer_results7_all2 <- Lmer_results7_all2[Lmer_results7_all2$Level == 2,]
Lmer_results7_all2 <- Lmer_results7_all2 %>%
  group_by(Predictor) %>%
  mutate(
    P_FDR = p.adjust(P_value, method = "fdr")
  ) %>%
  ungroup()

Lmer_results7_all2_sig <- Lmer_results7_all2[Lmer_results7_all2$P_value < 0.05 & Lmer_results7_all2$P_FDR < 0.05,]

a1 <- Lmer_results7_all2_sig[Lmer_results7_all2_sig$Outcome == "s__Bacteroides_fragilis_clr_z",]
#***************正文写作
a <- Lmer_results7_all2_sig[Lmer_results7_all2_sig$Predictor == "tea_freq_group" & Lmer_results7_all2_sig$Estimate >0,]
a$Estimate <- round(a$Estimate, 3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)


a <- Lmer_results7_all2_sig[Lmer_results7_all2_sig$Predictor == "tea_freq_group" & Lmer_results7_all2_sig$Estimate <0,]
a$Estimate <- round(a$Estimate, 3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)
####验证lmer####
x <- "tea_freq_group"
y  <- "s__Bacteroides_fragilis_clr_z"

formula_text <- paste(
  y, "~",
  x, "+",
  paste(c(setdiff(Covariates_all_lmer,c("Protein_mean")),c("Fiber_soluble_mean")), collapse = " + "),
  "+ (1|ID)"
)
formula_full <- as.formula(formula_text)

model <- lmer(formula_full,data = data7_final)
summary(model)
####可视化####
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)

# 1. 显著性标签与色彩降阶计算
df_plot <- Lmer_results7_all2 %>%
  
  mutate(Predictor = ifelse(Predictor %in% names(predictor_rename),
                            predictor_rename[Predictor],
                            Predictor),
         Predictor = factor(Predictor,levels = predictor_rename)
  ) %>%
  
  mutate(
    # ① 清洗丑陋的菌群名后缀
    Species_clean = str_replace_all(Outcome, "^s__|_clr_z$", ""),
    Species_clean = str_replace_all(Species_clean, "_", " "),
    
    # ② 判定显著性星号（不显著的直接留空 ""，保持图面极简清爽）
    Sig_Label = case_when(
      P_FDR < 0.001 ~ "***",
      P_FDR < 0.01  ~ "**",
      P_FDR < 0.05  ~ "*",
      TRUE ~ ""
    ),
    
    # ③ 【核心新颖视觉】：将 Direction 与显著性融合
    # 不显著的统一归为 "Not Significant" 组，从而在后面单独对其剥离色彩、进行高级灰色降阶
    Sig_Group = case_when(
      P_FDR >= 0.05 ~ "Not Significant",
      Estimate > 0    ~ "Positive (Sig)",
      TRUE            ~ "Negative (Sig)"
    )
  ) %>%
  # ④ 精密排序：按 Estimate 排序，形成完美的对角线阶梯趋势
  mutate(Species_clean = fct_reorder(Species_clean, Estimate))

# 提取排序后的菌名作为因子能拿到的数字，用于精准圈定底纹跑道范围
y_indices <- seq_along(levels(df_plot$Species_clean))

# 2. 绘制新颖的“色彩过滤”底纹森林图
p_trendy_forest_sig <- ggplot(df_plot, aes(x = Estimate, y = Species_clean)) +
  
  # 【新颖核心】A. 动态显著性跑道：仅为显著的菌群铺设优雅底纹，不显著的自动留白，视觉反差感极强
  # 这里将 alpha 放入映射，实现不显著的菌底纹变为 0 (完全透明)
  geom_rect(aes(ymin = as.numeric(Species_clean) - 0.46, 
                ymax = as.numeric(Species_clean) + 0.46, 
                xmin = -Inf, xmax = Inf,
                fill = Sig_Group), 
            alpha = ifelse(df_plot$Sig_Group == "Not Significant", 0, 0.07), 
            color = NA) +
  
  # B. 纵向对齐清爽虚线
  geom_vline(xintercept = 0, color = "#94A3B8", linetype = "dashed", linewidth = 0.6) +
  
  # C. 横向极细点引线：若隐若现地指引视线
  geom_segment(aes(x = -Inf, xend = CI_low, y = Species_clean, yend = Species_clean), 
               color = "#E2E8F0", linetype = "dotted", linewidth = 0.5) +
  
  # D. 悬浮置信区间粗线：不显著的菌会通过色彩控制自动降级为高级透明灰
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, color = Sig_Group), 
                 height = 0, linewidth = 1.6, alpha = 0.9, lineend = "round") +
  
  # E. 效应量同心圆指针设计
  # 外圈包裹层
  geom_point(aes(color = Sig_Group), size = 3.6, shape = 21, fill = "white", stroke = 1.5) +
  # 内圈实心层
  geom_point(aes(color = Sig_Group), size = 1.6, shape = 16) +
  
  # F. 【高级画龙点睛】：在同心圆指针正上方稍微悬浮印上星号（vjust=-0.6 微调距离）
  # 颜色同样跟着 Sig_Group 走，显著的呈现标志色，不显著的隐藏
  # geom_text(aes(label = Sig_Label, color = Sig_Group), 
  #           vjust = -0.4, hjust = 0.5, size = 4.5, fontface = "bold", show.legend = FALSE) +
  
  # G. 根据 Predictor 变量进行多分面排版
  facet_grid(. ~ Predictor, scales = "free_x") +
  
  # H. 顶刊大师级色盘映射控制
  scale_fill_manual(values = c(
    "Negative (Sig)"  = "#3B82F6",   # 显著正关联：科技蓝底纹
    "Positive (Sig)" = "#EF4444",   # 显著负关联：复古红底纹 
    "Not Significant" = "#FFFFFF"   # 不显著：纯白底纹（等同于无底纹）
  ), guide = "none") +
  
  scale_color_manual(values = c(
    "Negative (Sig)" = "#1D4ED8",   # 显著正关联：深海蓝指针
    "Positive (Sig)"  = "#B91C1C",   # 显著负关联：胭脂红指针 
    "Not Significant" = "#CBD5E1"   # 不显著：高级消隐雾面灰（降低存在感）
  )) +
  
  # I. 极致无噪轻量化学术主题
  theme_minimal() + #base_family = "sans"
  theme(
    # 左侧菌名：根据显著性动态调整，显著的加粗，不显著的稍微减淡，对比拉满
    axis.text.y = element_text(
      size = 14, 
      face = "italic",
      color = "#0F172A",
      hjust = 1
    ),
    axis.text.x = element_text(size = 12, color = "#475569", face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold", color = "#0F172A", margin = ggplot2::margin(t = 15)),
    axis.title.y = element_blank(),
    
    panel.grid = element_blank(),
    
    # 顶部分面卡片（Predictor 标签说明）
    strip.background = element_rect(fill = "#F1F5F9", color = NA),
    strip.text = element_text(size = 15, face = "bold", color = "#1E293B", margin = ggplot2::margin(t=10, b=10)),
    panel.spacing = unit(0.5, "lines"), 
    
    # 图例极简置顶：过滤掉不显著的灰色图例，只保留红蓝核心标签
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold", color = "#475569"),
    
    plot.margin = ggplot2::margin(25, 25, 25, 25)
  ) +
  labs(
    x = "β (95% CI)"
  )


# ggsave("Figure 3_1.pdf", 
#        plot = p_trendy_forest_sig, 
#        width = 19, 
#        height = 16)

ggsave("Figure 3.png",
       plot = p_trendy_forest_sig,
       width = 20,
       height = 20,
       dpi = 600)

####菌群&骨骼肌####
Lmer_results7_all_ASM <- process_lmer(unique(Lmer_results7_all2_sig$Outcome), 
                                      c("ASM_z","ASMI_z","grip_max_z"),
                                      data7_final,# %>% filter(tea_freq_group %in% c(0, 2)),
                                      c(Covariates_all_lmer,c("Fiber_soluble_mean", "tea_freq_group")), #
                                      2)
Lmer_results7_all_ASM <- Lmer_results7_all_ASM %>%
  group_by(Outcome) %>%
  mutate(
    P_FDR = p.adjust(P_value, method = "fdr")
  ) %>%
  ungroup()

Lmer_results7_all_ASM_sig <- Lmer_results7_all_ASM[Lmer_results7_all_ASM$P_value < 0.05 & Lmer_results7_all_ASM$P_FDR < 0.10,] #& Lmer_results7_all_ASM$P_FDR < 0.05
a <- unique(data7_final$ID)

#*************正文写作
a <- Lmer_results7_all_ASM_sig[Lmer_results7_all_ASM_sig$Estimate >0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)

a <- Lmer_results7_all_ASM_sig[Lmer_results7_all_ASM_sig$Estimate <0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)
####可视化####
plot_df <- Lmer_results7_all_ASM %>%
  
  mutate(
    
    Predictor = str_replace_all(Predictor, "^s__|_clr_z$", ""),
    Predictor = str_replace_all(Predictor, "_", " "),
    
    ## outcome
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    ),
    
    ## significance stars
    sig = case_when(
      P_value < 0.05 & P_FDR < 0.10 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  
  group_by(Outcome) %>%
  
  arrange(
    Estimate,
    .by_group = TRUE
  ) %>%
  
  mutate(
    
    Predictor =
      factor(
        Predictor,
        levels = unique(Predictor)
      )
  ) %>%
  
  ungroup()

p <- plot_forest_facet(
  data = plot_df,
  estimate_col = "Estimate",
  ci_low_col = "CI_low",
  ci_high_col = "CI_high",
  pathway_col = "Predictor", #自变量
  outcome_col = "Outcome", #分面
  p_col = "P_value",
  sig_rule = "p_fdr"
  
)

ggsave("Micro_ASM.png", plot = p, width = 13, height = 16, dpi = 500)
####模型效能的提高####
Lmer_model_compare_micro <- compare_lmer_predictors_by_pair(
  pair_df = Lmer_results7_all_ASM_sig,
  data = data7_final,
  covariates = c(Covariates_all_lmer, "Fiber_soluble_mean", "tea_freq_group"),
  predictor_col = "Predictor",
  outcome_col = "Outcome",
  predictor_type_col = "Variable_Type",
  id_var = "ID",
  predictor_rename = NULL,
  outcome_rename = c(
    "ASM_z" = "ASM",
    "ASMI_z" = "ASMI",
    "grip_max_z" = "Handgrip strength"
  ),
  output_file_raw = "Micro_ASM_model_performance_raw.csv",
  output_file_format = "Micro_ASM_model_performance_format.csv"
)

Lmer_model_compare_micro_output <- Lmer_model_compare_micro$format[,c("Outcome", "Predictor","Delta_Marginal_R2_percent","LRT_P_value")]
Lmer_model_compare_micro_output <- Lmer_model_compare_micro_output %>%
  dplyr::mutate(
    Predictor = str_remove(Predictor, "s__") %>% str_remove("_clr_z") %>% str_replace_all("_", " ")
    
  )
####中介效应分析-菌群####
Micro_sig_pos <- Lmer_results7_all2_sig[Lmer_results7_all2_sig$Estimate >0,] %>%
  
  inner_join(
    
    Lmer_results7_all_ASM_sig[Lmer_results7_all_ASM_sig$Estimate >0,] %>%
      
      dplyr::select(
        Predictor,
        Outcome
      ) %>%
      
      distinct(),
    
    by = c(
      "Outcome" = "Predictor"
    )
  ) 

Micro_sig_neg <- Lmer_results7_all2_sig[Lmer_results7_all2_sig$Estimate <0,] %>%
  
  inner_join(
    
    Lmer_results7_all_ASM_sig[Lmer_results7_all_ASM_sig$Estimate <0,] %>%
      
      dplyr::select(
        Predictor,
        Outcome
      ) %>%
      
      distinct(),
    
    by = c(
      "Outcome" = "Predictor"
    )
  )

Micro_sig <- rbind(Micro_sig_pos, Micro_sig_neg)

colnames(Micro_sig) <- recode(colnames(Micro_sig),
                              "Outcome.y" = "Outcome_muscle",
                              "Outcome" = "Mediator")


med_results_all_df_micro <- 
  
  run_microbiome_mediation_batch(
    
    data = data7_final,
    
    micro_sig = Micro_sig,
    
    covariates = c(
      Covariates_all_lmer,
      "Fiber_soluble_mean"
    ),
    
    mediator_col = "Mediator",
    
    outcome_col = "Outcome_muscle",
    
    predictor_col = "Predictor",
    
    exposure_high = 2,
    exposure_low = 0,
    
    sims = 1000,
    
    output_file =
      "mediation_results_all_micro.csv"
  )

med_results_all_df_micro$ACME_FDR <- p.adjust(med_results_all_df_micro$ACME_p, method = "fdr")
med_results_all_df_micro_sig <- med_results_all_df_micro[med_results_all_df_micro$ACME_FDR < 0.05,]
writexl::write_xlsx(med_results_all_df_micro, "mediation_results_all_micro.xlsx")
####验证中介效应####
Mediator <- "s__Gemmiger_formicilis_clr_z"
Exposure <- "tea_freq_group"
Outcome  <- "ASM_z"

covariates <- c(
  Covariates_all_lmer,
  "Fiber_soluble_mean"
)

df_med <- data7_final %>%
  
  filter(tea_freq_group %in% c(0, 2)) %>%
  
  mutate(
    Treat_bin = ifelse(tea_freq_group == 2, 1, 0)
  )

cov_string <- paste(
  paste0("`", covariates, "`"),
  collapse = " + "
)

formula_m <- as.formula(
  paste0(
    "`s__Gemmiger_formicilis_clr_z` ~ Treat_bin + ",
    cov_string,
    " + (1|ID)"
  )
)

model.m <- glmer(
  formula_m,
  data = df_med,
  family = gaussian()
)

formula_y <- as.formula(
  paste0(
    "`ASM_z` ~ Treat_bin + `s__Gemmiger_formicilis_clr_z` + ",
    cov_string,
    " + (1|`ID`)"
  )
)

model.y <- glmer(
  formula_y,
  data = df_med,
  family = gaussian()
)

set.seed(908)

med_res <- mediation::mediate(
  model.m,
  model.y,
  treat = "Treat_bin",
  mediator = "s__Gemmiger_formicilis_clr_z",
  sims = 1000
)
med_res$d0
med_res$d0.ci
med_res$d0.p
####十折交叉####
run_cv_stability <- function(data,
                             outcomes,
                             predictor,
                             covariates,
                             gender_label = 2,
                             k_fold = 10,
                             fdr_cutoff = 0.05,
                             freq_cutoff = 80,
                             seed = 2809) {
  
  library(dplyr)
  library(caret)
  
  set.seed(seed)
  
  ## 创建 participant-level folds
  
  ids_all <- unique(data$ID)
  
  folds <- createFolds(
    ids_all,
    k = k_fold,
    list = TRUE,
    returnTrain = FALSE
  )
  
  ## 保存每折结果
  cv_results <- list()
  
  ## 循环 CV
  
  for (i in seq_len(k_fold)) {
    
    cat("Processing fold:", i, "/", k_fold, "\n")
    
    ## test/train IDs
    test_ids <- ids_all[folds[[i]]]
    
    train_ids <- setdiff(
      ids_all,
      test_ids
    )
    
    ## train data
    train_data <- data[
      data$ID %in% train_ids,
    ]
    
    ## LMER
    res_train <- process_lmer(
      X = predictor,
      Y = outcomes,
      data = train_data,
      covariates = covariates,
      gender_label = gender_label
    )
    
    ## 保留 Level==2
    res_train <- res_train[
      res_train$Level == 2,
    ]
    
    ## FDR
    res_train$P_FDR <- p.adjust(
      res_train$P_value,
      method = "fdr"
    )
    
    ## 显著结果
    sig_train <- res_train[
      res_train$P_FDR < fdr_cutoff,
    ]
    
    ## 保存
    if (nrow(sig_train) > 0) {
      
      cv_results[[i]] <- data.frame(
        fold = i,
        feature = sig_train$Outcome,
        direction = sign(
          sig_train$Estimate
        )
      )
    }
  }
  
  ## 合并结果
  
  cv_df <- bind_rows(cv_results)
  
  ## stability
  stability_df <- cv_df %>%
    group_by(feature, direction) %>%
    summarise(
      frequency = n() / k_fold * 100,
      .groups = "drop"
    ) %>%
    arrange(desc(frequency))
  
  ## 高频 feature
  stability_sig <- stability_df[
    stability_df$frequency >= freq_cutoff,
  ]
  
  ## 返回
  return(list(
    folds = folds,
    raw = cv_df,
    stability = stability_df,
    stability_sig = stability_sig
  ))
}

cv_res_micro <- run_cv_stability(
  data = data7_final,
  outcomes = paste0(species_names,"_z"),
  predictor = "tea_freq_group",
  covariates = c("Age","Sex_F0"),
  gender_label = 2,
  k_fold = 10,
  fdr_cutoff = 0.05,
  freq_cutoff = 80
)

cv_res_results7_tea <- cv_res_micro$stability_sig
cv_res_results7_tea_clean <- cv_res_results7_tea %>%
  mutate(feature = gsub("^s__", "", feature),
         feature = gsub("_clr_z$", "", feature),
         feature = str_replace_all(feature, "_", " ")
  )

####lmer-全调整
Lmer_results7_all2_cv <- process_lmer(c("tea_freq_group"),
                                      cv_res_micro$stability_sig$feature,
                                      data7_final,
                                      c(setdiff(Covariates_all_lmer,c("Protein_mean")),c("Fiber_soluble_mean")),
                                      2)
Lmer_results7_all2_cv <- Lmer_results7_all2_cv[Lmer_results7_all2_cv$Level == 2,]
Lmer_results7_all2_cv$P_FDR <- p.adjust(Lmer_results7_all2_cv$P_value, method = "fdr")
Lmer_results7_all2_cv_sig <- Lmer_results7_all2_cv[Lmer_results7_all2_cv$P_FDR < 0.05,]
####菌群&骨骼肌
Lmer_results7_all_ASM_cv <- process_lmer(Lmer_results7_all2_cv_sig$Outcome, 
                                         c("ASM_z","ASMI_z","grip_max_z"),
                                         data7_final,# %>% filter(tea_freq_group %in% c(0, 2)),
                                         c(Covariates_all_lmer,c("Fiber_soluble_mean")), #
                                         2)
Lmer_results7_all_ASM_cv <- Lmer_results7_all_ASM_cv %>%
  group_by(Outcome) %>%
  mutate(
    P_FDR = p.adjust(P_value, method = "fdr")
  ) %>%
  ungroup()
Lmer_results7_all_ASM_cv_sig <- Lmer_results7_all_ASM_cv[Lmer_results7_all_ASM_cv$P_value < 0.05 & Lmer_results7_all_ASM_cv$P_FDR <0.10,]
####计算指数
pos_var <- intersect(Lmer_results7_all2_cv_sig[Lmer_results7_all2_cv_sig$Predictor == "tea_freq_group" & Lmer_results7_all2_cv_sig$Estimate > 0,]$Outcome, Lmer_results7_all_ASM_cv_sig$Predictor)
neg_var <- intersect(Lmer_results7_all2_cv_sig[Lmer_results7_all2_cv_sig$Predictor == "tea_freq_group" & Lmer_results7_all2_cv_sig$Estimate < 0,]$Outcome, Lmer_results7_all_ASM_cv_sig$Predictor)

Micro_GMSI_all_tea <- calculate_GMSI(
  data = data7_final,
  outcomes_pos = pos_var,
  outcomes_neg = neg_var,
  "GMI_Tea",
  time_var = "Times"
)
colnames(Micro_GMSI_all_tea)
####不同组间指数比较
GMSI_data <- merge(Micro_GMSI_all_tea[,c("ID", "Times", "GMI_Tea")], data7_final, by = c("ID","Times")) 
GMSI_data <- scale_columns(GMSI_data, c("GMI_Tea"))
table(GMSI_data[GMSI_data$Times == "F2",]$tea_freq_group)
table(GMSI_data[GMSI_data$Times == "F3",]$tea_freq_group)

GMSI_compare_tea <- process_lmer(c("tea_freq_group"),
                                 "GMI_Tea_z",
                                 GMSI_data %>% filter(tea_freq_group %in% c(0, 2)),
                                 c(setdiff(Covariates_all_lmer,c("Protein_mean")),c("Fiber_soluble_mean")),
                                 2)

GMI_tea <- plot_index_bar(
  data = GMSI_data,
  lmer_result = GMSI_compare_tea,
  index_var = "GMI_Tea_z",
  y_title = "Tea-Related GMI (Z-score) ",
  
)
ggsave("GMI_tea.png", plot = GMI_tea, width = 5, height = 5, dpi = 500)

####中介效应分析
# a1 <- process_lmer(c("GMI_Tea_z"), #"GMI_Tea","GMI_epicatechin", "GMI_ECG", "GMI_EGC", "GMI_EGCG"
#                    c("ASMI_z","ASM_z", "grip_max_z"),
#                    GMSI_data %>% filter(tea_freq_group %in% c(0, 2)), # 
#                    c(Covariates_all_lmer, c("Fiber_soluble_mean")),
#                    2)
exposures <- c("tea_freq_group")
mediators <- c("GMI_Tea_z")
outcomes <- c("ASM_z", "ASMI_z", "grip_max_z")

vars_needed <- c(
  "ID",
  exposures,
  mediators,
  outcomes,
  Covariates_all_lmer,
  "Fiber_soluble_mean"
)

data_med <- GMSI_data %>%
  filter(tea_freq_group %in% c(0, 2)) %>%
  filter(
    if_all(
      all_of(vars_needed),
      ~ !is.na(.)
    )
  )

# 使用 map_df 自动合并结果
med_results_all <- map_df(outcomes, function(outcome_var) {
  cat("Processing outcome:", outcome_var, "\n")
  
  result <- run_mediation_glmer_batch(
    data        = data_med,
    id_var      = "ID",
    exposures   = exposures,
    mediators   = mediators,
    outcomes    = outcome_var,
    covariates  = c(Covariates_all_lmer, "Fiber_soluble_mean"),
    exposure_high = 2,
    exposure_low  = 0,
    sims = 1000
  )
  
  return(result)
}, .id = "Outcome_ID")  # 自动添加 ID 列
####***有时间顺序的中介效应分析***####
####F024####
data_list <- list(Cov_F0_final, Tea_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], ASM_Height_long[ASM_Height_long$Followup == "F4",], Micro_selected_filter_clr[Micro_selected_filter_clr$Times == "F2",]) #
Micro_med_data1_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
Micro_names <- grep("^s__.*_clr$", colnames(Micro_med_data1_v1), value = TRUE)

# 合并握力等
Micro_med_data1_v2 <- left_join(Micro_med_data1_v1,SPPB_grip_long[SPPB_grip_long$Followup == "F4",],by = "ID")
Micro_med_data1_v2 <- left_join(Micro_med_data1_v2, Catechin, by = "ID")
#Micro_med_data1_v2 <- Micro_med_data1_v2[complete.cases(Micro_med_data1_v2),]


Micro_med_data1_v3 <- group_variables3(
  data = Micro_med_data1_v2,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(Micro_med_data1_v3)%in% factor_name)
for(i in idx ){
  Micro_med_data1_v3[[i]] <-  as.factor(Micro_med_data1_v3[[i]])
}
#***************ASMI根据性别Z分数，菌群Z分数
Micro_med_data1_v4 <- scale_columns_group(Micro_med_data1_v3, c("ASM", "ASMI","grip_max"), c("Sex_F0"))
Micro_med_data1_final_F024 <- scale_columns(Micro_med_data1_v4, Micro_names)
####中介效应分析####
Micro_med_results2_F024 <- run_microbiome_mediation_batch_glm(
  data = Micro_med_data1_final_F024,
  micro_sig = med_results_all_df_micro_sig,
  covariates =  c(Covariates_all,"Fiber_soluble_mean"),
  predictor_col = "Exposure",
  mediator_col = "Mediator",
  outcome_col = "Outcome",
  exposure_high = 2,
  exposure_low = 0,
  sims = 1000,
  output_file = "microbiome_mediation_results_time.csv"
)

Micro_med_results2_F024_sig <- Micro_med_results2_F024[Micro_med_results2_F024$ACME_P_value < 0.05,]
Micro_med_results2_F024 <- format_mediation_table(Micro_med_results2_F024, include_fdr = FALSE)
Micro_med_results2_F024 <- Micro_med_results2_F024 %>%
  mutate(
    
    Exposure = recode(
      Exposure,
      !!!predictor_rename
    ),
    
    Mediator = str_remove(Mediator, "s__") %>% str_remove("_clr_z") %>% str_replace_all("_", " "),
    
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    )
  )
####F023####
data_list <- list(Cov_F0_final, Tea_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], ASM_Height_long[ASM_Height_long$Followup == "F3",], Micro_selected_filter_clr[Micro_selected_filter_clr$Times == "F2",]) #
Micro_med_data1_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
Micro_names <- grep("^s__.*_clr$", colnames(Micro_med_data1_v1), value = TRUE)

# 合并握力等
Micro_med_data1_v2 <- left_join(Micro_med_data1_v1,SPPB_grip_long[SPPB_grip_long$Followup == "F3",],by = "ID")
Micro_med_data1_v2 <- left_join(Micro_med_data1_v2, Catechin, by = "ID")
#Micro_med_data1_v2 <- Micro_med_data1_v2[complete.cases(Micro_med_data1_v2),]


Micro_med_data1_v3 <- group_variables3(
  data = Micro_med_data1_v2,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(Micro_med_data1_v3)%in% factor_name)
for(i in idx ){
  Micro_med_data1_v3[[i]] <-  as.factor(Micro_med_data1_v3[[i]])
}
#***************ASMI根据性别Z分数，菌群Z分数
Micro_med_data1_v4 <- scale_columns_group(Micro_med_data1_v3, c("ASM", "ASMI","grip_max"), c("Sex_F0"))
Micro_med_data1_final_F023 <- scale_columns(Micro_med_data1_v4, Micro_names)
####中介效应分析####
Micro_med_results2_F023 <- run_microbiome_mediation_batch_glm(
  data = Micro_med_data1_final_F023,
  micro_sig = med_results_all_df_micro_sig,
  covariates =  c(Covariates_all,"Fiber_soluble_mean"),
  predictor_col = "Exposure",
  mediator_col = "Mediator",
  outcome_col = "Outcome",
  exposure_high = 2,
  exposure_low = 0,
  sims = 1000,
  output_file = "microbiome_mediation_results_time.csv"
)

Micro_med_results2_F023_sig <- Micro_med_results2_F023[Micro_med_results2_F023$ACME_P_value < 0.05,]
Micro_med_results2_F023 <- format_mediation_table(Micro_med_results2_F023, include_fdr = FALSE)
Micro_med_results2_F023 <- Micro_med_results2_F023 %>%
  mutate(
    
    Exposure = recode(
      Exposure,
      !!!predictor_rename
    ),
    
    Mediator = str_remove(Mediator, "s__") %>% str_remove("_clr_z") %>% str_replace_all("_", " "),
    
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    )
  )
####*****************************data8-饮茶&功能&骨骼肌*****************************####
####儿茶素数据集处理合并####
Micro_ASMI <- merge(Micro_func_filter_clr, ASM_Height_long, by.x = c("ID","Times"), by.y = c("ID","Followup"))
table(Micro_ASMI$Times)

data_list <- list(Cov_F0_final, Tea_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], Micro_ASMI) #
data8_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
Func_names <- grep("_clr$", colnames(data8_v1), value = TRUE)

# 合并握力等
data8_v2 <- left_join(data8_v1,SPPB_grip_long,by = c("ID" = "ID","Times" = "Followup"))
data8_v2 <- left_join(data8_v2, Catechin, by = c("ID" = "ID"))

data8_v3 <- group_variables3(
  data = data8_v2,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(data8_v3)%in% factor_name)
for(i in idx ){
  data8_v3[[i]] <-  as.factor(data8_v3[[i]])
}
#***************ASMI根据性别Z分数，菌群Z分数
data8_v4 <- scale_columns_group(data8_v3, c("ASM", "ASMI","grip_max"), c("Sex_F0"))
data8_final <- scale_columns(data8_v4, Func_names)

table(data8_final$Times)
a <- unique(data8_final$ID)
a <- unique(data8_final[data8_final$tea_freq_group %in% c(0,2),]$ID)
a <- unique(data8_final$ID)
####lmer-简单调整####
Lmer_results8_all <- process_lmer(c("tea_freq_group"), 
                                  paste0(Func_names,"_z"), 
                                  data8_final, 
                                  c("Age","Sex_F0"), # 
                                  2)

Lmer_results8_all <- Lmer_results8_all[Lmer_results8_all$Level == 2,]
Lmer_results8_all$P_FDR <- p.adjust(Lmer_results8_all$P_value, method = "fdr")

Lmer_results8_all_sig <- Lmer_results8_all[Lmer_results8_all$P_value <0.05 & Lmer_results8_all$P_FDR < 0.10,] #& Lmer_results7_all_level$P_FDR < 0.05
#*************正文写作
a <- Lmer_results8_all_sig[Lmer_results8_all_sig$Estimate >0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)

a <- Lmer_results8_all_sig[Lmer_results8_all_sig$Estimate <0,]
round(a$Estimate,3)
round(a$CI_low,3)
round(a$CI_high,3)
#*******************数据清洗导出
# 数据清洗和转换（不使用函数，直接在 mutate 中判断）
Lmer_results8_all_output <- Lmer_results8_all %>%
  # 只保留需要的列
  dplyr::select(Outcome, Predictor, Level, Estimate, CI_low, CI_high, P_value, P_FDR) %>%
  # 清洗 Outcome 列
  mutate(
    Outcome = str_remove_all(Outcome, "_clr_z"),

    # 根据原始 Predictor 名称转换 Level 标签
    Level_Label = case_when(
      # tea_freq_group 的处理
      Predictor == "tea_freq_group" ~ "≥ 7 times/week vs. Non-drinker",
      # epicatechin 和 EGC 的处理
      Predictor %in% c("serum_I_epicatechin_F0_T", "serum_I_EGC_F0_T") ~ "High vs. Undetectable",
      # catechin, EGCG, ECG 的处理
      Predictor %in% c("serum_I_catechin_F0_T", "serum_I_EGCG_F0_T", "serum_I_ECG_F0_T") ~ "Tertile 3 vs. Tertile 1",
      TRUE ~ as.character(Level)
    ),
    
    # 重命名 Predictor（在 Level 转换之后进行）
    Predictor = recode(Predictor, !!!predictor_rename),
    
    # 构建 β (95% CI)
    `β (95% CI)` = sprintf("%.3f (%.3f, %.3f)", Estimate, CI_low, CI_high),
    
    # 格式化 P 值
    P_value_formatted = case_when(
      P_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", P_value)
    ),
    P_FDR_formatted = case_when(
      P_FDR < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", P_FDR)
    )
  ) %>%
  # 重新排列列顺序
  dplyr::select(Outcome, Predictor, Level = Level_Label, `β (95% CI)`, 
                p = P_value_formatted, 
                p_FDR = P_FDR_formatted) %>%
  # 按 Predictor 和 Outcome 排序
  arrange(Predictor, Outcome)
####可视化####
library(dplyr)
library(ggplot2)
library(stringr)
library(ggnewscale)

## 数据整理

plot_df <- Lmer_results8_all_sig %>%
  
  mutate(
    
    ## pathway名称
    Pathway = Outcome,
    
    Pathway = str_remove(
      Pathway,
      "_clr_z$"
    ),
    
    ## 自动换行
    Pathway = str_wrap(
      Pathway,
      width = 25
    ),
    
    ## 方向
    Direction = ifelse(
      Estimate > 0,
      "Positive",
      "Negative"
    )
  ) %>%
  
  arrange(Estimate) %>%
  
  mutate(
    id = 1:n()
  )

## 环形文字角度

number_of_bar <- nrow(plot_df)

angle <- 90 - 360 * (plot_df$id - 0.5) / number_of_bar

plot_df$hjust <- ifelse(
  angle < -90,
  1,
  0
)

plot_df$angle <- ifelse(
  angle < -90,
  angle + 180,
  angle
)

## 圆环上下限

inner_radius <- -0.22

outer_limit <- max(abs(plot_df$Estimate)) + 0.12

## 作图

p <- ggplot(
  plot_df,
  aes(
    x = factor(id),
    y = abs(Estimate)
  )
) +
  
  ## 浅灰背景环
  
  geom_col(
    aes(y = outer_limit),
    fill = "#F3F3F3",
    width = 0.92,
    alpha = 1
  ) +
  
  ## 主圆环
  
  geom_col(
    aes(
      fill = Estimate
    ),
    width = 0.92,
    alpha = 0.96,
    color = "white",
    linewidth = 0.5
  ) +
  
  ## 外圈点缀
  
  geom_point(
    aes(
      y = abs(Estimate) + 0.015,
      color = Estimate
    ),
    size = 2.8
  ) +
  
  ## pathway文字
  
  geom_text(
    aes(
      y = abs(Estimate) + 0.05,
      label = Pathway,
      hjust = hjust
    ),
    angle = plot_df$angle,
    size = 5.5,
    lineheight = 0.92,
    color = "black",
    family = "sans"
  ) +
  
  ## 中心文字
  
  annotate(
    "text",
    x = 0,
    y = inner_radius + 0.03,
    label = "Tea\nconsumption",
    fontface = "bold",
    size = 8,
    lineheight = 0.9,
    color = "black"
  ) +
  
  ## 副标题
  annotate(
    "text",
    x = 0,
    y = inner_radius - 0.02,
    label = "Gut microbial pathways",
    size = 4.2,
    color = "grey40"
  ) +
  
  ## 极坐标
  
  coord_polar(start = 0) +
  
  ## 中间空洞
  ylim(inner_radius, outer_limit + 0.12) +
  
  ## 高级配色
  
  scale_fill_gradient2(
    low = "#3D7FB6",
    mid = "#F7F7F7",
    high = "#D95F5F",
    midpoint = 0,
    name = "Estimate"
  ) +
  
  scale_color_gradient2(
    low = "#3D7FB6",
    mid = "#F7F7F7",
    high = "#D95F5F",
    midpoint = 0
  ) +
  
  ## 去掉重复legend
  guides(
    color = "none"
  ) +
  
  ## 极简主题
  
  theme_void(base_size = 14) +
  
  theme(
    
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 14
    ),
    
    legend.text = element_text(
      size = 12
    ),
    
    plot.margin = ggplot2::margin(
      15,
      15,
      15,
      15
    )
  )

# 4. 导出完美比例的长方形超清图（因为文字和表格是横向展开的，宽度必须大一点）
ggsave("Forest_Plot_Ultra_Clear.png", plot = p, width = 15, height = 15, dpi = 500)
####lmer-全调整####
Lmer_results8_all2 <- process_lmer(c("tea_freq_group", "serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T"),
                                   unique(Lmer_results8_all_sig$Outcome),
                                   data8_final,
                                   c(setdiff(Covariates_all_lmer,c("Protein_mean")),c("Fiber_soluble_mean")),
                                   2)
Lmer_results8_all2 <- Lmer_results8_all2[Lmer_results8_all2$Level == 2,]
Lmer_results8_all2 <- Lmer_results8_all2 %>%
  group_by(Predictor) %>%
  mutate(
    P_FDR = p.adjust(P_value, method = "fdr")
  ) %>%
  ungroup()

Lmer_results8_all2_sig <- Lmer_results8_all2[Lmer_results8_all2$P_value < 0.05 & Lmer_results8_all2$P_FDR < 0.10,]

#*************正文写作
a <- Lmer_results8_all2_sig[Lmer_results8_all2_sig$Estimate >0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)

a <- Lmer_results8_all2_sig[Lmer_results8_all2_sig$Estimate <0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)
####可视化####
# 1. 数据深度预处理（添加动态星号逻辑）
df_forest <- Lmer_results8_all2  %>%
  mutate(Predictor = ifelse(Predictor %in% names(predictor_rename),
                            predictor_rename[Predictor],
                            Predictor),
         Predictor = factor(Predictor, levels = predictor_rename)
  )

df_forest$Outcome_clean <- str_replace(df_forest$Outcome, "_clr_z$", "")
df_forest$Outcome_clean <- str_wrap(df_forest$Outcome_clean, width = 45)

# 核心阶梯状排序
df_forest$Outcome_clean <- fct_reorder(df_forest$Outcome_clean, df_forest$Estimate)

# 细化显著性分类（用于图形颜色/形状映射）
df_forest$Is_Significant <- ifelse(df_forest$P_value <0.05 & df_forest$P_FDR < 0.10, "Significant", "Non-significant")

# 【核心加分项】定义严格的学术星号，并计算安全悬浮的 X 轴坐标
df_forest <- df_forest %>%
  mutate(
    # 按标准划分星号：P < 0.05 (*), P < 0.01 (**), P < 0.001 (***)
    Significance_Star = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE            ~ "" # 不显著则留空
    ),
    # 动态避让：让星号落在置信区间右端点（CI_high）再往右侧延伸全图总跨度的 2% 处
    Star_X = CI_high + (max(CI_high) - min(CI_low)) * 0.02
  )

# 2. 斑马纹背景数据精准构建（保障分面平铺）
levels_y <- levels(df_forest$Outcome_clean)
bg_intervals <- data.frame(
  Outcome_clean = factor(levels_y, levels = levels_y),
  y_index = seq_along(levels_y)
) %>% 
  filter(y_index %% 2 == 1)

# 3. 绘制带星号的至臻森林图
p_refined <- ggplot(df_forest, aes(x = Estimate, y = Outcome_clean)) +
  
  # A. 【底层】全图层斑马纹阴影
  geom_rect(data = bg_intervals,
            aes(ymin = y_index - 0.5, ymax = y_index + 0.5, xmin = -Inf, xmax = Inf),
            fill = "#F1F5F9", color = NA, alpha = 0.6, inherit.aes = FALSE) +
  
  # B. 无效应参考虚线
  geom_vline(xintercept = 0, color = "#94A3B8", linetype = "twodash", linewidth = 0.7) +
  
  # C. 绘制置信区间（毛毛虫线）
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, color = Is_Significant), 
                 height = 0.22, linewidth = 0.95, alpha = 0.9) +
  
  # D. 绘制 Estimate 核心点（虚实结合高级质感）
  geom_point(aes(color = Is_Significant, fill = Is_Significant, size = Is_Significant), 
             shape = 21, stroke = 1.2) +
  
  # E. 【核心精修图层】动态标注显著性星号
  # vjust = 0.38 让星号在垂直方向上完美对齐线条中轴，不会飘高或下沉
  # geom_text(aes(x = Star_X, label = Significance_Star),
  #           color = "#C94A29", size = 4.5, fontface = "bold", hjust = 0, vjust = 0.38) +
  
  # F. 分面系统调整（左右并排对比）
  facet_grid(. ~ Predictor, scales = "free_x") +
  
  # G. 顶级学术色彩与大小映射
  scale_color_manual(values = c("Significant" = "#C94A29", "Non-significant" = "#64748B")) +
  scale_fill_manual(values = c("Significant" = "#C94A29", "Non-significant" = "#FFFFFF")) + 
  scale_size_manual(values = c("Significant" = 3.2, "Non-significant" = 2.0)) +
  
  # H. 适当微调右侧边界扩充，确保最长那根线的三个星号不会被画布边缘切掉
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08))) + 
  
  # 4. 极致美化主题（学术开放式版面）
  theme_minimal(base_size = 13, base_family = "sans") +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    
    # 分面条带高级雕刻
    strip.background = element_rect(fill = "#F8FAFC", color = "#E2E8F0", linewidth = 1),
    strip.text = element_text(size = 14, face = "bold", color = "#1E293B", margin = ggplot2::margin(t=10, b=10)),
    panel.spacing = unit(0.5, "lines"), # 稍微拉开分面间距，图表更有呼吸感
    
    # 坐标轴精致美化
    axis.line.x = element_line(color = "#475569", linewidth = 0.6), 
    axis.ticks.x = element_line(color = "#475569", linewidth = 0.6),
    axis.text.x = element_text(size = 11, color = "#475569", margin = ggplot2::margin(t = 5)),
    
    # Y轴结局变量排版调整
    axis.text.y = element_text(size = 11.5, face = "italic", color = "#0F172A", lineheight = 1.05),
    axis.title.x = element_text(size = 13, face = "bold", color = "#0F172A", margin = ggplot2::margin(t = 12)),
    
    legend.position = "top",
    plot.margin = ggplot2::margin(20, 20, 20, 20)
  ) +
  labs(
    x = "β (95% CI)",
    y = NULL,
    color = NULL,
    fill = NULL,
    size = NULL
  )

# 5. 出版级 PDF 高清渲染输出
# cairo_pdf("Figure 4.pdf", width = 18, height = 12)
# print(p_refined)
# dev.off()

ggsave("Figure 4.png",
       plot = p_refined,
       width = 18,
       height = 12,
       dpi = 600)
####功能&骨骼肌####
Lmer_results8_all_ASM <- process_lmer(Lmer_results8_all2_sig$Outcome, 
                                      c("ASM_z","ASMI_z","grip_max_z"),
                                      data8_final, #%>% filter(tea_freq_group %in% c(0, 2)),
                                      c(Covariates_all_lmer,c("Fiber_soluble_mean","tea_freq_group")), #
                                      2)

Lmer_results8_all_ASM <- Lmer_results8_all_ASM %>%
  group_by(Outcome) %>%
  mutate(
    P_FDR = p.adjust(P_value, method = "fdr")
  ) %>%
  ungroup()

Lmer_results8_all_ASM_sig <- Lmer_results8_all_ASM[Lmer_results8_all_ASM$P_value < 0.05 & Lmer_results8_all_ASM$P_FDR <0.10,]
#*************正文写作
a <- Lmer_results8_all_ASM_sig[Lmer_results8_all_ASM_sig$Estimate >0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)

a <- Lmer_results8_all_ASM_sig[Lmer_results8_all_ASM_sig$Estimate <0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)
####可视化####
library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)

plot_df <- Lmer_results8_all_ASM %>%
  
  mutate(
    
    ## pathway名称
    Pathway = Predictor,
    
    Pathway = str_remove(
      Pathway,
      "_clr_z$"
    ),
    
    Pathway = str_wrap(
      Pathway,
      width = 60
    ),
    
    ## outcome
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    ),
    
    ## significance stars
    sig = case_when(
      P_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  
  group_by(Outcome) %>%
  
  arrange(Estimate) %>%
  
  mutate(
    Pathway = factor(
      Pathway,
      levels = unique(Pathway)
    )
  ) %>%
  
  ungroup()


p <- plot_forest_facet(
  data = plot_df,
  estimate_col = "Estimate",
  ci_low_col = "CI_low",
  ci_high_col = "CI_high",
  pathway_col = "Pathway", #因变量
  outcome_col = "Outcome", #分面
  p_col = "P_value",
  sig_rule = "p_fdr"
)

ggsave("Func_ASM.png", plot = p, width = 13, height = 11, dpi = 500)
# pdf("Func_ASM.pdf", width = 12, height = 11)
# print(p)
# dev.off()
####模型效能的提高####
Lmer_model_compare_func <- compare_lmer_predictors_by_pair(
  pair_df = Lmer_results8_all_ASM_sig,
  data = data8_final,
  covariates = c(Covariates_all_lmer, "Fiber_soluble_mean", "tea_freq_group"),
  predictor_col = "Predictor",
  outcome_col = "Outcome",
  predictor_type_col = "Variable_Type",
  id_var = "ID",
  predictor_rename = NULL,
  outcome_rename = c(
    "ASM_z" = "ASM",
    "ASMI_z" = "ASMI",
    "grip_max_z" = "Handgrip strength"
  ),
  output_file_raw = "Func_ASM_model_performance_raw.csv",
  output_file_format = "Func_ASM_model_performance_format.csv"
)

Lmer_model_compare_func_output <- Lmer_model_compare_func$format[,c("Outcome", "Predictor","Delta_Marginal_R2_percent","LRT_P_value")]

Lmer_model_compare_func_output <- Lmer_model_compare_func_output %>%
  dplyr::mutate(
    Predictor = str_remove(Predictor, "_clr_z") 
    
  )

####中介效应分析-功能####
# Micro_sig <- Lmer_results7_all2_sig[Lmer_results7_all2_sig$Outcome %in% Lmer_results7_all4_sig$Predictor, ]
# Micro_ASM_sig <- Lmer_results7_all4_sig
# table(Micro_sig$Predictor)

Func_sig_pos <- Lmer_results8_all2_sig[Lmer_results8_all2_sig$Estimate >0,] %>%
  
  inner_join(
    
    Lmer_results8_all_ASM_sig[Lmer_results8_all_ASM_sig$Estimate >0,] %>%
      
      dplyr::select(
        Predictor,
        Outcome
      ) %>%
      
      distinct(),
    
    by = c(
      "Outcome" = "Predictor"
    )
  ) 

Func_sig_neg <- Lmer_results8_all2_sig[Lmer_results8_all2_sig$Estimate <0,] %>%
  
  inner_join(
    
    Lmer_results8_all_ASM_sig[Lmer_results8_all_ASM_sig$Estimate <0,] %>%
      
      dplyr::select(
        Predictor,
        Outcome
      ) %>%
      
      distinct(),
    
    by = c(
      "Outcome" = "Predictor"
    )
  )

Func_sig <- rbind(Func_sig_pos, Func_sig_neg)

colnames(Func_sig) <- recode(colnames(Func_sig),
                              "Outcome.y" = "Outcome_muscle",
                              "Outcome" = "Mediator")


med_results_all_df_Func <- 
  
  run_microbiome_mediation_batch(
    
    data = data8_final,
    
    micro_sig = Func_sig,
    
    covariates = c(
      Covariates_all_lmer,
      "Fiber_soluble_mean"
    ),
    
    mediator_col = "Mediator",
    
    outcome_col = "Outcome_muscle",
    
    predictor_col = "Predictor",
    
    exposure_high = 2,
    exposure_low = 0,
    
    sims = 1000,
    
    output_file =
      "mediation_results_all_Func1.csv"
  )

med_results_all_df_Func <- med_results_all_df_Func %>%
  group_by(Exposure) %>%
  mutate(
    ACME_FDR = p.adjust(ACME_p, method = "fdr")
  ) %>%
  ungroup()
med_results_all_df_Func_sig <- med_results_all_df_Func[med_results_all_df_Func$ACME_p < 0.05 & med_results_all_df_Func$ACME_FDR < 0.10,]
writexl::write_xlsx(med_results_all_df_Func, "med_results_all_df_Func.xlsx")
####***有时间顺序的中介效应***####
####F024####
data_list <- list(Cov_F0_final, Tea_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], ASM_Height_long[ASM_Height_long$Followup == "F4",], Micro_func_filter_clr[Micro_func_filter_clr$Times == "F2",]) #
Func_med_data1_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
Func_names <- grep("_clr$", colnames(Func_med_data1_v1), value = TRUE)

# 合并握力等
Func_med_data1_v2 <- left_join(Func_med_data1_v1,SPPB_grip_long[SPPB_grip_long$Followup == "F4",],by = "ID")
Func_med_data1_v2 <- left_join(Func_med_data1_v2, Catechin, by = "ID")

Func_med_data1_v3 <- group_variables3(
  data = Func_med_data1_v2,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(Func_med_data1_v3)%in% factor_name)
for(i in idx ){
  Func_med_data1_v3[[i]] <-  as.factor(Func_med_data1_v3[[i]])
}
#***************ASMI根据性别Z分数，菌群Z分数
Func_med_data1_v4 <- scale_columns_group(Func_med_data1_v3, c("ASM", "ASMI","grip_max"), c("Sex_F0"))
Func_med_data1_final_F024 <- scale_columns(Func_med_data1_v4, Func_names)
####中介效应分析####
Fun_med_results2_F024 <- run_microbiome_mediation_batch_glm(
  data = Func_med_data1_final_F024,
  micro_sig = med_results_all_df_Func_sig,
  covariates =  c(Covariates_all_lmer,"Fiber_soluble_mean"),
  predictor_col = "Exposure",
  mediator_col = "Mediator",
  outcome_col = "Outcome",
  exposure_high = 2,
  exposure_low = 0,
  sims = 1000,
  output_file = "Func_mediation_results.csv"
)

Fun_med_results2_F024_sig <- Fun_med_results2_F024[Fun_med_results2_F024$ACME_P_value < 0.05,]
Fun_med_results2_F024 <- format_mediation_table(Fun_med_results2_F024, include_fdr = FALSE)
Fun_med_results2_F024 <- Fun_med_results2_F024 %>%
  mutate(
    
    Exposure = recode(
      Exposure,
      !!!predictor_rename
    ),
    
    Mediator = str_remove(Mediator, "s__") %>% str_remove("_clr_z") %>% str_replace_all("_", " "),
    
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    )
    
  )
####F023####
data_list <- list(Cov_F0_final, Tea_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], ASM_Height_long[ASM_Height_long$Followup == "F3",], Micro_func_filter_clr[Micro_func_filter_clr$Times == "F2",]) #
Func_med_data1_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))
Func_names <- grep("_clr$", colnames(Func_med_data1_v1), value = TRUE)

# 合并握力等
Func_med_data1_v2 <- left_join(Func_med_data1_v1,SPPB_grip_long[SPPB_grip_long$Followup == "F3",],by = "ID")
Func_med_data1_v2 <- left_join(Func_med_data1_v2, Catechin, by = "ID")

Func_med_data1_v3 <- group_variables3(
  data = Func_med_data1_v2,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(Func_med_data1_v3)%in% factor_name)
for(i in idx ){
  Func_med_data1_v3[[i]] <-  as.factor(Func_med_data1_v3[[i]])
}
#***************ASMI根据性别Z分数，菌群Z分数
Func_med_data1_v4 <- scale_columns_group(Func_med_data1_v3, c("ASM", "ASMI","grip_max"), c("Sex_F0"))
Func_med_data1_final_F023 <- scale_columns(Func_med_data1_v4, Func_names)
####中介效应分析####
Fun_med_results2_F023 <- run_microbiome_mediation_batch_glm(
  data = Func_med_data1_final_F023,
  micro_sig = med_results_all_df_Func_sig,
  covariates =  c(Covariates_all_lmer,"Fiber_soluble_mean"),
  predictor_col = "Exposure",
  mediator_col = "Mediator",
  outcome_col = "Outcome",
  exposure_high = 2,
  exposure_low = 0,
  sims = 1000,
  output_file = "Func_mediation_results.csv"
)

Fun_med_results2_F023_sig <- Fun_med_results2_F023[Fun_med_results2_F023$ACME_P_value < 0.05,]
Fun_med_results2_F023 <- format_mediation_table(Fun_med_results2_F023, include_fdr = FALSE)
Fun_med_results2_F023 <- Fun_med_results2_F023 %>%
  mutate(
    
    Exposure = recode(
      Exposure,
      !!!predictor_rename
    ),
    
    Mediator = str_remove(Mediator, "s__") %>% str_remove("_clr_z") %>% str_replace_all("_", " "),
    
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    )
    
  )
####验证非重复测量中介效应分析####
a <- run_mediation_lm_once(Func_med_data1_final, 
                           exposure = "tea_freq_group",
                           mediator = "ARGSYN-PWY: L-arginine biosynthesis I (via L-ornithine)_clr_z",
                           outcome = "ASMI_z",
                           covariates = c(
                             Covariates_all_lmer,
                             "Fiber_soluble_mean"
                           )
)
####菌群对功能通路的贡献####
# 2. 提取目标功能通路
# 来源：Lmer_results8_all_ASM_sig$Predictor

target_pathways <- unique(
  gsub("_clr_z$", "", Lmer_results8_all_ASM_sig$Predictor)
)

target_pathway_ids <- sub(":.*$", "", target_pathways)

target_pathway_df <- tibble::tibble(
  target_pathway = target_pathways,
  target_pathway_id = target_pathway_ids
)

cat("目标功能通路数量：", length(target_pathways), "\n")
print(target_pathway_df)

target_taxa_clean <- unique(gsub("_clr_z$", "", Lmer_results7_all_ASM_sig$Predictor))

# 3. 整理 Micro_func_stra 的列名信息

all_cols <- colnames(Micro_func_stra)

meta_cols <- intersect(
  c("...1", "id", "ID", "sample_id", "SampleID", "batch", "tag"),
  all_cols
)

col_info <- tibble::tibble(
  colname = all_cols
) %>%
  mutate(
    # 如果是 pathway|species，只保留 pathway 部分
    pathway_name = str_replace(colname, "\\|.*$", ""),
    
    # 提取 pathway ID，例如 ARGSYN-PWY
    pathway_id = str_replace(pathway_name, ":.*$", ""),
    
    # 是否是带物种的分层通路
    is_stratified = str_detect(colname, "\\|"),
    
    # 是否是 UNINTEGRATED
    is_unintegrated = pathway_id == "UNINTEGRATED",
    
    # 是否为基本信息列
    is_meta = colname %in% meta_cols,
    
    # 是否匹配目标功能通路
    match_by_full_name = pathway_name %in% target_pathways,
    match_by_pathway_id = pathway_id %in% target_pathway_ids,
    is_target_pathway = match_by_full_name | match_by_pathway_id
  )

a <- col_info[col_info$pathway_name == "POLYISOPRENSYN-PWY: polyisoprenoid biosynthesis (E. coli)",]
a <- col_info[col_info$pathway_name == "PWY-7400: L-arginine biosynthesis IV (archaebacteria)",]

# 4. 筛选目标功能通路对应的列

target_col_info <- col_info %>%
  filter(
    is_target_pathway,
    !is_meta,
    !is_unintegrated
  ) %>%
  arrange(
    pathway_id,
    desc(is_stratified),
    colname
  )


# 5. 检查每个目标通路匹配到了多少列

pathway_match_check <- target_col_info %>%
  group_by(pathway_id, pathway_name) %>%
  summarise(
    n_total_cols = n(),
    n_stratified_cols = sum(is_stratified),
    n_unstratified_cols = sum(!is_stratified),
    .groups = "drop"
  ) %>%
  arrange(pathway_id)


# 6. 生成只包含目标功能通路的宽表

target_pathway_cols <- target_col_info$colname

Micro_func_stra_target_pathways <- Micro_func_stra %>%
  dplyr::select(
    all_of(meta_cols),
    all_of(target_pathway_cols)
  )


# 7. 单独生成 species-stratified 目标通路表
#    后续计算“哪些菌贡献这些通路”主要用这个表

target_stratified_col_info <- target_col_info %>%
  filter(is_stratified)

target_stratified_cols <- target_stratified_col_info$colname

Micro_func_stra_target_pathways_stratified <- Micro_func_stra %>%
  dplyr::select(
    all_of(meta_cols),
    all_of(target_stratified_cols)
  )


# 8. 如果没有匹配到，输出诊断信息

if (nrow(target_col_info) == 0) {
  
  cat("\n没有匹配到目标功能通路。\n")
  cat("请检查目标通路 ID 和 Micro_func_stra 中的 pathway ID 是否一致。\n")
  
  cat("\n目标通路 ID：\n")
  print(target_pathway_ids)
  
  cat("\nMicro_func_stra 中前50个 pathway ID：\n")
  print(
    col_info %>%
      filter(!is_meta) %>%
      distinct(pathway_id) %>%
      head(50)
  )
}

# Species-stratified pathway contribution analysis
# 计算目标菌对每一条目标功能通路的贡献

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(forcats)

# 1. 输入目标菌

target_taxa_clean <- c(
  "s__Erysipelatoclostridium_ramosum",
  "s__Actinomyces_SGB17132",
  "s__Actinomyces_sp_ph3",
  "s__Tyzzerella_nexilis",
  "s__Gemmiger_formicilis",
  "s__Phocaeicola_massiliensis",
  "s__Bacteroides_cellulosilyticus"
)

# 如果你还想把审稿人点名的物种也一起检查，可以取消下面注释
# target_taxa_clean <- unique(c(
#   target_taxa_clean,
#   "s__Bacteroides_fragilis",
#   "s__Faecalibacterium_prausnitzii",
#   "s__Coprococcus_comes"
# ))


# 2. 物种名标准化函数
#    目的是让：
#    s__Tyzzerella_nexilis
#    和
#    g__Tyzzerella.s__Tyzzerella_nexilis
#    能够匹配上

make_taxon_key <- function(x) {
  
  x <- as.character(x)
  
  # 如果是 pathway|taxon 格式，只保留 taxon
  x <- sub("^.*\\|", "", x)
  
  # 如果是 g__Genus.s__Species，保留 Species
  x <- sub("^.*\\.s__", "", x)
  
  # 去掉 s__
  x <- gsub("^s__", "", x)
  x <- gsub("^g__", "", x)
  
  # 统一格式
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x <- tolower(x)
  
  return(x)
}

make_taxon_show <- function(x) {
  
  x <- as.character(x)
  x <- sub("^.*\\|", "", x)
  x <- sub("^.*\\.s__", "", x)
  x <- gsub("^s__", "", x)
  x <- gsub("^g__", "", x)
  x <- gsub("_", " ", x)
  
  return(x)
}

target_taxa_df <- tibble::tibble(
  Target_Taxon_raw = target_taxa_clean,
  Taxon_key = make_taxon_key(target_taxa_clean),
  Target_Taxon_show = make_taxon_show(target_taxa_clean)
) %>%
  distinct(Taxon_key, .keep_all = TRUE)

# 3. 解析目标功能通路 stratified 表
#    这里使用你已经筛选好的：
#    Micro_func_stra_target_pathways_stratified

dat_stra <- Micro_func_stra_target_pathways_stratified

all_cols <- colnames(dat_stra)

meta_cols <- all_cols[!grepl("\\|", all_cols)]
func_cols <- all_cols[grepl("\\|", all_cols)]

if (length(func_cols) == 0) {
  stop("没有检测到包含 | 的 species-stratified pathway 列。")
}

col_info <- tibble::tibble(
  Colname = func_cols,
  Pathway = sub("\\|.*$", "", func_cols),
  Taxon_raw = sub("^.*\\|", "", func_cols)
) %>%
  mutate(
    Pathway_ID = sub(":.*$", "", Pathway),
    Taxon_key = make_taxon_key(Taxon_raw),
    Taxon_show = make_taxon_show(Taxon_raw)
  )

# 4. 转换为数值矩阵

X_df <- dat_stra %>%
  dplyr::select(all_of(func_cols))

X_df <- as.data.frame(X_df, check.names = FALSE)

X_df[] <- lapply(X_df, function(x) {
  suppressWarnings(as.numeric(as.character(x)))
})

X_mat <- as.matrix(X_df)
colnames(X_mat) <- func_cols

if (any(X_mat < 0, na.rm = TRUE)) {
  stop("检测到负值。贡献比例不能用 CLR 或 z-score 后的数据计算，请使用原始 species-stratified pathway abundance 表。")
}


# 5. 计算所有物种对每条通路的 pooled contribution
#    pooled contribution =
#    该菌该通路在所有样本中的总丰度 /
#    该通路所有菌在所有样本中的总丰度 × 100

col_sum_abundance <- colSums(X_mat, na.rm = TRUE)
col_mean_abundance <- colMeans(X_mat, na.rm = TRUE)
col_prevalence <- colMeans(X_mat > 0, na.rm = TRUE) * 100
col_n_present <- colSums(X_mat > 0, na.rm = TRUE)

all_species_contribution <- col_info %>%
  mutate(
    sum_abundance = col_sum_abundance[Colname],
    mean_abundance = col_mean_abundance[Colname],
    prevalence_percent = col_prevalence[Colname],
    n_present = col_n_present[Colname],
    n_samples = nrow(X_mat)
  ) %>%
  group_by(Pathway, Pathway_ID, Taxon_key, Taxon_show) %>%
  summarise(
    sum_abundance = sum(sum_abundance, na.rm = TRUE),
    mean_abundance = sum(mean_abundance, na.rm = TRUE),
    prevalence_percent = max(prevalence_percent, na.rm = TRUE),
    n_present = max(n_present, na.rm = TRUE),
    n_samples = max(n_samples, na.rm = TRUE),
    matched_cols = paste(Colname, collapse = "; "),
    .groups = "drop"
  ) %>%
  group_by(Pathway, Pathway_ID) %>%
  mutate(
    total_pathway_abundance = sum(sum_abundance, na.rm = TRUE),
    pooled_contribution_percent = ifelse(
      total_pathway_abundance > 0,
      sum_abundance / total_pathway_abundance * 100,
      NA_real_
    ),
    contribution_rank = dense_rank(rank(-pooled_contribution_percent))
  ) %>%
  ungroup() %>%
  arrange(Pathway_ID, contribution_rank)


# 6. 计算目标菌对每条通路的贡献
#    同时计算 pooled contribution 和 sample-level contribution

pathways <- unique(col_info$Pathway)

target_contribution_list <- list()
counter <- 1

for (pw in pathways) {
  
  pw_cols <- col_info %>%
    filter(Pathway == pw) %>%
    pull(Colname)
  
  pw_id <- unique(col_info$Pathway_ID[col_info$Pathway == pw])
  
  # 每个样本中，这条通路所有物种的总丰度
  denom <- rowSums(X_mat[, pw_cols, drop = FALSE], na.rm = TRUE)
  
  for (i in seq_len(nrow(target_taxa_df))) {
    
    taxon_key_i <- target_taxa_df$Taxon_key[i]
    taxon_raw_i <- target_taxa_df$Target_Taxon_raw[i]
    taxon_show_i <- target_taxa_df$Target_Taxon_show[i]
    
    matched_cols_i <- col_info %>%
      filter(
        Pathway == pw,
        Taxon_key == taxon_key_i
      ) %>%
      pull(Colname)
    
    if (length(matched_cols_i) > 0) {
      numer <- rowSums(X_mat[, matched_cols_i, drop = FALSE], na.rm = TRUE)
    } else {
      numer <- rep(0, nrow(X_mat))
    }
    
    valid <- denom > 0
    
    sample_contribution <- rep(NA_real_, length(denom))
    sample_contribution[valid] <- numer[valid] / denom[valid] * 100
    
    pooled_contribution <- ifelse(
      sum(denom, na.rm = TRUE) > 0,
      sum(numer, na.rm = TRUE) / sum(denom, na.rm = TRUE) * 100,
      NA_real_
    )
    
    rank_i <- all_species_contribution %>%
      filter(
        Pathway == pw,
        Taxon_key == taxon_key_i
      ) %>%
      dplyr::select(contribution_rank) %>%
      pull()
    
    if (length(rank_i) == 0) {
      rank_i <- NA
    }
    
    target_contribution_list[[counter]] <- tibble::tibble(
      Pathway = pw,
      Pathway_ID = pw_id,
      Target_Taxon_raw = taxon_raw_i,
      Target_Taxon_show = taxon_show_i,
      Taxon_key = taxon_key_i,
      
      Detected_in_this_pathway = length(matched_cols_i) > 0,
      Matched_column_number = length(matched_cols_i),
      Matched_columns = ifelse(
        length(matched_cols_i) > 0,
        paste(matched_cols_i, collapse = "; "),
        NA_character_
      ),
      
      pooled_contribution_percent = pooled_contribution,
      mean_sample_contribution_percent = mean(sample_contribution, na.rm = TRUE),
      median_sample_contribution_percent = median(sample_contribution, na.rm = TRUE),
      prevalence_percent = mean(numer > 0, na.rm = TRUE) * 100,
      n_present = sum(numer > 0, na.rm = TRUE),
      n_samples = length(numer),
      contribution_rank = rank_i,
      
      contribution_level = case_when(
        is.na(pooled_contribution) ~ "Not detected",
        pooled_contribution >= 5 ~ "Major contributor",
        pooled_contribution >= 1 ~ "Moderate contributor",
        pooled_contribution > 0 ~ "Low contributor",
        TRUE ~ "Not detected"
      )
    )
    
    counter <- counter + 1
  }
}

target_taxa_contribution <- bind_rows(target_contribution_list) %>%
  arrange(Pathway_ID, desc(pooled_contribution_percent))

target_taxa_contribution_match <- target_taxa_contribution[target_taxa_contribution$Matched_column_number == 1,]
colnames(target_taxa_contribution_match)
target_taxa_contribution_match_output <- target_taxa_contribution_match[,c("Pathway","Target_Taxon_show","pooled_contribution_percent", "prevalence_percent")]


writexl::write_xlsx(target_taxa_contribution_match_output, "target_taxa_contribution_match_output.xlsx")
####*****************************data9-饮茶&蛋白质&骨骼肌*****************************####
Protein_Age_all <- merge(Protein_long_all_fill, ASM_Height_long, by.x = c("ID","followup"), by.y = c("ID","Followup"))
Protein_Age_all <- Protein_Age_all %>%
  mutate(
    across(
      where(is.numeric),
      ~ replace_na(., 0)
    )
  )

data_list <- list(Cov_F0_final, Tea_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], Protein_Age_all) 
data9_v1<- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))

# 合并握力等
data9_v2 <- left_join(data9_v1, SPPB_grip_long, by = c("ID" = "ID","followup" = "Followup"))
data9_v2 <- left_join(data9_v2, Catechin, by = c("ID" = "ID"))

Protein_names_var <- grep("_fill_log2$", colnames(data9_v2), value = TRUE) 

data9_v3 <- group_variables3(
  data = data9_v2,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(data9_v3)%in% factor_name)
for(i in idx ){
  data9_v3[[i]] <-  as.factor(data9_v3[[i]])
}
####蛋白质数据转换####
data9_v4 <- data9_v3 %>%
  mutate(
    Phase_group = ifelse(
      Phase %in% c("phase1", "phase2", "phase3"),
      "phase1_3",
      "phase4"
    )
  )
data9_v5 <- scale_columns_group(data9_v4, Protein_names_var, "Phase_group")
data9_final <- scale_columns_group(data9_v5, c("ASM", "ASMI","grip_max"), c("Sex_F0"))
table(data9_final$followup)
a <- unique(data9_final[data9_final$tea_freq_group %in% c(0,2),]$ID)
a <- unique(data9_final$ID)
####简单调整####
Lmer_results9_all <- process_lmer(c("tea_freq_group"), 
                                  paste0(Protein_names_var,"_z"),
                                  data9_final,
                                  c("Age","Sex_F0","Phase_group"), #
                                  2)
Lmer_results9_all <- Lmer_results9_all[Lmer_results9_all$Level == 2,]
Lmer_results9_all$P_FDR <- p.adjust(Lmer_results9_all$P_value, method = "fdr")
Lmer_results9_all_FDR_sig <- Lmer_results9_all[Lmer_results9_all$P_FDR <0.05,]

Lmer_results9_all$Outcome_name <- sub("_fill_log2_z", "", Lmer_results9_all$Outcome)
Lmer_results9_all$Outcome_name <- toupper(Lmer_results9_all$Outcome_name)
Lmer_results9_all <- merge(Lmer_results9_all, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Outcome_name", by.y = "Protein_names")
#*************正文写作
a <- Lmer_results9_all_FDR_sig[Lmer_results9_all_FDR_sig$Estimate >0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)

a <- Lmer_results9_all_FDR_sig[Lmer_results9_all_FDR_sig$Estimate <0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)
# Lmer_results9_all <- Lmer_results9_all %>%
#   group_by(Predictor) %>%
#   mutate(
#     P_FDR = p.adjust(P_value, method = "fdr")
#   ) %>%
#   ungroup()
# Lmer_results9_all$P_FDR <- p.adjust(Lmer_results9_all$P_value, method = "fdr")
# Lmer_results9_all_sig <- Lmer_results9_all[Lmer_results9_all$P_value < 0.05 & Lmer_results9_all$P_FDR < 0.05,]

Lmer_results9_all_output <- Lmer_results9_all %>%
  # 只保留需要的列
  dplyr::select(Outcome, hgnc_symbol, Predictor, Level, Estimate, CI_low, CI_high, P_value, P_FDR) %>%
  # 清洗 Outcome 列
  mutate(
    Outcome = str_remove_all(Outcome, "_fill_log2_z"),
    
    # 根据原始 Predictor 名称转换 Level 标签
    Level_Label = case_when(
      # tea_freq_group 的处理
      Predictor == "tea_freq_group" ~ "≥ 7 times/week vs. Non-drinker",
      # epicatechin 和 EGC 的处理
      Predictor %in% c("serum_I_epicatechin_F0_T", "serum_I_EGC_F0_T") ~ "High vs. Undetectable",
      # catechin, EGCG, ECG 的处理
      Predictor %in% c("serum_I_catechin_F0_T", "serum_I_EGCG_F0_T", "serum_I_ECG_F0_T") ~ "Tertile 3 vs. Tertile 1",
      TRUE ~ as.character(Level)
    ),
    
    # 重命名 Predictor（在 Level 转换之后进行）
    Predictor = recode(Predictor, !!!predictor_rename),
    
    # 构建 β (95% CI)
    `β (95% CI)` = sprintf("%.3f (%.3f, %.3f)", Estimate, CI_low, CI_high),
    
    # 格式化 P 值
    P_value_formatted = case_when(
      P_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", P_value)
    ),
    P_FDR_formatted = case_when(
      P_FDR < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", P_FDR)
    )
  ) %>%
  # 重新排列列顺序
  dplyr::select(Outcome, hgnc_symbol, Predictor, Level = Level_Label, `β (95% CI)`, 
                p = P_value_formatted, 
                p_FDR = P_FDR_formatted) %>%
  # 按 Predictor 和 Outcome 排序
  arrange(Predictor, Outcome)

colnames(Lmer_results9_all_output) <- recode(colnames(Lmer_results9_all_output),
                                             "hgnc_symbol" = "Outcome (HGNC Symbol)")
Lmer_results9_all_output_sig <- Lmer_results9_all_output[Lmer_results9_all_output$P_FDR <0.05,]
####可视化#####
Lmer_results9_all_adjusted <- Lmer_results9_all
# 为了区分VTN和CFI
Lmer_results9_all_adjusted[Lmer_results9_all_adjusted$hgnc_symbol == "VTN",]$Estimate <- 0.215
# 数据整理
df_plot <- Lmer_results9_all_adjusted %>%
  mutate(
    
    # y轴
    log_p = -log10(P_FDR),
    
    # 分组
    Group = case_when(
      P_FDR < 0.05 & Estimate > 0 ~ "FDR Sig (Up)",
      P_FDR < 0.05 & Estimate < 0 ~ "FDR Sig (Down)",
      TRUE ~ "Not Significant"
    )
    
  )

# 火山图
p_volcano <- ggplot(
  df_plot,
  aes(x = Estimate, y = log_p)
) +
  
  # FDR阈值线
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "#CBD5E1",
    linewidth = 0.6
  ) +
  
  # Estimate=0
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    color = "#94A3B8",
    linewidth = 0.6
  ) +
  
  # 所有点
  geom_point(
    aes(fill = Group, size = Group),
    shape = 21,
    color = "white",
    stroke = 0.4,
    alpha = 0.85
  ) +
  
  # 显著基因加黑圈
  geom_point(
    data = filter(df_plot, P_FDR < 0.05),
    shape = 21,
    color = "#0F172A",
    fill = NA,
    size = 4.8,
    stroke = 1.2
  ) +
  
  # 显著基因标签
  geom_text_repel(
    data = filter(df_plot, P_FDR < 0.05),
    aes(label = hgnc_symbol),
    
    size = 4,
    fontface = "bold",
    color = "#0F172A",
    
    box_padding = 0.6,
    point_padding = 0.3,
    
    segment.color = "#475569",
    segment.size = 0.5,
    
    min.segment.length = 0,
    max.overlaps = Inf,
    
    force = 5,
    force_pull = 1
  ) +
  
  # 配色
  scale_fill_manual(
    values = c(
      "FDR Sig (Up)" = "#E07A5F",
      "FDR Sig (Down)" = "#2A6F97",
      "Not Significant" = "#F1F5F9"
    )
  ) +
  
  # 点大小
  scale_size_manual(
    values = c(
      "FDR Sig (Up)" = 4.2,
      "FDR Sig (Down)" = 4.2,
      "Not Significant" = 2.5
    )
  ) +
  
  # 主题
  theme_bw() +
  
  theme(
    
    axis.text = element_text(
      size = 11,
      color = "#475569"
    ),
    
    axis.title = element_text(
      size = 13,
      face = "bold",
      color = "#0F172A"
    ),
    
    panel.grid.major = element_line(
      color = "#F8F9FA",
      linewidth = 0.5
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "top",
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 10,
      face = "bold",
      color = "#334155"
    ),
    
    plot.margin = ggplot2::margin(
      20,
      20,
      20,
      20
    )
  ) +
  
  labs(
    x = "β",
    y = expression(bold(-log[10](P[FDR]))),
    subtitle = NULL
  )

# 导出
ggsave(
  "Hgnc_Volcano_Plot.png",
  plot = p_volcano,
  width = 7.5,
  height = 6.5,
  dpi = 500
)

####lmer-全调整####
Lmer_results9_all2 <- process_lmer(c("tea_freq_group", "serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T"),
                                   unique(Lmer_results9_all_FDR_sig$Outcome),
                                   data9_final,
                                   c(setdiff(Covariates_all_lmer,c("Protein_mean")), "Phase_group"),
                                   2)
Lmer_results9_all2 <- Lmer_results9_all2[Lmer_results9_all2$Level == 2,]
Lmer_results9_all2 <- Lmer_results9_all2 %>%
  group_by(Predictor) %>%
  mutate(
    P_FDR = p.adjust(P_value, method = "fdr")
  ) %>%
  ungroup()
Lmer_results9_all2_sig <- Lmer_results9_all2[Lmer_results9_all2$P_FDR < 0.05,]

Lmer_results9_all2$Outcome_name <- sub("_fill_log2_z", "", Lmer_results9_all2$Outcome)
Lmer_results9_all2$Outcome_name <- toupper(Lmer_results9_all2$Outcome_name)
Lmer_results9_all2 <- merge(Lmer_results9_all2, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Outcome_name", by.y = "Protein_names")
a <- Lmer_results9_all2[Lmer_results9_all2$P_FDR < 0.05 & Lmer_results9_all2$Predictor == "tea_freq_group",]
writexl::write_xlsx(a, "a.xlsx")
#*************正文写作
a <- Lmer_results9_all2_sig[Lmer_results9_all2_sig$Estimate >0 & Lmer_results9_all2_sig$Predictor == "tea_freq_group",]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)

a <- Lmer_results9_all2_sig[Lmer_results9_all2_sig$Estimate <0 & Lmer_results9_all2_sig$Predictor == "tea_freq_group",]



a <- Lmer_results9_all2_sig[Lmer_results9_all2_sig$Estimate >0  & Lmer_results9_all2_sig$Predictor != "tea_freq_group",]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)
####可视化####
df_plot <- Lmer_results9_all2 %>%
  
  mutate(
    
    ## Predictor 重命名
    Predictor_clean = ifelse(
      Predictor %in% names(predictor_rename),
      predictor_rename[Predictor],
      Predictor
    ),
    
    ## Predictor 顺序
    Predictor_clean = factor(
      Predictor_clean,
      levels = predictor_rename
    ),
    
    ## 显著性方向
    Direction = case_when(
      P_value < 0.05 & P_FDR <0.05 & Estimate > 0 ~ "Positive Sig",
      P_value < 0.05 & P_FDR <0.05 & Estimate < 0 ~ "Negative Sig",
      TRUE ~ "Not Significant"
    )
  )

## Matrix Facet Forest Plot

p_matrix <- ggplot(
  df_plot,
  aes(
    x = Estimate,
    y = hgnc_symbol
  )
) +
  
  ## 背景
  geom_rect(
    aes(
      ymin = -Inf,
      ymax = Inf,
      xmin = -Inf,
      xmax = Inf
    ),
    fill = "#F8F9FA",
    color = NA,
    alpha = 0.55
  ) +
  
  ## 0线
  geom_vline(
    xintercept = 0,
    color = "#CBD5E1",
    linetype = "dashed",
    linewidth = 0.65
  ) +
  
  ## CI
  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high,
      color = Direction,
      linewidth = Direction
    ),
    height = 0.16,
    alpha = 0.9
  ) +
  
  ## 点（固定大小）
  geom_point(
    aes(
      fill = Direction
    ),
    shape = 21,
    size = 4.2,
    color = "white",
    stroke = 0.75
  ) +
  
  ## 分面
  facet_grid(
    . ~ Predictor_clean,
    scales = "free_x"
  ) +
  
  ## 配色
  scale_color_manual(
    values = c(
      "Positive Sig" = "#E07A5F",
      "Negative Sig" = "#2A6F97",
      "Not Significant" = "#94A3B8"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Positive Sig" = "#E07A5F",
      "Negative Sig" = "#2A6F97",
      "Not Significant" = "#E2E8F0"
    )
  ) +
  
  ## CI粗细
  scale_linewidth_manual(
    values = c(
      "Positive Sig" = 1.15,
      "Negative Sig" = 1.15,
      "Not Significant" = 0.55
    ),
    guide = "none"
  ) +
  
  ## 图例标题
  guides(
    fill = guide_legend(
      title = NULL
    ),
    color = "none"
  ) +
  
  ## 主题
  theme_minimal(base_family = "sans") +
  
  theme(
    
    axis.text.y = element_text(
      size = 12,
      face = "bold",
      color = "#0F172A"
    ),
    
    axis.text.x = element_text(
      size = 9,
      color = "#475569"
    ),
    
    axis.title.x = element_text(
      size = 12,
      face = "bold",
      color = "#0F172A",
      margin = ggplot2::margin(t = 12)
    ),
    
    axis.title.y = element_blank(),
    
    panel.grid = element_blank(),
    
    strip.background = element_rect(
      fill = "#F1F5F9",
      color = NA
    ),
    
    strip.text = element_text(
      size = 10,
      face = "bold",
      color = "#1E293B",
      margin = ggplot2::margin(
        t = 8,
        b = 8
      )
    ),
    
    panel.spacing = unit(
      1.25,
      "lines"
    ),
    
    legend.position = "top",
    
    legend.title = element_text(
      size = 10,
      face = "bold",
      color = "#334155"
    ),
    
    legend.text = element_text(
      size = 9,
      color = "#475569"
    ),
    
    plot.margin = ggplot2::margin(
      20,
      20,
      20,
      20
    )
  ) +
  
  labs(
    x = "β (95% CI)",
    fill = NULL
  )

## 保存
ggsave(
  "Figure 5.png",
  plot = p_matrix,
  width = 13,
  height = 5.5,
  dpi = 600
)
# pdf("Figure 5.pdf", width = 13, height = 5.5)
# print(p_matrix)
# dev.off()
####蛋白质&骨骼肌####
Lmer_results9_all_ASM <- process_lmer(unique(Lmer_results9_all2_sig$Outcome), 
                                      c("ASM_z","ASMI_z","grip_max_z"),
                                      data9_final, #%>% filter(tea_freq_group %in% c(0, 2)),
                                      c(Covariates_all_lmer,"tea_freq_group"), #
                                      2)
Lmer_results9_all_ASM <- Lmer_results9_all_ASM %>%
  group_by(Outcome) %>%
  mutate(
    P_FDR = p.adjust(P_value, method = "fdr")
  ) %>%
  ungroup()

Lmer_results9_all_ASM_sig <- Lmer_results9_all_ASM[Lmer_results9_all_ASM$P_value < 0.05 & Lmer_results9_all_ASM$P_FDR <0.10,]

Lmer_results9_all_ASM$Predictor_name <- sub("_fill_log2_z", "", Lmer_results9_all_ASM$Predictor)
Lmer_results9_all_ASM$Predictor_name <- toupper(Lmer_results9_all_ASM$Predictor_name)
Lmer_results9_all_ASM <- merge(Lmer_results9_all_ASM, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Predictor_name", by.y = "Protein_names")

#*************正文写作
a <- Lmer_results9_all_ASM_sig[Lmer_results9_all_ASM_sig$Estimate >0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)

a <- Lmer_results9_all_ASM_sig[Lmer_results9_all_ASM_sig$Estimate <0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)
####可视化####
plot_df <- Lmer_results9_all_ASM %>%
  
  mutate(
    ## outcome
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    ),
    
    ## significance stars
    sig = case_when(
      P_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) 

p <- plot_forest_facet(
  data = plot_df,
  estimate_col = "Estimate",
  ci_low_col = "CI_low",
  ci_high_col = "CI_high",
  pathway_col = "hgnc_symbol", #因变量
  outcome_col = "Outcome", #分面
  p_col = "P_value",
  y_text_face = "plain",
  sig_rule = "p_fdr"
  
)

ggsave("Pro_ASM.png", plot = p, width = 13, height = 6, dpi = 500)
# pdf("Func_ASM.pdf", width = 12, height = 11)
# print(p)
# dev.off()
####模型效能的提高####
Lmer_model_compare_pro <- compare_lmer_predictors_by_pair(
  pair_df = Lmer_results9_all_ASM_sig,
  data = data9_final,
  covariates = c(Covariates_all_lmer,"tea_freq_group"),
  predictor_col = "Predictor",
  outcome_col = "Outcome",
  predictor_type_col = "Variable_Type",
  id_var = "ID",
  predictor_rename = NULL,
  outcome_rename = c(
    "ASM_z" = "ASM",
    "ASMI_z" = "ASMI",
    "grip_max_z" = "Handgrip strength"
  ),
  output_file_raw = "Pro_ASM_model_performance_raw.csv",
  output_file_format = "Pro_ASM_model_performance_format.csv"
)

Lmer_model_compare_pro_output <- Lmer_model_compare_pro$format[,c("Outcome", "Predictor","Delta_Marginal_R2_percent","LRT_P_value")]

Lmer_model_compare_pro_output <- Lmer_model_compare_pro_output %>%
  dplyr::mutate(
    Predictor = str_remove(Predictor, "_fill_log2_z") 
  )

Lmer_model_compare_pro_output$Predictor_name <- toupper(Lmer_model_compare_pro_output$Predictor)
Lmer_model_compare_pro_output <- merge(Lmer_model_compare_pro_output, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Predictor_name", by.y = "Protein_names")

Lmer_model_compare_pro_output <- dplyr::select(Lmer_model_compare_pro_output, -c("Predictor_name","Predictor"))
colnames(Lmer_model_compare_pro_output) <- recode(colnames(Lmer_model_compare_pro_output),
                                                  "hgnc_symbol" = "Predictor")

####中介效应分析-蛋白质####
Pro_sig_pos <- Lmer_results9_all2_sig[Lmer_results9_all2_sig$Estimate >0,] %>%
  
  inner_join(
    
    Lmer_results9_all_ASM_sig[Lmer_results9_all_ASM_sig$Estimate >0,] %>%
      
      dplyr::select(
        Predictor,
        Outcome
      ) %>%
      
      distinct(),
    
    by = c(
      "Outcome" = "Predictor"
    )
  ) 

Pro_sig_neg <- Lmer_results9_all2_sig[Lmer_results9_all2_sig$Estimate <0,] %>%
  
  inner_join(
    
    Lmer_results9_all_ASM_sig[Lmer_results9_all_ASM_sig$Estimate <0,] %>%
      
      dplyr::select(
        Predictor,
        Outcome
      ) %>%
      
      distinct(),
    
    by = c(
      "Outcome" = "Predictor"
    )
  )

Pro_sig <- rbind(Pro_sig_pos, Pro_sig_neg)

colnames(Pro_sig) <- recode(colnames(Pro_sig),
                              "Outcome.y" = "Outcome_muscle",
                              "Outcome" = "Mediator")


med_results_all_df_Pro <- 
  
  run_microbiome_mediation_batch(
    
    data = data9_final,
    
    micro_sig = Pro_sig,
    
    covariates = c(Covariates_all_lmer, "Phase_group"),

    mediator_col = "Mediator",
    
    outcome_col = "Outcome_muscle",
    
    predictor_col = "Predictor",
    
    exposure_high = 2,
    exposure_low = 0,
    
    sims = 1000,
    
    output_file =
      "mediation_results_all_Pro1.csv"
  )

# med_results_all_df_Pro <- med_results_all_df_Pro %>%
#   group_by(Exposure) %>%
#   mutate(
#     ACME_P_FDR = p.adjust(ACME_p, method = "fdr")
#   ) %>%
#   ungroup()
med_results_all_df_Pro$ACME_FDR <- p.adjust(med_results_all_df_Pro$ACME_p, method = "fdr")
med_results_all_df_Pro_sig <- med_results_all_df_Pro[med_results_all_df_Pro$ACME_FDR < 0.05,]
writexl::write_xlsx(med_results_all_df_Pro, "med_results_all_df_Pro.xlsx")
####中介效应一起可视化####
mediation_results_all_micro_final <- read_excel("D:/OneDrive/Papers/9 饮茶频率&骨骼肌/Output/mediation_results_all_micro.xlsx", sheet=1)
plot_data_species <- mediation_results_all_micro_final %>%
  filter(
    ACME_FDR < 0.05,
    ACME > 0
  )
plot_data_species$hgnc_symbol = str_remove(plot_data_species$Mediator, "s__") %>% str_remove("_clr_z") %>% str_replace_all("_", " ")


mediation_results_all_Func <- read_excel("D:/OneDrive/Papers/9 饮茶频率&骨骼肌/Output/med_results_all_df_Func.xlsx",sheet = 1)
plot_data_func <- mediation_results_all_Func %>%
  filter(
    ACME_p < 0.05,
    ACME_FDR < 0.10,
    ACME > 0
  ) 
plot_data_func$hgnc_symbol = str_remove(plot_data_func$Mediator,"_clr_z")

mediation_results_all_pro <- read_excel("D:/OneDrive/Papers/9 饮茶频率&骨骼肌/Output/med_results_all_df_Pro.xlsx",sheet = 1)
plot_data_pro <- mediation_results_all_pro %>%
  filter(
    ACME_FDR < 0.05,
    ACME > 0
  )
plot_data_pro$Mediator = str_remove(plot_data_pro$Mediator,"_fill_log2_z")
plot_data_pro$Mediator_name <- toupper(plot_data_pro$Mediator)
plot_data_pro <- merge(plot_data_pro, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Mediator_name", by.y = "Protein_names")
plot_data_pro <- dplyr::select(plot_data_pro,-Mediator_name)


plot_data_mediation <- rbind(plot_data_species, plot_data_func, plot_data_pro)
plot_data_mediation$PropMediated <- plot_data_mediation$PropMediated*100

max(plot_data_mediation$PropMediated)
min(plot_data_mediation$PropMediated)

plot_data_mediation <- plot_data_mediation %>%
  mutate(
    
    Exposure = recode(
      Exposure,
      !!!predictor_rename
    ),
    
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    )
    
  )

plot_data_mediation <- plot_mediation_alluvial(
  data = plot_data_mediation,
  title = NULL,
  fill_low = "#E6D8F5",
  fill_high = "#6A51A3" ,
  mediator_col = "hgnc_symbol"
)


ggsave(
  "Figure 6.png",
  plot = plot_data_mediation,
  width = 11,
  height = 12,
  dpi = 600
)

#***********正文写作
round(min(plot_data_species$ACME),3)
round(max(plot_data_species$ACME),3)
round(min(plot_data_species$ACME_low),3)
round(max(plot_data_species$ACME_high),3)
round(min(plot_data_species$PropMediated*100),1)
round(max(plot_data_species$PropMediated*100),1)

round(min(plot_data_func$ACME),3)
round(max(plot_data_func$ACME),3)
round(min(plot_data_func$ACME_low),3)
round(max(plot_data_func$ACME_high),3)
round(min(plot_data_func$PropMediated*100),1)
round(max(plot_data_func$PropMediated*100),1)

round(min(plot_data_pro$ACME),3)
round(max(plot_data_pro$ACME),3)
round(min(plot_data_pro$ACME_low),3)
round(max(plot_data_pro$ACME_high),3)
round(min(plot_data_pro$PropMediated*100),1)
round(max(plot_data_pro$PropMediated*100),1)
####中介效应正文写作####
a <- mediation_results_all_pro[mediation_results_all_pro$ACME_p < 0.05,]
a$Mediator <- str_remove(a$Mediator,"_fill_log2_z")
a$Mediator <- toupper(a$Mediator)
a <- merge(a, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Mediator", by.y = "Protein_names")

a$ACME <- round(a$ACME, 3)
a$ACME_low <- round(as.numeric(a$ACME_low), 3)
a$ACME_high <- round(as.numeric(a$ACME_high), 3)
a$PropMediated <- round(a$PropMediated, 1)

a1 <- a[a$hgnc_symbol == "VTN",]
min(a1$ACME)
max(a1$ACME)
min(a1$ACME_low)
max(a1$ACME_high)
min(a1$PropMediated)
max(a1$PropMediated)

a1 <- a[a$hgnc_symbol == "CFI",]
min(a1$ACME)
max(a1$ACME)
min(a1$ACME_low)
max(a1$ACME_high)
min(a1$PropMediated)
max(a1$PropMediated)

a1 <- a[a$hgnc_symbol == "ITIH4",]
min(a1$ACME)
max(a1$ACME)
min(a1$ACME_low)
max(a1$ACME_high)
min(a1$PropMediated)
max(a1$PropMediated)

a1 <- a[a$hgnc_symbol == "APOF",]
min(a1$ACME)
max(a1$ACME)
min(a1$ACME_low)
max(a1$ACME_high)
min(a1$PropMediated)
max(a1$PropMediated)

a1 <- a[a$hgnc_symbol == "CNDP1",]
min(a1$ACME)
max(a1$ACME)
min(a1$ACME_low)
max(a1$ACME_high)
min(a1$PropMediated)
max(a1$PropMediated)
####***有时间顺序的中介效应***####
####F0-F2-F3####
Pro_names <- grep("_fill_log2$", colnames(Protein_long_all_fill), value = TRUE) 
data_list <- list(Cov_F0_final, Tea_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], ASM_Height_long[ASM_Height_long$Followup == "F3",], Protein_long_all_fill[Protein_long_all_fill$followup == "F2",c("ID","followup","Phase",Pro_names)]) #
Pro_med_data1_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))

# 合并握力等
Pro_med_data1_v2 <- left_join(Pro_med_data1_v1,SPPB_grip_long[SPPB_grip_long$Followup == "F3",],by = "ID")
Pro_med_data1_v2 <- left_join(Pro_med_data1_v2, Catechin, by = "ID")

Pro_med_data1_v3 <- group_variables3(
  data = Pro_med_data1_v2,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(Pro_med_data1_v3)%in% factor_name)
for(i in idx ){
  Pro_med_data1_v3[[i]] <-  as.factor(Pro_med_data1_v3[[i]])
}
#***************ASMI根据性别Z分数，菌群Z分数
Pro_med_data1_v4 <- scale_columns_group(Pro_med_data1_v3, c("ASM", "ASMI","grip_max"), c("Sex_F0"))
Pro_med_data1_v4 <- Pro_med_data1_v4 %>%
  mutate(
    Phase_group = ifelse(
      Phase %in% c("phase1", "phase2", "phase3"),
      "phase1_3",
      "phase4"
    )
  )
Pro_med_data1_final_F023 <- scale_columns_group(Pro_med_data1_v4, Pro_names, "Phase_group")
####中介效应####
Pro_med_results2_F023 <- run_microbiome_mediation_batch_glm(
  data = Pro_med_data1_final_F023,
  micro_sig = med_results_all_df_Pro_sig,
  covariates =  c(Covariates_all, "Phase_group"),
  predictor_col = "Exposure",
  mediator_col = "Mediator",
  outcome_col = "Outcome",
  exposure_high = 2,
  exposure_low = 0,
  sims = 1000,
  output_file = "Proteins_mediation_results_time.csv"
)

Pro_med_results2_F023_sig <- Pro_med_results2_F023[Pro_med_results2_F023$ACME_P_value < 0.05,]

Pro_med_results2_F023 <- format_mediation_table(Pro_med_results2_F023, include_fdr = FALSE)

Pro_med_results2_F023 <- Pro_med_results2_F023 %>%
  mutate(
    
    Exposure = recode(
      Exposure,
      !!!predictor_rename
    ),
    
    Mediator = str_remove(
      Mediator,
      "_fill_log2_z"
    ),
    
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    )
    
  )
Pro_med_results2_F023$Mediator <- toupper(Pro_med_results2_F023$Mediator)
Pro_med_results2_F023 <- merge(Pro_med_results2_F023, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Mediator", by.y = "Protein_names")
Pro_med_results2_F023 <- dplyr::select(Pro_med_results2_F023, -Mediator)
colnames(Pro_med_results2_F023) <- recode(colnames(Pro_med_results2_F023), "hgnc_symbol" = "Mediator")
####F0-F2-F4####
Pro_names <- grep("_fill_log2$", colnames(Protein_long_all_fill), value = TRUE) 
data_list <- list(Cov_F0_final, Tea_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], ASM_Height_long[ASM_Height_long$Followup == "F4",], Protein_long_all_fill[Protein_long_all_fill$followup == "F2",c("ID","followup","Phase",Pro_names)]) #
Pro_med_data1_v1 <- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))

# 合并握力等
Pro_med_data1_v2 <- left_join(Pro_med_data1_v1,SPPB_grip_long[SPPB_grip_long$Followup == "F4",],by = "ID")
Pro_med_data1_v2 <- left_join(Pro_med_data1_v2, Catechin, by = "ID")

Pro_med_data1_v3 <- group_variables3(
  data = Pro_med_data1_v2,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(Pro_med_data1_v3)%in% factor_name)
for(i in idx ){
  Pro_med_data1_v3[[i]] <-  as.factor(Pro_med_data1_v3[[i]])
}
#***************ASMI根据性别Z分数，菌群Z分数
Pro_med_data1_v4 <- scale_columns_group(Pro_med_data1_v3, c("ASM", "ASMI","grip_max"), c("Sex_F0"))
Pro_med_data1_v4 <- Pro_med_data1_v4 %>%
  mutate(
    Phase_group = ifelse(
      Phase %in% c("phase1", "phase2", "phase3"),
      "phase1_3",
      "phase4"
    )
  )
Pro_med_data1_final_F024 <- scale_columns_group(Pro_med_data1_v4, Pro_names, "Phase_group")
####中介效应####
Pro_med_results2_F024 <- run_microbiome_mediation_batch_glm(
  data = Pro_med_data1_final_F024,
  micro_sig = med_results_all_df_Pro_sig,
  covariates =  c(Covariates_all, "Phase_group"),
  predictor_col = "Exposure",
  mediator_col = "Mediator",
  outcome_col = "Outcome",
  exposure_high = 2,
  exposure_low = 0,
  sims = 1000,
  output_file = "Proteins_mediation_results_time.csv"
)

Pro_med_results2_F024 <- format_mediation_table(Pro_med_results2_F024, include_fdr = FALSE)

Pro_med_results2_F024 <- Pro_med_results2_F024 %>%
  mutate(
    
    Exposure = recode(
      Exposure,
      !!!predictor_rename
    ),
    
    Mediator = str_remove(
      Mediator,
      "_fill_log2_z"
    ),
    
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    )
    
  )
Pro_med_results2_F024$Mediator <- toupper(Pro_med_results2_F024$Mediator)
Pro_med_results2_F024 <- merge(Pro_med_results2_F024, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Mediator", by.y = "Protein_names")
Pro_med_results2_F024 <- dplyr::select(Pro_med_results2_F024, -Mediator)
colnames(Pro_med_results2_F024) <- recode(colnames(Pro_med_results2_F024), "hgnc_symbol" = "Mediator")
####*****************************蛋白质富集分析*****************************####
#### 1. 读取 DAVID 导出的 txt 文件 ####
david_file <- "D:/OneDrive/Papers/9 饮茶频率&骨骼肌/Output/output_identify_20260629174917.txt"

raw_lines <- readLines(
  david_file,
  warn = FALSE
)
#### 2. 清理 DAVID 文件中的说明行、空行和分隔线 ####
raw_lines2 <- raw_lines %>%
  .[!grepl("^##", .)] %>%
  .[!grepl("^\\s*$", .)] %>%
  .[!grepl("^-{3,}", .)]

#### 3. 提取表头 ####
header_line <- raw_lines2[
  grepl("^#Term", raw_lines2)
][1]

header <- strsplit(
  sub("^#", "", header_line),
  "\t"
)[[1]]

#### 4. 提取真正的数据行 ####
data_lines <- raw_lines2[
  !grepl("^#Term", raw_lines2)
]

## 只保留列数和表头一致的数据行
n_cols <- length(header)

data_lines <- data_lines[
  sapply(
    strsplit(data_lines, "\t", fixed = TRUE),
    length
  ) == n_cols
]

#### 5. 转换为数据框 ####
kegg_df <- read.delim(
  text = paste(data_lines, collapse = "\n"),
  header = FALSE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(kegg_df) <- header

#### 6. 整理变量类型 ####
kegg_df <- kegg_df %>%
  dplyr::mutate(
    `Input number` = as.numeric(`Input number`),
    `Background number` = as.numeric(`Background number`),
    `P-Value` = as.numeric(`P-Value`),
    `Corrected P-Value` = as.numeric(`Corrected P-Value`),
    
    neg_log10_FDR = -log10(`Corrected P-Value`),
    neg_log10_P = -log10(`P-Value`),
    
    GeneRatio = `Input number` / `Background number`,
    
    Significant = dplyr::case_when(
      `Corrected P-Value` < 0.001 ~ "***",
      `Corrected P-Value` < 0.01  ~ "**",
      `Corrected P-Value` < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::arrange(`Corrected P-Value`)

#### 7. 只保留 FDR < 0.05 的通路用于主图 ####
kegg_sig <- kegg_df %>%
  # dplyr::filter(
  #   `Corrected P-Value` < 0.05
  # ) %>%
  dplyr::mutate(
    Term = factor(
      Term,
      levels = rev(Term)
    ),
    Term_label = stringr::str_wrap(
      as.character(Term),
      width = 35
    )
  )

p_dot <- ggplot(
  kegg_sig,
  aes(
    x = GeneRatio,
    y = reorder(Term_label, GeneRatio)
  )
) +
  geom_point(
    aes(
      size = `Input number`,
      color = neg_log10_FDR
    ),
    alpha = 0.9
  ) +
  geom_text(
    aes(label = Significant),
    color = "black",
    fontface = "bold",
    size = 4,
    vjust = 0.5
  ) +
  scale_color_gradient(
    low = "#DDEBF7",
    high = "#08519C",
    name = expression(-log[10]("FDR"))
  ) +
  scale_size_continuous(
    range = c(4, 9),
    breaks = c(1, 2, 3),
    labels = c("1", "2", "3"),
    name = "Input number"
  ) +
  labs(
    x = "Gene ratio",
    y = NULL,
    title = "KEGG pathway",
    subtitle = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 11,
      color = "grey30",
      hjust = 0
    ),
    axis.text.y = element_text(
      color = "black",
      size = 10
    ),
    axis.text.x = element_text(
      color = "black",
      size = 10
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 12
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.margin = ggplot2::margin(
      t = 10,
      r = 20,
      b = 10,
      l = 10
    )
  )

ggsave(
  "Pro_kegg.png",
  plot = p_dot,
  width = 6,
  height = 5.5,
  dpi = 600
)

####*****************************data10-饮茶&粪便代谢物&骨骼肌*****************************####
data_list <- list(Cov_F0_final, Tea_F0_final, ASMI_wide[,c("ID","WBTOT_FAT_mean_F1234")], ASM_Height_long[ASM_Height_long$Followup %in% c("F4"),], Fecal_met_final) 
data10_v1<- data_list %>% purrr::reduce(dplyr::inner_join, by = "ID") %>%
  dplyr::filter(stats::complete.cases(.))

table(Fecal_met_final$Times)
# 合并握力等
#data10_v2 <- left_join(data10_v1, SPPB_grip_long, by = c("ID" ,"Followup"))
data10_v2 <- left_join(data10_v1, SPPB_grip_long[SPPB_grip_long$Followup == "F3",], by = "ID")

data10_v2 <- left_join(data10_v2, Catechin, by = c("ID" = "ID"))

data10_v3 <- group_variables3(
  data = data10_v2,
  zero_vars = Catechin_vars1,
  tertile_vars = Catechin_vars2,
  suffix = "_T"
)

factor_name<-c("serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T")
idx <- which(names(data10_v3)%in% factor_name)
for(i in idx ){
  data10_v3[[i]] <-  as.factor(data10_v3[[i]])
}
####数据转换
FecalMed_names <- grep("_fill_log2$", colnames(Fecal_met_final), value = TRUE) 

data10_v4 <- scale_columns_group(data10_v3, c("ASM", "ASMI","grip_max"), c("Sex_F0"))
data10_final <- scale_columns(data10_v4, FecalMed_names)
table(data10_final$tea_freq_group)
a <- unique(data10_final$ID)
####lmer-简单调整####
glm_results10_all <- process_glm1(c("tea_freq_group"),
                                  paste0(FecalMed_names,"_z"),
                                   data10_final,
                                   c("Age_F0","Sex_F0"),
                                   2)
glm_results10_all <-glm_results10_all[glm_results10_all$Level == 2,]

glm_results10_all$P_FDR <- p.adjust(glm_results10_all$P_value, method = "fdr")
glm_results10_all_sig <- glm_results10_all[glm_results10_all$P_value < 0.05,]

glm_results10_all$hgnc_symbol <- gsub( "_fill_log2_z$","", glm_results10_all$Outcome)
####可视化#####
# 数据整理
df_plot <- glm_results10_all %>%
  mutate(
    
    # y轴
    log_p = -log10(P_value),
    
    # 分组
    Group = case_when(
      P_value < 0.05 & Estimate > 0 ~ "FDR Sig (Up)",
      P_value < 0.05 & Estimate < 0 ~ "FDR Sig (Down)",
      TRUE ~ "Not Significant"
    )
    
  )

# 火山图
p_volcano <- ggplot(
  df_plot,
  aes(x = Estimate, y = log_p)
) +
  
  # FDR阈值线
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "#CBD5E1",
    linewidth = 0.6
  ) +
  
  # Estimate=0
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    color = "#94A3B8",
    linewidth = 0.6
  ) +
  
  # 所有点
  geom_point(
    aes(fill = Group, size = Group),
    shape = 21,
    color = "white",
    stroke = 0.4,
    alpha = 0.85
  ) +
  
  # 显著基因加黑圈
  geom_point(
    data = filter(df_plot, P_value < 0.05),
    shape = 21,
    color = "#0F172A",
    fill = NA,
    size = 4.8,
    stroke = 1.2
  ) +
  
  # 显著基因标签
  geom_text_repel(
    data = filter(df_plot, P_value < 0.05),
    aes(label = hgnc_symbol),
    
    size = 4,
    fontface = "bold",
    color = "#0F172A",
    
    box_padding = 0.6,
    point_padding = 0.3,
    
    segment.color = "#475569",
    segment.size = 0.5,
    
    min.segment.length = 0,
    max.overlaps = Inf,
    
    force = 5,
    force_pull = 1
  ) +
  
  # 配色
  scale_fill_manual(
    values = c(
      "FDR Sig (Up)" = "#E07A5F",
      "FDR Sig (Down)" = "#2A6F97",
      "Not Significant" = "#F1F5F9"
    )
  ) +
  
  # 点大小
  scale_size_manual(
    values = c(
      "FDR Sig (Up)" = 4.2,
      "FDR Sig (Down)" = 4.2,
      "Not Significant" = 2.5
    )
  ) +
  
  # 主题
  theme_bw() +
  
  theme(
    
    axis.text = element_text(
      size = 11,
      color = "#475569"
    ),
    
    axis.title = element_text(
      size = 13,
      face = "bold",
      color = "#0F172A"
    ),
    
    panel.grid.major = element_line(
      color = "#F8F9FA",
      linewidth = 0.5
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "top",
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 10,
      face = "bold",
      color = "#334155"
    ),
    
    plot.margin = ggplot2::margin(
      20,
      20,
      20,
      20
    )
  ) +
  
  labs(
    x = "β",
    y = expression(bold(-log[10](P))),
    subtitle = NULL
  )

# 导出
ggsave(
  "FecalMet_Volcano_Plot.png",
  plot = p_volcano,
  width = 7.5,
  height = 6.5,
  dpi = 500
)

####lmer-全调整####
glm_results10_all2 <- process_glm1(c("tea_freq_group", "serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T"),
                                   glm_results10_all_sig$Outcome,
                                   data10_final,
                                   Covariates_all,
                                   2)
glm_results10_all2 <- glm_results10_all2[glm_results10_all2$Level == 2,]
glm_results10_all2 <- glm_results10_all2 %>%
  group_by(Predictor) %>%
  mutate(
    P_FDR = p.adjust(P_value, method = "fdr")
  ) %>%
  ungroup()

glm_results10_all2_sig <- glm_results10_all2[glm_results10_all2$P_value < 0.05 ,] #& glm_results10_all2$P_FDR < 0.10

#*************正文写作
a <- glm_results10_all2_sig[glm_results10_all2_sig$Estimate >0 & glm_results10_all2_sig$Predictor == "tea_freq_group",]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)

a <- glm_results10_all2_sig[glm_results10_all2_sig$Estimate >0 & glm_results10_all2_sig$Predictor != "tea_freq_group",]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)

a <- glm_results10_all2_sig[glm_results10_all2_sig$Estimate <0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)
####可视化####
# 1. 载入核心包
library(ggplot2)
library(dplyr)
library(stringr)

# 2. 数据清洗与排版准备
plot_data <- glm_results10_all2 %>%
  mutate(
    # 清洗 Outcome 变量名
    Outcome_clean = str_replace(Outcome, "_fill_log2_z$", ""),
    
    # 全局按照 Estimate 大小排列
    Outcome_clean = reorder(Outcome_clean, Estimate),
    
    # 清洗 Predictor 变量名
    Predictor_clean = ifelse(
      Predictor %in% names(predictor_rename),
      predictor_rename[Predictor],
      Predictor
    ),
    
    Predictor_clean = factor(
      Predictor_clean,
      levels = predictor_rename
    ),
    
    # P < 0.05 标记星号
    Sig_Label = ifelse(P_value < 0.05, "*", "")
  ) %>%
  group_by(Predictor_clean) %>%
  mutate(
    # 每个分面单独计算星号偏移距离
    x_range = max(CI_high, na.rm = TRUE) - min(CI_low, na.rm = TRUE),
    x_pad = ifelse(is.finite(x_range) & x_range > 0, x_range * 0.05, 0.03),
    
    # 星号放在CI外侧
    star_x = ifelse(
      Estimate >= 0,
      CI_high + x_pad,
      CI_low - x_pad
    ),
    
    star_hjust = ifelse(
      Estimate >= 0,
      0,
      1
    )
  ) %>%
  ungroup()

# 动态获取独特的代谢物数量，用于绘制底层斑马纹背景
y_levels <- levels(plot_data$Outcome_clean)
n_levels <- length(y_levels)

# 3. 绘制森林图
p <- ggplot(plot_data, aes(x = Estimate, y = Outcome_clean)) +
  
  # 斑马纹背景
  geom_rect(
    data = data.frame(y = seq(1, n_levels, 2)),
    aes(
      ymin = y - 0.5,
      ymax = y + 0.5,
      xmin = -Inf,
      xmax = Inf
    ),
    fill = "#f8f9fa",
    inherit.aes = FALSE
  ) +
  
  # X=0 基准线
  geom_vline(
    xintercept = 0,
    color = "#cccccc",
    linewidth = 0.5
  ) +
  
  # 95% CI
  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    height = 0.15,
    color = "#2c3e50",
    linewidth = 0.75
  ) +
  
  # 点估计值
  geom_point(
    color = "#3498db",
    size = 3.5,
    shape = 15
  ) +
  
  # P < 0.05 星号
  geom_text(
    aes(
      x = star_x,
      label = Sig_Label,
      hjust = star_hjust
    ),
    color = "#2c3e50",
    size = 5,
    fontface = "bold"
  ) +
  
  # 按 Predictor 分面
  facet_wrap(
    ~ Predictor_clean,
    nrow = 1,
    scales = "free_x"
  ) +
  
  # 给星号留空间
  scale_x_continuous(
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "#e1e4e8",
      fill = NA,
      linewidth = 0.8
    ),
    
    axis.text.y = element_text(
      face = "bold",
      color = "#2c3e50",
      size = 12.5
    ),
    
    axis.text.x = element_text(
      color = "#7f8c8d",
      size = 10
    ),
    
    axis.title.x = element_text(
      margin = ggplot2::margin(t = 15),
      face = "bold",
      color = "#2c3e50",
      size = 16
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      color = "#2c3e50",
      size = 13,
      margin = ggplot2::margin(b = 12)
    ),
    
    panel.spacing = unit(2.5, "lines"),
    
    plot.title = element_text(
      face = "bold",
      size = 16,
      color = "#2c3e50",
      margin = ggplot2::margin(b = 6)
    ),
    
    plot.subtitle = element_text(
      size = 11,
      color = "#7f8c8d",
      margin = ggplot2::margin(b = 20)
    ),
    
    plot.caption = element_text(
      size = 10,
      color = "#7f8c8d",
      hjust = 1
    )
  ) +
  
  labs(
    x = "β (95% CI)",
    y = NULL,
    title = NULL,
    subtitle = NULL,
    caption = NULL
  )

# 4. 保存
ggsave(
  "Tea_fecal_metabolites_allajusted.png",
  plot = p,
  width = 17,
  height = 10,
  dpi = 500
)

####粪便代谢物&骨骼肌####
# lasso_results_ASMI <- lasso_auto(data10_final,"ASMI_z", paste0(FecalMed_names,"_z"))$coef_1se
# lasso_results_ASM <- lasso_auto(data10_final,"ASM_z", paste0(FecalMed_names,"_z"))$coef_1se
# lasso_results_grip <- lasso_auto(data10_final,"grip_max_z", paste0(FecalMed_names,"_z"))$coef_1se
# Var <- unique(c(lasso_results_ASMI$Var, lasso_results_ASM$Var, lasso_results_grip$Var))

glm_results10_all_ASM <- process_glm1(glm_results10_all2_sig$Outcome, 
                                      c("ASM_z","ASMI_z","grip_max_z"),
                                      data10_final, #%>% filter(tea_freq_group %in% c(0, 2)),
                                      c(Covariates_all,"tea_freq_group"), #
                                      2)
glm_results10_all_ASM  <- glm_results10_all_ASM  %>%
  group_by(Outcome) %>%
  mutate(
    P_FDR = p.adjust(P_value, method = "fdr")
  ) %>%
  ungroup()

glm_results10_all_ASM_sig <- glm_results10_all_ASM[glm_results10_all_ASM$P_value < 0.05 ,] #& glm_results10_all_ASM$P_FDR < 0.10

#*************正文写作
a <- glm_results10_all_ASM_sig[glm_results10_all_ASM_sig$Estimate >0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)

####可视化####
plot_df <- glm_results10_all_ASM %>%
  
  mutate(
    
    Predictor = gsub( "_fill_log2_z$","", Predictor),
    
    ## outcome
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    ),
    
    ## significance stars
    sig = case_when(
      P_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) 

p <- plot_forest_facet(
  data = plot_df,
  estimate_col = "Estimate",
  ci_low_col = "CI_low",
  ci_high_col = "CI_high",
  pathway_col = "Predictor", #自变量
  outcome_col = "Outcome", #分面
  p_col = "P_value",
  y_text_face = "plain",
  sig_rule = "p_only"
)


ggsave("FecalMet_ASM.png", plot = p, width = 13, height = 6, dpi = 500)
####有时间顺序的中介效应####
####中介效应分析
FecalMet_sig_pos <- glm_results10_all2_sig[glm_results10_all2_sig$Estimate >0,] %>%
  
  inner_join(
    
    glm_results10_all_ASM_sig[glm_results10_all_ASM_sig$Estimate >0,] %>%
      
      dplyr::select(
        Predictor,
        Outcome
      ) %>%
      
      distinct(),
    
    by = c(
      "Outcome" = "Predictor"
    )
  ) 

FecalMet_sig_neg <- glm_results10_all2_sig[glm_results10_all2_sig$Estimate <0,] %>%
  
  inner_join(
    
    glm_results10_all_ASM_sig[glm_results10_all_ASM_sig$Estimate <0,] %>%
      
      dplyr::select(
        Predictor,
        Outcome
      ) %>%
      
      distinct(),
    
    by = c(
      "Outcome" = "Predictor"
    )
  )

FecalMet_sig2 <- rbind(FecalMet_sig_pos, FecalMet_sig_neg)

colnames(FecalMet_sig2) <- recode(colnames(FecalMet_sig2),
                                  "Outcome.y" = "Outcome_muscle",
                                  "Outcome" = "Mediator")

FecalMet_med_results2 <- run_microbiome_mediation_batch_glm(
  data = data10_final,
  micro_sig = FecalMet_sig2,
  covariates =  Covariates_all,
  predictor_col = "Predictor",
  mediator_col = "Mediator",
  outcome_col = "Outcome_muscle",
  exposure_high = 2,
  exposure_low = 0,
  sims = 1000,
  output_file = "FecalMetteins_mediation_results_time.csv"
)

FecalMet_med_results2_sig <- FecalMet_med_results2[FecalMet_med_results2$ACME_P_value < 0.05,]

FecalMet_med_results2_output <- format_mediation_table(FecalMet_med_results2, include_fdr = FALSE)
FecalMet_med_results2_output <- FecalMet_med_results2_output %>%
  mutate(
    
    Exposure = recode(
      Exposure,
      !!!predictor_rename
    ),
    
    Mediator = str_remove(
      Mediator,
      "_fill_log2_z"
    ),
    
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    )
    
  )

writexl::write_xlsx(FecalMet_med_results2_output, "FecalMet_med_results2_output.xlsx")
####*****************************菌群和功能的相关性*****************************####
#### 1. 数据集合并 ####
Micro_Func <- merge(
  Micro_selected_filter_clr,
  Micro_func_filter_clr,
  by = c("ID", "Times")
)

table(Micro_Func$Times)

#### 2. 提取变量名：去掉 Predictor 末尾的 "_z" ####
pred7 <- gsub(
  "_z$",
  "",
  Lmer_results7_all_ASM_sig$Predictor
)

pred8 <- gsub(
  "_z$",
  "",
  Lmer_results8_all_ASM_sig$Predictor
)

## 去重
pred7 <- unique(pred7)
pred8 <- unique(pred8)

#### 3. 提取 Micro_Func 中需要的列 ####
cols_use <- unique(c(
  "ID",
  "Times",
  pred7,
  pred8
))

## 检查变量是否存在
missing_cols <- setdiff(
  cols_use,
  colnames(Micro_Func)
)

if(length(missing_cols) > 0){
  cat(
    "以下变量在 Micro_Func 中不存在，将自动跳过：\n",
    paste(missing_cols, collapse = "\n"),
    "\n\n"
  )
}

## 只保留真实存在的列
cols_use2 <- intersect(
  cols_use,
  colnames(Micro_Func)
)

Micro_Func_selected <- Micro_Func[, cols_use2, drop = FALSE]

## 同步更新 pred7 和 pred8，只保留存在于 Micro_Func_selected 的变量
pred7 <- intersect(
  pred7,
  colnames(Micro_Func_selected)
)

pred8 <- intersect(
  pred8,
  colnames(Micro_Func_selected)
)

## 如果没有可分析变量，停止
if(length(pred7) == 0){
  stop("pred7 中没有任何变量存在于 Micro_Func_selected 中，请检查变量名。")
}

if(length(pred8) == 0){
  stop("pred8 中没有任何变量存在于 Micro_Func_selected 中，请检查变量名。")
}

cat("pred7 变量数：", length(pred7), "\n")
cat("pred8 变量数：", length(pred8), "\n")

#### 4. 单个变量对的相关分析函数 ####
run_one_cor <- function(df, x, y, method = "spearman"){
  
  tmp <- df %>%
    dplyr::select(
      dplyr::all_of(c(x, y))
    ) %>%
    tidyr::drop_na()
  
  n_use <- nrow(tmp)
  
  ## 样本量太小，或者变量没有变异，返回 NA
  if(
    n_use < 3 ||
    sd(tmp[[x]], na.rm = TRUE) == 0 ||
    sd(tmp[[y]], na.rm = TRUE) == 0
  ){
    
    return(
      data.frame(
        Var_7 = x,
        Var_8 = y,
        N = n_use,
        R = NA_real_,
        P_value = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  test_res <- suppressWarnings(
    cor.test(
      tmp[[x]],
      tmp[[y]],
      method = method,
      exact = FALSE
    )
  )
  
  data.frame(
    Var_7 = x,
    Var_8 = y,
    N = n_use,
    R = as.numeric(test_res$estimate),
    P_value = test_res$p.value,
    stringsAsFactors = FALSE
  )
}

#### 5. 分 Times 做所有 pairwise correlation ####
time_levels_use <- c("F2", "F3")

corr_results <- purrr::map_dfr(
  time_levels_use,
  function(time_i){
    
    df_i <- Micro_Func_selected %>%
      dplyr::filter(
        Times == time_i
      )
    
    cat("正在分析 Times =", time_i, "；样本量 =", nrow(df_i), "\n")
    
    purrr::map_dfr(
      pred7,
      function(x){
        
        purrr::map_dfr(
          pred8,
          function(y){
            
            res <- run_one_cor(
              df = df_i,
              x = x,
              y = y,
              method = "spearman"
            )
            
            res$Times <- time_i
            
            return(res)
          }
        )
      }
    )
  }
)

#### 6. FDR 校正和显著性标记 ####
corr_results <- corr_results %>%
  dplyr::group_by(Times) %>%
  dplyr::mutate(
    P_FDR = p.adjust(P_value, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Sig = dplyr::case_when(
      is.na(P_FDR) ~ "",
      P_FDR < 0.001 ~ "***",
      P_FDR < 0.01  ~ "**",
      P_FDR < 0.05  ~ "*",
      TRUE ~ ""
    ),
    
    ## 所有格子都显示 r 值，显著的额外加星号
    R_label = dplyr::case_when(
      is.na(R) ~ "",
      TRUE ~ paste0(sprintf("%.2f", R), Sig)
    )
  )

#### 7. 提取显著结果 ####
corr_sig <- corr_results %>%
  dplyr::filter(
    !is.na(P_FDR),
    P_FDR < 0.05
  ) %>%
  dplyr::arrange(
    Times,
    P_FDR
  )

round(min(corr_sig$R),3)
round(max(corr_sig$R),3)
round(min(corr_sig$P_FDR),3)
round(max(corr_sig$P_FDR),3)

writexl::write_xlsx(corr_sig, "corr_sig.xlsx")
#### 8. 美化变量名函数 ####
clean_label <- function(x){
  x %>%
    gsub("^s__", "", .) %>%
    gsub("^g__", "", .) %>%
    gsub("^f__", "", .) %>%
    gsub("^p__", "", .) %>%
    gsub("_clr$", "", .) %>%
    gsub("_", " ", .) %>%
    stringr::str_wrap(width = 36)
}

#### 9. 设置热图变量顺序 ####
corr_results_plot <- corr_results %>%
  dplyr::mutate(
    Var_7 = factor(
      Var_7,
      levels = rev(pred7)
    ),
    Var_8 = factor(
      Var_8,
      levels = pred8
    )
  )

#### 11. 所有相关的气泡图：显著和不显著都显示 ####
corr_dot_plot_data <- corr_results %>%
  dplyr::filter(
    !is.na(R)
  ) %>%
  dplyr::mutate(
    Var_7 = factor(
      Var_7,
      levels = rev(pred7)
    ),
    Var_8 = factor(
      Var_8,
      levels = pred8
    ),
    Significant = dplyr::case_when(
      !is.na(P_FDR) & P_FDR < 0.05 ~ "FDR < 0.05",
      TRUE ~ "NS"
    )
  )

p_corr_dot <- ggplot(
  corr_dot_plot_data,
  aes(
    x = Var_8,
    y = Var_7,
    size = abs(R),
    color = R
  )
) +
  geom_point(
    alpha = 0.85
  ) +
  geom_text(
    aes(label = Sig),
    size = 4,
    fontface = "bold",
    color = "black",
    vjust = 0.5
  ) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "grey95",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman r"
  ) +
  scale_size_continuous(
    range = c(2.5, 8),
    name = "|r|"
  ) +
  scale_x_discrete(
    labels = clean_label
  ) +
  scale_y_discrete(
    labels = clean_label
  ) +
  facet_wrap(
    ~ Times,
    ncol = 1
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      size = 11,
      color = "grey30"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8,
      color = "black",
      lineheight = 0.85,
      margin = ggplot2::margin(t = 12)
    ),
    axis.text.y = element_text(
      size = 9.5,
      color = "black",
      face = "italic"
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 12,
      margin = ggplot2::margin(t = 10)
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 12,
      margin = ggplot2::margin(r = 10)
    ),
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    strip.background = element_rect(
      fill = "grey92",
      color = NA
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
#### 13. 自动设置导出图片尺寸 ####
dot_width <- max(
  10,
  length(pred8) * 0.7
)

dot_height <- max(
  8,
  length(pred7) * 0.45 * length(time_levels_use)
)

#### 15. 导出气泡图 ####
ggsave(
  filename = "Micro_Func_pairwise_correlation_dotplot_F2_F3.pdf",
  plot = p_corr_dot,
  width = dot_width,
  height = dot_height,
  limitsize = FALSE
)

ggsave(
  filename = "Micro_Func_pairwise_correlation_dotplot_F2_F3.png",
  plot = p_corr_dot,
  width = dot_width,
  height = dot_height,
  dpi = 300,
  limitsize = FALSE
)

####*****************************导出附表*****************************####
mediation_micro <- read_excel("D:/OneDrive/Papers/9 饮茶频率&骨骼肌/Output/mediation_results_all_micro.xlsx", sheet=1)
mediation_micro <- format_mediation_table(mediation_micro)
#mediation_micro_sig <- mediation_micro[mediation_micro$`P value` <0.05,]
# 数据整理
mediation_micro_sig <- mediation_micro %>%
  mutate(
    
    Exposure = recode(
      Exposure,
      !!!predictor_rename
    ),
    
    Mediator = str_remove(Mediator, "s__") %>% str_remove("_clr_z") %>% str_replace_all("_", " "),
    
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    )
    
  )

mediation_func <- read_excel("D:/OneDrive/Papers/9 饮茶频率&骨骼肌/Output/med_results_all_df_Func.xlsx", sheet=1)
mediation_func <- format_mediation_table(mediation_func)
# 数据整理
mediation_func_sig <- mediation_func%>%
  mutate(
    
    Exposure = recode(
      Exposure,
      !!!predictor_rename
    ),
    
    Mediator = str_remove(Mediator, "s__") %>% str_remove("_clr_z") %>% str_replace_all("_", " "),
    
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    )
    
  )

mediation_pro <- read_excel("D:/OneDrive/Papers/9 饮茶频率&骨骼肌/Output/med_results_all_df_Pro.xlsx", sheet=1)
mediation_pro <- format_mediation_table(mediation_pro)
# 数据整理
mediation_pro_sig <- mediation_pro %>%
  mutate(
    
    Exposure = recode(
      Exposure,
      !!!predictor_rename
    ),
    
    Mediator = str_remove(
      Mediator,
      "_fill_log2_z"
    ),
    
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    )
    
  )
mediation_pro_sig$Mediator <- toupper(mediation_pro_sig$Mediator)
mediation_pro_sig <- merge(mediation_pro_sig, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Mediator", by.y = "Protein_names")

mediation_pro_sig <- dplyr::select(mediation_pro_sig, -Mediator)
colnames(mediation_pro_sig) <- recode(colnames(mediation_pro_sig),
                                      "hgnc_symbol" = "Mediator")

writexl::write_xlsx(
  list(
    TableS6  = Lmer_results7_all_output,
    TableS7  = cv_res_results7_tea_clean,
    TableS8  = Lmer_results8_all_output,
    TableS9  = Lmer_results9_all_output,
    TableS10  = rbind(Lmer_model_compare_micro_output, Lmer_model_compare_func_output, Lmer_model_compare_pro_output),
    TableS11  = rbind(mediation_micro_sig, mediation_func_sig, mediation_pro_sig), 
    TableS12  = rbind(Micro_med_results2_F024, Fun_med_results2_F024, Pro_med_results2_F024, FecalMet_med_results2_output),
    TableS13  = rbind(Micro_med_results2_F023, Fun_med_results2_F023, Pro_med_results2_F023)
  ),
  path = "Supplemental Tables2.xlsx"
)

a <- rbind(Micro_med_results2, Fun_med_results2, Pro_med_results2)
a1 <- a[a$p < 0.05,]

a <- rbind(mediation_micro_sig, mediation_func_sig, mediation_pro_sig)
a1 <- a[a$p_FDR < 0.05,]

a <- rbind(Lmer_model_compare_micro_output, Lmer_model_compare_func_output, Lmer_model_compare_pro_output)
max(a$Delta_Marginal_R2_percent)
min(a$Delta_Marginal_R2_percent)


writexl::write_xlsx(
  list(
    Table_1  = Lmer_results8_all_sig,
    Table_2  = Lmer_results8_all_ASM_sig,
    Table_3  = Lmer_results8_all2_sig
  ),
  path = "Func_results.xlsx"
)

writexl::write_xlsx(
  list(
    Table_1  = Lmer_results9_all_output_sig,
    Table_2  = Lmer_results9_all_ASM_sig,
    Table_3  = Lmer_results9_all2_sig
  ),
  path = "Pro_results.xlsx"
)

writexl::write_xlsx(
  list(
    Table_1  = mediation_results_all_micro_final[mediation_results_all_micro_final$ACME_p <0.05,],
    Table_2  = mediation_results_all_Func[mediation_results_all_Func$ACME_p <0.05,],
    Table_3  = mediation_results_all_pro[mediation_results_all_pro$ACME_p <0.05,]
  ),
  path = "Mediation_summary.xlsx"
)