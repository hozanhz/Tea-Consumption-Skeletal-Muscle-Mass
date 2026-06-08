
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
    p_col = "P_value",
    sig_level1 = 0.05,
    sig_level2 = 0.01,
    sig_level3 = 0.001,
    ncol = 3,
    point_size = 5.8,
    ci_width = 1.4,
    base_size = 15,
    low_color = "#4C78A8",
    mid_color = "white",
    high_color = "#D65F5F",
    y_text_face = "italic"   # 新增参数
) {
  
  library(ggplot2)
  library(dplyr)
  library(rlang)
  library(grid)
  
  ## 1. 非标准求值
  
  estimate_sym <- sym(estimate_col)
  ci_low_sym   <- sym(ci_low_col)
  ci_high_sym  <- sym(ci_high_col)
  pathway_sym  <- sym(pathway_col)
  outcome_sym  <- sym(outcome_col)
  p_sym        <- sym(p_col)
  
  ## 2. 显著性星号
  
  plot_df <- data %>%
    
    mutate(
      
      sig = case_when(
        
        !!p_sym < sig_level3 ~ "***",
        !!p_sym < sig_level2 ~ "**",
        !!p_sym < sig_level1 ~ "*",
        TRUE ~ ""
      )
    )
  
  ## 3. 绘图
  
  p <- ggplot(
    plot_df,
    aes(
      x = !!estimate_sym,
      y = !!pathway_sym
    )
  ) +
    
    ## CI
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
    
    ## 点
    geom_point(
      aes(
        fill = !!estimate_sym
      ),
      shape = 21,
      size = point_size,
      color = "black",
      stroke = 0.25
    ) +
    
    ## 星号
    geom_text(
      aes(label = sig),
      color = "black",
      size = 4.5,
      fontface = "bold"
    ) +
    
    ## 0线
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey70",
      linewidth = 0.5
    ) +
    
    ## facet
    facet_wrap(
      vars(!!outcome_sym),
      ncol = ncol
    ) +
    
    ## fill颜色
    scale_fill_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = 0,
      name = "Estimate"
    ) +
    
    ## line颜色
    scale_color_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = 0
    ) +
    
    ## 删除重复legend
    guides(
      color = "none"
    ) +
    
    ## 标签
    labs(
      x = "β (95% CI)",
      y = NULL
    ) +
    
    ## 主题
    theme_minimal(
      base_size = base_size
    ) +
    
    theme(
      
      panel.grid.major.y =
        element_blank(),
      
      panel.grid.minor =
        element_blank(),
      
      axis.text.y =
        element_text(
          size = 12,
          color = "black",
          face = y_text_face   # 参数化
        ),
      
      axis.text.x =
        element_text(
          size = 13,
          color = "black"
        ),
      
      axis.title =
        element_text(
          face = "bold",
          size = 16
        ),
      
      strip.text =
        element_text(
          face = "bold",
          size = 15
        ),
      
      strip.background =
        element_rect(
          fill = "#F3F3F3",
          color = NA
        ),
      
      legend.position =
        "right",
      
      panel.spacing =
        unit(1.5, "lines")
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
####变量名转换####
predictor_rename <- c(
  "tea_freq_group" = "Tea consumption",
  "serum_I_catechin_F0_T" = "Catechin",
  "serum_I_EGC_F0_T" = "Epigallocatechin",
  "serum_I_EGCG_F0_T" = "Epigallocatechin gallate",
  "serum_I_epicatechin_F0_T" = "Epicatechin",
  "serum_I_ECG_F0_T" = "Epicatechin gallate"
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
colnames(Tea)
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
Tea_F2_filter <- Tea %>%
  
  mutate(
    
    每周泡茶次数_F2 =
      coalesce(每周冲茶次数_F2,
               平均每周泡茶次数_F2),
    
    过去一年茶叶总量_斤_F2 =
      if_else(
        过去一年是否经常喝茶_F2 == 0,
        0,
        coalesce(过去一年茶叶量斤_F2, 0) +
          coalesce(过去一年茶叶量两_F2, 0) / 10
      )
  ) %>%
  
  rowwise() %>%
  mutate(
    常饮用茶种类_F2 =
      str_c(
        na.omit(c(
          ifelse(绿茶_F2 == 1, "绿茶", NA),
          ifelse(红茶_F2 == 1, "红茶", NA),
          ifelse(乌龙_F2 == 1, "乌龙茶", NA),
          ifelse(其他茶_F2 == 1, 其他茶叶名称_F2, NA)
        )),
        collapse = "；"
      )
  ) %>%
  ungroup() %>%
  
  dplyr::select(
    CODE_F2,
    过去一年是否经常喝茶_F2,
    每周泡茶次数_F2,
    喝茶浓淡_F2,
    过去一年茶叶总量_斤_F2,
    常饮用茶种类_F2
  )
head(Tea_F1_filter$CODE_F1)
colnames(Tea_F2_filter)

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

# 转换为长数据
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
a <- Micro_selected_filter$filtered_data
#****************clr转换，因为clr是按行转换，所以跟整体分布无关
Micro_selected_filter_clr <- perform_micro_clr(Micro_selected_filter$filtered_data, exclude_cols = c("ID", "Times"))
Micro_selected_filter_clr <-Micro_selected_filter_clr[grepl("^NL", Micro_selected_filter_clr$ID),]

a <- Micro_selected_filter_clr[Micro_selected_filter_clr$Times %in% c("F2","F3"),]
table(a$Times)
a$identified <- paste0(a$ID,a$Times)
a1 <- Micro_func_filter_clr
a1$identified <- paste0(a1$ID, a1$Times)
setdiff(unique(a$identified ), unique(a1$identified))
table(a$Times)
####多样性####
Diversity <- as.data.frame(read.csv("D:/OneDrive/Papers/2.2 中大队列/database/大队列数据库整理/原始数据库/GNHS_metaphlan3_diversity_2664_221028_VF.csv",header = TRUE))
Diversity <- Diversity %>%
  mutate(
    Times = substr(id, 1, 2),
    ID    = substr(id, 3, nchar(id))
  ) %>%
  relocate(ID, Times, .after = 1)
####微生物功能####
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
####*********混杂定义*********####
Covariates_all = c("Age_F0","Sex_F0","WBTOT_FAT_mean_F1234","Total_score","Income_F0","Education_F0","coffee","Met_mean","Protein_mean","Alcohol_F0","Smoke_F0","Calcium_F0","Vitamin_F0","Disease_F0","Estrogen_F0","Menopause_F0","Oil_sup_follow")
Covariates_female = c("Age_F0","WBTOT_FAT_mean_F1234","Total_score","Income_F0","Education_F0","coffee","Met_mean","Protein_mean","Alcohol_F0","Calcium_F0","Vitamin_F0","Disease_F0","Estrogen_F0","Menopause_F0","Oil_sup_follow")
Covariates_male = c("Age_F0","WBTOT_FAT_mean_F1234","Total_score","Income_F0","Education_F0","coffee","Met_mean","Protein_mean","Alcohol_F0","Smoke_F0","Calcium_F0","Vitamin_F0","Disease_F0","Oil_sup_follow")

# 线形混合效应混杂
Covariates_all_lmer = c("Age","Sex_F0","WBTOT_FAT","Total_score","Income_F0","Education_F0","coffee","Met_mean","Protein_mean","Alcohol_F0","Smoke_F0","Calcium_F0","Vitamin_F0","Disease_F0","Estrogen_F0","Menopause_F0","Oil_sup_follow")

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

####lmer####
Lmer_results1_all <- process_lmer(c("tea_freq_group","serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T"),
                                  c("ASM_z","ASMI_z", "grip_max_z","gait_speed_z"), #,"SPPB_total"
                                  data1_final, 
                                  Covariates_all_lmer, 
                                  2)
Lmer_results1_all_sig <- Lmer_results1_all[Lmer_results1_all$P_value < 0.05,]
table(Lmer_results1_all$Level)
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
  "p_interaction.png",
  plot = p,
  width = 12,
  height = 12,
  dpi = 500
)
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
#### 仅保留饮茶者
colnames(Tea_F0_filter)
tea_drinkers <- Tea_F0_filter %>%
  
  filter(
    是否喝茶_F0 == 1
  )
tea_drinkers <- tea_drinkers[complete.cases(tea_drinkers) & tea_drinkers$累计喝茶总年限_F0 != 0,]
min_amount <- min(
  tea_drinkers$过去一年茶叶总量_斤_F0,
  na.rm = TRUE
)

min_duration <- min(
  tea_drinkers$累计喝茶总年限_F0,
  na.rm = TRUE
)
#### A. Tea type distribution

tea_type_long <- tea_drinkers %>%
  dplyr::select(常饮用茶种类_F0) %>%
  mutate(ID = row_number()) %>%
  separate_rows(常饮用茶种类_F0, sep = "[；，、]") %>%
  mutate(
    常饮用茶种类_F0 = str_trim(常饮用茶种类_F0),
    常饮用茶种类_F0 = case_when(
      # 普洱类合并
      str_detect(常饮用茶种类_F0, "普|普洱|白茶普洱|花茶普洱|普洱花茶|普洱苦丁|普洱菊花|普洱龙井") ~ "普洱",
      # 花茶类合并
      str_detect(常饮用茶种类_F0, "花|花茶|茉|菊|菊花|菊花茶") ~ "花茶",
      # 苦丁茶合并
      str_detect(常饮用茶种类_F0, "苦|苦丁茶") ~ "苦丁茶",
      # 凤凰单丛合并
      str_detect(常饮用茶种类_F0, "凤凰|单") ~ "单丛",
      # 其他明确保留
      TRUE ~ 常饮用茶种类_F0
    )
  ) %>%
  filter(!is.na(常饮用茶种类_F0), 常饮用茶种类_F0 != "")


tea_type_summary <- tea_type_long %>%
  
  count(常饮用茶种类_F0) %>%
  
  mutate(
    Percent = 100 * n / sum(n),
    
    # 英文翻译 + 格式化百分比（保留两位小数）
    茶种类_EN = case_when(
      常饮用茶种类_F0 == "乌龙茶"   ~ "Oolong tea",
      常饮用茶种类_F0 == "决明子"   ~ "Cassia seed tea",
      常饮用茶种类_F0 == "减肥茶（"   ~ "Slimming tea",
      常饮用茶种类_F0 == "单丛"     ~ "Dancong tea",
      常饮用茶种类_F0 == "大麦茶"   ~ "Barley tea",
      常饮用茶种类_F0 == "山楂茶"   ~ "Hawthorn tea",
      常饮用茶种类_F0 == "普洱"     ~ "Pu-erh tea",
      常饮用茶种类_F0 == "红茶"     ~ "Black tea",
      常饮用茶种类_F0 == "绿茶"     ~ "Green tea",
      常饮用茶种类_F0 == "罗布麻茶" ~ "Luobuma tea",
      常饮用茶种类_F0 == "花茶"     ~ "Flower tea",
      常饮用茶种类_F0 == "苦丁茶"   ~ "Kuding tea",
      常饮用茶种类_F0 == "藏茶"     ~ "Tibetan tea",
      常饮用茶种类_F0 == "银杏叶"   ~ "Ginkgo leaf tea",
      TRUE ~ 常饮用茶种类_F0
    ),
    
    # 格式化百分比：保留两位小数
    Percent_label = sprintf("%.2f%%", Percent),
    
    # 用于排序（按百分比，保留中文原名或英文均可）
    常饮用茶种类_F0 = fct_reorder(常饮用茶种类_F0, Percent)
  )
table(tea_type_summary$常饮用茶种类_F0)

p1 <- ggplot(
  tea_type_summary,
  aes(x = Percent, y = 茶种类_EN)  # 这里改成英文
) +
  geom_col(fill = "#238B45", width = 0.75) +
  geom_text(
    aes(label = Percent_label),  # 使用格式化后的标签
    hjust = -0.2,
    size = 4
  ) +
  expand_limits(x = max(tea_type_summary$Percent) * 1.15) +
  labs(
    title = "A. Tea type",
    x = "Percentage (%)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

#### B. Tea amount

p2 <- ggplot(
  tea_drinkers,
  aes(
    x = 过去一年茶叶总量_斤_F0
  )
) +
  
  geom_density(
    linewidth = 1.2,
    color = "#E76F51",
    fill = "#E76F51",
    alpha = 0.25,
    trim = TRUE
  ) +
  
  coord_cartesian(
    xlim = c(
      min_amount,
      NA
    )
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  
  labs(
    title = "B. Tea amount",
    x = "Tea amount (jin/year)",
    y = "Density"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )
min(tea_drinkers$累计喝茶总年限_F0)

#### C. Brewing strength

p3 <- ggplot(
  tea_drinkers,
  aes(
    x =
      factor(
        喝茶浓淡_F0,
        levels = c(1,2,3,4,5),
        labels = c(
          "Strongest",
          "Strong",
          "Moderate",
          "Weak",
          "Weakest"
        )
      )
  )
) +
  
  geom_bar(
    fill = "#457B9D"
  ) +
  
  labs(
    title = "C. Brewing strength",
    x = NULL,
    y = "Count"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    axis.text.x =
      element_text(
        angle = 25,
        hjust = 1
      )
  )

#### D. Tea drinking duration

p4 <- ggplot(
  tea_drinkers,
  aes(
    x = 累计喝茶总年限_F0
  )
) +
  
  geom_density(
    linewidth = 1.2,
    color = "#6A4C93",
    fill = "#6A4C93",
    alpha = 0.25,
    trim = TRUE
  ) +
  
  coord_cartesian(
    xlim = c(
      min_duration,
      NA
    )
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  
  labs(
    title = "D. Tea drinking duration",
    x = "Duration (years)",
    y = "Density"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

#### Combine figure

p_sup_tea <-
  
  (p1 | p2) /
  (p3 | p4) +
  
  plot_annotation(
    title = NULL
  )

ggsave(
  "Supplementary_Figure_Tea_Characteristics.png",
  p_sup_tea,
  width = 12,
  height = 10,
  dpi = 600
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
      str_detect(tea, "普|耳|饵|菊普|生普|功夫茶普洱|普洱花茶|普洱菊花|普洱龙井|普洱毛尖|普洱学菊|普洱玫瑰花|普洱绞股蓝|普洱单木丛|普洱，|普洱/|普洱田七|小青柑|小金柑|白茶普洱|花茶普洱|普洱苦丁") ~ "普洱",
      
      # 乌龙茶类（不包括单丛）
      str_detect(tea, "乌龙茶|乌龙$|大红袍|凤凰茶|清茶|水仙") ~ "乌龙茶",
      
      # 单丛类
      str_detect(tea, "单丛|单松|单枞|枞茶|凤凰|单") ~ "单丛",
      
      # 绿茶类
      str_detect(tea, "绿茶|龙井") ~ "绿茶",
      
      # 红茶类
      str_detect(tea, "^红茶$") ~ "红茶",
      
      # 花茶类
      str_detect(tea, "花|花茶|菊花|茉莉|香花茶|花茶芦荟菊|菊花茶寿眉|昆仑雪菊|玫瑰花|枸杞红枣|陈皮茶") ~ "花茶",
      
      # 黑茶类
      str_detect(tea, "^黑茶$|藏茶") ~ "黑茶",
      
      # 苦丁茶类
      str_detect(tea, "苦丁茶|苦$") ~ "苦丁茶",
      
      # 罗布麻类
      str_detect(tea, "罗布麻") ~ "罗布麻茶",
      
      # 绞股蓝类
      str_detect(tea, "绞股蓝|绞胶蓝|搞股兰") ~ "绞股蓝茶",
      
      # 丹参类
      str_detect(tea, "丹参|丹心|丹参保心|天草丹参") ~ "丹参茶",
      
      # 灵芝类
      str_detect(tea, "灵芝") ~ "灵芝茶",
      
      # 谷物茶类
      str_detect(tea, "大麦茶|荞麦茶|苦荞麦|苦麦茶") ~ "谷物茶",
      
      # 保健茶类
      str_detect(tea, "保健茶|养生茶|降压茶") ~ "保健茶",
      
      # 牛蒡茶
      str_detect(tea, "牛蒡|牛膀") ~ "牛蒡茶",
      
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
    
    # 低 catechin
    "红茶"     = 1,
    "普洱"     = 1,
    "黑茶"     = 1,
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
    # ⭐ 剔除：茶种类为纯数字的行
    filter(
      !str_detect(as.character(.data[[type_col]]), "^[0-9]+$")
    ) %>%
    # ⭐ 剔除：喝茶但茶种类为空的行
    filter(
      !(as.numeric(.data[[drink_col]]) == 1 & 
          (is.na(.data[[type_col]]) | .data[[type_col]] == "" | str_trim(.data[[type_col]]) == ""))
    ) %>%
    transmute(
      ID = if(is.null(remove_prefix)) {
        .data[[code_col]]
      } else {
        sub(remove_prefix, "", .data[[code_col]])
      },
      visit = visit,
      tea_drink = as.numeric(.data[[drink_col]]),
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
        tea_drink == 1 ~ sapply(.data[[type_col]], get_tea_type_score),
        TRUE ~ NA_real_
      ),
      tea_duration = if(is.null(duration_col)) {
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
  
  strength_col = "喝茶浓淡_F2",
  
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

#### duration 只来自 baseline

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

#### 时间变量

Tea_long <- Tea_long %>%
  
  mutate(
    
    time = c(
      F0 = 0,
      F1 = 1,
      F2 = 2,
      F3 = 3
    )[visit]
  )

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
colnames(Tea_long)
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
         P_value = P_value_formatted, 
         P_FDR = P_FDR_formatted) %>%
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
####菌群&骨骼肌####
Lmer_results7_all_ASM <- process_lmer(Lmer_results7_all_FDR_sig$Outcome, 
                                      c("ASM_z","ASMI_z","grip_max_z"),
                                      data7_final,# %>% filter(tea_freq_group %in% c(0, 2)),
                                      c(Covariates_all_lmer,c("Fiber_soluble_mean", "tea_freq_group")), #
                                      2)
Lmer_results7_all_ASM_sig <- Lmer_results7_all_ASM[Lmer_results7_all_ASM$P_value < 0.05,]
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
      P_value < 0.05  ~ "*",
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
  p_col = "P_value"
)

ggsave("Micro_ASM.png", plot = p, width = 13, height = 16, dpi = 500)
####lmer-全调整####
Lmer_results7_all2 <- process_lmer(c("tea_freq_group", "serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T"),
                                   unique(Lmer_results7_all_ASM_sig$Predictor),
                                   data7_final,
                                   c(setdiff(Covariates_all_lmer,c("Protein_mean")),c("Fiber_soluble_mean")),
                                   2)
Lmer_results7_all2 <- Lmer_results7_all2[Lmer_results7_all2$Level == 2,]
Lmer_results7_all2_sig <- Lmer_results7_all2[Lmer_results7_all2$P_value < 0.05,]
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
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    
    # ③ 【核心新颖视觉】：将 Direction 与显著性融合
    # 不显著的统一归为 "Not Significant" 组，从而在后面单独对其剥离色彩、进行高级灰色降阶
    Sig_Group = case_when(
      P_value >= 0.05 ~ "Not Significant",
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
  geom_text(aes(label = Sig_Label, color = Sig_Group), 
            vjust = -0.4, hjust = 0.5, size = 4.5, fontface = "bold", show.legend = FALSE) +
  
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
       width = 19,
       height = 16,
       dpi = 600)

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

unique(Micro_sig$Mediator)
#Micro_sig_filter <- Micro_sig[!Micro_sig$Predictor %in% c("tea_freq_group"), ]
# Micro_sig_filter <- Micro_sig[Micro_sig$Mediator %in% c("s__Eubacterium_ventriosum_clr_z","s__Alistipes_shahii_clr_z",
#                                                         "s__Gemmiger_formicilis_clr_z", "s__Veillonella_parvula_clr_z",
#                                                         "s__Phocaeicola_massiliensis_clr_z", "s__Ruminococcus_bicirculans_clr_z"
#                                                      ), ]
# 
# Micro_sig_filter <- Micro_sig[Micro_sig$Mediator %in% c("s__Eubacterium_ventriosum_clr_z","s__Alistipes_shahii_clr_z",
#                                                         "s__Gemmiger_formicilis_clr_z", "s__Veillonella_parvula_clr_z"
# ), ]
# 
# Micro_sig_filter <- Micro_sig[Micro_sig$Mediator %in% c("s__Clostridium_leptum_clr_z", "s__Alistipes_shahii_clr_z", "s__Gemmiger_formicilis_clr_z",
#                                                         "s__Blautia_hansenii_clr_z"), ]

Micro_sig_filter <- Micro_sig[Micro_sig$Mediator == "s__Bacteroides_fragilis_clr_z" & Micro_sig$Predictor == "tea_freq_group" & Micro_sig$Outcome_muscle == "ASMI_z", ]


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


####可视化####
mediation_results_all_micro_final <- read_excel("D:/OneDrive/Papers/9 饮茶频率&骨骼肌/Output/mediation_results_all_micro.xlsx", sheet=1)
mediation_results_all_micro_final$PropMediated <- mediation_results_all_micro_final$PropMediated * 100
colnames(mediation_results_all_micro_final)
# 数据整理
plot_data_species <- mediation_results_all_micro_final %>%
  
  filter(
    ACME_p < 0.05,
    ACME > 0
  ) %>%
  
  dplyr::select(
    Exposure,
    Mediator,
    Outcome,
    ACME,
    ACME_low,
    ACME_high,
    ACME_p,
    PropMediated
  ) %>%
  
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

p_mediation_species <- plot_mediation_alluvial(
  data = plot_data_species,
  title = "Mediation Effects: Tea Consumption/Biomarkers → Gut Microbial Species → Muscle",
  fill_low = "#D8F0D2",
  fill_high = "#238B45" #
)

pdf(
  "Species_mediation.pdf",
  width = 11,
  height = 4
)
print(p_mediation_species)
dev.off()

#************正文写作
a <- mediation_results_all_micro_final[mediation_results_all_micro_final$ACME_p < 0.05,]
a$ACME <- round(a$ACME, 3)
a$ACME_low <- round(a$ACME_low, 3)
a$ACME_high <- round(a$ACME_high, 3)
a$PropMediated <- round(a$PropMediated, 1)

a1 <- a[a$Exposure == "tea_freq_group",]
min(a1$ACME)
max(a1$ACME)
min(a1$ACME_low)
max(a1$ACME_high)
min(a1$PropMediated)
max(a1$PropMediated)

a1 <- a[a$Exposure != "tea_freq_group",]
min(a1$ACME)
max(a1$ACME)
min(a1$ACME_low)
max(a1$ACME_high)
min(a1$PropMediated)
max(a1$PropMediated)
####验证中介效应####
Mediator <- "s__Bacteroides_fragilis_clr_z"
Exposure <- "tea_freq_group"
Outcome  <- "ASMI_z"

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
    "`s__Bacteroides_fragilis_clr_z` ~ Treat_bin + ",
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
    "`ASMI_z` ~ Treat_bin + `s__Bacteroides_fragilis_clr_z` + ",
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

med_res <- mediate(
  model.m,
  model.y,
  treat = "Treat_bin",
  mediator = "s__Bacteroides_fragilis_clr_z",
  sims = 1000
)
med_res$d0
med_res$d0.ci
####十折交叉验证####
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
####菌群&骨骼肌
Lmer_results7_all_ASM_cv <- process_lmer(cv_res_micro$stability_sig$feature, 
                                      c("ASM_z","ASMI_z","grip_max_z"),
                                      data7_final,# %>% filter(tea_freq_group %in% c(0, 2)),
                                      c(Covariates_all_lmer,c("Fiber_soluble_mean")), #
                                      2)
Lmer_results7_all_ASM_cv_sig <- Lmer_results7_all_ASM_cv[Lmer_results7_all_ASM_cv$P_value < 0.05,]

####lmer-全调整
Lmer_results7_all2_cv <- process_lmer(c("tea_freq_group"),
                                   unique(Lmer_results7_all_ASM_sig$Predictor),
                                   data7_final,
                                   c(setdiff(Covariates_all_lmer,c("Protein_mean")),c("Fiber_soluble_mean")),
                                   2)
Lmer_results7_all2_cv <- Lmer_results7_all2_cv[Lmer_results7_all2_cv$Level == 2,]
Lmer_results7_all2_cv_sig <- Lmer_results7_all2_cv[Lmer_results7_all2_cv$P_value < 0.05,]

####计算指数
Micro_GMSI_all_tea <- calculate_GMSI(
  data = data7_final,
  outcomes_pos = Lmer_results7_all2_cv_sig[Lmer_results7_all2_cv_sig$Predictor == "tea_freq_group" & Lmer_results7_all2_cv_sig$Estimate > 0,]$Outcome,
  outcomes_neg = Lmer_results7_all2_cv_sig[Lmer_results7_all2_cv_sig$Predictor == "tea_freq_group" & Lmer_results7_all2_cv_sig$Estimate < 0,]$Outcome,
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
                P_value = P_value_formatted, 
                P_FDR = P_FDR_formatted) %>%
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
####功能&骨骼肌####
Lmer_results8_all_ASM <- process_lmer(Lmer_results8_all_sig$Outcome, 
                                      c("ASM_z","ASMI_z","grip_max_z"),
                                      data8_final, #%>% filter(tea_freq_group %in% c(0, 2)),
                                      c(Covariates_all_lmer,c("Fiber_soluble_mean","tea_freq_group")), #
                                      2)
Lmer_results8_all_ASM_sig <- Lmer_results8_all_ASM[Lmer_results8_all_ASM$P_value < 0.05,]
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
  p_col = "P_value"
)

ggsave("Func_ASM.png", plot = p, width = 13, height = 11, dpi = 500)
# pdf("Func_ASM.pdf", width = 12, height = 11)
# print(p)
# dev.off()



####lmer-全调整####
Lmer_results8_all2 <- process_lmer(c("tea_freq_group", "serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T"),
                                   unique(Lmer_results8_all_ASM_sig$Predictor),
                                   data8_final,
                                   c(setdiff(Covariates_all_lmer,c("Protein_mean")),c("Fiber_soluble_mean")),
                                   2)
Lmer_results8_all2 <- Lmer_results8_all2[Lmer_results8_all2$Level == 2,]
Lmer_results8_all2_sig <- Lmer_results8_all2[Lmer_results8_all2$P_value < 0.05,]
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
df_forest$Is_Significant <- ifelse(df_forest$P_value < 0.05, "Significant", "Non-significant")

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
  geom_text(aes(x = Star_X, label = Significance_Star),
            color = "#C94A29", size = 4.5, fontface = "bold", hjust = 0, vjust = 0.38) +
  
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
    
    legend.position = "none",
    plot.margin = ggplot2::margin(20, 20, 20, 20)
  ) +
  labs(
    x = "β (95% CI)",
    y = NULL
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
####可视化####
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggalluvial)
library(scales)

# 读取数据
mediation_results_all_Func <- read_excel(
  "D:/OneDrive/Papers/9 饮茶频率&骨骼肌/Output/mediation_results_all_Func1.xlsx",
  sheet = 1
)
a <- mediation_results_all_Func[mediation_results_all_Func$ACME_p < 0.05,]
mediation_results_all_Func$PropMediated <- mediation_results_all_Func$PropMediated*100
# 数据整理
plot_data_func <- mediation_results_all_Func %>%
  
  filter(
    ACME_p < 0.05,
    ACME > 0
  ) %>%
  
  dplyr::select(
    Exposure,
    Mediator,
    Outcome,
    ACME,
    PropMediated
  ) %>%
  
  mutate(
    
    Exposure = recode(
      Exposure,
      !!!predictor_rename
    ),
    
    Mediator = str_remove(
      Mediator,
      "_clr_z"
    ),
    
    Outcome = recode(
      Outcome,
      "ASM_z" = "ASM",
      "ASMI_z" = "ASMI",
      "grip_max_z" = "Handgrip strength"
    )
    
  )

p_mediation_func <- plot_mediation_alluvial(
  data = plot_data_func,
  title = "Mediation Effects: Tea Consumption/Biomarkers → Microbial Functional Pathways → Muscle",
  fill_low = "#CFE8F3",
  fill_high = "#0072B2" #
)

pdf(
  "Func_mediation.pdf",
  width = 11,
  height = 8
)
print(p_mediation_func)
dev.off()
#************正文写作
a <- mediation_results_all_Func[mediation_results_all_Func$ACME_p < 0.05,]
a$ACME <- round(a$ACME, 3)
a$ACME_low <- round(as.numeric(a$ACME_low), 3)
a$ACME_high <- round(a$ACME_high, 3)
a$PropMediated <- round(a$PropMediated,1)

a1 <- a[a$Exposure == "tea_freq_group",]
min(a1$ACME)
max(a1$ACME)
min(a1$ACME_low)
max(a1$ACME_high)
min(a1$PropMediated)
max(a1$PropMediated)

a1 <- a[a$Exposure != "tea_freq_group",]
min(a1$ACME)
max(a1$ACME)
min(a1$ACME_low)
max(a1$ACME_high)
min(a1$PropMediated)
max(a1$PropMediated)
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
                P_value = P_value_formatted, 
                P_FDR = P_FDR_formatted) %>%
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

####蛋白质&骨骼肌####
Lmer_results9_all_ASM <- process_lmer(Lmer_results9_all_FDR_sig$Outcome, 
                                      c("ASM_z","ASMI_z","grip_max_z"),
                                      data9_final, #%>% filter(tea_freq_group %in% c(0, 2)),
                                      c(Covariates_all_lmer,"tea_freq_group"), #
                                      2)
Lmer_results9_all_ASM_sig <- Lmer_results9_all_ASM[Lmer_results9_all_ASM$P_value < 0.05,]

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
  y_text_face = "plain"
)

ggsave("Pro_ASM.png", plot = p, width = 13, height = 6, dpi = 500)
# pdf("Func_ASM.pdf", width = 12, height = 11)
# print(p)
# dev.off()
####lmer-全调整####
Lmer_results9_all2 <- process_lmer(c("tea_freq_group", "serum_I_catechin_F0_T","serum_I_epicatechin_F0_T","serum_I_EGC_F0_T","serum_I_EGCG_F0_T","serum_I_ECG_F0_T"),
                                   unique(Lmer_results9_all_ASM_sig$Predictor),
                                   data9_final,
                                   c(setdiff(Covariates_all_lmer,c("Protein_mean")), "Phase_group"),
                                   2)
Lmer_results9_all2 <- Lmer_results9_all2[Lmer_results9_all2$Level == 2,]
Lmer_results9_all2_sig <- Lmer_results9_all2[Lmer_results9_all2$P_value < 0.05,]
table(Lmer_results9_all2_sig$Predictor)

Lmer_results9_all2$Outcome_name <- sub("_fill_log2_z", "", Lmer_results9_all2$Outcome)
Lmer_results9_all2$Outcome_name <- toupper(Lmer_results9_all2$Outcome_name)
Lmer_results9_all2 <- merge(Lmer_results9_all2, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Outcome_name", by.y = "Protein_names")
#*************正文写作
a <- Lmer_results9_all2_sig[Lmer_results9_all2_sig$Estimate >0,]
a$Estimate <- round(a$Estimate,3)
a$CI_low <- round(a$CI_low, 3)
a$CI_high <- round(a$CI_high, 3)

min(a$Estimate)
max(a$Estimate)
min(a$CI_low)
max(a$CI_high)

a <- Lmer_results9_all2_sig[Lmer_results9_all2_sig$Estimate <0,]
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
      P_value < 0.05 & Estimate > 0 ~ "Positive Sig",
      P_value < 0.05 & Estimate < 0 ~ "Negative Sig",
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
    
    covariates = Covariates_all_lmer,

    mediator_col = "Mediator",
    
    outcome_col = "Outcome_muscle",
    
    predictor_col = "Predictor",
    
    exposure_high = 2,
    exposure_low = 0,
    
    sims = 1000,
    
    output_file =
      "mediation_results_all_Pro1.csv"
  )
####可视化####
# 读取数据
mediation_results_all_pro <- read_excel(
  "D:/OneDrive/Papers/9 饮茶频率&骨骼肌/Output/mediation_results_all_Pro1.xlsx",
  sheet = 1
)
a <- mediation_results_all_pro[mediation_results_all_pro$ACME_p < 0.05,]

mediation_results_all_pro$PropMediated <- mediation_results_all_pro$PropMediated*100
# 数据整理
plot_data_pro <- mediation_results_all_pro %>%
  
  filter(
    ACME_p < 0.05,
    ACME > 0
  ) %>%
  
  dplyr::select(
    Exposure,
    Mediator,
    Outcome,
    ACME,
    PropMediated
  ) %>%
  
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
plot_data_pro$Mediator_name <- toupper(plot_data_pro$Mediator)
plot_data_pro <- merge(plot_data_pro, Protein_names[,c("Protein_names","hgnc_symbol")], by.x = "Mediator_name", by.y = "Protein_names")



p_mediation_pro <- plot_mediation_alluvial(
  data = plot_data_pro,
  title = "Mediation Effects: Tea Consumption/Biomarkers → Serum Proteins → Muscle",
  fill_low = "#E6D8F5",
  fill_high = "#6A51A3" ,
  mediator_col = "hgnc_symbol"
)

pdf(
  "Pro_mediation.pdf",
  width = 11,
  height = 9
)
print(p_mediation_pro)
dev.off()
#************正文写作
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