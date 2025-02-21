


# GSE3701
load("~/R/RNA/data/NASH/GSE37031/Data_GSE37031.RData")
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE37031, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(Affymetrix_GSE37031[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_Affymetrix_GSE37031 <- sweep(Affymetrix_GSE37031, 1, nc_mean, FUN = "/")



load("~/R/RNA/data/NASH/GSE63067/Data_GSE63067.RData")
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE63067, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(Affymetrix_GSE63067[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_Affymetrix_GSE63067 <- sweep(Affymetrix_GSE63067, 1, nc_mean, FUN = "/")


# GSE61260
load("~/R/RNA/data/NASH/GSE61260/Data_GSE61260.RData")
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE61260, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(Affymetrix_GSE61260[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_Affymetrix_GSE61260 <- sweep(Affymetrix_GSE61260, 1, nc_mean, FUN = "/")


# 
load("~/R/RNA/data/NASH/GSE66676/Data_GSE66676.RData")
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE66676, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(Affymetrix_GSE66676[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_Affymetrix_GSE66676 <- sweep(Affymetrix_GSE66676, 1, nc_mean, FUN = "/")


#
load("~/R/RNA/data/NASH/GSE89632/Data_GSE89632.RData")
# 1. 提取表达矩阵
Illumina_GSE89632_matrix <- Illumina_GSE89632$E
# 2. 赋值行名（基因名）
rownames(Illumina_GSE89632_matrix) <- Illumina_GSE89632$genes$Symbol
# 3. 赋值列名（样本名）并清理路径信息
colnames(Illumina_GSE89632_matrix) <- basename(Illumina_GSE89632$targets$IDATfile)
Illumina_GSE89632_matrix  <- as.data.frame(Illumina_GSE89632_matrix)
Illumina_GSE89632_matrix$gene <- rownames(Illumina_GSE89632_matrix)
Illumina_GSE89632_matrix$gene <- sub("\\..*", "", Illumina_GSE89632_matrix$gene)
Illumina_GSE89632_matrix <- aggregate(. ~ gene, data = Illumina_GSE89632_matrix, FUN = mean)
rownames(Illumina_GSE89632_matrix) <- Illumina_GSE89632_matrix$gene
Illumina_GSE89632_matrix$gene <- NULL
colnames(Illumina_GSE89632_matrix) <- sub("_.*", "", colnames(Illumina_GSE89632_matrix))
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE89632, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(Illumina_GSE89632_matrix[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_Illumina_GSE89632 <- sweep(Illumina_GSE89632_matrix, 1, nc_mean, FUN = "/")


# 
load("~/R/RNA/data/NASH/GSE48452/Data_GSE48452.RData")
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE48452, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(Affymetrix_GSE48452[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_Illumina_GSE48452 <- sweep(Affymetrix_GSE48452, 1, nc_mean, FUN = "/")

# 
load("~/R/RNA/data/NASH/GSE83452/Data_GSE83452.RData")
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE83452, condition == "NC_baseline"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(Affymetrix_GSE83452[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_Affymetrix_GSE83452 <- sweep(Affymetrix_GSE83452, 1, nc_mean, FUN = "/")



# GSE106737
load("~/R/RNA/data/NASH/GSE106737/Data_GSE106737.RData")
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE106737, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(Affymetrix_GSE106737[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_Illumina_GSE106737 <- sweep(Affymetrix_GSE106737, 1, nc_mean, FUN = "/")



# 
load("~/R/RNA/data/NASH/GSE105127/Data_GSE105127.RData")
vsd_GSE105127_matrix <- assay(vsd_GSE105127)
vsd_GSE105127 <- as.data.frame(vsd_GSE105127_matrix)
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE105127, condition == "CV_NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(vsd_GSE105127[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_vsd_GSE105127 <- sweep(vsd_GSE105127, 1, nc_mean, FUN = "/")


# 
load("~/R/RNA/data/NASH/GSE115193/Data_GSE115193.RData")
vsd_GSE115193_matrix <- assay(vsd_GSE115193)
vsd_GSE115193 <- as.data.frame(vsd_GSE115193_matrix)
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE115193, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(vsd_GSE115193[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_vsd_GSE115193 <- sweep(vsd_GSE115193, 1, nc_mean, FUN = "/")


# 
load("~/R/RNA/data/NASH/GSE126848/Data_GSE126848.RData")
vsd_GSE126848_matrix <- assay(vsd_GSE126848)
vsd_GSE126848 <- as.data.frame(vsd_GSE126848_matrix)
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE126848, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(vsd_GSE126848[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_vsd_GSE126848 <- sweep(vsd_GSE126848, 1, nc_mean, FUN = "/")


# GSE134146
load("~/R/RNA/data/NASH/GSE134146/Data_GSE134146.RData")
Arraystar_GSE134146 <- log2(Arraystar_GSE134146 + 1)
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE134146, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(Arraystar_GSE134146[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_Arraystar_GSE134146 <- sweep(Arraystar_GSE134146, 1, nc_mean, FUN = "/")


# 
load("~/R/RNA/data/NASH/GSE135251/Data_GSE135251.RData")
vsd_level_GSE135251_matrix <- assay(vsd_level_GSE135251)
vsd_level_GSE135251 <- as.data.frame(vsd_level_GSE135251_matrix)
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE135251, level == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(vsd_level_GSE135251[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_vsd_level_GSE135251 <- sweep(vsd_level_GSE135251, 1, nc_mean, FUN = "/")



# GSE159676
load("~/R/RNA/data/NASH/GSE159676/Data_GSE159676.RData")
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE159676, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(Affymetrix__GSE159676[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_Affymetrix_GSE159676 <- sweep(Affymetrix__GSE159676, 1, nc_mean, FUN = "/")



# GSE147304
load("~/R/RNA/data/NASH/GSE147304/Data_GSE147304.RData")
vsd_GSE147304_matrix <- assay(vsd_GSE147304)
vsd_GSE147304 <- as.data.frame(vsd_GSE147304_matrix)
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE147304, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(vsd_GSE147304[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_vsd_GSE147304 <- sweep(vsd_GSE147304, 1, nc_mean, FUN = "/")



# 
load("~/R/RNA/data/NASH/GSE164760/Data_GSE164760.RData")
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE164760, level == "Healthy"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(affy_GSE164760[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_affy_GSE164760 <- sweep(affy_GSE164760, 1, nc_mean, FUN = "/")



# GSE173735
load("~/R/RNA/data/NASH/GSE173735/Data_GSE173735.RData")
vsd_GSE173735 <- as.data.frame(vsd_GSE173735)
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE173735, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(vsd_GSE173735[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_vsd_GSE173735 <- sweep(vsd_GSE173735, 1, nc_mean, FUN = "/")



# 
load("~/R/RNA/data/NASH/GSE185051/Data_GSE185051.RData")
vsd_level_GSE185051_matrix <- assay(vsd_level_GSE185051)
vsd_level_GSE185051 <- as.data.frame(vsd_level_GSE185051_matrix)
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE185051, level == "Normal"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(vsd_level_GSE185051[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_vsd_level_GSE185051 <- sweep(vsd_level_GSE185051, 1, nc_mean, FUN = "/")



# 加载数据
load("~/R/RNA/data/NASH/GSE207310/Data_GSE207310.RData")
vsd_level_GSE207310_matrix <- assay(vsd_level_GSE207310)
vsd_level_GSE207310 <- as.data.frame(vsd_level_GSE207310_matrix)
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE207310, level == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(vsd_level_GSE207310[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_vsd_level_GSE207310 <- sweep(vsd_level_GSE207310, 1, nc_mean, FUN = "/")



# 加载数据
load("~/R/RNA/data/NASH/GSE260666/Data_GSE260666.RData")
vsd_GSE260666_matrix <- assay(vsd_GSE260666)
vsd_GSE260666 <- as.data.frame(vsd_GSE260666_matrix)
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE260666, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(vsd_GSE260666[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_vsd_GSE260666 <- sweep(vsd_GSE260666, 1, nc_mean, FUN = "/")



#
load("~/R/RNA/data/NASH/GSE274114/Data_GSE274114.RData")
vsd_GSE274114_matrix <- assay(vsd_GSE274114)
vsd_GSE274114 <- as.data.frame(vsd_GSE274114_matrix)
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE274114, condition == "NC"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(vsd_GSE274114[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_vsd_GSE274114 <- sweep(vsd_GSE274114, 1, nc_mean, FUN = "/")



#
load("~/R/RNA/data/NASH/GSE163211/Data_GSE163211.RData")
vsd_GSE163211_matrix <- assay(vsd_GSE163211)
vsd_GSE163211 <- as.data.frame(vsd_GSE163211_matrix)
# 1. 确定 "NC" 组的样本
nc_samples <- rownames(subset(colData_GSE163211, level == "Normal"))
# 2. 计算 NC 组的基因平均表达值
nc_mean <- rowMeans(vsd_GSE163211[, nc_samples], na.rm = TRUE)
# 3. 计算比值归一化（所有样本的表达值 / NC 组均值）
normalized_vsd_GSE163211 <- sweep(vsd_GSE163211, 1, nc_mean, FUN = "/")

















# 1. 创建归一化数据的列表
normalized_data_list <- list(
  normalized_Affymetrix_GSE37031,
  normalized_Affymetrix_GSE63067,
  normalized_Affymetrix_GSE61260,
  normalized_Affymetrix_GSE66676,
  normalized_Illumina_GSE89632,
  normalized_Illumina_GSE48452,
  normalized_Affymetrix_GSE83452,
  normalized_Illumina_GSE106737,
  normalized_vsd_GSE105127,
  normalized_vsd_GSE115193,
  normalized_vsd_GSE126848,
  normalized_Arraystar_GSE134146,
  normalized_vsd_level_GSE135251,
  normalized_Affymetrix_GSE159676,
  normalized_vsd_GSE147304,
  normalized_affy_GSE164760,
  normalized_vsd_GSE173735,
  normalized_vsd_level_GSE185051,
  normalized_vsd_level_GSE207310,
  normalized_vsd_GSE260666,
  normalized_vsd_GSE274114,
  normalized_vsd_GSE163211
)

# 2. 为每个数据集添加基因列并转换为数据框
normalized_data_list <- lapply(normalized_data_list, function(df) {
  df <- as.data.frame(df)
  df$gene <- rownames(df)  # 将行名作为基因列
  return(df)
})

# 3. 使用 dplyr::reduce() 按基因合并数据集
combined_data <- Reduce(function(x, y) {
  left_join(x, y, by = "gene", suffix = c("_x", "_y"))
}, normalized_data_list)

# 4. 设置基因名为行名，删除gene列
rownames(combined_data) <- combined_data$gene
combined_data$gene <- NULL  # 删除gene列






# 1. 创建包含所有数据集的 colData 列表
condition_list <- list(
  colData_GSE37031,
  colData_GSE63067,
  colData_GSE61260,
  colData_GSE66676,
  colData_GSE89632,
  colData_GSE48452,
  colData_GSE83452,
  colData_GSE106737,
  colData_GSE105127,
  colData_GSE115193,
  colData_GSE126848,
  colData_GSE134146,
  colData_GSE135251,
  colData_GSE159676,
  colData_GSE147304,
  colData_GSE164760,
  colData_GSE173735,
  colData_GSE185051,
  colData_GSE207310,
  colData_GSE260666,
  colData_GSE274114,
  colData_GSE163211
)

# 2. 提取每个数据集的 'condition' 或 'level' 列，统一为 'condition' 列
condition_list <- lapply(condition_list, function(df) {
  if ("condition" %in% colnames(df)) {
    return(df[, "condition", drop = FALSE])  # 只保留 condition 列
  } else if ("level" %in% colnames(df)) {
    return(df[, "level", drop = FALSE])  # 只保留 level 列
  }
})

# 3. 确保每个数据集都是数据框并设置列名为 'condition'
condition_list <- lapply(condition_list, function(df) {
  if (is.vector(df)) {
    df <- as.data.frame(df)  # 将向量转换为数据框
  }
  colnames(df) <- "condition"  # 统一列名为 "condition"
  return(df)
})

# 5. 将所有数据按行名（样本名）纵向叠加
combined_colData <- do.call(rbind, condition_list)

# 使用 gsub 进行字符串替换
combined_colData$condition <- gsub("Definite_NASH|Borderline_NASH|NASH_baseline|PP_EARLY_NASH|CV_EARLY_NASH|IZ_EARLY_NASH|NASH_F1_F4|NASH_F0|NASH_HCC_notumor|Borderline", "NASH", combined_colData$condition)
combined_colData$condition <- gsub("CV_NC|PP_NC|IZ_NC|NC_baseline|Healthy|Normal", "NC", combined_colData$condition)
combined_colData$condition <- gsub("IZ_STEA|CV_STEA|PP_STEA|Steatosis", "Obese", combined_colData$condition)
combined_colData$condition <- gsub("CV_HO|IZ_HO|PP_HO", "NAFLD", combined_colData$condition)
# 使用 ifelse 进行条件替换
combined_colData$condition <- ifelse(combined_colData$condition == "NAFL", 
                                     "NAFLD", 
                                     combined_colData$condition)

combined_colData$condition <- factor(combined_colData$condition)



save(combined_data,combined_colData, file = "~/R/RNA/data/NASH/Summary/Data_Summary_Human.RData")


load("~/R/RNA/data/NASH/Summary/Data_Summary_Human.RData")

# 去除极值
analyze_gene_expression_combined <- function(GeneX, combined_data, combined_colData) {
  # Step 1: 检查基因是否存在于表达矩阵
  if (!(GeneX %in% rownames(combined_data))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 获取目标基因的表达值
  gene_expression <- combined_data[GeneX, ]
  
  # Step 3: 将基因表达值和条件信息结合，确保是数据框格式
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),  # 确保基因表达值为数值类型
    Condition = combined_colData$condition
  )
  
  # Step 4: 可选：根据需要只保留特定的条件组（例如“NC”，“Obese”，“NAFLD”，“NASH”）
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "Obese", "NAFLD", "NASH"), ]
  
  # Step 5: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "Obese", "NAFLD", "NASH"))
  
  # Step 6: 按因子水平排序数据
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 7: 去除每组5%的极端值
  gene_df_clean <- gene_df %>%
    group_by(Condition) %>%
    mutate(
      lower_quantile = quantile(Expression, 0.01, na.rm = TRUE),  # 计算下2.5%分位数
      upper_quantile = quantile(Expression, 0.99, na.rm = TRUE)   # 计算上2.5%分位数
    ) %>%
    filter(Expression >= lower_quantile & Expression <= upper_quantile)  # 保留中间的90%
  
  # Step 8: 统计分析（NC vs. Obese, NC vs. NAFLD, NC vs. NASH）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "Obese"),
                                                     c("NC", "NAFLD"),
                                                     c("NC", "NASH"),
                                                     c("Obese", "NAFLD"),
                                                     c("Obese", "NASH"),
                                                     c("NAFLD", "NASH")
                                                     ), 
                                  method = "t.test", label = "p.signif")
  
  # Step 9: 输出每个组的非空样本数
  non_na_counts <- table(gene_df_clean$Condition)
  
  # 创建图例标签
  legend_labels <- paste(names(non_na_counts), " (n=", non_na_counts, ")", sep = "")
  
  # Step 10: 自定义颜色（动态生成）
  unique_conditions <- levels(gene_df_clean$Condition)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  
  # Step 11: 绘制箱线图，并在图例中显示每个组的非空样本数
  plot <- ggplot(gene_df_clean, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    geom_jitter(aes(color = Condition), width = 0.2, size = 0.5, alpha = 0.6) + 
    scale_fill_manual(values = custom_colors_fill, labels = legend_labels) +  
    scale_color_manual(values = custom_colors_color, labels = legend_labels) +  
    stat_test +  
    labs(
      title = paste("Expression of", GeneX),
      x = NULL,
      y = "Expression Level",
      fill = "Condition",   # 图例标题
      color = "Condition"   # 图例标题
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "right"
    )
  
  # 返回绘图
  return(plot)
}

# 调用函数进行绘制
analyze_gene_expression_combined("SREBF1", combined_data, combined_colData)

