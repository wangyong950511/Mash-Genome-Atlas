# library(sva)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(dplyr)
library(data.table)
library(AnnotationDbi)
library(babelgene)
library(patchwork)
library(myplotColors)
library(cowplot)
library(myplotColors)


# 函数集合 ----------

####总结数据
load("~/R/RNA/data/NASH/Summary/Data_Summary_Human.RData")


###Fibrosis
# 
load("~/R/RNA/data/NASH/GSE130970/Data_GSE130970.RData")
analyze_gene_expression_GSE130970 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nafld' 或 'fibrosis'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nafld" = c("nafld_0", "nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_6"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # **确保所有组都显示，即使某些组样本数 <2**
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # **确定可以进行统计检验的组（样本数 >= 2）**
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE130970"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE89632/Data_GSE89632.RData")
analyze_gene_expression_GSE89632 <- function(GeneX, norm_data, colData, group) {
  # Step 1: 检查并匹配目标基因在 norm_data$genes$Symbol 中
  if (!("genes" %in% names(norm_data))) {
    stop("norm_data 对象中没有 'genes' 信息。")
  }
  if (!("Symbol" %in% colnames(norm_data$genes))) {
    stop("norm_data$genes 中不包含 'Symbol' 列。")
  }
  idx <- which(norm_data$genes$Symbol == GeneX)
  if (length(idx) == 0) {
    stop(paste("目标基因", GeneX, "在 norm_data$genes$Symbol 中未找到。"))
  }
  # 如果有多个匹配，则取第一个（也可以根据需要修改策略）
  gene_expression <- norm_data$E[idx[1], ]
  
  # Step 2: 检查分组变量是否存在于 colData 中
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'condition', 'fibrosis' 或 'nafld'。")
  }
  
  # Step 3: 构建绘图数据框，提取目标基因表达和对应分组信息
  gene_df <- data.frame(Expression = gene_expression,
                        Condition  = colData[[group]])
  
  # 筛除 Condition 列中为 NA 或字符串 "NA" 的样本
  gene_df <- subset(gene_df, !is.na(Condition) & Condition != "NA")
  
  # Step 4: 根据不同分组设定固定的分组顺序（所有组都显示）
  group_order <- switch(group,
                        "condition" = c("NC", "Obese", "NASH"),
                        "fibrosis"  = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4", "NA"),
                        "nafld"     = c("nafld_0", "nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_6", "nafld_8", "NA"),
                        unique(gene_df$Condition)  # 默认按原顺序
  )
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算各组样本数量，并确定可以进行统计检验的组（样本数 >= 2）
  group_counts <- table(gene_df$Condition)
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: 自动生成 comparisons 列表，仅比较有效组中第一个分组与其他分组
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()
    warning("没有足够的样本进行统计检验，仅显示箱线图。")
  }
  
  # Step 7: 设置颜色映射
  # 填充色设为白色，边框颜色采用 myplotColors 函数生成的调色板
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color)
  
  # Step 9: 添加统计检验（仅在有效组数 >= 2 时）
  if (length(comparisons_list) > 0) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list,
                                      method = "t.test",
                                      label = "p.signif")
  }
  
  # Step 10: 添加标题和图形主题
  plot <- plot + labs(
    title = paste("GSE89632"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE48452/Data_GSE48452.RData")
analyze_gene_expression_GSE48452 <- function(GeneX, expr_data, meta_data, group) {
  # Step 1: 确保基因存在
  if (!(GeneX %in% rownames(expr_data))) {
    stop(paste("目标基因", GeneX, "不存在于表达矩阵中，请检查数据。"))
  }
  
  # Step 2: 提取目标基因的表达数据
  gene_expression <- expr_data[GeneX, ]
  
  # Step 3: 确保分组变量存在
  if (!(group %in% colnames(meta_data))) {
    stop(paste("分组变量", group, "不存在，请检查。"))
  }
  
  # Step 4: 创建数据框
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Condition  = meta_data[[group]]
  )
  
  # Step 5: 移除 NA 值
  gene_df <- subset(gene_df, !is.na(Condition) & Condition != "fibrosis_NA")
  
  # Step 6: 设定分组顺序
  group_order <- switch(group,
                        "condition" = c("NC", "Obese", "Steatosis", "NASH"),
                        "fibrosis"  = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        "nas"       = c("nas_0", "nas_1", "nas_2", "nas_3", "nas_4", "nas_5", "nas_6", "nas_7"),
                        unique(gene_df$Condition))  # 默认按数据顺序
  
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 7: 统计学比较（基准为分组的第一个）
  group_counts <- table(gene_df$Condition)
  valid_groups <- names(group_counts[group_counts >= 2])
  
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()
    warning("没有足够的样本进行统计检验，仅显示箱线图。")
  }
  
  # Step 8: 设置颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 9: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(width = 0.2, size = 1, alpha = 0.6)+
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color)
  
  # Step 10: 添加统计检验（仅在有效组数 >= 2 时）
  if (length(comparisons_list) > 0) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list,
                                      method = "t.test",
                                      label = "p.signif")
  }
  
  # Step 11: 添加标题和图形主题
  plot <- plot + labs(
    title = paste("GSE48452"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE135251/Data_GSE135251.RData")
analyze_gene_expression_GSE135251 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nas', 'fibrosis' 或 'steatohep'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  
  # Step 5: 设定分组顺序
  group_order <- switch(group,
                        "nas" = c("nas_0", "nas_1", "nas_2", "nas_3", "nas_4", "nas_5", "nas_6", "nas_7", "nas_8"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        "level" = c("NC", "NAFL", "NASH_F0_F1", "NASH_F2", "NASH_F3", "NASH_F4"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # 确保分组顺序只包含数据中的分组
  group_order <- intersect(group_order, unique(gene_df$Condition))
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 6: **自动生成 comparisons**
  condition_levels <- levels(gene_df$Condition)
  comparisons_list <- lapply(condition_levels[-1], function(x) c(condition_levels[1], x))
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(unique(gene_df$Condition))), unique(gene_df$Condition))
  custom_colors_color <- setNames(getplotColors(length(unique(gene_df$Condition))), unique(gene_df$Condition))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) + 
    stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif") +  # 统计显著性
    labs(
      title = paste("GSE135251"),
      x = group,
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE162694/Data_GSE162694.RData")
analyze_gene_expression_GSE162694 <- function(GeneX, vsd, colData) {
  # 检测基因是否包含小写字母（判定为小鼠基因）
  if (grepl("[a-z]", GeneX)) {
    # 使用 babelgene::orthologs() 查询同源人类基因
    result <- orthologs(genes = GeneX, species = "mouse", human = FALSE)
    # 提取人类基因名
    GeneX <- result$human_symbol
  } 
  # Step 1: 检查目标基因是否存在于表达矩阵
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因未找到，请检查基因名是否正确或是否存在于表达矩阵中。")
  }
  
  # Step 2: 获取目标基因的表达值
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition
  )
  
  # Step 3: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "0", "1", "2", "3", "4"))
  
  # 按 Condition 因子水平排序数据框
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 4: 绘制箱线图
  custom_colors_fill <- c("NC" = "white", "0" = "white", "1" = "white", "2" = "white", "3" = "white", "4" = "white")
  custom_colors_color <- c("NC" = "#98DF8A", "0" = "#FFB347", "1" = "#FF8C42", "2" = "#FF6F31", "3"  = "#D62728", "4" = "#A51220")
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  # 边框颜色根据 Condition 设置
    scale_fill_manual(values = custom_colors_fill) +  # 自定义填充颜色
    scale_color_manual(values = custom_colors_color) + 
    stat_compare_means(comparisons = list(c("NC", "0"), 
                                          c("NC", "1"), 
                                          c("NC", "2"), 
                                          c("NC", "3"), 
                                          c("NC", "4")), 
                       method = "t.test",
                       label = "p.signif") + # 两两比较，显示显著性标记
    labs(
      title = paste("GSE162694"),
      x = "Condition",
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE174478/Data_GSE174478.RData")
analyze_gene_expression_GSE174478 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nafld' 或 'fibrosis'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nafld" = c("nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_6", "nafld_7", "nafld_8"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # **确保所有组都显示，即使某些组样本数 <2**
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # **确定可以进行统计检验的组（样本数 >= 2）**
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE174478"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE185051/Data_GSE185051.RData")
analyze_gene_expression_GSE185051 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nas', 'fibrosis' 或 'level'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nas" = c("nas_0", "nas_1", "nas_2", "nas_3", "nas_4", "nas_5", "nas_6", "nas_7"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3"),
                        "level" = c("Normal", "Steatosis", "Borderline", "NASH"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # 确保所有组都显示，即使某些组样本数 <2
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # 确定可以进行统计检验的组（样本数 >= 2）
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE185051"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE193066/Data_GSE193066.RData")
analyze_gene_expression_GSE193066 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nafld' 或 'fibrosis'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nafld" = c("nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_6", "nafld_7", "nafld_8"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # **确保所有组都显示，即使某些组样本数 <2**
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # **确定可以进行统计检验的组（样本数 >= 2）**
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE193066"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE193080/Data_GSE193080.RData")
analyze_gene_expression_GSE193080 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nafld' 或 'fibrosis'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nafld" = c("nafld_0", "nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_6", "nafld_7", "nafld_8"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # **确保所有组都显示，即使某些组样本数 <2**
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # **确定可以进行统计检验的组（样本数 >= 2）**
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE193080"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE207310/Data_GSE207310.RData")
analyze_gene_expression_GSE207310 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nafld', 'fibrosis' 或 'level'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nafld" = c("nafld_0", "nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_7", "nafld_8"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_1A", "fibrosis_1B", "fibrosis_2", "fibrosis_3"),
                        "level" = c("NC", "NAFL", "NASH"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # **确保所有组都显示，即使某些组样本数 <2**
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # **确定可以进行统计检验的组（样本数 >= 2）**
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE207310"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE246221/Data_GSE246221.RData")
analyze_gene_expression_GSE246221 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nafld' 或 'fibrosis'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nafld" = c("NC", "NAFLD_2_3", "NAFLD_5_6", "NAFLD_4_7"),
                        "fibrosis" = c("NC", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # **确保所有组都显示，即使某些组样本数 <2**
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # **确定可以进行统计检验的组（样本数 >= 2）**
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE246221"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE225740/Data_GSE225740.RData")
analyze_gene_expression_GSE225740 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nas', 'fibrosis' 或 'steatohep'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 处理 NA 值 **(去除 NA 组)**
  gene_df <- gene_df[!is.na(gene_df$Condition) & gene_df$Condition != "fibrosis_NA" & gene_df$Condition != "steatohep_NA", ]
  
  # Step 5: 设定分组顺序
  group_order <- switch(group,
                        "nas" = c("nas_0", "nas_1", "nas_2", "nas_3", "nas_4", "nas_5", "nas_6", "nas_7"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        "steatohep" = c("steatohep_0", "steatohep_1", "steatohep_2"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # 确保分组顺序只包含数据中的分组
  group_order <- intersect(group_order, unique(gene_df$Condition))
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 6: **自动生成 comparisons**
  condition_levels <- levels(gene_df$Condition)
  comparisons_list <- lapply(condition_levels[-1], function(x) c(condition_levels[1], x))
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(unique(gene_df$Condition))), unique(gene_df$Condition))
  custom_colors_color <- setNames(getplotColors(length(unique(gene_df$Condition))), unique(gene_df$Condition))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) + 
    stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif") +  # 统计显著性
    labs(
      title = paste("GSE225740"),
      x = group,
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE240729/Data_GSE240729.RData")
analyze_gene_expression_GSE240729 <- function(GeneX, vsd, colData) {
  # Step 1: 检查基因是否存在于表达矩阵
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 获取目标基因的表达值
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition
  )
  
  # Step 3: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("fibrosisscore: F0", "fibrosisscore: F1", "fibrosisscore: F2", "fibrosisscore: F3", "fibrosisscore: F4"))
  
  # 按 Condition 因子水平排序数据框
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 4: 绘制箱线图
  custom_colors_fill <- c("fibrosisscore: F0" = "white", "fibrosisscore: F1" = "white", "fibrosisscore: F2" = "white", "fibrosisscore: F3" = "white", "fibrosisscore: F4" = "white")
  custom_colors_color <- c("fibrosisscore: F0" = "#2ca02c", "fibrosisscore: F1" = "#FFC07C", "fibrosisscore: F2" = "#ff7f0e", "fibrosisscore: F3" = "#d62728", "fibrosisscore: F4" = "#d62728")
  
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  # 边框颜色根据 Condition 设置
    scale_fill_manual(values = custom_colors_fill) +  # 自定义填充颜色
    scale_color_manual(values = custom_colors_color) +  # 自定义边框颜色
    stat_compare_means(comparisons = list(c("fibrosisscore: F0", "fibrosisscore: F1"), 
                                          c("fibrosisscore: F0", "fibrosisscore: F2"), 
                                          c("fibrosisscore: F0", "fibrosisscore: F3"), 
                                          c("fibrosisscore: F0", "fibrosisscore: F4")),
                       method = "t.test",
                       label = "p.signif") + # 两两比较，显示显著性标记
    labs(
      title = paste("GSE240729"),
      x = "Condition",
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


###NAFLD
# GSE130970
load("~/R/RNA/data/NASH/GSE130970/Data_GSE130970.RData")
analyze_gene_expression_GSE130970 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nafld' 或 'fibrosis'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nafld" = c("nafld_0", "nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_6"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # **确保所有组都显示，即使某些组样本数 <2**
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # **确定可以进行统计检验的组（样本数 >= 2）**
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE130970"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE89632/Data_GSE89632.RData")
analyze_gene_expression_GSE89632 <- function(GeneX, norm_data, colData, group) {
  # Step 1: 检查并匹配目标基因在 norm_data$genes$Symbol 中
  if (!("genes" %in% names(norm_data))) {
    stop("norm_data 对象中没有 'genes' 信息。")
  }
  if (!("Symbol" %in% colnames(norm_data$genes))) {
    stop("norm_data$genes 中不包含 'Symbol' 列。")
  }
  idx <- which(norm_data$genes$Symbol == GeneX)
  if (length(idx) == 0) {
    stop(paste("目标基因", GeneX, "在 norm_data$genes$Symbol 中未找到。"))
  }
  # 如果有多个匹配，则取第一个（也可以根据需要修改策略）
  gene_expression <- norm_data$E[idx[1], ]
  
  # Step 2: 检查分组变量是否存在于 colData 中
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'condition', 'fibrosis' 或 'nafld'。")
  }
  
  # Step 3: 构建绘图数据框，提取目标基因表达和对应分组信息
  gene_df <- data.frame(Expression = gene_expression,
                        Condition  = colData[[group]])
  
  # 筛除 Condition 列中为 NA 或字符串 "NA" 的样本
  gene_df <- subset(gene_df, !is.na(Condition) & Condition != "NA")
  
  # Step 4: 根据不同分组设定固定的分组顺序（所有组都显示）
  group_order <- switch(group,
                        "condition" = c("NC", "Obese", "NASH"),
                        "fibrosis"  = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4", "NA"),
                        "nafld"     = c("nafld_0", "nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_6", "nafld_8", "NA"),
                        unique(gene_df$Condition)  # 默认按原顺序
  )
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算各组样本数量，并确定可以进行统计检验的组（样本数 >= 2）
  group_counts <- table(gene_df$Condition)
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: 自动生成 comparisons 列表，仅比较有效组中第一个分组与其他分组
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()
    warning("没有足够的样本进行统计检验，仅显示箱线图。")
  }
  
  # Step 7: 设置颜色映射
  # 填充色设为白色，边框颜色采用 myplotColors 函数生成的调色板
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color)
  
  # Step 9: 添加统计检验（仅在有效组数 >= 2 时）
  if (length(comparisons_list) > 0) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list,
                                      method = "t.test",
                                      label = "p.signif")
  }
  
  # Step 10: 添加标题和图形主题
  plot <- plot + labs(
    title = paste("GSE89632"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
# 
load("~/R/RNA/data/NASH/GSE48452/Data_GSE48452.RData")
analyze_gene_expression_GSE48452 <- function(GeneX, expr_data, meta_data, group) {
  # Step 1: 确保基因存在
  if (!(GeneX %in% rownames(expr_data))) {
    stop(paste("目标基因", GeneX, "不存在于表达矩阵中，请检查数据。"))
  }
  
  # Step 2: 提取目标基因的表达数据
  gene_expression <- expr_data[GeneX, ]
  
  # Step 3: 确保分组变量存在
  if (!(group %in% colnames(meta_data))) {
    stop(paste("分组变量", group, "不存在，请检查。"))
  }
  
  # Step 4: 创建数据框
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Condition  = meta_data[[group]]
  )
  
  # Step 5: 移除 NA 值
  gene_df <- subset(gene_df, !is.na(Condition) & Condition != "fibrosis_NA")
  
  # Step 6: 设定分组顺序
  group_order <- switch(group,
                        "condition" = c("NC", "Obese", "Steatosis", "NASH"),
                        "fibrosis"  = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        "nas"       = c("nas_0", "nas_1", "nas_2", "nas_3", "nas_4", "nas_5", "nas_6", "nas_7"),
                        unique(gene_df$Condition))  # 默认按数据顺序
  
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 7: 统计学比较（基准为分组的第一个）
  group_counts <- table(gene_df$Condition)
  valid_groups <- names(group_counts[group_counts >= 2])
  
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()
    warning("没有足够的样本进行统计检验，仅显示箱线图。")
  }
  
  # Step 8: 设置颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 9: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(width = 0.2, size = 1, alpha = 0.6)+
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color)
  
  # Step 10: 添加统计检验（仅在有效组数 >= 2 时）
  if (length(comparisons_list) > 0) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list,
                                      method = "t.test",
                                      label = "p.signif")
  }
  
  # Step 11: 添加标题和图形主题
  plot <- plot + labs(
    title = paste("GSE89632"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
# 
load("~/R/RNA/data/NASH/GSE135251/Data_GSE135251.RData")
analyze_gene_expression_GSE135251 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nas', 'fibrosis' 或 'steatohep'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  
  # Step 5: 设定分组顺序
  group_order <- switch(group,
                        "nas" = c("nas_0", "nas_1", "nas_2", "nas_3", "nas_4", "nas_5", "nas_6", "nas_7", "nas_8"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        "level" = c("NC", "NAFL", "NASH_F0_F1", "NASH_F2", "NASH_F3", "NASH_F4"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # 确保分组顺序只包含数据中的分组
  group_order <- intersect(group_order, unique(gene_df$Condition))
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 6: **自动生成 comparisons**
  condition_levels <- levels(gene_df$Condition)
  comparisons_list <- lapply(condition_levels[-1], function(x) c(condition_levels[1], x))
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(unique(gene_df$Condition))), unique(gene_df$Condition))
  custom_colors_color <- setNames(getplotColors(length(unique(gene_df$Condition))), unique(gene_df$Condition))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) + 
    stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif") +  # 统计显著性
    labs(
      title = paste("GSE135251"),
      x = group,
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
# 
load("~/R/RNA/data/NASH/GSE174478/Data_GSE174478.RData")
analyze_gene_expression_GSE174478 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nafld' 或 'fibrosis'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nafld" = c("nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_6", "nafld_7", "nafld_8"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # **确保所有组都显示，即使某些组样本数 <2**
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # **确定可以进行统计检验的组（样本数 >= 2）**
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE174478"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
# 
load("~/R/RNA/data/NASH/GSE185051/Data_GSE185051.RData")
analyze_gene_expression_GSE185051 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nas', 'fibrosis' 或 'level'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nas" = c("nas_0", "nas_1", "nas_2", "nas_3", "nas_4", "nas_5", "nas_6", "nas_7"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3"),
                        "level" = c("Normal", "Steatosis", "Borderline", "NASH"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # 确保所有组都显示，即使某些组样本数 <2
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # 确定可以进行统计检验的组（样本数 >= 2）
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE185051"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
# 
load("~/R/RNA/data/NASH/GSE193066/Data_GSE193066.RData")
analyze_gene_expression_GSE193066 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nafld' 或 'fibrosis'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nafld" = c("nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_6", "nafld_7", "nafld_8"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # **确保所有组都显示，即使某些组样本数 <2**
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # **确定可以进行统计检验的组（样本数 >= 2）**
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE193066"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
# 
load("~/R/RNA/data/NASH/GSE193080/Data_GSE193080.RData")
analyze_gene_expression_GSE193080 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nafld' 或 'fibrosis'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nafld" = c("nafld_0", "nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_6", "nafld_7", "nafld_8"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # **确保所有组都显示，即使某些组样本数 <2**
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # **确定可以进行统计检验的组（样本数 >= 2）**
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE193080"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
# 
load("~/R/RNA/data/NASH/GSE207310/Data_GSE207310.RData")
analyze_gene_expression_GSE207310 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nafld', 'fibrosis' 或 'level'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nafld" = c("nafld_0", "nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_7", "nafld_8"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_1A", "fibrosis_1B", "fibrosis_2", "fibrosis_3"),
                        "level" = c("NC", "NAFL", "NASH"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # **确保所有组都显示，即使某些组样本数 <2**
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # **确定可以进行统计检验的组（样本数 >= 2）**
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE207310"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE246221/Data_GSE246221.RData")
analyze_gene_expression_GSE246221 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nafld' 或 'fibrosis'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nafld" = c("NC", "NAFLD_2_3", "NAFLD_5_6", "NAFLD_4_7"),
                        "fibrosis" = c("NC", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # **确保所有组都显示，即使某些组样本数 <2**
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # **确定可以进行统计检验的组（样本数 >= 2）**
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE246221"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}
#
load("~/R/RNA/data/NASH/GSE225740/Data_GSE225740.RData")
analyze_gene_expression_GSE225740 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nas', 'fibrosis' 或 'steatohep'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 处理 NA 值 **(去除 NA 组)**
  gene_df <- gene_df[!is.na(gene_df$Condition) & gene_df$Condition != "fibrosis_NA" & gene_df$Condition != "steatohep_NA", ]
  
  # Step 5: 设定分组顺序
  group_order <- switch(group,
                        "nas" = c("nas_0", "nas_1", "nas_2", "nas_3", "nas_4", "nas_5", "nas_6", "nas_7"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        "steatohep" = c("steatohep_0", "steatohep_1", "steatohep_2"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # 确保分组顺序只包含数据中的分组
  group_order <- intersect(group_order, unique(gene_df$Condition))
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 6: **自动生成 comparisons**
  condition_levels <- levels(gene_df$Condition)
  comparisons_list <- lapply(condition_levels[-1], function(x) c(condition_levels[1], x))
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(unique(gene_df$Condition))), unique(gene_df$Condition))
  custom_colors_color <- setNames(getplotColors(length(unique(gene_df$Condition))), unique(gene_df$Condition))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) + 
    stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif") +  # 统计显著性
    labs(
      title = paste("GSE225740"),
      x = group,
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}



####NCvsNASH
# GSE24807
load("~/R/RNA/data/NASH/GSE24807/Data_GSE24807.RData")
analyze_gene_expression_GSE24807 <- function(GeneX, exprs_matrix, colData) {
  if (!(GeneX %in% rownames(exprs_matrix))) {
    stop(paste("目标基因", GeneX, "不存在于表达矩阵中，请检查数据是否包含该基因。"))
  }
  
  gene_expression <- exprs_matrix[GeneX, ]
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Sample     = names(gene_expression),
    Condition  = colData$condition
  )
  
  group_order <- c("NC", "NASH")
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  if (all(c("NC", "NASH") %in% gene_df$Condition)) {
    comparisons_list <- list(c("NC", "NASH"))
  } else {
    comparisons_list <- list()
    warning("NC 和 NASH 组中有一组样本数不足，无法进行统计检验，仅显示箱线图。")
  }
  
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  p <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color) +
    stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  
  p <- p + labs(
    title = "GSE24807",  # 仅显示GSE号
    x = "",             # 去掉 x 轴标签
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(p)
}

# GSE37031
load("~/R/RNA/data/NASH/GSE37031/Data_GSE37031.RData")
analyze_gene_expression_GSE37031 <- function(GeneX, exprs_matrix, colData) {
  group_var <- "condition"
  if (!(GeneX %in% rownames(exprs_matrix))) {
    stop(paste("目标基因", GeneX, "不存在于表达矩阵中，请检查数据是否包含该基因。"))
  }
  if (!(group_var %in% colnames(colData))) {
    stop(paste("指定的分组变量", group_var, "不存在于 colData 中，请检查。"))
  }
  
  gene_expression <- exprs_matrix[GeneX, ]
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Sample     = names(gene_expression),
    Condition  = colData[[group_var]]
  )
  
  group_order <- c("NC", "NASH")
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  if (all(c("NC", "NASH") %in% gene_df$Condition)) {
    comparisons_list <- list(c("NC", "NASH"))
  } else {
    comparisons_list <- list()
    warning("NC 和 NASH 组中有一组样本数不足，无法进行统计检验，仅显示箱线图。")
  }
  
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  p <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color)
  
  if (length(comparisons_list) > 0) {
    p <- p + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  }
  
  p <- p + labs(
    title = "GSE37031",
    x = "",
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(p)
}

# GSE106737
load("~/R/RNA/data/NASH/GSE106737/Data_GSE106737.RData")
analyze_gene_expression_GSE106737 <- function(GeneX, exprs_matrix, colData) {
  group_var <- "condition"
  if (!(GeneX %in% rownames(exprs_matrix))) {
    stop(paste("目标基因", GeneX, "不存在于表达矩阵中，请检查数据是否包含该基因。"))
  }
  if (!(group_var %in% colnames(colData))) {
    stop(paste("指定的分组变量", group_var, "不存在于 colData 中，请检查。"))
  }
  
  gene_expression <- exprs_matrix[GeneX, ]
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Sample     = names(gene_expression),
    Condition  = colData[[group_var]]
  )
  
  default_group_order <- c("NC", "NASH")
  gene_df$Condition <- factor(gene_df$Condition, levels = default_group_order)
  
  group_counts <- table(gene_df$Condition)
  valid_groups <- names(group_counts[group_counts >= 2])
  
  if (length(valid_groups) > 1) {
    comparisons_list <- combn(valid_groups, 2, simplify = FALSE)
  } else {
    comparisons_list <- list()
  }
  
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  p <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color)
  
  if (length(comparisons_list) > 0) {
    p <- p + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的组或有效组小于2个，无法进行统计检验。仅显示箱线图。")
  }
  
  p <- p + labs(
    title = "GSE106737",
    x = "",
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(p)
}

# GSE134146
load("~/R/RNA/data/NASH/GSE134146/Data_GSE134146.RData")
analyze_gene_expression_GSE134146 <- function(GeneX, vsd, colData) {
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因未找到，请检查基因名是否正确或是否存在于表达矩阵中。")
  }
  
  gene_expression <- vsd[GeneX, ]
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition
  )
  
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NASH"))
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color) +
    stat_compare_means(comparisons = list(c("NC", "NASH")), method = "t.test", label = "p.signif")
  
  plot <- plot + labs(
    title = "GSE134146",
    x = "",
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# GSE159676
load("~/R/RNA/data/NASH/GSE159676/Data_GSE159676.RData")
analyze_gene_expression_GSE159676 <- function(GeneX, countData, colData) {
  if (!(GeneX %in% rownames(countData))) {
    stop("目标基因不存在于表达矩阵中，请检查基因名是否正确。")
  }
  
  gene_expression <- countData[GeneX, ]
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Condition = colData$condition
  )
  
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NASH"))
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  custom_colors_fill <- c("NC" = "white", "NASH" = "white")
  custom_colors_color <- c("NC" = "#2ca02c", "NASH" = "#d62728")
  
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color) +
    stat_compare_means(comparisons = list(c("NC", "NASH")), method = "t.test", label = "p.signif")
  
  plot <- plot + labs(
    title = "GSE159676",
    x = "",
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# GSE147304
load("~/R/RNA/data/NASH/GSE147304/Data_GSE147304.RData")
analyze_gene_expression_GSE147304 <- function(GeneX, vsd, colData) {
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition
  )
  
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NASH"))
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color) +
    stat_compare_means(comparisons = list(c("NC", "NASH")), method = "t.test", label = "p.signif")
  
  plot <- plot + labs(
    title = "GSE147304",
    x = "",
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# GSE173735
load("~/R/RNA/data/NASH/GSE173735/Data_GSE173735.RData")
analyze_gene_expression_GSE173735 <- function(GeneX, vsd, colData) {
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因未找到，请检查基因名是否正确或是否存在于表达矩阵中。")
  }
  
  gene_expression <- vsd[GeneX, ]
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition
  )
  
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NASH"))
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color) +
    stat_compare_means(comparisons = list(c("NC", "NASH")), method = "t.test", label = "p.signif")
  
  plot <- plot + labs(
    title = "GSE173735",
    x = "",
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# GSE63067
load("~/R/RNA/data/NASH/GSE63067/Data_GSE63067.RData")
analyze_gene_expression_GSE63067 <- function(GeneX, exprs_matrix, colData) {
  # 固定分组变量名称
  group_var <- "condition"
  
  # Step 1: 检查目标基因是否存在于表达矩阵中
  if (!(GeneX %in% rownames(exprs_matrix))) {
    stop(paste("目标基因", GeneX, "不存在于表达矩阵中，请检查数据是否包含该基因。"))
  }
  
  # Step 2: 检查分组变量是否存在于 colData 中
  if (!(group_var %in% colnames(colData))) {
    stop(paste("指定的分组变量", group_var, "不存在于 colData 中，请检查。"))
  }
  
  # Step 3: 提取目标基因表达数据
  gene_expression <- exprs_matrix[GeneX, ]
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Sample     = names(gene_expression),
    Condition  = colData[[group_var]]
  )
  
  # Step 4: 保留所有组，不过滤（以便绘图显示3组）
  # Step 5: 固定分组顺序：这里设定顺序为 "NC", "NASH", "Obese"
  group_order <- c("NC", "Obese", "NASH")
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 6: 统计比较：仅比较 "NC" 与 "NASH"
  if (all(c("NC", "NASH") %in% gene_df$Condition)) {
    comparisons_list <- list(c("NC", "NASH"))
  } else {
    comparisons_list <- list()
    warning("NC 和 NASH 组中有一组样本数不足，无法进行统计检验，仅显示箱线图。")
  }
  
  # Step 7: 设置颜色映射  
  # 填充色统一为白色；边框颜色使用全局函数 myplotColors 返回的颜色集合
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: 绘制箱线图
  p <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color)
  
  # Step 9: 添加统计检验（仅对 NC 与 NASH 进行 t 检验）
  if (length(comparisons_list) > 0) {
    p <- p + stat_compare_means(
      comparisons = comparisons_list,
      method = "t.test",
      label = "p.signif"
    )
  }
  
  # Step 10: 添加标题和主题
  p <- p + labs(
    title = paste("GSE63067"),
    x = group_var,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(p)
}

# GSE61260
load("~/R/RNA/data/NASH/GSE61260/Data_GSE61260.RData")
analyze_gene_expression_GSE61260 <- function(GeneX, vsd, colData) {
  # Step 1: 检查基因是否存在于表达矩阵
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 获取目标基因的表达值
  gene_expression <- vsd[GeneX, ]
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Condition = colData$condition
  )
  
  # Step 3: 仅保留 NC, Obese, NAFLD, NASH
  valid_conditions <- c("NC", "Obese", "NAFLD", "NASH")
  gene_df <- subset(gene_df, Condition %in% valid_conditions)
  
  # Step 4: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = valid_conditions)
  
  # Step 5: 统计学比较（NC 作为基准）
  comparisons_list <- list(
    c("NC", "Obese"),
    c("NC", "NAFLD"),
    c("NC", "NASH")
  )
  
  # Step 6: 自定义颜色
  custom_colors_fill <- c("NC" = "white", "Obese" = "white", "NAFLD" = "white", "NASH" = "white")
  custom_colors_color <- setNames(getplotColors(length(valid_conditions)), valid_conditions)
  
  # Step 7: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_compare_means(comparisons = comparisons_list,
                       method = "t.test",
                       label = "p.signif") +  
    labs(
      title = paste("GSE61260"),
      x = "Condition",
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE66676/Data_GSE66676.RData")
analyze_gene_expression_GSE66676 <- function(GeneX, vsd, colData) {
  # Step 1: 检查基因是否存在于表达矩阵
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 获取目标基因的表达值
  gene_expression <- vsd[GeneX, ]  # vsd 现在是数据框
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Condition = colData$condition
  )
  
  # Step 3: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NAFLD", "Borderline_NASH", "Definite_NASH"))
  
  # 按 Condition 因子水平排序数据框
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 4: 自定义颜色
  custom_colors_fill <- c("NC" = "white", "NAFLD" = "white", "Borderline_NASH" = "white", "Definite_NASH" = "white")
  custom_colors_color <- getplotColors(length(levels(gene_df$Condition)))
  # Step 5: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  # 边框颜色根据 Condition 设置
    scale_fill_manual(values = custom_colors_fill) +  # 自定义填充颜色
    scale_color_manual(values = custom_colors_color) +  # 自定义边框颜色
    stat_compare_means(comparisons = list(c("NC", "NAFLD"), 
                                          c("NC", "Borderline_NASH"),
                                          c("NC", "Definite_NASH")),
                       method = "t.test",
                       label = "p.signif") + # 两两比较，显示显著性标记
    labs(
      title = paste("GSE66676"),
      x = "Condition",
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

#
load("~/R/RNA/data/NASH/GSE89632/Data_GSE89632.RData")
analyze_gene_expression_GSE89632 <- function(GeneX, norm_data, colData, group) {
  # Step 1: 检查并匹配目标基因在 norm_data$genes$Symbol 中
  if (!("genes" %in% names(norm_data))) {
    stop("norm_data 对象中没有 'genes' 信息。")
  }
  if (!("Symbol" %in% colnames(norm_data$genes))) {
    stop("norm_data$genes 中不包含 'Symbol' 列。")
  }
  idx <- which(norm_data$genes$Symbol == GeneX)
  if (length(idx) == 0) {
    stop(paste("目标基因", GeneX, "在 norm_data$genes$Symbol 中未找到。"))
  }
  # 如果有多个匹配，则取第一个（也可以根据需要修改策略）
  gene_expression <- norm_data$E[idx[1], ]
  
  # Step 2: 检查分组变量是否存在于 colData 中
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'condition', 'fibrosis' 或 'nafld'。")
  }
  
  # Step 3: 构建绘图数据框，提取目标基因表达和对应分组信息
  gene_df <- data.frame(Expression = gene_expression,
                        Condition  = colData[[group]])
  
  # 筛除 Condition 列中为 NA 或字符串 "NA" 的样本
  gene_df <- subset(gene_df, !is.na(Condition) & Condition != "NA")
  
  # Step 4: 根据不同分组设定固定的分组顺序（所有组都显示）
  group_order <- switch(group,
                        "condition" = c("NC", "Obese", "NASH"),
                        "fibrosis"  = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4", "NA"),
                        "nafld"     = c("nafld_0", "nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_6", "nafld_8", "NA"),
                        unique(gene_df$Condition)  # 默认按原顺序
  )
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算各组样本数量，并确定可以进行统计检验的组（样本数 >= 2）
  group_counts <- table(gene_df$Condition)
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: 自动生成 comparisons 列表，仅比较有效组中第一个分组与其他分组
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()
    warning("没有足够的样本进行统计检验，仅显示箱线图。")
  }
  
  # Step 7: 设置颜色映射
  # 填充色设为白色，边框颜色采用 myplotColors 函数生成的调色板
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color)
  
  # Step 9: 添加统计检验（仅在有效组数 >= 2 时）
  if (length(comparisons_list) > 0) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list,
                                      method = "t.test",
                                      label = "p.signif")
  }
  
  # Step 10: 添加标题和图形主题
  plot <- plot + labs(
    title = paste("GSE89632"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE48452/Data_GSE48452.RData")
analyze_gene_expression_GSE48452 <- function(GeneX, expr_data, meta_data, group) {
  # Step 1: 确保基因存在
  if (!(GeneX %in% rownames(expr_data))) {
    stop(paste("目标基因", GeneX, "不存在于表达矩阵中，请检查数据。"))
  }
  
  # Step 2: 提取目标基因的表达数据
  gene_expression <- expr_data[GeneX, ]
  
  # Step 3: 确保分组变量存在
  if (!(group %in% colnames(meta_data))) {
    stop(paste("分组变量", group, "不存在，请检查。"))
  }
  
  # Step 4: 创建数据框
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Condition  = meta_data[[group]]
  )
  
  # Step 5: 移除 NA 值
  gene_df <- subset(gene_df, !is.na(Condition) & Condition != "fibrosis_NA")
  
  # Step 6: 设定分组顺序
  group_order <- switch(group,
                        "condition" = c("NC", "Obese", "Steatosis", "NASH"),
                        "fibrosis"  = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        "nas"       = c("nas_0", "nas_1", "nas_2", "nas_3", "nas_4", "nas_5", "nas_6", "nas_7"),
                        unique(gene_df$Condition))  # 默认按数据顺序
  
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 7: 统计学比较（基准为分组的第一个）
  group_counts <- table(gene_df$Condition)
  valid_groups <- names(group_counts[group_counts >= 2])
  
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()
    warning("没有足够的样本进行统计检验，仅显示箱线图。")
  }
  
  # Step 8: 设置颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 9: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(width = 0.2, size = 1, alpha = 0.6)+
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color)
  
  # Step 10: 添加统计检验（仅在有效组数 >= 2 时）
  if (length(comparisons_list) > 0) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list,
                                      method = "t.test",
                                      label = "p.signif")
  }
  
  # Step 11: 添加标题和图形主题
  plot <- plot + labs(
    title = paste("GSE89632"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE83452/Data_GSE83452.RData")
analyze_gene_expression_GSE83452 <- function(GeneX,  exprs_matrix, colData  ) {
  # 固定分组变量
  group_var <- "condition"
  
  # Step 1: 检查目标基因是否存在
  if (!(GeneX %in% rownames(exprs_matrix))) {
    stop(paste("目标基因", GeneX, "不存在于表达矩阵中，请检查数据是否包含该基因。"))
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group_var %in% colnames(colData))) {
    stop(paste("指定的分组变量", group_var, "不存在于 colData 中，请检查。"))
  }
  
  # Step 3: 提取目标基因表达
  gene_expression <- exprs_matrix[GeneX, ]
  
  # 构建绘图数据框
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Sample     = names(gene_expression),
    Condition  = colData[[group_var]]
  )
  
  # Step 4: 过滤掉 undefined 组
  valid_conditions <- c("NC_baseline", "NASH_baseline")
  gene_df <- subset(gene_df, Condition %in% valid_conditions)
  
  # Step 5: 固定分组顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = valid_conditions)
  
  # Step 6: 统计各组样本数
  group_counts <- table(gene_df$Condition)
  
  # 仅保留样本数 >= 2 的组用于统计检验
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 7: comparisons 列表（基于 `valid_groups`）
  comparisons_list <- list(c("NC_baseline", "NASH_baseline"))
  
  # Step 8: 自定义颜色
  custom_colors_fill <- c("NC_baseline" = "white", "NASH_baseline" = "white")
  custom_colors_color <- setNames(getplotColors(length(valid_conditions)), valid_conditions)
  
  # Step 9: 绘制箱线图
  p <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color)
  
  # Step 10: 统计检验
  if (length(comparisons_list) > 0) {
    p <- p + stat_compare_means(
      comparisons = comparisons_list,
      method = "t.test",
      label = "p.signif"
    )
  } else {
    warning("没有足够的组或有效组小于2个，无法进行统计检验。仅显示箱线图。")
  }
  
  # Step 11: 添加标题和主题
  p <- p + labs(
    title = paste("GSE83452"),
    x = "Condition",
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(p)
}

# 
load("~/R/RNA/data/NASH/GSE105127/Data_GSE105127.RData")
analyze_gene_expression_GSE105127 <- function(GeneX, vsd, colData) {
  # Step 1: 检查基因是否存在于表达矩阵
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 获取目标基因的表达值
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition
  )
  
  # Step 3: 设置 Condition 的因子顺序
  # 按分组名称排序（CV -> IZ -> PP，NC -> STEA -> HO -> EARLY_NASH）
  condition_levels <- c(
    "CV_NC", "CV_STEA", "CV_HO", "CV_EARLY_NASH",
    "IZ_NC", "IZ_STEA", "IZ_HO", "IZ_EARLY_NASH",
    "PP_NC", "PP_STEA", "PP_HO", "PP_EARLY_NASH"
  )
  gene_df$Condition <- factor(gene_df$Condition, levels = condition_levels)
  
  # Step 4: 按 Condition 因子水平排序数据框
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 5: 定义颜色
  # 按分组前缀（CV/IZ/PP）和后缀（NC/STEA/HO/EARLY_NASH）定义颜色
  custom_colors_fill <- c(
    "CV_NC" = "#FFFFFF", "CV_STEA" = "#FFFFFF", "CV_HO" = "#FFFFFF", "CV_EARLY_NASH" = "#FFFFFF",
    "IZ_NC" = "#FFFFFF", "IZ_STEA" = "#FFFFFF", "IZ_HO" = "#FFFFFF", "IZ_EARLY_NASH" = "#FFFFFF",
    "PP_NC" = "#FFFFFF", "PP_STEA" = "#FFFFFF", "PP_HO" = "#FFFFFF", "PP_EARLY_NASH" = "#FFFFFF"
  )
  
  custom_colors_color <- c(
    "CV_NC" = "#2ca02c", "CV_STEA" = "#ff7f0e", "CV_HO" = "#d62728", "CV_EARLY_NASH" = "#9467bd",
    "IZ_NC" = "#2ca02c", "IZ_STEA" = "#ff7f0e", "IZ_HO" = "#d62728", "IZ_EARLY_NASH" = "#9467bd",
    "PP_NC" = "#2ca02c", "PP_STEA" = "#ff7f0e", "PP_HO" = "#d62728", "PP_EARLY_NASH" = "#9467bd"
  )
  
  # Step 6: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  # 边框颜色根据 Condition 设置
    scale_fill_manual(values = custom_colors_fill) +  # 自定义填充颜色
    scale_color_manual(values = custom_colors_color) +  # 自定义边框颜色
    stat_compare_means(
      comparisons = list(
        c("CV_NC", "CV_STEA"), c("CV_NC", "CV_HO"), c("CV_NC", "CV_EARLY_NASH"),
        c("IZ_NC", "IZ_STEA"), c("IZ_NC", "IZ_HO"), c("IZ_NC", "IZ_EARLY_NASH"),
        c("PP_NC", "PP_STEA"), c("PP_NC", "PP_HO"), c("PP_NC", "PP_EARLY_NASH")
      ),
      method = "t.test",
      label = "p.signif"
    ) + # 两两比较，显示显著性标记
    labs(
      title = paste("GSE105127"),
      x = "Condition",
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE115193/Data_GSE115193.RData")
analyze_gene_expression_GSE115193 <- function(GeneX, vsd, colData) {
  # Step 1: 检查基因是否存在于表达矩阵
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 获取目标基因的表达值
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition
  )
  
  # Step 3: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NAFLD", "NASH"))
  
  # 按 Condition 因子水平排序数据框
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 4: 绘制箱线图
  custom_colors_fill <- c("NC" = "white", "NAFLD" = "white", "NASH" = "white")
  custom_colors_color <- c("NC" = "#2ca02c", "NAFLD" = "#ff7f0e", "NASH" = "#d62728")
  
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  # 边框颜色根据 Condition 设置
    scale_fill_manual(values = custom_colors_fill) +  # 自定义填充颜色
    scale_color_manual(values = custom_colors_color) +  # 自定义边框颜色
    stat_compare_means(comparisons = list(c("NC", "NAFLD"), 
                                          c("NC", "NASH"), 
                                          c("NASH", "NAFLD")),
                       method = "t.test",
                       label = "p.signif") + # 两两比较，显示显著性标记
    labs(
      title = paste("GSE115193"),
      x = "Condition",
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE126848/Data_GSE126848.RData")
analyze_gene_expression_GSE126848 <- function(GeneX, vsd, colData) {
  # Step 1: 检查基因是否存在于表达矩阵
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 获取目标基因的表达值
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition
  )
  
  # Step 3: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "Obese", "NAFLD", "NASH"))
  
  # 按 Condition 因子水平排序数据框
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 4: 绘制箱线图
  custom_colors_fill <- c("NC" = "white", "Obese" = "white", "NAFLD" = "white", "NASH" = "white")
  custom_colors_color <- c("NC" = "#2ca02c", "Obese" = "#FFC07C", "NAFLD" = "#ff7f0e", "NASH" = "#d62728")
  
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  # 边框颜色根据 Condition 设置
    scale_fill_manual(values = custom_colors_fill) +  # 自定义填充颜色
    scale_color_manual(values = custom_colors_color) +  # 自定义边框颜色
    stat_compare_means(comparisons = list(c("NC", "Obese"), 
                                          c("NC", "NAFLD"), 
                                          c("NC", "NASH"), 
                                          c("Obese", "NASH"), 
                                          c("Obese", "NAFLD"), 
                                          c("NASH", "NAFLD")),
                       method = "t.test",
                       label = "p.signif") + # 两两比较，显示显著性标记
    labs(
      title = paste("GSE126848"),
      x = "Condition",
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE135251/Data_GSE135251.RData")
analyze_gene_expression_GSE135251 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nas', 'fibrosis' 或 'steatohep'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  
  # Step 5: 设定分组顺序
  group_order <- switch(group,
                        "nas" = c("nas_0", "nas_1", "nas_2", "nas_3", "nas_4", "nas_5", "nas_6", "nas_7", "nas_8"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3", "fibrosis_4"),
                        "level" = c("NC", "NAFL", "NASH"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # 确保分组顺序只包含数据中的分组
  group_order <- intersect(group_order, unique(gene_df$Condition))
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 6: **自动生成 comparisons**
  condition_levels <- levels(gene_df$Condition)
  comparisons_list <- lapply(condition_levels[-1], function(x) c(condition_levels[1], x))
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(unique(gene_df$Condition))), unique(gene_df$Condition))
  custom_colors_color <- setNames(getplotColors(length(unique(gene_df$Condition))), unique(gene_df$Condition))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) + 
    stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif") +  # 统计显著性
    labs(
      title = paste("GSE135251"),
      x = group,
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE164760/Data_GSE164760.RData")
analyze_gene_expression_GSE164760 <- function(GeneX,exprs_matrix, colData ) {
  # 固定分组变量
  group_var <- "level"
  
  # Step 1: 检查目标基因是否存在
  if (!(GeneX %in% rownames(exprs_matrix))) {
    stop(paste("目标基因", GeneX, "不存在于表达矩阵中，请检查数据是否包含该基因。"))
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group_var %in% colnames(colData))) {
    stop(paste("指定的分组变量", group_var, "不存在于 colData 中，请检查。"))
  }
  
  # Step 3: 提取目标基因表达
  gene_expression <- exprs_matrix[GeneX, ]
  
  # 构建绘图数据框
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Sample     = names(gene_expression),
    Condition  = colData[[group_var]]
  )
  
  # Step 4: 如果需要固定分组顺序（自行修改）
  default_group_order <- c("Healthy", "NASH", "Cirrhotic", "NASH_HCC_notumor", "NASH_HCC_tumor")
  gene_df$Condition <- factor(gene_df$Condition, levels = default_group_order)
  
  # Step 5: 统计各组样本数
  group_counts <- table(gene_df$Condition)
  
  # 仅保留样本数 >= 2 的组用于统计检验
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: comparisons 列表
  # 使用 combn() 对所有有效组进行两两比较
  if (length(valid_groups) > 1) {
    comparisons_list <- combn(valid_groups, 2, simplify = FALSE)
  } else {
    comparisons_list <- list()
  }
  
  # Step 7: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 8: 绘制箱线图
  p <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_jitter(width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color)
  
  # Step 9: 统计检验
  if (length(comparisons_list) > 0) {
    p <- p + stat_compare_means(
      comparisons = comparisons_list,
      method = "t.test",
      label = "p.signif"
    )
  } else {
    warning("没有足够的组或有效组小于2个，无法进行统计检验。仅显示箱线图。")
  }
  
  # Step 10: 添加标题和主题
  p <- p + labs(
    title = paste("GSE164760"),
    x = group_var,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(p)
}

# 
load("~/R/RNA/data/NASH/GSE167523/Data_GSE167523.RData")
analyze_gene_expression_GSE167523 <- function(GeneX, vsd, colData) {
  # Step 1: 检查目标基因是否存在于表达矩阵
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因未找到，请检查基因名是否正确或是否存在于表达矩阵中。")
  }
  
  # Step 2: 获取目标基因的表达值
  gene_expression <- assay(vsd)[GeneX, ]  # 使用 assay() 提取表达矩阵
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$level  # 使用 level 列作为分组信息
  )
  
  # Step 3: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NAFL", "NASH"))
  
  # 按 Condition 因子水平排序数据框
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 4: 绘制箱线图
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition), width = 0.6) +  # 设置箱线图宽度
    # geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, aes(color = Condition)) +  # 添加散点
    scale_fill_manual(values = custom_colors_fill) +  # 自定义填充颜色
    scale_color_manual(values = custom_colors_color) +  # 自定义边框和散点颜色
    stat_compare_means(
      comparisons = list(c("NAFL", "NASH")), 
      method = "t.test",  # 使用 t 检验
      label = "p.signif",  # 显示显著性标记
      label.y = max(gene_df$Expression) * 1.1  # 将显著性标记放在图的上方
    ) +
    labs(
      title = paste("GSE167523"),
      x = "Condition",
      y = "Expression Level (VST)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # 调整 x 轴标签
      plot.title = element_text(hjust = 0.5),  # 调整标题
      legend.position = "none"  # 移除图例
    )
  
  # 返回绘图
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE185051/Data_GSE185051.RData")
analyze_gene_expression_GSE185051 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nas', 'fibrosis' 或 'level'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nas" = c("nas_0", "nas_1", "nas_2", "nas_3", "nas_4", "nas_5", "nas_6", "nas_7"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_2", "fibrosis_3"),
                        "level" = c("Normal", "Steatosis", "Borderline", "NASH"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # 确保所有组都显示，即使某些组样本数 <2
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # 确定可以进行统计检验的组（样本数 >= 2）
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE185051"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE207310/Data_GSE207310.RData")
analyze_gene_expression_GSE207310 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否为 'nafld', 'fibrosis' 或 'level'。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序（所有组都要显示）
  group_order <- switch(group,
                        "nafld" = c("nafld_0", "nafld_1", "nafld_2", "nafld_3", "nafld_4", "nafld_5", "nafld_7", "nafld_8"),
                        "fibrosis" = c("fibrosis_0", "fibrosis_1", "fibrosis_1A", "fibrosis_1B", "fibrosis_2", "fibrosis_3"),
                        "level" = c("NC", "NAFL", "NASH"),
                        unique(gene_df$Condition))  # 默认按原顺序
  
  # **确保所有组都显示，即使某些组样本数 <2**
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # **确定可以进行统计检验的组（样本数 >= 2）**
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: **自动生成 comparisons（仅比较有效组）**
  condition_levels <- levels(gene_df$Condition)
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c(valid_groups[1], x))
  } else {
    comparisons_list <- list()  # 避免只有一个分组时报错
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 8: **绘制箱线图（所有组都要显示）**
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
    scale_fill_manual(values = custom_colors_fill) + 
    scale_color_manual(values = custom_colors_color) 
  
  # **仅在组数 >=2 且样本足够时进行统计检验**
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("GSE207310"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE260666/Data_GSE260666.RData")
analyze_gene_expression_GSE260666 <- function(GeneX, vsd, colData) {
  # Step 1: 检查基因是否存在于表达矩阵
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 获取目标基因的表达值
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition
  )
  
  # Step 3: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NAFLD", "NASH"))
  
  # 按 Condition 因子水平排序数据框
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  # Step 4: 绘制箱线图
  custom_colors_fill <- c("NC" = "white", "NAFLD" = "white", "NASH" = "white")
  custom_colors_color <- c("NC" = "#2ca02c", "NAFLD" = "#ff7f0e", "NASH" = "#d62728")
  
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  # 边框颜色根据 Condition 设置
    scale_fill_manual(values = custom_colors_fill) +  # 自定义填充颜色
    scale_color_manual(values = custom_colors_color) +  # 自定义边框颜色
    stat_compare_means(comparisons = list(c("NC", "NAFLD"), 
                                          c("NC", "NASH"), 
                                          c("NASH", "NAFLD")),
                       method = "t.test",
                       label = "p.signif") + # 两两比较，显示显著性标记
    labs(
      title = paste("GSE260666"),
      x = "Condition",
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

#
load("~/R/RNA/data/NASH/GSE274114/Data_GSE274114.RData")
analyze_gene_expression_GSE274114 <- function(GeneX, vsd_object, colData) {
  # 固定分组变量名称
  group_var <- "condition"
  
  # Step 1: 检查目标基因是否存在于表达矩阵中
  # 从 S4 对象中提取表达矩阵
  if (is(vsd_object, "DESeqTransform") || is(vsd_object, "SummarizedExperiment")) {
    exprs_matrix <- assay(vsd_object)  # 提取表达矩阵
  } else {
    stop("输入的表达数据必须是 DESeqTransform 或 SummarizedExperiment 对象。")
  }
  
  if (!(GeneX %in% rownames(exprs_matrix))) {
    stop(paste("目标基因", GeneX, "不存在于表达矩阵中，请检查数据是否包含该基因。"))
  }
  
  # Step 2: 检查分组变量是否存在于 colData 中
  if (!(group_var %in% colnames(colData))) {
    stop(paste("指定的分组变量", group_var, "不存在于 colData 中，请检查。"))
  }
  
  # Step 3: 提取目标基因表达数据
  gene_expression <- exprs_matrix[GeneX, ]
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Sample     = names(gene_expression),
    Condition  = colData[[group_var]]
  )
  
  # Step 4: 保留所有组，不过滤（以便绘图显示4组）
  # 修改点：按照数据实际因子顺序设置分组顺序
  group_order <- c("NC", "ENEG", "NASH", "ENEG_NASH")  # 直接使用数据中的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 统计比较：仅比较 "NC" vs "ENEG" 和 "NASH" vs "ENEG_NASH"
  comparisons_list <- list()
  
  if (all(c("NC", "ENEG") %in% gene_df$Condition)) {
    comparisons_list <- c(comparisons_list, list(c("NC", "ENEG")))
  } else {
    warning("NC 和 ENEG 组中有一组样本数不足，无法进行统计检验。")
  }
  
  if (all(c("NASH", "ENEG_NASH") %in% gene_df$Condition)) {
    comparisons_list <- c(comparisons_list, list(c("NASH", "ENEG_NASH")))
  } else {
    warning("NASH 和 ENEG_NASH 组中有一组样本数不足，无法进行统计检验。")
  }
  
  # Step 6: 设置颜色映射  
  # 填充色统一为白色；边框颜色使用全局函数 myplotColors 返回的颜色集合
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  
  # Step 7: 绘制箱线图
  p <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition, color = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color)
  
  # Step 8: 添加统计检验（仅对 "NC" vs "ENEG" 和 "NASH" vs "ENEG_NASH" 进行 t 检验）
  if (length(comparisons_list) > 0) {
    p <- p + stat_compare_means(
      comparisons = comparisons_list,
      method = "t.test",
      label = "p.signif"
    )
  }
  
  # Step 9: 添加标题和主题
  p <- p + labs(
    title = paste("GSE274114"),
    x = group_var,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(p)
}

# 
load("~/R/RNA/data/NASH/GSE288077/Human/Data_GSE288077.RData")
analyze_gene_expression_GSE288077_Human <- function(GeneX, vsd, colData) {
  # Step 1: 直接检查基因是否存在于表达矩阵
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 获取目标基因的表达值
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition  # 直接使用 colData$condition
  )
  # Step 5: 固定分组顺序：这里设定顺序为 "NC", "NASH", "Obese"
  group_order <- c("ENEG", "ENEG_NASH")
  # Step 3: 绘制箱线图
  # 填充色统一为白色；边框颜色使用全局函数 myplotColors 返回的颜色集合
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(getplotColors(length(group_order)), group_order)
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  # 边框颜色根据 Condition 设置
    scale_fill_manual(values = custom_colors_fill) +  # 自定义填充颜色
    scale_color_manual(values = custom_colors_color) +  # 自定义边框颜色
    stat_compare_means(comparisons = list(c("ENEG", "ENEG_NASH")), 
                       method = "t.test", 
                       label = "p.signif") +  # 仅比较 ENEG 与 ENEG_NASH
    labs(
      title = paste("GSE288077"),
      x = "Condition",
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

#
load("~/R/RNA/data/NASH/GSE163211/Data_GSE163211.RData")
analyze_gene_expression_GSE163211 <- function(GeneX, vsd, colData, group) {
  # Step 1: 检查基因是否存在
  if (!(GeneX %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 2: 检查分组变量是否存在
  if (!(group %in% colnames(colData))) {
    stop("指定的分组变量不存在，请检查是否正确。")
  }
  
  # Step 3: 提取目标基因的表达数据
  gene_expression <- assay(vsd)[GeneX, ]
  gene_df <- data.frame(Expression = gene_expression, Condition = colData[[group]])
  
  # Step 4: 设定分组顺序
  group_order <- c("Normal", "Steatosis", "NASH_F0", "NASH_F1_F4")
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # Step 5: 计算每个组的样本数量
  group_counts <- table(gene_df$Condition)
  
  # 仅保留样本数 >= 2 的组
  valid_groups <- names(group_counts[group_counts >= 2])
  
  # Step 6: 生成 comparisons
  if (length(valid_groups) > 1) {
    comparisons_list <- lapply(valid_groups[-1], function(x) c("Normal", x))
  } else {
    comparisons_list <- list()  # 避免错误
  }
  
  # Step 7: 颜色映射
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(rainbow(length(group_order)), group_order)
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color)
  
  # Step 9: 添加统计检验
  if (length(valid_groups) > 1) {
    plot <- plot + stat_compare_means(comparisons = comparisons_list, method = "t.test", label = "p.signif")
  } else {
    warning("没有足够的样本进行统计检验，仅显示所有组")
  }
  
  # 添加标题
  plot <- plot + labs(
    title = paste("Expression of", GeneX, "in GSE163211"),
    x = group,
    y = "Expression Level"
  ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE151158/Data_GSE151158.RData")
analyze_gene_expression_GSE151158 <- function(GeneX, countData, colData) {
  # Step 1: 检查基因是否存在于表达矩阵
  if (!(GeneX %in% rownames(countData))) {
    stop("目标基因不存在于表达矩阵中，请检查基因名是否正确。")
  }
  
  # Step 2: 提取目标基因的表达值
  gene_expression <- countData[GeneX, ]
  gene_df <- data.frame(
    Expression = as.numeric(gene_expression),
    Condition = colData$nas  # 使用 nas 作为分组变量
  )
  
  # Step 3: 设置 Condition 的因子顺序
  group_order <- c("nas_0", "nas_2", "nas_3", "nas_5", "nas_6")  # 指定分组顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = group_order)
  
  # 按 Condition 因子水平排序数据框
  gene_df <- gene_df[order(gene_df$Condition), ]
  
  
  # Step 5: 绘制箱线图
  custom_colors_fill <- setNames(rep("white", length(group_order)), group_order)
  custom_colors_color <- setNames(rainbow(length(group_order)), group_order)
  
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  # 边框颜色根据 Condition 设置
    scale_fill_manual(values = custom_colors_fill) +  # 自定义填充颜色
    scale_color_manual(values = custom_colors_color) +  # 自定义边框颜色
    stat_compare_means(comparisons = list(c("nas_0", "nas_3"), c("nas_0", "nas_5"), c("nas_0", "nas_6")),  # 进行 t 检验，比较 nas_0 和其他组
                       method = "t.test",
                       label = "p.signif") +  # 显示显著性标记
    labs(
      title = paste("Expression of", GeneX, "in GSE151158"),
      x = "Condition",
      y = "Expression Level"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  # 返回绘图
  return(plot)
}






# 网页端设置 ----------
# 预设的用户名和密码
USER <- "admin"
PASSWORD <- "1234"

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      /* 主容器样式 */
      body {
        position: relative;
        margin: 0;
        height: 100vh;
      }
      
      /* 登录界面居中样式 */
      #login_panel {
        position: absolute;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        background: #f8f9fa;
        padding: 30px;
        border-radius: 10px;
        box-shadow: 0 0 15px rgba(0,0,0,0.1);
        z-index: 1000;
      }
      
      /* 侧边栏样式 */
      #sidebar_panel {
        position: fixed;
        top: 0;
        left: 0;
        width: 300px;
        height: 100vh;
        background: #f0f0f0;
        padding: 20px;
        border-right: 1px solid #ddd;
        z-index: 999;
      }
      
      /* 主内容区样式 */
      .main-content {
        margin-left: 320px;
        padding: 20px;
        height: 100vh;
        overflow-y: auto;
      }
      
     /* 自定义 titlePanel 样式 */
      #sidebar_panel .title {
        font-size: 18px; /* 调整字体大小 */
        white-space: nowrap; /* 防止换行 */
        overflow: hidden; /* 超出部分隐藏 */
        text-overflow: ellipsis; /* 超出部分显示省略号 */
        text-align: center;  /* 居中文本 */
        width: 100%; /* 确保元素的宽度为 100% */
        padding: 10px 0; /* 使用 padding 来增加高度 */
        line-height: 1.4; /* 确保文本垂直居中 */
        margin-bottom: 20px; /* 增加 titlePanel 和下面内容的间距 */
      }
      
    "))
  ),
  
  # 登录面板
  conditionalPanel(
    condition = "!output.user_authenticated",
    div(id = "login_panel",
        h3("用户登录", style = "text-align: center;"),
        textInput("user", "用户名", value = "admin", width = "100%"),
        passwordInput("password", "密码", value = "1234", width = "100%"),
        div(actionButton("login", "登录", class = "btn-primary"), 
            style = "text-align: center; margin-top: 15px;"),
        textOutput("login_status")
    )
  ),
  
  # 登录后的界面
  conditionalPanel(
    condition = "output.user_authenticated",
    div(id = "sidebar_panel",
        div(class = "title", HTML("MGA <br>（Mash Genome Atlas）")),  # 自定义 titlePanel
        radioButtons("species", "物种", 
                     choices = c("Human","Mouse（建设中）")),  
        
        conditionalPanel(
          condition = "input.species == 'Human'",
          radioButtons("analysis_type", "选择分析类型", 
                       choices = c("Summary","Classification","NAFLD Score", "Fibrosis Stages" ,"scRNAseq"))  # 根据选择显示子分析选项
        ),
        conditionalPanel(
          condition = "input.species == 'Mouse（建设中）'",
          radioButtons("feed_type", "造模方式", 
                       choices = c("WD","HFHC", "MCD", "HFD+CDA"))  # 根据选择显示子分析选项
        ),
        # 当选择 "单细胞转录组" 时显示细胞类型选择框
        conditionalPanel(
          condition = "input.analysis_type == 'scRNAseq'", 
          selectInput("cell_type", "选择细胞类型", 
                      choices = c("Live cell", "Macrophages", "CD4+Tcell", "CD8+Tcell", "B cell", "NK cell", "Endothelial cell"))  # 你可以根据需求修改这些细胞类型
        ),
        textInput("gene", "输入基因名称", value = "COL1A1"),
        actionButton("update", "更新图表")
    ),
    
    div(class = "main-content",
        uiOutput("dynamic_plot")
    )
    
    
    # div(class = "main-content",
    #     plotOutput("combined_plot", width = "100%", height = "1000px")
    # )
  )
)


# 服务器端逻辑
server <- function(input, output, session) {
  
  # 用户身份验证状态
  user_auth <- reactiveVal(FALSE)
  
  observeEvent(input$login, {
    if (input$user == USER && input$password == PASSWORD) {
      user_auth(TRUE)
    } else {
      user_auth(FALSE)
    }
  })
  
  output$login_status <- renderText({
    if (user_auth()) {
      "✅ 登录成功！"
    } else if (input$login > 0) {
      "❌ 登录失败，请检查用户名或密码。"
    }
  })
  
  output$user_authenticated <- reactive({
    user_auth()
  })
  outputOptions(output, "user_authenticated", suspendWhenHidden = FALSE)
  
  # 定义分析函数
  analyze_nafld_singlegene <- function(Gene) {
    
    # 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))+ labs(x = NULL)
    }
    
    # 数据集列表
    datasets <- list(
      GSE130970 = list(data = vsd_nafld_GSE130970, colData = colData_GSE130970, func = analyze_gene_expression_GSE130970, group = "nafld"),
      GSE89632  = list(data = Illumina_GSE89632, colData = colData_GSE89632, func = analyze_gene_expression_GSE89632, group = "nafld"),
      GSE48452  = list(data = Affymetrix_GSE48452, colData = colData_GSE48452, func = analyze_gene_expression_GSE48452, group = "nas"),
      GSE135251 = list(data = vsd_nas_GSE135251, colData = colData_GSE135251, func = analyze_gene_expression_GSE135251, group = "nas"),
      GSE174478 = list(data = vsd_nafld_GSE174478, colData = colData_GSE174478, func = analyze_gene_expression_GSE174478, group = "nafld"),
      GSE185051 = list(data = vsd_nas_GSE185051, colData = colData_GSE185051, func = analyze_gene_expression_GSE185051, group = "nas"),
      GSE193066 = list(data = vsd_nafld_GSE193066, colData = colData_GSE193066, func = analyze_gene_expression_GSE193066, group = "nafld"),
      GSE193080 = list(data = vsd_nafld_GSE193080, colData = colData_GSE193080, func = analyze_gene_expression_GSE193080, group = "nafld"),
      GSE207310 = list(data = vsd_nafld_GSE207310, colData = colData_GSE207310, func = analyze_gene_expression_GSE207310, group = "nafld"),
      GSE246221 = list(data = vsd_nafld_GSE246221, colData = colData_GSE246221, func = analyze_gene_expression_GSE246221, group = "nafld"),
      GSE225740 = list(data = vsd_nas_GSE225740, colData = colData_GSE225740, func = analyze_gene_expression_GSE225740, group = "nas")
    )
    
    # 处理所有数据集，捕获可能的错误
    plots <- lapply(names(datasets), function(name) {
      dataset <- datasets[[name]]
      tryCatch({
        plot <- dataset$func(Gene, dataset$data, dataset$colData, dataset$group)
        modify_plot(plot) # 统一格式化
      }, error = function(e) {
        message(paste("跳过", name, "数据集，错误:", e$message))
        NULL  # 遇到错误时返回 NULL，避免影响整个流程
      })
    })
    
    # 过滤掉 NULL 结果
    plots <- Filter(Negate(is.null), plots)
    
    # 组合所有的图（4 列布局）
    combined_plot <- wrap_plots(plots, ncol = 4) + 
      plot_annotation(
        title = paste("Expression of", Gene, "in NAFLD Score"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20))
        )
      )
    
    return(combined_plot)
  }
  
  
  analyze_fibrosis_singlegene <- function(Gene) {
    
    # 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))+ labs(x = NULL)
    }
    
    # 数据集列表
    datasets <- list(
      GSE130970 = list(data = vsd_fibrosis_GSE130970, colData = colData_GSE130970, func = analyze_gene_expression_GSE130970, group = "fibrosis"),
      GSE89632  = list(data = Illumina_GSE89632, colData = colData_GSE89632, func = analyze_gene_expression_GSE89632, group = "fibrosis"),
      GSE48452  = list(data = Affymetrix_GSE48452, colData = colData_GSE48452, func = analyze_gene_expression_GSE48452, group = "fibrosis"),
      GSE135251 = list(data = vsd_fibrosis_GSE135251, colData = colData_GSE135251, func = analyze_gene_expression_GSE135251, group = "fibrosis"),
      GSE174478 = list(data = vsd_fibrosis_GSE174478, colData = colData_GSE174478, func = analyze_gene_expression_GSE174478, group = "fibrosis"),
      GSE185051 = list(data = vsd_fibrosis_GSE185051, colData = colData_GSE185051, func = analyze_gene_expression_GSE185051, group = "fibrosis"),
      GSE193066 = list(data = vsd_fibrosis_GSE193066, colData = colData_GSE193066, func = analyze_gene_expression_GSE193066, group = "fibrosis"),
      GSE193080 = list(data = vsd_fibrosis_GSE193080, colData = colData_GSE193080, func = analyze_gene_expression_GSE193080, group = "fibrosis"),
      GSE207310 = list(data = vsd_fibrosis_GSE207310, colData = colData_GSE207310, func = analyze_gene_expression_GSE207310, group = "fibrosis"),
      GSE246221 = list(data = vsd_fibrosis_GSE246221, colData = colData_GSE246221, func = analyze_gene_expression_GSE246221, group = "fibrosis"),
      GSE225740 = list(data = vsd_fibrosis_GSE225740, colData = colData_GSE225740, func = analyze_gene_expression_GSE225740, group = "fibrosis")
    )
    
    # 处理所有数据集，捕获可能的错误
    plots <- lapply(names(datasets), function(name) {
      dataset <- datasets[[name]]
      tryCatch({
        plot <- dataset$func(Gene, dataset$data, dataset$colData, dataset$group)
        modify_plot(plot) # 统一格式化
      }, error = function(e) {
        message(paste("跳过", name, "数据集，错误:", e$message))
        NULL  # 遇到错误时返回 NULL，避免影响整个流程
      })
    })
    
    # 过滤掉 NULL 结果
    plots <- Filter(Negate(is.null), plots)
    
    # 组合所有的图（4 列布局）
    combined_plot <- wrap_plots(plots, ncol = 4) + 
      plot_annotation(
        title = paste("Expression of", Gene, "in Fibrosis Stages"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20))
        )
      )
    
    return(combined_plot)
  }
  

  analyze_classification_singlegene <- function(Gene) {
    
    # 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))+ labs(x = NULL)
    }
    
    # 数据集列表
    datasets <- list(
      GSE24807 = list(data = codelink_GSE24807, colData = colData_GSE24807, func = analyze_gene_expression_GSE24807),
      GSE37031 = list(data = Affymetrix_GSE37031, colData = colData_GSE37031, func = analyze_gene_expression_GSE37031),
      GSE106737 = list(data = Affymetrix_GSE106737, colData = colData_GSE106737, func = analyze_gene_expression_GSE106737),
      GSE134146 = list(data = Arraystar_GSE134146, colData = colData_GSE134146, func = analyze_gene_expression_GSE134146),
      GSE159676 = list(data = Affymetrix__GSE159676, colData = colData_GSE159676, func = analyze_gene_expression_GSE159676),
      GSE147304 = list(data = vsd_GSE147304, colData = colData_GSE147304, func = analyze_gene_expression_GSE147304),
      GSE173735 = list(data = vsd_GSE173735, colData = colData_GSE173735, func = analyze_gene_expression_GSE173735),
      GSE83452 = list(data = Affymetrix_GSE83452, colData = colData_GSE83452, func = analyze_gene_expression_GSE83452),
      GSE167523 = list(data = vsd_GSE167523, colData = colData_GSE167523, func = analyze_gene_expression_GSE167523),
      GSE288077 = list(data = vsd_GSE288077_Human, colData = colData_GSE288077_Human, func = analyze_gene_expression_GSE288077_Human),
      GSE89632 = list(data = Illumina_GSE89632, colData = colData_GSE89632, func = analyze_gene_expression_GSE89632, group = "condition"),
      GSE63067 = list(data = Affymetrix_GSE63067, colData = colData_GSE63067, func = analyze_gene_expression_GSE63067),
      GSE115193 = list(data = vsd_GSE115193, colData = colData_GSE115193, func = analyze_gene_expression_GSE115193),
      GSE260666 = list(data = vsd_GSE260666, colData = colData_GSE260666, func = analyze_gene_expression_GSE260666),
      GSE135251 = list(data = vsd_level_GSE135251, colData = colData_GSE135251, func = analyze_gene_expression_GSE135251, group = "level"),
      GSE207310 = list(data = vsd_level_GSE207310, colData = colData_GSE207310, func = analyze_gene_expression_GSE207310, group = "level"),
      GSE61260 = list(data = Affymetrix_GSE61260, colData = colData_GSE61260, func = analyze_gene_expression_GSE61260),
      GSE66676 = list(data = Affymetrix_GSE66676, colData = colData_GSE66676, func = analyze_gene_expression_GSE66676),
      GSE126848 = list(data = vsd_GSE126848, colData = colData_GSE126848, func = analyze_gene_expression_GSE126848),
      GSE274114 = list(data = vsd_GSE274114, colData = colData_GSE274114, func = analyze_gene_expression_GSE274114),
      GSE48452 = list(data = Affymetrix_GSE48452, colData = colData_GSE48452, func = analyze_gene_expression_GSE48452, group = "condition"),
      GSE185051 = list(data = vsd_level_GSE185051, colData = colData_GSE185051, func = analyze_gene_expression_GSE185051, group = "level"),
      GSE164760 = list(data = affy_GSE164760, colData = colData_GSE164760, func = analyze_gene_expression_GSE164760),
      GSE105127 = list(data = vsd_GSE105127, colData = colData_GSE105127, func = analyze_gene_expression_GSE105127),
      GSE163211 = list(data = vsd_GSE163211, colData = colData_GSE163211, func = analyze_gene_expression_GSE163211, group = "level"),
      GSE151158 = list(data = RCC_GSE15115, colData = colData_GSE151158, func = analyze_gene_expression_GSE151158)
      
    )
    
    # 处理所有数据集，捕获可能的错误
    plots <- lapply(names(datasets), function(name) {
      dataset <- datasets[[name]]
      tryCatch({
        plot <- if ("group" %in% names(dataset)) {
          dataset$func(Gene, dataset$data, dataset$colData, dataset$group)
        } else {
          dataset$func(Gene, dataset$data, dataset$colData)
        }
        modify_plot(plot) # 统一格式化
      }, error = function(e) {
        message(paste("跳过", name, "数据集，错误:", e$message))
        NULL  # 遇到错误时返回 NULL，避免影响整个流程
      })
    })
    
    # 过滤掉 NULL 结果
    plots <- Filter(Negate(is.null), plots)
    
    # 定义布局
    design <- c(
      area(1, 1), area(1, 2), area(1, 3), area(1, 4), area(1, 5), area(1, 6), area(1, 7), area(1, 8), area(1, 9), area(1, 10), # 第一行：10 张图
      area(2, 1, 2, 2), area(2, 3, 2, 4), area(2, 5, 2, 6), area(2, 7, 2, 8), area(2, 9, 2, 10), # 第二行：5 张图
      area(3, 1, 3, 2), area(3, 3, 3, 4), area(3, 5, 3, 6), area(3, 7, 3, 8), area(3, 9, 3, 10), # 第三行：5 张图
      area(4, 1, 4, 2), area(4, 3, 4, 4), area(4, 5, 4, 7), area(4, 8, 4, 10), # 第四行：4 张图
      area(5, 1, 5, 2), area(5, 3, 5, 4)
    )
    
    # 组合所有的图
    combined_plot <- wrap_plots(plots) + plot_layout(design = design)
    
    # 添加标题
    combined_plot <- combined_plot + 
      plot_annotation(
        title = paste("Expression of", Gene, "in Classification"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20))
        )
      )
    
    return(combined_plot)
  }
  

  analyze_summary_singlegene <- function(Gene) {
    # Step 1: 检查基因是否存在于表达矩阵
    if (!(Gene %in% rownames(combined_data))) {
      stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
    }
    
    # Step 2: 获取目标基因的表达值
    gene_expression <- combined_data[Gene, ]
    
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
        lower_quantile = quantile(Expression, 0.01, na.rm = TRUE),  # 计算下1%分位数
        upper_quantile = quantile(Expression, 0.99, na.rm = TRUE)   # 计算上1%分位数
      ) %>%
      filter(Expression >= lower_quantile & Expression <= upper_quantile)  # 保留中间的98%
    
    # Step 8: 统计分析（NC vs. Obese, NC vs. NAFLD, NC vs. NASH）
    stat_test <- stat_compare_means(comparisons = list(c("NC", "Obese"),
                                                       c("NC", "NAFLD"),
                                                       c("NC", "NASH"),
                                                       c("Obese", "NAFLD"),
                                                       c("Obese", "NASH"),
                                                       c("NAFLD", "NASH")
    ), 
    method = "t.test", label = "p.signif")
    
    # Step 9: 计算所有选择组的总样本数
    total_samples <- sum(grepl("NC|Obese|NAFLD|NASH", gene_df_clean$Condition))
    
    # 创建图例标签
    non_na_counts <- table(gene_df_clean$Condition)
    legend_labels <- paste(names(non_na_counts), "(n=", non_na_counts, ")", sep = "")
    
    # Step 10: 自定义颜色（动态生成）
    unique_conditions <- levels(gene_df_clean$Condition)
    custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
    custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
    
    # Step 11: 创建标题并添加每个条件组的样本总数
    plot_title <- paste("Expression of ", Gene, 
                        "   (n=", total_samples, ")", sep="")
    
    # Step 12: 绘制箱线图，并在图例中显示每个组的非空样本数
    plot <- ggplot(gene_df_clean, aes(x = Condition, y = Expression, fill = Condition)) +
      geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
      geom_jitter(aes(color = Condition), width = 0.2, size = 0.5, alpha = 0.6) + 
      scale_fill_manual(values = custom_colors_fill, labels = legend_labels) +  
      scale_color_manual(values = custom_colors_color, labels = legend_labels) +  
      stat_test +  
      labs(
        title = plot_title,  # 使用更新后的标题
        x = NULL,
        y = "Expression Level",
        fill = "Condition",   # 图例标题
        color = "Condition"   # 图例标题
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),   # 设置x轴字体大小
        axis.text.y = element_text(size = 12),                            # 设置y轴字体大小
        legend.title = element_text(size = 12),                           # 设置图例标题字体大小
        legend.text = element_text(size = 10),                            # 设置图例标签字体大小
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # 标题居中
        legend.position = "right"
      )
    
    # 返回绘图
    return(plot)
  }
  
  # 监听基因输入并更新图表
  observeEvent(input$update, {
    req(user_auth())  # 确保只有在登录后才能更新
    gene <- input$gene  # 获取用户输入的基因名称
    analysis_type <- input$analysis_type  # 获取分析类型
    
    # 根据分析类型选择不同的绘图函数
    if (analysis_type == "NAFLD Score") {
      combined_plot <- analyze_nafld_singlegene(gene)
      plot_width <- 1400
      plot_height <- 1000
    } else if (analysis_type == "Classification") {
      combined_plot <- analyze_classification_singlegene(gene)
      plot_width <- 1400  # 由于图较多，设置较大尺寸
      plot_height <- 1600
    } else if (analysis_type == "Summary") {
      combined_plot <- analyze_summary_singlegene(gene)
      plot_width <- 700 # 由于图较多，设置较大尺寸
      plot_height <- 530
    } else {
      combined_plot <- analyze_fibrosis_singlegene(gene)
      plot_width <- 1400
      plot_height <- 1000
    }
    # 动态渲染 `plotOutput()`
    output$dynamic_plot <- renderUI({
      plotOutput("combined_plot", width = paste0(plot_width, "px"), height = paste0(plot_height, "px"))
    })
    # 生成图像
    output$combined_plot <- renderPlot({
      combined_plot
    })
    
  })
}

# 运行 Shiny 应用
shinyApp(ui = ui, server = server)
