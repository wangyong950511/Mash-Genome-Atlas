library(sva)
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


# 函数集合 ----------
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


###Fibrosis
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
      legend.position = "right",
      text = element_text(size = 11)
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
      legend.position = "right",
      text = element_text(size = 11)
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
  custom_colors_color <- setNames(rainbow(length(unique_conditions)), unique_conditions)
  
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
      legend.position = "right",
      text = element_text(size = 11)
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
  
  custom_colors_fill <- c("NC" = "white", "NASH" = "white")
  custom_colors_color <- c("NC" = "#98DF8A", "NASH" = "#FFB347")
  
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
      legend.position = "right",
      text = element_text(size = 11)
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
      legend.position = "right",
      text = element_text(size = 11)
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
  
  custom_colors_fill <- c("NC" = "white", "NASH" = "white")
  custom_colors_color <- c("NC" = "#98DF8A", "NASH" = "#FFB347")
  
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
      legend.position = "right",
      text = element_text(size = 11)
    )
  
  return(plot)
}
# GSE147304
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
  
  custom_colors_fill <- c("NC" = "white", "NASH" = "white")
  custom_colors_color <- c("NC" = "#98DF8A", "NASH" = "#FFB347")
  
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
      legend.position = "right",
      text = element_text(size = 11)
    )
  
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
                       choices = c("Classification","NAFLD Score", "Fibrosis Stages" ,"scRNAseq"))  # 根据选择显示子分析选项
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
        plotOutput("combined_plot", width = "100%", height = "1000px")
    )
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
    modify_plot <- function(plot) {
      plot + labs(x = NULL) + 
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    }
    
    plots <- list(
      plot_GSE130970 = analyze_gene_expression_GSE130970(Gene, vsd_nafld_GSE130970, colData_GSE130970, "nafld"),
      plot_GSE89632 = analyze_gene_expression_GSE89632(Gene, Illumina_GSE89632, colData_GSE89632, "nafld"),
      plot_GSE48452 = analyze_gene_expression_GSE48452(Gene, Affymetrix_GSE48452, colData_GSE48452, "nas"),
      plot_GSE135251 = analyze_gene_expression_GSE135251(Gene, vsd_nas_GSE135251, colData_GSE135251, "nas"),
      plot_GSE174478 = analyze_gene_expression_GSE174478(Gene, vsd_nafld_GSE174478, colData_GSE174478, "nafld"),
      plot_GSE185051 = analyze_gene_expression_GSE185051(Gene, vsd_nas_GSE185051, colData_GSE185051, "nas"),
      plot_GSE193066 = analyze_gene_expression_GSE193066(Gene, vsd_nafld_GSE193066, colData_GSE193066, "nafld"),
      plot_GSE193080 = analyze_gene_expression_GSE193080(Gene, vsd_nafld_GSE193080, colData_GSE193080, "nafld"),
      plot_GSE207310 = analyze_gene_expression_GSE207310(Gene, vsd_nafld_GSE207310, colData_GSE207310, "nafld"),
      plot_GSE246221 = analyze_gene_expression_GSE246221(Gene, vsd_nafld_GSE246221, colData_GSE246221, "nafld"),
      plot_GSE225740 = analyze_gene_expression_GSE225740(Gene, vsd_nas_GSE225740, colData_GSE225740, "nas")
    )
    
    # 使用 patchwork 为所有组合图像添加大标题
    combined_plot <- plot_grid(plotlist = plots, ncol = 4) + 
      plot_annotation(title = paste("Expression of the ",Gene, " in NAFLD Score"), 
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 30, face = "bold",margin = margin(b = 20))))
    
    return(combined_plot)
  }
  
  analyze_fibrosis_singlegene <- function(Gene) {
    modify_plot <- function(plot) {
      plot + labs(x = NULL) + 
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    }
    plots <- list(
      plot_GSE130970 = analyze_gene_expression_GSE130970(Gene, vsd_fibrosis_GSE130970, colData_GSE130970, "fibrosis"),
      plot_GSE89632 = analyze_gene_expression_GSE89632(Gene, Illumina_GSE89632, colData_GSE89632, "fibrosis"),
      plot_GSE48452 = analyze_gene_expression_GSE48452(Gene, Affymetrix_GSE48452, colData_GSE48452, "fibrosis"),
      plot_GSE135251 = analyze_gene_expression_GSE135251(Gene, vsd_fibrosis_GSE135251, colData_GSE135251, "fibrosis"),
      plot_GSE174478 = analyze_gene_expression_GSE174478(Gene, vsd_fibrosis_GSE174478, colData_GSE174478, "fibrosis"),
      plot_GSE185051 = analyze_gene_expression_GSE185051(Gene, vsd_fibrosis_GSE185051, colData_GSE185051, "fibrosis"),
      plot_GSE193066 = analyze_gene_expression_GSE193066(Gene, vsd_fibrosis_GSE193066, colData_GSE193066, "fibrosis"),
      plot_GSE193080 = analyze_gene_expression_GSE193080(Gene, vsd_fibrosis_GSE193080, colData_GSE193080, "fibrosis"),
      plot_GSE207310 = analyze_gene_expression_GSE207310(Gene, vsd_fibrosis_GSE207310, colData_GSE207310, "fibrosis"),
      plot_GSE246221 = analyze_gene_expression_GSE246221(Gene, vsd_fibrosis_GSE246221, colData_GSE246221, "fibrosis"),
      plot_GSE225740 = analyze_gene_expression_GSE225740(Gene, vsd_fibrosis_GSE225740, colData_GSE225740, "fibrosis")
    )
    
    # 使用 patchwork 为所有组合图像添加大标题
    combined_plot <- plot_grid(plotlist = plots, ncol = 4) + 
      plot_annotation(title = paste("Expression of the ",Gene, " in Fibrosis Stages"), 
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 30, face = "bold",margin = margin(b = 20))))
    
    return(combined_plot)
  }
  
  analyze_classification_singlegene <- function(Gene) {
    modify_plot <- function(plot) {
      plot + labs(x = NULL) + 
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    }
    
    # 创建一个空列表用于保存图形
    plots <- list()
    
    # 使用 tryCatch 捕获每个数据集的错误
    datasets <- list(
      GSE24807 = list(data = codelink_GSE24807, colData = colData_GSE24807, func = analyze_gene_expression_GSE24807),
      GSE37031 = list(data = Affymetrix_GSE37031, colData = colData_GSE37031, func = analyze_gene_expression_GSE37031),
      GSE106737 = list(data = Affymetrix_GSE106737, colData = colData_GSE106737, func = analyze_gene_expression_GSE106737),
      GSE134146 = list(data = Arraystar_GSE134146, colData = colData_GSE134146, func = analyze_gene_expression_GSE134146),
      GSE159676 = list(data = Affymetrix__GSE159676, colData = colData_GSE159676, func = analyze_gene_expression_GSE159676),
      GSE147304 = list(data = vsd_GSE147304, colData = colData_GSE147304, func = analyze_gene_expression_GSE147304),
      GSE173735 = list(data = vsd_GSE173735, colData = colData_GSE173735, func = analyze_gene_expression_GSE173735)
    )
    
    for (dataset in names(datasets)) {
      tryCatch({
        # 尝试调用分析函数，如果成功则保存图形
        plot <- datasets[[dataset]]$func(Gene, datasets[[dataset]]$data, datasets[[dataset]]$colData)
        plots[[paste0("plot_", dataset)]] <- modify_plot(plot)
      }, error = function(e) {
        # 如果遇到错误，打印错误信息并跳过当前图
        message(paste("跳过", dataset, "数据集，错误:", e$message))
      })
    }
    
    # 使用 patchwork 为所有组合图像添加大标题
    combined_plot <- plot_grid(plotlist = plots, ncol = 3) + 
      plot_annotation(title = paste("Expression of the ", Gene, " in Classification"), 
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20))))
    
    return(combined_plot)
  }
  
  # 监听基因输入并更新图表
  observeEvent(input$update, {
    req(user_auth())  # 确保只有在登录后才能更新
    gene <- input$gene  # 获取用户输入的基因名称
    analysis_type <- input$analysis_type  # 获取分析类型
    
    if (analysis_type == "NAFLD Score") {
      combined_plot <- analyze_nafld_singlegene(gene)
    } else if (analysis_type == "Classification") {
      combined_plot <- analyze_classification_singlegene(gene)  # 替换为你想要的函数
    } else {
      combined_plot <- analyze_fibrosis_singlegene(gene)
    }
    
    output$combined_plot <- renderPlot({
      combined_plot
    }, height = 1000 , width = 1200 )
  })
}

# 运行 Shiny 应用
shinyApp(ui = ui, server = server)
