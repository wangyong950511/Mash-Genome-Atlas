library(sva)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(dplyr)
library(data.table)
library(AnnotationDbi)
library(babelgene)
library(myplotColors)
library(cowplot)
library(myplotColors)
library(babelgene)
library(parallel)
library(limma)
library(DESeq2)
library(shinycssloaders)
library(shinyjs)
library(patchwork)


# 函数集合 ----------

####总结数据
load("~/R/RNA/data/NASH/Summary/Data_Summary_Human.RData")
load("~/R/RNA/data/NASH/Summary/WGCNA_Summary_Human.RData")

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
# 函数集合 
# GSE24807
load("~/R/RNA/data/NASH/GSE24807/Data_GSE24807.RData")
analyze_gene_expression_GSE24807 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 构建分组设计矩阵
  group <- factor(colData$condition, levels = c("NC", "NASH"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # Step 2: 线性模型拟合
  fit <- lmFit(expr_matrix, design)
  contrast.matrix <- makeContrasts(NASH - NC, levels = design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  
  # Step 3: 获取差异表达基因
  results <- topTable(fit, adjust.method = "fdr", number = Inf)
  
  # Step 4: 检查目标基因是否存在
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取目标基因的表达值
  gene_expression <- as.numeric(expr_matrix[GeneX, ])  # 确保为数值向量
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 6: 仅保留 "NC" 和 "NASH" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "NASH"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NASH"))
  
  # Step 8: 获取 logFC 和 p.adj 值
  if (GeneX %in% rownames(results)) {
    logFC_value <- results[GeneX, "logFC"]
    logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  } else {
    logFC_label <- "logFC = NA"
  }
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 11: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE24807"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# GSE37031
load("~/R/RNA/data/NASH/GSE37031/Data_GSE37031.RData")
analyze_gene_expression_GSE37031 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionNASH", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "NASH" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "NASH"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NASH"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE37031"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# GSE106737
load("~/R/RNA/data/NASH/GSE106737/Data_GSE106737.RData")
analyze_gene_expression_GSE106737 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionNASH", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "NASH" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "NASH"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NASH"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE106737"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# GSE134146
load("~/R/RNA/data/NASH/GSE134146/Data_GSE134146.RData")
analyze_gene_expression_GSE134146 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  # 防止空值
  offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionNASH", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "NASH" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "NASH"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NASH"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE134146"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
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
analyze_gene_expression_GSE159676 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionNASH", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "NASH" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "NASH"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NASH"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE159676"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
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
analyze_gene_expression_GSE147304 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "NASH" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "NASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NASH"))
  
  # Step 5: 统计分析（仅 NC vs. NASH）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "NASH", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE147304"),
      x = "Condition",
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    labs(x = NULL)+
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# GSE173735
load("~/R/RNA/data/NASH/GSE173735/Data_GSE173735.RData")
analyze_gene_expression_GSE173735 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "NASH" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "NASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NASH"))
  
  # Step 5: 统计分析（仅 NC vs. NASH）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "NASH", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE173735"),
      x = "Condition",
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    labs(x = NULL)+
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# GSE63067
load("~/R/RNA/data/NASH/GSE63067/Data_GSE63067.RData")
analyze_gene_expression_GSE63067 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionNASH", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "Obese"),c("NC", "NASH"),c("Obese", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "NASH" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "Obese", "NASH"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "Obese", "NASH"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE63067"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# GSE61260
load("~/R/RNA/data/NASH/GSE61260/Data_GSE61260.RData")
analyze_gene_expression_GSE61260 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionNASH", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "Obese"),c("NC", "NAFLD"),c("NC", "NASH"),c("Obese", "NAFLD"),c("Obese", "NASH"),c("NAFLD", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "NASH" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "Obese", "NAFLD", "NASH"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "Obese", "NAFLD", "NASH"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE61260"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE66676/Data_GSE66676.RData")
analyze_gene_expression_GSE66676 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionDefinite_NASH", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NAFLD"),c("NC", "Borderline_NASH"),c("NC", "Definite_NASH"),c("NAFLD", "Borderline_NASH"),c("NAFLD", "Definite_NASH"),c("Borderline_NASH", "Definite_NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "NASH" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "NAFLD", "Borderline_NASH", "Definite_NASH"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NAFLD", "Borderline_NASH", "Definite_NASH"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE66676"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

#
load("~/R/RNA/data/NASH/GSE89632/Data_GSE89632_classification.RData")
analyze_gene_expression_GSE89632_classification <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionNASH", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "Obese"),c("NC", "NASH"),c("Obese", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "NASH" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "Obese", "NASH"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "Obese", "NASH"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE89632"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE48452/Data_GSE48452_classification.RData")
analyze_gene_expression_GSE48452_classification <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionNASH", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "Obese"),c("NC", "Steatosis"),c("NC", "NASH"),c("Obese", "Steatosis"),c("Obese", "NASH"),c("Steatosis", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "NASH" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "Obese", "Steatosis", "NASH"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "Obese", "Steatosis", "NASH"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE48452"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
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
analyze_gene_expression_GSE83452 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC_baseline` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC_baseline")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionNASH_baseline", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC_baseline", "NASH_baseline")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC_baseline" 和 "NASH_baseline" 组
  gene_df <- subset(gene_df, Condition %in% c("NC_baseline", "NASH_baseline"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC_baseline", "NASH_baseline"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE83452"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE105127/Data_GSE105127.RData")
analyze_gene_expression_GSE105127 <- function(GeneX, dds, colData) {
  vsd <- vst(dds, blind = FALSE)
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
      x = NULL,
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
analyze_gene_expression_GSE115193 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "NASH" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "NAFLD", "NASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NAFLD", "NASH"))
  
  # Step 5: 统计分析（仅 NC vs. NASH）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NAFLD"),c("NC", "NASH"),c("NAFLD", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "NASH", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE115193"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
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
analyze_gene_expression_GSE126848 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "NASH" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "Obese", "NAFLD", "NASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "Obese", "NAFLD", "NASH"))
  
  # Step 5: 统计分析（仅 NC vs. NASH）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "Obese"), 
                                                     c("NC", "NAFLD"), 
                                                     c("NC", "NASH"), 
                                                     c("Obese", "NASH"), 
                                                     c("Obese", "NAFLD"), 
                                                     c("NASH", "NAFLD")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "NASH", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE126848"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
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
load("~/R/RNA/data/NASH/GSE135251/Data_GSE135251_classification.RData")
analyze_gene_expression_GSE135251_classification <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "NASH" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "NAFL", "NASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NAFL", "NASH"))
  
  # Step 5: 统计分析（仅 NC vs. NASH）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NAFL"),c("NC", "NASH"),c("NAFL", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "NASH", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE135251"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
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
load("~/R/RNA/data/NASH/GSE164760/Data_GSE164760.RData")
analyze_gene_expression_GSE164760 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionNASH", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NASH"),c("NC", "Cirrhotic"),c("NC", "NASH"),c("NC", "NASH_HCC_notumor"),c("NC", "NASH_HCC_tumor")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "NASH" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "NASH", "Cirrhotic", "NASH_HCC_notumor", "NASH_HCC_tumor"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NASH", "Cirrhotic", "NASH_HCC_notumor", "NASH_HCC_tumor"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE164760"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE167523/Data_GSE167523.RData")
analyze_gene_expression_GSE167523 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NAFL" 和 "NASH" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NAFL", "NASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NAFL", "NASH"))
  
  # Step 5: 统计分析（仅 NAFL vs. NASH）
  stat_test <- stat_compare_means(comparisons = list(c("NAFL", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "NASH", "NAFL"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE167523"),
      x = "Condition",
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    labs(x = NULL)+
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 
load("~/R/RNA/data/NASH/GSE185051/Data_GSE185051_classification.RData")
analyze_gene_expression_GSE185051_classification <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "NASH" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "Steatosis", "Borderline", "NASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "Steatosis", "Borderline", "NASH"))
  
  # Step 5: 统计分析（仅 NC vs. NASH）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "Steatosis"),c("NC", "Borderline"),c("NC", "NASH"),c("Steatosis", "NASH"),c("Steatosis", "Borderline"),c("Borderline", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "NASH", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE185051"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE207310/Data_GSE207310_classification.RData")
analyze_gene_expression_GSE207310_classification <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "NASH" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "NAFL", "NASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NAFL", "NASH"))
  
  # Step 5: 统计分析（仅 NC vs. NASH）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NAFL"),c("NC", "NASH"),c("NAFL", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "NASH", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE207310"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE260666/Data_GSE260666.RData")
analyze_gene_expression_GSE260666 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "NASH" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "NAFLD", "NASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NAFLD", "NASH"))
  
  # Step 5: 统计分析（仅 NC vs. NASH）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NAFLD"),c("NC", "NASH"),c("NAFLD", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "NASH", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE260666"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
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
analyze_gene_expression_GSE274114 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "NASH" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "NASH","ENEG", "ENEG_NASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "ENEG", "NASH","ENEG_NASH"))
  
  # Step 5: 统计分析（仅 NC vs. NASH）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "ENEG"),c("NASH","ENEG_NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "NASH", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE274114"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
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
load("~/R/RNA/data/NASH/GSE288077/Human/Data_GSE288077.RData")
analyze_gene_expression_GSE288077_Human <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "CDAHFD_7w" 组
  gene_df <- gene_df[gene_df$Condition %in% c("ENEG", "ENEG_NASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("ENEG", "ENEG_NASH"))
  
  # Step 5: 统计分析（仅 NC vs. CDAHFD_7w）
  stat_test <- stat_compare_means(comparisons = list(c("ENEG", "ENEG_NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "ENEG_NASH", "ENEG"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE288077"),
      x = "Condition",
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    labs(x = NULL)+
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE59045/Data_GSE59045.RData")
analyze_gene_expression_GSE59045 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionNASH", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "NAFLD"),c("NC", "NASH"),c("NAFLD", "NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "NASH" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "NAFLD", "NASH"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "NAFLD", "NASH"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE59045"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE49541/Data_GSE49541.RData")
analyze_gene_expression_GSE49541 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `Fibrosis_0_1` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "Fibrosis_0_1")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionFibrosis_3_4", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("Fibrosis_0_1", "Fibrosis_3_4")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "Fibrosis_0_1" 和 "Fibrosis_3_4" 组
  gene_df <- subset(gene_df, Condition %in% c("Fibrosis_0_1", "Fibrosis_3_4"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("Fibrosis_0_1", "Fibrosis_3_4"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE49541"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}





####mouse
##STAM
# 加载数据
load("~/R/RNA/data/NASH/GSE236832/Data_GSE236832.RData")
analyze_gene_expression_GSE236832 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "STAM_4W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "STAM_4W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "STAM_4W"))
  
  # Step 5: 统计分析（仅 NC vs. STAM_4W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "STAM_4W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "STAM_4W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE236832"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE210517/Data_GSE210517.RData")
analyze_gene_expression_GSE210517 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  # 先进行 log2 变换
  expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionSTAM", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "STAM")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "STAM" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "STAM"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "STAM"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE210517"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE114261/Data_GSE114261.RData")
analyze_gene_expression_GSE114261 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "STAM_8W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "STAM_8W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "STAM_8W"))
  
  # Step 5: 统计分析（仅 NC vs. STAM_8W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "STAM_8W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "STAM_8W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE114261"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE83596/Data_GSE83596.RData")
analyze_gene_expression_GSE83596 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  # log2 变换
  expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionSTAM_20w_NonTumor", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "STAM_6w"),c("NC", "STAM_8w"),c("NC", "STAM_12w"),c("NC", "STAM_20w_NonTumor"),
                                                     c("NC", "STAM_20w_Tumor")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "STAM_20w_NonTumor" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "STAM_6w", "STAM_8w", "STAM_12w", "STAM_20w_NonTumor", "STAM_20w_Tumor"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "STAM_6w", "STAM_8w", "STAM_12w", "STAM_20w_NonTumor", "STAM_20w_Tumor"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE83596"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE246221/Data_GSE246221_Mouse.RData")
analyze_gene_expression_GSE246221_Mouse <- function(GeneX, dds, colData) {
  
  vsd <- vst(dds, blind = FALSE)
  
  # 直接指定 selected_conditions
  selected_conditions <- c("SCDonly20w", "STZSCD08w", "STZSCD20w", "HFDonly20w", "STZHFD14w", 
                           "STZHFD20w", "STZHFD32w", "STZHFD44w", "STZHFD50w", "STZHFD56w")
  
  # 提取 ENSEMBL
  gene_id <- GeneX
  
  # Step 2: 检查基因是否存在于表达矩阵
  if (!(gene_id %in% rownames(vsd))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 3: 获取目标基因的表达值
  gene_expression <- assay(vsd)[gene_id, ]
  gene_df <- data.frame(
    Expression = gene_expression,
    Condition = colData$condition
  )
  
  # Step 4: 筛选特定分组
  gene_df <- gene_df[gene_df$Condition %in% selected_conditions, ]
  gene_df$Condition <- factor(gene_df$Condition, levels = selected_conditions)
  
  # Step 9: 自定义颜色
  custom_colors_fill <- setNames(rep("white", length(gene_df$Condition)), gene_df$Condition)
  custom_colors_color <- setNames(rainbow(length(gene_df$Condition)), gene_df$Condition)
  
  # Step 5: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  # 边框颜色根据 Condition 设置
    scale_fill_manual(values = custom_colors_fill) +
    scale_color_manual(values = custom_colors_color) +
    labs(
      title = paste("GSE246221"),
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

# 加载数据
load("~/R/RNA/data/NASH/GSE162863/Data_STAM_GSE162863.RData")
analyze_gene_expression_STAM_GSE162863 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "STAM_10W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "STAM_4W", "STAM_10W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "STAM_4W", "STAM_10W"))
  
  # Step 5: 统计分析（仅 NC vs. STAM_10W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "STAM_4W"),c("NC", "STAM_10W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "STAM_10W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE162863"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

##MCD
# 加载数据
load("~/R/RNA/data/NASH/GSE156918/Data_GSE156918.RData")
analyze_gene_expression_GSE156918 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "MCD_12W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "MCD_12W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "MCD_12W"))
  
  # Step 5: 统计分析（仅 NC vs. MCD_12W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "MCD_12W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "MCD_12W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE156918"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE144443/Data_GSE144443.RData")
analyze_gene_expression_GSE144443 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "MCD_21d" 组
  gene_df <- gene_df[gene_df$Condition %in% c("MCD_3d", "MCD_6d", "MCD_12d", "MCD_21d"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("MCD_3d", "MCD_6d", "MCD_12d", "MCD_21d"))
  
  # Step 5: 统计分析（仅 NC vs. MCD_21d）
  stat_test <- stat_compare_means(comparisons = list(c("MCD_3d", "MCD_6d"),c("MCD_3d", "MCD_12d"),c("MCD_3d", "MCD_21d")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "MCD_21d", "MCD_3d"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE144443"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE274949/Data_GSE274949.RData")
analyze_gene_expression_GSE274949 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "MCD_4W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "MCD_4W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "MCD_4W"))
  
  # Step 5: 统计分析（仅 NC vs. MCD_4W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "MCD_4W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "MCD_4W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE274949"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE273291/Data_GSE273291.RData")
analyze_gene_expression_GSE273291 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "MCD_8W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "MCD_8W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "MCD_8W"))
  
  # Step 5: 统计分析（仅 NC vs. MCD_8W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "MCD_8W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "MCD_8W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE273291"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE218174/Data_GSE218174.RData")
analyze_gene_expression_GSE218174 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "MCD_5W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "MCD_5W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "MCD_5W"))
  
  # Step 5: 统计分析（仅 NC vs. MCD_5W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "MCD_5W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "MCD_5W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE218174"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE165174/Data_GSE165174.RData")
analyze_gene_expression_GSE165174 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "MCD_4W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "MCD_4W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "MCD_4W"))
  
  # Step 5: 统计分析（仅 NC vs. MCD_4W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "MCD_4W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "MCD_4W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE165174"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE113843/Data_GSE113843.RData")
analyze_gene_expression_GSE113843 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1.5  # 计算偏移量，使最小值变为 1
  expr_matrix <- log2(expr_matrix + offset)
  # 先进行 log2 变换
  expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionMCD_8W", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "MCD_8W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "MCD_8W" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "MCD_8W"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "MCD_8W"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE113843"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE205974/Data_GSE205974.RData")
analyze_gene_expression_GSE205974 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "MCD_4W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "MCD_4W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "MCD_4W"))
  
  # Step 5: 统计分析（仅 NC vs. MCD_4W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "MCD_4W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "MCD_4W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE205974"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE162863/Data_GSE162863.RData")
analyze_gene_expression_GSE162863 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "MCD_5W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "MCD_5W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "MCD_5W"))
  
  # Step 5: 统计分析（仅 NC vs. MCD_5W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "MCD_5W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "MCD_5W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE162863"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

###CDAHFD
# 加载数据
load("~/R/RNA/data/NASH/GSE269493/Data_GSE269493.RData")
analyze_gene_expression_GSE269493 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "CDAHFD_20w" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "CDAHFD_3W", "CDAHFD_6W", "CDAHFD_12W", "CDAHFD_20W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "CDAHFD_3W", "CDAHFD_6W", "CDAHFD_12W", "CDAHFD_20W"))
  
  # Step 5: 统计分析（仅 NC vs. CDAHFD_20w）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "CDAHFD_3W"),c("NC", "CDAHFD_6W"),c("NC", "CDAHFD_12W"),c("NC", "CDAHFD_20W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "CDAHFD_20W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE269493"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 获取数据
load("~/R/RNA/data/NASH/GSE263975/Data_GSE263975.RData")
analyze_gene_expression_GSE263975 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "CDAHFD_7w" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "CDAHFD_7w"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "CDAHFD_7w"))
  
  # Step 5: 统计分析（仅 NC vs. CDAHFD_7w）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "CDAHFD_7w")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "CDAHFD_7w", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE263975"),
      x = "Condition",
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    labs(x = NULL)+
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE239861/Data_GSE239861.RData")
analyze_gene_expression_GSE239861 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "CDAHFD_7w" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "CDAHFD_8W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "CDAHFD_8W"))
  
  # Step 5: 统计分析（仅 NC vs. CDAHFD_7w）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "CDAHFD_8W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "CDAHFD_8W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE239861"),
      x = "Condition",
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    labs(x = NULL)+
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE233051/Data_GSE233051.RData")
analyze_gene_expression_GSE233051 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "CDAHFD_12W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "CDAHFD_12W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "CDAHFD_12W"))
  
  # Step 5: 统计分析（仅 NC vs. CDAHFD_12W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "CDAHFD_12W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "CDAHFD_12W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE233051"),
      x = "Condition",
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    labs(x = NULL)+
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE207856/Data_GSE207856.RData")
analyze_gene_expression_GSE207856 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "CDAHFD_16W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "CDAHFD_8W", "CDAHFD_16W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "CDAHFD_8W", "CDAHFD_16W"))
  
  # Step 5: 统计分析（仅 NC vs. CDAHFD_16W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "CDAHFD_8W"),c("NC", "CDAHFD_16W"),c("CDAHFD_8W", "CDAHFD_16W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "CDAHFD_16W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE207856"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE200409/Data_GSE200409.RData")
analyze_gene_expression_GSE200409 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  expr_matrix <- log2(expr_matrix + 1)
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionCDAHFD_6W", number = Inf, adjust.method = "BH")
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "CDAHFD_1W"),c("NC", "CDAHFD_2W"),c("NC", "CDAHFD_6W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "CDAHFD_12W" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "CDAHFD_1W", "CDAHFD_2W", "CDAHFD_6W"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "CDAHFD_1W", "CDAHFD_2W", "CDAHFD_6W"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE200409"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}

# 加载数据
load("~/R/RNA/data/NASH/GSE197322/Data_GSE197322.RData")
analyze_gene_expression_GSE197322 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "CDAHFD_12W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "CDAHFD_12W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "CDAHFD_12W"))
  
  # Step 5: 统计分析（仅 NC vs. CDAHFD_12W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "CDAHFD_12W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "CDAHFD_12W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE197322"),
      x = "Condition",
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    labs(x = NULL)+
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# GSE174543
load("~/R/RNA/data/NASH/GSE174543/Data_GSE174543.RData")
analyze_gene_expression_GSE174543 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "CDAHFD_16W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "CDAHFD_16W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "CDAHFD_16W"))
  
  # Step 5: 统计分析（仅 NC vs. CDAHFD_16W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "CDAHFD_16W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "CDAHFD_16W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE174543"),
      x = "Condition",
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    labs(x = NULL)+
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE154892/Data_GSE154892.RData")
analyze_gene_expression_GSE154892 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionCDAHFD_12W", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "CDAHFD_12W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "CDAHFD_12W" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "CDAHFD_12W"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "CDAHFD_12W"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE154892"),
      x = "Condition",
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE120977/Data_GSE120977.RData")
analyze_gene_expression_GSE120977 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "CDAHFD_12W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "CDAHFD_12W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "CDAHFD_12W"))
  
  # Step 5: 统计分析（仅 NC vs. CDAHFD_12W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "CDAHFD_12W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "CDAHFD_12W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE120977"),
      x = "Condition",
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    labs(x = NULL)+
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}



# 加载数据
load("~/R/RNA/data/NASH/GSE35961/Data_GSE35961.RData")
analyze_gene_expression_GSE35961 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionCDAHFD_8W", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "CDAHFD_8W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "CDAHFD_8W" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "CDAHFD_8W"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "CDAHFD_8W"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE35961"),
      x = "Condition",
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}


###HFD
# 加载数据
load("~/R/RNA/data/NASH/GSE94790/Data_GSE94790.RData")
analyze_gene_expression_GSE94790 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionHFD_17W", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFD_17W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "HFD_17W" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "HFD_17W"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFD_17W"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE94790"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}


# 获取数据
load("~/R/RNA/data/NASH/GSE165855/Data_GSE165855.RData")
analyze_gene_expression_GSE165855 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFD_18W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFD_18W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFD_18W"))
  
  # Step 5: 统计分析（仅 NC vs. HFD_18W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFD_18W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFD_18W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE165855"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE242881/Data_GSE242881.RData")
analyze_gene_expression_GSE242881 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFD_26W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFD_26W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFD_26W"))
  
  # Step 5: 统计分析（仅 NC vs. HFD_26W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFD_26W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFD_26W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE242881"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE189066/Data_GSE189066.RData")
analyze_gene_expression_GSE189066 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFD_12W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFD_12W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFD_12W"))
  
  # Step 5: 统计分析（仅 NC vs. HFD_12W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFD_12W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFD_12W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE189066"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE205846/Data_GSE205846.RData")
analyze_gene_expression_GSE205846 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFD_15W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFD_15W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFD_15W"))
  
  # Step 5: 统计分析（仅 NC vs. HFD_15W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFD_15W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFD_15W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE205846"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE218025/Data_GSE218025.RData")
analyze_gene_expression_GSE218025 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFD_16W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFD_16W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFD_16W"))
  
  # Step 5: 统计分析（仅 NC vs. HFD_16W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFD_16W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFD_16W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE218025"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE186116/Data_GSE186116.RData")
analyze_gene_expression_GSE186116 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFD_16W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFD_16W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFD_16W"))
  
  # Step 5: 统计分析（仅 NC vs. HFD_16W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFD_16W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFD_16W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE186116"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE126790/Data_GSE126790.RData")
analyze_gene_expression_GSE126790 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "HFD_21W_Steatosis" 和 "HFD_21W_NASH" 组
  gene_df <- gene_df[gene_df$Condition %in% c("HFD_21W_Steatosis", "HFD_21W_NASH"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("HFD_21W_Steatosis", "HFD_21W_NASH"))
  
  # Step 5: 统计分析（仅 HFD_21W_Steatosis vs. HFD_21W_NASH）
  stat_test <- stat_compare_means(comparisons = list(c("HFD_21W_Steatosis", "HFD_21W_NASH")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFD_21W_NASH", "HFD_21W_Steatosis"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE126790"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE59042/Data_GSE59042.RData")
analyze_gene_expression_GSE59042 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC_36W")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionHFD_36W", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  # stat_test <- stat_compare_means(comparisons = list(c("NC", "HFD_36W")), 
  #                                 method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "HFD_36W" 组
  gene_df <- subset(gene_df, Condition %in% c("NC_12W", "NC_36W", "HFD_12W", "HFD_36W"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC_12W", "NC_36W", "HFD_12W", "HFD_36W"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    # stat_test +  
    labs(
      title = paste("GSE59042"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE145665/Data_GSE145665.RData")
analyze_gene_expression_GSE145665 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFD_24W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFD_24W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFD_24W"))
  
  # Step 5: 统计分析（仅 NC vs. HFD_24W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFD_24W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFD_24W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE145665"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE129306/Data_GSE129306.RData")
analyze_gene_expression_GSE129306 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFD_24W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFD_24W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFD_24W"))
  
  # Step 5: 统计分析（仅 NC vs. HFD_24W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFD_24W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFD_24W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE129306"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE24031/Data_GSE24031.RData")
analyze_gene_expression_GSE24031 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `LFD_L` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "LFD_21W_L")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionHFD_21W_H", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("LFD_21W_L", "LFD_21W_H"),
                                                     c("HFD_21W_L", "HFD_21W_H")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "LFD_L" 和 "HFD_H" 组
  gene_df <- subset(gene_df, Condition %in% c("LFD_21W_L", "LFD_21W_H","HFD_21W_L", "HFD_21W_H"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("LFD_21W_L", "LFD_21W_H","HFD_21W_L", "HFD_21W_H"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE24031"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE43106/Data_GSE43106.RData")
analyze_gene_expression_GSE43106 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `B6J_NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "B6J_NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionB6J_HFD_3W", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("129_NC", "129_HFD_3W"),
                                                     c("B6J_NC", "B6J_HFD_3W"),
                                                     c("B6N_NC","B6N_HFD_3W"),
                                                     c("C3H_NC","C3H_HFD_3W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "B6J_NC" 和 "B6J_HFD_3W" 组
  gene_df <- subset(gene_df, Condition %in% c("129_NC", "129_HFD_3W", "B6J_NC", "B6J_HFD_3W",
                                              "B6N_NC","B6N_HFD_3W","C3H_NC","C3H_HFD_3W"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("129_NC", "129_HFD_3W", "B6J_NC", "B6J_HFD_3W",
                                                            "B6N_NC","B6N_HFD_3W","C3H_NC","C3H_HFD_3W"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE43106"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE179394/Data_GSE179394.RData")
analyze_gene_expression_GSE179394 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFD_50W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFD_50W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFD_50W"))
  
  # Step 5: 统计分析（仅 NC vs. HFD_50W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFD_50W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFD_50W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE179394"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE288077/Data_GSE288077_Mouse.RData")
analyze_gene_expression_GSE288077_Mouse <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFD_16W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFD_16W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFD_16W"))
  
  # Step 5: 统计分析（仅 NC vs. HFD_16W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFD_16W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFD_16W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE288077"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE109345/Data_GSE109345.RData")
analyze_gene_expression_GSE109345 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFD_24W" 组
  # gene_df <- gene_df[gene_df$Condition %in% c("NC_24W", "HFD_24W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC_00W","NC_06W", "HFD_06W","NC_12W", "HFD_12W", "NC_18W", "HFD_18W", "NC_24W", "HFD_24W"))
  
  # Step 5: 统计分析（仅 NC vs. HFD_24W）
  stat_test <- stat_compare_means(comparisons = list(c("NC_06W", "HFD_06W"),c("NC_12W", "HFD_12W"),c("NC_18W", "HFD_18W"),c("NC_24W", "HFD_24W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFD_24W", "NC_24W"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE109345"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE201819/Data_GSE201819.RData")
analyze_gene_expression_GSE201819 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC_C57BL","WD_17W_C57BL","NC_129S1","WD_17W_129S1","NC_A_J","WD_17W_A_J","NC_CAST","WD_17W_CAST","NC_DBA","WD_17W_DBA","NC_PWK","WD_17W_PWK", "NC_WSB",
                                                            "WD_17W_WSB"))
  
  # Step 5: 统计分析（仅 NC_C57BL vs. WD_17W_C57BL）
  stat_test <- stat_compare_means(comparisons = list(c("NC_C57BL", "WD_17W_C57BL"),c("NC_129S1","WD_17W_129S1"),c("NC_A_J","WD_17W_A_J"),c("NC_CAST","WD_17W_CAST"),c("NC_DBA","WD_17W_DBA"),c("NC_PWK","WD_17W_PWK"),c("NC_WSB",
                                                                                                                                                                                                                        "WD_17W_WSB")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "WD_17W_C57BL", "NC_C57BL"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE201819"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 获取数据
load("~/R/RNA/data/NASH/GSE195798/Data_GSE195798.RData")
analyze_gene_expression_GSE195798 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "FFD_38W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "FFD_38W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "FFD_38W"))
  
  # Step 5: 统计分析（仅 NC vs. FFD_38W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "FFD_38W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "FFD_38W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE195798"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE226496/Data_GSE226496.RData")
analyze_gene_expression_GSE226496 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "FFD_37W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "FFD_25W", "FFD_37W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "FFD_25W", "FFD_37W"))
  
  # Step 5: 统计分析（仅 NC vs. FFD_37W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "FFD_25W"),c("NC", "FFD_37W"),c("FFD_25W", "FFD_37W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "FFD_37W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE226496"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


###HFHC
# 加载数据
load("~/R/RNA/data/NASH/GSE53381/Data_GSE53381.RData")
analyze_gene_expression_GSE53381 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionHFHC_4W", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHC_4W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "HFHC_4W" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "HFHC_4W"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHC_4W"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE53381"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE229188/Data_GSE229188.RData")
analyze_gene_expression_GSE229188 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHC_32W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHC_32W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHC_32W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHC_32W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHC_32W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHC_32W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE229188"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE51432/Data_GSE51432.RData")
analyze_gene_expression_GSE51432 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionHFHC_12W", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHC_12W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "HFHC_12W" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "HFHC_12W"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHC_12W"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE51432"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}



###HFHS
# 加载数据
load("~/R/RNA/data/NASH/GSE154426/Data_GSE154426.RData")
analyze_gene_expression_GSE154426 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHS_8W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHS_8W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHS_8W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHS_8W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHS_8W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHS_8W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE154426"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE67678/Data_GSE67678.RData")
analyze_gene_expression_GSE67678 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  expr_matrix <- log2(expr_matrix + 1)
  expr_matrix <- normalizeBetweenArrays(expr_matrix, method = "quantile")
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionWDF_8W", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "WDF_8W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "WDF_8W" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "WDF_8W"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "WDF_8W"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE67678"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE94593/Data_GSE94593.RData")
analyze_gene_expression_GSE94593 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "WDF_14W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "WDF_14W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "WDF_14W"))
  
  # Step 5: 统计分析（仅 NC vs. WDF_14W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "WDF_14W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "WDF_14W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE94593"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE197884/Data_GSE197884.RData")
analyze_gene_expression_GSE197884 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "WDF_54W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "WDF_54W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "WDF_54W"))
  
  # Step 5: 统计分析（仅 NC vs. WDF_54W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "WDF_54W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "WDF_54W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE197884"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE215225/Data_GSE215225.RData")
analyze_gene_expression_GSE215225 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "WDF_16W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "WDF_16W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "WDF_16W"))
  
  # Step 5: 统计分析（仅 NC vs. WDF_16W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "WDF_16W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "WDF_16W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE215225"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE162211/Data_GSE162211.RData")
analyze_gene_expression_GSE162211 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "WDF_24W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "WDF_24W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "WDF_24W"))
  
  # Step 5: 统计分析（仅 NC vs. WDF_24W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "WDF_24W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "WDF_24W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE162211"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE253217/Data_GSE253217.RData")
analyze_gene_expression_GSE253217 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "WDF_14W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "WDF_14W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "WDF_14W"))
  
  # Step 5: 统计分析（仅 NC vs. WDF_14W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "WDF_14W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "WDF_14W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE253217"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE273290/Data_GSE273290.RData")
analyze_gene_expression_GSE273290 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHS_16W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHS_16W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHS_16W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHS_16W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHS_16W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHS_16W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE273290"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE230745/Data_GSE230745.RData")
analyze_gene_expression_GSE230745 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "WDF_24W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "WDF_24W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "WDF_24W"))
  
  # Step 5: 统计分析（仅 NC vs. WDF_24W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "WDF_24W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "WDF_24W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE230745"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE110404/Data_GSE110404.RData")
analyze_gene_expression_GSE110404 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHS_32W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHS_32W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHS_32W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHS_32W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHS_32W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHS_32W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE110404"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


###HFHSHC
# 获取数据
load("~/R/RNA/data/NASH/GSE190982/Data_GSE190982.RData")
analyze_gene_expression_GSE190982 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # # Step 3: 仅保留 "LF_24W" 和 "HFHSHC_24W" 组
  # gene_df <- gene_df[gene_df$Condition %in% c("LF_24W", "HFHSHC_24W"), ]
  # 
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("LF_8W","HF_8W","FFHC_8W", "HFHSHC_8W",
                                                            "LF_24W","HF_24W","FFHC_24W", "HFHSHC_24W"))
  
  # Step 5: 统计分析（仅 LF_24W vs. HFHSHC_24W）
  stat_test <- stat_compare_means(comparisons = list(c("LF_24W", "HFHSHC_24W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_24W", "LF_24W"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE190982"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE52748/Data_GSE52748.RData")
analyze_gene_expression_GSE52748 <- function(GeneX, expr_matrix, colData) {
  # Step 1: 确保表达矩阵是数值型
  expr_matrix <- as.matrix(expr_matrix)
  ## 防止空值
  # offset <- abs(min(expr_matrix, na.rm = TRUE)) + 1  # 计算偏移量，使最小值变为 1
  # expr_matrix <- log2(expr_matrix + offset)
  ## log2 变换
  # expr_matrix <- log2(expr_matrix + 1)
  
  # Step 2: 确保 `NC` 为基准组
  colData$condition <- relevel(factor(colData$condition), ref = "NC")
  
  # Step 3: 运行 limma 进行差异表达分析
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "conditionHFHSHC_12W", number = Inf, adjust.method = "BH")
  
  # Step 4: 检查基因是否在表达矩阵中
  if (!(GeneX %in% rownames(expr_matrix))) {
    stop("目标基因不存在于表达矩阵中，请检查数据是否包含该基因。")
  }
  
  # Step 5: 获取基因表达值
  gene_expression <- expr_matrix[GeneX, ]
  gene_df <- data.frame(
    Sample = rownames(colData),
    Expression = gene_expression,
    Condition = colData$condition,
    stringsAsFactors = FALSE
  )
  
  # Step 5: 统计分析（t-test）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_12W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 仅保留 "NC" 和 "HFHSHC_12W" 组
  gene_df <- subset(gene_df, Condition %in% c("NC", "HFHSHC_12W"))
  
  # Step 7: 设置 Condition 因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_12W"))
  
  # Step 8: 获取 logFC 和 padj 值
  logFC_value <- results[GeneX, "logFC"]
  
  # 生成 logFC 显示文本
  logFC_label <- paste0("logFC = ", round(logFC_value, 3))
  
  # Step 9: 颜色设置
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 10: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE52748"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(plot)
}


# 获取数据
load("~/R/RNA/data/NASH/GSE168069/Data_GSE168069.RData")
analyze_gene_expression_GSE168069 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_32W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHSHC_32W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_32W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_32W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_32W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_32W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE168069"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE200482/Data_GSE200482.RData")
analyze_gene_expression_GSE200482 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_16W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHSHC_16W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_16W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_16W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_16W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_16W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE200482"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 获取数据
load("~/R/RNA/data/NASH/GSE133566/Data_GSE133566.RData")
analyze_gene_expression_GSE133566 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_12W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHSHC_12W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_12W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_12W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_12W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_12W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE133566"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE160292/Data_GSE160292.RData")
analyze_gene_expression_GSE160292 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_18W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHSHC_18W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_18W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_18W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_18W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_18W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE160292"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 获取数据
load("~/R/RNA/data/NASH/GSE195483/Data_GSE195483.RData")
analyze_gene_expression_GSE195483 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_36W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHSHC_36W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_36W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_36W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_36W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_36W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE195483"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE196908/Data_GSE196908.RData")
analyze_gene_expression_GSE196908 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_12W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC","HFHSHC_8W", "HFHSHC_12W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC","HFHSHC_8W", "HFHSHC_12W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_12W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_8W"),c("NC", "HFHSHC_12W"),c("HFHSHC_8W", "HFHSHC_12W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_12W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE196908"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE162869/Data_GSE162869.RData")
analyze_gene_expression_GSE162869 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "FPC_24W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "FPC_24W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "FPC_24W"))
  
  # Step 5: 统计分析（仅 NC vs. FPC_24W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "FPC_24W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "FPC_24W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE162869"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE256501/Data_GSE256501.RData")
analyze_gene_expression_GSE256501 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_16W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHSHC_16W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_16W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_16W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_16W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_16W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE256501"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE126204/Data_GSE126204.RData")
analyze_gene_expression_GSE126204 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_24W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHSHC_24W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_24W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_24W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_24W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_24W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE126204"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}



# 加载数据
load("~/R/RNA/data/NASH/GSE164084/Data_GSE164084.RData")
analyze_gene_expression_GSE164084 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_30W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHSHC_30W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_30W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_30W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_30W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_30W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE164084"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE190140/Data_GSE190140.RData")
analyze_gene_expression_GSE190140 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_21W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHSHC_21W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_21W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_21W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_21W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_21W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE190140"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE220519/Data_GSE220519.RData")
analyze_gene_expression_GSE220519 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_18W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHSHC_18W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_18W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_18W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_18W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_18W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE220519"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 获取数据
load("~/R/RNA/data/NASH/GSE243976/Data_GSE243976.RData")
analyze_gene_expression_GSE243976 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_72W_NonTumor" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHSHC_72W_NonTumor","HFHSHC_72W_Tumor"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_72W_NonTumor","HFHSHC_72W_Tumor"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_72W_NonTumor）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_72W_NonTumor"),c("NC", "HFHSHC_72W_Tumor"),c("HFHSHC_72W_Tumor", "HFHSHC_72W_NonTumor")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_72W_NonTumor", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE243976"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 获取数据
load("~/R/RNA/data/NASH/GSE233767/Data_GSE233767.RData")
analyze_gene_expression_GSE233767 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_28W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHSHC_28W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_28W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_28W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_28W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_28W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE233767"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}


# 加载数据
load("~/R/RNA/data/NASH/GSE256063/Data_GSE256063.RData")
analyze_gene_expression_GSE256063 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "HFHSHC_12W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "HFHSHC_12W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "HFHSHC_12W"))
  
  # Step 5: 统计分析（仅 NC vs. HFHSHC_12W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "HFHSHC_12W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "HFHSHC_12W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE256063"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}



####WDCCL4
# 加载数据
load("~/R/RNA/data/NASH/GSE184256/Data_GSE184256.RData")
analyze_gene_expression_GSE184256 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "WDCCL4_24W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "WDCCL4_12W", "WDCCL4_24W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "WDCCL4_12W", "WDCCL4_24W"))
  
  # Step 5: 统计分析（仅 NC vs. WDCCL4_24W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "WDCCL4_24W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "WDCCL4_24W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE184256"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}

# 分析数据
load("~/R/RNA/data/NASH/GSE287943/Data_GSE287943.RData")
analyze_gene_expression_GSE287943 <- function(GeneX, dds,colData) {
  vsd <- vst(dds, blind = FALSE)
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
  
  # Step 3: 仅保留 "NC" 和 "WDCCL4_12W" 组
  gene_df <- gene_df[gene_df$Condition %in% c("NC", "WDCCL4_12W"), ]
  
  # Step 4: 设置 Condition 的因子顺序
  gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "WDCCL4_12W"))
  
  # Step 5: 统计分析（仅 NC vs. WDCCL4_12W）
  stat_test <- stat_compare_means(comparisons = list(c("NC", "WDCCL4_12W")), 
                                  method = "t.test", label = "p.signif")
  
  # Step 6: 自定义颜色
  unique_conditions <- levels(gene_df$Condition)
  custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
  custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
  
  # Step 7: 获取 logFC 值
  res <- results(dds, contrast = c("condition", "WDCCL4_12W", "NC"))
  logFC_value <- res[GeneX, "log2FoldChange"]
  logFC_label <- paste("logFC =", round(logFC_value, 3))
  
  # Step 8: 绘制箱线图
  plot <- ggplot(gene_df, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, aes(color = Condition)) +  
    scale_fill_manual(values = custom_colors_fill) +  
    scale_color_manual(values = custom_colors_color) +  
    stat_test +  
    labs(
      title = paste("GSE287943"),
      x = NULL,
      y = "Expression Level",
      caption = logFC_label  # 在图的底部添加 logFC 值
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # 返回绘图
  return(plot)
}





# 网页端设置 ----------


ui <- fluidPage(
  shinyjs::useShinyjs(),
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
      
      /* 白色遮罩，只覆盖 .main-content 内容区域 */
      #plot_mask {
        position: absolute;
        top: 0;
        left: 0;
        width: calc(100% - 320px); /* 侧边栏宽度为 300px + margin */
        height: 100%;
        background: white;
        z-index: 10000;
        display: none;
        margin-left: 320px; /* 对齐 main-content 的起点 */
      }

     /* 自定义 titlePanel 样式 */
      #sidebar_panel .title {
        font-size: 18px; 
        white-space: nowrap; 
        overflow: hidden; 
        text-overflow: ellipsis; 
        text-align: center;  
        width: 100%; 
        padding: 10px 0; 
        line-height: 1.4; 
        margin-bottom: 20px; 
      }
    "))
  ),
  
  # 登录面板
  conditionalPanel(
    condition = "!output.user_authenticated",
    div(id = "login_panel",
        h3("用户登录", style = "text-align: center;"),
        textInput("user", "用户名", value = "test", width = "100%"),
        passwordInput("password", "密码", value = "", width = "100%"),
        div(actionButton("login", "登录", class = "btn-primary"), 
            style = "text-align: center; margin-top: 15px;"),
        textOutput("login_status")
    ),
    
    div(
      id = "login_spinner",
      img(src = "https://i.gifer.com/ZKZg.gif", height = "25px"),
      style = "
        position: absolute;
        top: 50%;
        left: calc(100% + 20px);
        transform: translateY(-50%);
        display: none;
      "
    )
  ),
  
  # 登录后的界面
  conditionalPanel(
    condition = "output.user_authenticated",
    div(id = "sidebar_panel",
        div(class = "title", HTML("MGA <br>（Mash Genome Atlas）")),
        
        # 物种选择 (Human / Mouse)
        radioButtons("species", "物种", 
                     choices = c("Human", "Mouse")),  
        
        # 仅当选择 Human 时，显示分析类型
        conditionalPanel(
          condition = "input.species == 'Human'",
          radioButtons("analysis_type", "选择分析类型", 
                       choices = c("Summary", "Classification", "NAFLD Score", "Fibrosis Stages",  "多基因分析"))
        ),
        
        # 仅当选择 Mouse 时，显示造模方式
        conditionalPanel(
          condition = "input.species == 'Mouse'",
          radioButtons("feed_type", "造模方式", 
                       choices = c("HFD", "HFHS", "HFHC", "HFHSHC", "WDCCL4", "STAM", "MCD", "CDAHFD"))
        ),
        
        # "多基因分析" 需要额外的子选项
        conditionalPanel(
          condition = "input.species == 'Human' && input.analysis_type == '多基因分析'",
          radioButtons("multi_gene_analysis", "分析方式", 
                       choices = c("寻找相关基因", "双基因相关性"),
                       selected = "寻找相关基因")  # 默认选择 "寻找相关基因"
        ),
        
        # 基因输入框 (适用于所有分析类型，但“多基因分析”有特殊情况)
        conditionalPanel(
          # 让 Mouse 物种也显示基因输入框
          condition = "(input.species == 'Human' && (input.analysis_type != '多基因分析' || input.multi_gene_analysis == '寻找相关基因')) || input.species == 'Mouse'",
          textInput("gene", "输入基因名称  (标准基因Symbol)", value = "COL1A1")
        ),
        
        # "双基因相关性" -> 显示两个基因输入框
        conditionalPanel(
          condition = "input.species == 'Human' && input.analysis_type == '多基因分析' && input.multi_gene_analysis == '双基因相关性'",
          textInput("gene1", "输入基因1", value = "COL1A1"),
          textInput("gene2", "输入基因2", value = "TP53")
        ),
        
        fluidRow(
          column(6, actionButton("update", "更新图表", class = "btn-primary")),
          column(6, div(id = "spinner_area", img(src = "https://i.gifer.com/ZKZg.gif", height = "25px"), style = "display: none;"))
        )
    ),
    
    div(class = "main-content",
        div(id = "plot_mask"),  # 遮罩层
        uiOutput("dynamic_plot") # 图像输出
    )
  )
)




# 服务器端逻辑
server <- function(input, output, session) {
  #初始化设置
  options(future.globals.maxSize = 3 * 1024^3) # 3GB
  
  # 用户身份验证状态
  user_authenticated <- reactiveVal(FALSE)
  
  # 读取密码
  users_db <- read.csv("~/R/RNA/data/NASH/网页客户端/users.csv", stringsAsFactors = FALSE)
  
  observeEvent(input$login, {
    req(input$user, input$password)
    # 确保用户名和密码在同一行匹配
    user_match <- users_db[users_db$username == input$user & users_db$password == input$password, ]
    
    if (nrow(user_match) == 1) {
      user_authenticated(TRUE)
      # 控制台打印登录成功 + 时间
      login_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      message(paste(login_time,"用户", input$user, "登录成功"))
      output$login_status <- renderText("")  # 清空错误信息
    } else {
      user_authenticated(FALSE)
      output$login_status <- renderText("❌ 用户名或密码错误")
      # 控制台打印登录成功 + 时间
      login_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      message(paste(login_time,"用户", input$user, "登录失败"))
    }
    
  })
  
  output$user_authenticated <- reactive({
    user_authenticated()
  })
  outputOptions(output, "user_authenticated", suspendWhenHidden = FALSE)
  
  # ✅ 物种切换时重置相关选项（这里放非常合适）
  observeEvent(input$species, {
    if (input$species == "Mouse") {
      updateRadioButtons(session, "analysis_type", selected = character(0))
      updateRadioButtons(session, "multi_gene_analysis", selected = "寻找相关基因")
    } else {
      updateRadioButtons(session, "analysis_type", selected = "Summary")
    }
  })
  
  
  # 定义分析函数
  analyze_nafld_singlegene <- function(Gene) {
    
    # 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))+ labs(x = NULL)
    }
    
    # 定义一个占位图（空白图）
    empty_plot <- ggplot() + 
      theme_void() + 
      ggtitle("Data Unavailable")
    
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
        empty_plot  # 失败时返回空白图，保持顺序不变
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
    
    # 定义一个占位图（空白图）
    empty_plot <- ggplot() + 
      theme_void() + 
      ggtitle("Data Unavailable")
    
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
        empty_plot  # 失败时返回空白图，保持顺序不变
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
    
    # ✅ 定义数据集
    datasets <- list(
      GSE24807 = list(data = codelink_GSE24807, colData = colData_GSE24807, func = analyze_gene_expression_GSE24807),
      GSE37031 = list(data = Affymetrix_GSE37031, colData = colData_GSE37031, func = analyze_gene_expression_GSE37031),
      GSE106737 = list(data = Affymetrix_GSE106737, colData = colData_GSE106737, func = analyze_gene_expression_GSE106737),
      GSE134146 = list(data = Arraystar_GSE134146, colData = colData_GSE134146, func = analyze_gene_expression_GSE134146),
      GSE159676 = list(data = Affymetrix_GSE159676, colData = colData_GSE159676, func = analyze_gene_expression_GSE159676),
      GSE147304 = list(data = dds_GSE147304, colData = colData_GSE147304, func = analyze_gene_expression_GSE147304),
      GSE173735 = list(data = dds_GSE173735, colData = colData_GSE173735, func = analyze_gene_expression_GSE173735),
      GSE83452 = list(data = Affymetrix_GSE83452, colData = colData_GSE83452, func = analyze_gene_expression_GSE83452),
      GSE167523 = list(data = dds_GSE167523, colData = colData_GSE167523, func = analyze_gene_expression_GSE167523),
      GSE288077 = list(data = dds_GSE288077_Human, colData = colData_GSE288077_Human, func = analyze_gene_expression_GSE288077_Human),
      GSE49541 = list(data = Affymetrix_GSE49541, colData = colData_GSE49541, func = analyze_gene_expression_GSE49541),
      
      GSE63067 = list(data = Affymetrix_GSE63067, colData = colData_GSE63067, func = analyze_gene_expression_GSE63067),
      GSE89632 = list(data = Illumina_GSE89632_classification, colData = colData_GSE89632_classification, func = analyze_gene_expression_GSE89632_classification),
      GSE260666 = list(data = dds_GSE260666, colData = colData_GSE260666, func = analyze_gene_expression_GSE260666),
      GSE135251 = list(data = dds_GSE135251_classification, colData = colData_GSE135251_classification, func = analyze_gene_expression_GSE135251_classification),
      GSE207310 = list(data = dds_GSE207310_classification, colData = colData_GSE207310_classification, func = analyze_gene_expression_GSE207310_classification),
      GSE59045 = list(data = Affymetrix_GSE59045, colData = colData_GSE59045, func = analyze_gene_expression_GSE59045),
      GSE115193 = list(data = dds_GSE115193, colData = colData_GSE115193, func = analyze_gene_expression_GSE115193),
      
      GSE61260 = list(data = Affymetrix_GSE61260, colData = colData_GSE61260, func = analyze_gene_expression_GSE61260),
      GSE66676 = list(data = Affymetrix_GSE66676, colData = colData_GSE66676, func = analyze_gene_expression_GSE66676),
      GSE48452 = list(data = Affymetrix_GSE48452, colData = colData_GSE48452_classification, func = analyze_gene_expression_GSE48452_classification),
      GSE126848 = list(data = dds_GSE126848, colData = colData_GSE126848, func = analyze_gene_expression_GSE126848),
      
      GSE185051 = list(data = dds_GSE185051_classification, colData = colData_GSE185051_classification, func = analyze_gene_expression_GSE185051_classification),
      GSE274114 = list(data = dds_GSE274114, colData = colData_GSE274114, func = analyze_gene_expression_GSE274114),
      GSE164760 = list(data = affy_GSE164760, colData = colData_GSE164760, func = analyze_gene_expression_GSE164760),
      GSE105127 = list(data = dds_GSE105127, colData = colData_GSE105127, func = analyze_gene_expression_GSE105127)
    )
    
    
    # ✅ 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + labs(x = NULL)
    }
    
    # ✅ 定义一个占位图（空白图）
    empty_plot <- ggplot() + 
      theme_void() + 
      ggtitle("Data Unavailable")
    
    # ✅ 处理所有数据集，捕获可能的错误
    plots <- lapply(names(datasets), function(name) {
      dataset <- datasets[[name]]
      tryCatch({
        # ✅ 判断是否有 group 变量
        if ("group" %in% names(dataset)) {
          plot <- dataset$func(Gene, dataset$data, dataset$colData, dataset$group)
        } else {
          plot <- dataset$func(Gene, dataset$data, dataset$colData)
        }
        modify_plot(plot)  # 统一格式化
      }, error = function(e) {
        message(paste("跳过", name, "数据集，错误:", e$message))
        empty_plot  # 失败时返回空白图，保持顺序不变
      })
    })
    
    # ✅ 定义布局
    design <- c(
      area(1, 1), area(1, 2), area(1, 3), area(1, 4), area(1, 5), area(1, 6), area(1, 7), area(1, 8), area(1, 9),area(1, 10), area(1, 11),# 第一行：图1-4
      area(2, 1, 2,2),  area(2, 3, 2,4),area(2, 5, 2,6),area(2, 7, 2,8),area(2, 9,2,10),area(2, 11,2,12), # 第二行：图5-8
      area(3, 1,3,2),  area(3, 3,3,4),area(3, 5,3,6),area(3, 7, 3,8),area(3, 9,3,10), area(3, 11,3,12),  # 第三行：图9-12
      area(4, 1,4,2),area(4, 3,4,5),area(4, 6, 4,8)
    )
    
    # ✅ 组合所有的图
    combined_plot <- wrap_plots(plots) + plot_layout(design = design)
    
    # ✅ 添加标题
    combined_plot <- combined_plot + 
      plot_annotation(
        title = paste("Expression of", Gene, "in multiple datasets"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20))
        )
      )
    
    return(combined_plot)
  }
  

  analyze_summary_singlegene <- function(Gene) {
  # 使用 tryCatch 处理异常，防止出错时中断
  tryCatch({
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
    
    # Step 4: 只保留特定的条件组
    gene_df <- gene_df[gene_df$Condition %in% c("NC", "Obese", "NAFLD", "NASH"), ]
    
    # Step 5: 设置 Condition 的因子顺序
    gene_df$Condition <- factor(gene_df$Condition, levels = c("NC", "Obese", "NAFLD", "NASH"))
    
    # Step 6: 按因子水平排序数据
    gene_df <- gene_df[order(gene_df$Condition), ]
    
    # Step 7: 去除每组5%的极端值
    gene_df_clean <- gene_df %>%
      group_by(Condition) %>%
      mutate(
        lower_quantile = quantile(Expression, 0.01, na.rm = TRUE),  
        upper_quantile = quantile(Expression, 0.99, na.rm = TRUE)   
      ) %>%
      filter(Expression >= lower_quantile & Expression <= upper_quantile)  
    
    # Step 8: 统计分析
    stat_test <- stat_compare_means(comparisons = list(c("NC", "Obese"),
                                                       c("NC", "NAFLD"),
                                                       c("NC", "NASH"),
                                                       c("Obese", "NAFLD"),
                                                       c("Obese", "NASH"),
                                                       c("NAFLD", "NASH")),
                                    method = "t.test", label = "p.signif")
    
    # Step 9: 计算所有选择组的总样本数
    total_samples <- sum(grepl("NC|Obese|NAFLD|NASH", gene_df_clean$Condition))
    
    # Step 10: 生成图例标签
    non_na_counts <- table(gene_df_clean$Condition)
    legend_labels <- paste(names(non_na_counts), "(n=", non_na_counts, ")", sep = "")
    
    # Step 11: 自定义颜色
    unique_conditions <- levels(gene_df_clean$Condition)
    custom_colors_color <- setNames(getplotColors(length(unique_conditions)), unique_conditions)
    custom_colors_fill <- setNames(rep("white", length(unique_conditions)), unique_conditions)
    
    # Step 12: 创建标题并添加样本数
    plot_title <- paste("Expression of ", Gene, "   (n=", total_samples, ")", sep="")
    
    # Step 13: 绘制箱线图
    plot <- ggplot(gene_df_clean, aes(x = Condition, y = Expression, fill = Condition)) +
      geom_boxplot(outlier.shape = NA, aes(color = Condition)) + 
      geom_jitter(aes(color = Condition), width = 0.2, size = 0.5, alpha = 0.6) + 
      scale_fill_manual(values = custom_colors_fill, labels = legend_labels) +  
      scale_color_manual(values = custom_colors_color, labels = legend_labels) +  
      stat_test +  
      labs(
        title = plot_title,
        x = NULL,
        y = "Expression Level",
        fill = "Condition",
        color = "Condition"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.position = "right"
      )
    
    return(plot)
    
  }, error = function(e) {
    message("Error: ", e$message)  # 打印错误信息
    return(NULL)  # 返回 NULL 以防止程序崩溃
  })
}


  analyze_correlation_doublegene <- function(Gene, Gene2) {
    
    # Step 1: 检查基因是否存在于表达矩阵
    if (!(Gene %in% rownames( WGCNA_data))) {
      stop(paste("目标基因", Gene, "不存在于表达矩阵中，请检查数据是否包含该基因。"))
    }
    if (!(Gene2 %in% rownames( WGCNA_data))) {
      stop(paste("目标基因", Gene2, "不存在于表达矩阵中，请检查数据是否包含该基因。"))
    }
    group="NASH"
    
    # Step 2: 确保 group 是合法的
    valid_groups <- c("NC", "Obese", "NAFLD", "NASH")
    if (!(group %in% valid_groups)) {
      stop("提供的 group 不在有效范围内，请使用 'NC', 'Obese', 'NAFLD', 或 'NASH'。")
    }
    
    # Step 3: 获取目标基因的表达值
    Gene_expression <-  WGCNA_data[Gene, ]
    Gene2_expression <-  WGCNA_data[Gene2, ]
    
    # Step 4: 组合数据并筛选指定 group
    gene_df <- data.frame(
      Gene_Expression = as.numeric(Gene_expression),
      Gene2_Expression = as.numeric(Gene2_expression),
      Condition = combined_colData$condition
    ) %>% 
      filter(Condition == group)  # 只保留指定的 group 数据
    
    # Step 5: 去除 5% 极端值
    gene_df_clean <- gene_df %>%
      mutate(
        lower_quantile_x = quantile(Gene_Expression, 0.01, na.rm = TRUE),  
        upper_quantile_x = quantile(Gene_Expression, 0.99, na.rm = TRUE),  
        lower_quantile_y = quantile(Gene2_Expression, 0.01, na.rm = TRUE),  
        upper_quantile_y = quantile(Gene2_Expression, 0.99, na.rm = TRUE)   
      ) %>%
      filter(Gene_Expression >= lower_quantile_x & Gene_Expression <= upper_quantile_x & 
               Gene2_Expression >= lower_quantile_y & Gene2_Expression <= upper_quantile_y) 
    
    # Step 6: 计算 Spearman 相关性
    cor_value <- cor(gene_df_clean$Gene_Expression, gene_df_clean$Gene2_Expression, method = "spearman", use = "complete.obs")
    p_value <- cor.test(gene_df_clean$Gene_Expression, gene_df_clean$Gene2_Expression, method = "spearman", use = "complete.obs")$p.value
    
    # Step 7: 绘制相关性散点图
    plot <- ggplot(gene_df_clean, aes(x = Gene_Expression, y = Gene2_Expression)) +
      geom_point(alpha = 0.6, size = 2, color = "blue") +  
      geom_smooth(method = "lm", se = FALSE, formula = y ~ x, linetype = "dashed", color = "red", size = 1) +
      labs(
        title = paste("Correlation between", Gene, "and", Gene2, "in", group),
        subtitle = paste("Spearman Correlation:", round(cor_value, 3), " (p =", signif(p_value, 3), ")"),
        x = paste("Expression of", Gene),
        y = paste("Expression of", Gene2)
      ) +
      theme_minimal() +
      theme(
        legend.position = "none"
      )
    
    # Step 8: 输出相关性结果并返回绘图
    message(paste("Group:", group, "| Spearman Correlation:", round(cor_value, 3), "| p-value:", signif(p_value, 3)))
    
    return(plot)
  }
  
  
  analyze_find_correlationGene <- function(target_gene) {
    top_n <- 30  # 选择 top 30 相关基因
    
    # 确保数据是数值型矩阵，并转置（样本为行，基因为列）
    WGCNA_data <- as.matrix(WGCNA_data)  # 确保是矩阵
    WGCNA_data <- t(WGCNA_data)  # 转置，使得基因作为列，样本作为行
    
    # 检查基因是否在数据集中
    if (!(target_gene %in% colnames(WGCNA_data))) {
      stop(paste("目标基因", target_gene, "不存在于表达矩阵中，请检查数据是否包含该基因。"))
    }
    
    # 目标基因表达数据
    target_data <- WGCNA_data[, target_gene]
    
    # **使用 mclapply 进行并行计算**
    num_cores <- detectCores() - 2  # 获取可用核心，保留 2 个给系统
    cor_values <- unlist(mclapply(colnames(WGCNA_data), function(gene) {
      cor(target_data, WGCNA_data[, gene], method = "spearman", use = "pairwise.complete.obs")
    }, mc.cores = num_cores))  # 并行计算相关性
    
    # 转换为数据框
    cor_df <- data.frame(Gene = colnames(WGCNA_data), Correlation = cor_values)
    
    # 按相关性排序（绝对值最大排前面），去掉自身
    top_correlated_genes <- cor_df %>%
      filter(Gene != target_gene) %>%
      arrange(desc(abs(Correlation))) %>%
      head(top_n)
    
    # 绘制相关性条形图
    p <- ggplot(top_correlated_genes, aes(x = reorder(Gene, Correlation), y = Correlation, fill = Correlation)) +
      geom_bar(stat = "identity", width = 0.7) +
      coord_flip() +  # 旋转坐标轴
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      labs(title = paste("Top", top_n, "Genes Correlated with", target_gene),
           x = "Gene", y = "Spearman Correlation") +
      theme_minimal()
    
    return(p)  # 返回绘图对象
  }
  
  
  analyze_Mouse_STAM_singlegene <- function(Gene) {
    
    # 数据集列表
    datasets <- list(
      GSE162863 = list(data = dds_STAM_GSE162863, colData = colData_STAM_GSE162863, func = analyze_gene_expression_STAM_GSE162863),
      GSE114261 = list(data = dds_GSE114261, colData = colData_GSE114261, func = analyze_gene_expression_GSE114261),
      GSE210517 = list(data = Agilent_GSE210517, colData = colData_GSE210517, func = analyze_gene_expression_GSE210517),
      GSE236832 = list(data = dds_GSE236832, colData = colData_GSE236832, func = analyze_gene_expression_GSE236832),
      GSE83596 = list(data = Agilent_GSE83596, colData = colData_GSE83596, func = analyze_gene_expression_GSE83596),
      GSE246221 = list(data = dds_GSE246221_Mouse, colData = colData_GSE246221_Mouse, func = analyze_gene_expression_GSE246221_Mouse)
    )
    
    # 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + labs(x = NULL)
    }
    
    # 定义一个占位图（空白图）
    empty_plot <- ggplot() + 
      theme_void() + 
      ggtitle("Data Unavailable")
    
    # 处理所有数据集，捕获可能的错误
    plots <- lapply(names(datasets), function(name) {
      dataset <- datasets[[name]]
      tryCatch({
        plot <- dataset$func(Gene, dataset$data, dataset$colData)
        modify_plot(plot)  # 统一格式化
      }, error = function(e) {
        message(paste("跳过", name, "数据集，错误:", e$message))
        empty_plot  # 失败时返回空白图，保持顺序不变
      })
    })
    
    # 定义布局
    design <- c(
      area(1, 1), area(1, 2), area(1, 3),area(1, 4) , # 第一行：图1-4
      area(2, 1, 2, 2), area(2, 3, 2, 5)  # 第二行：图5
    )
    
    # 组合所有的图
    combined_plot <- wrap_plots(plots) + plot_layout(design = design)
    
    # 添加标题
    combined_plot <- combined_plot + 
      plot_annotation(
        title = paste("Expression of", Gene, "in STAM mouse model"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20))
        )
      )
    
    return(combined_plot)
  }
  
  
  analyze_Mouse_MCD_singlegene <- function(Gene) {
    
    # 数据集列表
    datasets <- list(
      GSE274949 = list(data = dds_GSE274949, colData = colData_GSE274949, func = analyze_gene_expression_GSE274949),
      GSE165174 = list(data = dds_GSE165174, colData = colData_GSE165174, func = analyze_gene_expression_GSE165174),
      GSE205974 = list(data = dds_GSE205974, colData = colData_GSE205974, func = analyze_gene_expression_GSE205974),
      GSE162863 = list(data = dds_GSE162863, colData = colData_GSE162863, func = analyze_gene_expression_GSE162863),
      GSE218174 = list(data = dds_GSE218174, colData = colData_GSE218174, func = analyze_gene_expression_GSE218174),
      GSE273291 = list(data = dds_GSE273291, colData = colData_GSE273291, func = analyze_gene_expression_GSE273291),
      GSE156918 = list(data = dds_GSE156918, colData = colData_GSE156918, func = analyze_gene_expression_GSE156918),      
      GSE144443 = list(data = dds_GSE144443, colData = colData_GSE144443, func = analyze_gene_expression_GSE144443),
      GSE113843 = list(data = Phalanx_GSE113843, colData = colData_GSE113843, func = analyze_gene_expression_GSE113843)
    )
    
    # 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + labs(x = NULL)
    }
    
    # 定义一个占位图（空白图）
    empty_plot <- ggplot() + 
      theme_void() + 
      ggtitle("Data Unavailable")
    
    # 处理所有数据集，捕获可能的错误
    plots <- lapply(names(datasets), function(name) {
      dataset <- datasets[[name]]
      tryCatch({
        plot <- dataset$func(Gene, dataset$data, dataset$colData)
        modify_plot(plot)  # 统一格式化
      }, error = function(e) {
        message(paste("跳过", name, "数据集，错误:", e$message))
        empty_plot  # 失败时返回空白图，保持顺序不变
      })
    })
    
    # 定义布局
    design <- c(
      area(1, 1), area(1, 2), area(1, 3), area(1, 4), area(1, 5),  # 第一行：图1-5
      area(2, 1), area(2, 2),area(2, 3, 2, 4), area(2, 5)  # 第二行：图6-8
    )
    
    # 组合所有的图
    combined_plot <- wrap_plots(plots) + plot_layout(design = design)
    
    # 添加标题
    combined_plot <- combined_plot + 
      plot_annotation(
        title = paste("Expression of", Gene, "in MCD mouse model"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20))
        )
      )
    
    return(combined_plot)
  }
  
  
  analyze_Mouse_CDAHFD_singlegene <- function(Gene) {
    
    # 数据集列表
    datasets <- list(
      
      GSE263975 = list(data = dds_GSE263975, colData = colData_GSE263975, func = analyze_gene_expression_GSE263975),
      GSE239861 = list(data = dds_GSE239861, colData = colData_GSE239861, func = analyze_gene_expression_GSE239861),
      GSE233051 = list(data = dds_GSE233051, colData = colData_GSE233051, func = analyze_gene_expression_GSE233051),
      GSE197322 = list(data = dds_GSE197322, colData = colData_GSE197322, func = analyze_gene_expression_GSE197322),
      GSE154892 = list(data = Affymetrix_GSE154892, colData = colData_GSE154892, func = analyze_gene_expression_GSE154892),
      GSE120977 = list(data = dds_GSE120977, colData = colData_GSE120977, func = analyze_gene_expression_GSE120977),
      GSE35961  = list(data = Affymetrix_GSE35961, colData = colData_GSE35961, func = analyze_gene_expression_GSE35961),
      GSE174543  = list(data = dds_GSE174543, colData = colData_GSE174543, func = analyze_gene_expression_GSE174543),
      
      
      GSE207856 = list(data = dds_GSE207856, colData = colData_GSE207856, func = analyze_gene_expression_GSE207856),  
      GSE200409 = list(data = Agilent_GSE200409, colData = colData_GSE200409, func = analyze_gene_expression_GSE200409),
      GSE269493 = list(data = dds_GSE269493, colData = colData_GSE269493, func = analyze_gene_expression_GSE269493)
      
      
    )
    
    # 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + labs(x = NULL)
    }
    # 定义一个占位图（空白图）
    empty_plot <- ggplot() + 
      theme_void() + 
      ggtitle("Data Unavailable")
    
    # 处理所有数据集，捕获可能的错误
    plots <- lapply(names(datasets), function(name) {
      dataset <- datasets[[name]]
      tryCatch({
        plot <- dataset$func(Gene, dataset$data, dataset$colData)
        modify_plot(plot)  # 统一格式化
      }, error = function(e) {
        message(paste("跳过", name, "数据集，错误:", e$message))
        empty_plot  # 失败时返回空白图，保持顺序不变
      })
    })
    
    # 定义布局
    design <- c(
      area(1, 1), area(1, 2), area(1, 3), area(1, 4), area(1, 5), area(1, 6), area(1, 7),area(1, 8),# 第一行：图1-5
      area(2, 1,2, 2), area(2, 3, 2, 4), area(2, 5, 2, 7)  # 第二行：图6-8
    )
    
    # 组合所有的图
    combined_plot <- wrap_plots(plots) + plot_layout(design = design)
    
    # 添加标题
    combined_plot <- combined_plot + 
      plot_annotation(
        title = paste("Expression of", Gene, "in CDAHFD mouse model"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20))
        )
      )
    
    return(combined_plot)
    
  }
  
  
  analyze_Mouse_HFD_singlegene <- function(Gene) {
    
    # 数据集列表
    datasets <- list(
      GSE189066 = list(data = dds_GSE189066, colData = colData_GSE189066, func = analyze_gene_expression_GSE189066),
      GSE205846 = list(data = dds_GSE205846, colData = colData_GSE205846, func = analyze_gene_expression_GSE205846),
      GSE218025 = list(data = dds_GSE218025, colData = colData_GSE218025, func = analyze_gene_expression_GSE218025),
      GSE288077 = list(data = dds_GSE288077_Mouse, colData = colData_GSE288077_Mouse, func = analyze_gene_expression_GSE288077_Mouse),
      GSE186116 = list(data = dds_GSE186116, colData = colData_GSE186116, func = analyze_gene_expression_GSE186116),
      GSE94790 = list(data = Agilent_GSE94790, colData = colData_GSE94790, func = analyze_gene_expression_GSE94790),
      GSE165855 = list(data = dds_GSE165855, colData = colData_GSE165855, func = analyze_gene_expression_GSE165855),
      GSE126790 = list(data = dds_GSE126790, colData = colData_GSE126790, func = analyze_gene_expression_GSE126790),
      GSE145665 = list(data = dds_GSE145665, colData = colData_GSE145665, func = analyze_gene_expression_GSE145665),
      GSE129306 = list(data = dds_GSE129306, colData = colData_GSE129306, func = analyze_gene_expression_GSE129306),
      GSE242881 = list(data = dds_GSE242881, colData = colData_GSE242881, func = analyze_gene_expression_GSE242881),
      GSE226496 = list(data = dds_GSE226496, colData = colData_GSE226496, func = analyze_gene_expression_GSE226496),
      GSE195798 = list(data = dds_GSE195798, colData = colData_GSE195798, func = analyze_gene_expression_GSE195798),
      GSE179394 = list(data = dds_GSE179394, colData = colData_GSE179394, func = analyze_gene_expression_GSE179394),
      GSE24031 = list(data = Affymetrix_GSE24031, colData = colData_GSE24031, func = analyze_gene_expression_GSE24031),
      GSE59042 = list(data = Agilent_GSE59042, colData = colData_GSE59042, func = analyze_gene_expression_GSE59042),
      GSE43106 = list(data = Affymetrix_GSE43106, colData = colData_GSE43106, func = analyze_gene_expression_GSE43106),
      GSE109345 = list(data = dds_GSE109345, colData = colData_GSE109345, func = analyze_gene_expression_GSE109345),
      GSE201819 = list(data = dds_GSE201819, colData = colData_GSE201819, func = analyze_gene_expression_GSE201819)
    )
    
    # 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + labs(x = NULL)
    }
    
    # 定义一个占位图（空白图）
    empty_plot <- ggplot() + 
      theme_void() + 
      ggtitle("Data Unavailable")
    
    # 处理所有数据集，捕获可能的错误
    plots <- lapply(names(datasets), function(name) {
      dataset <- datasets[[name]]
      tryCatch({
        plot <- dataset$func(Gene, dataset$data, dataset$colData)
        modify_plot(plot)  # 统一格式化
      }, error = function(e) {
        message(paste("跳过", name, "数据集，错误:", e$message))
        empty_plot  # 失败时返回空白图，保持顺序不变
      })
    })
    
    # 定义布局
    design <- c(
      area(1, 1), area(1, 2), area(1, 3), area(1, 4), area(1, 5), area(1, 6), area(1, 7), area(1, 8), # 第一行：图1-5
      area(2, 1), area(2, 2), area(2, 3), area(2, 4), area(2, 5), area(2, 6), area(2, 7, 2, 8), 
      area(3, 1,3,2), area(3, 3,3,4), area(3, 5,3,6), area(3, 7,3,8) # 第二行：图6-8
    )
    
    # 组合所有的图
    combined_plot <- wrap_plots(plots) + plot_layout(design = design)
    
    # 添加标题
    combined_plot <- combined_plot + 
      plot_annotation(
        title = paste("Expression of", Gene, "in HFD mouse model"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20))
        )
      )
    
    return(combined_plot)
  }
  
  
  analyze_Mouse_HFHC_singlegene <- function(Gene) {
    
    # 数据集列表
    datasets <- list(
      GSE53381 = list(data = Agilent_GSE53381, colData = colData_GSE53381, func = analyze_gene_expression_GSE53381),
      GSE51432 = list(data = Affymetrix_GSE51432, colData = colData_GSE51432, func = analyze_gene_expression_GSE51432),
      GSE229188 = list(data = dds_GSE229188, colData = colData_GSE229188, func = analyze_gene_expression_GSE229188)
    )
    
    # 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + labs(x = NULL)
    }
    
    # 定义一个占位图（空白图）
    empty_plot <- ggplot() + 
      theme_void() + 
      ggtitle("Data Unavailable")
    
    # 处理所有数据集，捕获可能的错误
    plots <- lapply(names(datasets), function(name) {
      dataset <- datasets[[name]]
      tryCatch({
        plot <- dataset$func(Gene, dataset$data, dataset$colData)
        modify_plot(plot)  # 统一格式化
      }, error = function(e) {
        message(paste("跳过", name, "数据集，错误:", e$message))
        empty_plot  # 失败时返回空白图，保持顺序不变
      })
    })
    
    # 定义布局
    design <- c(
      area(1, 1), area(1, 2), area(1, 3)
      # area(2, 1), area(2, 2),area(2, 3, 2, 4), area(2, 5)  # 第二行：图6-8
    )
    
    # 组合所有的图
    combined_plot <- wrap_plots(plots) + plot_layout(design = design)
    
    # 添加标题
    combined_plot <- combined_plot + 
      plot_annotation(
        title = paste("Expression of", Gene, "in HFHC mouse model"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20))
        )
      )
    
    return(combined_plot)
  }
  
  
  analyze_Mouse_HFHS_singlegene <- function(Gene) {
    
    # 数据集列表
    datasets <- list(
      GSE154426 = list(data = dds_GSE154426, colData = colData_GSE154426, func = analyze_gene_expression_GSE154426),
      GSE67678 = list(data = Illumina_GSE67678, colData = colData_GSE67678, func = analyze_gene_expression_GSE67678),
      GSE94593 = list(data = dds_GSE94593, colData = colData_GSE94593, func = analyze_gene_expression_GSE94593),
      GSE253217 = list(data = dds_GSE253217, colData = colData_GSE253217, func = analyze_gene_expression_GSE253217),
      GSE215225 = list(data = dds_GSE215225, colData = colData_GSE215225, func = analyze_gene_expression_GSE215225),
      GSE273290 = list(data = dds_GSE273290, colData = colData_GSE273290, func = analyze_gene_expression_GSE273290),
      GSE162211 = list(data = dds_GSE162211, colData = colData_GSE162211, func = analyze_gene_expression_GSE162211),
      GSE230745 = list(data = dds_GSE230745, colData = colData_GSE230745, func = analyze_gene_expression_GSE230745),
      GSE110404 = list(data = dds_GSE110404, colData = colData_GSE110404, func = analyze_gene_expression_GSE110404),
      GSE197884 = list(data = dds_GSE197884, colData = colData_GSE197884, func = analyze_gene_expression_GSE197884)
    )
    
    # 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + labs(x = NULL)
    }
    
    # 定义一个占位图（空白图）
    empty_plot <- ggplot() + 
      theme_void() + 
      ggtitle("Data Unavailable")
    
    # 处理所有数据集，捕获可能的错误
    plots <- lapply(names(datasets), function(name) {
      dataset <- datasets[[name]]
      tryCatch({
        plot <- dataset$func(Gene, dataset$data, dataset$colData)
        modify_plot(plot)  # 统一格式化
      }, error = function(e) {
        message(paste("跳过", name, "数据集，错误:", e$message))
        empty_plot  # 失败时返回空白图，保持顺序不变
      })
    })
    
    # 定义布局
    design <- c(
      area(1, 1), area(1, 2), area(1, 3), area(1, 4), area(1, 5),  # 第一行：图1-5
      area(2, 1), area(2, 2),area(2, 3), area(2, 4),area(2, 5)  # 第二行：图6-8
    )
    
    # 组合所有的图
    combined_plot <- wrap_plots(plots) + plot_layout(design = design)
    
    # 添加标题
    combined_plot <- combined_plot + 
      plot_annotation(
        title = paste("Expression of", Gene, "in HFHS mouse model"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20))
        )
      )
    
    return(combined_plot)
  }
  
  
  analyze_Mouse_HFHSHC_singlegene <- function(Gene) {
    
    # 数据集列表
    datasets <- list(
      GSE196908 = list(data = dds_GSE196908, colData = colData_GSE196908, func = analyze_gene_expression_GSE196908),
      GSE256063 = list(data = dds_GSE256063, colData = colData_GSE256063, func = analyze_gene_expression_GSE256063),
      GSE52748 = list(data = Affymetrix_GSE52748, colData = colData_GSE52748, func = analyze_gene_expression_GSE52748),
      GSE133566 = list(data = dds_GSE133566, colData = colData_GSE133566, func = analyze_gene_expression_GSE133566),
      GSE256501 = list(data = dds_GSE256501, colData = colData_GSE256501, func = analyze_gene_expression_GSE256501),
      GSE200482 = list(data = dds_GSE200482, colData = colData_GSE200482, func = analyze_gene_expression_GSE200482),
      GSE160292 = list(data = dds_GSE160292, colData = colData_GSE160292, func = analyze_gene_expression_GSE160292),
      GSE220519 = list(data = dds_GSE220519, colData = colData_GSE220519, func = analyze_gene_expression_GSE220519),
      GSE190140 = list(data = dds_GSE190140, colData = colData_GSE190140, func = analyze_gene_expression_GSE190140),
      GSE162869 = list(data = dds_GSE162869, colData = colData_GSE162869, func = analyze_gene_expression_GSE162869),
      GSE190982 = list(data = dds_GSE190982, colData = colData_GSE190982, func = analyze_gene_expression_GSE190982),
      GSE126204 = list(data = dds_GSE126204, colData = colData_GSE126204, func = analyze_gene_expression_GSE126204),
      GSE164084 = list(data = dds_GSE164084, colData = colData_GSE164084, func = analyze_gene_expression_GSE164084),
      GSE233767 = list(data = dds_GSE233767, colData = colData_GSE233767, func = analyze_gene_expression_GSE233767),
      GSE168069 = list(data = dds_GSE168069, colData = colData_GSE168069, func = analyze_gene_expression_GSE168069),
      GSE195483 = list(data = dds_GSE195483, colData = colData_GSE195483, func = analyze_gene_expression_GSE195483),
      GSE243976 = list(data = dds_GSE243976, colData = colData_GSE243976, func = analyze_gene_expression_GSE243976)
    )
    
    # 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + labs(x = NULL)
    }
    
    # 定义一个占位图（空白图）
    empty_plot <- ggplot() + 
      theme_void() + 
      ggtitle("Data Unavailable")
    
    # 处理所有数据集，捕获可能的错误
    plots <- lapply(names(datasets), function(name) {
      dataset <- datasets[[name]]
      tryCatch({
        plot <- dataset$func(Gene, dataset$data, dataset$colData)
        modify_plot(plot)  # 统一格式化
      }, error = function(e) {
        message(paste("跳过", name, "数据集，错误:", e$message))
        empty_plot  # 失败时返回空白图，保持顺序不变
      })
    })
    
    # 定义布局
    design <- c(
      area(1, 1), area(1, 2), area(1, 3), area(1, 4), area(1, 5),area(1, 6),  # 第一行：图1-5
      area(2, 1), area(2, 2),area(2, 3),area(2, 4), area(2, 5, 2, 6),  # 第二行：图6-8
      area(3, 1), area(3, 2),area(3, 3),area(3, 4), area(3, 5), area(3, 6)
    )
    
    # 组合所有的图
    combined_plot <- wrap_plots(plots) + plot_layout(design = design)
    
    # 添加标题
    combined_plot <- combined_plot + 
      plot_annotation(
        title = paste("Expression of", Gene, "in HFHSHC mouse model"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20))
        )
      )
    
    return(combined_plot)
  }
  
  
  analyze_Mouse_WDCCL4_singlegene <- function(Gene) {
    
    # 数据集列表
    datasets <- list(
      GSE184256 = list(data = dds_GSE184256, colData = colData_GSE184256, func = analyze_gene_expression_GSE184256),
      GSE287943 = list(data = dds_GSE287943, colData = colData_GSE287943, func = analyze_gene_expression_GSE287943)
    )
    
    # 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + labs(x = NULL)
    }
    
    # 定义一个占位图（空白图）
    empty_plot <- ggplot() + 
      theme_void() + 
      ggtitle("Data Unavailable")
    
    # 处理所有数据集，捕获可能的错误
    plots <- lapply(names(datasets), function(name) {
      dataset <- datasets[[name]]
      tryCatch({
        plot <- dataset$func(Gene, dataset$data, dataset$colData)
        modify_plot(plot)  # 统一格式化
      }, error = function(e) {
        message(paste("跳过", name, "数据集，错误:", e$message))
        empty_plot  # 失败时返回空白图，保持顺序不变
      })
    })
    
    # 定义布局
    design <- c(
      area(1, 1), area(1, 2)
    )
    
    # 组合所有的图
    combined_plot <- wrap_plots(plots) + plot_layout(design = design)
    
    # 添加标题
    combined_plot <- combined_plot + 
      plot_annotation(
        title = paste("Expression of", Gene, "in WDCCL4 mouse model"),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20))
        )
      )
    
    return(combined_plot)
  }
  
  
  get_orthologs <- function(Gene) {
    tryCatch({
      # 检查基因名是否包含小写字母（如果包含，说明是小鼠基因）
      if (grepl("[a-z]", Gene)) {
        # 小鼠基因 -> 查询对应的人类基因
        result <- orthologs(genes = Gene, species = "mouse", human = FALSE)
        Gene_human <- result$human_symbol
        Gene_mouse <- Gene
      } else {
        # 人类基因 -> 查询对应的小鼠基因
        result <- orthologs(genes = Gene, species = "mouse")
        Gene_human <- Gene
        Gene_mouse <- result$symbol
      }
      
      # 如果查询结果为空，则返回输入基因作为默认值
      if (is.null(Gene_human) || is.null(Gene_mouse) ||
          length(Gene_human) == 0 || length(Gene_mouse) == 0) {
        Gene_human <- Gene
        Gene_mouse <- Gene
      }
      
      # 返回一个包含人类和小鼠基因的列表
      return(list(human = Gene_human, mouse = Gene_mouse))
      
    }, error = function(e) {
      message("Error: Unable to retrieve orthologs for gene: ", Gene, "\n", e$message)
      # 查询出错时返回输入基因名称
      return(list(human = Gene, mouse = Gene))
    })
  }
  
  # 监听基因输入并更新图表
  observeEvent(input$update, {
    req(user_authenticated())  # 确保登录后才能运行
    
    shinyjs::show("spinner_area")  # ✅ 开始时显示

    shinyjs::runjs("$('#plot_mask').show();")           # 显示遮罩
    
    # 使用 isolate 确保不会触发其他依赖
    isolate({
      
      species <- input$species  # 获取物种（Human 或 Mouse）
      analysis_type <- input$analysis_type  # 获取分析类型
      feed_type <- input$feed_type  # 获取 Mouse 造模方式
      multi_gene_analysis <- input$multi_gene_analysis  # 获取多基因分析子选项
      
      # 判断是 "多基因分析" 还是普通的单基因分析
      if (analysis_type == "多基因分析" && multi_gene_analysis == "双基因相关性") {
        # 获取两个基因输入
        GeneX <- get_orthologs(input$gene1)
        GeneY <- get_orthologs(input$gene2)
        
        # 选择相关性分析
        combined_plot <- analyze_correlation_doublegene(GeneX$human, GeneY$human)
        plot_width <- 800
        plot_height <- 600
        
      } else {
        # 处理单基因输入
        Gene <- get_orthologs(input$gene)
        
        # 根据分析类型选择不同的分析函数
        if (species == "Mouse") {
          if (feed_type == "STAM") {
            combined_plot <- analyze_Mouse_STAM_singlegene(Gene$mouse)
            plot_width <- 750
            plot_height <- 733
          } else if (feed_type == "MCD") {
            combined_plot <- analyze_Mouse_MCD_singlegene(Gene$mouse)
            plot_width <- 750
            plot_height <- 733
          } else if (feed_type == "CDAHFD") {
            combined_plot <- analyze_Mouse_CDAHFD_singlegene(Gene$mouse)
            plot_width <- 1200
            plot_height <- 733
          } else if (feed_type == "HFD") {
            combined_plot <- analyze_Mouse_HFD_singlegene(Gene$mouse)
            plot_width <- 1200
            plot_height <- 1100
          } else if (feed_type == "HFHC") {
            combined_plot <- analyze_Mouse_HFHC_singlegene(Gene$mouse)
            plot_width <- 450
            plot_height <- 366
          } else if (feed_type == "HFHS") {
            combined_plot <- analyze_Mouse_HFHS_singlegene(Gene$mouse)
            plot_width <- 750
            plot_height <- 733
          } else if (feed_type == "HFHSHC") {
            combined_plot <- analyze_Mouse_HFHSHC_singlegene(Gene$mouse)
            plot_width <- 900
            plot_height <- 1100
          } else if (feed_type == "WDCCL4") {
            combined_plot <- analyze_Mouse_WDCCL4_singlegene(Gene$mouse)
            plot_width <- 300
            plot_height <- 366
          } else {
            combined_plot <- analyze_Mouse_STAM_singlegene(Gene$mouse)
            plot_width <- 750
            plot_height <- 733
          }
        } else {
          if (analysis_type == "NAFLD Score") {
            combined_plot <- analyze_nafld_singlegene(Gene$human)
            plot_width <- 1200
            plot_height <- 1100
          } else if (analysis_type == "Classification") {
            combined_plot <- analyze_classification_singlegene(Gene$human)
            plot_width <- 1400
            plot_height <- 1466
          } else if (analysis_type == "Summary") {
            combined_plot <- analyze_summary_singlegene(Gene$human)
            plot_width <- 700
            plot_height <- 530
          } else if (analysis_type == "多基因分析" && multi_gene_analysis == "寻找相关基因") {
            combined_plot <- analyze_find_correlationGene(Gene$human)
            plot_width <- 700
            plot_height <- 530
          } else {
            combined_plot <- analyze_fibrosis_singlegene(Gene$human)
            plot_width <- 1200
            plot_height <- 1100
          }
        }
      }
    })
    
    # 动态渲染 `plotOutput()`
    output$dynamic_plot <- renderUI({
      plotOutput("combined_plot", width = paste0(plot_width, "px"), height = paste0(plot_height, "px"))
    })
    
    # 生成图像
    output$combined_plot <- renderPlot({
      combined_plot
    })
    
    shinyjs::delay(2000, shinyjs::runjs("$('#plot_mask').fadeOut();"))  # 延迟后隐藏
    
    shinyjs::delay(2000, shinyjs::hide("spinner_area"))
  })
}

# 运行 Shiny 应用
shinyApp(ui = ui, server = server)
