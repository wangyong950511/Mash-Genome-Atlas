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
library(bslib)
library(fontawesome)
library(shinyWidgets)


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




