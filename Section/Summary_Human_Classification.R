
# 库加载 ------------------------------------
library(DESeq2)
library(sva)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(dplyr)
library(data.table)
library(AnnotationDbi)
library(babelgene)
library(patchwork)



# 函数集合 ------------------------------------
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
      text = element_text(size = FIXED_TEXT_SIZE)
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
      text = element_text(size = FIXED_TEXT_SIZE)
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
      text = element_text(size = FIXED_TEXT_SIZE)
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
      text = element_text(size = FIXED_TEXT_SIZE)
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
      text = element_text(size = FIXED_TEXT_SIZE)
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
      text = element_text(size = FIXED_TEXT_SIZE)
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
      text = element_text(size = FIXED_TEXT_SIZE)
    )
  
  return(plot)
}



# 分析 ------------------------------------
# 固定文字大小常量
FIXED_TEXT_SIZE <- 10
# 单基因分析
Gene <- "LONP1"
# 检测基因是否包含小写字母（判定为小鼠基因）
if (grepl("[a-z]", Gene)) {
  # 使用 babelgene::orthologs() 查询同源人类基因
  result <- orthologs(genes = Gene, species = "mouse", human = FALSE)
  # 提取基因名
  Gene_hum <- result$human_symbol
  Gene_mum <- Gene
} else{
  # 使用 babelgene::orthologs() 查询同源人类基因
  result <- orthologs(genes = Gene, species = "mouse")
  Gene_hum <- Gene
  Gene_mum <- result$symbol
}


plot_GSE24807 <- analyze_gene_expression_GSE24807(Gene, codelink_GSE24807, colData_GSE24807)+ theme(legend.position = "none")
plot_GSE37031 <- analyze_gene_expression_GSE37031(Gene, Affymetrix_GSE37031, colData_GSE37031)+ theme(legend.position = "none")
plot_GSE106737 <- analyze_gene_expression_GSE106737(Gene,Affymetrix_GSE106737,colData_GSE106737)+ theme(legend.position = "none")
plot_GSE134146 <- analyze_gene_expression_GSE134146(Gene, Arraystar_GSE134146, colData_GSE134146)+ theme(legend.position = "none")
plot_GSE159676 <- analyze_gene_expression_GSE159676(Gene,Affymetrix__GSE159676,colData_GSE159676)+ theme(legend.position = "none")
plot_GSE147304 <- analyze_gene_expression_GSE147304(Gene, vsd_GSE147304, colData_GSE147304)+ theme(legend.position = "none")
plot_GSE173735 <- analyze_gene_expression_GSE173735(Gene, vsd_GSE173735, colData_GSE173735)+ theme(legend.position = "none")



combined_plot <- plot_GSE24807 +
  plot_GSE37031 +
  plot_GSE106737 +
  plot_GSE134146 +
  plot_GSE159676 +
  plot_GSE147304 +
  plot_GSE173735 +
  plot_layout(ncol = 4)  # 设置每行两张图

# 显示组合后的图形
print(combined_plot)

ggsave("combined_plot.eps", combined_plot, width = 12, height = 18)

