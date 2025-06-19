
# 网页端设置 ----------

# 服务器端逻辑
server <- function(input, output, session) {
  #初始化设置
  options(future.globals.maxSize = 3 * 1024^3) # 3GB
  
  #日志说明
  write_log <- function(user, action) {
    log_file <- "~/R/RNA/data/NASH/网页客户端/user_log.csv"
    time <- Sys.time()
    log_entry <- data.frame(Time = time, User = user, Action = action)
    if (!file.exists(log_file)) {
      write.csv(log_entry, log_file, row.names = FALSE)
    } else {
      write.table(log_entry, log_file, row.names = FALSE, col.names = FALSE,
                  sep = ",", append = TRUE)
    }
  }
  
  
  # ---- 输入合法性检查 -------------------------------------------------
  is_valid_gene <- function(symbol) {
    # 基本正则：字母/数字/下划线，不能全为空
    nzchar(symbol) && grepl("^[A-Za-z0-9._-]+$", symbol)
  }
  
  check_gene_input <- function(raw_symbol, reference_mat) {
    # raw_symbol  : 用户输入的字符
    # reference_mat: 你想用来判断“是否存在” 的表达矩阵
    if (!is_valid_gene(raw_symbol)) {
      # showNotification("⚠️ 请输入合法的基因符号！",
      #                  type = "error", duration = 5)
      return(FALSE)
    }
    if (!(raw_symbol %in% rownames(reference_mat))) {
      # showNotification(paste0("⚠️ 数据库中找不到基因：", raw_symbol),
      #                  type = "error", duration = 5)
      return(FALSE)
    }
    TRUE
  }
  
  
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
      write_log(input$user, "登录成功")
      output$login_status <- renderText("")  # 清空错误信息
    } else {
      user_authenticated(FALSE)
      output$login_status <- renderText("❌ 用户名或密码错误")
      # 控制台打印登录成功 + 时间
      login_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      write_log(input$user, "登录失败")
    }
    
  })
  

  # 读取使用说明文本
  usage_html <- paste(readLines("~/R/RNA/data/NASH/网页客户端/www/MGA_User_Guide.html", encoding = "UTF-8"), collapse = "\n")
  
  # Help 页中展示 HTML 格式内容
  output$help_text <- renderUI({
    HTML(usage_html)
  })
  
  # 首页显示静态 HTML 文件
  output$home_html <- renderUI({
    includeHTML("~/R/RNA/data/NASH/网页客户端/www/mash_intro.html")
  })
  
  
  # 在用户点击 Analyze 后淡出说明
  observeEvent(input$update, {
    shinyjs::hide("usage_section", anim = TRUE, animType = "fade", time = 1)
  })
  
  
  output$user_authenticated <- reactive({
    user_authenticated()
  })
  outputOptions(output, "user_authenticated", suspendWhenHidden = FALSE)
  
  # ✅ 物种切换时重置相关选项（这里放非常合适）
  observeEvent(input$species, {
    if (input$species == "Mouse") {
      updateRadioButtons(session, "analysis_type", selected = character(0))
      updateRadioButtons(session, "multi_gene_analysis", selected = "Find Correlated Genes")
    } else {
      updateRadioButtons(session, "analysis_type", selected = "Summary")
    }
  })
  
  # 定义图像格式
  PLOT_MARGIN_PT <- 20
  
  
  # 定义分析函数
  analyze_nafld_singlegene <- function(Gene) {
    
    # 统一格式化 ggplot 主题
    modify_plot <- function(plot) {
      plot +
        theme(
          plot.title  = element_text(hjust = 0.5, size = 14, face = "bold"),
          # 边距统一用全局常量
          plot.margin = margin(
            PLOT_MARGIN_PT, PLOT_MARGIN_PT,
            PLOT_MARGIN_PT, PLOT_MARGIN_PT, unit = "pt"
          )
        ) +
        labs(x = NULL)
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
      plot +
        theme(
          plot.title  = element_text(hjust = 0.5, size = 14, face = "bold"),
          # 边距统一用全局常量
          plot.margin = margin(
            PLOT_MARGIN_PT, PLOT_MARGIN_PT,
            PLOT_MARGIN_PT, PLOT_MARGIN_PT, unit = "pt"
          )
        ) +
        labs(x = NULL)
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
      plot +
        theme(
          plot.title  = element_text(hjust = 0.5, size = 14, face = "bold"),
          # 边距统一用全局常量
          plot.margin = margin(
            PLOT_MARGIN_PT, PLOT_MARGIN_PT,
            PLOT_MARGIN_PT, PLOT_MARGIN_PT, unit = "pt"
          )
        ) +
        labs(x = NULL)
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
    
    
    # ------ 修补行 ------
    analysis_type       <- input$analysis_type        # ← 新增
    multi_gene_analysis <- input$multi_gene_analysis  # ← 新增
    
    # ---------- ① 先做输入合法性检查 ----------
    if (analysis_type == "Multi-Gene Analysis" &&
        multi_gene_analysis == "Double Gene Correlation") {
      
      if (!check_gene_input(input$gene1, combined_data) ||
          !check_gene_input(input$gene2, combined_data)) {
        
        shinyWidgets::sendSweetAlert(
          session,
          title = "❌ Invalid input",
          text  = "Please enter two valid gene symbols that exist in the database!",
          type  = "error"
        )
        
        shinyjs::hide("spinner_area")
        shinyjs::runjs("$('#plot_mask').hide()")
        return()
      }
      
    } else {
      
      if (!check_gene_input(input$gene, combined_data)) {
        
        # --- 单基因检查不通过 -------------------------
        shinyWidgets::sendSweetAlert(
          session,
          title = "❌ Invalid input",
          text  = "Please enter a valid gene symbol that exists in the database!",
          type  = "error"
        )
        
        shinyjs::hide("spinner_area")
        shinyjs::runjs("$('#plot_mask').hide()")
        return()
      }
    }
    
    
    
    
    # 使用 isolate 确保不会触发其他依赖
    isolate({
      
      species <- input$species  # 获取物种（Human 或 Mouse）
      analysis_type <- input$analysis_type  # 获取分析类型
      feed_type <- input$feed_type  # 获取 Mouse 造模方式
      multi_gene_analysis <- input$multi_gene_analysis  # 获取Multi-Gene Analysis子选项
      
      # 判断是 "Multi-Gene Analysis" 还是普通的单基因分析
      if (analysis_type == "Multi-Gene Analysis" && multi_gene_analysis == "Double Gene Correlation") {
        # 获取两个基因输入
        GeneX <- get_orthologs(input$gene1)
        GeneY <- get_orthologs(input$gene2)
        
        # 选择相关性分析
        combined_plot <- analyze_correlation_doublegene(GeneX$human, GeneY$human)
        plot_width <- 800
        plot_height <- 600
        
        write_log(input$user, paste("运行Double Gene Correlation分析：", GeneX, " vs ", GeneY))
        
      } else {
        # 处理单基因输入
        Gene <- get_orthologs(input$gene)
        
        # 根据分析类型选择不同的分析函数
        if (species == "Mouse") {
          if (feed_type == "STAM") {
            combined_plot <- analyze_Mouse_STAM_singlegene(Gene$mouse)
            write_log(input$user, paste(Gene, "的STAM模型分析"))
            plot_width <- 750
            plot_height <- 733
          } else if (feed_type == "MCD") {
            combined_plot <- analyze_Mouse_MCD_singlegene(Gene$mouse)
            write_log(input$user, paste(Gene, "的MCD模型分析"))
            plot_width <- 750
            plot_height <- 733
          } else if (feed_type == "CDAHFD") {
            combined_plot <- analyze_Mouse_CDAHFD_singlegene(Gene$mouse)
            write_log(input$user, paste(Gene, "的CDAHFD模型分析"))
            plot_width <- 1200
            plot_height <- 733
          } else if (feed_type == "HFD") {
            combined_plot <- analyze_Mouse_HFD_singlegene(Gene$mouse)
            write_log(input$user, paste(Gene, "的HFD模型分析"))
            plot_width <- 1200
            plot_height <- 1100
          } else if (feed_type == "HFHC") {
            combined_plot <- analyze_Mouse_HFHC_singlegene(Gene$mouse)
            write_log(input$user, paste(Gene, "的HFHC模型分析"))
            plot_width <- 450
            plot_height <- 366
          } else if (feed_type == "HFHS") {
            combined_plot <- analyze_Mouse_HFHS_singlegene(Gene$mouse)
            write_log(input$user, paste(Gene, "的HFHS模型分析"))
            plot_width <- 750
            plot_height <- 733
          } else if (feed_type == "HFHSHC") {
            combined_plot <- analyze_Mouse_HFHSHC_singlegene(Gene$mouse)
            write_log(input$user, paste(Gene, "的HFHSHC模型分析"))
            plot_width <- 900
            plot_height <- 1100
          } else if (feed_type == "WDCCL4") {
            combined_plot <- analyze_Mouse_WDCCL4_singlegene(Gene$mouse)
            write_log(input$user, paste(Gene, "的WDCCL4模型分析"))
            plot_width <- 300
            plot_height <- 366
          } else {
            combined_plot <- analyze_Mouse_STAM_singlegene(Gene$mouse)
            write_log(input$user, paste(Gene, "的STAM模型分析"))
            plot_width <- 750
            plot_height <- 733
          }
        } else {
          if (analysis_type == "NAFLD Score") {
            combined_plot <- analyze_nafld_singlegene(Gene$human)
            write_log(input$user, paste(Gene, "的NAFLD Score分析"))
            plot_width <- 1200
            plot_height <- 1100
          } else if (analysis_type == "Classification") {
            combined_plot <- analyze_classification_singlegene(Gene$human)
            write_log(input$user, paste(Gene, "的Classification分析"))
            plot_width <- 1400
            plot_height <- 1466
          } else if (analysis_type == "Summary") {
            combined_plot <- analyze_summary_singlegene(Gene$human)
            write_log(input$user, paste(Gene, "的Summary分析"))
            plot_width <- 700
            plot_height <- 530
          } else if (analysis_type == "Multi-Gene Analysis" && multi_gene_analysis == "Find Correlated Genes") {
            combined_plot <- analyze_find_correlationGene(Gene$human)
            write_log(input$user, paste(Gene, "的Multi-Gene Analysis"))
            plot_width <- 700
            plot_height <- 530
          } else {
            combined_plot <- analyze_fibrosis_singlegene(Gene$human)
            write_log(input$user, paste(Gene, "的Fibrosis score分析"))
            plot_width <- 1200
            plot_height <- 1100
          }
        }
      }
    })
    
    # 动态渲染 `plotOutput()`
    output$dynamic_plot <- renderUI({
      div(
        style = "text-align: left;",
        div(
          style = "margin-left: 20px;",
          plotOutput("combined_plot", width = paste0(plot_width, "px"), height = paste0(plot_height, "px"))
        )
      )
    })
    
    # 生成图像
    output$combined_plot <- renderPlot({
      combined_plot
    })
    
    shinyjs::delay(2000, shinyjs::runjs("$('#plot_mask').fadeOut();"))  # 延迟后隐藏
    
    shinyjs::delay(2000, shinyjs::hide("spinner_area"))
  })
}