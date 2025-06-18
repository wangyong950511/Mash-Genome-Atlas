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
        h3("User Login", style = "text-align: center;"),
        textInput("user", "User", value = "test", width = "100%"),
        passwordInput("password", "Password", value = "", width = "100%"),
        div(actionButton("login", "Login", class = "btn-primary"), 
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
        radioButtons("species", "Species", 
                     choices = c("Human", "Mouse")),  
        
        # 仅当选择 Human 时，显示分析类型
        conditionalPanel(
          condition = "input.species == 'Human'",
          radioButtons("analysis_type", "Select Analysis Type", 
                       choices = c("Summary", "Classification", "NAFLD Score", "Fibrosis Stages",  "Multi-Gene Analysis"))
        ),
        
        # 仅当选择 Mouse 时，显示造模方式
        conditionalPanel(
          condition = "input.species == 'Mouse'",
          radioButtons("feed_type", "Select Feed Type", 
                       choices = c("HFD", "HFHS", "HFHC", "HFHSHC", "WDCCL4", "STAM", "MCD", "CDAHFD"))
        ),
        
        # "Multi-Gene Analysis" 需要额外的子选项
        conditionalPanel(
          condition = "input.species == 'Human' && input.analysis_type == 'Multi-Gene Analysis'",
          radioButtons("multi_gene_analysis", "Select Type", 
                       choices = c("Find Correlated Genes", "Double Gene Correlation"),
                       selected = "Find Correlated Genes")  # 默认选择 "Find Correlated Genes"
        ),
        
        # 基因输入框 (适用于所有分析类型，但“Multi-Gene Analysis”有特殊情况)
        conditionalPanel(
          # 让 Mouse 物种也显示基因输入框
          condition = "(input.species == 'Human' && (input.analysis_type != 'Multi-Gene Analysis' || input.multi_gene_analysis == 'Find Correlated Genes')) || input.species == 'Mouse'",
          textInput("gene", 
                    label = HTML("Enter Gene Name<br/><span style='font-size:smaller'>(Standard Gene Symbol, e.g., COL1A1, Col1a1)</span>"),
                    value = "")
        ),
        
        # "Double Gene Correlation" -> 显示两个基因输入框
        conditionalPanel(
          condition = "input.species == 'Human' && input.analysis_type == 'Multi-Gene Analysis' && input.multi_gene_analysis == 'Double Gene Correlation'",
          textInput("gene1", "Enter Gene 1", value = ""),
          textInput("gene2", "Enter Gene 2", value = "")
        ),
        
        fluidRow(
          column(6, actionButton("update", "Analyze", class = "btn-primary")),
          column(6, div(id = "spinner_area", img(src = "https://i.gifer.com/ZKZg.gif", height = "25px"), style = "display: none;"))
        )
    ),
    
    div(class = "main-content",
        uiOutput("usage_ui"),    # 使用说明，生成图像后可在 server 中隐藏
        div(id = "plot_mask"),   # 遮罩层
        uiOutput("dynamic_plot") # 图像输出
    )
  )
)
