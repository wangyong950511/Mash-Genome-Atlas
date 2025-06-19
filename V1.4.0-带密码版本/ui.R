

ui <- fluidPage(
  shinyjs::useShinyjs(),
  
  # ---------- global <head> ----------
  tags$head(
    tags$style(HTML("
      body{position:relative;margin:0;height:100vh;}

      #login_panel{position:absolute;top:50%;left:50%;
        transform:translate(-50%,-50%);background:#f8f9fa;padding:30px;
        border-radius:10px;box-shadow:0 0 15px rgba(0,0,0,.1);z-index:1000;}

      #sidebar_panel{width:300px;margin-left:0;}

      #plot_mask{
        position:absolute;      /* 相对最近的 position:relative; 的父元素(main_area_db) */
        top:0; left:0;          /* 完全贴合容器四边 */
        width:100%; height:100%;
        background:white;
        z-index:999;            /* 只需盖住图即可，不必 9999，这样不会压住导航栏 */
        display:none;           /* 需要时再通过 JS/shinyjs show() */
      }

      .navbar-custom{
        background:#e6f2ff;        /* 保持原浅蓝底 */
        color:#2c3e50;
      
        display:flex;              /* ① 让标题条本身变成 flex 容器 */
        justify-content:center;    /* ② 水平居中子元素 */
        align-items:center;        /* ③ 垂直居中子元素 */
      
        height:56px;               /* ④ 与 .navbar-default 高度一致 */
        margin-bottom:-10px;       /* ⑤ 向下偏移一半高度，使文字位于两条色带正中央 */
      
        padding:0 10px;            /* 去掉上下 padding，只保留左右 */
        font-size:30px;
        font-weight:bold;
        text-align:center;         /* 冗余安全设置，可保留 */
      }

      /* ---- tab-navbar ---- */
      .navbar-default{background:#e6f2ff!important;border-color:#e6f2ff!important;}
      .navbar-default .navbar-nav>li>a{
        position:relative;top:8px;background:#e6f2ff;color:#333;
        border:1px solid #cce0ff;margin:5px;border-radius:6px;padding:8px 15px;
        font-weight:500;transition:all .3s ease;}
      .navbar-default .navbar-nav>li>a:hover{background:#d4eaff;color:#000;border-color:#99ccff;}
      .navbar-default .navbar-nav>.active>a,
      .navbar-default .navbar-nav>.active>a:focus,
      .navbar-default .navbar-nav>.active>a:hover{
        background:#3399ff;color:#fff!important;border-color:#3399ff;}
      .navbar-default .navbar-nav>li>a>.fa{margin-right:6px;}

      /* main_area 额外左内边距（只给 Database 页加 class） */
      .with-sidebar{padding-left:300px;}
      
      .main_area{
      padding-top: 12px !important;   /* 原来 30px → 12px，可再自行微调 */
      }

      .main_area h1:first-child,
      .main_area h2:first-child{
        margin-top: 0 !important;
      }

      /* 顶栏贴满两侧 */
      .navbar-default .container,
      .navbar-default .container-fluid{padding-left:0!important;padding-right:0!important;}
      
      /* ------- 始终显示垂直滚动条，消除标签抖动 ------- */
      html{ overflow-y: scroll; }
  
      /* ---------- 让 main_area 内的文字行距和段距舒展一些 ---------- */
      .main_area p,
      .main_area li{
      line-height: 1.8;          /* 行距 (默认约 1.2) → 1.8 可再酌情调 */
      }
      .main_area p{
        margin: 0 0 14px 0;        /* 段后空 14px；保持段前 0，更紧凑时可减小 */
      }
      .main_area ul{
        margin: 0 0 18px 24px;     /* 整个列表与上下文的间隔 + 左缩进 */
      }
      .main_area h2{               /* 小节标题与正文的距离 */
        margin-top: 28px;
        margin-bottom: 14px;
      }
      .main_area h1{
        margin-bottom: 22px;
      }
  
      /* 若你的 HTML 里还有 <div class=\"section\"> 可以单独加：*/
      .main_area .section{
        margin-bottom: 32px;
      }
    
    
    "))
  ),
  
  # ---------- login panel ----------
  conditionalPanel(
    condition = "!output.user_authenticated",
    div(
      id = "login_panel",
      h3("User Login", style = "text-align:center;"),
      textInput("user","User","test",width="100%"),
      passwordInput("password","Password","",width="100%"),
      div(actionButton("login","Login",class="btn-primary"),
          style="text-align:center;margin-top:15px;"),
      textOutput("login_status")
    ),
    div(id="login_spinner",
        img(src="ZKZg.gif",height="25px"),
        style="position:absolute;top:50%;left:calc(100% + 20px);
               transform:translateY(-50%);display:none;")
  ),
  
  # ---------- main authenticated area ----------
  conditionalPanel(
    condition = "output.user_authenticated",
    tagList(
      div(class="navbar-custom","Mash Genome Atlas (MGA) Database"),
      
      navbarPage(
        title = NULL,
        fluid = TRUE,
        
        # ---------- Home ----------
        tabPanel(
          title = tagList(icon("house"),"Home"),
          div(
            id    = "main_area_home",
            class = "main_area",
            style = "
              max-width:1500px;padding:50px;margin:0;        /* 去掉 auto 以取消水平居中 */
              text-align:left;                               /* 明确文本左对齐 */
            ",
            htmlOutput("home_html")
          )
        ),
        
        # ---------- Database ----------
        tabPanel(
          title = tagList(icon("database"),"Database"),
          div(style="display:flex;flex-direction:row;",
              
              # --- sidebar ---
              div(
                id="sidebar_panel",
                style="
                width:300px;padding:20px;background:#e6f2ff;
                align-self:flex-start;
              ",
                
                radioButtons("species","Species",choices=c("Human","Mouse")),
                
                conditionalPanel(
                  condition="input.species == 'Human'",
                  radioButtons("analysis_type","Select Analysis Type",
                               choices=c("Summary","Classification","NAFLD Score",
                                         "Fibrosis Stages","Multi-Gene Analysis"))
                ),
                
                conditionalPanel(
                  condition="input.species == 'Mouse'",
                  radioButtons("feed_type","Select Feed Type",
                               choices=c("HFD","HFHS","HFHC","HFHSHC","WDCCL4",
                                         "STAM","MCD","CDAHFD"))
                ),
                
                conditionalPanel(
                  condition="input.species == 'Human' && input.analysis_type == 'Multi-Gene Analysis'",
                  radioButtons("multi_gene_analysis","Select Type",
                               choices=c("Find Correlated Genes","Double Gene Correlation"),
                               selected="Find Correlated Genes")
                ),
                
                conditionalPanel(
                  condition="(input.species == 'Human' && (input.analysis_type != 'Multi-Gene Analysis' ||
                           input.multi_gene_analysis == 'Find Correlated Genes')) ||
                           input.species == 'Mouse'",
                  textInput(inputId = "gene",
                            label   = HTML(
                              "Enter Gene Name<br/>
                               <span style='font-size:smaller'>
                                 (Standard Gene Symbol, e.g., <code>TP53</code>, <code>Tp53</code>)
                               </span>"
                            ),value="",
                            placeholder="")
                ),
                
                conditionalPanel(
                  condition="input.species == 'Human' && input.analysis_type == 'Multi-Gene Analysis' &&
                           input.multi_gene_analysis == 'Double Gene Correlation'",
                  textInput("gene1","Enter Gene 1",value=""),
                  textInput("gene2","Enter Gene 2",value="")
                ),
                
                div(style="margin-top:10px;margin-bottom:30px;",
                    fluidRow(
                      column(6,actionButton("update","Analyze",class="btn-primary")),
                      column(6,
                             div(id="spinner_area",
                                 img(src="ZKZg.gif",style="width:20px;height:20px;border:none;box-shadow:none;"),
                                 style="display:none;text-align:center;")
                      )
                    )
                )
              ),
              
              # --- main area ---
              div(
                id    = "main_area_db",
                class = "main_area with-sidebar",
                style = "
                position:relative;flex-grow:1;max-width:1500px;
                padding:30px;
              ",
                ## ① <<< 追加这段覆盖样式 >>>
                tags$style("#main_area_db img{border:none!important;box-shadow:none!important;}"),
                
                uiOutput("usage_ui"),
                uiOutput("dynamic_plot"),
                
                # plot mask
                div(id="plot_mask")
              )
          )
        ),
        
        # ---------- Help ----------
        tabPanel(
          title = tagList(icon("circle-question"),"Help"),
          div(
            id    = "main_area_help",
            class = "main_area",
            style = "
              max-width:1500px;padding:50px;margin:0;        /* 去掉 auto 以取消水平居中 */
              text-align:left;                               /* 明确文本左对齐 */
            ",
            htmlOutput("help_text")
          )
        )
      )
    )
  )
)