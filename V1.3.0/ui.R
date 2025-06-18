ui <- fluidPage(
  shinyjs::useShinyjs(),

  # ---------- global <head> ----------
  tags$head(
    tags$style(HTML(
      "
      body {
        position: relative;
        margin: 0;
        height: 100vh;
      }

      #login_panel {
        position: absolute;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        background: #f8f9fa;
        padding: 30px;
        border-radius: 10px;
        box-shadow: 0 0 15px rgba(0, 0, 0, 0.1);
        z-index: 1000;
      }

      #sidebar_panel {
        width: 300px;
        margin-left: 0;
      }

      .main-content {
        margin-left: 100px;
        padding: 20px;
        height: 100vh;
        overflow-y: auto;
        background-color: #e6f2ff;
      }

      #plot_mask {
        position: fixed;
        top: 60px;
        left: 310px;
        width: calc(100% - 300px);
        height: calc(100% - 60px);
        background: rgba(255, 255, 255, 0.8);
        z-index: 9999;
        display: none;
      }

      #sidebar_panel .title {
        font-size: 18px;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
        text-align: center;
        width: 100%;
        padding: 10px 10px;
        line-height: 1.4;
        margin-bottom: 20px;
      }

      .navbar-custom {
        background-color: #e6f2ff;
        color: #2c3e50;
        padding: 10px 10px;
        font-size: 30px;
        font-weight: bold;
        text-align: center;
      }

      /* ---- tab-navbar styling ---- */
      .navbar-default {
        background-color: #e6f2ff !important;
        border-color:    #e6f2ff !important;
      }

      .navbar-default .navbar-nav > li > a {
        position: relative;   /* makes top offset effective */
        top: 8px;             /* move the whole button downward */
        background-color: #e6f2ff;
        color: #333;
        border: 1px solid #cce0ff;
        margin: 5px 5px;
        border-radius: 6px;
        padding: 8px 15px;
        font-weight: 500;
        transition: all 0.3s ease;
      }

      .navbar-default .navbar-nav > li > a:hover {
        background-color: #d4eaff;
        color: #000;
        border-color: #99ccff;
      }

      .navbar-default .navbar-nav > .active > a,
      .navbar-default .navbar-nav > .active > a:focus,
      .navbar-default .navbar-nav > .active > a:hover {
        background-color: #3399ff;
        color: white !important;
        border-color: #3399ff;
      }

      .navbar-default .navbar-nav > li > a > .fa {
        margin-right: 6px;  /* icon–text gap */
      }
      "
    ))
  ),

  # ---------- login panel ----------
  conditionalPanel(
    condition = "!output.user_authenticated",
    div(
      id = "login_panel",
      h3("User Login", style = "text-align: center;"),
      textInput("user", "User", value = "test", width = "100%"),
      passwordInput("password", "Password", value = "", width = "100%"),
      div(
        actionButton("login", "Login", class = "btn-primary"),
        style = "text-align: center; margin-top: 15px;"
      ),
      textOutput("login_status")
    ),
    div(
      id = "login_spinner",
      img(src = "ZKZg.gif", height = "25px"),
      style = "
        position: absolute;
        top: 50%;
        left: calc(100% + 20px);
        transform: translateY(-50%);
        display: none;
      "
    )
  ),

  # ---------- main authenticated area ----------
  conditionalPanel(
    condition = "output.user_authenticated",
    tagList(
      div(class = "navbar-custom", "Mash Genome Atlas (MGA) Database"),

      navbarPage(
        title = NULL,

        # ----- Home -----
        tabPanel(
          title = tagList(icon("house"), "Home"),
          div(
            style = "display: flex; flex-direction: row;",

            # hidden placeholder sidebar
            div(
              id = "sidebar_panel",
              style = "
                width: 300px;
                padding: 20px;
                background-color: #f0f0f0;
                align-self: flex-start;
                visibility: hidden;
              "
            ),

            # main area wrapper
            div(
              id = "main_area",
              style = "
                position: relative;
                flex-grow: 1;
                margin-left: 10px;
                padding: 2px;
                background-color: transparent;
              ",
              div(
                id = "home_text_container",
                style = "
                  max-width: 1500px;
                  padding: 30px;
                  margin-left: -270px;
                  margin-top: -50px;
                ",
                htmlOutput("home_html")
              )
            )
          )
        ),

        # ----- Database -----
        tabPanel(
          title = tagList(icon("database"), "Database"),
          div(
            style = "display: flex; flex-direction: row;",

            # sidebar
            div(
              id = "sidebar_panel",
              style = "
                width: 300px;
                padding: 20px;
                background-color: #f0f0f0;
                align-self: flex-start;
              ",
              div(class = "title", HTML("MGA <br>（Mash Genome Atlas）")),

              radioButtons("species", "Species", choices = c("Human", "Mouse")),

              conditionalPanel(
                condition = "input.species == 'Human'",
                radioButtons(
                  "analysis_type",
                  "Select Analysis Type",
                  choices = c(
                    "Summary",
                    "Classification",
                    "NAFLD Score",
                    "Fibrosis Stages",
                    "Multi-Gene Analysis"
                  )
                )
              ),

              conditionalPanel(
                condition = "input.species == 'Mouse'",
                radioButtons(
                  "feed_type",
                  "Select Feed Type",
                  choices = c(
                    "HFD", "HFHS", "HFHC", "HFHSHC", "WDCCL4",
                    "STAM", "MCD", "CDAHFD"
                  )
                )
              ),

              conditionalPanel(
                condition = "input.species == 'Human' && input.analysis_type == 'Multi-Gene Analysis'",
                radioButtons(
                  "multi_gene_analysis",
                  "Select Type",
                  choices   = c("Find Correlated Genes", "Double Gene Correlation"),
                  selected  = "Find Correlated Genes"
                )
              ),

              conditionalPanel(
                condition =
                  "(input.species == 'Human' && (input.analysis_type != 'Multi-Gene Analysis' ||
                                                input.multi_gene_analysis == 'Find Correlated Genes')) ||
                   input.species == 'Mouse'",
                textInput(
                  "gene",
                  label = HTML(
                    'Enter Gene Name<br/><span style="font-size:smaller">(e.g., COL1A1, Col1a1)</span>'
                  ),
                  value = ""
                )
              ),

              conditionalPanel(
                condition = "input.species == 'Human' &&
                             input.analysis_type == 'Multi-Gene Analysis' &&
                             input.multi_gene_analysis == 'Double Gene Correlation'",
                textInput("gene1", "Enter Gene 1", value = ""),
                textInput("gene2", "Enter Gene 2", value = "")
              ),

              div(
                style = "margin-top:10px; margin-bottom:30px;",
                fluidRow(
                  column(
                    6,
                    actionButton("update", "Analyze", class = "btn-primary")
                  ),
                  column(
                    6,
                    div(
                      id = "spinner_area",
                      img(
                        src   = "ZKZg.gif",
                        style = "
                          width: 20px;     /* ❶ 缩小宽度，想多小写多小 */
                          height: 20px;    /* ❷ 保持比例也可用 height:auto */
                          border: none;    /* ❸ 去掉边框 */
                          box-shadow: none;
                        "
                      ),
                      style = "display:none; text-align:center;"  # ⬅️ 原来的隐藏 + 居中
                    )
                  )
                )
              )
            ),

            # main display
            div(
              id = "main_area",
              style = "
                position: relative;
                flex-grow: 1;
                margin-left: 10px;
                padding: 2px;
                background-color: transparent;
              ",
              uiOutput("usage_ui"),

              div(
                id = "plot_mask",
                style = "
                  position: absolute;
                  top: 0;
                  left: 0;
                  width: 100%;
                  height: 100%;
                  background: white;
                  z-index: 1000;
                  display: none;
                "
              ),

              div(
                style = "
                  max-width: 1500px;
                  padding: 15px;
                  margin-left: 0px;
                ",
                uiOutput("dynamic_plot")
              )
            )
          )
        ),

        # ----- Help -----
        tabPanel(
          title = tagList(icon("circle-question"), "Help"),
          div(
            style = "display: flex; flex-direction: row;",

            div(
              id = "sidebar_panel",
              style = "
                width: 300px;
                padding: 20px;
                background-color: #f0f0f0;
                align-self: flex-start;
                visibility: hidden;
              "
            ),

            div(
              id = "main_area",
              style = "
                position: relative;
                flex-grow: 1;
                margin-left: 10px;
                padding: 2px;
                background-color: transparent;
              ",
              div(
                id = "help_text_container",
                style = "
                  max-width: 1500px;
                  padding: 30px;
                  margin-left: -270px;
                  margin-top: -50px;
                ",
                uiOutput("help_text")
              )
            )
          )
        )
      )
    )
  )
)