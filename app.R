library(shiny)
library(shinydashboard)
library(tidyverse)
library(survival)
library(JMbayes2)
library(lubridate)
library(shinyjs)

source("IgANJM.R", encoding = "UTF-8")

## app.R ##
#load(".RData")
JM6A = read_rds("jointFit_clinical6A.RDS")

#FD = read_rds("FD.RDS")
ExampleFile = read_csv("du.csv")
sidebar <- function () {
  dashboardSidebar(
    disable = TRUE,
    sidebarMenu(
      id = "sidebar",
      # width  = 350,
      menuItem(
        "Dashboard",
        tabName = "dashboard",
        icon = icon("dashboard", verify_fa = FALSE)
      ),
      menuItem(
        "IgAN-JM",
        tabName = "IgANJM",
        icon = icon("book", verify_fa = FALSE)
      ),
      menuItem(
        "Source code",
        icon = icon("file-code-o", verify_fa = FALSE),
        href = "https://github.com/rstudio/shinydashboard/"
      ),
      textInput("pID", "patientId"),
      bookmarkButton()
    )
  )
}

body <-  function () {
  dashboardBody(fluidRow(useShinyjs(),
                         tabItems(
                           tabItem(
                             tabName = "dashboard",
                             h2("IgA肾病尿毒症风险预测", align = "center"),
                             h2("IgAN ESKD predictor", align = "center"),
                             box(
                               title = "背景与使用方法",
                               width = 12,
                               solidHeader = TRUE,
                               status = "primary",
                               tags$div(
                                 h4("背景：IgA肾病是最常见的肾小球肾炎。大概三分之一的病人会在20年内进入尿毒症。
                                 尿蛋白、血白蛋白和肾功能水平（血肌酐）是最重要的判断尿毒症风险的指标。
                                 我们联合使用多层线性模型和Cox模型（即联合模型Joint Modle）预测患者进入尿毒症的概率。
                                 相对于传统的横截面模型(Cox模型、AI模型等)，本模型能够利用患者的多次检验结果，降低预测误差。
                                 可以用于全病程任何时刻，不局限于肾活检时刻。本程序也是当前世界上第一个和唯一一个动态预测IgA肾病尿毒症风险的模型。"),
                                 h4("使用方法：请下载以下范例数据, 并将其中的数据更改为您自己的检验日期和结果。然后上传提交。推荐在电脑或笔记本上使用。
                                    微信中无法直接下载范例数据，请点击右上方在浏览器中使用。"),
                                  
                               ),
                               downloadLink("download", label = "范例数据下载(download the Example data)")
                             ),
                             box(
                               title = "数据输入(Data Input)",
                               width = 12,
                               solidHeader = TRUE,
                               status = "info",
                               radioButtons("sex", "第一步：请选择性别(Sex):", c("男 Male", "女 Female")),
                               sliderInput("age", "第二步：请选择年龄(Age):", 16, 80, 16),
                               strong("第三步：请下载范例数据,更改为自己的结果后上传(upload your data):"),
                               fileInput("upload",
                                 label = NULL,
                                 buttonLabel = "浏览Browse...",
                                 accept = c(".csv", ".tsv")
                               )
                             ),
                             box(
                               title = "数据预览(Data Preview)",
                               width = 12,
                               status = "primary",
                               solidHeader = TRUE,
                               tableOutput("files"),
                               ##数据提交前隐藏submit
                               shinyjs::hidden(
                                 actionButton("Submit", "提交Submit", class = "btn btn-primary btn-lg btn-block")
                               )
                             ),
                             box(
                               title = "结果(Result)",
                               width = 12,
                               solidHeader = TRUE,
                               status = "success",
                               hidden(
                                 plotOutput("IgANJM", height = "600px")
                                 )
                             ),
                             box(
                               title = "解读(Interpretation)",
                               width = 12,
                               solidHeader = TRUE,
                               status = "warning",
                               h4("本模型的预测结果是基于既往的临床经验，患者尿毒症预后取决于后期对治疗的反应和病情变化，仅供参考。"),
                               h4("本模型预测的起病开始到第10年是否进入尿毒症的准确性为(c-index) 为88% - 95%"),
                               h4("本模型结果由三幅图形构成，左侧两幅图分别为尿蛋白和肾功能(eGFR)随时间的变化,右侧大图为发生尿毒症概率风险预测。
                                  右图中红色实线是点估计，反应平均风险概率；浅红色为95%置信区间（即最坏情况和最好情况）。大多数情况下，预测的置信区间范围都是很大的。
                                  因为疾病的进展受很多因素影响"),
                               h4("中位肾脏生存期：进入尿毒症概率风险达到50%的时间，在这个时间点有50%的概率发生尿毒症。")
                             ),
                           )
                         )))
}

# Put them together into a dashboardPage
ui <- function(request) {
  dashboardPage(dashboardHeader(title = "IgAN", disable = TRUE),
                sidebar(),
                body())
}

server <- function(input, output, session) {
  ##下载范例数据
  output$download <- downloadHandler(
    filename = function() {
      "my.csv"
    },
    content = function(file) {
      write.csv(ExampleFile, file)
    }
  )
  ##上传数据处理
  data <- reactive({
    req(input$upload)
    ext <- tools::file_ext(input$upload$name)
    switch(
      ext,
      csv = vroom::vroom(input$upload$datapath, delim = ","),
      tsv = vroom::vroom(input$upload$datapath, delim = "\t"),
      validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })
  
  ##上传数据处理
  data_eGFR_age <- reactive({
    req(input$upload)
    data()  %>% mutate( patientId = as.factor(as.numeric(Sys.time())),
                        DATE_FU = as.Date(DATE_FU),
                        Age_Fu = input$age - as.numeric(as.Date(Sys.time()) - DATE_FU)/365,
                        CREA = if_else(CREA > 25, CREA/88.4,CREA),
                        SEX = if_else(input$sex =="男 Male",1,0),
                        eGFR = if_else(SEX == 1,
                                       142* (pmin(CREA/0.9,1)^-0.302) * (pmax(CREA/0.9,1)^-1.2) * 0.9938^Age_Fu,   #Male
                                       142* (pmin(CREA/0.7,1)^-0.241) * (pmax(CREA/0.7,1)^-1.2) * 0.9938^Age_Fu * 1.012),
                        PRO24H_log = log(PRO24H),
                        datefirst = min(DATE_FU),
                        ti = as.numeric(DATE_FU - datefirst)/365)
  }
  )
  
  output$files =  renderTable(data())
  #output$orders = renderDataTable(FD %>% filter (patientId == input$pID) %>% select(patientId,t,eGFR,PRO24H,ALB,UPCR))
  
  ##数据提交后显示submit
  observeEvent(input$upload, {
    shinyjs::show("Submit")
  })
  
  ##submit后开始计算
  observeEvent(input$Submit, {
    output$IgANJM = renderPlot(dynamicplotIgAN(
      JMmodel = JM6A,
      newD = isolate(data_eGFR_age()),
      years_pred = 10
    ))
    shinyjs::show("IgANJM")
  })
}



shinyApp(ui, server, enableBookmarking = "url")
