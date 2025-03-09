library(shiny)
library(SingleR)
library(celldex)
library(hdf5r)
library(Seurat)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(shinydashboard)
source("Global.R")
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = tags$span(
      style = "font-weight: bold; color: white;", 
      "FPCAM"
    ),
    titleWidth = 300
  ),
  dashboardSidebar(
    width = 250,
    sidebarMenu(
      id = 'tab',
      useShinyjs(),
      menuItem(
        "Home", 
        tabName = "home", 
        icon = icon("home"),
        selected = TRUE
      ),
      menuItem(
        "Cell Annotation", 
        tabName = "input", 
        icon = icon("microscope")
      ),
      menuItem(
        "UMAP Visualization", 
        tabName = "umap_", 
        icon = icon("eye")
      ),
      menuItem(
        "Score Heatmap", 
        tabName = "score_", 
        icon = icon("dna")
      ),
      menuItem(
        "Delta Distribution", 
        tabName = "delta_", 
        icon = icon("chart-bar")
      ),
      menuItem(
        "About", 
        tabName = "about", 
        icon = icon("info-circle")
      ),
      # 条件面板，用于文件上传区域
      conditionalPanel(
        condition = "input.tab == 'input'",
        div(
          style = "padding: 15px;",
          h4("File Upload"),
          fileInput(
            "file1", 
            label = tags$span(icon("file-excel"), "Upload File 1 (M)"),
            multiple = FALSE, 
            accept = c('.xlsx', '.xls'),
            width = "100%",
            placeholder = "Choose Excel file"
          ),
          fileInput(
            "file2", 
            label = tags$span(icon("file-excel"), "Upload File 2 (m)"),
            multiple = FALSE, 
            accept = c('.xlsx', '.xls'),
            width = "100%",
            placeholder = "Choose Excel file"
          ),
          div(
            style = "display: flex; flex-direction: column; gap: 10px; padding-top: 20px;",
            actionButton(
              "run", 
              "Run", 
              icon = icon("play"), 
              class = "btn-success btn-lg",
              style = "width: 75%; padding: 12px; border-radius: 5px;"
            ),
            actionButton(
              "reset", 
              "Reset", 
              icon = icon("undo"), 
              class = "btn-danger btn-lg",
              style = "width: 75%; padding: 12px; border-radius: 5px;"
            )
          )
        )
      )
    )
  ),
  dashboardBody(
    tabItems(
      # Home Tab
      tabItem(
        tabName = "home",
        box(
          width = 12,
          status = "primary",
          solidHeader = TRUE,
          title = "Welcome",
          tags$div(
            style = "text-align: center; padding: 20px;",
            # 添加图片
            tags$div(
              style = "display: flex; justify-content: center; align-items: center; margin: 30px 0;",
              tags$img(
                src = "cell_analysis.png",  # 图片需要放在 www 文件夹中
                width = "100%",
                style = "max-width: 1100px; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.1);"
              )
            ),
            tags$p(
              "Upload your Excel files and perform cell annotation analysis with ease.",
              style = "font-size: 16px; text-align: center;"
            )
          )
        )
      ),
      
      # Cell Annotation Tab
      tabItem(
        tabName = "input",
        tabsetPanel(
          id = 'main_tabs',
          tabPanel(
            "Instructions", 
            includeMarkdown("./markdown/instructions.md")
          ),
          tabPanel(
            "Result",
            tableOutput("resultTable"),
            downloadButton("downloadData", "Download Results")
          )
        )
      ),
      
      # Single Cell Analysis Tab
      tabItem(
        tabName = "umap_",
        fluidRow(
          column(
            12,
            box(
              width = NULL,
              status = "primary",
              title = "Analysis Input",
              fileInput(
                "h5_file", 
                "Choose H5 File",
                accept = c(".h5")
              ),
              actionButton(
                "run_umap_analysis", 
                "Run UMAP Analysis",
                class = "btn-primary",
                style = "width: 100%"
              ),
              # 显式添加下载按钮
              downloadButton("download_umap", "Download UMAP Plot"),
              hr(),
              helpText("Upload a 10X H5 file to start the analysis")
            )
          )
        ),
        fluidRow(
          column(
            12,
            tabBox(
              width = NULL,
              tabPanel(
                "UMAP Plot",
                plotOutput("umap_plot", height = "600px")
              )
            )
          )
        )
      ),
      tabItem(
        tabName = "score_",
        fluidRow(
          column(
            12,
            box(
              width = NULL,
              status = "primary",
              title = "Analysis Input",
              fileInput(
                "h5_file", 
                "Choose H5 File",
                accept = c(".h5")
              ),
              actionButton(
                "run_heatmap_analysis", 
                "Run Heatmap Analysis",
                class = "btn-primary",
                style = "width: 100%"
              ),
              # 显式添加下载按钮
              downloadButton("download_heatmap", "Download Heatmap"),
              hr(),
              helpText("Upload a 10X H5 file to start the analysis")
            )
          )
        ),
        fluidRow(
          column(
            12,
            tabBox(
              width = NULL,
              tabPanel(
                "Score Heatmap",
                plotOutput("score_heatmap", height = "600px")
              )
            )
          )
        )
      ),
      tabItem(
        tabName = "delta_",
        fluidRow(
          column(
            12,
            box(
              width = NULL,
              status = "primary",
              title = "Analysis Input",
              fileInput(
                "h5_file", 
                "Choose H5 File",
                accept = c(".h5")
              ),
              actionButton(
                "run_delta_analysis", 
                "Run Delta Analysis",
                class = "btn-primary",
                style = "width: 100%"
              ),
              # 显式添加下载按钮
              downloadButton("download_delta", "Download Delta Plot"),
              hr(),
              helpText("Upload a 10X H5 file to start the analysis")
            )
          )
        ),
        fluidRow(
          column(
            12,
            tabBox(
              width = NULL,
              tabPanel(
                "Delta Distribution",
                plotOutput("delta_plot", height = "600px")
              )
            )
          )
        )
      ),
      # About Tab 标签页
      tabItem(
        tabName = "about",
        fluidRow(
          column(
            12,
            box(
              width = NULL,
              status = "info",
              title = "About This Tool",
              tags$div(
                style = "text-align: center; padding: 20px;",
                tags$h3("FPCAM - Cell Annotation Tool"),
                tags$p(
                  "This tool provides a streamlined approach to annotating cells in single-cell RNA sequencing data. ",
                  style = "font-size: 16px; text-align: center;"
                ),
                tags$p(
                  "Developed by: Guo Jun Liu、Yan Shi",
                  style = "font-size: 14px; text-align: center;"
                ),
                tags$p(
                  "Contact: gjliu77@qq.com、1490261148@qq.com",
                  style = "font-size: 14px; text-align: center;"
                )
              )
            )
          )
        )
      )
    )
  )
)
# Server 函数
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 300 * 1024^2)  # 设置最大上传文件大小
  
  temp_folder <- tempdir()  # 定义临时文件夹路径
  cat("临时文件夹路径:", temp_folder, "\n")  # 打印路径以便调试
  
  # 定义响应式变量以存储 result_data
  result_data <- reactiveVal(NULL)
  
  
  observeEvent(input$run, {
    req(input$file1, input$file2)  # 确保两个文件都已上传
    
    # 显示进度条
    withProgress(message = 'Processing Files...', value = 0, {
      
      # 将文件保存到临时文件夹中，使用原始文件名
      file1_path <- file.path(temp_folder, input$file1$name)
      file2_path <- file.path(temp_folder, input$file2$name)
      
      # 复制上传的文件到临时文件夹
      file.copy(input$file1$datapath, file1_path, overwrite = TRUE)
      file.copy(input$file2$datapath, file2_path, overwrite = TRUE)
      
      # 检查文件是否成功复制
      if (file.exists(file1_path) && file.exists(file2_path)) {
        cat("文件已成功复制到临时文件夹。\n")
        print(list.files(temp_folder))  # 列出临时文件夹中的文件
        
        # 设置临时文件夹为工作目录（可选）
        setwd(temp_folder)
        rm(list = ls())
        Marker <- read_excel(file.path(temp_folder, "Marker_expression_test.xlsx"), 
                             sheet = 1, col_names = TRUE, col_types = NULL, na = "", skip = 0)
        # 初始化一个向量来存储每列的比值
        ratios <- list()
        # 循环遍历偶数索引的列
        sum(!is.na(Marker[, 3]))
        for (j in 1:sum(!is.na(Marker[, 3]))) {
          print(j)
        }
        for (i in seq(2, ncol(Marker), 2)) {
          # 计算列的总和
          column_sum <- sum(Marker[, i], na.rm = TRUE)
          # 计算每个值与列总和的比值
          column_ratios <- Marker[, i] / column_sum
          # 将比值向量添加到ratios列表中
          ratios[[i]] <- column_ratios
        }
        
        
        # ### 手动验证
        # col_name <- paste0("cluster","_", i/2-1, "_", "Markergene")
        # gene <- Marker[2,col_name]
        # column_sum <- sum(Marker[, 58], na.rm = TRUE)
        # sum(Marker$cluster_28_Express)
        # column_ratios <- Marker[, 58] / column_sum
        ###################################################################
        ## Script name:  计算列和基因列整合
        ###################################################################
        for (i in seq(1, ncol(Marker), 2)) {
          print(i)
          for (j in seq(2, ncol(Marker), 2)) {
            Express_var <- paste0("cluster","_", j/2-1, "_", "Express")
            # eval(Express_var)
            ratios[[j]][[eval(Express_var)]]
            # print(paste(i,j))
            # Marker[,i]
          }
        }
        # Marker[,1]
        # Express_var <- paste0("cluster","_", 4/2-1, "_", "Express")
        # eval(Express_var)
        # ratios[[4]][[eval(Express_var)]]
        # combined_list <- list(Marker[,1], ratios[[2]][[eval(Express_var)]])
        # combined_list <- list(Marker[,3], ratios[[4]][[eval(Express_var)]])
        # combined_list <- do.call(c, list(Marker[,1], ratios[[2]][[eval(Express_var)]]))
        # 初始化一个空列表
        inner_list <- list()
        # 外层循环
        for (i in seq(0, ncol(Marker)/2, 2)) {
          # 初始化内层列表
          inner_list <- list()
          # 内层循环
          for (j in seq(2, ncol(Marker), 2)) {
            # 为当前位置生成一个值
            Express_var <- paste0("cluster","_", j/2-1, "_", "Express")
            eval(Express_var)
            var <- list(ratios[[j]][[eval(Express_var)]])
            combined_list <- list(Marker[,j-1], ratios[[j]][[eval(Express_var)]])
            # combined_list <- paste("Item", i, "Subitem", j, sep = "_")
            # 将值添加到内层列表
            inner_list[[j]] <- combined_list
          }
          # 将内层列表添加到多维列表中
          # multi_dimensional_list[[i]] <- inner_list
        }
        
        
        # 测试从列表打印基因和表达量
        # print(inner_list)
        # inner_list[[2]][[1]][[1]][1]
        # inner_list[[2]][[1]][[1]]
        # # count the length of
        # inner_list[[2]][[2]][1]
        # inner_list[[58]][[1]][[1]][1]
        # inner_list[[58]][[2]][1]
        ###################################################################
        ## Script name: 比对
        ###################################################################
        db <- read_excel(file.path(temp_folder, "marker_database.xlsx"), 
                         sheet = 1, col_names = TRUE, col_types = NULL, na = "", skip = 0)
        a <- unlist(unique(db$cell_cluster))
        # 先创建列再创建行
        for (j in seq(2, ncol(Marker), 2)) {
          num =j/2-1
          tmp <- paste("Cluster","_",num,sep="")
          db[tmp] <- 0
          print(paste("Calculating",tmp))
          # j = 56
          my_list <- inner_list[[j]][[1]][[1]]
          for(m in 1:nrow(db)){
            # db[m,tmp] <- NA
            for (l in 1:length(na.omit(unlist(my_list)))) { # Calculate the length of list
              # print("The row",m)
              if (toupper(inner_list[[j]][[1]][[1]][l]) == toupper(db$marker_gene[m])) {
                db[m,tmp] <- inner_list[[j]][[2]][l]*db$percent[m]
              }
            }
          }
        }
        write.csv(db, file = file.path(temp_folder, "Cluster.csv"), row.names = FALSE)
        ###################################################################
        ## Script name: Calculate the cellular origin of clusters.
        ###################################################################
        result <- db %>%
          group_by(cell_cluster) %>%
          summarise(Total=sum(Cluster_1)) # set Cluster_0
        tmp <- "Cluster_1" #
        max_index <- which.max(result$Total)
        a<- result[max_index,1]
        print(paste0("The cellular origin of ",tmp," is ",a))
        
        ###################################################################
        ## Script name: 
        ###################################################################
        for (i in seq(2, ncol(Marker), 2)) {
          num =i/2-1
          num_2 = num+2
          tmp <- paste("Cluster",num,sep="_")
          # db[,tmp] <- as.numeric(db[,tmp])
          result <- db %>%
            group_by(cell_cluster) %>%
            summarise(across(starts_with("Cluster"),~sum(.,na.rm = TRUE)))
          max_index <- which.max(result[[num_2]])
          a<- result[max_index,1]
          print(paste0("The cellular origin of ",tmp," is ",a))
        }
        write.csv(result, file = file.path(temp_folder, "result.csv"), row.names = FALSE, quote = TRUE)
        
        ###################################################################
        ## Script name: 挑出最高和次高值
        ###################################################################
        # 读取CSV文件
        data <- read.csv(file.path(temp_folder, "result.csv"))
        # 初始化一个空的数据框，用于存储每一列的最大值、次高值及对应的细胞群
        result <- data.frame(Cluster = character(), Max_Cell_Cluster = character(), Max_Value = numeric(),
                             Second_Max_Cell_Cluster = character(), Second_Max_Value = numeric(), stringsAsFactors = FALSE)
        # 遍历每一列（从第二列到最后一列），找到每列的最大值、次高值及对应的细胞群
        for (j in 2:ncol(data)) {
          # 找到当前列的最大值
          max_value <- max(data[, j])
          # 找到最大值对应的细胞群
          max_cell_cluster <- data$cell_cluster[which.max(data[, j])]
          # 找到次高值（去除最大值后再找）
          temp <- data[, j]
          temp[which.max(temp)] <- -Inf  # 将最大值替换为负无穷大，确保不会再次被选中
          second_max_value <- max(temp)
          # 找到次高值对应的细胞群
          second_max_cell_cluster <- data$cell_cluster[which.max(temp)]
          # 将簇名、最大值对应的细胞群、次高值及其对应的细胞群存入结果数据框
          result <- rbind(result, data.frame(Cluster = names(data)[j], 
                                             Max_Cell_Cluster = max_cell_cluster, Max_Value = max_value, 
                                             Second_Max_Cell_Cluster = second_max_cell_cluster, Second_Max_Value = second_max_value))
        }
        # 将结果写入CSV文件
        write.csv(result, file = file.path(temp_folder, "result_max_second_values_per_column.csv"), row.names = FALSE)
        
        ###################################################################
        ## Script name: 计算delta制表
        ###################################################################
        # 读取CSV文件
        data <- read.csv(file.path(temp_folder, "result.csv"))
        # 读取result_max_second_values_per_column.csv文件
        max_data<- read.csv(file.path(temp_folder, "result_max_second_values_per_column.csv"), 
                            header = TRUE, stringsAsFactors = FALSE)
        # 创建一个向量用于存储每列的 delta1 值，长度为数据列数减1(去掉第一列)
        delta1_values <- numeric(ncol(data) - 1)
        # 遍历数据框中除第一列外的所有列
        for (i in 2:ncol(data)) {
          # 提取当前列的数据，drop=FALSE 保持数据框结构
          column_data <- data[, i, drop = FALSE]
          # 过滤数据：
          column_data <- column_data[!is.na(column_data) & column_data != max(column_data, na.rm = TRUE)]
          # 检查过滤后的数据是否还有值
          if (length(column_data) > 0) {
            # 计算过滤后数据的中位数
            median_value <- median(column_data, na.rm = TRUE)
            # 计算 delta1：原始数据的最大值减去过滤后数据的中位数
            delta1 <- max(data[, i], na.rm = TRUE) - median_value
            # 将计算得到的 delta1 存储在结果向量中
            # i-1 是因为我们的结果向量从1开始，而循环从2开始
            delta1_values[i - 1] <- delta1
          } else {
            delta1_values[i - 1] <- NA
          }
        }
        # 输出delta1值
        delta1_values
        # 第二部分：计算 λ2 (delta2)，加入阈值判断
        delta2_values <- numeric(ncol(data) - 1)
        for (i in 2:ncol(data)) {
          # 提取当前列的数据并去掉缺失值
          column_data <- na.omit(data[, i])
          
          # 计算最大值和次最大值
          sorted_values <- sort(unique(column_data), decreasing = TRUE)
          max_value <- sorted_values[1]
          second_max_value <- ifelse(length(sorted_values) > 1, sorted_values[2], sorted_values[1])
          
          # 检查最大值与次最大值的差是否大于等于0.00001
          if (max_value - second_max_value >= 0.00001) {
            # 只去掉次最大值
            column_data <- column_data[column_data != second_max_value]
            
            # 检查剩余数据是否足够计算中位数
            if (length(column_data) > 0) {
              median_value <- median(column_data)
              delta2 <- second_max_value - median_value
              delta2_values[i - 1] <- delta2
            } else {
              delta2_values[i - 1] <- NA
            }
          } else {
            # 如果差值不大于0.00001，存储NA
            delta2_values[i - 1] <- NA
          }
        }
        # 输出delta2值
        delta2_values
        # 检查长度匹配并更新max_data
        length_delta1 <- length(delta1_values)
        length_delta2 <- length(delta2_values)
        num_rows_max_data <- nrow(max_data)
        
        # 添加错误处理和输出长度信息
        print(paste("delta1_values length:", length_delta1))
        print(paste("delta2_values length:", length_delta2))
        print(paste("max_data rows:", num_rows_max_data))
        
        if (length_delta1 == num_rows_max_data && length_delta2 == num_rows_max_data) {
          # 将delta值添加到max_data中
          max_data$Delta1_Value <- delta1_values
          max_data$Delta2_Value <- delta2_values
          
          # 选择需要的列并创建新的数据框
          result_data <- max_data[, c("Cluster", "Max_Cell_Cluster", "Delta1_Value", 
                                      "Second_Max_Cell_Cluster", "Delta2_Value")]
          # 将结果写入一个新的CSV文件
          write.csv(result_data, file = file.path(temp_folder, "result_updated.csv"), row.names = FALSE)
          result_updated0<-read.csv(file.path(temp_folder, "result_updated.csv"))
          
          output$resultTable <- renderTable({
            result_updated0[] <- lapply(result_updated0, function(x) if(is.numeric(x)) format(x, scientific = FALSE) else x)
            result_updated0
          })                   
          
          
          
          
        } else {
          cat("文件复制失败，请检查。\n")
        }
      }
    })
    
    incProgress(1 / total_steps)
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("result_updated", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        # 确保文件路径正确
        temp_file <- file.path(temp_folder, "result_updated.csv")
        if (file.exists(temp_file)) {
          file.copy(temp_file, file)
        } else {
          stop("The result file does not exist.")
        }
      }
    )
    
    # 计算完成后关闭进度条，输出“计算完成”
    shinyjs::alert("Calculation Finished!")
  })
  
  
  # 设置上传文件大小限制
  options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB
  
  # 每个标签页的反应式值
  umap_values <- reactiveValues(seurat_obj = NULL)
  heatmap_values <- reactiveValues(singleR_pred = NULL)
  delta_values <- reactiveValues(singleR_pred = NULL)
  
  # UMAP标签页数据处理
  process_umap_data <- eventReactive(input$run_umap_analysis, {
    req(input$h5_file)  # 确保文件上传后才处理数据
    
    withProgress(message = 'Processing data...', value = 0, {
      # Read data
      incProgress(0.1, detail = "Reading H5 file")
      hdf5_obj <- Read10X_h5(filename = input$h5_file$datapath, use.names = TRUE, unique.features = TRUE)
      
      # Create Seurat object and preprocess
      incProgress(0.2, detail = "Creating Seurat object")
      pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)
      
      # QC and filtering
      incProgress(0.3, detail = "QC and filtering")
      pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
      pbmc.seurat.filtered <- subset(pbmc.seurat, 
                                     subset = nCount_RNA > 800 &
                                       nFeature_RNA > 500 &
                                       mitoPercent < 10)
      
      # Normalization and dimensionality reduction
      incProgress(0.4, detail = "Normalizing data")
      pbmc.seurat.filtered <- NormalizeData(pbmc.seurat.filtered)
      pbmc.seurat.filtered <- FindVariableFeatures(pbmc.seurat.filtered, selection.method = "vst", nfeatures = 2000)
      
      incProgress(0.5, detail = "Scaling data")
      pbmc.seurat.filtered <- ScaleData(pbmc.seurat.filtered, features = rownames(pbmc.seurat.filtered))
      
      incProgress(0.6, detail = "Running PCA")
      pbmc.seurat.filtered <- RunPCA(pbmc.seurat.filtered, features = VariableFeatures(pbmc.seurat.filtered))
      
      incProgress(0.7, detail = "Finding neighbors")
      pbmc.seurat.filtered <- FindNeighbors(pbmc.seurat.filtered, dims = 1:20)
      pbmc.seurat.filtered <- FindClusters(pbmc.seurat.filtered, resolution = 0.5)
      
      incProgress(0.8, detail = "Running UMAP")
      pbmc.seurat.filtered <- RunUMAP(pbmc.seurat.filtered, dims = 1:20)
      
      # Get reference data and run SingleR
      incProgress(0.9, detail = "Running SingleR")
      ref <- celldex::HumanPrimaryCellAtlasData()
      ref <- ref[,grepl('DC|B_cell|Neutrophils|T_cells|Monocyte|Erythroblast|Macrophage|NK_cell|Platelets|Myelocyte', ref$label.main)]
      
      pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')
      pred <- SingleR(test = pbmc_counts, ref = ref, labels = ref$label.main)
      
      # Add SingleR labels to Seurat object
      pbmc.seurat.filtered$singleR.labels <- pred$labels[match(rownames(pbmc.seurat.filtered@meta.data), rownames(pred))]
      
      # Store results
      umap_values$seurat_obj <- pbmc.seurat.filtered
      incProgress(1, detail = "Done!")
    })
  })
  
  # 热图标签页数据处理
  process_heatmap_data <- eventReactive(input$run_heatmap_analysis, {
    req(input$h5_file)  # 确保文件上传后才处理数据
    
    withProgress(message = 'Processing data...', value = 0, {
      # Read data
      incProgress(0.1, detail = "Reading H5 file")
      hdf5_obj <- Read10X_h5(filename = input$h5_file$datapath, use.names = TRUE, unique.features = TRUE)
      
      # Create Seurat object and preprocess
      incProgress(0.2, detail = "Creating Seurat object")
      pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)
      
      # QC and filtering
      incProgress(0.3, detail = "QC and filtering")
      pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
      pbmc.seurat.filtered <- subset(pbmc.seurat, 
                                     subset = nCount_RNA > 800 &
                                       nFeature_RNA > 500 &
                                       mitoPercent < 10)
      
      # Normalization
      incProgress(0.4, detail = "Normalizing data")
      pbmc.seurat.filtered <- NormalizeData(pbmc.seurat.filtered)
      
      # Get reference data and run SingleR
      incProgress(0.8, detail = "Running SingleR")
      ref <- celldex::HumanPrimaryCellAtlasData()
      ref <- ref[,grepl('DC|B_cell|Neutrophils|T_cells|Monocyte|Erythroblast|Macrophage|NK_cell|Platelets|Myelocyte', ref$label.main)]
      
      pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')
      pred <- SingleR(test = pbmc_counts, ref = ref, labels = ref$label.main)
      
      # Store results
      heatmap_values$singleR_pred <- pred
      incProgress(1, detail = "Done!")
    })
  })
  
  # Delta标签页数据处理
  process_delta_data <- eventReactive(input$run_delta_analysis, {
    req(input$h5_file)  # 确保文件上传后才处理数据
    
    withProgress(message = 'Processing data...', value = 0, {
      # Read data
      incProgress(0.1, detail = "Reading H5 file")
      hdf5_obj <- Read10X_h5(filename = input$h5_file$datapath, use.names = TRUE, unique.features = TRUE)
      
      # Create Seurat object and preprocess
      incProgress(0.2, detail = "Creating Seurat object")
      pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)
      
      # QC and filtering
      incProgress(0.3, detail = "QC and filtering")
      pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
      pbmc.seurat.filtered <- subset(pbmc.seurat, 
                                     subset = nCount_RNA > 800 &
                                       nFeature_RNA > 500 &
                                       mitoPercent < 10)
      
      # Normalization
      incProgress(0.4, detail = "Normalizing data")
      pbmc.seurat.filtered <- NormalizeData(pbmc.seurat.filtered)
      
      # Get reference data and run SingleR
      incProgress(0.8, detail = "Running SingleR")
      ref <- celldex::HumanPrimaryCellAtlasData()
      ref <- ref[,grepl('DC|B_cell|Neutrophils|T_cells|Monocyte|Erythroblast|Macrophage|NK_cell|Platelets|Myelocyte', ref$label.main)]
      
      pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')
      pred <- SingleR(test = pbmc_counts, ref = ref, labels = ref$label.main)
      
      # Store results
      delta_values$singleR_pred <- pred
      incProgress(1, detail = "Done!")
    })
  })
  
  # 生成UMAP图
  output$umap_plot <- renderPlot({
    req(umap_values$seurat_obj)
    DimPlot(umap_values$seurat_obj, 
            reduction = 'umap', 
            group.by = 'singleR.labels',
            label = TRUE,
            label.size = 4,
            repel = TRUE) +
      theme_classic() +
      ggtitle("UMAP visualization with SingleR annotations")
  })
  
  # 生成得分热图
  output$score_heatmap <- renderPlot({
    req(heatmap_values$singleR_pred)
    plotScoreHeatmap(heatmap_values$singleR_pred)
  })
  
  # 生成Delta分布图
  output$delta_plot <- renderPlot({
    req(delta_values$singleR_pred)
    plotDeltaDistribution(delta_values$singleR_pred)
  })
  
  # UMAP图下载功能
  output$download_umap <- downloadHandler(
    filename = function() {
      paste("umap-plot-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      req(umap_values$seurat_obj)
      
      # 创建UMAP图
      p <- DimPlot(umap_values$seurat_obj, 
                   reduction = 'umap', 
                   group.by = 'singleR.labels',
                   label = TRUE,
                   label.size = 4,
                   repel = TRUE) +
        theme_classic() +
        ggtitle("UMAP visualization with SingleR annotations")
      
      # 使用ggsave保存图像
      ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
    }
  )
  
  # 热图下载功能
  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste("score-heatmap-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      req(heatmap_values$singleR_pred)
      
      # 创建热图图
      p_heatmap <- plotScoreHeatmap(heatmap_values$singleR_pred)
      
      # 使用ggsave保存图像
      ggsave(file, plot = p_heatmap, width = 10, height = 8, dpi = 300)
    }
  )
  # Delta图下载功能
  output$download_delta <- downloadHandler(
    filename = function() {
      paste("delta-distribution-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      req(delta_values$singleR_pred)
      
      # 打开图形设备
      png(file, width = 3000, height = 2400, res = 300)
      
      # 绘制Delta分布图
      plotDeltaDistribution(delta_values$singleR_pred)
      
      # 关闭图形设备
      dev.off()
    }
  )
  
  # 触发数据处理
  observeEvent(input$run_umap_analysis, {
    process_umap_data()
  })
  
  observeEvent(input$run_heatmap_analysis, {
    process_heatmap_data()
  })
  
  observeEvent(input$run_delta_analysis, {
    process_delta_data()
  })
}



shinyApp(ui, server)