library(shiny)
library(Seurat)
library(plotly)
library(plyr)
library(dplyr)
library(varhandle)
library(DT)
library(rlist)
library(logging)
library(scales)
library(stringr)
library(Matrix)
source("utils.R")

#Start to read in the config file.
json_file <- rjson::fromJSON(file = "~/Downloads/shiny_cell_browser-master/data/pbmc_config.json")

json_data <- json_file$data
datasets <- 1:length(json_data)
dataset_names <- sapply(json_data, function(x) x$name)
dataset_selector <- as.list(c(datasets))
names(dataset_selector) <- c(dataset_names)

#Use only the first dataset in the config file
# dataset_name = dataset_names[[1]]
# dataset = datasets[[1]]

#Read the config data
config <- json_file$config

IsSeurat2 <- function() {
  return (packageVersion("Seurat") < 3)
}

SetAllIdent <- function(object, ids) {
  if (IsSeurat2()) {
    return(Seurat::SetAllIdent(object, ids))
  }
  Idents(object) <- ids
  return(object)
}

GetClusters <- function(object) {
  if (IsSeurat2()) {
    return(Seurat::GetClusters(object))
  }
  clusters <- data.frame(cell.name = names(object@active.ident), cluster = as.character(object@active.ident))
  rownames(clusters) <- NULL
  clusters$cell.name <- as.character(clusters$cell.name)
  return(clusters)
}

GetDimReduction <- function(object, reduction.type = "umap", slot = "cell.embeddings") {
  if (IsSeurat2()) {
    return(Seurat::GetDimReduction(object, reduction.type = reduction.type, slot = slot))
  }
  reduction <- object[[reduction.type]]
  return(eval(expr = parse(text = paste0("reduction", "@", slot))))
}

GetCellNames <- function(object) {
  if (IsSeurat2()) {
    return(object@cell.names)
  } else {
    return(colnames(object))
  }
}

GetActiveIdent <- function(object) {
  if (IsSeurat2()) {
    return(object@ident)
  } else {
    return(object@active.ident)
  }
}

GetAssayData <- function(object) {
  if (IsSeurat2()) {
    return(object@data)
  } else {
    return(Seurat::GetAssayData(object))
  }
}

calc_pt_size <- function(n) { 25 / n ^ 0.33 }

#Now read in the data
read_data <- function(x) {
  # load data and metadata specified by the JSON string.
  # x: individual json string, with [name, file, clusters embedding]
  seurat_data <- readRDS(x$file)
  seurat_data <- SetAllIdent(seurat_data, x$cluster)
  ncells <- length(GetCellNames(seurat_data))
  
  pt_size <- calc_pt_size(ncells)
  if (!is.null(x$pt_size)) {
    pt_size <- x$pt_size
  }
  font_scale <- 1
  if (!is.null(x$font_scale)) {
    font_scale <- x$font_scale
  }
  colors <- seurat_data@misc[[sprintf("%s_colors", x$cluster)]]
  if (is.null(colors)) {
    set.seed(2)
    colors <- sample(hue_pal()(n_distinct(GetActiveIdent(seurat_data))))
  }
  genes <- sort(rownames(GetAssayData(seurat_data)))
  

  
  #Parser additions
  #full_embedding <- as.data.frame(GetDimReduction(seurat_data, reduction.type = x$embedding, slot = "cell.embeddings"))
  #assign_clust <- as.data.frame(GetClusters(seurat_data))
 # assign_cell<-as.data.frame(GetCellInfoMk2(seurat_data, x$cell_details))#cell type minor
  #seurat_data<-SetAllIdent(seurat_data, x$cluster)
  #if (is.factor(seurat_data@meta.data[[x$cluster]])) {
   # assign_clust[, 2] <- factor(assign_clust[, 2], levels = seurat_data@meta.data[[x$cluster]] %>% levels)
  #}
  #if(is.factor(seurat_data@meta.data[[x$cell_details]])){
  #  assign_clust[, 3]<-factor(assing_clust[,3], levels = seurat_data@meta.data[[x$cell_details]] %>% levels)
  #}
  
  #colorVec = mapvalues(as.integer(assign_clust[, 2]), from = 1:length(colors), to = toupper(colors)) #1:length(colors)
  #df_plot = cbind(full_embedding, assign_clust[, 2], colorVec) #added
  #colnames(df_plot) = c("dim1", "dim2", "cluster", "colorVec") #added
  #y_range = max(full_embedding[, 2]) - min(full_embedding[, 2])
  #x_domain = max(full_embedding[, 1]) - min(full_embedding[, 1])
  #xScaleRatio_clusterPlot = y_range / x_domain
  #yScaleRatio_clusterPlot = x_domain / y_range
  #coords_title = group_by(df_plot, cluster) %>% dplyr::summarize(x_center = mean(dim1), y_center = mean(dim2))
  #if (!is.null(x$label_coordinates)) {
   # coords_title <- dplyr::bind_rows(x$label_coordinates)
    #colnames(coords_title) <- c("cluster", "x_center", "y_center")
  #}
  
  #Parser additions (modified) fetch UMAP coordinate, cluster, detailed cell type, and barcodes
  #if(is.null(x$cell_details)){
   # seurat_data$temp<-"1"
  #}
  full_embedding <- as.data.frame(GetDimReduction(seurat_data, reduction.type = x$embedding, slot = "cell.embeddings"))
  clust<-data.frame(cell.name=names(seurat_data@active.ident), cluster=seurat_data@active.ident)
  seurat_data<-SetAllIdent(seurat_data, x$cell_details) #reset ident to retrieve detailed cell types
  cell<-data.frame(cell=seurat_data@active.ident) 
  #subtype<-data.frame(subtype=x$subtype)
  assign_clust<-cbind(clust, cell)
  #combine all necessary information for hoverinfo
  assign_clust$info<-str_c(assign_clust$cluster, "\n", assign_clust$cell, "\n", assign_clust$cell.name) 
  seurat_data<-SetAllIdent(seurat_data, x$cluster)

  
  #assign levels to cluster/cell information
  if (is.factor(seurat_data@meta.data[[x$cluster]])) {
    assign_clust[, 2] <- factor(assign_clust[, 2], levels = seurat_data@meta.data[[x$cluster]] %>% levels)
  }
  if(is.factor(seurat_data@meta.data[[x$cell_details]])){
    assign_clust[, 3]<-factor(assing_clust[,3], levels = seurat_data@meta.data[[x$cell_details]] %>% levels)
  }
  
  #create the data frame for plot on single cell browser
  colorVec = mapvalues(as.integer(assign_clust[, 2]), from = 1:length(colors), to = toupper(colors)) #1:length(colors)
  df_plot = cbind(full_embedding, assign_clust[,1:4], colorVec) #added
  colnames(df_plot) = c("dim1", "dim2", "barcodes", "cluster", "cell_details", "info", "colorVec") #barcodes added
  y_range = max(full_embedding[, 2]) - min(full_embedding[, 2])
  x_domain = max(full_embedding[, 1]) - min(full_embedding[, 1])
  xScaleRatio_clusterPlot = y_range / x_domain
  yScaleRatio_clusterPlot = x_domain / y_range
  coords_title = group_by(df_plot, cluster) %>% dplyr::summarize(x_center = mean(dim1), y_center = mean(dim2))
  if (!is.null(x$label_coordinates)) {
    coords_title <- dplyr::bind_rows(x$label_coordinates)
    colnames(coords_title) <- c("cluster", "x_center", "y_center")
  }
  
  #Add the full description name on mouse over
  #adding information on mouse over
  if (is.null(x$cluster_name_mapping)) {
    cluster_names <- GetActiveIdent(seurat_data) %>% levels()
    names(cluster_names) <- cluster_names
    x$cluster_name_mapping <- as.list(cluster_names)
  }
  #if (is.null(x$cell_desc)){
   # cell_details<-GetCellInfoMk2(seurat_data, x$cell_details) %>% levels()
    #names(cell_details)<-cell_details
    #x$cell_desc<-as.list(cell_details)
  #}
  
  desc_df = list.flatten(x$cluster_name_mapping)
  source_abbv = names(desc_df)
  dest_desc = as.character(list.flatten(x$cluster_name_mapping))
  
  #cell_df=list.flatten(x$cell_desc)
  #cell_src=names(cell_df)
  #cell_dest=as.character(list.flatten(x$cell_desc))
  
  df_plot$cluster_description = as.character(mapvalues(df_plot$cluster, from = source_abbv, to = dest_desc))
  #df_plot$cell_description=as.character(mapvalues(df_plot$cell_details, from= cell_src, to=cell_dest))
  
  #Differential expression data
  differential_expression = read.csv(file = x$diff_ex_file, header = TRUE, sep = ",")
  plot_tab <- differential_expression # %>% select(-c("id")) #%>% select(-c("id","cluster","is_max_pct","p_val","myAUC","power"))
  
  if (!is.null(x$diff_ex_cluster) && x$cluster != x$diff_ex_cluster) {
    seurat_data2 <- SetAllIdent(seurat_data, x$diff_ex_cluster)
    assign_clust2 <- as.data.frame(GetClusters(seurat_data2))
    merged = dplyr::left_join(assign_clust, assign_clust2, by = "cell.name")
    keyMap = distinct(merged %>% select(cluster.x, cluster.y))
    
    plot_tab$cluster = as.character(mapvalues(as.character(plot_tab$cluster), from = as.character(keyMap$cluster.y), to = as.character(keyMap$cluster.x)))
    seurat_data <- SetAllIdent(seurat_data, x$cluster)
  }
  
  return(
    list(
      name = x$name,
      seurat_data = seurat_data,
      ncells = ncells,
      pt_size = pt_size,
      font_scale = font_scale,
      embedding = x$embedding,
      colors = colors,
      genes = genes,
     # mtx = mtx,
      
      #Parser additions
      plot_df = df_plot,
      x_scale_ratio_clusterPlot = xScaleRatio_clusterPlot,
      y_scale_ratio_clusterPlot = yScaleRatio_clusterPlot,
      title_coords = coords_title,
      diff_eq_table = plot_tab,
      cluster_name_mapping = x$cluster_name_mapping
      #cell information mapping
      
    ))
}

# code to load all data (may slow down app startup)
# data_list <- lapply(json_data, read_data)

logging::loginfo("loading data...")
data_list <- rep(list(NULL), length(json_data))
data_list[[1]] <- read_data(json_data[[1]])
logging::loginfo("loaded dataset #1.")

#OLD WAY TO UPDATE EXPRESSION PLOT VIA PLOTLY UPDATE
#updateExpressionPlot <- function(input, output, session, inputGene)
#{
#  updateTextInput(session, "hidden_selected_gene", value = inputGene)

#new_plot_data = GetPlotData(organoid,inputGene)
#plotlyProxy("expression_plot", session) %>% plotlyProxyInvoke("addTraces",list(type="scattergl",mode="markers",hoverinfo="text",text=as.double(unlist(select(new_plot_data,"gene"))),marker=list(size=2,colors=c("grey90", "red"),color=as.double(unlist(select(new_plot_data,"gene")))),x=as.double(unlist(select(new_plot_data,"dim1"))),y=as.double(unlist(select(new_plot_data,"dim2")))))
#plotlyProxy("expression_plot", session) %>% plotlyProxyInvoke("deleteTraces",list(0))
#plotlyProxy("expression_plot", session) %>% plotlyProxyInvoke("relayout",list(title=inputGene))
#}

server <- function(input, output, session) {
  
  updateSelectInput(session, "selected_dataset", choices = dataset_names, selected = dataset_names[[1]])
  
  #Updates dataset index on selection and updates gene list
  current_dataset_index <- eventReactive({ input$selected_dataset }, {
    current_index <- dataset_selector[[input$selected_dataset]]
    return(current_index)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  observeEvent({ current_dataset_index() }, {
    current_index <- current_dataset_index()
    if (is.null(data_list[[current_index]])) {
      # Use <<- to modify global variable (shared across sessions)
      data_list[[current_index]] <<- read_data(json_data[[current_index]])
      logging::loginfo("loaded dataset #%s.", current_index)
    }
  }, ignoreInit = TRUE, ignoreNULL = TRUE, priority = 100)
  
  #Return current organoid and update values
  organoid <- eventReactive({ current_dataset_index() }, {
    return(data_list[[current_dataset_index()]])
  })
  
  #Update the gene list on change
  observeEvent({ organoid() }, {
    updateSelectizeInput(session, 'selected_gene', choices = organoid()$genes, server = TRUE)
  })
  
  #Logging
  observeEvent({ input$client }, {
    logging::loginfo("New client with ip: %s", input$client$ip)
  },
  ignoreNULL = TRUE, ignoreInit = TRUE)
  
  #Update expression plot on click
  observeEvent({
    s <- event_data("plotly_click", source = "plot_dot")
    !is.null(s$y)
  }, {
    s <- event_data("plotly_click", source = "plot_dot")
    updateTextInput(session, "hidden_selected_gene", value = s$y)
  },
  ignoreNULL = TRUE, ignoreInit = TRUE)
  

  
  #popup information on click
  observeEvent({s<-event_data("plotly_click", source="plot_expression")},
               {s<-event_data("plotly_click", source="plot_expression")
               #process the expression matrix
               inputDataObj = data_list[[current_dataset_index()]]
               mtx<-as.matrix(GetAssayData(inputDataObj$seurat_data))
               colIndex<-which(s$key == colnames(mtx))
               expList<-mtx[,colIndex]
               rankExp<-sort(expList,decreasing = TRUE)
               rankExp<-head(rankExp, 10)
               rankExp<-as.data.frame(rankExp)
                 showModal(
                 modalDialog(
                   title = paste0("Top 10 expressed genes in ", s$key, "(by expression value)"), 
                   as.list(rownames(rankExp))
                 )
                 #create the popup
                 )}
  )
  
  #Update expression plot from selectize input
  observeEvent({ input$selected_gene }, {
    updateTextInput(session, "hidden_selected_gene", value = input$selected_gene)
    logging::loginfo("Gene selection from text input: %s", input$selected_gene)
  },
  ignoreNULL = TRUE, ignoreInit = TRUE)
  
  #Update expression plot on table row click 
  observeEvent({ input$cluster_gene_table_rows_selected }, {
    rowid <- input$cluster_gene_table_rows_selected
    gene_selected <- current_table()[rowid, 'gene']
    updateTextInput(session, "hidden_selected_gene", value = gene_selected)
  },
  ignoreNULL = TRUE, ignoreInit = TRUE)
  
  #Get plot window width using the cluster plot as a reference
  plot_window_width = eventReactive({ session$clientData$output_cluster_plot_width }, {
    return(session$clientData$output_cluster_plot_width)
  })
  
  #Get plot window height using the cluster plot as a reference (force height = width)
  plot_window_height = eventReactive({ session$clientData$output_cluster_plot_width }, {
    return(session$clientData$output_cluster_plot_width)
  })
  
  #Generate the current table based on the current hidden selected cluster
  current_table <- eventReactive({
    input$hidden_selected_cluster
    current_dataset_index()
  }, {
    if (as.character(input$hidden_selected_cluster) == "") {
      return(organoid()$diff_eq_table)
    }
    else {
      subTable = filter(organoid()$diff_eq_table, cluster == input$hidden_selected_cluster)
      return(subTable)
    }
  })
  
  #Monitor cluster plot for changes and update hidden_selected_cluster field
  observeEvent({
    s <- event_data("plotly_click", source = "plot_cluster")
    !is.null(s)
  }, {
    s <- event_data("plotly_click", source = "plot_cluster")
    if (!is.null(s)) {
      updateTextInput(session, "hidden_selected_cluster", value = s$key)
    }
  })
  
  #Set the hidden_selected_cluster field to nothing when the reset button is clicked
  observeEvent(eventExpr = { input$reset_table }, handlerExpr = {
    updateTextInput(session, "hidden_selected_cluster", value = "")
  })
  
  #Set the hidden_selected_cluster field to nothing when the the dataset is changed 
  observeEvent(eventExpr = { current_dataset_index() }, handlerExpr = {
    updateTextInput(session, "hidden_selected_cluster", value = "")
  })
  
  #Update the gene table when current_table() changes
  observeEvent({ current_table() }, {
    dataTableProxy("cluster_gene_table", session, deferUntilFlush = TRUE) %>% replaceData(current_table(), rownames = FALSE)
  })
  
  #Update the dot plot with new gene list
  observeEvent(c({ input$gene_list_submit }, { current_dataset_index() }), {
    gene_listy <- trimws(strsplit(input$gene_list, '\n')[[1]])
    filtered_gene_list <- get_shared_genes(gene_listy, organoid()$genes, 10)
    updateTextAreaInput(session, "hidden_gene_list", value = paste(filtered_gene_list, collapse = ","))
  })
  
  current_gene_list <- eventReactive(c({ input$hidden_gene_list }, { current_dataset_index() }), {
    gene_listy = strsplit(paste(input$hidden_gene_list, collapse = ","), split = ",")[[1]]
    return(gene_listy)
  })
  
  ##GRAPHIC OUTPUTS
  output$cluster_plot <- renderPlotly({
    GetClusterPlot(data_list, current_dataset_index(), plot_window_width(), plot_window_height())
  }
  )
  output$expression_plot <- renderPlotly({
    GetExpressionPlot(data_list, current_dataset_index(), input$hidden_selected_gene, plot_window_width(), plot_window_height())
  }
  )
  output$dot_plot <- renderPlotly({
    GetDotPlot(data_list, current_dataset_index(), current_gene_list(), plot_window_width(), plot_window_height())
  }
  )
  
  clusterString <- eventReactive({ input$hidden_selected_cluster }, {
    baseString = "all clusters"
    if (input$hidden_selected_cluster != "") {
      baseString = organoid()$cluster_name_mapping[input$hidden_selected_cluster]
    }
    return(sprintf("Genes differentially expressed in %s", baseString))
  })
  
  #Export dot plot with new gene list
  observeEvent(c({ input$gene_list_submit }, { current_dataset_index() }), {
    gene_listy <- trimws(strsplit(input$gene_list, '\n')[[1]])
    filtered_gene_list <- get_shared_genes(gene_listy, organoid()$genes, 10)
    updateTextAreaInput(session, "hidden_gene_list", value = paste(filtered_gene_list, collapse = ","))
  })
  
  current_gene_list <- eventReactive(c({ input$hidden_gene_list }, { current_dataset_index() }), {
    gene_listy = strsplit(paste(input$hidden_gene_list, collapse = ","), split = ",")[[1]]
    return(gene_listy)
  })
  
  #Format the cluster gene table and add links to Addgene and ENSEMBL
  
  decimal_columns <- c('avg_logFC', 'p_val', 'p_val_adj', 'avg_diff')
  important_columns <- c('gene', 'cluster', 'p_val','p_val_adj')
  
  output$cluster_gene_table_title <- renderText({ clusterString() })
  output$cluster_gene_table <-
    DT::renderDT({
      datatable(organoid()$diff_eq_table,
                rownames = FALSE,
                extensions = c('Responsive'),
                selection = 'single',
                options =
                  list(
                    columnDefs =
                      list(
                        list(responsivePriority = 1, targets = important_columns),
                        list(
                          render = DT::JS(
                            "function(data, type, row, meta) {",
                            "return type === 'display'?",
                            "'<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=' + data + '\">' + data + '</a>' : data;",
                            "}"), targets = c(0)) #,
                      )
                  )
      ) %>%
        formatSignif(decimal_columns[decimal_columns %in% colnames(organoid()$diff_eq_table)], 3)
    },
    server = TRUE
    )
  #Download DEG Table
  output$TableExport<-downloadHandler(
    filename<-function() {
      paste("name","_DEG.csv",sep="")
      },
    content<-function(file){
     write.csv(organoid()$diff_eq_table,file) 
    }
    )
  
    #Download UMAP file
    output$UMAPExport<-downloadHandler(
      filename<-function() {
        paste("name","_UMAP.pdf",sep="")
      },
      content<-function(file){
        pdf(file)
        inputDataObj = data_list[[current_dataset_index()]]
        print(DimPlot(inputDataObj$seurat_data,reduction = "umap",label = TRUE,label.size = 6))
        dev.off()
        
      }
    )
    #Download Dotplot file
    output$DotplotExport<-downloadHandler(
      filename<-function() {
        paste("name","_Dotplot.pdf",sep="")
      },
      content<-function(file){
        pdf(file)
        inputDataObj = data_list[[current_dataset_index()]]
        print(DotPlot(inputDataObj$seurat_data, features = current_gene_list() # should not include
        ), dot.scale = 10) + 
          theme_classic() + coord_flip() + RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "white", high = "red")
        dev.off()      
        }
    )
  
}

