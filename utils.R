library(Seurat)
library(dplyr)
library(ggplot2)
library(plotly)
library(plyr)
library(dplyr)
library(varhandle)
library(reshape2)


##Helper calculation and data functions

#All plots were desigend around a width of 330 pixels, so we scale around that for different screen sizes
scaleRatio <- function(inputWidth) {
  return(inputWidth / 330)
}

PercentAbove <- function(x) {
  return(length(x = x[x > 0]) / length(x = x))
}

MaxMutate <- function(x) {
  return(x / max(x))
}

IsSeurat2 <- function() {
  return (packageVersion("Seurat") < 3)
}

get_shared_genes <- function(inputGeneList1, inputGeneList2, topN) {
  gene_list1 = dplyr::distinct(as.data.frame(inputGeneList1) %>% mutate_if(is.factor, as.character))
  gene_list2 = dplyr::distinct(as.data.frame(inputGeneList2) %>% mutate_if(is.factor, as.character))
  colnames(gene_list1) = c("gene")
  colnames(gene_list2) = c("gene")
  shared <- dplyr::semi_join(as.data.frame(gene_list1), as.data.frame(gene_list2), by = "gene")
  return(dplyr::top_n(shared, topN)$gene)
}

FetchGenes <- function(
  object,
  vars.all = NULL
) {
  if (IsSeurat2()) {
    return(FetchData(object, vars.all))
  }
  cells.use <- colnames(object)
  gene.check <- vars.all %in% rownames(GetAssayData(object))
  if (!all(gene.check)) {
    for (my.var in vars.all) {
      if (!(my.var %in% rownames(GetAssayData(object)))) {
        stop(paste("Error:", my.var, "not found"))
      }
    }
  }
  data.expression <- GetAssayData(object)
  data.expression <- t(as.matrix(data.expression[vars.all[gene.check], cells.use, drop=FALSE]))
  return(as.matrix(x = data.expression))
}



##Helper plotting functions

GetClusterPlot <- function(inputDataList, inputDataIndex, inputWidth, inputHeight) {

  inputDataObj = inputDataList[[inputDataIndex]]

  x_ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = TRUE,
    scaleanchor = 'y',
    scaleratio = inputDataObj$x_scale_ratio_clusterPlot
  )

  y_ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = TRUE
  #scaleratio = inputDataObj$y_scale_ratio_clusterPlot
  )

  #p <- plot_ly(inputDataObj$plot_df,source="plot_cluster",hoverinfo="skip",x=~dim1,y=~dim2,type="scattergl",width=300,mode="markers",marker=list(color=~colorVec,size=2)) %>%
  #  layout(
  #    #autosize = TRUE,
  #    title=inputDataObj$name,
  #    xaxis = x_ax,
  #    yaxis= y_ax
  #  ) %>% config(displayModeBar = F)

  p <- plot_ly(inputDataObj$plot_df, source = "plot_cluster", width = inputWidth, height = inputHeight) %>%
    add_trace(
      x = ~dim1,
      y = ~dim2,
      hoverinfo = "text",
      type = "scattergl",
      mode = "markers",
      text = ~info,#cell_description
      key = ~cluster,
      marker = list(size = 2 * scaleRatio(inputWidth) * inputDataObj$pt_size, color = ~colorVec),
      opacity = 0.5
    ) %>%
    
    add_annotations(
      x = inputDataObj$title_coords$x_center,
      y = inputDataObj$title_coords$y_center,
      text = sprintf("<b>%s</b>", inputDataObj$title_coords$cluster),
      showarrow = FALSE,
      font = list(size = 11 * scaleRatio(inputWidth) * inputDataObj$font_scale)
    ) %>%
    layout(
      autosize = TRUE,
      title = inputDataObj$name,
      xaxis = x_ax,
      yaxis = y_ax
    ) %>%
    hide_colorbar() %>%
    config(displayModeBar = F)

  return(p)
}

GetPlotData <- function(inputDataObj, inputGene, quantile=0.99) {
  data <- as.numeric(FetchGenes(inputDataObj$seurat_data, inputGene))
  cutoff <- quantile(x = data[data > 0], probs = quantile) 
  data[data > cutoff] <- cutoff
  single_gene <- mutate(inputDataObj$plot_df[, 1:3], inputDataObj$plot_df[,6], gene = data) %>% arrange(gene)#barcodes added
  colnames(single_gene) = c("dim1", "dim2","barcodes","info", "gene") #add cell metadata into the data frame for expression plot
  single_gene$infoExp<-str_c(single_gene$info, "\n", signif(single_gene$gene, digits = 4)) #add a column that combines cell information with exp value
  return(single_gene)
}

GetExpressionPlot <- function(inputDataList, inputDataIndex, inputGene, inputWidth, inputHeight) {

  inputDataObj = inputDataList[[inputDataIndex]]

  #On initialization, check if the inputGene is not defined
  if (inputGene == "") {
    return(NULL)
  }
  else {
    x_ax <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE,
      scaleanchor = 'y',
      scaleratio = inputDataObj$x_scale_ratio_clusterPlot
    )

    y_ax <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE
    #scaleratio = inputDataObj$y_scale_ratio_clusterPlot
    )

    single_gene <- GetPlotData(inputDataObj, inputGene)
    p <- plot_ly(
      single_gene, 
      source = "plot_expression", 
      x = ~dim1, 
      y = ~dim2, 
      type = "scattergl", 
      width = inputWidth, 
      height = inputHeight, 
      mode = "markers", 
      text = ~infoExp, #expression value with cell information
      color = ~gene, 
      marker = list(size = 2 * scaleRatio(inputWidth) * inputDataObj$pt_size), 
      hoverinfo = "text", 
      name = inputGene, 
      key=~barcodes,
      colors = c("grey90", "red") 
    )%>%
    
    layout(
    #autosize = TRUE,
      title = inputGene,
      xaxis = x_ax,
      yaxis = y_ax
    ) %>% 
    hide_colorbar() %>% 
    config(displayModeBar = F)

    return(p)
  }
}

GetDotPlot <- function(inputDataList, inputDataIndex, inputGeneList, inputWidth, inputHeight) {

  inputDataObj = inputDataList[[inputDataIndex]]

  if (length(inputGeneList) == 0) {
    return(NULL)
  }
  else {
    #Get the gene expression values and scale them so the max value for each gene is 1
    gene_exp = FetchGenes(inputDataObj$seurat_data, inputGeneList)
    #Combine the cluster assignments with the gene expression data
    multiple_genes <- as.data.frame(cbind(cluster = as.character(inputDataObj$plot_df$cluster), as.data.frame(gene_exp)))

    #Calculate the average expression per gene per cluster
    avgs <- multiple_genes %>% group_by(cluster) %>% dplyr::summarise_all(funs(mean))
    #Normalize so max is 1, melt the dataframe so we can plot it, and max sure the clusters are factors for proper plotting
    avgs <- melt(cbind(cluster = avgs$cluster, avgs %>% select(-cluster)), id.vars = c("cluster"))
    colnames(avgs) = c("cluster", "gene", "average_expression")
    avgs <- mutate(avgs, average_expression_relative = MaxMutate(average_expression))
    #avgs$cluster = as.factor(avgs$cluster)

    #Calculate the percent of cells that each gene was detected in per cluster
    p_above <- melt(multiple_genes %>% group_by(cluster) %>% dplyr::summarise_all(funs(PercentAbove)), id.vars = c("cluster"))

    #Combine the calculations
    combined = cbind(avgs, percent_above = 100 * p_above[, 3])
    combined = cbind(combined, percent_above_scaled = 100 * p_above[, 3] * scaleRatio(inputWidth))
    #Reverse the row order so it plots correctly - from https://stat.ethz.ch/pipermail/r-help/2008-September/175012.html
    rev_combined <- combined[rev(rownames(combined)),]
    #Add the hover text
    rev_combined <- mutate(rev_combined, hover_text = sprintf("Cluster: %s\nAvg. Expression: %0.3f\nPercent Cells: %0.2f%%", cluster, average_expression, percent_above))

    y_ax <- list(
      title = "",
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = TRUE,
      showgrid = FALSE,
      categoryorder = "trace"
    )

    x_ax <- list(
      title = "",
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = TRUE,
      showgrid = FALSE,
      categoryorder = "array",
      categoryarray = inputDataObj$category_order
    )

    t <- list(
      size = 12 * scaleRatio(inputWidth))

    #colorbar=list(title='Avg. expr.')
    p <- plot_ly(
      rev_combined, 
      source = "plot_dot", 
      x = ~cluster, 
      y = ~gene, 
      type = "scattergl", 
      mode = "markers", 
      width = inputWidth, 
      height = inputHeight, 
      text = ~hover_text, 
      hoverinfo = "text", 
      marker = list(symbol = "circle", size = rev_combined$percent_above_scaled, sizemode = "area", color = ~average_expression_relative)
    ) %>%
    layout(
      title = 'Dot Plot',
  #autosize = TRUE,
      showlegend = FALSE,
      yaxis = y_ax,
      xaxis = x_ax,
      font = t
    ) #%>% config(displayModeBar = F)

    return(p)
  }
}