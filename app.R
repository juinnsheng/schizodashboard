#The RShiny dashboard contains 3 elements, ui, server and app initialization.

# This Youtube video comes in handy to install the packages : https://www.youtube.com/watch?v=vRr78s37CI4
#Use biocmanager to install the remaining packages
library(shiny)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(BiocManager)
options(repos = BiocManager::repositories())
# Load the organism database
org <- org.Hs.eg.db #Human database

# Read the data from a CSV file
df <- read.csv("data/schizo_genes.csv", header = TRUE) #header = TRUE indicates first row of csv contains the column names
new_colnames <- c("Gene", "Symbol", "Species", "Experiment_accession", "Comparison", "log2foldchange", "Adjustedpvalue")
colnames(df) <- new_colnames #Renaming the dataframe since there are a lot of blank spaces in the names of columns
df <- df[-1, ] #Removing the first row containing the old column names
df$Adjustedpvalue <- as.numeric(as.character(df$Adjustedpvalue), na.rm = TRUE) #Volcano plot takes in numeric values only

# Convert log2foldchange column to numeric with NA handling
df$log2foldchange <- as.numeric(as.character(df$log2foldchange)) #Volcano plot takes numeric values only
df$log2foldchange[is.na(df$log2foldchange)] <- NA 

# Define the UI
ui <- fluidPage(
  # Application title
  titlePanel("Differential Gene Expression Analysis for Schizophrenia"), #Setting the title page
  
  # Sidebar layout
  sidebarLayout(
    sidebarPanel(
      selectizeInput(
        inputId = "geneSelector", # Just specifying one of the inputs is associated with gene selection. This is important especially where there are multiple input elements. In server, I will use input$geneSelector to update its behaviour
        label = "Select Gene", # Label of the input widget 
        choices = unique(df$Symbol), #Choices from the dataframe consisting gene names from Ensembl
        multiple = FALSE, # If multiple selection is allowed. If FALSE, multiple selection is not allowed
        options = list(
          placeholder = "Select a gene...", #"Select a gene" text will be temporarily displayed when there is no input
          onInitialize = I('function() { this.setValue(""); }') #Setting the select input to be null. This is to allow to visualize the volcano plot for all the genes else only one gene's volcano plot can be viewed at one time
        )
      )
    ),
    
    mainPanel(
      # Output plots
      plotOutput("volcanoPlot"), #Creates a placeholder for the volcano plot to be displayed, to be accessed via Output$volcanoPlot in server side later with the code to show the volcanoplot
      plotOutput("dotPlot"), #Creates a placeholder for the dotplot to be displayed, to be accessed via Output$dotplot later in server side with the corresponding code to show the dotplot
      
      # Output references
      verbatimTextOutput("references") #creates a placeholder for displaying references, to be accessed via Output$references later with code to print the references
    )
  )
)

# Define the server logic
shinyServer <- function(input, output) {
  
  # Volcano plot
  output$volcanoPlot <- renderPlot({
    EnhancedVolcano(df, lab = df$Symbol, x = 'log2foldchange', y = 'Adjustedpvalue', pCutoff = 1e-4, FCcutoff = 1)
  }) #Code to plot the volcano plot based on the log2foldchange and p.adjusted
  
  # Dot plot
  output$dotPlot <- renderPlot({
    enrich_result <- enrichGO(gene = df$Gene,
                              OrgDb = org,
                              keyType = "ENSEMBL",
                              ont = "MF",
                              pvalueCutoff = 0.2,
                              qvalueCutoff = 0.2) #Code to conduct Gene Ontology Enrichment Analysis to obtain the molecular functions for each genes in the Ensemble database, a less stringent p-value is used to conduct gene ontology (GO) enrichment analysis 
    par(cex.lab = 0.5) #Minimizing the font size at y-axis or else it gets messy
    dotplot(enrich_result, showCategory = 15) #Creates a dotplot to study the top 15 molecular functions of genes that may associated with schizophrenia
  })
  
  # Render the references
  output$references <- renderPrint({
    cat("References:\n")
    cat("1. Hu, J., Xu, J., Pang, L., Zhao, H., Li, F., Deng, Y., Liu, L., Lan, Y., Zhang, X., Zhao, T., Xu, C., Xu, C., Xiao, Y., & Li, X. (2016). Systematically characterizing dysfunctional long intergenic non-coding RNAs in multiple brain regions of major psychosis. Oncotarget, 7(44), 71087-71098. [DOI: 10.18632/oncotarget.12122]\n")
    cat("2. Blighe, K., Rana, S., Turkes, E., Ostendorf, B., Grioni, A., Lewis, M. (2021). EnhancedVolcano (Version 1    .12.0). Bioconductor. [DOI: 10.18129/B9.bioc.EnhancedVolcano]\n")
    cat("3. Yu, G., Wang, L., Hu, E., Luo, X., Chen, M., Dall'Olio, G., Wei, W., Gao, C. (2021). clusterProfiler (Version 3.20.1). Bioconductor. [DOI: 10.18129/B9.bioc.clusterProfiler]\n")
  }) #A code to print the references, '\n' to push the text down
  
  # Produces the volcano plot and dot plot based on gene selected
  #observe event produces a response based on a conditional loop. If something happens, execute this for example
  #!is.null(input$geneSelector indicates that a value has been selected and it is not null,
  # input$geneSelector != "" acts as a validation check to ensure that the value is not null
  observeEvent(input$geneSelector, {
    if (!is.null(input$geneSelector) && input$geneSelector != "") { 
      selected_gene <- input$geneSelector
      gene_df <- df[df$Symbol == selected_gene, ]
      
      # Render the volcano plot for the selected gene
      output$volcanoPlot <- renderPlot({
        EnhancedVolcano(gene_df, lab = gene_df$Symbol, x = 'log2foldchange', y = 'Adjustedpvalue', pCutoff = 1e-4, FCcutoff = 1)
      })
      
      # Render the dot plot for the selected gene
      output$dotPlot <- renderPlot({
        enrich_result <- enrichGO(gene = gene_df$Gene,
                                  OrgDb = org,
                                  keyType = "ENSEMBL",
                                  ont = "MF",
                                  pvalueCutoff = 0.2,
                                  qvalueCutoff = 0.2)
        par(cex.lab = 0.5)
        dotplot(enrich_result, showCategory = 15) #All gene contains ENSEMBL identifiers
      })
    }
  })
}

# Run shinyApp
shinyApp(ui = ui, server = shinyServer)



