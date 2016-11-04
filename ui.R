library(shiny)
library(shinyjs)
library(shinythemes)
# library(shinyTable)
# library(shinyFiles)
# library(matrixStats)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(reshape2)

library(plotly)

panel_width = 3
result_with = 10

labelMandatory <- function(label) {
  tagList(label,span("*", class = "mandatory_star"))
}

appCSS <- ".mandatory_star { color: red; }"

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(#theme = shinytheme("cerulean"),
            # theme='cloneR.css',
    navbarPage("CloneR",

             tabPanel("Clone composition",

                      titlePanel("CloneR"),
                      h4("Assessing the clone compositon of tumours"),
                      shinyjs::useShinyjs(),
                      shinyjs::inlineCSS(appCSS),

                      sidebarLayout(
                      # fluidRow(
                        # SIDEBAR
                        # column(panel_width,
                        sidebarPanel(
                               # div(
                                 textInput("title",                 "Analysis", ""),
                                 textInput("outdir",                labelMandatory("Result directory"), ""),
                                 fileInput('sample_file',           labelMandatory('Samples'), accept = c('text/tab-separated-values', '.tsv' )),
                                 fileInput('mutation_file',         labelMandatory('Mutations'), accept = c('text/tab-separated-values', '.tsv' )),
                                 fileInput('cnv_file',              labelMandatory('CNV regions'), accept = c('text/tab-separated-values', '.tsv' )),
                                 # textInput("snp_folder",            labelMandatory("Path to germline heterozygous SNP folder"), ""),
                                 fileInput('list_gene_file',       'Genes of interest', accept = c('text/tab-separated-values', '.tsv' )),
                                 fileInput('gene_coordinate_file', 'Gene coordinates', accept = c('text/tab-separated-values', '.tsv' )),
                                 p(labelMandatory("Mandatory fields are marked with ")),
                                 actionButton("submit", "Run CloneR", class = "btn-primary")
                               # )
                               ),

                        # MAIN PANEL
                        # column(result_with,
                        mainPanel(
                               # shinyjs::hidden(
                                 # div(
                          tabsetPanel(
                            tabPanel("Clone Composition",      plotlyOutput("global_report")),
                            tabPanel('Summary',                DT::dataTableOutput('composition')),
                            tabPanel('Alteration Clonalities', DT::dataTableOutput('dataset'))
                          )
                                 # )
                               # )
                               )
                       )
             ),

             tabPanel("Help",
                      titlePanel("Help page")
                      )

    )
  )
)
