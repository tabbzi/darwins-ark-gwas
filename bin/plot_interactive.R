library(shiny)
library(tidyverse)

# Load the GWAS data from a .tsv file
gwas_data <- read_tsv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-gwas/test.mlma") %>%
  mutate(phe = "Q121",
         n = 2115) %>%
  select(phe,
         chr = Chr,
         pos = bp,
         snp = SNP,
         ref = A1,
         alt = A2,
         freq = Freq,
         p,
         b,
         n)

# Define the user interface
ui <- fluidPage(
  
  # Add a sidebar
  sidebarLayout(
    sidebarPanel(
      
      # Input for y-axis
      selectInput("yaxis", "Y-Axis",
                  choices = c("significance" = "p",
                              "effect" = "b"),
                  selected = "p")),
    
    # Add the main panel
    mainPanel(
      
      # Output for the GWAS manhattan plot
      plotOutput("manhattan"),
      
      # Output for the GWAS data table
      tableOutput("gwas")
    )
  )
)

# Define the server logic
server <- function(input, output) {
  
  # Filter the GWAS data based on user input
  filtered_data <- reactive({
    data <- gwas_data
    return(data)
  })
  
  # Create the GWAS manhattan plot
  output$manhattan <- renderPlot({
    data <- filtered_data()
    
    # Create a new variable for the cumulative position per chromosome
    data$cumulative_pos <- ave(data$pos, data$chr, FUN = cumsum)
    
    ggplot(data, aes(x = cumulative_pos,
                     y = -log10(input$yaxis))) +
      geom_point(aes(color = factor(chr))) +
      geom_text(data = subset(data, -log10(input$yaxis) > 2),
                aes(label = snp), hjust = 0, vjust = 0, size = 3, alpha = 0.5) +
      scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a",
                                    "#33a02c", "#fb9a99", "#e31a1c",
                                    "#fdbf6f", "#ff7f00", "#cab2d6",
                                    "#6a3d9a", "#ffff99", "#b15928"))
  })
  
  
  # Create the GWAS data table
  output$gwas <- renderTable({
    filtered_data()
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
