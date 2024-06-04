library(shiny)
library(shinyalert)
library(shinydisconnect)
library(colourpicker)
library(bslib)
library(DT)
library(RSQLite)

source("./functionsScrapping.R")
source("./mutationsFunctions.R")

# The sidebar
sidebarSelect <- sidebar(
  "Select all necessary parameters before analyzing",
  
  # A card containing the element to input a file
  card(
    card_header("Select a fasta file"),
    fileInput(inputId = "file", label = "Select a file", accept = ".fna")
  ),
  
  # The card to decide whether or not you want to use the web scrapping
  card(
    card_header("Scrap website data"),
    radioButtons(
      "scrap",
      "Select option",
      choices = list("If not done recently" = 1, "Yes" = 2, "No" = 3),
      selected = 1
    )
  ),
  
  # The type of hiv (hiv1 or hiv2) you want to use
  card(
    card_header("Hiv type"),
    selectInput(
      "HivType",
      "Select option",
      choices = list("Hiv-1" = 1, "Hiv-2" = 2),
      selected = 1
    )
  ),
  
  # The start button
  card(
    card_header("If everything has been selected, Start analyzing"),
    actionButton("analyze", "Analyze")
  ),
)

# The main page
ui <- page_sidebar(
  
  title = "Hiv mutation sequence analyser",
  sidebar = sidebarSelect,
  
  # The two tables containing the data
  card(
    card_header("Mutations associated with resistances"),
    dataTableOutput("MutationsFound"),
  ),
  card(
    card_header("Mutations possibly associated with resistances"),
    dataTableOutput("PossibleMutationsFound")
  ),
)


server <- function(input, output) {
  
  # This code starts when the analyse button is pressed
  observeEvent(input$analyze,{
    
    # If no file is selected, an alert is shown
    if(is.null(input$file$datapath)){
      shinyalert(
        title = "Please select a fna file",
        text = "This file is necessary for analysis",
        size = "s", 
        closeOnEsc = TRUE,
        closeOnClickOutside = FALSE,
        html = FALSE,
        type = "warning",
        showConfirmButton = TRUE,
        showCancelButton = FALSE,
        confirmButtonText = "OK",
        confirmButtonCol = "#AEDEF4",
        timer = 0,
        imageUrl = "",
        animation = TRUE
      )
    }else{
    # Message to make the user know that the app did not just crash
    shinyalert(
      title = "Please wait for the sequences analysis",
      text = "This may take a few minutes",
      size = "s", 
      closeOnEsc = TRUE,
      closeOnClickOutside = FALSE,
      html = FALSE,
      type = "info",
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "OK",
      confirmButtonCol = "#AEDEF4",
      timer = 0,
      imageUrl = "",
      animation = TRUE
    )
    
      
    # To scrap if not done recently
    if(input$scrap == 1){
      if(!is.na(file.info("./scriptsR/DatabaseHivResistance")$mtime) && difftime(Sys.time(),file.info("./scriptsR/DatabaseHivResistance")$mtime, units = "days") > 30){
        scrap("https://hivfrenchresistance.org/hiv-french-resitance-tables-of-rules/")
      }
    }
    # To force the scrapping (even if not done recently)
    else if(input$scrap == 2){
      scrap("https://hivfrenchresistance.org/hiv-french-resitance-tables-of-rules/")
    }
      
    # Connect to the database and grab the data
    conn <- dbConnect(SQLite(), "../data/DatabaseHivResistance")
    data <- lapply(setNames(nm = dbListTables(conn)), dbReadTable, conn= conn)
    dbDisconnect(conn)
    
    # Get the mutations and possible mutations
    datapath <- input$file$datapath
    if(input$HivType == 1){
      
      # Analyse and get the data to generate a report
      mutationsFound <-  suppressMessages(analyseSequences(datapath, data[[1]], "../data/genomeHIV.fna"))
      possibleMutationFound <- suppressMessages(analyseSequences(datapath, data[[2]], "../data/genomeHIV.fna"))
      suppressMessages(generateReport(mutationsFound, possibleMutationFound, outputFile = "../Report_Mutations.pdf"))
    }else{
      # Change the reference sequence's path
      referencePath <- "../data/genomeHIV2.fna"
      
      # Select the Hiv-2 data
      data[[1]] <- data[[1]][(data[[1]][,5] %in% "Hiv-2"),]
      data[[1]] <- data[[1]][!(data[[1]][2] == ""),]
      data[[2]] <- data[[2]][(data[[2]][,5] %in% "Hiv-2"),]
      data[[2]] <- data[[2]][!(data[[2]][2] == ""),]
      data[[3]] <- data[[3]][!(data[[3]][,2] %in% "Hiv-2"),]
      
      # Associate Hiv-2 data to the target treatment
      for(i in 1:nrow(data[[1]])){
        data[[1]][i,5] <- data[[3]][grepl(data[[1]][i,2], data[[3]][,1]),][1,2]
      }
      for(i in 1:nrow(data[[2]])){
        data[[2]][i,5] <- data[[3]][grepl(data[[2]][i,2], data[[3]][,1]),][1,2]
      }
      
      # Analyse and get the data to generate a report
      mutationsFound <-  suppressMessages(analyseSequencesHiv2(datapath, data[[1]], "../data/genomeHIV2.fna"))
      possibleMutationFound <- suppressMessages(analyseSequencesHiv2(datapath, data[[2]], "../data/genomeHIV2.fna"))
      suppressMessages(generateReport(mutationsFound, possibleMutationFound, outputFile = "../Report_Mutations_Hiv2.pdf"))
    }
    
    # Add the tables to their respective cards
    output$MutationsFound <- renderDataTable(DT::datatable(mutationsFound))
    output$PossibleMutationsFound <- renderDataTable(DT::datatable(possibleMutationFound))
    }
  })
  
}

shinyApp(ui = ui, server = server )