suppressMessages(library(rvest))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(RSQLite))
suppressMessages(library(readxl))
suppressMessages(library(readr))



scrap <- function(mainLink){
  
  # Get all the relevant links
  links <- generateListLink(mainLink)
  # Get that data associated with each links
  donneesCompletes <- getData(links)
  # Parse the raw data to get well structured table
  listTab <- parseData(donneesCompletes)
  # Write the data as csv files
  writeCsv(listTab)
  # Write the data in a database
  writeDatabase(listTab)
}

writeCsv <- function(list){
  #' This function write the list of dataframe into csv files
  #' @param list list of tables
  
  # Write all csv files
  write.csv(list[[1]],"../data/TablesCsvVih/RelationMutation.csv")
  write.csv(list[[2]],"../data/TablesCsvVih/RelationMutationPossible.csv")
  write.csv(list[[3]],"../data/TablesCsvVih/At_Least_N.csv")
  write.csv(list[[4]],"../data/TablesCsvVih/At_Least_N_Possible.csv")
  write.csv(list[[5]],"../data/TablesCsvVih/Treatment.csv")
}

writeDatabase <- function(list){
  #' Writes the table in the database
  #' @param list list of tables
  
  # Connect to database
  conn <- dbConnect(dbConnect(SQLite(), "../data/data/DatabaseHivResistance"))
  
  # Drop tables if they exist already
  dbSendStatement(conn,"DROP TABLE IF EXISTS RelationMutation")
  dbSendStatement(conn,"DROP TABLE IF EXISTS At_Least_N")
  dbSendStatement(conn,"DROP TABLE IF EXISTS RelationMutationPossible")
  dbSendStatement(conn,"DROP TABLE IF EXISTS At_Least_N_Possible")
  dbSendStatement(conn,"DROP TABLE IF EXISTS Treatment")
  
  # Create the table from the data
  dbWriteTable(conn, "RelationMutation", list[[1]])
  dbWriteTable(conn, "RelationMutationPossible", list[[3]])
  dbWriteTable(conn, "Treatment", list[[5]])
  
  # Disconnect to avoid conficts
  dbDisconnect(conn)
}

parseData <- function(donnees){
  #' This function creates all table in the form of a dataframe, refines the raw data and adds it into them
  #' @param donnees The raw data in a dataframe
  #' @return A list of tables with the refined data in it 
  
  # RelationMutation: table linking the relationID, the treatmentName, the mutation(s), if it require multiple ones and the treatment target
  RelationMutation <- data.frame(matrix(nrow = 0, ncol = 5)) # An empty dataframe
  x <- c("RelationID", "TreatmentName", "Mutation", "At_Least", "Treatment_Target")  # A vector containing column names
  colnames(RelationMutation) <- x
  RelationMutation[]<- lapply(RelationMutation, as.character)
  
  # RelationMutationPossible: table linking the relationID, the treatmentName, the possible mutation(s), if it require multiple ones and the treatment target
  RelationMutationPossible <- data.frame(matrix(nrow = 0, ncol = 5)) # An empty dataframe
  x <- c("PossibleRelationID", "TreatmentName", "Mutation", "At_Least", "Treatment_Target") # A vector containing column names
  colnames(RelationMutationPossible) <- x
  RelationMutationPossible[]<- lapply(RelationMutationPossible, as.character)
  
  # At_Least_N: table linking the relationID, the mutation, the minimum number of it
  At_Least_N <- data.frame(matrix(nrow = 0, ncol = 3)) # An empty dataframe
  x <- c("RelationID", "Mutation","Number") # A vector containing column names
  colnames(At_Least_N) <- x
  At_Least_N[]<- lapply(At_Least_N, as.character) 
  
  # At_Least_N_Possible: table linking the possiblerelationID, the possible mutation, the minimum number of it
  At_Least_N_Possible <- data.frame(matrix(nrow = 0, ncol = 3)) # An empty dataframe
  x <- c("PossibleRelationID", "Mutation","Number") # A vector containing column names
  colnames(At_Least_N_Possible)<- x
  At_Least_N_Possible[]<- lapply(At_Least_N_Possible, as.character) 
  
  # Treatment: table linking a treatment with what it target
  Treatment <- data.frame(matrix(nrow = 0, ncol = 2)) # An empty dataframe
  x <- c("TreatmentName", "TreatmentTarget") # A vector containing column names
  colnames(Treatment)<- x
  Treatment[]<- lapply(Treatment, as.character) 
  
  
  #This loop parse data in donnees to use other functions to format and structure them properly
  for(row in 1:nrow(donnees)){ 
    
    donnees[row,1] <- renameTreatment(donnees,row)  # Treatments are renamed
    donnees[row,4] <- urlTraduction(donnees, row)
    Treatment <- rbind(Treatment,
                       c(donnees[row,1],  # Treatments are renamed
                         donnees[row,4])) # The target of the treatment is renamed depending of its url
    
    # Create the list of tables
    listTab <- addMutations(donnees[row,], 
                            RelationMutation, 
                            At_Least_N, 
                            RelationMutationPossible, 
                            At_Least_N_Possible)
    
    # Each table is added to the main ones
    RelationMutation <- listTab[[1]]
    At_Least_N <- listTab[[2]]
    RelationMutationPossible <- listTab[[3]]
    At_Least_N_Possible <- listTab[[4]]
  }
  
  return(list(RelationMutation, At_Least_N, RelationMutationPossible, At_Least_N_Possible, Treatment))
}

generateListLink <- function(mainLink){
  #' This function gets the links containing the data tables
  #' @param mainLink The link containing all the other links we need
  #' @return The list of links containing data table
  
  #Read the main page
  mainPage <- rvest::read_html(mainLink)
  # Extract links from page
  links <- mainPage %>%
    rvest::html_nodes("a") %>%  # Select elements containing links
    rvest::html_attr("href") %>% # Get link content
    unique() # Removes every duplicate link
  # First filtering to remove url that aren't interesting to us
  links <- links[str_detect(links, "^https://hivfrenchresistance.org/hiv-french-resitance-")]
  links <- links[grepl(pattern = "protease|transcriptase|attachment|hiv-2|fusion|integrase|capsid",links,ignore.case = T)]  # Second filtering where unnecessary links are removed (hardcoded)
  return(links)
}

getData <- function(links){
  #' This function gets the data tables of multiple links
  #' @param links a list of links
  #' @return A dataframe containing the full data' 
  
  # Create a dataframe and define its column names
  donnees <- data.frame(matrix(nrow = 0, ncol = 4)) 
  x <- c("Treatment", "Mutation with associated resistance", "Mutation with possible resistance", "Treatment target") 
  colnames(donnees) <- x 
  donnees[]<- lapply(donnees, as.character) # Turns donnees' data into char to avoid conflicts with bind_rows()
  # for each link selected, tables from links are added to donnees
  for(link in links){
    donnees <- dplyr::bind_rows( # Add content of each table to donnees
      donnees, functionScrapping(link)) # Use functionScrapping on each link to get its table
  }
  return(donnees)
}

functionScrapping <- function(url){
  #' This function take tables from an url, put it in a dataframe, removes the first row and changes the column names to "Treatment", "Mutation with associated resistance", "Mutation with possible resistance"
  #' @param url The url of the scrapped page (urls from https://hivfrenchresistance.org/ for this script)
  #' @return a dataframe with the data we want
  table <- read_html(url) %>% html_elements("table")
  #Create the dataframe
  dataframe <- data.frame(matrix(nrow = 0, ncol = 4))
  # extract every line from the table
  lignes <- table %>% html_elements("tr")
  # remove the first line (contains no data)
  lignes <- lignes[-1]
  for (ligne in lignes) {
    # Extract each cell from each line
    cellules <- ligne %>% html_elements("td")
    # Add each line to a vector of line
    cleanedData <- cleanData(cellules)
    if(length(cleanedData)==2){
      cleanedData[length(cleanedData) +1] <- NA
    }
    cleanedData[length(cleanedData)+1] <- url
    dataframe <- rbind(dataframe, cleanedData)
  }
  #rename the dataframe to match each table
  colnames(dataframe) <- c("Treatment", "Mutation with associated resistance", "Mutation with possible resistance", "Treatment target") 
  return(dataframe) # Return the data frame
}

functionAutoIncrement <- function(tableName){
  #' This function takes a table with an id
  #' @param tableName the name of the table we want to auto increment
  #' @return This function returns the value of the next id
  return(nrow(tableName)+1)
}

addMutations <- function(line, tabMutat, n_mutat ,tabMutatPossible, n_possible){
  #' This function turn a line of raw data into tables linking mutations, treatments their target 
  #' and if there is multiple mutations linked to the resistance
  #' @param line the line of the raw data
  #' @param tabMutat The data table linking mutation, treatments and their target
  #' @param n_mutat The data table linking multiple mutations to a mutationID
  #' @param tabMutatPossible The data table linking possible mutation, treatments and their target
  #' @param n_possible The data table linking multiple possible mutations to a PossibleMutationID
  #' @return return a list of dataframe corresponding
  
  # Store the treatment and it's target
  treatment <- as.character(line[1])
  treatmentTarget <- as.character(line[4])
  if(is.na(line[2])) line[2] <- ""
  if(is.na(line[3])) line[3] <- ""
  # Split each different mutation and format them
  mutations <- str_split_1(as.character(line[2]), "\\\n")
  mutations <- lapply(mutations, mutationAmong)
  
  # Split each different possible mutation and format them
  mutationsPossibles <- str_split_1(as.character(line[3]), "\\\n")
  mutationsPossibles <- lapply(mutationsPossibles, mutationAmong)
  
  # Parse each possible mutation
  for (mut in mutations) {
    
    at_least <- substring(mut,1,1)# To check if the first character is a number
    # If it's one, it's stored in a variable and the mutation and split
    if(grepl("[0-9]", at_least)){ 
      number <- at_least
      at_least <- TRUE
      mutat_list <- str_split_1(mut,",")
      for(m in mutat_list){
        # Removes unnecessary characters
        m <- gsub("[0-9]\\(","",m)
        m <- gsub(")","",m, fixed = TRUE)
        
        # Add each mutation associated with the mutation ID and the minimum number of it
        vectorAtLeast <- c(as.character(functionAutoIncrement(tabMutat)), m, number)
        vectorAtLeast <- as.data.frame(t(vectorAtLeast))
        colnames(vectorAtLeast) <- colnames(n_mutat)
        n_mutat <- rbind(n_mutat,vectorAtLeast)
      }
    }
    else{
      at_least <- FALSE
    }
    # Create a vector with the data to add to the table
    finalVector <- c(as.character(functionAutoIncrement(tabMutat)), treatment, mut, at_least, treatmentTarget)
    finalVector <- as.data.frame(t(finalVector))
    colnames(finalVector) <- colnames(tabMutat)
    
    #Add the data to the table
    tabMutat<- rbind(tabMutat,finalVector)
  }
  
  # Parse each possible mutation
  for (mut in mutationsPossibles) {
    at_least <- substring(mut,1,1)# To check if the first character is a number
    # If it's one, it's stored in a variable and the mutation and split
    if(grepl("[0-9]", at_least)){
      number <- at_least
      at_least = TRUE
      mutat_list <- str_split_1(mut,",")
      for(m in mutat_list){
        # Removes unnecessary characters
        m <- gsub("[0-9]\\(","",m)
        m <- gsub(")","",m, fixed = TRUE)
        
        # Add each mutation associated with the mutation ID and the minimum number of it
        vectorAtLeast <- c(as.character(functionAutoIncrement(tabMutatPossible)), m,number)
        vectorAtLeast <- as.data.frame(t(vectorAtLeast))
        colnames(vectorAtLeast) <- colnames(n_possible)
        n_possible <- rbind(n_possible,vectorAtLeast)
      }
    }else{
      at_least <- FALSE
    }
    # Create a vector with the data to add to the table
    finalVector <- c(as.character(functionAutoIncrement(tabMutatPossible)), treatment, mut, at_least, treatmentTarget)
    finalVector <- as.data.frame(t(finalVector))
    colnames(finalVector) <- colnames(tabMutatPossible)
    
    #Add the data to the table
    tabMutatPossible<- rbind(tabMutatPossible,finalVector)
  }
  resList <- list(tabMutat, n_mutat ,tabMutatPossible, n_possible)
  return(resList)
}

addTreatment <- function(treatments, t){
  #' This function add a string containing treatments separated by / to a list of treatments
  #' @param treatments The list of treatments
  #' @param t The treatment we want to add
  #' @return a dataframe with Treatments added
  t <- strsplit(t,"/") # string separator
  t[]<- lapply(t, as.character) # t content is typed as character
  t <- data.frame(t)
  colnames(t) <- c("TreatmentName") # give it the same column names for biding
  colnames(treatments) <- c("TreatmentName") # give it the same column names for biding
  treatments <- rbind(treatments,t)
  return(treatments) # Return the data frame
}

cleanData <- function(data){
  #' This function clean the data of a row from remaining CSS or HTML
  #' @param line The row we want to clean
  #' 
  data <- html_text2(data)
  for(i in 1:length(data)){
    data[i] <- gsub(" ", "", data[i],perl = TRUE)   # delete remaining html without <br>
    data[i] <- gsub("@.*", "", data[i], perl = TRUE)  # delete remaining CSS 
    data[i] <- gsub("\\[(.*?)\\]", "", data[i]) ## delete unnecessary brackets
  }
  return(data)
}

renameTreatment <- function(donnees, line){
  #' This renames the abreviation of a treatment with it's real name
  #' @param donnees The dataframe we want to modify
  #' @param line The row whose treatments we want to rename
  donnees[line,1] <- gsub("3TC", "Lamivudine", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("FTC", "Emtricitabine", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("ddl", "Didanosine", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("d4T", "Stavudine", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("ABC", "Abacavir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("TDF", "Tenofovir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("TAF", "Tenofovir alafenamide", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("ZDV", "Zidovudine", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("ISL", "Islatravir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("LPV", "Lopinavir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("ATV", "Atazanavir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("DRV", "Darunavir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("RTV", "Ritonavir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("EFV", "Efavirenz", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("NVP", "Nevirapine", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("ETR", "Etravirine", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("RPV", "Rilpivirine", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("DOR", "Doravirine", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("ENFT20", "Enfurvitide", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("FTR", "Fostemsavir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("RAL", "Raltegravir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("EVG", "Elvitegravir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("DTG", "Dolutegravir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("CAB", "Cabotegravir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("BIC", "Bictegravir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("RAL", "Raltegravir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("LEN", "Lenacapavir", donnees[line,1],fixed = TRUE)
  donnees[line,1] <- gsub("\\*", "", donnees[line,1])
  donnees[line,1] <- gsub("<.*>", "", donnees[line,1])
  
  return(donnees[line,1])
}

urlTraduction <- function(donnees, line){
  #' This renames the urls into their corresponding effects
  #' @param donnees The dataframe we want to modify
  #' @param line The row whose url we want to rename
  donnees[line,4] <- gsub("https://hivfrenchresistance.org/hiv-french-resitance-nucleoside-and-nucleotide-reverse-transcriptase-inhibitors/", "Nucleoside and nucleotide reverse transcriptase inhibitors", donnees[line,4],fixed = TRUE)
  donnees[line,4] <- gsub("https://hivfrenchresistance.org/hiv-french-resitance-protease-inhibitors/", "Protease inhibitors", donnees[line,4],fixed = TRUE)
  donnees[line,4] <- gsub("https://hivfrenchresistance.org/hiv-french-resitance-non-nucleoside-reverse-transcriptase-inhibitors/", "Non nucleoside reverse transcriptase inhibitors", donnees[line,4],fixed = TRUE)
  donnees[line,4] <- gsub("https://hivfrenchresistance.org/hiv-french-resitance-fusion-inhibitor/", "Fusion inhibitor", donnees[line,4],fixed = TRUE)
  donnees[line,4] <- gsub("https://hivfrenchresistance.org/hiv-french-resitance-attachment-inhibitor/", "Resistance attachment inhibitor", donnees[line,4],fixed = TRUE)
  donnees[line,4] <- gsub("https://hivfrenchresistance.org/hiv-french-resitance-integrase-strand-transfer-inhibitors/", "Integrase strand transfer inhibitor", donnees[line,4],fixed = TRUE)
  donnees[line,4] <- gsub("https://hivfrenchresistance.org/hiv-french-resitance-capsid-inhibitors/", "Capsid inhibitors", donnees[line,4],fixed = TRUE)
  donnees[line,4] <- gsub("https://hivfrenchresistance.org/hiv-french-resitance-hiv-2/", "Hiv-2", donnees[line,4],fixed = TRUE)
  
}

mutationAmong <- function(mutations){
  #' This functions removes "+" and "N mutations among" to format it all in a similar way (N(mutation1, mutation2...))
  #' @param mutations
  #' @return mutations in a more formal notation
  #' 
  temp <- unlist(str_split(mutations,"\\+")) # Split the mutation if it contains "+"
  temp <- unlist(lapply(temp, str_trim)) # removes unecessary spaces
  # This block takes mut1 + mut2 + mut3 +...+ mutN and makes it into N(mut1, mut2, mut3, ..., mutN)
  if(length(temp) > 1){   # check if it was split
    # this encloses the mutations between parenthesis
    temp[1] <- paste("(", temp[1],sep = "") 
    temp[1] <- paste(length(temp), temp[1],sep = "")
    temp[length(temp)] <- paste(temp[length(temp)], ")",sep = "") 
    mutations <- paste(temp, collapse = ",") # this add all split parts and put a "," between them
  }
  
  # This block takes "At least N mutations among: mut1, mut2 +...+ mutM and makes it into N(mut1, mut2, mut3, ..., mutM)
  temp <- str_split_1(mutations,"(mutation.?among.?:|mutation.?|among.?:)")
  temp <- temp[temp != ""] # Removes empty strings
  if(length(temp) > 1){ # Check if it was split
    # Removes everything not useful that remains
    temp <- lapply(temp, function(x) gsub("Atleast","", x, ignore.case = T))
    temp <- lapply(temp, function(x) gsub(":","", x, ignore.case = T))
    temp <- lapply(temp, function(x) gsub("or",",", x))
    temp <- lapply(temp, function(x) gsub("among","", x, ignore.case = T))
    temp <- lapply(temp, function(x) gsub("mutation.?","", x, ignore.case = T))
    temp <- unlist(temp)
    
    # replace some written numbers by their numeric counterpart
    temp[1] <- gsub("one", "1", temp[1])
    temp[1] <- gsub("two", "2", temp[1])
    temp[1] <- gsub("three", "3", temp[1])
    temp[1] <- gsub("four", "4", temp[1])
    temp[1] <- str_trim(temp[1])
    temp[2] <- str_trim(temp[2])
    #Enclose content between parethesis
    temp[1] <- paste(temp[1], "(",sep = "") 
    temp[2] <- gsub(pattern = " ", replacement = "", temp[2])
    temp[length(temp)] <- paste(temp[length(temp)], ")",sep = "")
    mutations <- paste(temp, collapse = ",")# Put back each part together with "," inbetween
  }
  mutations <- str_trim(gsub("\\(,", "\\(", mutations))
  return(mutations) 
}