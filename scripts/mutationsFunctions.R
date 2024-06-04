suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(Biostrings))
suppressMessages(library(dplyr))
suppressMessages(library(RSQLite))
suppressMessages(library(seqinr))
suppressMessages(library(readxl))
suppressMessages(library(read.gb))
suppressMessages(library(readr))
suppressMessages(library(rmarkdown))
suppressMessages(library(utils))

# Global variables (to avoid leading them multiple times)
geneGagPol <- read.gb("../data/HIVgagpol.gbk")
geneEnv <- read.gb("../data/envGene.gbk")


main <- function(content){
  # Connect to the database and grab the data
  conn <- dbConnect(dbConnect(SQLite(), "../data/DatabaseHivResistance"))
  data <- lapply(setNames(nm = dbListTables(conn)), dbReadTable, conn= conn)
  dbDisconnect(conn)
  
  if(content[[2]]) print("Data successfully collected")
  
  filePath <- content[[1]]
  referencePath <- "../data/genomeHIV.fna"
  
  if(content[[3]]){
    
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
    
    # Get the mutations and possible mutations
    mutationsFound <-  suppressMessages(analyseSequencesHiv2(filePath, data[[1]],referencePath))
    possibleMutationFound <-  suppressMessages(analyseSequencesHiv2(filePath, data[[2]],referencePath))
    if(content[[2]]) print("Data analyzed")
    
    # Generate the report
    suppressMessages(generateReport(mutationsFound, possibleMutationFound, outputFile = "../Report_Mutations_Hiv2.pdf"))
    if(content[[2]]) print("Report generated")
  }
  else{
    # Get the mutations and possible mutations
    mutationsFound <-  suppressMessages(analyseSequences(filePath, data[[1]],referencePath))
    possibleMutationFound <-  suppressMessages(analyseSequences(filePath, data[[2]],referencePath))
    if(content[[2]]) print("Data analyzed")
    
    # Generate the report
    suppressMessages(generateReport(mutationsFound, possibleMutationFound))
    if(content[[2]]) print("Report generated")
  }
}

generateReport <- function(possibleResistances, resistances, outputFile = "../Report_Mutations.pdf"){
  #' This function generate a pdf report from the data found
  #' @param possibleResistances A dataframe containing all possible resistances found
  #' @param resistances A dataframe containing all resistances found
  #' @param outputFile The path of the output pdf file

  render("./reportMaker.Rmd", params = list(allPossibleResistancesFound =possibleResistances, allResistancesFound=resistances), output_file = outputFile)
}

analyseSequencesHiv2 <- function(fastaFile, donneesMutation, fastaRef){
  #' this function analyses hiv-2 sequences so find mutations associated with resistances
  #' @param fastaFile A fasta file containing the sequence whose mutations we want to find
  #' @param donneesMutation A dataframe containing all mutations
  #' @param fastaRef A file containing the fasta reference sequence
  #' @return A dataframe containing all the associated treatment resistance along with the mutations
  
  # Create a nested list of all readable mutations
  readableMutations <- list()
  mutations <- donneesMutation
  for(i in 1:nrow(mutations)){
    mutationList <- getReadableMutation(mutations[i,3])
    readableMutations[[i]] <- list(mutationList, mutations[i,2], mutations[i,5])
  }
  
  # The hiv sequences whose mutation we want to detect
  seqList <- read.fasta(fastaFile, as.string = TRUE)
  
  # Initialization of the list of resistances we found
  resistancesFound <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(resistancesFound) <- c("SequenceNumber", "Mutation", "Target", "Treatment")
  
  # Get the reference sequence
  hivReferenceSequence <- read.fasta(fastaRef, as.string = TRUE, seqonly = TRUE)[[1]]
  
  # Progress bar
  pb = txtProgressBar(min = 0, max = 100, initial = 0,title = "Sequence analysis in progress", width = 50)
  
  # Parse all the sequences 
  for(i in 1:length(seqList)){
    
    # Get the alignment, mutation and the list that are not relevant to the sequence
    mutationsProtein <- getProteinMutationsHiv2(seqList[[i]][1], hivReferenceSequence)
    targetReadableMutations <- removeUnecessaryLists(readableMutations, mutationsProtein[[2]])
    # if the sequence contains possible relevant mutations
    if(!isEmpty(targetReadableMutations)){
      
      # Parse the treatments associated to the mutation and look for the mutation list
      for(j in 1:length(targetReadableMutations)){
        depthSearch <- depthSearchMutation(
          targetReadableMutations[[j]][[1]], 
          mutationsProtein[[1]])
        
        # Si des mutation correspondantes sont trouvées, elles sont ajoutées
        if(!is.null(depthSearch) && !isEmpty(depthSearch)){
          depthSearch <- paste0(lapply(depthSearch, unlist, use.names=FALSE),collapse = ",")
          resistancesFound[nrow(resistancesFound)+1, ] <- list(i, depthSearch,
                                                               targetReadableMutations[[j]][[2]],
                                                               targetReadableMutations[[j]][[3]])
        }
      }
    }
    setTxtProgressBar(pb, 100*i/length(seqList))
  }
  close(pb)
  return(resistancesFound)
}

analyseSequences <- function(fastaFile, donneesMutation, fastaRef){
  #' this function analyses hiv sequences so find mutations associated with resistances
  #' @param fastaFile A fasta file containing the sequence whose mutations we want to find
  #' @param donneesMutation A dataframe containing all mutations
  #' @param fastaRef A file containing the fasta reference sequence
  #' @return A dataframe containing all the associated treatment resistance along with the mutations
  
  # Create a nested list of all readable mutations
  readableMutations <- list()
  mutations <- donneesMutation
  for(i in 1:nrow(mutations)){
    mutationList <- getReadableMutation(mutations[i,3])
    readableMutations[[i]] <- list(mutationList, mutations[i,2], mutations[i,5])
  }
  
  # The hiv sequences whose mutation we want to detect
  seqList <- read.fasta(fastaFile, as.string = TRUE)
  
  # Initialization of the list of resistances we found
  resistancesFound <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(resistancesFound) <- c("SequenceNumber", "Mutation", "Target", "Treatment")
  
  # Get the reference sequence
  hivReferenceSequence <- read.fasta(fastaRef, as.string = TRUE, seqonly = TRUE)[[1]]
  
  # Progress bar
  pb = txtProgressBar(min = 0, max = 100, initial = 0,title = "Sequence analysis in progress", width = 50)
  
  # Parse all the sequences 
  for(i in 1:length(seqList)){
    
    # Get the alignment, mutation and the list that are not relevant to the sequence
    mutationsProtein <- getProteinMutations(seqList[[i]][1], hivReferenceSequence)
    targetReadableMutations <- removeUnecessaryLists(readableMutations, mutationsProtein[[2]])
    
    # if the sequence contains possible relevant mutations
    if(!isEmpty(targetReadableMutations)){
      
      # Parse the treatments associated to the mutation and look for the mutation list
      for(j in 1:length(targetReadableMutations)){
        depthSearch <- depthSearchMutation(
          targetReadableMutations[[j]][[1]], 
          mutationsProtein[[1]])
        
        # Si des mutation correspondantes sont trouvées, elles sont ajoutées
        if(!is.null(depthSearch) && !isEmpty(depthSearch)){
          depthSearch <- paste0(lapply(depthSearch, unlist, use.names=FALSE),collapse = ",")
          resistancesFound[nrow(resistancesFound)+1, ] <- list(i, depthSearch,
                                                               targetReadableMutations[[j]][[2]],
                                                               targetReadableMutations[[j]][[3]])
        }
      }
    }
    setTxtProgressBar(pb, 100*i/length(seqList))
  }
  close(pb)
  return(resistancesFound)
}

getProteinMutationsHiv2 <- function(sequence, referenceSequence){
  #' This function takes a sequence and it's reference sequence, translate them into protein and returns the mutations with their locations
  #' @param sequence The sequence in input
  #' @param referenceSequence The reference sequence that will be used to find mutations
  #' @return a list with the mutations, their location and on which protein they are
  
  # Align the sequence on its reference
  alignmentSequence <- getAlignmentDna(toupper(sequence), toupper(referenceSequence))
  insertionAt <- insertion(alignmentSequence)
  deletionAt <- deletion(alignmentSequence)
  inDels <- c()
  for(i in 1:length(insertionAt)){
    if(!is.na(start(insertionAt)[[1]][i])) inDels <- c(inDels, paste0(c("Insertion at codon", start(insertionAt)[[1]][i]%/%3, "-", end(insertionAt)[[1]][i]%/%3), collapse = ""))
  }
  
  # Extract aligned sequence and reference
  alignedSequence <- as.character(aligned(pattern(alignmentSequence)))
  alignedReference <- as.character(aligned(subject(alignmentSequence)))
  
  # Get sequence location on the dna to get the sequence location of the protein
  sequenceLocation <- start(subject(alignmentSequence))
  
  # Get the protein names associated to the sequence and their associated shortened sequences
  sequenceWithProteinLocation <- positionToProteinsHiv2(sequenceLocation, sequence)
  
  mutationsProtein <- list(list(),list())
  # This loop get every elements from this list and includes it in the associated protein sequence
  for(elem in sequenceWithProteinLocation){
    # Get the dna from the reference associated with the protein
    referenceProtein <- substr(referenceSequence, elem[5], elem[6])
    
    # Create a version of the reference alignment sequence without indels
    alignedReference2 <- gsub("-", "", alignedReference, fixed = T)
    
    # Use this version to replace the protein w indels
    referenceProtein <- gsub(pattern = alignedReference2, replacement = alignedReference, x = referenceProtein)
    
    # Add replace the sequence in the reference genome
    referenceWithSequence <- gsub(pattern = alignedReference,replacement = alignedSequence, x = referenceProtein)
    # First the sequences are translated into proteins
    if(!is.na(referenceWithSequence)&&!is.na(referenceProtein)){
      proteinSequence <- paste0(translateDna(referenceWithSequence), collapse = "")
      proteinReference <- paste0(translateDna(referenceProtein), collapse = "")
    
      # Get the locations of mutations
      mutationsProtein[[1]] <- append(mutationsProtein[[1]],geneComparison(proteinReference, proteinSequence))
      mutationsProtein[[2]] <- append(mutationsProtein[[2]], elem[1])
    }
  }
  mutationsProtein[[1]] <- append(mutationsProtein[[1]], inDels)
  return(mutationsProtein)
}

getProteinMutations <- function(sequence, referenceSequence){
  #' This function takes a sequence and it's reference sequence, translate them into protein and returns the mutations with their locations
  #' @param sequence The sequence in input
  #' @param referenceSequence The reference sequence that will be used to find mutations
  #' @return a list with the mutations, their location and on which protein they are
  
  # Align the sequence on its reference
  alignmentSequence <- getAlignmentDna(toupper(sequence), toupper(referenceSequence))
  insertionAt <- insertion(alignmentSequence)
  deletionAt <- deletion(alignmentSequence)
  inDels <- c()
  for(i in 1:length(insertionAt)){
    if(!is.na(start(insertionAt)[[1]][i])) inDels <- c(inDels, paste0(c("Insertion at codon", start(insertionAt)[[1]][i]%/%3, "-", end(insertionAt)[[1]][i]%/%3), collapse = ""))
  }

  # Extract aligned sequence and reference
  alignedSequence <- as.character(aligned(pattern(alignmentSequence)))
  alignedReference <- as.character(aligned(subject(alignmentSequence)))
  
  # Get sequence location on the dna to get the sequence location of the protein
  sequenceLocation <- start(subject(alignmentSequence))
  
  # Get the protein names associated to the sequence and their associated shortened sequences
  sequenceWithProteinLocation <- positionToProteins(sequenceLocation, sequence)
  
  mutationsProtein <- list(list(),list())
  # This loop get every elements from this list and includes it in the associated protein sequence
  for(elem in sequenceWithProteinLocation){

    # Get the dna from the reference associated with the protein
    referenceProtein <- substr(referenceSequence, elem[5], elem[6])
    
    # Create a version of the reference alignment sequence without indels
    alignedReference2 <- gsub("-", "", alignedReference, fixed = T)
    
    # Use this version to replace the protein w indels
    referenceProtein <- gsub(pattern = alignedReference2, replacement = alignedReference, x = referenceProtein)
    
    # Add replace the sequence in the reference genome
    referenceWithSequence <- gsub(pattern = alignedReference,replacement = alignedSequence, x = referenceProtein)
    # First the sequences are translated into proteins
    proteinSequence <- paste0(translateDna(referenceWithSequence), collapse = "")
    proteinReference <- paste0(translateDna(referenceProtein), collapse = "")
    
    # Get the locations of mutations
    mutationsProtein[[1]] <- append(mutationsProtein[[1]],geneComparison(proteinReference, proteinSequence))
    mutationsProtein[[2]] <- append(mutationsProtein[[2]], elem[1])
  }
  mutationsProtein[[1]] <- append(mutationsProtein[[1]], inDels)
  return(mutationsProtein)
}

geneComparison <- function(referenceGene, seqGene){
  #' This function compares the reference gene to the sequenced one, turns it into a protein
  #' detects its mutation and compare each one to the ones in the database
  #' @param referenceGene The reference gene from nbci database
  #' @param seqGene The sequence gene that we want to analyze
  #' @return A list of mutations
  
  ref <- unlist(strsplit(referenceGene, ""))
  seq <- unlist(strsplit(seqGene, ""))
  iMut <- (1:length(ref))[ref != seq]
  return(paste0(ref[iMut], iMut, seq[iMut]))
}

removeUnecessaryLists <- function(listSearched, targets){
  #' Create a list without unnecessary sublists (which does not contain the right targets)
  #' @param listSearched The list of mutations associated with resistance
  #' @param target The target we want to look for
  #' @return A list containing only the relevant sublists (the ones with the right target)
  
  # The result list we want to get
  res <- list()
  # Check if each list contains the right target and add it to the result list if it is
  for(i in 1:length(listSearched)){
    for(t in targets){
      switch (t,   # Check for all cases to make it match with the database
              "p66" = t <-"transcriptase",
              "gp120" = t <-"attachment",
              "gp41" = t <-"fusion",
              "integrase" = t <- "strand")
      
      # If it matches it, add it
      if(grepl(t, listSearched[[i]][[3]], ignore.case = TRUE)){
        res[[length(res)+1]] <- listSearched[[i]]
      }
    }
  }
  return(res)
}

searchMutations <- function(element, mutationList){
  #' Look for a mutation for a reference compared to mesured mutations
  #' @param element The reference mutation
  #' @param mutationList The mutations whose presence we want to check
  #' @return A matching mutation or NULL if there is none
  
  # We compare the mutation to each detected one (until one is found)
  for(i in 1:length(mutationList)){
    
    # The location of the mutations are stored for later
    mutationLocation <- parse_number(mutationList[[i]])
    splitMut <- str_split_1(mutationList[[i]], "[0-9]")
    mutatedProtein <- splitMut[nzchar(splitMut)][2]
    referenceLocation <- as.numeric(element[4])
    
    
    # In case of an insertion
    if((length(element) == 1) && grepl("Insertion", mutationList[[i]])&& grepl("Insertion", element)){
      # Get the start and end of the detected insertion
      splitMutation <- str_split_1(mutationList[[i]],"-")
      insertionStartEnd <- c(parse_number(splitMutation[1]), parse_number(splitMutation[2]), parse_number(element))
      # If it matches returns it
      if((insertionStartEnd[1] <= insertionStartEnd[3]) && (insertionStartEnd[3] <= insertionStartEnd[2])){
        return(mutationList[[i]])
      }
    }
    # If there is one, add it to the result list and returns it
    else if((mutationLocation == referenceLocation) && grepl(mutatedProtein, element[3], fixed = TRUE)){
      return(mutationList[[i]])
    }
  }
  return(NULL)
}

depthSearchMutation <- function(listSearched, mutationList){
  #' Look for the presence of mutation in a vector matching the ones in a nested list related to a target
  #' @param listSearched The nested list containing mutations and the target linked to them
  #' @param mutationList The mutations whose presence we want to check
  #' @return A vector containing all the mutations with their target and the treatment related to them
  
  # If there's an empty list or if the list is null return null
  if(is.null(listSearched)||isEmpty(listSearched)) return(NULL)
  
  # Value initialized at -1 for mutation like At least N mutations among:
  at_Least <- -1
  
  # The return value is initialized
  result <- list()
  
  # If a character vector is directly inserted check the mutation presence in list search and add the matching ones to result
  if(class(listSearched) == "character") result <-append(result, searchMutations(listSearched, mutationList)) 
  else { 
    # Else parse all elements in the list
    for(element in listSearched){
      
      # If it's a sublist
      if(class(element) == "list"){
        # In case of "at least N among"
        if(at_Least > 0){
          # Initialise the sublist
          AtLeastList <- list()
          AtLeastList<- depthSearchMutation(element, mutationList)
          # If there is enough mutations add it to the results
          if(!length(AtLeastList)< at_Least){
            result <- append(result, AtLeastList)
          }
          # Block the if as long as a number isn't found
          at_Least <- -1
        }
        else{
          result <-append(result, depthSearchMutation(element, mutationList))
        }
      }
      else  # If it's a string/character vector    
        if(class(element) == "character"){
          # If the character is a single number (for 2(A111B, C222D...) cases) 
          if((length(element) == 1)&&grepl("^[0-9]+$", element)){
            at_Least <- parse_number(element)
          }else{
            # Else look for matching content and adds it to the result
            result <-append(result, searchMutations(element, mutationList)) 
          }
        }
    }
  }
  if(isEmpty(result)) return(NULL)
  return(result)
}

positionToProteinsHiv2 <- function(position, sequence){
  #' Take a sequence and it's position and returns proteins associated with it
  #' @param position the location of the sequence of the HIV reference genome
  #' @param sequence The sequence we want to translate
  #' @return A vector that contains the name of the protein with its position in the genome
  
  geneList <- list()
  
  # Look if the sequence is in gag-pol
  if(1103 <= position+length(sequence) &&position <= 5754){
    # If the sequence aligns on gag-pol, extract it 
    geneList[[length(geneList)+1]] <- c("gag-pol", substr(sequence, max(1103-position,1), 4642-position),position-1102)
  }
  if(isEmpty(geneList)) return(c("any", sequence, position, nchar(sequence), 1, 10359))
  
  # The result list is initialized
  resList <- list()
  
  # parse each gene found
  for(elem in geneList){
    
    # If the gene is gag-pol
    if(elem[1] == "gag-pol"){
      
      # Extract the target protein that are interesting from a gag pol gbk file
      peptides <- extract.gb(geneGagPol, "mat_peptide")$NC_001802
      protease <- as.numeric(str_split_1(peptides[[3]]$mat_peptide[1,2], "\\.\\."))
      p66 <- as.numeric(str_split_1(peptides[[4]]$mat_peptide[1,2], "\\.\\."))
      integrase <- as.numeric(str_split_1(peptides[[6]]$mat_peptide[1,2], "\\.\\."))
      capsid <- as.numeric(str_split_1(peptides[[8]]$mat_peptide[1,2], "\\.\\."))
      
      # Set the relative location of the sequence on the gene
      location <- as.numeric(elem[3])
      
      # If the sequence covers the protease protein
      if(1536 <= location+length(elem[2])&& location <= 1830){
        
        # Get the start and end of the protein on the whole reference genome
        proteinStart <- protease[1]+1102
        proteinEnd <- protease[2]+1102
        
        # Extract the part of the sequence associated with that protein and add it to the result
        reducedSequence <- substr(elem[2],max(protease[1]-location,1) , protease[2] - location)
        # The protein name, sequence location, size, and the start and end of the protein is also returned
        resList[[length(resList)+1]] <- c("protease", reducedSequence, location, nchar(reducedSequence), proteinStart, proteinEnd)
      }
      
      # If the sequence covers the p66 (or p55) protein
      if(1833 <= location+length(elem[2]) && location <= 3507){
        
        # Get the start and end of the protein on the whole reference genome
        proteinStart <- p66[1]+1102
        proteinEnd <- p66[2]+1102
        
        # Extract the part of the sequence associated with that protein and add it to the result
        reducedSequence <- substr(elem[2],max(p66[1]-location,1) , p66[2] - location)
        # The protein name, sequence location, size, and the start and end of the protein is also returned
        resList[[length(resList)+1]] <- c("p66", reducedSequence, location, nchar(reducedSequence), proteinStart, proteinEnd)
      }
      
      # If the sequence covers the integrase protein
      if(3507  <= location+length(elem[2]) && location <= 4395){
        
        # Get the start and end of the protein on the whole reference genome
        proteinStart <- integrase[1]+1102
        proteinEnd <- integrase[2]+1102
        
        # Extract the part of the sequence associated with that protein and add it to the result
        reducedSequence <- substr(elem[2],max(integrase[1]-location,1) , integrase[2] - location)
        # The protein name, sequence location, size, and the start and end of the protein is also returned
        resList[[length(resList)+1]] <- c("integrase", reducedSequence, location, nchar(reducedSequence), proteinStart, proteinEnd)
      }
      
      # If the sequence covers the integrase protein
      if(390 <= location+length(elem[2]) && location <= 1077){
        
        # Get the start and end of the protein on the whole reference genome
        proteinStart <- capsid[1]+1102
        proteinEnd <- capsid[2]+1102
        
        # Extract the part of the sequence associated with that protein and add it to the result
        reducedSequence <- substr(elem[2],max(capsid[1]-location,1) , capsid[2] - location)
        # The protein name, sequence location, size, and the start and end of the protein is also returned
        resList[[length(resList)+1]] <- c("capsid", reducedSequence, location, nchar(reducedSequence), proteinStart, proteinEnd)
      }
    }
  }
  # If no result has been found
  if(isEmpty(resList)) return(c("any", sequence, position, length(sequence), 1, 10359))
  return(resList)
}

positionToProteins <- function(position, sequence){
  #' Take a sequence and it's position and returns proteins associated with it
  #' @param position the location of the sequence of the HIV reference genome
  #' @param sequence The sequence we want to translate
  #' @return A vector that contains the name of the protein with its position in the genome
  
  #Create the list of substrings depending of the gene
  geneList <- list()
  
  # Look if the sequence is in gag-pol
  if(336 <= position+length(sequence) &&position <= 4642){
    # If the sequence aligns on gag-pol, extract it 
    geneList[[length(geneList)+1]] <- c("gag-pol", substr(sequence, max(336-position,1), 4642-position),position-335)
  }
  # Look if the sequence covers env
  if(5771 <= position+length(sequence) && position <= 8341){
    geneList[[length(geneList)+1]] <- c("env", substr(sequence, max(5771-position,1), 8341-position), position-5770)
  }
  # If it doesn't cover either, return any protein
  if(isEmpty(geneList)) return(c("any", sequence, position, nchar(sequence), 1, 9181))
  
  # The result list is initialized
  resList <- list()
  
  # parse each gene found
  for(elem in geneList){
    
    # If the gene is gag-pol
    if(elem[1] == "gag-pol"){
      
      # Extract the target protein that are interesting from a gag pol gbk file
      peptides <- extract.gb(geneGagPol, "mat_peptide")$NC_001802
      protease <- as.numeric(str_split_1(peptides[[3]]$mat_peptide[1,2], "\\.\\."))
      p66 <- as.numeric(str_split_1(peptides[[4]]$mat_peptide[1,2], "\\.\\."))
      integrase <- as.numeric(str_split_1(peptides[[6]]$mat_peptide[1,2], "\\.\\."))
      capsid <- as.numeric(str_split_1(peptides[[8]]$mat_peptide[1,2], "\\.\\."))
      
      # Set the relative location of the sequence on the gene
      location <- as.numeric(elem[3])
      
      # If the sequence covers the protease protein
      if(protease[1] <= location+length(elem[2])&& location <= protease[2]){
        
        # Get the start and end of the protein on the whole reference genome
        proteinStart <- protease[1]+335
        proteinEnd <- protease[2]+335
        
        # Extract the part of the sequence associated with that protein and add it to the result
        reducedSequence <- substr(elem[2],max(protease[1]-location,1) , protease[2] - location)
        # The protein name, sequence location, size, and the start and end of the protein is also returned
        resList[[length(resList)+1]] <- c("protease", reducedSequence, location, nchar(reducedSequence), proteinStart, proteinEnd)
      }
      
      # If the sequence covers the p66 (or p55) protein
      if(p66[1] <= location+length(elem[2]) && location <= p66[2]){
        
        # Get the start and end of the protein on the whole reference genome
        proteinStart <- p66[1]+335
        proteinEnd <- p66[2]+335
        
        # Extract the part of the sequence associated with that protein and add it to the result
        reducedSequence <- substr(elem[2],max(p66[1]-location,1) , p66[2] - location)
        # The protein name, sequence location, size, and the start and end of the protein is also returned
        resList[[length(resList)+1]] <- c("p66", reducedSequence, location, nchar(reducedSequence), proteinStart, proteinEnd)
      }
      
      # If the sequence covers the integrase protein
      if(integrase[1]  <= location+length(elem[2]) && location <= integrase[2]){
        
        # Get the start and end of the protein on the whole reference genome
        proteinStart <- integrase[1]+335
        proteinEnd <- integrase[2]+335
        
        # Extract the part of the sequence associated with that protein and add it to the result
        reducedSequence <- substr(elem[2],max(integrase[1]-location,1) , integrase[2] - location)
        # The protein name, sequence location, size, and the start and end of the protein is also returned
        resList[[length(resList)+1]] <- c("integrase", reducedSequence, location, nchar(reducedSequence), proteinStart, proteinEnd)
      }
      
      # If the sequence covers the integrase protein
      if(capsid[1]  <= location+length(elem[2]) && location <= capsid[2]){
        
        # Get the start and end of the protein on the whole reference genome
        proteinStart <- capsid[1]+335
        proteinEnd <- capsid[2]+335
        
        # Extract the part of the sequence associated with that protein and add it to the result
        reducedSequence <- substr(elem[2],max(capsid[1]-location,1) , capsid[2] - location)
        # The protein name, sequence location, size, and the start and end of the protein is also returned
        resList[[length(resList)+1]] <- c("capsid", reducedSequence, location, nchar(reducedSequence), proteinStart, proteinEnd)
      }
    }
    # If the gene is env
    else{
      
      # Extract the target protein that are interesting from a env gbk file
      peptides <- extract.gb(geneEnv,"mat_peptide")$NC_001802
      gp120 <- str_split_1(peptides[[1]]$mat_peptide[1,2], "..")
      gp41 <- str_split_1(peptides[[2]]$mat_peptide[1,2],"..")
      
      # If the sequence covers the gp120 protein
      if(gp120[1]  <= location+length(elem[2]) && location <= gp120[2]){
        
        # Get the start and end of the protein on the whole reference genome
        proteinStart <- gp120[1]+5770
        proteinEnd <- gp120[2]+5770
        
        # Extract the part of the sequence associated with that protein and add it to the result
        reducedSequence <- substr(elem[2], max(gp120[1]-location,1), gp120[2]-location)
        # The protein name, sequence location, size, and the start and end of the protein is also returned
        resList[[length(resList)+1]] <- c("gp120", reducedSequence, location, nchar(reducedSequence), proteinStart, proteinEnd)
      }
      
      # If the sequence covers the gp41 protein
      if(gp41[1] <= location+length(elem[2]) && location <= gp41[2]){
        
        # Get the start and end of the protein on the whole reference genome
        proteinStart <- gp41[1]+5770
        proteinEnd <- gp41[2]+5770
        
        # Extract the part of the sequence associated with that protein and add it to the result
        reducedSequence <- substr(elem[2], max(gp41[1]-location,1), gp41[2]-location)
        # The protein name, sequence location, size, and the start and end of the protein is also returned
        resList[[length(resList)+1]] <- c("gp41", reducedSequence, location, nchar(reducedSequence), proteinStart, proteinEnd)
      }
    }
  }
  # If no result has been found
  if(isEmpty(resList)) return(c("any", sequence, position, length(sequence), 1, 9181))
  return(resList)
}

executionTime <- function(time = NA){
  #' This function is executed twice, the first time to get the 
  #' current time (before execution), the second after execution (with the first time in parameter)
  #' @param time The time given at the first execution
  #' @return The current time if it exists
  
  if(is.na(time)) return(Sys.time())
  else return(Sys.time() - time)
}

getAlignmentDna <- function(dnaSeq, referenceSequence){
  #' This function tries to align a dna sequence with the main one
  #' @param dnaSeq The sequence we cant to align
  #' @param mainSequence The main sequence we want to compare with
  #' @return The position of the dna sequence on the hiv dna
  
  matrixSub <- substitutionMatrixBuilder(3,1,-3)
  return(pairwiseAlignment(pattern = dnaSeq, subject = referenceSequence, substitutionMatrix= matrixSub, type="local", gapOpening=-3, gapExtension=-1))
}

substitutionMatrixBuilder <- function(weightMatch, weightMissmatchTransition, weightMissmatchTranslation){
  #' This function returns a subsitution matrix based on the weight of a match, 
  #' a missmatch with a translation and a missmatch with a transition
  #' @param weightMatch The weight of a match
  #' @param weightMissmatchTransition The weight of a missmatch with a transition
  #' @param weightMissmatchTranslation The weight of a missmatch with a translation
  #' @return The subsitution matrix
  
  
  # Bases are listed based on IUPAC codes and the matrix is instanciated
  bases <- c("A","T","C","G","R","Y","S","W","K","M","B","D","H","V","N")
  res <- matrix(nrow = length(bases), ncol = length(bases))
  
  # Columns and rows are renamed according to the bases
  colnames(res) <- bases
  rownames(res) <- bases
  
  # A double loop to match each character with each other
  for(i in 1:length(bases)){
    for(j in i:length(bases)){# Because the matrix is symetrical only half of it is needed (and its diagonal)
      
      # This case if for match
      if((bases[i] == bases[j])&&grepl(pattern = "A|C|G|T", x= bases[i])){
        res[i,j] <- weightMatch 
      }
      
      # This case is for transition
      else if((grepl(pattern = "A|G|R", x = bases[i])&&grepl(pattern = "A|G|R", x = bases[j]))||(grepl(pattern = "C|T|Y", x = bases[i])&&grepl(pattern = "C|T|Y", x = bases[j]))){
        res[i,j] <- weightMissmatchTransition
        res[j,i] <- weightMissmatchTransition
      }
      
      # All of the other cases are translation
      else {
        res[i,j] <- weightMissmatchTranslation
        res[j,i] <- weightMissmatchTranslation
      }
    }
  }
  return(res)
}

translateDna <- function(strDna){
  #' This function transform a dna sequence into a vector of their corresponding proteins
  #' @param strDna The dna sequence string
  #' @return A vector of proteins (IUPAC code)
  charset <- unlist(strsplit(strDna, split = ""))
  res <- translate(charset)
  return(res)
}

getReadableMutation <- function(mutation){ 
  #' This function get the location of a given mutation associate it with its transformation, 
  #' treatment and the target of the treatment
  #' @param row The row of the mutation whose location we want to get
  #' @return A vector or sublist containing the data
  
  # Check if the resistance depends of a single mutation (ex: A111B)
  if(grepl(pattern = "^[A-Z][0-9]+[A-Z](\\/[A-Z])*(,[A-Z][0-9]+[A-Z](\\/[A-Z])*)*$", x = mutation)){
    
    # Separate each part of the mutation
    res <- gsub("[A-Z]|\\/",replacement = "",x= mutation)
    firstBase <- substr(mutation, 1, 1)
    lastBase <- gsub(firstBase, "",mutation)
    lastBase <- gsub(res, "" ,lastBase)
    replacedBases <- paste0(str_split_1(lastBase, "\\/"),sep = "" ,collapse = "")
    
    # Create a vector with the data in it
    mutationList<- c(mutation, firstBase, replacedBases, res)
    
    # Return the vector
    return(mutationList)
  }
  
  # Check if the resistance depends of multiple mutation (ex: 2(A111B/C/E, C222D,...))
  else if(grepl(pattern = "[0-9]*\\(([A-Z][0-9]+[A-Z](\\/[A-Z])*)(,?.(([A-Z][0-9]+[A-Z](\\/[A-Z])*)|([0-9]*\\(([A-Z][0-9]+[A-Z](\\/[A-Z])*)(,?.[A-Z][0-9]+[A-Z](\\/[A-Z])*)*\\))))*\\)", x = mutation)){
    
    # Take the number of required mutations for the resistance to be active and add it to a list another list of mutation
    if(grepl(pattern = "^[0-9]", x = mutation)){
      content <- substr(mutation, 3, nchar(mutation)-1)
      numberMut <-  substr(mutation, 1, 1)
      
      mutationList <- list()
      mutationList[[1]] <- numberMut
    }else{
      content <- mutation
      mutationList <- list()
    }
    
    # If there is not multiple required mutations in the required mutations
    if(!grepl(",[0-9]",content)){
      # Split the content of the parenthesis and add sublists for each mutation
      mutation <- unlist(strsplit(content, ","))
      subListMutation <- list()
      for(mut in mutation){
        subListMutation[[length(subListMutation)+1]] <- getReadableMutation(mut)
      }
      mutationList[[2]] <- subListMutation
      
    }
    # If the mutation looks like: 2(A123B, C456D, 2(E123F,...))
    else if(grepl(pattern = "^[A-Z]", x = mutation)){
      # While the pattern is still detected, remove the first mutation and add it as a  new sublist
      while(grepl(pattern= "^[A-Z]", x = mutation)){
        mutation <- unlist(regmatches(mutation, regexpr(",", mutation), invert = TRUE))
        mutationList[[length(mutationList)+1]] <- getReadableMutation(mutation[1])
        mutation <- mutation[2]
      }
      mutationList[[length(mutationList)+1]] <- getReadableMutation(mutation)
    }# In any other cases
    else{
      mutation <- content
      subtree <- getReadableMutation(mutation)
      mutationList[[length(mutationList)+1]] <- subtree
    }
    return(mutationList)
  }
  
  # Else if the resistance is due to an insertion/deletion of codon
  else{
    return(mutation)
  }
}

getReadableMutation2 <- function(mutation){
  
  
}