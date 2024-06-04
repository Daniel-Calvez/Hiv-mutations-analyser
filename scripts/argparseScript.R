suppressMessages(source("./functionsScrapping.R"))
suppressMessages(library(argparse))
mainScriptParser <- function(){
  parser <- ArgumentParser()
  parser$add_argument("file", nargs=1, help= "File to be displayed")
  parser$add_argument("-n", "--noscrap", action ="store_true", help= "If you prefer not using the web scrapping", default = FALSE)
  parser$add_argument("-v", "--verbose", action="store_true" , help= "If you want more text to be displayed", default = FALSE)
  parser$add_argument("-f", "--forcescrap", action = "store_true", help = "Force the web scrapping", default = FALSE)
  parser$add_argument("-s", "--switchHiv2", action = "store_true", help = "To analyse hiv2 instead of hiv1", default = FALSE)
  
  args <- parser$parse_args()
  file <- args$file
  noscrap <- args$noscrap
  forcescrap <- args$forcescrap
  verbose <- args$verbose
  Hiv2 <- args$switchHiv2
  
  if(file.access(file) == -1){
    stop(sprintf("Specified file ( %s ) does not exist", file))
  }
  if(forcescrap){
    scrap("https://hivfrenchresistance.org/hiv-french-resitance-tables-of-rules/")
    if(verbose) print("Data sucessfully scrapped")
  }else
  if(!noscrap){
    if(!is.na(file.info("../data/DatabaseHivResistance")$mtime) && difftime(Sys.time(),file.info("../data/DatabaseHivResistance")$mtime, units = "days") > 30){
      scrap("https://hivfrenchresistance.org/hiv-french-resitance-tables-of-rules/")
      if(verbose) print("Data sucessfully scrapped")
    }
  }
  return(list(file, verbose,Hiv2))
}