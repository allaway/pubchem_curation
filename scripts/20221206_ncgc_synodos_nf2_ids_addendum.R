library(reticulate)
use_condaenv("synapse250", required = T)
synapseclient <- reticulate::import('synapseclient')
syn <- synapseclient$Synapse()
library(tidyverse)
syn$login()
library(dtexbuilder)

convert_name_to_cid <- function (input_id, id_type = c("name"),
                                 output_type = c("parent")) {
  Sys.sleep(0.25)
  
  input <- URLencode(input_id)
  
  statement <- glue::glue("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{id_type}/{input}/cids/XML?cids_type=parent")
  
  res <- httr::with_config(httr::config(http_version = 0), httr::GET(statement), override = T)
  
  if (res$status_code == 200) {
    res_2 <- XML::xmlToList(rawToChar(res$content))
    response <- res_2$CID 
    if (is.null(response)) {
      response <- NA
    }
  } else {
    message(glue::glue("input \"{input_id}\" appears to be invalid"))
    response <- NA
  }
  response
}

map_submit <- tibble(
  'a' = "NCGC00384479-01", 
  'b' = "E7107",
  'PUBCHEM_EXT_DATASOURCE_CID' = convert_name_to_cid("E7107","name","parent"),
  'PUBCHEM_EXT_DATASOURCE_REGID' = "NF-OSI_72")

colnames(map_submit)[colnames(map_submit)%in%c('a', 'b')] <- 'PUBCHEM_SUBSTANCE_SYNONYM'

write_csv(map_submit, "submission_files/pubchem_submission_2_12062022.csv", na = "")






