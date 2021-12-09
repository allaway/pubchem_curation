library(reticulate)
use_condaenv("synapse250", required = T)
synapseclient <- reticulate::import('synapseclient')
syn <- synapseclient$Synapse()
library(tidyverse)
syn$login()
library(dtexbuilder)


##to avoid Error in the HTTP2 framing layer issue
httr::set_config(httr::config(http_version = 0))

convert_xref_to_sid <- function(input_id, id_type = c("RegistryID"),
                                output_type = c("original")) {
  Sys.sleep(0.25)
  
  input <- URLencode(input_id, repeated = T, reserved = T)
  
  statement <- glue::glue("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/xref/{id_type}/{input}/sids/XML")
  
  res <- httr::with_config(httr::config(http_version = 0), httr::GET(statement), override = T)
  
  if (res$status_code == 200) {
    res_2 <- XML::xmlToList(rawToChar(res$content))
    response <- res_2$SID
    if (is.null(response)) {
      response <- NA
    }
  } else {
    message(glue::glue("input \"{input_id}\" appears to be invalid"))
    response <- NA
  }
  response
}
convert_sid_to_cid <- function(input_id, id_type = c("sid"),
                                output_type = c('all', "standardized")) {
  Sys.sleep(0.25)
  
  input <- URLencode(input_id, repeated = T, reserved = T)
  
  statement <- glue::glue("https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/{id_type}/{input}/cids/XML?cids_type={output_type}")
  
  res <- httr::with_config(httr::config(http_version = 0), httr::GET(statement), override = T)
  
  if (res$status_code == 200) {
    res_2 <- XML::xmlToList(rawToChar(res$content))
    response <- res_2$Information$CID 
    if (is.null(response)) {
      response <- NA
    }
  } else {
    message(glue::glue("input \"{input_id}\" appears to be invalid"))
    response <- NA
  }
  response
}



convert_name_to_cid <- function(input_id, id_type = c("name"),
                                 output_type = c("original", 'parent')) {
  Sys.sleep(0.25)
  
  input <- URLencode(input_id, repeated = T, reserved = T)
  
  statement <- glue::glue("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{id_type}/{input}/cids/XML?cids_type=original")
  
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


convert_structure_to_cid <- function(input_id) {
  Sys.sleep(0.25)
  
  input <- URLencode(input_id, repeated = T, reserved = T)
  
  statement <- glue::glue("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/identity/smiles/{input}/XML")
  
  res <- httr::with_config(httr::config(http_version = 0), httr::GET(statement), override = T)
  
  listkey <- XML::xmlToList(rawToChar(res$content))$ListKey 
  
  statement_2 <- glue::glue("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{listkey}/cids/XML")
  
  res <- httr::with_config(httr::config(http_version = 0), httr::GET(statement_2), override = T)

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

secondary_screen <- readr::read_csv(syn$get('syn12292395')$path, col_types = "c")

id_to_name_map <- readr::read_csv(syn$get('syn11527695')$path, col_types = "c")

cmpd_names <- id_to_name_map %>% 
  select(`Compound ID`, `Drug Name`, Name, `Supplier ID`) %>% 
  distinct() %>% 
  filter(`Compound ID` %in% secondary_screen$compoundName) 

cmpd_names_1 <- cmpd_names %>% 
  mutate(compoundName = Name) %>% 
  filter(!is.na(compoundName))

cids_1 <- sapply(cmpd_names_1$`compoundName`, function(x){
  httr::set_config(httr::config(http_version = 0))
  foo <- convert_name_to_cid(x, id_type = 'name')
})

cmpd_names_1$cid <- cids_1
cmpd_names_1_filt <- cmpd_names_1 %>% filter(!is.na(cid))

cmpd_names_2 <- cmpd_names %>% 
  mutate(compoundName = `Drug Name`) %>% 
  filter(!is.na(compoundName)) %>% 
  filter(!compoundName %in% cmpd_names_1_filt$`Drug Name`)

cids_2 <- sapply(cmpd_names_2$compoundName, function(x){
  httr::set_config(httr::config(http_version = 0))
  foo <- convert_name_to_cid(x, id_type = 'name')
})

cmpd_names_2$cid <- cids_2

cmpd_names_2_filt <- cmpd_names_2 %>% filter(!is.na(cid))

cmpd_names_3 <- cmpd_names %>% 
  mutate(compoundName = `Supplier ID`) %>% 
  filter(!is.na(compoundName)) %>% 
  filter(!Name %in% cmpd_names_1_filt$compoundName) %>% 
  filter(!`Drug Name` %in% cmpd_names_2_filt$compoundName)

cids_3 <- sapply(cmpd_names_3$compoundName, function(x){
  httr::set_config(httr::config(http_version = 0))
  foo <- convert_name_to_cid(x, id_type = 'name')
})

cmpd_names_3$cid <- cids_3
cmpd_names_3_filt <- cmpd_names_3 %>% filter(!is.na(cid))

cmpd_names_4 <- cmpd_names %>% 
  filter(!Name %in% cmpd_names_1_filt$compoundName) %>% 
  filter(!`Drug Name` %in% cmpd_names_2_filt$compoundName) %>% 
  filter(!`Supplier ID` %in% cmpd_names_3_filt$`Supplier ID`) %>% 
  filter(!is.na(`Supplier ID`))

cids_4 <- sapply(cmpd_names_4$`Supplier ID`, function(x){
  httr::set_config(httr::config(http_version = 0))
  foo <- convert_xref_to_sid(x, id_type = 'RegistryID')
  if(!is.na(foo)){
    foo <-  convert_sid_to_cid(as.character(foo), id_type = 'sid', output_type = 'standardized')
  }
  foo
})

cmpd_names_4$cid <- cids_4

cmpd_names_4_filt <- cmpd_names_4 %>% filter(!is.na(cid)) %>% 
  filter(!`Supplier ID` %in% c(1,01503223,21683,43308,401005))

cmpd_names_5 <- cmpd_names %>% 
  filter(!Name %in% cmpd_names_1_filt$compoundName) %>% 
  filter(!`Drug Name` %in% cmpd_names_2_filt$compoundName) %>% 
  filter(!`Supplier ID` %in% cmpd_names_3_filt$`Supplier ID`) %>% 
  filter(!`Supplier ID` %in% cmpd_names_4_filt$`Supplier ID`) %>% 
  filter(!is.na(`Name`))

#use SMILES 

cids_5 <- sapply(cmpd_names_5$Name, function(x){
  httr::set_config(httr::config(http_version = 0))
  convert_structure_to_cid(x)
})


c1 <- select(cmpd_names_1_filt, `Compound ID`,cid)

c2 <- select(cmpd_names_2_filt, `Compound ID`, cid)

c3 <- select(cmpd_names_3_filt, `Compound ID`,cid) 

c4 <- select(cmpd_names_4_filt, `Compound ID`, cid) 

all_c <- bind_rows(c1,c2,c3,c4) %>% unique()


convert_cid_to_title <- function (input_id, output_type = c("Title")) {
  Sys.sleep(0.25)
  
  input <- URLencode(input_id)
  
  statement <- glue::glue("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{input_id}/property/{output_type}/XML")
  
  res <- httr::with_config(httr::config(http_version = 0), httr::GET(statement), override = T)
  
  if (res$status_code == 200) {
    res_2 <- XML::xmlToList(rawToChar(res$content))
    response <- res_2$Properties$Title
    if (is.null(response)) {
      response <- NA
    }
  } else {
    message(glue::glue("input \"{input_id}\" appears to be invalid"))
    response <- NA
  }
  response
}

sanitycheck <- all_c %>% 
  filter(!is.na(cid)) %>% 
  mutate(pubchem_title = sapply(cid, function(x){
    convert_cid_to_title(x,  output_type = 'Title')
  })) 

all_c$pctitle <- sanitycheck

simple_cmpd <- cmpd_names %>% group_by(`Compound ID`) %>% slice(1)
mapped <- inner_join(simple_cmpd, all_c)

##Spot check suggests the first match is always (or nearly always, though I didn't find any examples otherwise) t
# the better match if there's more than one match, so...
mapped <- inner_join(simple_cmpd, all_c) %>% group_by(`Compound ID`) %>% slice(1) %>% 
  mutate(synonym = case_when(!is.na(Name) ~ Name,
                             is.na(Name) & !is.na(`Drug Name`) ~ `Drug Name`,
                             is.na(Name) & is.na(`Drug Name`) & !is.na(`Supplier ID`) ~ `Supplier ID`, 
                             !is.na(Name) & !is.na(`Drug Name`) & !is.na(`Supplier ID`) ~ `Compound ID`)) %>% 
  select(-Name, -`Drug Name`, -`Supplier ID`,-starts_with('pc')) %>% 
  rename(screen_id = `Compound ID`)

httr::set_config(httr::config(http_version = 0))

map_for_synapse <- mapped %>% 
  filter(!is.na(cid)) %>% 
  mutate(inchikey = sapply(cid, function(x){
    dtexbuilder::.convert_id_to_structure_pubchem(x, id_type = 'cid', output_type = 'InChIKey')
  })) 

write_csv(map_for_synapse,'synapse_files/synnf1_minnesota_screening_structures.csv')
syn$store(synapseclient$File('synapse_files/synnf1_minnesota_screening_structures.csv', parentId= 'syn26532680',
                             used = c('syn12292395','syn11527695'), executed = ))


# 
# 
# map_submit <- map_1 %>% 
#   filter(!is.na(cid)) %>% 
#   mutate(PUBCHEM_EXT_DATASOURCE_CID = as.numeric(cid), .keep = 'unused') %>%
#   group_by(PUBCHEM_EXT_DATASOURCE_CID) %>% 
#   pivot_longer(!PUBCHEM_EXT_DATASOURCE_CID, names_to = "trash", values_to = "synonyms") %>% 
#   select(-trash) %>% 
#   distinct() %>% 
#   mutate(PUBCHEM_EXT_DATASOURCE_REGID = stringr::str_c("NF-OSI_", cur_group_id())) %>% 
#   mutate(grpid = row_number()) %>% 
#   pivot_wider(names_from = 'grpid', values_from = 'synonyms') %>% 
#   ungroup() %>% 
#   distinct() %>%
#   arrange(PUBCHEM_EXT_DATASOURCE_CID)
# 
# #do i really have to format it like this??
# colnames(map_submit)[colnames(map_submit)%in%c('1', '2', '3', '4')] <- 'PUBCHEM_SUBSTANCE_SYNONYM'
# 
# write_csv(map_submit, "submission_files/pubchem_submission_12062022.csv", na = "")

##map file for synapse/my use
# 
# write_csv(map_for_synapse, 'synapse_files/synodos_nf2_screen_structures.csv')
# 





