library(reticulate)
use_condaenv("synapse250", required = T)
synapseclient <- reticulate::import('synapseclient')
syn <- synapseclient$Synapse()
library(tidyverse)
library(pbapply)
syn$login()
library(dtexbuilder)


ids <- c(
  "syn12292393",  ##NTAP pNF NCATS single agent
  "syn12292601",  ##NF2 Synodos 10x10 NCATS combination
  "syn12293222",  ##NF2 Synodos 6x6 NCATS combination
  "syn12296219",  ##NF2 Synodos NCATS single agent MIPE 4.0
  "syn12297785")  ##NF2 Synodos NCATS single agent MIPE 1.0 Syn5,Syn1 only

foo <- map(ids, function(x){
  readr::read_csv(syn$get(x)$path, col_types = "c")
}) %>% set_names(ids)

data <- dplyr::bind_rows(foo)  

cmpd_names <- data %>% 
  select(compoundName) %>% 
  distinct()

##to avoid Error in the HTTP2 framing layer issue
httr::set_config(httr::config(http_version = 0))

inchikeys <- sapply(cmpd_names$compoundName, function(x){
  foo <- dtexbuilder::.convert_id_to_structure_pubchem(x, id_type = 'name', output_type = 'InChIKey')
})

cmpd_names$inchikey <- inchikeys

##### missing ncgc names 
missing_ncgc <- dplyr::filter(cmpd_names, is.na(inchikey)) %>% 
  dplyr::filter(!grepl('GPHR.+',compoundName)) %>% 
  rename(id = compoundName) %>% 
  select(-inchikey)

#none from these
# screening_ids <- c("syn5522642",
#                    "syn5522643",
#                    "syn5522644",
#                    "syn5522645",
#                    "syn5522646",
#                    "syn5522647",
#                    "syn5522649",
#                    "syn5522648",
#                    "syn8556314"
# )
#
# screening_ids <- c("syn5917344", "syn5917348")
# 
# dat <- lapply(screening_ids, function(x){
# 
#   foo <- read.csv(syn$get(x)$path, header = T, colClasses = "character")
# 
#   foo <- foo %>%
#     select(contains("ID"),matches("name")) %>%
#     set_names(c("id","compoundName"))
#   foo
# })
# dat <- bind_rows(dat) %>% unique()


#Synodos NF2 HTS 6x6
dat_1 <- syn$get('syn18457540')$path %>%
  readxl::read_excel() %>%
  select(id = 'RowSid',
         compoundName = 'RowName') %>%
  distinct

dat_2 <- syn$get('syn18457540')$path %>%
  readxl::read_excel() %>%
  select(id = 'ColSid',
         compoundName = 'ColName') %>%
 distinct


#Synodos NF2 HTS 10x10
dat_3 <- syn$get('syn18457543')$path %>%
  readxl::read_excel() %>%
  select(id = 'RowSid',
         compoundName = 'RowName') %>%
  distinct

dat_4 <- syn$get('syn18457543')$path %>%
  readxl::read_excel() %>%
  select(id = 'ColSid',
         compoundName = 'ColName') %>%
  distinct

#Synodos NF2 MIPE 1.0
screening_ids <- c("syn18457501", "syn18457503")
dat <- lapply(screening_ids, function(x){
  
  foo <- read.csv(syn$get(x)$path, header = T, colClasses = "character")
  
  foo <- foo %>%
    select(contains("ID"),matches("name")) %>%
    set_names(c("id","compoundName"))
  foo
})
dat_5 <- bind_rows(dat) %>% unique()

#Synodos NF2 single agent mipe 4.0
dat_6 <- syn$get('syn18457533')$path %>%
  read_tsv %>%
  select(id = 'Sample ID',
         compoundName = 'Sample Name') %>%
  distinct


dat <- bind_rows(dat_1, dat_2, dat_3, dat_4, dat_5, dat_6) %>% 
  distinct()
dat$id[duplicated(dat$id)]

dat_single <- dat %>% group_by(id) %>% slice(1)
dat_single$id[duplicated(dat_single$id)]

map_1 <- left_join(missing_ncgc, dat_single)

##remap the one "null" string...
map_1$compoundName[map_1$id=='NCGC00384479-01'] <- 'NCGC00384479-01'


convert_name_to_cid <- function (input_id, id_type = c("name"),
                                 output_type = c("original", 'parent')) {
  Sys.sleep(0.25)
  
  input <- URLencode(input_id)

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

cids <- sapply(map_1$compoundName, function(x){
  convert_name_to_cid(x, 'name', 'parent')
})

map_1$cid <- cids

map_1$cid[map_1$compoundName=="17-DMAG (Alvespimycin)"] <- convert_name_to_cid("Alvespimycin", "name", "parent")
map_1$cid[map_1$compoundName=="NCGC00384479-01"] <- NA ##may be able to fix this later if we can find a real name
map_1$cid[map_1$compoundName=="MLN-7243"] <- convert_name_to_cid("MLN7243", "name", "parent")

map_submit <- map_1 %>% 
  filter(!is.na(cid)) %>% 
  mutate(PUBCHEM_EXT_DATASOURCE_CID = as.numeric(cid), .keep = 'unused') %>%
  group_by(PUBCHEM_EXT_DATASOURCE_CID) %>% 
  pivot_longer(!PUBCHEM_EXT_DATASOURCE_CID, names_to = "trash", values_to = "synonyms") %>% 
  select(-trash) %>% 
  distinct() %>% 
  mutate(PUBCHEM_EXT_DATASOURCE_REGID = stringr::str_c("NF-OSI_", cur_group_id())) %>% 
  mutate(grpid = row_number()) %>% 
  pivot_wider(names_from = 'grpid', values_from = 'synonyms') %>% 
  ungroup() %>% 
  distinct() %>%
  arrange(PUBCHEM_EXT_DATASOURCE_CID)

#do i really have to format it like this??
colnames(map_submit)[colnames(map_submit)%in%c('1', '2', '3', '4')] <- 'PUBCHEM_SUBSTANCE_SYNONYM'

write_csv(map_submit, "submission_files/pubchem_submission_12062022.csv", na = "")

##map file for synapse/my use
map_for_synapse <- map_1 %>% 
  filter(!is.na(cid)) %>% 
  rename(screen_id = id, synonym = compoundName) %>% 
  mutate(inchikey = sapply(cid, function(x){
    dtexbuilder::.convert_id_to_structure_pubchem(x, id_type = 'cid', output_type = 'InChIKey')
  })) 

write_csv(map_for_synapse, 'synapse_files/synodos_nf2_screen_structures.csv')


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


sanitycheck <- map_for_synapse %>% 
  filter(!is.na(cid)) %>% 
  mutate(pubchem_title = sapply(cid, function(x){
    convert_cid_to_title(x,  output_type = 'Title')
  })) 
  

