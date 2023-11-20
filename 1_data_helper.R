library(tidyverse)
library(lubridate)
library(stringi)
library(readxl)

# Functions for transforming/cleaning data
na_if_blank <- function(x) {
  if_else(x == "", as.character(NA), x)
}

blank_if_na <- function(x) {
  if_else(is.na(x), "", x)
}

false_if_na <- function(x) {
  if_else(is.na(x), FALSE, x)
}

normalize_str <- function(x) {
  x %>% 
    stri_trans_general("Latin-ASCII") %>% 
    toupper() %>% 
    str_squish()
}

# Functions for reading data

# Get Location Histories
get_latest_location_histories <- function() {
  location_hist_folder <- "[REDACTED]" # cannot be shared
  location_hist_files <- list.files(location_hist_folder) %>% 
    as.data.frame() %>% 
    dplyr::rename("file_path" = ".") %>% 
    mutate(
      string_date = str_extract(file_path, "\\d\\d\\d\\d-\\d\\d-\\d\\d"),
      date = parse_date_time(string_date, orders=c("ymd"))
    ) %>% 
    arrange(desc(date))
  loc_hists <- read_rds(paste0(location_hist_folder, location_hist_files$file_path[1]))
  cat("Loaded location history file created on", location_hist_files$string_date[1], "\n")
  return(loc_hists)
}

# Get SCC Employee Data
get_scc_employees <- function() {
  read_csv("[REDACTED]") # cannot be shared
}

# Get latest CalCONNECT
get_latest_calconnect <- function() {
  cc_folder <- "[REDACTED]" # cannot be shared
  cc_files <- list.files(cc_folder) %>% 
    as.data.frame() %>% 
    dplyr::rename("file_path" = ".") %>% 
    mutate(
      string_date = str_extract(file_path, "\\d\\d_\\d\\d_\\d\\d\\d\\d"),
      date = parse_date_time(string_date, orders=c("mdy"))
    ) %>% 
    arrange(desc(date))
  cc <- read_csv(paste0(cc_folder, cc_files$file_path[1]), guess_max = 200000)
  print(names(cc))
  cc <- cc %>% 
    transmute(
      Date.Time.Opened = `Date/Time Opened`,
      Record.Date.Time.Last.Modified = strptime(`Record Date/Time Last Modified`, "%m/%d/%Y, %I:%M %p"),
      Record.Number = `Record Number`,
      Parent.Record.Number = `Parent Record Number`,
      Record.Type = `Record Type`,
      CalREDIE.Incident.ID.. = `CalREDIE Incident ID +`,
      Record.Owner = `Record Owner`,
      Initial.Interview.Outcome = `Initial Interview Outcome`,
      Person.Account..First.Name.. = `Person Account: First Name +`,
      Person.Account..Last.Name.. = `Person Account: Last Name +`,
      Person.Account..Birthdate.. = `Person Account: Birthdate +`,
      Person.Account..Home.Phone.. = `Person Account: Home Phone +`,
      Person.Account..Mobile.. = `Person Account: Mobile +`,
      Household,
      Person.Account..Mailing.Street.. = `Person Account: Mailing Street +`,
      Apartment.Number.. = `Mailing Apartment Number +`,
      Person.Account..Mailing.Zip.Postal.Code.. = `Person Account: Mailing Zip/Postal Code +`,
      Person.Account..Language.. = `Person Account: Language +`,
      Status,
      Date.Time.Closed = `Date/Time Closed`, 
      Closed.Reason = `Closed Reason`,
      Closed.Reason..If.Other..Please.Specify = `Closed Reason: If Other, Please Specify`,
      Interview.Completed.Date = `Interview Completed Date`,
      DO.YOU.HAVE.CLOSE.CONTACTS. = `DO YOU HAVE CLOSE CONTACTS?`, 
      If.yes..total.number.of.close.contacts. = `If yes, total number of close contacts?`,
      Person.Account..Gender.. = `Person Account: Gender +`,
      Ethnicity.. = `Ethnicity +`,
      Race.. = `Race +`,
      Asian.Specify = `Asian Specify`,
      Account.ID = `Account ID`,
      Record.ID = `Record ID`,
      Congregate.Setting.Type.. = `Congregate Setting Type +`, 
      RESIDENT.STAFF.IN.CONGREGATE.SETTING.. = `RESIDENT/STAFF IN CONGREGATE SETTING +`,
      CalREDIE.Person.ID.. = `CalREDIE Person ID +`,
      Symptomatic.Status.. = `Symptomatic Status +`,
      Symptom.Onset.Date.. = `Symptom Onset Date +`,
      Symptom.Resolution.Date.. = `Symptom Resolution Date +`,
      Person..Occupation.. = `Person: Occupation +`,
      Person..Describe.Specify.Occupation.. = `Person: Describe/Specify Occupation +`,
      Person..Occupation.Setting.. = `Person: Occupation Setting +`,
      Person..Occupation.Setting.Describe.Sett = `Person: Occupation Setting Describe/Sett`,
      Patient.affiliated.with.a.school... = `Patient affiliated with a school? +`,
      What.type.of.school... = `What type of school? +`,
      Grade,
      What.is.the.name.of.the.school... = `What is the name of the school? +`,
      Specify.school.location..city... = `Specify school location (city) +`,
      What.type.of.close.contact... = `What type of close contact? +`,
      Exposure.Source.Details = `Exposure Source Details`,
      Healthcare.work.location.. = `Healthcare work location +`,
      Congregate.Setting.Location.. = `Congregate Setting Location +`,
      Name.of.Congregate.Setting.. = `Name of Congregate Setting +`,
      Congregate.Setting.Location.2.. = `Congregate Setting Location 2 +`,
      Name.of.Congregate.Setting.2.. = `Name of Congregate Setting 2 +`,
      Exposure.to.a.variant.of.concern = `Exposure to a variant of concern`,
      Person..Employer.Name = `Person: Employer Name`,
      Originating.SPOT.Portal = `Originating SPOT Portal`,
      Date.First.Positive = `LR: Date of First Positive Result`,
      Date.Specimen.Collected = `LR: Date Specimen Collected (Recent Pos)`,
      Date.Specimen.Collected.Neg = `LR: (Recent Neg) Date Specimen Collected`,
      Variant.Identified = `LR: Variant Identified`,
      Variant.Identified.Other = `LR: Variant Identified - Other`,
      Supervisor.Name...Phone.Number. = `Supervisor Name & Phone Number+`,
      Supervisor.Name...Phone.Number2 = `Supervisor Name & Phone Number2`,
      Occupation.Location = `Occupation Location`,
      Supervisor.Email.Address = `Supervisor Email Address`,
      Supervisor.Email.Address.2 = `Supervisor Email Address 2`,
      CalCONNECT.Episode.Date = `CalCONNECT Episode Date`,
      SPOT.Intake.Form.. = `SPOT Intake Form #`
    )
  cat("Loaded CalCONNECT export from", cc_files$string_date[1], "\n")
  return(cc)
}

# Get latest CalCONNECT
get_jan12_calconnect <- function() {
  cc_folder <- "[REDACTED]" # cannot be shared
  cc_files <- list.files(cc_folder) %>% 
    as.data.frame() %>% 
    dplyr::rename("file_path" = ".") %>% 
    mutate(
      string_date = str_extract(file_path, "\\d\\d_\\d\\d_\\d\\d\\d\\d"),
      date = parse_date_time(string_date, orders=c("mdy"))
    ) %>% 
    filter(string_date == "01_12_2022") %>% 
    arrange(desc(date))
  cc <- read_csv(paste0(cc_folder, cc_files$file_path[1]), guess_max = 200000)
  print(names(cc))
  cc <- cc %>% 
    transmute(
      Date.Time.Opened = `Date/Time Opened`,
      Record.Date.Time.Last.Modified = strptime(`Record Date/Time Last Modified`, "%m/%d/%Y, %I:%M %p"),
      Record.Number = `Record Number`,
      Parent.Record.Number = `Parent Record Number`,
      Record.Type = `Record Type`,
      CalREDIE.Incident.ID.. = `CalREDIE Incident ID +`,
      Record.Owner = `Record Owner`,
      Initial.Interview.Outcome = `Initial Interview Outcome`,
      Person.Account..First.Name.. = `Person Account: First Name +`,
      Person.Account..Last.Name.. = `Person Account: Last Name +`,
      Person.Account..Birthdate.. = `Person Account: Birthdate +`,
      Person.Account..Home.Phone.. = `Person Account: Home Phone +`,
      Person.Account..Mobile.. = `Person Account: Mobile +`,
      Household,
      Person.Account..Mailing.Street.. = `Person Account: Mailing Street +`,
      Apartment.Number.. = `Mailing Apartment Number +`,
      Person.Account..Mailing.Zip.Postal.Code.. = `Person Account: Mailing Zip/Postal Code +`,
      Person.Account..Language.. = `Person Account: Language +`,
      Status,
      Date.Time.Closed = `Date/Time Closed`, 
      Closed.Reason = `Closed Reason`,
      Closed.Reason..If.Other..Please.Specify = `Closed Reason: If Other, Please Specify`,
      Interview.Completed.Date = `Interview Completed Date`,
      DO.YOU.HAVE.CLOSE.CONTACTS. = `DO YOU HAVE CLOSE CONTACTS?`, 
      If.yes..total.number.of.close.contacts. = `If yes, total number of close contacts?`,
      Person.Account..Gender.. = `Person Account: Gender +`,
      Ethnicity.. = `Ethnicity +`,
      Race.. = `Race +`,
      Asian.Specify = `Asian Specify`,
      Account.ID = `Account ID`,
      Record.ID = `Record ID`,
      Congregate.Setting.Type.. = `Congregate Setting Type +`, 
      RESIDENT.STAFF.IN.CONGREGATE.SETTING.. = `RESIDENT/STAFF IN CONGREGATE SETTING +`,
      CalREDIE.Person.ID.. = `CalREDIE Person ID +`,
      Symptomatic.Status.. = `Symptomatic Status +`,
      Symptom.Onset.Date.. = `Symptom Onset Date +`,
      Symptom.Resolution.Date.. = `Symptom Resolution Date +`,
      Person..Occupation.. = `Person: Occupation +`,
      Person..Describe.Specify.Occupation.. = `Person: Describe/Specify Occupation +`,
      Person..Occupation.Setting.. = `Person: Occupation Setting +`,
      Person..Occupation.Setting.Describe.Sett = `Person: Occupation Setting Describe/Sett`,
      Patient.affiliated.with.a.school... = `Patient affiliated with a school? +`,
      What.type.of.school... = `What type of school? +`,
      Grade,
      What.is.the.name.of.the.school... = `What is the name of the school? +`,
      Specify.school.location..city... = `Specify school location (city) +`,
      What.type.of.close.contact... = `What type of close contact? +`,
      Exposure.Source.Details = `Exposure Source Details`,
      Healthcare.work.location.. = `Healthcare work location +`,
      Congregate.Setting.Location.. = `Congregate Setting Location +`,
      Name.of.Congregate.Setting.. = `Name of Congregate Setting +`,
      Congregate.Setting.Location.2.. = `Congregate Setting Location 2 +`,
      Name.of.Congregate.Setting.2.. = `Name of Congregate Setting 2 +`,
      Exposure.to.a.variant.of.concern = `Exposure to a variant of concern`,
      Person..Employer.Name = `Person: Employer Name`,
      Originating.SPOT.Portal = `Originating SPOT Portal`,
      Date.First.Positive = `LR: Date of First Positive Result`,
      Date.Specimen.Collected = `LR: Date Specimen Collected (Recent Pos)`,
      Date.Specimen.Collected.Neg = `LR: (Recent Neg) Date Specimen Collected`,
      Variant.Identified = `LR: Variant Identified`,
      Variant.Identified.Other = `LR: Variant Identified - Other`,
      Supervisor.Name...Phone.Number. = `Supervisor Name & Phone Number+`,
      CalCONNECT.Episode.Date = `CalCONNECT Episode Date`,
      SPOT.Intake.Form.. = `SPOT Intake Form #`
    )
  cat("Loaded CalCONNECT export from", cc_files$string_date[1], "\n")
  return(cc)
}

# Get latest location account data
get_latest_location_accounts <- function() {
  location_acct_folder <- "[REDACTED]" # cannot be shared
  location_acct_files <- list.files(location_acct_folder) %>% 
    as.data.frame() %>% 
    dplyr::rename("file_path" = ".") %>% 
    mutate(
      string_date = str_extract(file_path, "\\d\\d\\d\\d-\\d\\d-\\d\\d"),
      date = parse_date_time(string_date, orders=c("ymd"))
    ) %>% 
    arrange(desc(date))
  loc_accts <- read_rds(paste0(location_acct_folder, location_acct_files$file_path[1]))
  cat("Loaded location accounts file created on", location_acct_files$string_date[1], "\n")
  return(loc_accts)
}

# Get latest PHAGE
get_latest_phage <- function() {
  phage.folder <- "[REDACTED]" # cannot be shared
  phage.files <- list.files(phage.folder) %>%
    as.data.frame() %>%
    dplyr::rename("file_path" = ".") %>%
    mutate(
      string_date = str_extract(file_path, "\\d\\d\\d\\d-\\d\\d-\\d\\d"),
      date = parse_date_time(string_date, orders=c("ymd"))
    ) %>% 
    arrange(desc(date))
  phage.data <- read_excel(paste0(phage.folder, phage.files$file_path[1]), guess_max=10500, sheet = "All Runs") %>% 
    mutate(
      `External ID` = ifelse(is.na(`External ID`), `Sequencing ID`, `External ID`)
    )
  cat("Loaded PHAGE file created on", phage.files$string_date[1], "\n")
  return(phage.data)
}

# Get latest Terra metadata
get_latest_terra_folder <- function() {
  terra.parent.folder <- "[REDACTED]" # cannot be shared
  terra.folders <- list.files(terra.parent.folder) %>% 
    as.data.frame %>% 
    dplyr::rename("file_path" = ".") %>% 
    mutate(
      string_date = str_extract(file_path, "\\d\\d\\d\\d\\d\\d\\d\\d"),
      date = parse_date_time(string_date, orders=c("mdy"))
    ) %>% 
    filter(!is.na(date)) %>% 
    arrange(desc(date))
  cat("Found Terra data from", terra.folders$string_date[1], "... \n")
  terra.folders$file_path[1] %>% paste0(terra.parent.folder, ., "/")
}

get_latest_terra_metadata <- function() {
  terra_folder <- get_latest_terra_folder()
  metadata <- read_tsv(paste0(terra_folder, "[REDACTED]"), guess_max=6000) # cannot be shared
  cat("Loaded Terra metadata file.", "\n")
  return(metadata)
}

# Get latest Nextstrain metadata
get_latest_nextstrain_metadata <- function() {
  terra_folder <- get_latest_terra_folder()
  nextstrain_metadata <- read_tsv(paste0(terra_folder, '[REDACTED]'), guess_max=6000) # cannot be shared
}

# Get latest genome aligned data
get_latest_genome_aligned_data <- function() {
  terra_folder <- get_latest_terra_folder()
  genome.aligned.data <- Biostrings::readDNAStringSet(paste0(terra_folder,'[REDACTED]')) # cannot be shared
}

# Get latest address map
get_latest_address_map <- function() {
  address_map_folder <- "[REDACTED]" # cannot be shared
  address.map.files <- list.files(address_map_folder) %>%
    as.data.frame() %>%
    dplyr::rename("file_path" = ".") %>%
    mutate(
      string_date = str_extract(file_path, "\\d\\d\\d\\d-\\d\\d-\\d\\d"),
      date = parse_date_time(string_date, orders="ymd")
    ) %>% 
    filter(!is.na(date)) %>% 
    arrange(desc(date))
  address_map <- read_rds(paste0(address_map_folder, address.map.files$file_path[1]))
  cat("Loaded address map created on", address.map.files$string_date[1], "\n")
  return(address_map)
}

# Get latest exposure file
get_latest_exposure <- function() {
  exposure_folder <- "[REDACTED]" # cannot be shared
  exposure.files <- list.files(exposure_folder) %>% 
    as.data.frame() %>% 
    dplyr::rename("file_path" = ".") %>% 
    mutate(
      string_date = str_extract(file_path, "\\d\\d\\d\\d-\\d\\d-\\d\\d"),
      date = parse_date_time(string_date, orders=c("ymd"))
    ) %>% 
    arrange(desc(date))
  exposure <- read_rds(paste0(exposure_folder, exposure.files$file_path[1]))
  cat("Loaded exposure file created on", exposure.files$string_date[1], "\n")
  return(exposure)
}

# Get schools master file
get_schools <- function() {
  clean_school_name <- function(x) {
    x %>% 
      stri_trans_general(id = "Latin-ASCII") %>% 
      gsub("preschool", "", ., ignore.case=T) %>% 
      gsub("[[:punct:]]", "", .) %>%
      gsub("elementary", "", ., ignore.case=T) %>% 
      gsub("remote", "", ., ignore.case=T) %>% 
      gsub(" classes", "", ., ignore.case=T) %>%
      gsub(" only", "", ., ignore.case=T) %>%
      gsub(" zoom", "", ., ignore.case=T) %>% 
      gsub("online", "", ., ignore.case=T) %>%
      gsub(" learning", "", ., ignore.case=T) %>% 
      gsub("preparatory", "PREP", ., ignore.case=T) %>% 
      gsub("middle", "", ., ignore.case=T) %>% 
      gsub(" high$| high ", "", ., ignore.case=T) %>% 
      gsub("school", "", ., ignore.case=T) %>% 
      gsub(" HS", "", ., ignore.case=T) %>% 
      toupper() %>% 
      str_squish() %>% 
      gsub("^LINCOLN", "ABRAHAM LINCOLN", .)
  }
  school_list <- read_csv("[REDACTED]") %>% # cannot be shared
    mutate(
      School.Name = clean_school_name(school_name),
      School.Name = case_when(
        School.Name == "LAURELWOOD" & city == "San Jose" ~ "LAURELWOOD SJ",
        School.Name == "MT PLEASANT" & grade == "HS" ~ "MT PLEASANT HS",
        School.Name == "CUPERTINO" & grade == "HS" ~ "CUPERTINO HS",
        School.Name == "CHRISTOPHER" & grade == "HS" ~ "CHRISTOPHER HS",
        School.Name == "SANTA TERESA" & grade == "HS" ~ "SANTA TERESA HS",
        School.Name == "ABRAHAM LINCOLN" & city == "San Jose" ~ "ABRAHAM LINCOLN HS",
        School.Name == "HERBERT HOOVER" & city == "San Jose" ~ "HERBERT HOOVER MS",
        School.Name == "JOHN MUIR" & city == "San Jose" ~ "JOHN MUIR MS",
        School.Name == "WILLOW GLEN" & grade == "H.S."  ~ "WILLOW GLEN HS",
        School.Name == "WILLOW GLEN" & grade == "Middle"  ~ "WILLOW GLEN MS",
        School.Name == "FOOTHILL" & city == "San Jose" ~ "FOOTHILL HS",
        School.Name == "SARATOGA" & city == "Los Gatos" ~ "SARATOGA HS",
        School.Name == "ALTA VISTA" & city == "M.V." ~ "ALTA VISTA HS",
        School.Name == "STRATFORD SAN JOSE" & grade == "MS" ~ "STRATFORD MS SAN JOSE",
        TRUE ~ School.Name
      )
    )
  school_list %>% select(School.Name, everything())
}

# Get congregate settings
get_congregate_list <- function() {
  read_rds("[REDACTED]") %>% # cannot be shared
    select(
      FACILITY_NAME = Facility_Name,
      FACILITY_TYPE = Facility_Type,
      FACILITY_ZIP = ZIP
    ) %>% filter(!duplicated(.))
}


# Get Melissa addresses
get_melissa <- function() {
  melissa <- read_rds("[REDACTED]") # cannot be shared
  return(melissa)
}

match_to_schools <- function(calconnect) {
  school_list <- get_schools()
  childcare_prek_family_schools <- school_list %>%
    filter(grade %in% c("Preschool", "Family", "Birth-MS"))
  
  elementary_schools <- school_list %>% 
    filter(grade %in% c("Elementary", "K-8th", "K-5th", "K-12th", "K-12", "K-8", "K", "K-6th", "TK-8", "K-4th", "K-7th", "K-2nd", "1st-12", "TK-*", "K-MS", "K-11th", "K-10h", "Birth-MS", "3rd-6th", "3rd-10th", "2nd-8th", "2nd-12th", "1st-8th"))
  
  middle_schools <- school_list %>% 
    filter(grade %in% c("Middle", "K-8th", "K-12th", "K-12", "MS", "K-8", "MS & HS", "K-6th", "6th-12th", "TK-8", "M.S. & H.S.", "K-7th", "1st-12", "TK-*", "MS &HS", "K-MS", "K-11th", "K-10h", "Birth-MS", "3rd-6th", "3rd-10th", "2nd-8th", "2nd-12th", "1st-8th"))
  
  high_schools <- school_list %>% 
    filter(grade %in% c("HS", "H.S.", "K-12th", "K-12", "MS & HS", "6th-12th", "M.S. & H.S.", "1st-12", "TK-*", "MS &HS", "K-11th", "K-10h", "3rd-10th", "2nd-12th"))
  
  postsecondary_adult_schools <- school_list %>% 
    filter(grade %in% c("Adult"))
  
  school_cases <- calconnect %>% 
    mutate(
      School.Affiliation = case_when(
        Patient.affiliated.with.a.school... == "Yes, student" ~ "STUDENT",
        Patient.affiliated.with.a.school... == "Yes, staff/faculty or volunteer" ~ "STAFF/FACULTY",
        TRUE ~ "NONE"
      ),
      Original.School.Name = What.is.the.name.of.the.school...,
      School.Name = clean_school_name(What.is.the.name.of.the.school...),
      School.Type = case_when(
        What.type.of.school... == "Elementary School" ~ "Elementary",
        What.type.of.school... == "Middle School" ~ "Middle",
        What.type.of.school... == "High School" ~ "High",
        What.type.of.school... %in% c("Vocational School or Community College", "College/University (Undergraduate or Graduate School)") ~ "Postsecondary",
        What.type.of.school... %in% c("Daycare/Preschool", "Other Childcare/After-school Care") ~ "Preschool_Daycare_Childcare",
        TRUE ~ "Unknown"
      ),
      School.Location = Specify.school.location..city... %>% normalize_str()
    ) %>% 
    filter(!is.na(School.Name))
  
  do_join <- function(x, y, threshold) {
    joined <- stringdist_left_join(x, y, method="jw", max_dist=threshold, by="School.Name") %>%
      mutate(dist = stringdist(School.Name.x, School.Name.y, method="jw")) %>% 
      group_by(Record.Number) %>% 
      arrange(dist) %>% 
      dplyr::slice(1) %>% ungroup()
  }
  
  calc_match_rate <- function(df) {
    match_rate <- 1 - (df$School.Name.y %>% is.na() %>% mean())
    cat("Match rate: ", match_rate, "\n")
    df
  }
  
  groups <- split(school_cases, school_cases$School.Type)
  school_threshold <- 0.15
  joined_child <- do_join(groups$Preschool_Daycare_Childcare, childcare_prek_family_schools, threshold=school_threshold) %>% calc_match_rate()
  joined_elementary <- do_join(groups$Elementary, elementary_schools, threshold=school_threshold) %>% calc_match_rate()
  joined_middle <- do_join(groups$Middle, middle_schools, threshold=school_threshold) %>% calc_match_rate()
  joined_high <- do_join(groups$High, high_schools, threshold=school_threshold) %>% calc_match_rate()
  joined_adult <- do_join(groups$Postsecondary, postsecondary_adult_schools, threshold=school_threshold) %>% calc_match_rate()
  joined_final <- rbind(
    joined_child,
    joined_elementary,
    joined_middle,
    joined_high,
    joined_adult
  ) %>% 
    calc_match_rate() %>% 
    mutate(
      Matched.School.Name = School.Name.y,
      Matched.School.Name = case_when(
        Matched.School.Name %in% c("ACE ESPERANZA","ACE EMPOWER ACADEMY") ~ ifelse(grepl("ESPERANZA",toupper(Original.School.Name)), "ACE ESPERANZA", "ACE EMPOWER ACADEMY"),
        Matched.School.Name %in% c("ACHIEVE KIDS SAN JOSE","ACHIEVEKIDS PALO ALTO") ~ ifelse(grepl("ALTO",toupper(School.Location)), "ACHIEVEKIDS PALO ALTO", "ACHIEVE KIDS SAN JOSE"),
        Matched.School.Name %in% c("ALMADEN","ALMOND") ~ ifelse(grepl("ALTOS",toupper(School.Location)), "ALMOND", "ALMADEN"),
        Matched.School.Name %in% c("BLOSSOM HILL","BLUE HILLS") ~ ifelse(grepl("GATOS",toupper(School.Location)), "BLOSSOM HILL", "BLUE HILLS"),
        Matched.School.Name %in% c("CASTILLEJA","CASTILLERO") ~ ifelse(grepl("ALTO",toupper(School.Location)), "CASTILLEJA", "CASTILLERO"),
        Matched.School.Name %in% c("CASTILLEJA","CASTILLERO") ~ ifelse(grepl("MILLER",toupper(Original.School.Name)), "MILLER", Matched.School.Name),
        Matched.School.Name %in% c("DISCOVERY CHARTER","DISCOVERY CHARTER II") ~ ifelse(grepl("II",toupper(Original.School.Name)), "DISCOVERY CHARTER II", "DISCOVERY CHARTER"),
        Matched.School.Name %in% c("ESTHER B CLARK PALO ALTO","ESTHER B CLARK S SOUTH BAY") ~ ifelse(grepl("ALTO",toupper(School.Location)), "ESTHER B CLARK PALO ALTO", "ESTHER B CLARK S SOUTH BAY"),
        Matched.School.Name %in% c("FAIRMEADOW","FAIRWOOD") ~ ifelse(grepl("ALTO",toupper(School.Location)), "FAIRMEADOW", "FAIRWOOD"),
        Matched.School.Name %in% c("FUSION ACADEMY PALO ALTO","FUSION ACADEMY SAN JOSE") ~ ifelse(grepl("ALTO",toupper(School.Location)), "FUSION ACADEMY PALO ALTO", "FUSION ACADEMY SAN JOSE"),
        Matched.School.Name %in% c("GARDNER BULLIS","GARDNER") ~ ifelse(grepl("ALTOS",toupper(School.Location)), "GARDNER BULLIS", "GARDNER"),
        Matched.School.Name %in% c("GEORGE C PAYNE","GEORGE MAYNE") ~ ifelse(grepl("JOSE",toupper(School.Location)), "GEORGE C PAYNE", "GEORGE MAYNE"),
        Matched.School.Name %in% c("GILROY","GILROY PREP A NAVIGATOR") ~ ifelse(grepl("PREP",toupper(Original.School.Name)), "GILROY PREP A NAVIGATOR", "GILROY"),
        Matched.School.Name %in% c("HORACE CURETON","HORACE MANN") ~ ifelse(grepl("CUR",toupper(Original.School.Name)), "HORACE CURETON", "HORACE MANN"),
        Matched.School.Name %in% c("LAKESIDE","LAKEWOOD") ~ ifelse(grepl("GATOS",toupper(School.Location)), "LAKESIDE","LAKEWOOD"),
        Matched.School.Name %in% c("LAURELWOOD") ~ ifelse(grepl("JOSE",toupper(School.Location)), "LAURELWOOD SJ","LAURELWOOD"),
        Matched.School.Name %in% c("LYNDALE","LYNHAVEN") ~ ifelse(grepl("D",toupper(Original.School.Name)), "LYNDALE","LYNHAVEN"),
        Matched.School.Name %in% c("MARSHALL LANE","MARSHALL POMEROY") ~ ifelse(grepl("JOSE",toupper(School.Location)), "MARSHALL LANE","MARSHALL POMEROY"),
        Matched.School.Name %in% c("ROBERT RANDALL","ROBERT SANDERS") ~ ifelse(grepl("JOSE",toupper(School.Location)), "ROBERT SANDERS","ROBERT RANDALL"),
        Matched.School.Name %in% c("RUCKER","RUSKIN") ~ ifelse(grepl("JOSE",toupper(School.Location)), "RUSKIN","RUCKER"),
        Matched.School.Name %in% c("SANTA CLARA COUNTY COMMUNITY","SANTA CLARA COUNTY COURT","SANTA CLARA COUNTY ROPSOUTH","SANTA CLARA COUNTY SPECIAL EDUCATION") ~ ifelse(grepl("GILROY",toupper(School.Location)), "SANTA CLARA COUNTY ROPSOUTH", ifelse(grepl("SPECIAL",Original.School.Name),"SANTA CLARA COUNTY SPECIAL EDUCATION", ifelse(grepl("COURT",Original.School.Name),"SANTA CLARA COUNTY COURT","SANTA CLARA COUNTY COMMUNITY"))),
        Matched.School.Name %in% c("WASHINGTON","WASHINGTON OPEN") ~ ifelse(grepl("JOSE",toupper(School.Location)), "WASHINGTON","WASHINGTON OPEN"),
        TRUE ~ Matched.School.Name
      )
    )
  
  calconnect <- calconnect %>% 
    left_join(joined_final %>% select(Record.Number, JOINED_SCHOOL = Matched.School.Name)
    )
  return(calconnect)
}

match_to_employers <- function(calconnect) {
  valid_employers <- get_latest_location_accounts() %>% 
    filter(grepl("axle", `Physical Description`, ignore.case=T)) %>% 
    transmute(LOCATION = `Name +` %>% normalize_str()) %>% 
    pull(LOCATION)
  calconnect <- calconnect %>% 
    mutate(
      EMPLOYER_NAME = if_else(normalize_str(Person..Employer.Name) %in% valid_employers,
                              normalize_str(Person..Employer.Name), as.character(NA))
    )
  return(calconnect)
}

match_to_congregate <- function(calconnect) {
  
  congregate_list <- get_congregate_list() %>% 
    rownames_to_column(var = "CONGREGATE_ID") %>% 
    select(CONGREGATE_ID, everything())
  cong_setting_1 <- calconnect %>% transmute(
    CC_RECORD_NUM = Record.Number,
    FACILITY_NAME = Name.of.Congregate.Setting.. %>% normalize_str(),
    FACILITY_LOCATION = Congregate.Setting.Location.. %>% normalize_str(),
    ZIP = Person.Account..Mailing.Zip.Postal.Code..
  ) %>% filter(!is.na(FACILITY_NAME))
  cong_setting_2 <- calconnect %>% transmute(
    CC_RECORD_NUM = Record.Number,
    FACILITY_NAME = Name.of.Congregate.Setting.2.. %>% normalize_str(),
    FACILITY_LOCATION = Congregate.Setting.Location.2.. %>% normalize_str(),
    ZIP = Person.Account..Mailing.Zip.Postal.Code..
  ) %>% filter(!is.na(FACILITY_NAME))
  
  joined_1 <- stringdist_left_join(cong_setting_1, congregate_list, by="FACILITY_NAME",
                                   method = "jw", max_dist = 0.15) %>% 
    filter(str_detect(ZIP, FACILITY_ZIP) | str_detect(FACILITY_LOCATION, FACILITY_ZIP)) %>% 
    mutate(score = stringdist(FACILITY_NAME.x, FACILITY_NAME.y, method="jw")) %>% 
    group_by(CC_RECORD_NUM) %>% 
    arrange(score) %>% 
    dplyr::slice(1) %>% 
    ungroup() %>% 
    select(CC_RECORD_NUM, MATCHED_CONG1 = CONGREGATE_ID)
  
  joined_2 <- stringdist_left_join(cong_setting_2, congregate_list, by="FACILITY_NAME",
                                   method = "jw", max_dist = 0.15) %>% 
    filter(str_detect(ZIP, FACILITY_ZIP) | str_detect(FACILITY_LOCATION, FACILITY_ZIP)) %>% 
    mutate(score = stringdist(FACILITY_NAME.x, FACILITY_NAME.y, method="jw")) %>% 
    group_by(CC_RECORD_NUM) %>% 
    arrange(score) %>% 
    dplyr::slice(1) %>% 
    ungroup() %>% 
    select(CC_RECORD_NUM, MATCHED_CONG2 = CONGREGATE_ID)
  
  calconnect %>% 
    left_join(joined_1, by=c("Record.Number" = "CC_RECORD_NUM")) %>% 
    left_join(joined_2, by=c("Record.Number" = "CC_RECORD_NUM"))
}