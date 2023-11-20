source("2_create_datasets_with_strains.R")

# Functions to limit analysis based on source of genomic data
subset_strains <- function(df, genome_source) {
  if (genome_source == "All") {return(df)}
  df %>% mutate(
    strain = if_else(source == genome_source, strain, as.character(NA))
  ) 
}

subset_paired_strains <- function(df, genome_source) {
  if (genome_source == "All") {return(df)}
  df %>% mutate(
    strain.x = if_else(source.x == "Terra/PHAGE", strain.x, as.character(NA)),
    strain.y = if_else(source.y == "Terra/PHAGE", strain.y, as.character(NA))
  )
}

# Helper function: Given df and grouping variables, get all proposed edges.
# Gives both unsequenced edges (as pairs of CalCONNECT IDs), 
# and sequenced edges (as pairs of strains). 
get_edges <- function(df, grouping_vars) {
  filtered <- df %>% 
    filter_at(grouping_vars, all_vars(!is.na(.)))
  
  with_strain <- filtered %>%
    filter(!is.na(strain))
  
  self_joined_unsequenced <- filtered %>% 
    left_join(filtered, by=grouping_vars) %>% 
    filter(CC_RECORD_NUM.x < CC_RECORD_NUM.y)
  
  self_joined_sequenced <- with_strain %>% 
    left_join(with_strain, by=grouping_vars) %>% 
    filter(strain.x < strain.y)
  
  if ("EPISODE_DATE.x" %in% names(self_joined_sequenced)) {
    self_joined_unsequenced <- self_joined_unsequenced %>% 
      select(from = CC_RECORD_NUM.x, to = CC_RECORD_NUM.y,
             from_episode_date = EPISODE_DATE.x, to_episode_date = EPISODE_DATE.y)
    self_joined_sequenced <- self_joined_sequenced %>% 
      select(from = strain.x, to = strain.y,
             from_episode_date = EPISODE_DATE.x, to_episode_date = EPISODE_DATE.y)
  } else {
    self_joined_unsequenced <- self_joined_unsequenced %>% 
      select(from = CC_RECORD_NUM.x, to = CC_RECORD_NUM.y,
             from_episode_date = EPISODE_DATE, to_episode_date = EPISODE_DATE)
    self_joined_sequenced <- self_joined_sequenced %>% 
      select(from = strain.x, to = strain.y,
             from_episode_date = EPISODE_DATE, to_episode_date = EPISODE_DATE)
  }
  
  self_joined_unsequenced <- self_joined_unsequenced %>% 
    mutate(
      edgename = paste0(from, "->", to),
      time_diff = as.numeric(abs(from_episode_date - to_episode_date))
    )
  
  self_joined_sequenced <- self_joined_sequenced %>% 
    mutate(
      edgename = paste0(from, "->", to),
      time_diff = as.numeric(abs(from_episode_date - to_episode_date))
    )
  
  all_unsequenced_edges <- self_joined_unsequenced %>% 
    filter(time_diff < 15) %>% 
    group_by(edgename) %>% 
    dplyr::slice(1) %>% 
    ungroup()
  
  all_sequenced_edges <- self_joined_sequenced %>% 
    filter(time_diff < 15) %>% 
    group_by(edgename) %>% 
    dplyr::slice(1) %>% ungroup()
  
  result <- list()
  result$cases_assigned <- unique(c(all_unsequenced_edges$from, all_unsequenced_edges$to))
  result$cases_with_strain <- df %>% 
    filter(CC_RECORD_NUM %in% result$cases_assigned) %>% 
    pull(strain) %>% unique()
  result$all_unsequenced_edges <- all_unsequenced_edges
  result$all_sequenced_edges <- all_sequenced_edges
  
  result
}

# Helper to combine two tables of edges, then deduplicate by edgename
concat_and_dedup <- function(x, y) {
  bind_rows(x, y) %>% 
    group_by(edgename) %>% 
    dplyr::slice(1) %>% 
    ungroup()
}

# Combine two surveillance methods (union, intersection, or set difference)
combine_methods <- function(result1, result2, how="union") {
  result <- list()
  if (how == "union") {
    result$cases_assigned <- unique(c(result1$cases_assigned, result2$cases_assigned))
    result$cases_with_strain <- unique(c(result1$cases_with_strain, result2$cases_with_strain))
    result$all_unsequenced_edges <- concat_and_dedup(result1$all_unsequenced_edges, result2$all_unsequenced_edges)
    result$all_sequenced_edges <- concat_and_dedup(result1$all_sequenced_edges, result2$all_sequenced_edges)
  } else if (how == "intersection") {
    result$cases_assigned <- intersect(result1$cases_assigned, result2$cases_assigned)
    result$cases_with_strain <- intersect(result1$cases_with_strain, result2$cases_with_strain)
    result$all_unsequenced_edges <- result1$all_unsequenced_edges %>% 
      filter(edgename %in% result2$all_unsequenced_edges$edgename)
    result$all_sequenced_edges <- result1$all_sequenced_edges %>% 
      filter(edgename %in% result2$all_sequenced_edges$edgename)
  } else if (how == "setdiff") {
    result$cases_assigned <- setdiff(result1$cases_assigned, result2$cases_assigned)
    result$cases_with_strain <- setdiff(result1$cases_with_strain, result2$cases_with_strain)
    result$all_unsequenced_edges <- result1$all_unsequenced_edges %>% 
      filter(!(edgename %in% result2$all_unsequenced_edges$edgename))
    result$all_sequenced_edges <- result1$all_sequenced_edges %>% 
      filter(!(edgename %in% result2$all_sequenced_edges$edgename))
  }
  
  result  
}

# One of 'Terra/PHAGE', 'GISAID/Fulgent', or 'All'
GENOME_SOURCE <- "All"

# CalCONNECT case records
calconnect_file <- "[REDACTED]" # cannot be shared
calconnect_dataset <- read_csv(calconnect_file, guess_max=200000) %>% 
  mutate(EPISODE_WEEK = floor_date(EPISODE_DATE, unit="week")) %>% 
  subset_strains(GENOME_SOURCE)
cc_supplement <- read_csv("[REDACTED]") %>% # cannot be shared
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
  ) %>% 
  select(Record.Number, CONGREGATE_RESIDENT_OR_STAFF = `RESIDENT.STAFF.IN.CONGREGATE.SETTING..`)
calconnect_dataset <- calconnect_dataset %>% 
  left_join(cc_supplement, by=c("CC_RECORD_NUM" = "Record.Number"))
cat("Sequencing rate:", mean(!is.na(calconnect_dataset$strain)), "\n")
cat("Number of strains:", sum(!is.na(calconnect_dataset$strain)), "\n")

# Exposure events
generic_exposure_ids <- "[REDACTED]" # cannot be shared
exposure <- read_csv("[REDACTED]", guess_max=100000) %>% # cannot be shared
  subset_strains(GENOME_SOURCE) %>% 
  filter(!(as.numeric(EXPOSURE_RECORD_NUM) %in% generic_exposure_ids))

# Case-contact links
case_contact <- read_csv("[REDACTED]", guess_max = 100000) %>% # cannot be shared
  subset_paired_strains(GENOME_SOURCE)

# Location Histories
specific_location_histories <- read_csv("[REDACTED]", guess_max=100000) %>% # cannot be shared
  subset_paired_strains(GENOME_SOURCE)

# Get all types of proposed edges
# Baseline: Same CBG
same_cbg <- get_edges(calconnect_dataset, c("CBG"))

# Concept: People who live in the same household or nearby
household <- get_edges(calconnect_dataset,  c("HOUSEHOLD"))
same_address <- calconnect_dataset %>% filter(PRECISION != "NEARBY ADDRESS") %>% 
  get_edges(.,  c("MELISSA_ID"))
parcel <- get_edges(calconnect_dataset,  c("PARCEL_ID"))

spec_longer_loc_hist <- specific_location_histories %>% 
  select(SPECIFIC_LOCATION_HISTORY_PAIR_ID,
         CC_RECORD_NUM = CC_RECORD_NUM.x,
         strain = strain.x,
         EPISODE_DATE = EPISODE_DATE.x) %>% 
  bind_rows(specific_location_histories %>% select(
    SPECIFIC_LOCATION_HISTORY_PAIR_ID,
    CC_RECORD_NUM = CC_RECORD_NUM.y,
    strain = strain.y,
    EPISODE_DATE = EPISODE_DATE.y))
location_history <- get_edges(spec_longer_loc_hist,  c("SPECIFIC_LOCATION_HISTORY_PAIR_ID"))

# Concept: People who work together
work_exposure <- exposure %>% 
  filter(EXPOSURE_TYPE == "work") 
work_exposure <- get_edges(work_exposure,  c("EXPOSURE_ID"))

same_employer_name <- calconnect_dataset %>% 
  filter(
    !str_detect(MESSY_EMPLOYER_NAME, "DECLINE"),
    !str_detect(MESSY_EMPLOYER_NAME, "DISCLOSE"),
    !str_detect(MESSY_EMPLOYER_NAME, "SELF"),
    !str_detect(MESSY_EMPLOYER_NAME, "EMPLOYED"),
    !(MESSY_EMPLOYER_NAME %in% c(
      "N/A",  "UNKNOWN", "RETIRED", "NOT APPLICABLE",
      "NOT PROVIDED", "REFUSED", "NONE", "SELF"))
  )
same_employer_name <- get_edges(same_employer_name,  c("MESSY_EMPLOYER_NAME"))

stacked_supervisor_name <- bind_rows(
  calconnect_dataset %>% select(CC_RECORD_NUM, EPISODE_DATE, SUPERVISOR_NAME_PHONE, strain),
  calconnect_dataset %>% select(CC_RECORD_NUM, EPISODE_DATE, SUPERVISOR_NAME_PHONE = SUPERVISOR_NAME_PHONE_2, strain)
) %>% filter(!is.na(SUPERVISOR_NAME_PHONE))

same_supervisor_name <- stacked_supervisor_name %>% 
  filter(
    !str_detect(SUPERVISOR_NAME_PHONE, "DISCLOSE"),
    !str_detect(SUPERVISOR_NAME_PHONE, "APPLICABLE"),
    !str_detect(SUPERVISOR_NAME_PHONE, "UNEMPLOYED"),
    !str_detect(SUPERVISOR_NAME_PHONE, "AVAILABLE"),
    !str_detect(SUPERVISOR_NAME_PHONE, "DIDNT SHARE"),
    !str_detect(SUPERVISOR_NAME_PHONE, "RETIRED"),
    !str_detect(SUPERVISOR_NAME_PHONE, "DIDNT WANT TO SHARE")
  )
same_supervisor_name <- get_edges(same_supervisor_name,  c("SUPERVISOR_NAME_PHONE"))

stacked_supervisor_phone <- bind_rows(
  calconnect_dataset %>% select(CC_RECORD_NUM, EPISODE_DATE, SUPERVISOR_PHONE, strain),
  calconnect_dataset %>% select(CC_RECORD_NUM, EPISODE_DATE, SUPERVISOR_PHONE = SUPERVISOR_PHONE_2, strain)
) %>% filter(!is.na(SUPERVISOR_PHONE))
same_supervisor_phone <- get_edges(stacked_supervisor_phone,  c("SUPERVISOR_PHONE"))

stacked_supervisor_email <- bind_rows(
  calconnect_dataset %>% select(CC_RECORD_NUM, EPISODE_DATE, SUPERVISOR_EMAIL, strain),
  calconnect_dataset %>% select(CC_RECORD_NUM, EPISODE_DATE, SUPERVISOR_EMAIL = SUPERVISOR_EMAIL_2, strain)
) %>% 
  filter(!is.na(SUPERVISOR_EMAIL), str_detect(SUPERVISOR_EMAIL, ".@."))
same_supervisor_email <- get_edges(stacked_supervisor_email,  c("SUPERVISOR_EMAIL"))

same_supervisor <- combine_methods(same_supervisor_name, same_supervisor_phone) %>% 
  combine_methods(same_supervisor_email)

occ_location <- get_edges(calconnect_dataset,  c("OCCUPATION_LOCATION"))

# Concept: People who go to school together
school_exposure <- exposure %>% 
  filter(EXPOSURE_TYPE == "school") 
school_exposure <- get_edges(school_exposure,  c("EXPOSURE_ID"))

same_messy_school <- calconnect_dataset %>% 
  filter(SCHOOL_TYPE %in% c("ELEMENTARY SCHOOL", "MIDDLE SCHOOL", "HIGH SCHOOL")) %>% 
  filter(
    !str_detect(SCHOOL_NAME, "UNKNOWN"), 
    !str_detect(SCHOOL_NAME, "CANT FIND"),
    !str_detect(SCHOOL_NAME, "DISCLOSED"),
    !str_detect(SCHOOL_NAME, "N/A")) %>% 
  filter(
    !str_detect(SCHOOL_CITY, "ONLINE"),
    !str_detect(SCHOOL_CITY, "REMOTE")) %>% 
  mutate(
    SCHOOL_CITY = str_replace_all(SCHOOL_CITY, "[[:punct:]]", "")
  )
same_messy_school <- get_edges(same_messy_school,  c("SCHOOL_NAME", "SCHOOL_TYPE", "SCHOOL_CITY"))

# Concept: Same congregate settings
all_congregate_exposure <- exposure %>% 
  filter(EXPOSURE_TYPE %in% c("acute", "jail", "LTCF", "unhoused")) %>% 
  get_edges( c("EXPOSURE_ID"))

jails_exposure <- exposure %>% 
  filter(EXPOSURE_TYPE == "jail") %>% 
  get_edges( c("EXPOSURE_ID"))

ltcf_exposure <- exposure %>% 
  filter(EXPOSURE_TYPE == "LTCF") %>% 
  get_edges( c("EXPOSURE_ID"))

unhoused_exposure <- exposure %>% 
  filter(EXPOSURE_TYPE == "unhoused") %>% 
  get_edges( c("EXPOSURE_ID"))

acute_exposure <- exposure %>% 
  filter(EXPOSURE_TYPE == "acute") %>% 
  get_edges( c("EXPOSURE_ID"))

# Concept: Contact tracing
longer_case_contact <- case_contact %>% 
  select(CASE_CONTACT_PAIR,
         CC_RECORD_NUM = Parent.Case,
         strain = strain.x,
         EPISODE_DATE = Parent.Episode.Date) %>% 
  bind_rows(case_contact %>% select(
    CASE_CONTACT_PAIR,
    CC_RECORD_NUM = Linked.Case,
    strain = strain.y,
    EPISODE_DATE = Linked.Episode.Date
  ))
case_contact_edges <- get_edges(longer_case_contact,  c("CASE_CONTACT_PAIR"))

# Concept: Spatiotemporal Clustering
satscan_precise <- read_json("[REDACTED]", simplifyDataFrame=T) %>% # cannot be shared
  as_tibble() %>% 
  unnest_longer(col=points) %>% 
  unnest_wider(col=points) %>%
  dplyr::rename(`point_lng` = `...1`, `point_lat` = `...2`) %>% 
  mutate(
    start_date = as.Date(start_date, format="%Y/%m/%d"),
    end_date = as.Date(end_date, format="%Y/%m/%d")
  ) %>% 
  left_join(calconnect_dataset, by=c("point_lat" = "LAT", "point_lng" = "LON")) %>% 
  filter(EPISODE_DATE >= start_date, EPISODE_DATE <= end_date) %>% 
  get_edges(.,  c("id"))

# Other methods
ct_excluding_same_address <- combine_methods(case_contact_edges, same_address, how="setdiff")
satscan_excluding_same_address <- combine_methods(satscan_precise, same_address, how="setdiff")
multibuilding <- combine_methods(parcel, same_address, how="setdiff")
satscan_multibuilding <- combine_methods(satscan_precise, multibuilding, how="intersection")

ltcf_staff <- exposure %>% 
  left_join(calconnect_dataset %>% select(CC_RECORD_NUM, CONGREGATE_RESIDENT_OR_STAFF), by=c("CC_RECORD_NUM")) %>% 
  filter(EXPOSURE_TYPE == "LTCF") %>% 
  filter(CONGREGATE_RESIDENT_OR_STAFF %in% c("Yes, staff", "Yes, both")) %>% 
  get_edges(c("EXPOSURE_ID"))
ltcf_residents <- exposure %>% 
  left_join(calconnect_dataset %>% select(CC_RECORD_NUM, CONGREGATE_RESIDENT_OR_STAFF), by=c("CC_RECORD_NUM")) %>% 
  filter(EXPOSURE_TYPE == "LTCF") %>% 
  filter(CONGREGATE_RESIDENT_OR_STAFF %in% c("Yes, resident", "Yes, both")) %>% 
  get_edges(c("EXPOSURE_ID"))
ltcf_same_address <- combine_methods(ltcf_exposure, same_address, how="intersection")
ltcf_occ_location <- combine_methods(ltcf_exposure, occ_location, how="intersection")
ltcf_staff_combined <- combine_methods(ltcf_occ_location, ltcf_staff, how="union")
ltcf_residents_combined <- combine_methods(ltcf_same_address, ltcf_residents, how="union")

jails_staff <- exposure %>% 
  left_join(calconnect_dataset %>% select(CC_RECORD_NUM, CONGREGATE_RESIDENT_OR_STAFF), by=c("CC_RECORD_NUM")) %>% 
  filter(EXPOSURE_TYPE == "jail") %>% 
  filter(CONGREGATE_RESIDENT_OR_STAFF %in% c("Yes, staff", "Yes, both")) %>% 
  get_edges(c("EXPOSURE_ID"))
jails_residents <- exposure %>% 
  left_join(calconnect_dataset %>% select(CC_RECORD_NUM, CONGREGATE_RESIDENT_OR_STAFF), by=c("CC_RECORD_NUM")) %>% 
  filter(EXPOSURE_TYPE == "jail") %>% 
  filter(CONGREGATE_RESIDENT_OR_STAFF %in% c("Yes, resident", "Yes, both")) %>% 
  get_edges(c("EXPOSURE_ID"))
jails_same_address <- combine_methods(jails_exposure, same_address, how="intersection")
jails_occ_location <- combine_methods(jails_exposure, occ_location, how="intersection")
jails_staff_combined <- combine_methods(jails_occ_location, jails_staff, how="union")
jails_residents_combined <- combine_methods(jails_same_address, jails_residents, how="union")
jails_other <- combine_methods(jails_exposure, jails_staff_combined, how="setdiff") %>% 
  combine_methods(., jails_residents_combined, how="setdiff")

existing_methods <- c(
  "case_contact_edges", "school_exposure", "all_congregate_exposure", "jails_exposure",
  "ltcf_exposure", "acute_exposure", "unhoused_exposure", "work_exposure", "same_address", "same_cbg")
potential_methods <- c(
  "satscan_precise", "same_supervisor", "occ_location", "same_messy_school",
  "same_employer_name", "location_history", "same_cbg")
boutique_methods <- c(
  "case_contact_edges", "ct_excluding_same_address",
  "satscan_precise", "satscan_excluding_same_address", "multibuilding", "satscan_multibuilding",
  "ltcf_exposure", "ltcf_residents_combined", "ltcf_staff_combined",
  "jails_exposure", "jails_residents_combined", "jails_staff_combined", "jails_other",
  "same_cbg")

# All methods
list_of_methods <- c(
  "same_cbg", "household", "same_address", "parcel", "location_history",
  "work_exposure", "same_employer_name", "same_supervisor", "occ_location",
  "school_exposure", "same_messy_school", "all_congregate_exposure",
  "jails_exposure", "ltcf_exposure", "unhoused_exposure",
  "acute_exposure", "case_contact_edges", "satscan_precise",
  "ct_excluding_same_address", "satscan_excluding_same_address", "multibuilding", 
  "satscan_multibuilding", "ltcf_residents_combined", "ltcf_staff_combined",
  "jails_residents_combined", "jails_staff_combined", "jails_other")

# All proposed edges
all_proposed_edges <- bind_rows(
  case_contact_edges$all_sequenced_edges,
  school_exposure$all_sequenced_edges,
  all_congregate_exposure$all_sequenced_edges,
  jails_exposure$all_sequenced_edges,
  ltcf_exposure$all_sequenced_edges,
  acute_exposure$all_sequenced_edges,
  unhoused_exposure$all_sequenced_edges,
  work_exposure$all_sequenced_edges,
  same_address$all_sequenced_edges,
  same_cbg$all_sequenced_edges,
  satscan_precise$all_sequenced_edges,
  same_supervisor$all_sequenced_edges,
  occ_location$all_sequenced_edges,
  same_messy_school$all_sequenced_edges, 
  same_employer_name$all_sequenced_edges,
  location_history$all_sequenced_edges
) %>% 
  group_by(edgename) %>% 
  dplyr::slice(1) %>% 
  ungroup()
