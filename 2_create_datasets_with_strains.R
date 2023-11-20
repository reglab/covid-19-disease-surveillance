library(tidyverse)
library(lubridate)
library(progress)
library(stringdist)
library(forcats)
library(pkgcond)
library(sf)
source("1_data_helper.R")

genomic_edges <- read_csv("[REDACTED]") # cannot be shared
cc_to_strains <- read_csv("[REDACTED]") # cannot be shared

# Create CalCONNECT file with strains + cluster information
address_map <- get_latest_address_map() %>% 
  select(-ZIP) %>% st_set_geometry(NULL)
melissa <- get_melissa()
parcels <- read_rds("[REDACTED]") %>%  # cannot be shared
  st_set_geometry(NULL) %>% 
  transmute(
    MELISSA_ID = ID,
    STRING_ADDRESS = str_squish(paste(STREET_NUM_NUMERIC, STREET_NUM_NON_NUMERIC, ADDR, UNIT, ZIP)),
    PARCEL_ID = PARCEL_ID
  )
melissa <- melissa %>% 
  left_join(parcels, by=c("MAK" = "MELISSA_ID")) %>% 
  select(
    MELISSA_ID = MAK,
    APN,
    PARCEL_ID
  )

calconnect_cases <- get_jan12_calconnect() %>% 
  filter(
    Record.Type == "COVID-19 Case"
  ) %>% mutate(
      first.name = Person.Account..First.Name.. %>% normalize_str(),
      last.name = Person.Account..Last.Name.. %>% normalize_str(),
      dob = Person.Account..Birthdate.. %>% parse_date_time(orders=c("mdy")),
      coll.date = `Date.Specimen.Collected` %>% parse_date_time(orders=c("mdy")),
      first.pos = Date.First.Positive %>% parse_date_time(orders=c("mdy")),
      zip = Person.Account..Mailing.Zip.Postal.Code.. %>% str_extract("^[0-9]{5}")
  )
calconnect_supplement <- get_latest_calconnect() %>% 
  select(Record.Number, Supervisor.Name...Phone.Number2, Occupation.Location, 
         Supervisor.Email.Address, Supervisor.Email.Address.2)

preprocessed_cc <- calconnect_cases %>% 
  left_join(calconnect_supplement, by=c("Record.Number")) 
preprocessed_cc <- preprocessed_cc %>% 
  match_to_schools() %>% 
  match_to_employers() %>% 
  match_to_congregate()
preprocessed_cc <- preprocessed_cc %>% 
  transmute(
    CC_RECORD_NUM = Record.Number,
    CALREDIE_ID = as.character(CalREDIE.Incident.ID..),
    PARENT_RECORD_NUM = Parent.Record.Number,
    RECORD_TYPE = Record.Type,
    OPENED_DATE = Date.Time.Opened %>% 
      str_extract("^.*,") %>% 
      as.Date(format="%m/%d/%Y,"),
    EPISODE_DATE = as.Date(CalCONNECT.Episode.Date, format="%m/%d/%Y"),
    FIRST_NAME = Person.Account..First.Name.. %>% normalize_str(),
    LAST_NAME = Person.Account..Last.Name.. %>% normalize_str(),
    DOB = Person.Account..Birthdate..,
    ZIP = Person.Account..Mailing.Zip.Postal.Code.. %>% substr(1, 5),
    PHONE = Person.Account..Home.Phone.. %>% str_replace_all("[^0-9]", ""),
    INTERVIEW_OUTCOME = Initial.Interview.Outcome,
    INTERVIEW_OUTCOME_WITH_DATE = paste0(
      Initial.Interview.Outcome, 
      if_else(is.na(Date.Time.Closed), "", paste0(", Closed ", Date.Time.Closed))) %>% str_replace("^NA$", ""),
    HOUSEHOLD = Household %>% normalize_str() %>% 
      str_replace_all("ASKED NOT TO CALL", "") %>% str_squish(),    
    CANDIDATE_HH_ID = str_extract(HOUSEHOLD, "#[0-9]{4,}$") %>% str_replace_all("#", ""),
    CONGREGATE_SETTING = Name.of.Congregate.Setting.. %>% normalize_str(),
    CONGREGATE_SETTING_2 = Name.of.Congregate.Setting.2.. %>% normalize_str(),
    MATCHED_CONG1 = MATCHED_CONG1,
    MATCHED_CONG2 = MATCHED_CONG2,
    SCHOOL_NAME = What.is.the.name.of.the.school... %>% normalize_str(),
    SCHOOL_CITY = Specify.school.location..city... %>% normalize_str() %>% 
      str_replace_all(", CA|UNITED STATES", "") %>% str_squish(),
    SCHOOL_TYPE = What.type.of.school... %>% normalize_str(),
    JOINED_SCHOOL = JOINED_SCHOOL,
    EMPLOYER_NAME = EMPLOYER_NAME,
    MESSY_EMPLOYER_NAME = Person..Employer.Name %>% normalize_str(),
    SUPERVISOR_NAME_PHONE = Supervisor.Name...Phone.Number. %>% normalize_str() %>% 
      str_replace_all("&", " ") %>% str_replace_all("\\(|\\)", "") %>% 
      str_replace_all("(?<=\\d) |[[:punct:]](?=\\d)", "") %>% 
      str_replace_all("[[:punct:]]", "") %>% str_squish(),
    SUPERVISOR_NAME_PHONE_2 = Supervisor.Name...Phone.Number2 %>% normalize_str() %>% 
      str_replace_all("&", " ") %>% str_replace_all("\\(|\\)", "") %>% 
      str_replace_all("(?<=\\d) |[[:punct:]](?=\\d)", "") %>% 
      str_replace_all("[[:punct:]]", "") %>% str_squish(),
    SUPERVISOR_EMAIL = Supervisor.Email.Address %>% normalize_str() %>% 
      if_else(str_detect(., "@"), ., as.character(NA)),
    SUPERVISOR_EMAIL_2 = Supervisor.Email.Address.2 %>% normalize_str() %>% 
      if_else(str_detect(., "@"), ., as.character(NA)),
    OCCUPATION_LOCATION = Occupation.Location %>% str_replace("\\[HT.*\\]", "") %>% 
      str_squish() %>% if_else(str_length(.) < 15, as.character(NA), .) %>% 
      if_else(str_detect(., "\\[HT"), as.character(NA), .) %>% 
      if_else(!str_detect(., "[0-9]"), as.character(NA), .) %>% normalize_str(),
    SPOT_FORM_ID = SPOT.Intake.Form..,
    GENDER = Person.Account..Gender..,
    RACE = Race..,
    ETHNICITY = Ethnicity..
  ) %>% 
  mutate(
    PHONE = if_else(str_detect(PHONE, "^[0-9]{10}$"), PHONE, as.character(NA)),
    HOUSEHOLD = case_when(
      str_length(CANDIDATE_HH_ID) > 1 & !is.na(CANDIDATE_HH_ID) ~ CANDIDATE_HH_ID,
      str_length(HOUSEHOLD) > 7 ~ HOUSEHOLD,
      TRUE ~ as.character(NA)
    ),
    SUPERVISOR_NAME_PHONE = case_when(
      grepl("UNABLE|UNWILLING|DECLINE|UNKNOWN|SELF|REFUSED|^NONE|NOT PROVIDED|DISCLOSED$|NOT PROVIDE", SUPERVISOR_NAME_PHONE) ~ as.character(NA),
      str_length(SUPERVISOR_NAME_PHONE) < 8 ~ as.character(NA),
      TRUE ~ SUPERVISOR_NAME_PHONE
    ),
    SUPERVISOR_NAME_PHONE_2 = case_when(
      grepl("UNABLE|UNWILLING|DECLINE|UNKNOWN|SELF|REFUSED|^NONE|NOT PROVIDED|DISCLOSED$|NOT PROVIDE", SUPERVISOR_NAME_PHONE_2) ~ as.character(NA),
      str_length(SUPERVISOR_NAME_PHONE) < 8 ~ as.character(NA),
      TRUE ~ SUPERVISOR_NAME_PHONE_2
    ),
    SUPERVISOR_PHONE = str_extract(SUPERVISOR_NAME_PHONE, "[\\d]{10,}") %>% 
      str_replace("^[01]", ""),
    SUPERVISOR_PHONE_2 = str_extract(SUPERVISOR_NAME_PHONE_2, "[\\d]{10,}") %>% 
      str_replace("^[01]", ""),
  ) %>% 
  mutate(
    CAL_ID = if_else(is.na(CALREDIE_ID), paste0("R", CC_RECORD_NUM), as.character(CALREDIE_ID))
  ) %>% 
  left_join(address_map, by=c("CAL_ID")) %>% 
  left_join(melissa, by=c("MELISSA_ID"))

with_strains <- preprocessed_cc %>% 
  left_join(cc_to_strains, by=c("CC_RECORD_NUM"))

cc_to_save_full <- with_strains %>% 
  filter(EPISODE_DATE >= as.Date("05/01/2021", format="%m/%d/%Y"), 
         EPISODE_DATE <= as.Date("12/31/2021", format="%m/%d/%Y"))

duplicated_record_nums <- cc_to_save_full %>% 
  group_by(CC_RECORD_NUM) %>%
  count() %>% 
  filter(n > 1) %>% 
  pull(CC_RECORD_NUM)
duplicated_records <- cc_to_save_full %>% 
  filter(CC_RECORD_NUM %in% duplicated_record_nums)
cc_to_save_full <- cc_to_save_full %>% 
  group_by(CC_RECORD_NUM) %>% 
  dplyr::slice(1) %>% ungroup()

cat("Sequencing rate for entire study period:", 1 - mean(is.na(cc_to_save_full$strain)), "\n")
cat("Sequences for entire study period:", sum(!is.na(cc_to_save_full$strain)), "\n")
write_csv(cc_to_save_full, "[REDACTED]")

## Gather and save exposure events
labeled_location_types <- read_csv("[REDACTED]") # cannot be shared
labeled_locations <- read_csv("[REDACTED]") # cannot be shared
exposure <- read_rds("[REDACTED]") %>% # cannot be shared
  select(
    CC_RECORD_NUM = `Record: Record Number`,
    EXPOSURE_ID = `Exposure Event: ID`,
    EXPOSURE_RECORD_NUM = `Exposure Event: Exposure Event Record #`,
    LOCATION = `Location`,
    LOCATION_TYPE = `Location Type`
  ) %>% 
  left_join(labeled_location_types, by=c("LOCATION_TYPE" = "Location.Type")) %>% 
  dplyr::rename(Location.Type.Label = Specific.Label) %>% 
  left_join(labeled_locations, by=c("LOCATION" = "Location")) %>% 
  dplyr::rename(Location.Label = Specific.Label) %>% 
  transmute(
    CC_RECORD_NUM = CC_RECORD_NUM,
    EXPOSURE_ID = EXPOSURE_ID,
    EXPOSURE_RECORD_NUM = EXPOSURE_RECORD_NUM,
    EXPOSURE_TYPE = if_else(is.na(Location.Type.Label), Location.Label, Location.Type.Label)
  ) %>% 
  left_join(preprocessed_cc, by=c("CC_RECORD_NUM")) %>% 
  select(CC_RECORD_NUM, EXPOSURE_ID, EXPOSURE_RECORD_NUM, EXPOSURE_TYPE, EPISODE_DATE) %>% 
  filter(EPISODE_DATE >= as.Date("05/01/2021", format="%m/%d/%Y"), 
         EPISODE_DATE <= as.Date("12/31/2021", format="%m/%d/%Y")) %>% 
  left_join(cc_to_strains, by=c("CC_RECORD_NUM"))
write_csv(exposure, "[REDACTED]")

# Gather and save location histories
calconnect_episode_dates <- cc_to_save_full %>% 
  select(CC_RECORD_NUM, EPISODE_DATE)
location_histories <- get_latest_location_histories() %>% 
  select(
    CC_RECORD_NUM = `Case Record ID`,
    LOCATION_NAME = `Location: Name +`,
    SPECIFIC_PLACE = `Specific Place in the Location`,
    LOCATION_ID = `Location: Account ID`,
    START = start_date,
    END = end_date
  ) %>%
  filter(
    START <= parse_date_time("12-31-2021", orders=c("mdy")),
    END >= parse_date_time("05-01-2021", orders = c("mdy"))
  ) %>% 
  mutate(
    START = difftime(START, parse_date_time("01/01/2021", orders=c("mdy")), units="days"),
    END = difftime(END, parse_date_time("01/01/2021", orders=c("mdy")), units="days")
  ) %>% 
  left_join(cc_to_strains, by=c("CC_RECORD_NUM")) %>% 
  left_join(calconnect_episode_dates, by=c("CC_RECORD_NUM")) %>% 
  filter(!is.na(EPISODE_DATE))

THRESHOLD <- 7
location_history_pairs <- location_histories %>%
  inner_join(location_histories, by=c("LOCATION_NAME")) %>% 
  filter(CC_RECORD_NUM.x != CC_RECORD_NUM.y) %>% 
  filter((START.x >= (START.y - THRESHOLD) & START.x <= (END.y + THRESHOLD)) | (START.y >= (START.x - THRESHOLD) & START.y <= (END.x + THRESHOLD))) %>% 
  mutate(time_diff = abs(as.numeric(difftime(EPISODE_DATE.x, EPISODE_DATE.y, units="days")))) %>% 
  filter(time_diff <= 14) %>% 
  rownames_to_column(var = "LOCATION_HISTORY_PAIR_ID") %>% 
  select(
    LOCATION_HISTORY_PAIR_ID,
    LOCATION_NAME,
    CC_RECORD_NUM.x,
    CC_RECORD_NUM.y,
    EPISODE_DATE.x,
    EPISODE_DATE.y,
    strain.x,
    strain.y,
    source.x,
    source.y
  )
write_csv(location_history_pairs, "[REDACTED]")

# Location history pairs -- requiring the "specific location" to match
location_history_specific_pairs <- location_histories %>%
  inner_join(location_histories, by=c("LOCATION_NAME", "SPECIFIC_PLACE")) %>% 
  filter(CC_RECORD_NUM.x != CC_RECORD_NUM.y) %>% 
  filter((START.x >= (START.y - THRESHOLD) & START.x <= (END.y + THRESHOLD)) | (START.y >= (START.x - THRESHOLD) & START.y <= (END.x + THRESHOLD))) %>% 
  mutate(time_diff = abs(as.numeric(difftime(EPISODE_DATE.x, EPISODE_DATE.y, units="days")))) %>% 
  filter(time_diff <= 14) %>% 
  rownames_to_column(var = "SPECIFIC_LOCATION_HISTORY_PAIR_ID") %>% 
  select(
    SPECIFIC_LOCATION_HISTORY_PAIR_ID,
    LOCATION_NAME,
    SPECIFIC_PLACE,
    CC_RECORD_NUM.x,
    CC_RECORD_NUM.y,
    EPISODE_DATE.x,
    EPISODE_DATE.y,
    strain.x,
    strain.y,
    source.x,
    source.y
  )
write_csv(location_history_specific_pairs, "[REDACTED]")

# Get case-contact links using Account ID
calconnect_contacts <- get_jan12_calconnect() %>% 
  filter(Record.Type == "COVID-19 Contact")
linked_cases <- calconnect_cases %>% # All cases that share an acct. ID with a contact
  filter(Account.ID %in% calconnect_contacts$Account.ID)
contact_to_case_map <- calconnect_contacts %>% # Join case to contact by Account ID
  transmute(
    Contact.Opened.Date = as.Date(Date.Time.Opened, "%m/%d/%Y"),
    Contact.Name = paste(Person.Account..First.Name.., Person.Account..Last.Name.., Person.Account..Birthdate..) %>% normalize_str(),
    Contact.Record = Record.Number,
    Account.ID
  ) %>% 
  left_join(
    linked_cases %>% 
      transmute(
        Case.Opened.Date = as.Date(Date.Time.Opened, "%m/%d/%Y"),
        Case.Name = paste(Person.Account..First.Name.., Person.Account..Last.Name.., Person.Account..Birthdate..),
        Case.Record = Record.Number,
        Account.ID
      )
  ) %>% 
  filter(
    Case.Opened.Date > Contact.Opened.Date,
    Case.Opened.Date - Contact.Opened.Date <= 21
  ) %>% 
  select(Contact.Record, Case.Record)

# Method based on name, DOB, ZIP, and phone
contact_to_case_map_2 <- calconnect_cases %>% 
  mutate(
    Person.Account..First.Name.. = Person.Account..First.Name.. %>% normalize_str(),
    Person.Account..Last.Name.. = Person.Account..Last.Name.. %>% normalize_str(),
    phone = if_else(is.na(Person.Account..Home.Phone..), Person.Account..Mobile.., Person.Account..Home.Phone..) %>% 
      str_replace_all("[^0-9]", "") %>% if_else(str_detect(., "^[0-9]{10}$"), ., as.character(NA))
  ) %>% 
  filter(
    !is.na(Person.Account..First.Name..), 
    !is.na(Person.Account..Last.Name..),
    !is.na(Person.Account..Birthdate..), 
    !is.na(Person.Account..Mailing.Zip.Postal.Code..),
    !is.na(phone)
  ) %>% 
  inner_join(calconnect_contacts %>% 
               mutate(
                 Person.Account..First.Name.. = Person.Account..First.Name.. %>% normalize_str(),
                 Person.Account..Last.Name.. = Person.Account..Last.Name.. %>% normalize_str(),
                 phone = if_else(is.na(Person.Account..Home.Phone..), Person.Account..Mobile.., Person.Account..Home.Phone..)
               ) %>% 
               filter(!is.na(Person.Account..First.Name..), !is.na(Person.Account..Last.Name..),
                      !is.na(Person.Account..Birthdate..), !is.na(Person.Account..Mailing.Zip.Postal.Code..),
                      !is.na(phone)), 
             by = c("Person.Account..First.Name..", "Person.Account..Last.Name..", "Person.Account..Birthdate..",
                    "phone", "Person.Account..Mailing.Zip.Postal.Code..")) %>%
  transmute(
    Name = paste(Person.Account..First.Name.., Person.Account..Last.Name.., Person.Account..Birthdate..),
    Case.Opened.Date = as.Date(Date.Time.Opened.x, "%m/%d/%Y"),
    Case.Record = Record.Number.x,
    Contact.Opened.Date = as.Date(Date.Time.Opened.y, "%m/%d/%Y"),
    Contact.Record = Record.Number.y,
  ) %>% 
  filter(
    Case.Opened.Date > Contact.Opened.Date,
    Case.Opened.Date - Contact.Opened.Date <= 21
  ) %>% 
  arrange(str_length(Name)) %>% 
  select(Contact.Record, Case.Record)

combined_contact_case_map <- bind_rows(contact_to_case_map, contact_to_case_map_2) %>% 
  mutate(edgename = paste(Contact.Record, Case.Record)) %>% 
  filter(!duplicated(edgename)) %>% 
  select(-edgename)

# Now get all case-contact pairs and swap out the contact ID for the case ID if it exists
calconnect_episode_dates <- calconnect_cases %>%
  filter(Record.Type == "COVID-19 Case") %>% 
  transmute(
    Record.Number = Record.Number,
    Episode.Date = as.Date(CalCONNECT.Episode.Date, "%m/%d/%Y")
  )
cases_linked_by_contacts <- calconnect_contacts %>% 
  select(Record.Number, Parent.Record.Number) %>% 
  filter(!is.na(Parent.Record.Number), Parent.Record.Number != 0) %>% 
  left_join(combined_contact_case_map, by=c("Record.Number" = "Contact.Record")) %>% 
  filter(!is.na(Case.Record)) %>% 
  select(
    Parent.Case = Parent.Record.Number,
    Linked.Case = Case.Record
  ) %>% rownames_to_column(var="CASE_CONTACT_PAIR") %>% 
  left_join(calconnect_episode_dates, by=c("Parent.Case" = "Record.Number")) %>% 
  dplyr::rename(Parent.Episode.Date = Episode.Date) %>% 
  left_join(calconnect_episode_dates, by=c("Linked.Case" = "Record.Number")) %>% 
  dplyr::rename(Linked.Episode.Date = Episode.Date) %>% 
  mutate(time_diff = abs(Parent.Episode.Date - Linked.Episode.Date)) %>% 
  filter(as.numeric(time_diff) < 15) %>% 
  left_join(cc_to_strains, by=c("Parent.Case" = "CC_RECORD_NUM")) %>% 
  left_join(cc_to_strains, by=c("Linked.Case" = "CC_RECORD_NUM")) %>% 
  filter(Parent.Episode.Date >= as.Date("05/01/2021", format="%m/%d/%Y")) %>% 
  filter(Parent.Episode.Date <= as.Date("12/31/2021", format="%m/%d/%Y")) %>% 
  filter(Linked.Episode.Date >= as.Date("05/01/2021", format="%m/%d/%Y")) %>% 
  filter(Linked.Episode.Date <= as.Date("12/31/2021", format="%m/%d/%Y"))
write_csv(cases_linked_by_contacts, "[REDACTED]")
