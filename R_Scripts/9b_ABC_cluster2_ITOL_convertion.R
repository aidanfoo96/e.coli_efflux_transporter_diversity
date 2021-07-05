library(tidyverse)
ABC_pfam_tbl <- read_tsv("pfam_annotations/run2/ABC_transporter_tbl2.txt")
ABC_pfam_tbl_tidy <- ABC_pfam_tbl %>%
  group_by(query_name) %>%
  mutate(target_name2 = make.unique(as.character(target_name))) %>%
  select(query_name, target_name, target_name2, tlen, qlen, env_from, env_to, e_value) 
ABC_pfam_tbl_tidy$query_name <- substr(ABC_pfam_tbl_tidy$query_name, 1, nchar(ABC_pfam_tbl_tidy$query_name)-2)

# import BLATS hits
ABC_pfam_tbl_tidy <- ABC_pfam_tbl_tidy %>%
  filter(query_name %in% cluster2_ABC_efflux$name)

ABC_pfam_tbl_tidy %>%
  group_by(target_name2) %>%
  summarise(n())


ABC_pfam_tbl_tidy$query_name <- sub("$", ",", ABC_pfam_tbl_tidy$query_name)
ABC_pfam_tbl_tidy$qlen <- sub("$", ",", ABC_pfam_tbl_tidy$qlen)
ABC_pfam_tbl_tidy$env_from <- sub("$", "|", ABC_pfam_tbl_tidy$env_from)
ABC_pfam_tbl_tidy$env_to <- sub("$", "|", ABC_pfam_tbl_tidy$env_to)
ABC_pfam_tbl_tidy$target_name2 <- sub("$", ",", ABC_pfam_tbl_tidy$target_name2)

ABC_pfam_tbl_tidy$target_name3 <- ABC_pfam_tbl_tidy$target_name

ABC_pfam_tbl_tidy2 <- ABC_pfam_tbl_tidy %>% 
  mutate(target_name3 = fct_recode(target_name3, 
                                   "RE|" = "ABC_tran",
                                   "HH|" = "ABC_membrane",
                                   "EL|" = "Peptidase_C39"))

ABC_pfam_tbl_tidy2 <- ABC_pfam_tbl_tidy %>% 
  mutate(target_name3 = fct_recode(target_name3, 
                                   "RE|" = "ABC_tran",
                                   "RE|" = "ABC_tran.1",
                                   "HH|" = "MacB_PCD",
                                   "HH|" = "MacB_PCD.1",
                                   "HH|" = "MacB_PCD.2"
                                   "EL|" = "OB_MalK", 
                                   "EL|" = "OB_MalK.1", 
                                   "HV|" = "AAA_21", 
                                   "HV|" = "AAA_21.1", 
                                   "HV|" = "AAA_21.2", 
                                   "HV|" = "AAA_21.3", 
                                   "DI|" = "NIL", 
                                   "PL|" = "oligo_HPY", 
                                   "PU|" = "SMC_N", 
                                   "PU|" = "SMC_N.1", 
                                   "OV|" = "TOBE_2"))



ABC_pfam_tbl_tidy2$col <- ABC_pfam_tbl_tidy$target_name

ABC_pfam_tbl_tidy3 <- ABC_pfam_tbl_tidy2 %>% 
  mutate(col = fct_recode(col,
                          "#ff0000|" = "ABC_tran",
                          "#0000ff|" = "ABC_membrane",
                          "#00ff00|" = "Peptidase_C39")) %>%
  select(!target_name )

ABC_pfam_tbl_tidy4 <- ABC_pfam_tbl_tidy3 %>%
  pivot_wider(names_from = target_name2, values_from = c(target_name2, tlen, env_from, env_to, target_name3, col))

ABC_pfam_tbl_tidy5 <- ABC_pfam_tbl_tidy4 %>%
  select(query_name, qlen, `target_name3_ABC_tran,`, `env_from_ABC_tran,`, `env_to_ABC_tran,`, `col_ABC_tran,`, `target_name2_ABC_tran,`,
         `target_name3_ABC_membrane,`, `env_from_ABC_membrane,`, `env_to_ABC_membrane,`, `col_ABC_membrane,`, `target_name2_ABC_membrane,`,
         `target_name3_Peptidase_C39,`, `env_from_Peptidase_C39,`, `env_to_Peptidase_C39,`, `col_Peptidase_C39,`, `target_name2_Peptidase_C39,`)

ACR_pfam_tbl_tidy5.2 <- ACR_pfam_tbl_tidy4 %>%
  select(query_name, qlen, `target_name3_ACR_tran,`, `env_from_ACR_tran,`, `env_to_ACR_tran,`, `col_ACR_tran,`, `target_name2_ACR_tran,`)


write.table(ABC_pfam_tbl_tidy5, 'pfam_annotations/run2/ABC_pfam_tbl_tidy5.txt', sep = "", row.names = FALSE, 
            na = "", col.names = FALSE, quote = FALSE)                                   



