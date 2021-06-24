# 9_ABC_ITOL_Conversion

library(tidyverse)

ABC_pfam_tbl <- read_tsv("pfam_annotations/ABC_pfam_tbl.txt")

ABC_pfam_tbl_tidy <- ABC_pfam_tbl %>%
  group_by(query_name) %>%
  mutate(target_name2 = make.unique(as.character(target_name))) %>%
  select(query_name, target_name, target_name2, tlen, qlen, env_from, env_to, E_value) %>%
  filter(E_value < 1.5e-12)

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
                                   "HH|" = "MacB_PCD",
                                   "EL|" = "OB_MalK", 
                                   "HV|" = "AAA_21", 
                                   "DI|" = "NIL", 
                                   "PL|" = "oligo_HPY", 
                                   "PU|" = "SMC_N", 
                                   "OV|" = "TOBE_2"))

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
                          "#0000ff|" = "MacB_PCD",
                          "#00ff00|" = "OB_MalK", 
                          "#ff0f87|" = "AAA_21", 
                          "#000024|" = "NIL", 
                          "#5c5cff|" = "oligo_HPY",
                          "#006100|" = "SMC_N", 
                          "#ffff1f|" = "TOBE_2")) %>%
  select(!target_name )

ABC_pfam_tbl_tidy4 <- ABC_pfam_tbl_tidy3 %>%
  pivot_wider(names_from = target_name2, values_from = c(target_name2, tlen, env_from, env_to, target_name3, col))

ABC_pfam_tbl_tidy5 <- ABC_pfam_tbl_tidy4 %>%
  select(query_name, qlen, `target_name3_ABC_tran,`, `env_from_ABC_tran,`, `env_to_ABC_tran,`, `col_ABC_tran,`, `target_name2_ABC_tran,`,
         `target_name3_ABC_tran.1,`, `env_from_ABC_tran.1,`, `env_to_ABC_tran.1,`, `col_ABC_tran.1,`, `target_name2_ABC_tran.1,`,
         `target_name3_AAA_21,`, `env_from_AAA_21,`, `env_to_AAA_21,`, `col_AAA_21,`, `target_name2_AAA_21,`, 
         `target_name3_AAA_21.1,`, `env_from_AAA_21.1,`, `env_to_AAA_21.1,`, `col_AAA_21.1,`, `target_name2_AAA_21.1,`, 
         `target_name3_AAA_21.2,`, `env_from_AAA_21.2,`, `env_to_AAA_21.2,`, `col_AAA_21.2,`, `target_name2_AAA_21.2,`, 
         `target_name3_AAA_21.3,`, `env_from_AAA_21.3,`, `env_to_AAA_21.3,`, `col_AAA_21.3,`, `target_name2_AAA_21.3,`, 
         `target_name3_MacB_PCD,`, `env_from_MacB_PCD,`, `env_to_MacB_PCD,`, `col_MacB_PCD,`, `target_name2_MacB_PCD,`, 
         `target_name3_MacB_PCD.1,`, `env_from_MacB_PCD.1,`, `env_to_MacB_PCD.1,`, `col_MacB_PCD.1,`, `target_name2_MacB_PCD.1,`, 
         `target_name3_MacB_PCD.2,`, `env_from_MacB_PCD.2,`, `env_to_MacB_PCD.2,`, `col_MacB_PCD.2,`, `target_name2_MacB_PCD.2,`, 
         `target_name3_NIL,`, `env_from_NIL,`, `env_to_NIL,`, `col_NIL,`, `target_name2_NIL,`, 
         `target_name3_OB_MalK,`, `env_from_OB_MalK,`, `env_to_OB_MalK,`, `col_OB_MalK,`, `target_name2_OB_MalK,`, 
         `target_name3_OB_MalK.1,`, `env_from_OB_MalK.1,`, `env_to_OB_MalK.1,`, `col_OB_MalK.1,`, `target_name2_OB_MalK.1,`,  
         `target_name3_oligo_HPY,`, `env_from_oligo_HPY,`, `env_to_oligo_HPY,`, `col_oligo_HPY,`, `target_name2_oligo_HPY,`, 
         `target_name3_SMC_N,`, `env_from_SMC_N,`, `env_to_SMC_N,`, `col_SMC_N,`, `target_name2_SMC_N,`, 
         `target_name3_SMC_N.1,`, `env_from_SMC_N.1,`, `env_to_SMC_N.1,`, `col_SMC_N.1,`, `target_name2_SMC_N.1,`, 
         `target_name3_TOBE_2,`, `env_from_TOBE_2,`, `env_to_TOBE_2,`, `col_TOBE_2,`, `target_name2_TOBE_2,`) 

ACR_pfam_tbl_tidy5.2 <- ACR_pfam_tbl_tidy4 %>%
  select(query_name, qlen, `target_name3_ACR_tran,`, `env_from_ACR_tran,`, `env_to_ACR_tran,`, `col_ACR_tran,`, `target_name2_ACR_tran,`)


write.table(ABC_pfam_tbl_tidy5, 'pfam_annotations/ABC_pfam_tbl_tidy5.txt', sep = "", row.names = FALSE, 
            na = "", col.names = FALSE, quote = FALSE)                                   



