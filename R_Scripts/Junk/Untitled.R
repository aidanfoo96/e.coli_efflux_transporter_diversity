library(tidyverse)

EamA_pfam_tbl <- read_tsv("pfam_annotations/EamA_pfam_tbl.txt")

EamA_pfam_tbl_tidy <- EamA_pfam_tbl %>%
  group_by(query_name) %>%
  mutate(target_name2 = make.unique(as.character(target_name))) %>%
  select(query_name, target_name, target_name2, tlen, qlen, env_from, env_to) 

EamA_pfam_tbl_tidy$query_name <- sub("$", ",", EamA_pfam_tbl_tidy$query_name)
EamA_pfam_tbl_tidy$qlen <- sub("$", ",", EamA_pfam_tbl_tidy$qlen)
EamA_pfam_tbl_tidy$env_from <- sub("$", "|", EamA_pfam_tbl_tidy$env_from)
EamA_pfam_tbl_tidy$env_to <- sub("$", "|", EamA_pfam_tbl_tidy$env_to)
EamA_pfam_tbl_tidy$target_name2 <- sub("$", ",", EamA_pfam_tbl_tidy$target_name2)

EamA_pfam_tbl_tidy$target_name3 <- EamA_pfam_tbl_tidy$target_name

EamA_pfam_tbl_tidy2 <- EamA_pfam_tbl_tidy %>% 
  mutate(target_name3 = fct_recode(target_name3, 
                                   "RE|" = "EamA",
                                   "HH|" = "Multi_Drug_Res"))



EamA_pfam_tbl_tidy2$col <- EamA_pfam_tbl_tidy$target_name

EamA_pfam_tbl_tidy3 <- EamA_pfam_tbl_tidy2 %>% 
  mutate(col = fct_recode(col,
                          "#ff0000|" = "EamA",
                          "#0000ff|" = "Multi_Drug_Res")) %>%
  select(!target_name )

EmaA_pfam_tbl_tidy4 <- EamA_pfam_tbl_tidy3 %>%
  pivot_wider(names_from = target_name2, values_from = c(target_name2, tlen, env_from, env_to, target_name3, col))

EamA_pfam_tbl_tidy5 <- EmaA_pfam_tbl_tidy4 %>%
  select(query_name, qlen, `target_name3_EamA,`, `env_from_EamA,`, `env_to_EamA,`, `col_EamA,`, `target_name2_EamA,`,
         `target_name3_EamA.1,`, `env_from_EamA.1,`, `env_to_EamA.1,`, `col_EamA.1,`, `target_name2_EamA.1,`, 
         `target_name3_Multi_Drug_Res,`, `env_from_Multi_Drug_Res,`, `env_to_Multi_Drug_Res,`, `col_Multi_Drug_Res,`, `target_name2_Multi_Drug_Res,`) 
         
write.table(EamA_pfam_tbl_tidy5, 'pfam_annotations/EamA_pfam_tbl_tidy5.txt', sep = "", row.names = FALSE, 
            na = "", col.names = FALSE, quote = FALSE)                                   
