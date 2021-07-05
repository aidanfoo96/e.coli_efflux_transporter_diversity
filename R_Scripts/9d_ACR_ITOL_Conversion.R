library(tidyverse)

ACR_pfam_tbl <- read_tsv("pfam_annotations/ACR_pfam_tbl.txt")

ACR_pfam_tbl_tidy <- ACR_pfam_tbl %>%
  group_by(query_name) %>%
  mutate(target_name2 = make.unique(as.character(target_name))) %>%
  select(query_name, target_name, target_name2, tlen, qlen, env_from, env_to) %>%
  fil

ACR_pfam_tbl_tidy$query_name <- sub("$", ",", ACR_pfam_tbl_tidy$query_name)
ACR_pfam_tbl_tidy$qlen <- sub("$", ",", ACR_pfam_tbl_tidy$qlen)
ACR_pfam_tbl_tidy$env_from <- sub("$", "|", ACR_pfam_tbl_tidy$env_from)
ACR_pfam_tbl_tidy$env_to <- sub("$", "|", ACR_pfam_tbl_tidy$env_to)
ACR_pfam_tbl_tidy$target_name2 <- sub("$", ",", ACR_pfam_tbl_tidy$target_name2)

ACR_pfam_tbl_tidy$target_name3 <- ACR_pfam_tbl_tidy$target_name

ACR_pfam_tbl_tidy2 <- ACR_pfam_tbl_tidy %>% 
  mutate(target_name3 = fct_recode(target_name3, 
                                   "RE|" = "ACR_tran",
                                   "HH|" = "SecD_SecF",
                                   "EL|" = "MMPL"))



ACR_pfam_tbl_tidy2$col <- ACR_pfam_tbl_tidy$target_name

ACR_pfam_tbl_tidy3 <- ACR_pfam_tbl_tidy2 %>% 
  mutate(col = fct_recode(col,
                          "#ff0000|" = "ACR_tran",
                          "#0000ff|" = "SecD_SecF",
                          "#00ff00|" = "MMPL")) %>%
  select(!target_name )

ACR_pfam_tbl_tidy4 <- ACR_pfam_tbl_tidy3 %>%
  pivot_wider(names_from = target_name2, values_from = c(target_name2, tlen, env_from, env_to, target_name3, col))

ACR_pfam_tbl_tidy5 <- ACR_pfam_tbl_tidy4 %>%
  select(query_name, qlen, `target_name3_ACR_tran,`, `env_from_ACR_tran,`, `env_to_ACR_tran,`, `col_ACR_tran,`, `target_name2_ACR_tran,`,
         `target_name3_MMPL,`, `env_from_MMPL,`, `env_to_MMPL,`, `col_ACR_tran,`, `target_name2_MMPL,`, 
         `target_name3_MMPL.1,`, `env_from_MMPL.1,`, `env_to_MMPL.1,`, `col_MMPL.1,`, `target_name2_MMPL.1,`, 
         `target_name3_MMPL.2,`, `env_from_MMPL.2,`, `env_to_MMPL.2,`, `col_MMPL.2,`, `target_name2_MMPL.2,`, 
         `target_name3_MMPL.3,`, `env_from_MMPL.3,`, `env_to_MMPL.3,`, `col_MMPL.3,`, `target_name2_MMPL.3,`, 
         `target_name3_SecD_SecF,`, `env_from_SecD_SecF,`, `env_to_SecD_SecF,`, `col_SecD_SecF,`, `target_name2_SecD_SecF,`, 
         `target_name3_SecD_SecF.1,`, `env_from_SecD_SecF.1,`, `env_to_SecD_SecF.1,`, `col_SecD_SecF.1,`, `target_name2_SecD_SecF.1,`) 

ACR_pfam_tbl_tidy5.2 <- ACR_pfam_tbl_tidy4 %>%
  select(query_name, qlen, `target_name3_ACR_tran,`, `env_from_ACR_tran,`, `env_to_ACR_tran,`, `col_ACR_tran,`, `target_name2_ACR_tran,`)
         
write.csv(ACR_pfam_tbl_tidy5, 'pfam_annotations/ACR_pfam_tbl_tidy5.csv', na = "")

write.table(ACR_pfam_tbl_tidy5, 'pfam_annotations/ACR_pfam_tbl_tidy5.txt', sep = "", row.names = FALSE, 
            na = "", col.names = FALSE, quote = FALSE)                                   
write.table(ACR_pfam_tbl_tidy5.2, 'pfam_annotations/ACR_pfam_tbl_tidy5.2.txt', sep = "", row.names = FALSE, 
            na = "", col.names = FALSE, quote = FALSE)                                   



