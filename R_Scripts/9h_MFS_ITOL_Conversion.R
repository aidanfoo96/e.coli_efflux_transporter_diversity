# 9_MFS_ITOL_Conversion
library(tidyverse)

MFS_pfam_tbl <- read_tsv("pfam_annotations/MFS_pfam_tbl.txt")

MFS_pfam_tbl_tidy <- MFS_pfam_tbl %>%
  group_by(query_name) %>%
  mutate(target_name2 = make.unique(as.character(target_name))) %>%
  select(query_name, target_name, target_name2, tlen, qlen, env_from, env_to) %>%
  filter(target_name == "Sugar_tr")

MFS_pfam_tbl_tidy$query_name <- sub("$", ",", MFS_pfam_tbl_tidy$query_name)
MFS_pfam_tbl_tidy$qlen <- sub("$", ",", MFS_pfam_tbl_tidy$qlen)
MFS_pfam_tbl_tidy$env_from <- sub("$", "|", MFS_pfam_tbl_tidy$env_from)
MFS_pfam_tbl_tidy$env_to <- sub("$", "|", MFS_pfam_tbl_tidy$env_to)
MFS_pfam_tbl_tidy$target_name2 <- sub("$", ",", MFS_pfam_tbl_tidy$target_name2)

MFS_pfam_tbl_tidy$target_name3 <- MFS_pfam_tbl_tidy$target_name

MFS_pfam_tbl_tidy2 <- MFS_pfam_tbl_tidy %>% 
  mutate(target_name3 = fct_recode(target_name3, 
                                   "RE|" = "MFS_1",
                                   "HH|" = "MFS_2",
                                   "EL|" = "MFS_4", 
                                   "TR|" = "Sugar_tr"))



MFS_pfam_tbl_tidy2$col <- MFS_pfam_tbl_tidy$target_name

MFS_pfam_tbl_tidy3 <- MFS_pfam_tbl_tidy2 %>% 
  mutate(col = fct_recode(col,
                          "#ff0000|" = "MFS_1",
                          "#0000ff|" = "MFS_2",
                          "#00ff00|" = "MFS_4", 
                          "#ffff33|" = "Sugar_tr")) %>%
  select(!target_name )

MFS_pfam_tbl_tidy4 <- MFS_pfam_tbl_tidy3 %>%
  pivot_wider(names_from = target_name2, values_from = c(target_name2, tlen, env_from, env_to, target_name3, col))

MFS_pfam_tbl_tidy5 <- MFS_pfam_tbl_tidy4 %>%
  select(query_name, qlen, `target_name3_MFS_1,`, `env_from_MFS_1,`, `env_to_MFS_1,`, `col_MFS_1,`, `target_name2_MFS_1,`,
         `target_name3_MFS_2,`, `env_from_MFS_2,`, `env_to_MFS_2,`, `col_MFS_2,`, `target_name2_MFS_2,`, 
         `target_name3_MFS_4,`, `env_from_MFS_4,`, `env_to_MFS_4,`, `col_MFS_4,`, `target_name2_MFS_4,`, 
         `target_name3_Sugar_tr,`, `env_from_Sugar_tr,`, `env_to_Sugar_tr,`, `col_Sugar_tr,`, `target_name2_Sugar_tr,`)



write.table(MFS_pfam_tbl_tidy5, 'pfam_annotations/MFS_pfam_tbl_tidy5.txt', sep = "", row.names = FALSE, 
            na = "", col.names = FALSE, quote = FALSE)                                   
           


