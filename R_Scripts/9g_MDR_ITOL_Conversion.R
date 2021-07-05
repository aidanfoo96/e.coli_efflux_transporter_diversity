MDR_pfam_tbl <- read_tsv("pfam_annotations/MDR_pfam_tbl.txt")

MDR_tidy_table <- MDR_pfam_tbl %>%
  group_by(query_name) %>%
  mutate(target_name2 = make.unique(as.character(target_name))) %>%
  select(query_name, target_name, target_name2, tlen, qlen, env_from, env_to) %>%
  filter(target_name == "Multi_Drug_Res")

MDR_tidy_table$query_name <- sub("$", ",", MDR_tidy_table$query_name)
MDR_tidy_table$qlen <- sub("$", ",", MDR_tidy_table$qlen)
MDR_tidy_table$env_from <- sub("$", "|", MDR_tidy_table$env_from)
MDR_tidy_table$env_to <- sub("$", "|", MDR_tidy_table$env_to)
MDR_tidy_table$target_name2 <- sub("$", ",", MDR_tidy_table$target_name2)

MDR_tidy_table$target_name3 <- MDR_tidy_table$target_name

MDR_tidy_table3 <- MDR_tidy_table %>% 
  mutate(target_name3 = fct_recode(target_name3, 
                                   "RE|" = "Multi_Drug_Res"))

MDR_tidy_table3$col <- MDR_tidy_table$target_name

MDR_tidy_table3 <- MDR_tidy_table3 %>% 
  mutate(col = fct_recode(col,
                          "#ff0000|" = "Multi_Drug_Res")) %>%
  select(!target_name )

MDR_tidy_table4 <- MDR_tidy_table3 %>%
  pivot_wider(names_from = target_name2, values_from = c(target_name2, tlen, env_from, env_to, target_name3, col))

MDR_tidy_table5 <- MDR_tidy_table4 %>%
  select(query_name, qlen, `target_name3_Multi_Drug_Res,`, `env_from_Multi_Drug_Res,`, `env_to_Multi_Drug_Res,`, `col_Multi_Drug_Res,`, `target_name2_Multi_Drug_Res,`)

write.table(MDR_tidy_table5, 'pfam_annotations/MDR_pfam_tbl_tidy5.txt', sep = "", row.names = FALSE, 
            na = "", col.names = FALSE, quote = FALSE)                                   




