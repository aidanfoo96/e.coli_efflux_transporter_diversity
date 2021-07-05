Mate_pfam_tbl <- read_tsv("pfam_annotations/MatE_pfam_tbl.txt")

MatE_tidy_table <- Mate_pfam_tbl %>%
  group_by(query_name) %>%
  mutate(target_name2 = make.unique(as.character(target_name))) %>%
  select(query_name, target_name, target_name2, tlen, qlen, env_from, env_to) %>%
  filter(target_name == "MatE")

MatE_tidy_table$query_name <- sub("$", ",", MatE_tidy_table$query_name)
MatE_tidy_table$qlen <- sub("$", ",", MatE_tidy_table$qlen)
MatE_tidy_table$env_from <- sub("$", "|", MatE_tidy_table$env_from)
MatE_tidy_table$env_to <- sub("$", "|", MatE_tidy_table$env_to)
MatE_tidy_table$target_name2 <- sub("$", ",", MatE_tidy_table$target_name2)

MatE_tidy_table$target_name3 <- MatE_tidy_table$target_name

MatE_tidy_table3 <- MatE_tidy_table %>% 
  mutate(target_name3 = fct_recode(target_name3, 
                                   "RE|" = "MatE"))



MatE_tidy_table3$col <- MatE_tidy_table$target_name

MatE_tidy_table3 <- MatE_tidy_table3 %>% 
  mutate(col = fct_recode(col,
                          "#ff0000|" = "MatE")) %>%
  select(!target_name )

MatE_tidy_table4 <- MatE_tidy_table3 %>%
  pivot_wider(names_from = target_name2, values_from = c(target_name2, tlen, env_from, env_to, target_name3, col))

MatE_tidy_table5 <- MatE_tidy_table4 %>%
  select(query_name, qlen, `target_name3_MatE,`, `env_from_MatE,`, `env_to_MatE,`, `col_MatE,`, `target_name2_MatE,`,
         `target_name3_MatE.1,`, `env_from_MatE.1,`, `env_to_MatE.1,`, `col_MatE.1,`, `target_name2_MatE.1,`, 
         `target_name3_MatE.2,`, `env_from_MatE.2,`, `env_to_MatE.2,`, `col_MatE.2,`, `target_name2_MatE.2,`, 
         `target_name3_MatE.3,`, `env_from_MatE.3,`, `env_to_MatE.3,`, `col_MatE.3,`, `target_name2_MatE.3,`) 

write.table(MatE_tidy_table5, 'pfam_annotations/MatE_pfam_tbl_tidy5.txt', sep = "", row.names = FALSE, 
            na = "", col.names = FALSE, quote = FALSE)                                   
                                 



