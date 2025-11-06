library(fs)

allele_df_files = dir_ls("/dataVolume/storage/single_cell_projects/resources/yang_et_al_proj/output/numbat", glob = "*_allele_counts.tsv.gz") %>% 
    identity()
              
test0 <- 
    allele_df_files %>% 
    map(read_and_check_allele_df)


read_and_check_allele_df <- 
    function(allele_df_file){
        allele_df_table <- data.table::fread(allele_df_file)
        
        table(allele_df_table$CHROM)
    }
