zcat /dfs9/ucightf-lab/kstachel/TOOLS/REFERENCES/mouse/parse_split_pipe/Mus_musculus.GRCm38.102.gtf.gz data/ParseWPRE_BFP2_NC_mT_mGMIN_redoGestalt_Celltag_postWPREv5_staticPLpgk0823Josie.gtf.gz  > data/merged_mm10_reference.gtf

zcat /dfs9/ucightf-lab/kstachel/TOOLS/REFERENCES/mouse/parse_split_pipe/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz data/WPREDoxMinmTmGminREDOGestalt_Celltag_newV5postWPREstaticPLpgk0722.fa.gz > data/merged_mm10_reference.fa

split-pipe --mode mkref --fasta merged_mm10_reference.fa --genes merged_mm10_reference.gtf --genome_name GRCm38 --nthreads 8 --output_dir merged_mm10_reference