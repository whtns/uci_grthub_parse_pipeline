# SCENIC rule for running the SCENIC analysis pipeline
# This rule uses Nextflow to run the SCENIC protocol with specified parameters
# The input is a loom file, and the output is a final loom file after processing
rule runscenic:
    input:
        loom_file = f"{OUTPUT_DIR}/loom/{{sample}}.loom"
    output:
        scenic_loom = f"{OUTPUT_DIR}/scenic/{{sample}}.loom"
    resources:
        mem_mb = 16000,  # 16GB for SCENIC
        cpus = 4,
        partition = "standard"
    log:
        f"{OUTPUT_DIR}/Rout/{{sample}}/scenic.Rout"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/{{sample}}_scenic.txt"
    params:
        TFs = config['TFs'],
        motifs = config['motifs'],
        feather_db = config['feather_db']
    shell:
        """
        module load singularity/3.11.3
        nextflow run aertslab/SCENICprotocol -profile singularity \
        --loom_input {input.loom_file} --loom_output {output.scenic_loom} \
        --TFs {params.TFs} --motifs {params.motifs} --db {params.feather_db} \
        --thr_min_genes 1
        module unload singularity/3.11.3
        """

