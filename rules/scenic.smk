
# Rule to subsample loom file to 10% of cells using loompy (Python)
rule subsample_loom:
    input:
        loom_file = f"{OUTPUT_DIR}/loom/{{sample}}.loom"
    output:
        subsampled_loom = f"{OUTPUT_DIR}/loom/{{sample}}_subsampled.loom"
    resources:
        mem_mb = 4000,
        cpus = 2
    log:
        f"{OUTPUT_DIR}/Rout/{{sample}}/subsample_loom.log"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/{{sample}}_subsample_loom.txt"
    threads: 2
    shell:
        """
        python scripts/subsample_loom.py \
            {input.loom_file} \
            {output.subsampled_loom} \
            --percentage 0.3
        """

# SCENIC rule for running the SCENIC analysis pipeline
rule runscenic:
    input:
        loom_file = f"{OUTPUT_DIR}/loom/{{sample}}.loom"
    output:
        loom_output_path = f"{OUTPUT_DIR}/scenic/{{sample}}.loom",
        loom_filtered_path = f"{OUTPUT_DIR}/scenic/{{sample}}_filtered.loom",
        report = f"{OUTPUT_DIR}/reports/{{sample}}_scenic_report.html"
    resources:
        mem_mb = 16000,  # 16GB for SCENIC
        partition = "standard"
    threads: 24
    log:
        f"{OUTPUT_DIR}/Rout/{{sample}}/scenic.Rout"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/{{sample}}_scenic.txt"
    params:
        TFs = config['TFs'],
        motifs = config['motifs'],
        feather_db = config['feather_db'],
        loom_filtered = f"{{sample}}_filtered.loom",
        loom_output = f"{{sample}}.loom"
    shell:
        """
        module load singularity/3.11.3
        nextflow run aertslab/SCENICprotocol -profile singularity \
        --loom_input {input.loom_file} \
        --loom_filtered_path {output.loom_filtered_path} \
        --loom_filtered {params.loom_filtered} \
        --loom_output_path {output.loom_output_path} \
        --loom_output {params.loom_output} \
        --TFs {params.TFs} --motifs {params.motifs} --db {params.feather_db} \
        --thr_min_genes 200 \
        --thr_min_cells 3 \
        --thr_n_genes 20000 \
        --thr_pct_mito 0.25 \
        --threads {threads} \
        -with-report {output.report}
        module unload singularity/3.11.3
        """

