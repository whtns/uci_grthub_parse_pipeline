# SCENIC rule for running the SCENIC analysis pipeline
# This rule uses Nextflow to run the SCENIC protocol with specified parameters
# The input is a loom file, and the output is a final loom file after processing
rule runscenic:
	input:
		loom_file = outputdir + "scenic/unfiltered.loom"
	output:
		final_loom = outputdir + "scenic/unfiltered-final.loom"
	log:
		outputdir + "Rout/scenic.Rout"
	benchmark:
		outputdir + "benchmarks/scenic.txt"
	params:
		TFs = config['TFs'],
		motifs = config['motifs'],
		feather_db = config['feather_db']
	shell:
	  "nextflow run aertslab/SCENICprotocol -profile docker "
	  "--loom_input {input.loom_file} --loom_output {output.final_loom} "
	  "--TFs {params.TFs} --motifs {params.motifs} --db {params.feather_db} "
	  "--thr_min_genes 1"
