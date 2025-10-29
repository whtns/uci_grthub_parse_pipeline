#!/bin/bash
#SBATCH --job-name=integrate_%j
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=6G
#SBATCH --time=08:00:00
#SBATCH --error=snakemake_master.%j.err
#SBATCH --output=snakemake_master.%j.out

set -euo pipefail

# Resolve the directory this script lives in so paths are stable regardless of current working directory
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default parameters
MIN_GENES=200
RUN_INTEGRATION=true

# Simple argument parsing: allow optional --min_genes N followed by one or more output_prefix arguments
while [[ $# -gt 0 ]]; do
	case "$1" in
		--min_genes|--min-genes)
			if [[ -z "${2:-}" ]]; then
				echo "Error: value required for $1" >&2
				exit 2
			fi
			MIN_GENES="$2"
			shift 2
			;;
		--no-integration)
			# Do not run the python integration step; continue with downstream steps
			RUN_INTEGRATION=false
			shift
			;;
		--help|-h)
			echo "Usage: $0 [--min_genes N] output_prefix [output_prefix ...]"
			echo
			echo "Runs src/harmony_integration.py for each provided output_prefix."
			echo
			echo "Example prefixes:" 
			echo "  output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017/control"
			echo "  output/scanpy/output-XETG00221__0069979__Adult_cohort__20251003__181017/heat_stress"
			echo "  output/scanpy/output-XETG00221__0069982__Aged_cohort__20251003__181018/control"
			echo "  output/scanpy/output-XETG00221__0069982__Aged_cohort__20251003__181018/heat_stress"
			exit 0
			;;
		-* )
			echo "Unknown option: $1" >&2
			exit 2
			;;
		* )
			break
			;;
	esac
done

if [[ $# -lt 1 ]]; then
	echo "Error: at least one output_prefix is required." >&2
	echo "Usage: $0 [--min_genes N] output_prefix [output_prefix ...]" >&2
	exit 2
fi

. ~/.mymambainit-24.3.0
mamba activate scvi-tools

# Loop over all provided prefixes and run the integration script once per prefix.
for output_prefix in "$@"; do
	echo "---- Running integration for: $output_prefix (min_genes=$MIN_GENES) ----"
		if [[ "$RUN_INTEGRATION" == true ]]; then
			python src/harmony_integration.py \
				--output_prefix "$output_prefix" \
				--min_genes "$MIN_GENES"
		else
			echo "Skipping python step for: $output_prefix (--no-python was set)"
		fi
    
    slug=$(basename "$output_prefix")
    DIRNAME=$(dirname "$output_prefix")

    papermill "src/inspect_integrated_anndata.ipynb" \
        "$DIRNAME/inspect_integrated_anndata_${slug}.ipynb" \
        -p combined_adata_path "$DIRNAME/${slug}_harmony_integrated.h5ad"
	
	jupyter nbconvert --to pdf --no-input \
	    "$DIRNAME/inspect_integrated_anndata_${slug}.ipynb"

	echo "---- Finished: $output_prefix ----"
done
