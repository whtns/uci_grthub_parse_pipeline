#!/bin/bash
#SBATCH --job-name=integrate_%j
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=08:00:00
#SBATCH --error=snakemake_master.%j.err
#SBATCH --output=snakemake_master.%j.out

set -euo pipefail

# Resolve the directory this script lives in so paths are stable regardless of current working directory
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default parameters
MIN_GENES=200
RUN_INTEGRATION=true
# Defaults for papermill inspection parameters
GROUPBY_VAR="condition"
GROUP1_VALUE="differentiated"
GROUP2_VALUE="ctrl"
NOTEBOOK="src/inspect_integrated_anndata.ipynb"

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
			--groupby_var|--groupby-var)
				if [[ -z "${2:-}" ]]; then
					echo "Error: value required for $1" >&2
					exit 2
				fi
				GROUPBY_VAR="$2"
				shift 2
				;;
			--group1_value|--group1-value)
				if [[ -z "${2:-}" ]]; then
					echo "Error: value required for $1" >&2
					exit 2
				fi
				GROUP1_VALUE="$2"
				shift 2
				;;
			--group2_value|--group2-value)
				if [[ -z "${2:-}" ]]; then
					echo "Error: value required for $1" >&2
					exit 2
				fi
				GROUP2_VALUE="$2"
				shift 2
				;;
			--notebook|--notebook-file)
				if [[ -z "${2:-}" ]]; then
					echo "Error: value required for $1" >&2
					exit 2
				fi
				NOTEBOOK="$2"
				shift 2
				;;
		--no-integration)
			# Do not run the python integration step; continue with downstream steps
			RUN_INTEGRATION=false
			shift
			;;
		--help|-h)
			echo "Usage: $0 [--min_genes N] [--groupby_var VAR] [--group1_value VAL] [--group2_value VAL] output_prefix [output_prefix ...]"
			echo
			echo "Runs src/harmony_integration.py for each provided output_prefix and then runs a papermill notebook for inspection."
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
	echo "Usage: $0 [--min_genes N] [--groupby_var VAR] [--group1_value VAL] [--group2_value VAL] output_prefix [output_prefix ...]" >&2
	exit 2
fi

. ~/.mymambainit-24.3.0
mamba activate scvi-tools

# Loop over all provided prefixes and run the integration script once per prefix.
for output_prefix in "$@"; do
	echo "---- Running integration for: $output_prefix (min_genes=$MIN_GENES, groupby_var=$GROUPBY_VAR, group1_value=$GROUP1_VALUE, group2_value=$GROUP2_VALUE) ----"
		if [[ "$RUN_INTEGRATION" == true ]]; then
			python src/parse_harmony_integration.py \
				--input_dir output/parse \
				--output_prefix "$output_prefix" \
				--min_genes "$MIN_GENES"
		else
			echo "Skipping python step for: $output_prefix (--no-python was set)"
		fi
    
    slug=$(basename "$output_prefix")
    DIRNAME=$(dirname "$output_prefix")

	# Run papermill and nbconvert but don't let failures stop the entire job.
	# Use explicit if/else so the script continues even if one of these steps fails.
	if papermill "$NOTEBOOK" \
		"$DIRNAME/inspect_integrated_anndata_${slug}.ipynb" \
		-p combined_adata_path "$DIRNAME/${slug}_harmony_integrated.h5ad" \
		-p groupby_var "$GROUPBY_VAR" \
		-p group1_value "$GROUP1_VALUE" \
		-p group2_value "$GROUP2_VALUE"; then
		echo "Papermill succeeded for ${slug}"
	else
		echo "Warning: papermill failed for ${slug}" >&2
	fi

	if jupyter nbconvert --to webpdf --no-input \
		"$DIRNAME/inspect_integrated_anndata_${slug}.ipynb"; then
		echo "nbconvert succeeded for ${slug}"
	else
		echo "Warning: nbconvert failed for ${slug}" >&2
	fi

	echo "---- Finished: $output_prefix ----"
done
