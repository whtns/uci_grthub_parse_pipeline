#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") [-r] <checksums.md5> [out.md5]

Read a checksums file (lines like: <md5>  <filename>) and rewrite filenames
that match the Parse-style sample pattern by prefixing them with the sample
directory (the sample is taken as the part before the first "_S<digits>").

Example:
  Input filename:
    0825I_44_BreuJ_Sublibrary4_S4_L001_R1_001.fastq.gz
  Output filename (default behaviour):
    0825I_44_BreuJ_Sublibrary4/0825I_44_BreuJ_Sublibrary4_S4_L001_R1_001.fastq.gz

Options:
  -r    Also replace R1/R2 with I1/I2 in the filename portion (optional).
        Use this only if you intentionally want R -> I mapping.
  -h    Show this help

The script writes to stdout by default or to [out.md5] if provided.
EOF
}

REPLACE_R_WITH_I=false
while getopts ":rh" opt; do
  case ${opt} in
    r ) REPLACE_R_WITH_I=true ;;
    h ) usage; exit 0 ;;
    \? ) echo "Invalid Option: -$OPTARG" 1>&2; usage; exit 1 ;;
  esac
done
shift $((OPTIND -1))

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    usage
    exit 1
fi

INFILE="$1"
OUTFILE="${2:--}"

if [ ! -f "$INFILE" ]; then
    echo "Input file not found: $INFILE" >&2
    exit 2
fi

awk_prog='BEGIN{ OFS="  " }
{
  checksum=$1
  sub(/^\S+\s+/, "")
  fname=$0
  gsub(/^ +| +$/, "", fname)
  if (match(fname, /(.*)_S[0-9]+_L[0-9]{3}_[RI][12]_001.fastq.gz$/, arr)) {
    sample=arr[1]
    newfname = "FASTQ/" sample "/" fname
    if (replace == "true") {
      gsub(/_R1_/, "_I1_", newfname)
      gsub(/_R2_/, "_I2_", newfname)
    }
    print checksum, newfname
  } else {
    print checksum, fname
  }
}'

if [ "$OUTFILE" = "-" ]; then
  # stream to stdout
  awk -v replace="$REPLACE_R_WITH_I" "$awk_prog" "$INFILE"
else
  awk -v replace="$REPLACE_R_WITH_I" "$awk_prog" "$INFILE" > "$OUTFILE"
fi
