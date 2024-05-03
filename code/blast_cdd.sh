#!/bin/bash

# Use getopts to parse command line arguments
while getopts ":i:o:b:d:" opt; do
  case $opt in
    i)
      input_dir="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      ;;
    b)
      blastp="$OPTARG"
      ;;
    d)
      database="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "option -$OPTARG need a parameter." >&2
      exit 1
      ;;
  esac
done


#blastp batch comparison
for file in `ls $input_dir/*.fasta`; do
      val=$(echo "${file##*/}" | cut -d '.' -f 1)
      "$blastp" -db "$database" -query "$file" -max_target_seqs 3 -evalue 1e-5 -out "$output_dir/$val.blast" 
done

