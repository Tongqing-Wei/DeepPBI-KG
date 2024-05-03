#!/bin/bash

# Use getopts to parse command line arguments
while getopts ":i:o:p:" opt; do
  case $opt in
    i)
      input_dir="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      ;;
    p)
      prokka="$OPTARG"
      ;;
    \?)
      echo "invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "option -$OPTARG need a patameter." >&2
      exit 1
      ;;
  esac
done


#prokka batch annotation
for file in `ls $input_dir/*.fasta`; do
      val=$(echo "${file##*/}" | cut -d '.' -f 1)
      "$prokka" "$file" --outdir "$output_dir/prokka_$val" --prefix $val 
done

