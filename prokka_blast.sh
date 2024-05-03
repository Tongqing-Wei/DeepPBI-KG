#!/bin/bash

# use getopts Parse command line arguments
while getopts ":i:o:p:d:b:O:" opt; do
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
    d)
      db="$OPTARG"
      ;;
    b)
      blast="$OPTARG"
      ;;
    O)
      out="$OPTARG"
      ;;
    \?)
      echo "invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "option -$OPTARG need a parameter." >&2
      exit 1
      ;;
  esac
done


#prokka batch annotation
for file in `ls $input_dir/*.fasta`; do
      val=$(echo "${file##*/}" | cut -d '.' -f 1)
      "$prokka" "$file" --outdir "$output_dir/prokka_$val" --prefix $val --force
      "$blast" -query "$file" -db "$db" -outfmt 6 -max_target_seqs 1 -out "$out/$val.out" 
done

