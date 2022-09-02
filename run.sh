#!/bin/bash

function setup(){
  # sra-toolkit
  PATH="$HOME/bio-bin/sratoolkit.3.0.0-ubuntu64/bin/":$PATH
  hash prefetch > /dev/null || error "prefetch not found"
  hash fasterq-dump > /dev/null || error "fasterq-dump not found"

  # MFCompressC
  PATH="$HOME/bio-bin/MFCompress/":$PATH
  hash MFCompressC > /dev/null || error "MFCompressC not found"

  # prophasm
  PATH="$HOME/bio-bin/prophasm/":$PATH
  hash prophasm > /dev/null || error "prophasm not found"

  # metagraph
  PATH="$HOME/bio-bin/metagraph/metagraph/build/":$PATH
  hash metagraph > /dev/null || error "metagraph not found"

  # bcalm
  PATH="$HOME/bio-bin/bcalm/build/":$PATH
  hash bcalm > /dev/null || error "bcalm not found"

  # UST
  PATH="$HOME/bio-bin/UST/":$PATH
  hash ust > /dev/null || error "ust not found"

  # numbers of threads
  NTHREAD=4

  # k-mer size
  KMER_SIZES="31"

  # Accessions list(s)
  SEQUENCES=$(cat sequences-*.txt | grep -vE "#") # escape comments
  if [[ -z ${SEQUENCES} ]]; then
    error "there are no sequence!"
  fi

  # result file
  touch results.txt
  RESULTS=$(realpath results.txt)
}

function error(){
  echo "An error occurred: ${1:-"unknown"}"
  exit 1
}

compress_and_write_csv(){
  echo "** Compressing ${1}..."
  MFCompressC -t $NTHREAD -3 -o "${1}.mfc" "${1}"
  make_csv "${1}.mfc"
}

make_csv(){
  size=$(stat -c %s "${1}")
  printf "%s,%s\n" "${1}" "${size}" >> "$RESULTS"
}

launch_prophasm(){
  if [[ $2 == "--counts" ]];
  then
      error "can't use prophasm with count!"
  fi

  echo "*** Launching prophasm with accession ${1} and k=${K}"
  mkdir -p "prophasm"
  (
  cd prophasm
  infile="../${1}.fasta"
  outfile="${1}.pro.k${K}.fasta"
  outfile_stat="${1}.pro.k${K}.stat"
  prophasm -s ${outfile_stat} -k ${K} -i ${infile} -o ${outfile}
  make_csv ${outfile}
  compress_and_write_csv ${outfile}
  )
}

launch_metagraph(){
  echo "*** Launching metagraph with accession ${1} and k=${K}"
  mkdir -p metagraph
  (
  cd metagraph
  infile="../${1}.fasta"
  outfile="${1}.met.k${K}.fasta"
  metagraph
  )
}

launch_bcalm(){
  counts_param= # null string!
  counts_desc="" # empty string!
  if [[ $2 == "--counts" ]]; then
    counts_param="-all-abundance-counts"
    counts_desc="-counts"
  fi
  echo "*** Launching bcalm${counts_desc} with accession ${1} and k=${K}"
  mkdir -p bcalm
  (
  cd bcalm
  infile="../${1}.fasta"
  outfile_base="${1}.bca${counts_desc}.k${K}"
  outfile_extension=".unitigs.fa"

  bcalm -nb-cores ${NTHREAD} \
  -kmer-size $K ${counts_param}\
  -in ${infile} -out ${outfile_base}

  make_csv "${outfile_base}${outfile_extension}"
  compress_and_write_csv "${outfile_base}${outfile_extension}"
  )
}

launch_ust(){
  counts_param="0"
  counts_desc="" # empty string!
  if [[ $2 == "--counts" ]]; then
    counts_param="1"
    counts_desc="-counts"
  fi

  bcalm_file="bcalm/${1}.bca${counts_desc}.k${K}.unitigs.fa"
  if [[ ! -e "${bcalm_file}" ]]; then
    error "${bcalm_file} needed!"
  fi

  echo "*** Launching UST${counts_desc} with accession ${1} and k=${K}"
  mkdir -p ust
  (
  cd ust
  infile="../${bcalm_file}"
  outfile="${bcalm_file}.ust.fa"
  outfile_counts="${bcalm_file}.ust.counts" # TODO compress count file!
  ust -k "${K}" -i "${infile}" -a ${counts_param}
  make_csv "${outfile}"
  compress_and_write_csv "${outfile}"
  )
}

download_and_launch(){
  for S in ${SEQUENCES}; do
    echo -n "*** Downloading $S..."
    if [ -f "$S/$S.fasta" ]; then
      echo "Skipped."
    else
      echo
      prefetch --progress "$S"
    fi
    (
    cd "$S"
    echo -n "*** Converting $S to FASTA format..."
    if [ -f "$S.fasta" ]; then
      echo "Skipped."
    else
      echo
      # keep only biological reads which are present at least 2 times
      # produce a FASTA file
      fasterq-dump --threads $NTHREAD --progress \
      --skip-technical --split-3 --fasta \
      --outfile "${S}.fasta" "${S}.sra"

      # keep only the first file
      rm "${S}.fasta" "${S}_2.fasta"
      mv "${S}_1.fasta" "${S}.fasta"
    fi

    # iterating k-mer size
    for K in $KMER_SIZES; do
      #launch_prophasm "$S"
      #launch_bcalm "$S"
      #launch_bcalm "$S" "--counts"
      launch_ust "$S"
      launch_ust "$S" "--counts"
      #launch_metagraph "$S"
    done
    )
  done
}

setup
download_and_launch
exit 0