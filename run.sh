#!/bin/bash

function setup(){
  # sra-toolkit
  PATH="$HOME/bio-bin/sratoolkit.3.0.0-ubuntu64/bin/":$PATH
  # MFCompressC
  PATH="$HOME/bio-bin/MFCompress/":$PATH
  # prophasm
  PATH="$HOME/bio-bin/prophasm/":$PATH
  # metagraph
  PATH="$HOME/bio-bin/metagraph/metagraph/build/":$PATH
  # bcalm
  PATH="$HOME/bio-bin/bcalm/build/":$PATH
  # UST
  PATH="$HOME/bio-bin/UST/":$PATH

  # numbers of threads
  NTHREAD=4

  # k-mer size
  K=31

  # Accessions list(s)
  SEQUENCES=$(cat sequences-*.txt | grep -vE "#") # escape comments
  if [[ -z ${SEQUENCES} ]]; then
    error "there are no sequence!"
  fi
  touch results.txt
  RESULTS=$(realpath results.txt)
}

function error(){
  echo "An error occurred: ${1:="unknown"}"
  exit 1
}

compress_and_make_csv(){
  MFCompressC -t $NTHREAD -3 -o "${1}.mfc" "${1}"
  make_csv "${1}.mfc"
}

make_csv(){
  size=$(stat -c %s "${1}")
  printf "%s,%s\n" "${1}" "${size}" >> "$RESULTS"
}

launch_prophasm(){
  echo "*** Launching prophasm with accession ${1} and k=${K}"
  mkdir -p "prophasm"
  (
  cd prophasm
  infile="../${1}.fasta"
  outfile="${1}.pro.k${K}.fasta"
  outfile_stat="${1}.pro.k${K}.stat"
  prophasm -s ${outfile_stat} -k ${K} -i ${infile} -o ${outfile}
  make_csv ${outfile}
  compress_and_make_csv ${outfile}
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

launch_bacalm(){
  echo "*** Launching UST with accession ${1} and k=${K}"
  mkdir -p bcalm
  (
  cd ust
  infile="../${1}.fasta"
  bcalm
  )
}

launch_ust(){
  echo "*** Launching UST with accession ${1} and k=${K}"
  mkdir -p ust
  (
  cd ust
  infile="../${1}.fasta"
  ust
  )
}

download_and_launch(){
  for S in ${SEQUENCES};
  do
    echo "Downloading $S..."
    prefetch --progress "$S"
    (
    cd "$S"
    echo "Converting $S to FASTA format..."

    # remove all fasta files to avoid fasterq-dump errors
    rm -f ./*.fasta

    # keep only biological reads which are present at least 2 times
    # produce a FASTA file
    fasterq-dump --threads $NTHREAD --progress \
    --skip-technical --split-3 --fasta \
    --outfile "${S}.fasta" "${S}.sra"

    # keep only the first file
    rm "${S}.fasta" "${S}_2.fasta"
    mv "${S}_1.fasta" "${S}.fasta"

    launch_prophasm "$S"
    #launch_bcalm "$S"
    #launch_ust "$S"
    #launch_metagraph "$S"
    )
  done
}

setup
download_and_launch
exit 0