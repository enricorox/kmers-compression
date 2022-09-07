#!/bin/bash

setup(){
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

  # k-mer size (separated by spaces)
  KMER_SIZES=$(grep -vE "#" kmer-sizes.txt)

  # Accessions list(s)
  SEQUENCES=$(grep -vE "#" sequences*.txt) # escape comments
  if [[ -z ${SEQUENCES} ]]; then
    error "there are no sequence!"
  fi

  # result file
  mkdir -p results
  RESULTS="results/results_$(date +%F_%H.%M.%S).csv"
  touch "$RESULTS"
  RESULTS=$(realpath "$RESULTS")

  # print headers
  printf "sequence,method,counts,kmer_size,file_type,compression,size\n" >> "${RESULTS}"
}

function error(){
  echo "An error occurred: ${1:-"unknown"}"
  exit 1
}

compress_and_write_csv(){
  # $1: file to compress
  # $2: compression method
  echo -n "** Compressing ${1} with "
  case ${2} in
  "mfcompress")
    echo "MFCompressC..."
    if [[ $(mimetype -b "${1}") != "text/plain" || $(head -c 1 -z "${1}") != ">" ]]; then
      echo "Not a FASTA file. Ignored."
      return
    fi

    outfile="${1}.mfc"
    [[ -f $outfile ]] || MFCompressC -t $NTHREAD -3 -o "${outfile}" "${1}"
    ;;
  "lzma")
    echo "lzma..."
    outfile="${1}.lzma"
    [[ -f $outfile ]] || lzma --force --best --keep "${1}"
    ;;
  *)
    echo "gzip..."
    outfile="${1}.gz"
    [[ -f $outfile ]] || gzip --force --keep --best "${1}"
    ;;
  esac

  write_to_csv "${outfile}" "${2:-"gzip"}" "${1##*.}"
}

write_to_csv(){
  local size=$(stat -c %s "${1}")
  #printf "%s,%s\n" "${1}" "${size}" >> "$RESULTS"
  compression=${2:-"none"}
  filetype=${3:-${1##*.}}
  # headers: sequence,method,counts,kmer-size,file-type,compression,size
  printf "%s,%s,%s,%s,%s,%s,%s\n" "$S" "$method" "$counts" "$K" "$filetype" "$compression" "$size" >> "$RESULTS"
}

compress_all_and_write_csv(){
  tools="mfcompress lzma gzip"
  for t in $tools; do
    compress_and_write_csv "${1}" "${t}"
  done
}

launch_prophasm(){
  method="prophasm"; counts="no-counts"
  if [[ $2 == "--counts" ]];
  then
      error "can't use prophasm with count!"
  fi

  echo "*** Launching prophasm with accession ${1} and k=${K}"
  mkdir -p "prophasm"
  (
  cd prophasm || error "cannot cd prophasm"
  infile="../${1}.fasta"
  outfile="${1}.pro.k${K}.fasta"
  outfile_stat="${1}.pro.k${K}.stat"
  [[ -f $outfile ]] || prophasm -s "${outfile_stat}" -k "${K}" -i "${infile}" -o "${outfile}"
  write_to_csv "${outfile}"
  compress_all_and_write_csv "${outfile}"
  )
}

launch_metagraph(){
  method="metagraph"; counts="no-counts"
  counts_param= # null string!
  counts_desc="" # empty string!
  if [[ $2 == "--counts" ]]; then
    counts_desc="-counts"
    counts="counts"
  fi

  echo "*** Launching metagraph${counts_desc} with accession ${1} and k=${K}"
  mkdir -p metagraph
  (
  cd metagraph || error "can't cd to metagraph"
  infile="../${1}.fasta"
  outfile_base="${S}.met.k${K}"
  graph="${outfile_base}.dbg"
  annotations_binary="${outfile_base}.column.annodbg"
  annotations_coord="${outfile_base}.column.annodbg.coords"

  # use count_dbg commands (with coordinates, not counts!)
  if [[ ! -f "$graph" ]]; then
    metagraph build --parallel "${NTHREAD}" --kmer-length "${K}" -o "${graph}" "${infile}"
  else
    echo "Skipped."
  fi

  if [[ $2 == "--counts" && ! -f $annotations_binary && ! -f $annotations_coord ]]; then
    # produce .column.annodbg (binary graph annotation) and .column.annodbg.coords (kmers coordinates) files
    metagraph annotate --anno-filename --coordinates\
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${infile}"

    # transform coordinates into column_coord
    metagraph transform_anno --anno-type column_coord --coordinates \
    --parallel "${NTHREAD}" --outfile-base "${outfile_base}" "${annotations_binary}"

    # transform coordinates into row_diff (3 stages)
    for i in {0,1,2}; do
      metagraph transform_anno --anno-type row_diff --coordinates --row-diff-stage "$i"\
      --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${annotations_binary}"
    done

    # transform coordinates into row_diff_coord
    metagraph transform_anno --anno-type row_diff_coord\
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${annotations_binary}"

    # transform coordinates into row_diff_brwt_coord
    metagraph transform_anno --anno-type row_diff_brwt_coord --greedy --fast --subsample 1000000\
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${annotations_binary}"
  fi

  # clean temp files
  echo "Press ENTER..."
  #read
  for file in "${outfile_base}"*; do
    [[ $file == "$graph" || $file == "$annotations_binary" || $file == "$annotations_coord" ]] ||\
    [[ $file == *.lzma || $file == *.gz || $file == *.mfc ]] || rm "$file"
  done

  write_to_csv "${graph}"
  compress_all_and_write_csv "${graph}"

  if [[ $2 == "--counts" ]]; then
    write_to_csv "${annotations_binary}"
    write_to_csv "${annotations_coord}"
    compress_all_and_write_csv "${annotations_binary}"
    compress_all_and_write_csv "${annotations_coord}"
  fi
  )
}

launch_bcalm(){
  method="bcalm"; counts="no-counts"
  counts_param= # null string!
  counts_desc="" # empty string!
  if [[ $2 == "--counts" ]]; then
    counts_param="-all-abundance-counts"
    counts_desc="-counts"
    counts="counts"
  fi
  echo "*** Launching bcalm${counts_desc} with accession ${1} and k=${K}"
  mkdir -p bcalm
  (
  cd bcalm || error "cannot cd to bcalm"
  infile="../${1}.fasta"
  outfile_base="${1}.bca${counts_desc}.k${K}"
  outfile_extension=".unitigs.fa"
  outfile="${outfile_base}${outfile_extension}"
  outfile1="${outfile_base}.fasta"
  if [[ ! -f $outfile1 ]]; then
    bcalm -nb-cores ${NTHREAD} \
    -kmer-size "${K}" ${counts_param}\
    -in "${infile}" -out "${outfile_base}"
    # give better filename
    mv "$outfile" "$outfile1"
  fi
  outfile=$outfile1
  write_to_csv "${outfile}"
  compress_all_and_write_csv "${outfile}"
  )
}

launch_ust(){
  method="ust"; counts="no-counts"
  counts_param="0"
  counts_desc="" # empty string!
  if [[ $2 == "--counts" ]]; then
    counts_param="1"
    counts_desc="-counts"
    counts="counts"
  fi

  bcalm_file="${1}.bca${counts_desc}.k${K}.fasta"

  echo "*** Launching UST${counts_desc} with accession ${1} and k=${K}"
  mkdir -p ust
  (
  cd ust || error "cannot cd to ust"
  infile="../bcalm/${bcalm_file}"
  outfile="${bcalm_file}.ust.fa"
  outfile_counts="${bcalm_file}.ust.counts"
  outfile1="${1}.ust${counts_desc}.k${K}.fasta" # better name
  outfile1_counts="${1}.ust-counts.k${K}.counts"
  [[ -f $outfile1 ]] || ust -k "${K}" -i "${infile}" -a ${counts_param}

  # give better names to output files
  mv "$outfile" "$outfile1"
  [[ ${counts_param} == "0" ]] || mv "$outfile_counts" "$outfile1_counts"
  outfile=$outfile1
  outfile_counts=$outfile1_counts

  write_to_csv "${outfile}"
  [[ ${counts_param} == "0" ]] || write_to_csv "${outfile_counts}"
  compress_all_and_write_csv "${outfile}"
  [[ ${counts_param} == "0" ]] || compress_all_and_write_csv "${outfile_counts}"
  )
}

launch_metagraph_assemble(){
  method="assembly"; counts="no-counts"
  if [[ $2 == "--counts" ]]; then
   error "can't use counts!"
  fi

  metagraph_file="${1}.met.k${K}.dbg"

  echo "*** Launching metagraph assemble with accession ${1} and k=${K}"
  mkdir -p metagraph
  (
  cd metagraph || error "cannot cd to metagraph"
  infile="${metagraph_file}"
  outfile_base="${1}.assembled"
  outfile="${outfile_base}.fasta"
  if [[ ! -f $outfile ]]; then
    metagraph assemble -p ${NTHREAD} -o "$outfile_base" "$infile"
    zless "${outfile_base}.fasta.gz" > "${outfile}"
    # not --best
    rm "${outfile_base}.fasta.gz"
  fi

  write_to_csv "${outfile}"
  compress_all_and_write_csv "${outfile}"
  )
}

download_and_launch(){
  for S in ${SEQUENCES}; do
    echo -n "*** Downloading $S..."
    if [ -f "$S/$S.sra" ]; then
      echo "Skipped."
    else
      echo
      prefetch --progress "$S" || continue
    fi
    (
    cd "$S" || error "can't cd to ${S}"
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
      rm -f "${S}.fasta" "${S}_2.fasta"
      mv -f "${S}_1.fasta" "${S}.fasta"
    fi

    # iterating k-mer size
    for K in $KMER_SIZES; do
      (
      launch_prophasm "$S"
      launch_bcalm "$S"
      launch_bcalm "$S" "--counts"
      launch_ust "$S"
      launch_ust "$S" "--counts"
      # produce the same results as with counts
      launch_metagraph "$S"
      launch_metagraph "$S" "--counts"
      launch_metagraph_assemble "$S"
      )
    done
    )
  done
}

setup
download_and_launch
./plot.py "$RESULTS"
exit 0