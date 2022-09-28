#!/bin/bash

setup(){
  DEBUG=true
  PLOT=false
  CLEAN=false

  # ---- add binaries to PATH ----
  # sra-toolkit
  PATH="$HOME/bio-bin/sratoolkit.3.0.0-ubuntu64/bin":$PATH
  hash prefetch > /dev/null || error "prefetch not found"
  hash fasterq-dump > /dev/null || error "fasterq-dump not found"

  # MFCompressC
  PATH="$HOME/bio-bin/MFCompress":$PATH
  hash MFCompressC > /dev/null || error "MFCompressC not found"

  # prophasm
  PATH="$HOME/bio-bin/prophasm":$PATH
  hash prophasm > /dev/null || error "prophasm not found"

  # metagraph
  PATH="$HOME/bio-bin/metagraph/metagraph/build":$PATH
  hash metagraph > /dev/null || error "metagraph not found"

  # bcalm
  PATH="$HOME/bio-bin/bcalm/build":$PATH
  hash bcalm > /dev/null || error "bcalm not found"

  # UST
  PATH="$HOME/bio-bin/UST":$PATH
  hash ust > /dev/null || error "ust not found"
  # --------------------------------

  # ---- define other global variables ----
  # root dir
  ROOT=$(realpath .)

  # numbers of threads
  NTHREAD=4

  # k-mer size (separated by spaces)
  KMER_SIZES=$(grep -vE "#" kmer-sizes.txt) # escape comments

  # Accessions list(s)
  SEQUENCES=$(grep -vE "#" sequences*.txt) # escape comments
  if [[ -z ${SEQUENCES} ]]; then
    error "there are no sequence!"
  fi

  # result file
  mkdir -p results
  RESULTS="results.csv"
  touch "$RESULTS"
  RESULTS=$(realpath "$RESULTS")

  # print headers
  printf "sequence,method,counts,kmer_size,file_type,compression,size\n" > "${RESULTS}"

  if $DEBUG; then
    echo "Debug mode activated"
    KMER_SIZES="15"
    SEQUENCES="SRR000001"
  fi;
  # --------------------------------
}

clean(){
  for S in $SEQUENCES; do
    (
    cd "$S" || exit
    # remove all except *.sra and *.fasta
    for f in *; do
      if [[ "$f" != *.fasta && "$f" != *.sra ]]; then
        rm -rf "$f"
      fi
    done
    )
  done
}

finalize(){
  # copy results on backup dir
  cp "$RESULTS" "${ROOT}/results/results_$(date +%F_%H.%M.%S).csv"
}

error(){
  echo "An error occurred: ${1:-"unknown"}"
  exit 1
}

compress_and_write_csv(){
  local uncompressed=$1
  local compressor=$2

  echo -n "** Compressing ${uncompressed} with "
  case ${compressor} in
  "mfcompress")
    echo "MFCompressC..."
    # check mimetype and first character
    if [[ $(mimetype -b "${uncompressed}") != "text/plain" || $(head -c 1 -z "${uncompressed}") != ">" ]]; then
      echo "Not a FASTA file. Ignored."
      return
    fi

    local compressed="${uncompressed}.mfc"
    [[ -f $compressed ]] || MFCompressC -t $NTHREAD -3 -o "${compressed}" "${uncompressed}"
    ;;
  "lzma")
    echo "lzma..."
    local compressed="${uncompressed}.lzma"
    [[ -f $compressed ]] || lzma --force --best --keep "${uncompressed}"
    ;;
  "gzip")
    echo "gzip..."
    local compressed="${uncompressed}.gz"
    [[ -f $compressed ]] || gzip --force --keep --best "${uncompressed}"
    ;;
  *)
    error "compressor not found"
  esac

  write_to_csv "${compressed}" "${compressor:-"gzip"}" "${uncompressed##*.}"
}

write_to_csv(){
  # find file size
  # shellcheck disable=SC2155
  local size=$(stat -c %s "${1}")
  # compression type if given otherwise none
  local compression=${2:-"none"}
  # filetype if given otherwise found it
  local filetype=${3:-${1##*.}}
  # recall csv headers: sequence,method,counts,kmer-size,file-type,compression,size
  printf "%s,%s,%s,%s,%s,%s,%s\n" "$S" "$method" "$counts" "$K" "$filetype" "$compression" "$size" >> "$RESULTS"
}

compress_all_and_write_csv(){
  tools="mfcompress lzma gzip"
  for t in $tools; do
    compress_and_write_csv "${1}" "${t}"
  done
}

write_starting_size(){
  local method="none"
  local counts="none"
  local K="0"
  local filetype="fasta"
  local compression="none"
  write_to_csv "${S}.fasta"
}

launch_prophasm(){
  local method="prophasm";
  local root_dir=$method
  local counts="no-counts"
  local infile="../${S}.fasta"
  local outfile="${S}.pro.k${K}.fasta"
  local outfile_stat="${S}.pro.k${K}.stat"

  echo "*** Launching prophasm ($counts) with accession ${S} and k=${K}"
  mkdir -p $root_dir
  (
  cd $root_dir || error "cannot cd $root_dir"
  [[ -f $outfile ]] || prophasm -s "${outfile_stat}" -k "${K}" -i "${infile}" -o "${outfile}"
  write_to_csv "${outfile}"
  compress_all_and_write_csv "${outfile}"
  )
}

launch_bcalm(){
  local method="bcalm";
  local root_dir=$method
  if [[ $1 == "--counts" ]]; then
    local counts="counts"
    local counts_param="-all-abundance-counts"
  else
    local counts="no-counts"
    local counts_param= # null string!
  fi
  local infile="../${S}.fasta"
  local outfile_base="${S}.${method}-${counts}.k${K}"
  local outfile="${outfile_base}.unitigs.fa"
  local outfile_better_name="${outfile_base}.fasta"

  # adjust parameters if kmer-size is too low
  if [[ $K -le 10 ]]; then
    local bca_min_size=5
  else
    local bca_min_size=10
  fi

  echo "*** Launching bcalm (${counts}) with accession ${S} and k=${K}"
  mkdir -p $root_dir
  (
  cd $root_dir || error "cannot cd to bcalm"

  if [[ ! -f $outfile_better_name ]]; then
    # shellcheck disable=SC2086
    bcalm -nb-cores ${NTHREAD} \
    -kmer-size "${K}" -abundance-min 1 -minimizer-size $bca_min_size ${counts_param}\
    -in "${infile}" -out "${outfile_base}"
    # give better filename
    mv "$outfile" "$outfile_better_name"
    local outfile=$outfile_better_name
  fi

  write_to_csv "${outfile}"
  compress_all_and_write_csv "${outfile}"
  )
}

launch_ust(){
  local method="ust";
  local root_dir=$method
  if [[ $1 == "--counts" ]]; then
    local counts_param="1"
    local counts="counts"
  else
    local counts_param="0"
    local counts="no-counts"
  fi
  local bcalm_file="${S}.bcalm-${counts}.k${K}.fasta"
  local infile="../bcalm/${bcalm_file}"
  local outfile="${bcalm_file}.ust.fa"
  local outfile_counts="${bcalm_file}.ust.counts"
  # better names
  local outfile_better_name="${S}.${method}-${counts}.k${K}.fasta"
  local outfile_better_name_counts="${S}.${method}-${counts}.k${K}.counts"

  echo "*** Launching UST (${counts}) with accession ${S} and k=${K}"
  mkdir -p $root_dir
  (
  cd $root_dir || error "cannot cd to ust"

  # if there is no old file
  if [[ ! -f $outfile_better_name ]]; then
    ust -k "${K}" -i "${infile}" -a ${counts_param}

    # give better names to output files
    mv "$outfile" "$outfile_better_name"
    [[ ${counts_param} == "0" ]] || mv "$outfile_counts" "$outfile_better_name_counts"
  fi

  local outfile=$outfile_better_name
  local outfile_counts=$outfile_better_name_counts

  write_to_csv "${outfile}"
  [[ ${counts_param} == "0" ]] || write_to_csv "${outfile_counts}"
  compress_all_and_write_csv "${outfile}"
  [[ ${counts_param} == "0" ]] || compress_all_and_write_csv "${outfile_counts}"
  )
}

launch_metagraph_count(){
  local method="metagraph";
  local root_dir="${method}-k${K}"
  local infile="../${S}.fasta"
  local outfile_base="${S}.${method}.k${K}"
  local graph="${outfile_base}.dbg"
  local annotations_binary="${outfile_base}.column.annodbg"
  local annotations_counts="${outfile_base}.column.annodbg.counts"

  if [[ $1 == "--counts" ]]; then
    local counts="counts"
  else
    local counts="no-counts"
  fi

  echo "*** Launching ${method} (${counts}) with accession ${S} and k=${K}"

  mkdir -p "$root_dir"
  (
  cd "$root_dir" || error "can't cd to metagraph"
  mkdir -p swap

  # use count_dbg commands (with coordinates, not counts!)
  echo -n "Building DBG..."
  if [[ ! -f "$graph" ]]; then
    metagraph build --parallel "${NTHREAD}" --kmer-length "${K}" --count-kmers -o "${graph}" "${infile}"
    echo "Done."
  else
    echo "Skipped."
  fi

  # if need counts and there are no old files
  if [[ $counts == "counts" && ! -f $annotations_binary && ! -f $annotations_counts ]]; then
    echo "* Dumping contigs to fasta"
    metagraph clean -p ${NTHREAD} --to-fasta -o "${outfile_base}.contigs" "${graph}"

    echo "* Building $annotations_binary and $annotations_counts"
    metagraph annotate --anno-filename --count-kmers\
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${outfile_base}.contigs.fasta.gz"

    echo "* Transforming coordinates into row_diff (3 stages)"
    echo "- stage: 0"
    metagraph transform_anno --count-kmers --anno-type row_diff --row-diff-stage 0 --disk-swap swap\
    --parallel "${NTHREAD}" --infile-base "${graph}" -o "${outfile_base}.row_count" "${annotations_binary}"
    echo "- stage: 1"
    metagraph transform_anno --count-kmers --anno-type row_diff --row-diff-stage 1 --disk-swap swap\
    --parallel "${NTHREAD}" --infile-base "${graph}" -o "${outfile_base}.row_reduction" "${annotations_binary}"
    echo "- stage: 2"
    metagraph transform_anno --count-kmers --anno-type row_diff --row-diff-stage 2 --disk-swap swap\
    --parallel "${NTHREAD}" --infile-base "${graph}" -o "${outfile_base}.row_diff" "${annotations_binary}"

    echo "* Transforming coordinates into row_diff_int_brwt"
    metagraph transform_anno --anno-type row_diff_int_brwt --greedy --fast\
    --parallel "${NTHREAD}" --infile-base "${graph}" -o "${outfile_base}" "${annotations_binary}"

    echo "* Relaxing annotations"
    metagraph relax_brwt \
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${outfile_base}.row_diff_int_brwt.annodbg"

    # clean temp files
    rm -rf swap *.pred *.pred_boundary *.rd_succ *.succ *.succ_boundary *.linkage *.row_count *.row_reduction *.anchors
  fi

  write_to_csv "${graph}"
  compress_all_and_write_csv "${graph}"

  if [[ $counts == "counts" ]]; then
    write_to_csv "${annotations_binary}"
    write_to_csv "${annotations_counts}"
    compress_all_and_write_csv "${annotations_binary}"
    compress_all_and_write_csv "${annotations_counts}"
  fi
  )
}

# launch metagraph_count before this!
launch_metagraph_assemble(){
  local method="assembly";
  local root_dir=$method
  local metagraph_contigs="../metagraph-k${K}/${S}.metagraph.k${K}.contigs.fasta.gz"
  local metagraph_counts="../metagraph-k${K}/${S}.metagraph.k${K}.contigs.kmer_counts.gz"
  local outfile_base="${S}.assembled.k${K}"
  local outfile_contigs="${outfile_base}.fasta"
  local outfile_counts="${outfile_base}.kmer_counts"

  if [[ $1 == "--counts" ]]; then
    local counts="counts"
  else
    local counts="no-counts"
  fi

  echo "*** Launching metagraph assemble (${counts}) with accession ${1} and k=${K}"

  mkdir -p $root_dir
  (
  cd $root_dir || error "cannot cd to metagraph"

  if [[ ! -f "$metagraph_contigs" || ! -f "$metagraph_counts" ]]; then
    error "${metagraph_contigs} or ${metagraph_counts} not found!"
  fi

  if [[ ! -f $outfile_contigs ]]; then
    # not --best compression, decompress!
    zless "${metagraph_contigs}" > "${outfile_contigs}"
  fi

  if [[ $counts == "counts" && ! -f $outfile_counts ]]; then
     zless "${metagraph_counts}" > "${outfile_counts}"
  fi

  write_to_csv "${outfile_contigs}"
  compress_all_and_write_csv "${outfile_contigs}"

  if [[ $counts == "counts" ]]; then
    write_to_csv "${outfile_counts}"
    compress_all_and_write_csv "${outfile_counts}"
  fi
  )
}

download_and_launch(){
  mkdir -p sequences
  cd sequences || error

  for S in ${SEQUENCES}; do
    echo -n "*** Downloading $S..."
    if [[ -f "$S/$S.sra" || -f "$S/$S.fasta" ]]; then
      echo "Skipped."
    else
      echo
      # download or skip sequence if errors!
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
      --skip-technical --fasta-unsorted \
      --outfile "${S}.fasta" "${S}.sra"
    fi

    write_starting_size

    # iterating k-mer size
    for K in $KMER_SIZES; do
      launch_metagraph_count
      launch_metagraph_count --counts
      launch_metagraph_assemble
      launch_metagraph_assemble --counts
      launch_prophasm
      launch_bcalm
      launch_bcalm --counts
      launch_ust
      launch_ust --counts
    done
    )
  done
}

setup
if [ $CLEAN ]; then clean; fi
download_and_launch
finalize

if $PLOT; then
  # plot
  source venv/bin/activate
  ./plot.py "$RESULTS"
fi

if ! $DEBUG; then
  # power off the computer
  shutdown -h +3
fi

exit 0