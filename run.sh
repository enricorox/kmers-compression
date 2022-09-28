#!/bin/bash

setup(){
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

  # ---- define global variables ----
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

  # tex table file
  RESULTS_TEX="results/results_$(date +%F_%H.%M.%S).tex"
  touch "$RESULTS_TEX"
  RESULTS_TEX=$(realpath "$RESULTS_TEX")
  printf  "%s\n"\
          "\\begin{center}"\
          "\\begin{tabular}{ |c|c|c|c|c|c|c| }"\
          "\\hline"\
          "sequence & method & counts & kmer-size & file-type & compression & size\\\\\\"\
          "\\hline" >> "${RESULTS_TEX}"
  # --------------------------------
}

clean(){
  for S in $SEQUENCES; do
    (
    cd "$S" || exit
    rm -rf $(ls | grep -vE "${S}.fasta")
    )
  done
}

finalize(){
  # finalize latex table
  printf  "%s\n"\
          "\\hline"\
          "\\end{tabular}"\
          "\\end{center}" >> "$RESULTS_TEX"
  # copy results on root dir
  cp "$RESULTS" results.csv
}

error(){
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
  compression=${2:-"none"}
  filetype=${3:-${1##*.}}
  # headers: sequence,method,counts,kmer-size,file-type,compression,size
  printf "%s,%s,%s,%s,%s,%s,%s\n" "$S" "$method" "$counts" "$K" "$filetype" "$compression" "$size" >> "$RESULTS"
  printf "%s & %s & %s & %s & %s & %s & %s\\\\\\\\\n" "$S" "$method" "$counts" "$K" "$filetype" "$compression" "$size" >> "$RESULTS_TEX"
}

compress_all_and_write_csv(){
  tools="mfcompress lzma gzip"
  for t in $tools; do
    compress_and_write_csv "${1}" "${t}"
  done
}

launch_prophasm(){
  local method="prophasm";
  local counts="no-counts"
  if [[ $2 == "--counts" ]];
  then
      error "can't use prophasm with count!"
  fi

  echo "*** Launching prophasm with accession ${1} and k=${K}"
  mkdir -p $method
  (
  cd $method || error "cannot cd $method"
  local infile="../${1}.fasta"
  local outfile="${1}.pro.k${K}.fasta"
  local outfile_stat="${1}.pro.k${K}.stat"
  [[ -f $outfile ]] || prophasm -s "${outfile_stat}" -k "${K}" -i "${infile}" -o "${outfile}"
  write_to_csv "${outfile}"
  compress_all_and_write_csv "${outfile}"
  )
}

launch_bcalm(){
  local method="bcalm";

  if [[ $2 == "--counts" ]]; then
    local counts="counts"
    local counts_param="-all-abundance-counts"
    local counts_desc="-counts"
  else
    local counts="no-counts"
    local counts_param= # null string!
    local counts_desc="" # empty string!
  fi

  # adjust parameters if kmer-size is too low
  if [[ $K -le 10 ]]; then
    bca_min_size=5
  else
    bca_min_size=10
  fi

  echo "*** Launching bcalm${counts_desc} with accession ${1} and k=${K}"
  mkdir -p $method
  (
  cd $method || error "cannot cd to bcalm"
  local infile="../${1}.fasta"
  local outfile_base="${1}.bca${counts_desc}.k${K}"
  local outfile_extension=".unitigs.fa"
  local outfile="${outfile_base}${outfile_extension}"
  local outfile_better_name="${outfile_base}.fasta"

  if [[ ! -f $outfile_better_name ]]; then
    bcalm -nb-cores ${NTHREAD} \
    -kmer-size "${K}" -abundance-min 1 -minimizer-size $bca_min_size ${counts_param}\
    -in "${infile}" -out "${outfile_base}"
    # give better filename
    mv "$outfile" "$outfile_better_name"
  fi

  local outfile=$outfile_better_name
  write_to_csv "${outfile}"
  compress_all_and_write_csv "${outfile}"
  )
}

launch_ust(){
  local method="ust";

  if [[ $2 == "--counts" ]]; then
    counts_param="1"
    counts_desc="-counts"
    counts="counts"
  else
    local counts_param="0"
    local counts_desc="" # empty string!
    local counts="no-counts"
  fi

  local bcalm_file="${1}.bca${counts_desc}.k${K}.fasta"

  echo "*** Launching UST${counts_desc} with accession ${1} and k=${K}"
  mkdir -p ust
  (
  cd ust || error "cannot cd to ust"
  local infile="../bcalm/${bcalm_file}"
  local outfile="${bcalm_file}.ust.fa"
  local outfile_counts="${bcalm_file}.ust.counts"

  # better names
  local outfile_better_name="${1}.ust${counts_desc}.k${K}.fasta"
  local outfile_better_name_counts="${1}.ust-counts.k${K}.counts"
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

launch_metagraph_coord1(){
  local method="metagraph";
  local counts="no-counts"

  if [[ $2 == "--counts" ]]; then
    local counts_desc="-counts"
    local counts="counts"
  else
    local counts_desc="" # empty string!
    local counts_param= # null string!
  fi

  echo "*** Launching ${method}${counts_desc} with accession ${1} and k=${K}"
  local root_dir="${method}-k${K}"
  mkdir -p "$root_dir"
  (
  cd "$root_dir" || error "can't cd to metagraph"
  mkdir -p swap
  local infile="../${1}.fasta"
  local outfile_base="${S}.met.k${K}"
  local graph="${outfile_base}.dbg"
  local annotations_binary="${outfile_base}.column.annodbg"
  local annotations_coord="${outfile_base}.column.annodbg.coords"

  # use count_dbg commands (with coordinates, not counts!)
  echo -n "Building DBG..."
  if [[ ! -f "$graph" ]]; then
    metagraph build --parallel "${NTHREAD}" --kmer-length "${K}" -o "${graph}" "${infile}"
    echo "Done."
  else
    echo "Skipped."
  fi

  if [[ $2 == "--counts" && ! -f $annotations_binary && ! -f $annotations_coord ]]; then
    # produce .column.annodbg (binary graph annotation) and .column.annodbg.coords (kmers coordinates) files
    echo "* Building $annotations_binary and $annotations_coord"
    metagraph annotate --anno-filename --coordinates\
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${infile}"

    # transform coordinates into column_coord
    echo "* Transforming coordinates into column_coord"
    metagraph transform_anno --anno-type column_coord --coordinates \
    --parallel "${NTHREAD}" -o "${outfile_base}" "${annotations_binary}"

    # transform coordinates into row_diff (3 stages)
    echo "* Transforming coordinates into row_diff (3 stages)"
    echo "- stage: 0"
    metagraph transform_anno --anno-type row_diff --coordinates --row-diff-stage 0 --disk-swap swap\
    --parallel "${NTHREAD}" --infile-base "${graph}" -o "${outfile_base}.row_count" "${annotations_binary}"
    echo "- stage: 1"
    metagraph transform_anno --anno-type row_diff --coordinates --row-diff-stage 1 --disk-swap swap\
    --parallel "${NTHREAD}" --infile-base "${graph}" -o "${outfile_base}.row_reduction" "${annotations_binary}"
    echo "- stage: 2"
    metagraph transform_anno --anno-type row_diff --coordinates --row-diff-stage 2 --disk-swap swap\
    --parallel "${NTHREAD}" --infile-base "${graph}" -o "${outfile_base}.row_diff" "${annotations_binary}"

    # transform coordinates into row_diff_coord
    echo "* Transforming coordinates into row_diff_coord"
    metagraph transform_anno --anno-type row_diff_coord\
    --parallel "${NTHREAD}" --infile-base "${graph}" -o "${outfile_base}" "${annotations_binary}"

    # transform coordinates into row_diff_brwt_coord
    echo "* Producing brwt linkage"
    metagraph transform_anno --anno-type row_diff_brwt_coord --greedy --fast --linkage\
    --parallel "${NTHREAD}" --infile-base "${graph}" -o "${outfile_base}.linkage" "${annotations_binary}"

    echo "* Transforming coordinates into row_diff_brwt_coord"
    metagraph transform_anno --anno-type row_diff_brwt_coord --linkage-file "${outfile_base}.linkage"\
    --parallel "${NTHREAD}" --infile-base "${graph}" -o "${outfile_base}" "${annotations_binary}"

    # relax annotations
    echo "* Relaxing annotations"
    metagraph relax_brwt --relax-arity 32\
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${outfile_base}.row_diff_brwt_coord.annodbg"

    # clean temp files
    rm -rf swap *.pred *.pred_boundary *.rd_succ *.succ *.succ_boundary *.linkage *.row_count *.row_reduction *.anchors
  fi

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

launch_metagraph_count(){
  local method="metagraph2";
  local counts="no-counts"

  if [[ $2 == "--counts" ]]; then
    local counts_desc="-counts"
    local counts="counts"
  else
    local counts_desc="" # empty string!
    local counts_param= # null string!
  fi

  echo "*** Launching ${method}${counts_desc} with accession ${1} and k=${K}"
  local root_dir="${method}-k${K}"
  mkdir -p "$root_dir"
  (
  cd "$root_dir" || error "can't cd to metagraph"
  mkdir -p swap
  local infile="../${1}.fasta"
  local outfile_base="${S}.met.k${K}"
  local graph="${outfile_base}.dbg"
  local annotations_binary="${outfile_base}.column.annodbg"
  local annotations_counts="${outfile_base}.column.annodbg.counts"

  # use count_dbg commands (with coordinates, not counts!)
  echo -n "Building DBG..."
  if [[ ! -f "$graph" ]]; then
    metagraph build --parallel "${NTHREAD}" --kmer-length "${K}" --count-kmers -o "${graph}" "${infile}"
    echo "Done."
  else
    echo "Skipped."
  fi

  if [[ $2 == "--counts" && ! -f $annotations_binary && ! -f $annotations_counts ]]; then
    # dumping contigs to fasta
    echo "* Dumping contigs to fasta"
    metagraph clean -p ${NTHREAD} --to-fasta -o "${outfile_base}.contigs" "${graph}"

    # produce .column.annodbg (binary graph annotation) and .column.annodbg.coords (kmers coordinates) files
    echo "* Building $annotations_binary and $annotations_counts"
    metagraph annotate --anno-filename --count-kmers\
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${outfile_base}.contigs.fasta.gz"

    # transform coordinates into row_diff (3 stages)
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

    # transform annotations into row_diff_int_brwt
    echo "* Transforming coordinates into row_diff_int_brwt"
    metagraph transform_anno --anno-type row_diff_int_brwt --greedy --fast\
    --parallel "${NTHREAD}" --infile-base "${graph}" -o "${outfile_base}" "${annotations_binary}"

    # relax annotations
    echo "* Relaxing annotations"
    metagraph relax_brwt \
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${outfile_base}.row_diff_int_brwt.annodbg"

    # clean temp files
    rm -rf swap *.pred *.pred_boundary *.rd_succ *.succ *.succ_boundary *.linkage *.row_count *.row_reduction *.anchors
  fi

  write_to_csv "${graph}"
  compress_all_and_write_csv "${graph}"

  if [[ $2 == "--counts" ]]; then
    write_to_csv "${annotations_binary}"
    write_to_csv "${annotations_counts}"
    compress_all_and_write_csv "${annotations_binary}"
    compress_all_and_write_csv "${annotations_counts}"
  fi
  )
}

# don't use this!
launch_metagraph_coord2(){
  method="metagraph_coord2"; counts="no-counts"
  counts_param= # null string!
  counts_desc="" # empty string!
  if [[ $2 == "--counts" ]]; then
    counts_desc="-counts"
    counts="counts"
  fi

  echo "*** Launching metagraph${counts_desc}-coord2 with accession ${1} and k=${K}"
  mkdir -p metagraph
  (
  cd metagraph || error "can't cd to metagraph"
  infile="../${1}.fasta"
  outfile_base="${S}.met-coord2.k${K}"
  graph="${outfile_base}.dbg"
  annotations_binary="${outfile_base}.column.annodbg"
  annotations_reduced="${outfile_base}.row_diff_brwt.annodbg"
  annotations_coord="${outfile_base}.column.annodbg.coords"

  # use count_dbg commands (with coordinates, not counts!)
  if [[ ! -f "$graph" ]]; then
    echo "Building DBG..."
    metagraph build --parallel "${NTHREAD}" --kmer-length "${K}" -o "${graph}" "${infile}"
  else
    echo "Skipped."
  fi

  if [[ 1 == 1 || $2 == "--counts" && ! -f $annotations_binary && ! -f $annotations_reduced ]]; then
    # produce .column.annodbg (binary graph annotation) and .column.annodbg.coords (kmers coordinates) files
    echo "Building *.column.annodbg and *.column.annodbg.coords"
    metagraph annotate --anno-filename --coordinates\
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${infile}"

    # transform coordinates into row_diff (3 stages)
    echo "Transform coordinates into row_diff (3 stages)"
    echo "stage: 0"
    metagraph transform_anno --anno-type row_diff --coordinates --row-diff-stage 0\
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}.row_count" "${annotations_binary}"
    echo "stage: 1"
    metagraph transform_anno --anno-type row_diff --coordinates --row-diff-stage 1\
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}.row_reduction" "${annotations_binary}"
    echo "stage: 2"
    metagraph transform_anno --anno-type row_diff --coordinates --row-diff-stage 2\
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${annotations_binary}"

    # transform coordinates into row_diff_bwrt
    echo "Transforming coordinates into row_diff_bwrt"
    metagraph transform_anno --anno-type row_diff_brwt --greedy\
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${outfile_base}.row_diff.annodbg"

    # relax annotations
    echo "Relaxing annotations"
    metagraph relax_brwt --relax-arity 32\
    --parallel "${NTHREAD}" --infile-base "${graph}" --outfile-base "${outfile_base}" "${outfile_base}.row_diff_brwt.annodbg"
  fi

  # clean temp files
  rm *.pred *.pred_boundary *.rd_succ *.succ *.succ_boundary *.linkage *.row_count *.row_reduction *.anchors
  #for file in "${outfile_base}"*; do
  #  [[ $file == "$graph" || $file == "$annotations_binary" || $file == "$annotations_coord" ]] ||\
  #  [[ $file == *.lzma || $file == *.gz || $file == *.mfc ]] || rm "$file"
  #done

  write_to_csv "${graph}"
  compress_all_and_write_csv "${graph}"

  if [[ $2 == "--counts" ]]; then
    write_to_csv "${annotations_reduced}"
    write_to_csv "${annotations_coord}"
    compress_all_and_write_csv "${annotations_reduced}"
    compress_all_and_write_csv "${annotations_coord}"
  fi
  )
}

# launch metagraph2 before this!
launch_metagraph_assemble(){
  local method="assembly";
  local counts="no-counts"


  local metagraph_contigs="../metagraph2-k${K}/${1}.met.k${K}.contigs.fasta.gz"
  local metagraph_counts="../metagraph2-k${K}/${1}.met.k${K}.contigs.kmer_counts.gz"
  local outfile_base="${1}.assembled.k${K}"
  local outfile_contigs="${outfile_base}.fasta"
  local outfile_counts="${outfile_base}.kmer_counts"

  echo "*** Launching metagraph assemble with accession ${1} and k=${K}"
  mkdir -p $method
  (
  cd $method || error "cannot cd to metagraph"

  if [[ ! -f $outfile_contigs ]]; then
    # not --best compression, decompress!
    zless "${metagraph_contigs}" > "${outfile_contigs}"
  fi

  if [[ $2 == "--counts" && ! -f $outfile_counts ]]; then
     zless "${metagraph_counts}" > "${outfile_counts}"
  fi

  write_to_csv "${outfile_contigs}"
  compress_all_and_write_csv "${outfile_contigs}"

  if [[ $2 == "--counts" ]]; then
    write_to_csv "${outfile_counts}"
    compress_all_and_write_csv "${outfile_counts}"
  fi
  )
}

download_and_launch(){
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
      --skip-technical --split-3 --fasta \
      --outfile "${S}.fasta" "${S}.sra"

      # keep only the first file
      rm -f "${S}.fasta" "${S}_2.fasta"
      mv -f "${S}_1.fasta" "${S}.fasta"
    fi

    # iterating k-mer size
    for K in $KMER_SIZES; do
      #launch_metagraph_count "$S" "--counts"
      #launch_prophasm "$S"
      #launch_bcalm "$S"
      #launch_bcalm "$S" "--counts"
      #launch_ust "$S"
      #launch_ust "$S" "--counts"
      #launch_metagraph_coord1 "$S"
      #launch_metagraph_coord1 "$S" "--counts"
      launch_metagraph_assemble "$S" "--counts"
    done
    )
  done
  exit
}

setup
# clean
download_and_launch
finalize

# plot
source venv/bin/activate
./plot.py "$RESULTS"

# power off the computer
shutdown -h +3

exit 0