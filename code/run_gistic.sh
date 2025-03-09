#!/bin/bash

cd ..

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <results_dir> <data_dir> <log_dir>"
    exit 1
fi

export results_dir="$1"
export data_dir="$2"
export log_dir="$3"

container_names=("CONTAINER_NAME_1" "CONTAINER_NAME_2" "CONTAINER_NAME_3")
num_workers=${#container_names[@]}

# Ensure log directorie exist
mkdir -p "$log_dir"

# function which takes a seg file and a container name as parameters and runs runs GISTIC
run_GISTIC() {
    local seg_file="$1"
    local container="$2"
    echo "Processing $seg_file using container $container"
    
    dir=$(dirname "$seg_file")
    filename=$(basename "$seg_file")
    filename_no_ext="${filename%.*}"
    
    output_dir="$(dirname "$seg_file")/GISTIC/$filename_no_ext"
    
    parent=$(basename "$(dirname "$seg_file")")           # returns cancer type
    grandparent=$(basename "$(dirname "$(dirname "$seg_file")")")  # returns arm
    log_file="$log_dir/$grandparent-$parent-$filename_no_ext.log"

    if [ -d "$output_dir" ]; then
        echo "Output exists for $output_dir. Skipping."
        return
    fi
    mkdir -p "$output_dir"
    
    docker run --volume "$(pwd)":/home "$container" -b "/home/$output_dir" \
      -seg "/home/$seg_file" \
      -mk "/home/$data_dir/genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt" \
      -refgene "/home/$data_dir/hg19_GENCODE_v18_20140127.mat" \
      -cnv "/home/$data_dir/SNP6.merged.151117.hg19.CNV.txt" \
      -brlen 0.7 -conf 0.99 -armpeel 1 -maxseg 2000 -genegistic 1 \
      -ta 0.1 -td 0.1 -cap 1.5 -js 4 > "$log_file" 2>&1
}

export -f run_GISTIC
export log_dir


# Create a temporary file that pairs each seg file with a container name using round-robin assignment
temp_file=$(mktemp)
index=0
for seg_file in $(find $results_dir  -type f -name "*.seg")
do
    # Select container name by modulating index with number of containers
    container="${container_names[$((index % num_workers))]}"
    echo "$seg_file $container" >> "$temp_file"
    index=$((index + 1))
done

# Launch jobs in parallel. Each invocation of run_GISTIC receives two arguments:
#   $1: seg_file and $2: container name.
cat "$temp_file" | xargs -P "$num_workers" -n 2 bash -c 'run_GISTIC "$1" "$2"' _

#rm "$temp_file"

echo "All tasks are done!"
