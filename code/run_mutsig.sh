#!/bin/bash

cd ..

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <results_dir> <log_dir>"
    exit 1
fi

export results_dir="$1"
export log_dir="$2"

container_names=("mutsig_new" "mutsigcv2" "mutsigcv3" "mutsigcv4")
num_workers=${#container_names[@]}

# Ensure log directorie exist
mkdir -p "$log_dir"

# function which takes a mutations file and a container name as parameters and runs MutSig
run_MutSig() {
    local mutations_file="$1"
    local container="$2"
    echo "Processing $mutations_file using container $container"

    dir=$(dirname "$mutations_file")
    filename=$(basename "$mutations_file")
    filename_no_ext="${filename%.*}"
    
    output_dir="$(dirname "$mutations_file")/MutSig/$filename_no_ext"
    
    parent=$(basename "$(dirname "$mutations_file")")           # returns cancer type
    grandparent=$(basename "$(dirname "$(dirname "$mutations_file")")")  # returns arm
    log_file="$log_dir/$grandparent-$parent-$filename_no_ext.log"

    if [ -d "$output_dir" ]; then
        echo "Output exists for $output_dir. Skipping."
        return
    fi
    mkdir -p "$output_dir"
    
    docker run --volume "$(pwd)":/home \
    --workdir /home/mutsig2cv "$container" \
    ./MutSig2CV "/home/$mutations_file" "/home/$output_dir" > "$log_file" 2>&1
}

export -f run_MutSig
export log_dir


# Create a temporary file that pairs each seg file with a container name using round-robin assignment
temp_file=$(mktemp)
index=0
for mutations_file in $(find $results_dir  -type f -name "*aneuploid_mutations.tsv")
do
    # Select container name by modulating index with number of containers
    container="${container_names[$((index % num_workers))]}"
    echo "$mutations_file $container" >> "$temp_file"
    index=$((index + 1))
done

# Launch jobs in parallel. Each invocation of run_MutSig receives two arguments:
#   $1: mutations_file and $2: container name.
cat "$temp_file" | xargs -P "$num_workers" -n 2 bash -c 'run_MutSig "$1" "$2"' _

rm "$temp_file"

echo "All tasks are done!"
