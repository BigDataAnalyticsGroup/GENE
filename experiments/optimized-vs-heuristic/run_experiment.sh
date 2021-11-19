# Check number of command line arguments
if [ "$#" -ne 0 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

# Create folders
mkdir -p data
mkdir -p workloads
mkdir -p results

# Specify arguments
EXPERIMENTS_DIR=${PWD}
cd ../../
ROOT_DIR=${PWD}
RELEASE_BUILD="${ROOT_DIR}/build/release"

OUTPUT="${EXPERIMENTS_DIR}/results/optimized_vs_heuristic.out"
PROGRAM="${RELEASE_BUILD}/bin/optimized_vs_heuristic"

# Compile project in release mode
make -C $RELEASE_BUILD


DATASETS=("uni-dense" "osm" "books")
DATASIZE=100000000 # 100M
WKLSIZE=10000000   # 10M

# Write header to OUTPUT file
echo n,distribution,workload_type,workload_size,index_type,inner_slot_size,leaf_slot_size,num_run,time > ${OUTPUT}

for data in "${DATASETS[@]}"
do
    # Generate data
    if [ ! -f "${EXPERIMENTS_DIR}/data/${data}_100M.data" ]; then
        ${RELEASE_BUILD}/bin/generate_data ${data} ${DATASIZE} ${EXPERIMENTS_DIR}/data/${data}_100M.data
    fi
    # Generate workload
    if [ ! -f "${EXPERIMENTS_DIR}/workloads/${data}_10M.data" ]; then
        ${RELEASE_BUILD}/bin/generate_opt-vs-heuristic_workload ${EXPERIMENTS_DIR}/data/${data}_100M.data ${WKLSIZE} ${EXPERIMENTS_DIR}/workloads/${data}_10M.wkl
    fi
    # Execute program
    ${RELEASE_BUILD}/bin/optimized_vs_heuristic --dataset="${EXPERIMENTS_DIR}/data/${data}_100M.data" --workload="${EXPERIMENTS_DIR}/workloads/${data}_10M.wkl" --output=${OUTPUT}
done
