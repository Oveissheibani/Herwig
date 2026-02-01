#!/bin/bash
#SBATCH --job-name=parallel_subsample
#SBATCH --output=parallel_subsample_%A_%a.out
#SBATCH --error=parallel_subsample_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --partition=mdtp
#SBATCH --array=0-99
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@domain.com

module purge
module load gnu7/7.3.0
module load root/6.28.10

# === CONFIGURATION ===
NUM_JOBS=100
JOB_ID=${SLURM_ARRAY_TASK_ID}

# Compile
g++ -std=c++17 -O2 -Wall -o parallelsubsample parallelsubsample.cpp `root-config --cflags --libs`
if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

BASE_DIR="./output"
# Choose one set of plot names (uncomment desired one)
#PLOT_NAMES="h_final_pair_density_2212_2112, h_final_pair_density_2212_m211, ..."  # (your full list)
PLOT_NAMES="h_Bs_211_211_class1"

# Target subdirectory in input files
SUBDIR="Class1/SymmetricBalanceFunctions_Bs"

# Output settings for partial averages
OUTPUT_PATH="./partial_averages"
OUTPUT_FILENAME="partial_avg_${JOB_ID}.root"

mkdir -p "${OUTPUT_PATH}"

echo "=== Parallel Subsample Job ${JOB_ID}/${NUM_JOBS} ==="
echo "Base dir: ${BASE_DIR}"
echo "Subdir: ${SUBDIR}"
echo "Output: ${OUTPUT_PATH}/${OUTPUT_FILENAME}"
echo "Folders containing ROOT files:"
find "${BASE_DIR}" -type f -name "*.root" -maxdepth 2 -exec dirname {} \; | sort -u
echo "----------------------------------------"

echo "Starting parallel subsample..."
./parallelsubsample "${BASE_DIR}" "${PLOT_NAMES}" "${SUBDIR}" "${OUTPUT_PATH}" "${OUTPUT_FILENAME}" "${NUM_JOBS}" "${JOB_ID}"

echo "Job ${JOB_ID} finished."
