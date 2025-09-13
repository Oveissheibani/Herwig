#!/bin/bash
#SBATCH --job-name=herwig_pp13
#SBATCH --output=output/%a/herwig_pp13_%a.out
#SBATCH --error=output/%a/herwig_pp13_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:25:00
#SBATCH --mem=2G
#SBATCH --partition=mdtp
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@domain.com
#SBATCH --array=1-10

# ---- Environment setup ----
module purge
module load gnu7/7.3.0
module load cmake/3.21.1
module load gsl/2.5
module load root/6.28.10

# ---- Environment setup ----
export PATH=$HOME/Herwig/install/bin:$PATH
export LD_LIBRARY_PATH=$HOME/Herwig/install/lib:$HOME/PEG/install/lib:$HOME/PEG/install/lib/ThePEG:$LD_LIBRARY_PATH
export LHAPDF_DATA_PATH=$HOME/LHAPDF/install/share/LHAPDF

# ---- Compile density calculator (check if already compiled) ----
cd $HOME/Herwig
if [ ! -f density_calc ] || [ density_calculator.cpp -nt density_calc ]; then
    echo "Compiling density calculator..."
    g++ density_calculator.cpp -o density_calc `root-config --cflags --libs`
fi

# ---- Create output directory structure ----
mkdir -p output/${SLURM_ARRAY_TASK_ID}
cd output/${SLURM_ARRAY_TASK_ID}

# Copy input file to job directory
cp $HOME/Herwig/pp13TeV_basic.in .

# Prepare run file
Herwig read pp13TeV_basic.in

# Run the event generation with unique seed based on array task ID
Herwig run pp13TeV_basic_run.run -N 1000 -s ${SLURM_ARRAY_TASK_ID}

# ---- Run density calculator on the generated events ----
echo "Running density calculator on events.hepmc..."

# Copy the compiled executable to this directory
cp $HOME/Herwig/density_calc .

# Run density calculator
./density_calc events.hepmc output.root 211 321 2212 3122
#rm *.hepmc

echo "Job ${SLURM_ARRAY_TASK_ID} completed successfully!"
