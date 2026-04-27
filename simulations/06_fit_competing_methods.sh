#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --account=cdslab
#SBATCH --job-name=06_sim_compare
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --time=2:00:00
#SBATCH --output=out/%a
#SBATCH --error=err/%a
#SBATCH --array=119,123-125,149-154,157,169-178,186-193,201-239 

#9,20,23-34,57-87,91,93,101,104,108,111,118, 

module load R
echo $SLURM_ARRAY_TASK_ID

LINE_NUMBER=$((SLURM_ARRAY_TASK_ID))

n_clocks=$(awk -F',' "NR==${LINE_NUMBER} { print \$1; exit }" parameter_combinations.txt)
n_events=$(awk -F',' "NR==${LINE_NUMBER} { print \$2; exit }" parameter_combinations.txt)
purity=$(awk -F',' "NR==${LINE_NUMBER} { print \$3; exit }" parameter_combinations.txt)
coverage=$(awk -F',' "NR==${LINE_NUMBER} { print \$4; exit }" parameter_combinations.txt)
mutations_density=$(awk -F',' "NR==${LINE_NUMBER} { print \$5; exit }" parameter_combinations.txt)

# awk -F',' "NR==${LINE_NUMBER} { print ; exit }" parameter_combinations.txt >> tickTack_sim_${n_clocks}_${n_events}_${purity}_${coverage}/config


echo $n_clocks
echo $n_events
echo $purity
echo $coverage
echo $mutations_density

echo "NR==${LINE_NUMBER}"

Rscript 06_competing_single_segments_methods_comparison.R ${n_clocks} ${n_events} ${purity} ${coverage} ${mutations_density}
