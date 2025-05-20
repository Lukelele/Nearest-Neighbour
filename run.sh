#!/bin/bash
#SBATCH --account=chem030406
#SBATCH --job-name=cpp_compile_and_run
#SBATCH --partition=teach_cpu
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --time=100:00
#SBATCH --mem-per-cpu=100M

cd "${SLURM_SUBMIT_DIR}"


srun --mpi=pmi2 ./main argon120.xyz 0 9 > output120_0.txt
srun --mpi=pmi2 ./main argon120.xyz 1 9 > output120_1.txt
srun --mpi=pmi2 ./main argon120.xyz 2 9 > output120_2.txt
srun --mpi=pmi2 ./main argon120.xyz 3 9 > output120_3.txt
srun --mpi=pmi2 ./main argon120.xyz 4 9 > output120_4.txt
srun --mpi=pmi2 ./main argon120.xyz 5 9 > output120_5.txt
srun --mpi=pmi2 ./main argon120.xyz 6 9 > output120_6.txt
srun --mpi=pmi2 ./main argon120.xyz 7 9 > output120_7.txt
srun --mpi=pmi2 ./main argon120.xyz 8 9 > output120_8.txt
srun --mpi=pmi2 ./main argon120.xyz 9 9 > output120_9.txt

srun --mpi=pmi2 ./main argon10549.xyz 0 9 > output10549_0.txt
srun --mpi=pmi2 ./main argon10549.xyz 1 9 > output10549_1.txt
srun --mpi=pmi2 ./main argon10549.xyz 2 9 > output10549_2.txt
srun --mpi=pmi2 ./main argon10549.xyz 3 9 > output10549_3.txt
srun --mpi=pmi2 ./main argon10549.xyz 4 9 > output10549_4.txt
srun --mpi=pmi2 ./main argon10549.xyz 5 9 > output10549_5.txt
srun --mpi=pmi2 ./main argon10549.xyz 6 9 > output10549_6.txt
srun --mpi=pmi2 ./main argon10549.xyz 7 9 > output10549_7.txt
srun --mpi=pmi2 ./main argon10549.xyz 8 9 > output10549_8.txt
srun --mpi=pmi2 ./main argon10549.xyz 9 9 > output10549_9.txt

srun --mpi=pmi2 ./main argon147023.xyz 0 9 > output147023_0.txt
srun --mpi=pmi2 ./main argon147023.xyz 1 9 > output147023_1.txt
srun --mpi=pmi2 ./main argon147023.xyz 2 9 > output147023_2.txt
srun --mpi=pmi2 ./main argon147023.xyz 3 9 > output147023_3.txt
srun --mpi=pmi2 ./main argon147023.xyz 4 9 > output147023_4.txt
srun --mpi=pmi2 ./main argon147023.xyz 5 9 > output147023_5.txt
srun --mpi=pmi2 ./main argon147023.xyz 6 9 > output147023_6.txt
srun --mpi=pmi2 ./main argon147023.xyz 7 9 > output147023_7.txt
srun --mpi=pmi2 ./main argon147023.xyz 8 9 > output147023_8.txt
srun --mpi=pmi2 ./main argon147023.xyz 9 9 > output147023_9.txt

