# Mapping with Hisat2

```
#!/bin/bash -e
#SBATCH -p nbi-medium # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=EMAIL # send-to address
#SBATCH --mem=32000
#SBATCH --cpus-per-task=8
#SBATCH --job-name="hisat2"
#SBATCH --array=0-n

source hisat2-2.1.0
source samtools-0.1.19

ARRAY=()

hisat2 -p 8 -x ../TAIR10_wCM_rRNA_masked -U ./Trimmed_data/${ARRAY[$SLURM_ARRAY_TASK_ID]}_trimmed.fq.gz -S Hisat2_out.${ARRAY[$SLURM_ARRAY_TASK_ID]}.sam
```

options used are
-x the basename of the Hisat index
-U reads to be mapped (single end). For paired end its -1 A_1.fq -2 A_2.fq
-S output sam file

other options are defaults. Full options can be found at https://ccb.jhu.edu/software/hisat2/manual.shtml 
