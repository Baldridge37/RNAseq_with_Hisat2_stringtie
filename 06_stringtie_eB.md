# Stringtie -e -B

Here we run stringtie again but with slightly different options.

By using the -e option we limit stringtie to estimating abundances only at reference transcripts supplied by -G
We could have used this earlier too if we didn't want stringtie to find novel transcripts.

The -B options writes the output ready to be read by ballgown.
So here we are creating a directory called ballgown, with a directory for each sample.
Within that will be each sample's output at A_abund.gtf for example.

```
#!/bin/bash -e
#SBATCH -p nbi-medium # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=EMAIL # send-to address
#SBATCH --mem=32000
#SBATCH --cpus-per-task=8
#SBATCH --job-name="stringtie"
#SBATCH --array=0-n

source stringtie-1.3.5

ARRAY=()

annotation_file=All_samples_transcripts.gtf

stringtie -e -B -p 8 -o ./ballgown/Sample_${ARRAY[$SLURM_ARRAY_TASK_ID]}/All_${ARRAY[$SLURM_ARRAY_TASK_ID]}_abund.gtf -G $annotation_file Hisat2_out.${ARRAY[$SLURM_ARRAY_TASK_ID]}_sorted.bam
```

As inputs it takes the bam files and the combined annotation file.

Note that if we had made stringtie find annotated genes only we could have skipped the stringtie merge step and jumped to ballgown.
