# Stringtie

Stringtie takes the mapped reads and builds them into models of transcripts.

```
#!/bin/bash -e
#SBATCH -p nbi-medium # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=EMAIL # send-to address
#SBATCH --mem=32000
#SBATCH --cpus-per-task=8
#SBATCH --job-name="stringtie"
#SBATCH --array=0-n

source samtools-0.1.19
source stringtie-1.3.5

ARRAY=()
annotation_file=../TAIR10_GFF3_genes.gff

stringtie -p 8 -o ${ARRAY[$SLURM_ARRAY_TASK_ID]}_stringout.gtf -G $annotation_file -v Hisat2_out.${ARRAY[$SLURM_ARRAY_TASK_ID]}_sorted.bam
```

options used are
-o output file to write
-G reference annotation to use to guide transcript discovery (remember to use file with correct Chr/chr)
-v verbose mode

There are also strandedness options for libraries
--rf	Assumes a stranded library fr-firststrand.
--fr	Assumes a stranded library fr-secondstrand. 

-x <seqid_list> can be used to ignore mapping to certain reference features e.g. chromosomes. 
If you wanted to just map to the nuclear genome you could use -x 'chrM,chrC' for example
