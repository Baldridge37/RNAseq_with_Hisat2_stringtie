###Trimming the RNAseq data

Depending on how your library was made you might want to trim your data differently. 

I used trim_galore with the illumina settings.
This will remove adapters and low quality bases (score <20)

```
#!/bin/bash -e
###SBATCH -p normal # partition (queue)
#SBATCH -p nbi-medium
#SBATCH -t 0-48:00 # time (D-HH:MM)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=billy.aldridge@jic.ac.uk # send-to address
#SBATCH --array=0-n
#SBATCH --mem=4000

ARRAY=()

#running with illumina defaults.
srun ~/group-data/Martin_Vickers/bin/trim_galore --illumina ./Raw_data/${ARRAY[$SLURM_ARRAY_TASK_ID]}.fastq.gz
```

I performed this as an array to trim all files at once. The output will the (name in array)_trimmed.fq
There will also be a trimming report created.
