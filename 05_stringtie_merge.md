# Stringtie merge

This special use of Stringtie combines transcripts from different files into one master list of transcripts.
Rather than each sample having its own transcript names 'STRG.1, STRG.2' etc, they will be common to all samples.
The combined transcripts will be written to 'All_sample_transcripts.gtf'. We will then use this to calculate abundances

```
#!/bin/bash -e
#SBATCH -p nbi-medium # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=EMAIL # send-to address
#SBATCH --mem=16000
#SBATCH --cpus-per-task=8
#SBATCH --job-name="stringtie"

source stringtie-1.3.5

annotation_file=../TAIR10_GFF3_genes.gff

stringtie --merge  -p 8 *stringout.gtf -G $annotation_file -o All_sample_transcripts.gtf -v
```

other options for stringtie --merge are 
-G <guide_gff>	reference annotation to include in the merging (GTF/GFF3)
-o <out_gtf>	output file name for the merged transcripts GTF (default: stdout)
-m <min_len>	minimum input transcript length to include in the merge (default: 50)
-c <min_cov>	minimum input transcript coverage to include in the merge (default: 0)
-F <min_fpkm>	minimum input transcript FPKM to include in the merge (default: 0)
-T <min_tpm>	minimum input transcript TPM to include in the merge (default: 0)
-f <min_iso>	minimum isoform fraction (default: 0.01)
-i 	keep merged transcripts with retained introns (default: these are not kept unless there is strong evidence for them)
-l <label>	name prefix for output transcripts (default: MSTRG) 

https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
