# Building the Hisat2 reference

Before we start we need to build a reference genome to map our RNAseq reads to. 
This is pretty straightforward but there are a few steps where things can go wrong and cause problems later on.
For this analysis we will need a reference sequence (TAIR10_reference_wCM.fas)

```
>chr1
CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCAT
GAATCCCTAAATACCTAATTCCCTAAACCCGAAACCGGTTTCTCTGGTTGAAAATCATTGTGTATATAATGATAATTTT
ATCGTTTTTATGTAATTGCTTATTGTTGTGTGTAGATTTTTTAAAAATATCATTTGAGGTCAATACAAATCCTATTTCT
TGTGGTTTTCTTTCCTTCACTTAGCTATGGATGGTTTATCTTCATTTGTTATATTGGATACAAGCTTTGCTACGATCTA
CATTTGGGAATGTGAGTCTCTTATTGTAACCTTAGGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACTGTTTA
TGTTTGGACATTTATTGTCATTCTTACTCCTTTGTGGAAATGTTTGTTCTATCAATTTATCTTTTGTGGGAAAATTATT
TAGTTGTAGGGATGAAGTCTTTCTTCGTTGTTGTTACGCTTGTCATCTCATCTCTCAATGATATGGGATGGTCCTTTAG
CATTTATTCTGAAGTTCTTCTGCTTGATGATTTTATCCTTAGCCAAAAGGATTGGTGGTTTGAAGACACATCATATCAA
AAAAGCTATCGCCTCGACGATGCTCTATTTCTATCCTTGTAGCACACATTTTGGCACTCAAAAAAGTATTTTTAGATGT
```

a gene annotation (TAIR10_GFF3-genes.gff)

```
Chr1    TAIR10  chromosome      1       30427671        .       .       .       ID=chr1;Name=chr1
Chr1    TAIR10  gene    3631    5899    .       +       .       ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
Chr1    TAIR10  mRNA    3631    5899    .       +       .       ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1;Index=1
Chr1    TAIR10  protein 3760    5630    .       +       .       ID=AT1G01010.1-Protein;Name=AT1G01010.1;Derives_from=AT1G01010.1
Chr1    TAIR10  exon    3631    3913    .       +       .       Parent=AT1G01010.1
Chr1    TAIR10  five_prime_UTR  3631    3759    .       +       .       Parent=AT1G01010.1
Chr1    TAIR10  CDS     3760    3913    .       +       0       Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1    TAIR10  exon    3996    4276    .       +       .       Parent=AT1G01010.1
Chr1    TAIR10  CDS     3996    4276    .       +       2       Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1    TAIR10  exon    4486    4605    .       +       .       Parent=AT1G01010.1

```
And an annotation of what we don't want to map to (i.e. rRNA genes, TAIR10_GFF3-rRNA_XF.gff)

```
chr2    TAIR10  rRNA    3706    5513    .       +       .       ID=AT2G01010.1;Parent=AT2G01010;Name=AT2G01010.1;Index=1
chr2    TAIR10  rRNA    5782    5945    .       +       .       ID=AT2G01020.1;Parent=AT2G01020;Name=AT2G01020.1;Index=1
chr3    TAIR10  rRNA    14197677        14199484        .       +       .       ID=AT3G41768.1;Parent=AT3G41768;Name=AT3G41768.1;Index=1
chr3    TAIR10  rRNA    14199753        14199916        .       +       .       ID=AT3G41979.1;Parent=AT3G41979;Name=AT3G41979.1;Index=1
chrC    TAIR10  rRNA    101012  102502  .       +       .       ID=ATCG00920.1;Parent=ATCG00920;Name=ATCG00920.1;Index=1
chrC    TAIR10  rRNA    104691  107500  .       +       .       ID=ATCG00950.1;Parent=ATCG00950;Name=ATCG00950.1;Index=1
chrC    TAIR10  rRNA    107599  107701  .       +       .       ID=ATCG00960.1;Parent=ATCG00960;Name=ATCG00960.1;Index=1
chrC    TAIR10  rRNA    107949  108069  .       +       .       ID=ATCG00970.1;Parent=ATCG00970;Name=ATCG00970.1;Index=1
chrC    TAIR10  rRNA    130580  130700  .       -       .       ID=ATCG01160.1;Parent=ATCG01160;Name=ATCG01160.1;Index=1
chrC    TAIR10  rRNA    130948  131050  .       -       .       ID=ATCG01170.1;Parent=ATCG01170;Name=ATCG01170.1;Index=1
```

Now we have these three files you may notice a slight problem. 
In the .fas file and the rRNA GFF the chromosomes are labelled with lowercase letters, while in the gene file they are "Chr".

It's easier, and probably safer, to just change it in the gene folder using sed.

```
sed '{s/Chr/chr/g}' TAIR10_GFF3-genes.gff > TAIR10_GFF3-genes_correct.gff
```

Now we will be able to combine the gene annotations with the output of stringtie.


## Masking rRNA genes

When performing RNAseq analysis we don't want to map to rRNA genes.
There are different ways to remove them but for this analysis I think the easiest is to produce a copy of the genome with rRNA genes masked.

To do this we will use bedtools maskfasta (https://bedtools.readthedocs.io/en/latest/content/tools/maskfasta.html).

On an interactive node I performed.

```
bedtools-2.27.1
$ bedtools maskfasta [OPTIONS] -fi TAIR10_reference_wCM.fas -bed TAIR10_GFF3-rRNA_XF.gff -fo TAIR10_reference_wCM_rRNA_masked.fas
```

The locations of the rRNA genes stored in the GFF file will now be replaced by Ns in the fasta file.

## Building the hisat2 reference

We can then build the Hisat2 reference with our masked fasta file with a SLURM submission script.

```
#!/bin/bash -e
#SBATCH -p nbi-medium # partition (queue)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=EMAIL # send-to address
#SBATCH --mem=32000
#SBATCH --cpus-per-task=8
#SBATCH --job-name="hisat2"

source hisat2-2.1.0

hisat2-build TAIR10_reference_wCM_rRNA_masked.fas TAIR10_wCM_rRNA_masked
```
We will then have a series of .ht2 files starting with TAIR10_wCM_rRNA_masked that Hisat2 will map to.
And we're ready to move onto mapping

