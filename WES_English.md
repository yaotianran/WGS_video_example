# WES_English_example
For training projects

---
# WES sequencing data analysis tutorial

In this WES sequencing data analysis tutorial we will use a PE150 data sampled from human blood as an example to introduce how to discover mutations related to hunman genetic diease sickle cell anemia.
The workflow includes:

1. Preprocess
2. Quality control
3. Mapping / Alignment
4. Duplication removal
5. Variants calling
6. Variants annotation

---
## 1. Preprocess
Here we assume you have already had basic knowledge about Linux OS such as file manipulation, executing a program etc. 
We will use the following softwares and data in this tutorial.

Softwares:
- fastp
- bwa
- samtools
- sambamba
- strelka
- annovar

Data
- Two pair-end sequencing fils：raw_R1.fq.gz, raw_R2.fq.gz
- Human reference genome (hg19)：hg19.fasta

---
## 2. Quality control
Quality control is an essential step in total pipeline. We use fastp to evaluate two sequencing file and discard the low-quality data. By doing this
we assure that only high-quality will be used in alignment.

```
fastp \
	--in1 raw_R1.fq.gz \
	--in2 raw_R2.fq.gz \
	--out1 clean_R1.fq.gz \
	--out2 clean_R2.fq.gz \
	--json fastp.json \
	--html fastp.html \
	--detect_adapter_for_pe \
	--qualified_quality_phred 20 \
	--unqualified_percent_limit 40 \
	--n_base_limit 5 \
	--length_required 35 \
	--low_complexity_filter \
	--complexity_threshold 30 \
```
Output from fastp (partial)：
![Js2fCw7.png](https://iili.io/Js2fCw7.png)


fastp software manual：https://github.com/OpenGene/fastp?tab=readme-ov-file#all-options

---
## 3. Mapping / Alignment
Sequences mapping or Sequences alignment means the clean sequences from previous step are mapped onto reference genome so that we can locate each fragment in genome.
In this step we will use bwa, samtools and genome reference file hg19.fasta

When you are first time to use genome reference index files have to be made by the following commmand: 
```
bwa index hg19.fasta
```
After index files are generated 'bwa mem' command was executed to map sequences. Then results are sorted and indexed and are finally stored in BAM format.
```
bwa mem -M -R "@RG\tID:sampleid\tSM:sampleid" hg19.fasta clean_R1.fq.gz clean_R1.fq.gz |\
 samtools view -bS - |\
 samtools sort -o sorted.bam -

samtools index sorted.bam
```

bwa manual：https://bio-bwa.sourceforge.net/bwa.shtml

samtools manual：http://www.htslib.org/doc/samtools.html

Partial mapping results:
[![Jsa6GTb.png](https://iili.io/Jsa6GTb.png)](https://freeimage.host/)

---
## 4. Duplication removal
Duplicates can be problematic for variant calling as they can introduce systematic error (e.g. copying errors during PCR). Removing them is usually recommended for WGS WES libraries.
But it is generally not recommended for RNA sequencing data or ultra-high coverage gene panel data based on amplicon enrichment unless UMIs (Unique Molecular Identifier) were used.
Hereby we use "sambamba markdup" to mark duplicated sequences in the data so that variants caller programs will ignore them automatically.

```
sambamba markdup -r sorted.bam rmdup.sorted.bam --tmpdir tmp --overflow-list-size 1000000
```

After duplicated sequences were marked use 'sambamba flagstat' command to verfify the result.
```
sambamba flagstat rmdup.sorted.bam > rmdup.sorted.bam.flagstat
```

Sambamba manual: https://lomereiter.github.io/sambamba/docs/sambamba-view.html

Sambamba output (partial)：
![Jsa4aea.png](https://iili.io/Jsa4aea.png)

---
## 5. Variants calling
Variants calling refers to the use of high-throughput sequencing technology to search for DNA/RNA polymorphisms in species with known genome sequences, including single nucleotide polymorphisms (SNP), short indels (InDel), structural variations (SVs) and copy number variation (CNV), etc.
In this step we use manta to detect SVs and strelka to detect SNPs and InDels. Both software require two steps to run.


First， run 'config Manta' and 'configure Strelka Germline Workflow' commands to generate pipeline scripts. 
```
# Manta
configManta.py \
   --config configManta.py.ini \ 
   --bam rmdup.sorted.bam \
   --referenceFasta hg19.fasta \
   --runDir ./manta

# strelka
configureStrelkaGermlineWorkflow.py \
   --config configureStrelkaGermlineWorkflow.py.ini \
   --bam rmdup.sorted.bam \
   --referenceFasta hg19.fasta \
   --indelCandidates ./manta/results/variants/candidateSmallIndels.vcf.gz \
   --runDir ./strelka

```

Second，Run the scripts generateed in step 1.
```
# Manta
./manta/runWorkflow.py

# strelka
./strelka/runWorkflow.py
```
The results were stored in VCF format.

strelka manual：https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md

Variants calling results (partial)：
[![JsaiQhG.png](https://iili.io/JsaiQhG.png)](https://freeimage.host/)

---
## 6. Variant annotation
Variant annotation is the process by which variants and mutations in the DNA are assigned functional information and is a crucial process in genomic sequence analysis. 
The outcomes of such annotation are beneficial because they can directly influence the conclusions arrived at in disease studies.
There are many software and database that can help us search for the functional information。 Hereby we use ANNOVAR and the refGene database to annotate the mutation location and genes
```
table_annovar.pl \
   ./strelka/results/variants/variants.vcf.gz \
   humandb \
   --buildver hg19 \
   --remove \
   --protocol refGene \
   --operation g \
   --outfile annovar.txt \
   --nastring . \
   --polish \
   --vcfinput
```
Variant annotation results (partial)：
![Jsas5Cv.png](https://iili.io/Jsas5Cv.png)

ANNOVAR manual： https://annovar.openbioinformatics.org/en/latest/


