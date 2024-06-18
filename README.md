# WGS_video_example
录制伯杰教材项目的素材

---
# WGS数据处理全部步骤

在这个WGS测序数据分析的教程中，我们将以一组PE150的测序FASTQ数据（样本类型：1例人全血样本）为例介绍WGS测序数据的人类常规遗传性疾病检测分析流程。全部的分析步骤包括：

1. 前期准备
2. 数据质量控制
3. 序列比对
4. 去除重复序列
5. 变异检测
6. 变异注释

---
## 1. 前期准备
这里我们假定您已经具备了基础的Linux文件操作能力，我们将使用一组开源的生物信息学软件和数据文件来完成以下全基因组测序比对流程

需要下载安装的开源软件有：
- fastp
- bwa
- samtools
- sambamba
- strelka
- annovar

需要的数据文件有：
- 两个测序仪的输出文件：raw_R1.fq.gz, raw_R2.fq.gz
- hg19版本的人类基因组参考文件：hg19.fasta

---
## 2. 数据质量控制
数据质量控制是全基因组测序数据分析的关键步骤之一，这里使用软件fastp对两个测序仪的输出文件进行质量评估和过滤，以确保我们分析的数据质量良好，减少后续分析中的误差.

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
部分fastp输出结果：
![Js2fCw7.png](https://iili.io/Js2fCw7.png)


fastp软件文档：https://github.com/OpenGene/fastp?tab=readme-ov-file#all-options

---
## 3. 序列比对
序列比对是指将上一步质量控制的输出文件中的序列，比对到人类基因组上，此一步的目的是确定每一个序列在基因组的位置。
该步骤将用到开源软件bwa，samtools和人类参考基因组的fasta数据。
在第一次使用参考基因组fasta数据的时候，需要使用bwa index命令，创建参考基因组fasta数据的索引文件。
```
bwa index hg19.fasta
```
之后就可以使用bwa mem命令，比对的结果经过排序，以BAM格式储存在文件中
```
bwa mem -M -R "@RG\tID:sampleid\tSM:sampleid" hg19.fasta clean_R1.fq.gz clean_R1.fq.gz |\
 samtools view -bS - |\
 samtools sort -o sorted.bam -
```

最后不要忘了将排序后的BAM文件再次index
```
samtools index sorted.bam
```

bwa软件文档：https://bio-bwa.sourceforge.net/bwa.shtml

samtools软件文档：http://www.htslib.org/doc/samtools.html

部分比对后结果
[![Jsa6GTb.png](https://iili.io/Jsa6GTb.png)](https://freeimage.host/)

---
## 4. 去除重复序列
为了下游分析的准确，一般需要将，由PCR操作产生的重复序列去除或标记出来。
在这一步骤中，我们使用软件sambamba去掉比对结果中的重复序列。这里的输入是我们在上一步产生的比对后的BAM文件，sambamba markdup命令的输出就是去掉重复系列以后的BAM文件
```
sambamba markdup -r sorted.bam rmdup.sorted.bam --tmpdir tmp --overflow-list-size 1000000
```

去重之后我们可以使用sambamba flagstate命令查看一下效果
```
sambamba flagstat rmdup.sorted.bam > rmdup.sorted.bam.flagstat
```

Sambamba软件文档: https://lomereiter.github.io/sambamba/docs/sambamba-view.html

部分Sambamba输出结果：
![Jsa4aea.png](https://iili.io/Jsa4aea.png)

---
## 5. 变异检测
变异检测是指利用高通量测序技术对已知基因组序列的物种，寻找DNA水平的变异和多态性，包括单核苷酸多态性、短插入缺失、结构变异和拷贝数变异等
在这一步骤中我们使用manta检测结构变异，使用strelka检测单核苷酸多态性和短插入缺失。这两个软件运行都需要两个步骤。

第一步，运行config Manta命令和configure Strelka Germline Workflow命令，生成检测脚本。
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

第二步，运行第一步生成的检测脚本runWorkflow
```
# Manta
./manta/runWorkflow.py

# strelka
./strelka/runWorkflow.py
```
变异检测的结果以VCF格式储存在文件中，可以使用文本编辑器（记事本或写字板打开）

strelka软件文档：https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md

部分变异检测的结果：
[![JsaiQhG.png](https://iili.io/JsaiQhG.png)](https://freeimage.host/)

---
## 6. 变异注释
这里我们对上一步得到的变异检测结果进行注释。注释的目的，是获取变异的表型信息，例如变异属于哪个基因，在该基因的哪个位置等等。
因为上一步的变异结果以标准的VCF文件格式储存，所以可以使用多种公共开源软件进行注释。这里我们用ANNOVAR软件举例。
该软件安装完成后自带人类参考基因组hg19版本的refGene数据库，它可以注释变异位点位于哪个基因、位于基因的什么位置以及该变异是否引起蛋白质改变等信息
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
部分变异注释的结果：
![Jsas5Cv.png](https://iili.io/Jsas5Cv.png)

ANNOVAR软件文档： https://annovar.openbioinformatics.org/en/latest/


