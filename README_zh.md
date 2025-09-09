# 🧬HACGA1 软件介绍

🌐 [English](README.md) | [简体中文](README_zh.md)

HACGA1软件包通过Python流程将各个软件工具有效串联起来，通过函数调用和命令行执行等方式，实现了从数据预处理到基因组注释的自动化流程（**图1**）。

本手册包含以下内容：

- **原理说明**：详细介绍 HACGA1 的设计理念、工作机制以及在基因组注释中的应用原理。
- **安装方法**：指导用户在 Linux 系统环境下（推荐使用 **CentOS**）进行软件安装和流程部署。
- **依赖环境配置**：列出并解释所需的软件依赖、Python 包以及第三方工具的配置方法，确保流程能够顺利运行。
- **使用说明**：提供从输入数据准备到结果解析的完整操作流程示例，帮助用户快速上手。

通过阅读本手册，用户可以在 **Linux 系统（推荐 CentOS）** 上高效完成基因组注释分析，从而实现基因预测、结构注释与多维数据整合的自动化处理。

------



## 1. HACGA1 软件包原理介绍

HACGA1工具的开发逻辑首先从基因组数据的预处理开始，通过split_fasta_by_length()函数将大基因组按指定长度（如100kb）拆分为小片段，以便后续分析可以并行处理，提高计算效率。拆分后的小片段会作为输入传递给后续模块进行处理。在转录本组装过程中，HACGA1通过调用RNADenovo模块进行分析，该模块通过subprocess.run函数，使用StringTie工具，将BAM文件中的比对信息用于转录本组装，并输出为GTF格式文件，接着使用TransDecoder工具识别编码区，来识别转录本中的长开放阅读框。接着，在基因预测阶段，HACGA1同时使用AUGUSTUS、GlimmerHMM和GeneMark.ET这三个工具进行基因边界预测。通过调用augustus、glimmerhmm和genemarket函数，命令行调用相应工具可以从不同算法中获得多维度的基因预测结果，并结合RNA-seq数据进一步优化预测准确性。在RNA-seq数据的处理上，HACGA1通过RNADenovo模块生成初步转录本组装，再通过trinity_pasa模块进行转录本的校准，以确保最终基因结构的高精度预测，同样通过subprocess.run函数调用命令行执行操作。整个流程中的每一步都经过标准化设计，确保各模块之间的数据流畅传递并避免人工干预。在执行过程中，Python脚本通过subprocess模块调用外部工具，并通过错误捕获机制确保容错性，比如当StringTie失败时，脚本会记录错误并跳过当前步骤继续执行后续分析。这样即使某一环节出错，整个流程也能继续执行，避免了流程中断，提高了系统的鲁棒性。此外，HACGA1还结合了SLURM任务调度系统，通过分布式计算集群进一步提高分析效率。SLURM在任务调度中的应用使得HACGA1能够在多个计算节点上并行执行任务，尤其在处理大规模基因组数据时，能够有效地分配计算资源，减少任务执行时间。

具体执行过程中，对于RNA-seq数据，HACGA1使用Hisat2以及StringTie2软件进行比对和组装，得到初步的转录本组装注释结果，并利用程序进行过滤处理。接着，使用GeneMark.ET软件进行基因结构的从头预测，得到第一个基因预测结果。同时，RNA-seq数据还通过Trinity进行转录本组装，并进行后续的PASA分析。如果有全长转录组测序（ISO-seq）数据，可以通过SMRTLink直接生成高质量的全长转录本，用于验证和更新RNA-seq转录本结果。然后，基于整合后得到的高质量转录本信息，使用PASA进行初步的开放阅读框（ORF）预测，获得基因结构预测的注释文件。最后，基于PASA得到的注释结果，分别使用AUGUSTUS和GlimmerHMM进行基因结构预测，生成第二个和第三个预测结果。对于同源预测，选择来自OrthoDB v10数据库或者需要注释的物种对应的模式物种加近缘物种蛋白序列的集合获取高质量同源蛋白序列。使用Genemark.EP+将同源蛋白序列和参考基因组序列进行比对，获得外显子、内含子等边界信息，即第四个预测结果。

最后，EvidenceModler（EVM）将上述所有证据（包括从头预测、转录组证据和同源蛋白证据）按照不同权重值进行整合，生成最终的基因模型。使用PASA进一步处理和过滤这些基因模型，进行UTR注释和可变剪切注释，以确保整合得到的最终基因集具有高完整性。

<img src="pic/pic1.png" alt="pic1" />

**图1 HACGA1 结构注释流程示意图**

---

## 2. Getting Started

### 2.1 安装 HACGA1

```bash
# 获取代码
git clone https://github.com/liqianQ1/annotation_gene.git
cd annotation_gene

# 创建 Conda 环境并安装依赖
conda create -n ann python=3.12.0 -y
conda activate ann
```

### 2.2 Conda 安装依赖工具

| 工具                  | 功能                                         |
| --------------------- | -------------------------------------------- |
| hisat2                | RNA-seq 比对到基因组                         |
| stringtie             | 转录本组装和定量                             |
| transdecoder          | 从转录本预测蛋白编码序列                     |
| gffread               | GFF/GTF 文件操作、格式转换                   |
| pasa                  | 蛋白和转录本证据整合、基因结构更新           |
| augustus              | ab initio 基因预测                           |
| bam2hints             | 为 augustus 生成 hints 文件（augustus 自带） |
| glimmerhmm            | ab initio 基因预测                           |
| evidencemodeler (EVM) | 整合多种基因预测证据生成最终注释             |

安装示例：

```bash
conda install -c bioconda hisat2 -y
conda install -c bioconda stringtie -y
conda install -c bioconda transdecoder -y
conda install -c bioconda gffread -y
conda install -c bioconda pasa -y
conda install -c bioconda augustus -y
conda install -c bioconda glimmerhmm -y
conda install -c bioconda evidencemodeler -y
```

------

## 2.3 GeneMarkS / GeneMark-ES

该软件需从官方网站获取：

- 下载地址：[GeneMark 下载页面](https://genemark.bme.gatech.edu/GeneMark/license_download.cgi)
- 注册后下载并解压，配置环境变量：

```bash
export PATH=/path/to/gm_et_linux_64:$PATH
```

------

## 2.4 Slurm 配置（可选）

推荐使用 Slurm 调度系统进行多任务管理与加速：
参考 [Slurm 官方教程](https://www.schedmd.com/slurm/installation-tutorial/) 进行配置。

------

## 2.5 配置 HACGA1

配置文件路径：

```
config/config.txt
```

示例配置：

```txt
hisat2-build       = Anaconda3/envs/ann/bin/hisat2-build
python2            = /usr/bin/python2
slurm              = Pipeline/slurm/slurm_Duty.pl
queue              = q_all
transdecoder       = Anaconda3/envs/ann/opt/transdecoder
stringtie          = Anaconda3/envs/ann/bin/stringtie
gffread            = Anaconda3/envs/ann/bin/gffread
pasa               = Anaconda3/envs/ann/opt/pasa-2.5.2
augustus           = Anaconda3/envs/ann/bin/augustus
bam2hints          = Anaconda3/envs/ann/bin/bam2hints
augustus_extrinsic = Anaconda3/envs/ann/config/extrinsic/extrinsic.M.RM.E.W.cfg
gmes               = Software/gmes_linux_64_4
perl               = Anaconda3/envs/ann/bin/perl
trainGlimmerHMM    = Anaconda3/envs/ann/bin/trainGlimmerHMM
glimmerhmm         = Anaconda3/envs/ann/bin/glimmerhmm
evidencemodeler    = Anaconda3/envs/ann/opt/evidencemodeler-2.1.0
```

## 3. HACGA1 使用说明

### 3.1 输入文件说明

HACGA1 注释流程需要以下输入文件：

1. **输入序列（3类）：**
   - 基因组序列
   - 基因组重复注释序列
   - RNA 数据组装的 Trinity 序列

2. **辅助注释文件（5个）：**
   - **EST.tab**  
     - 两列，制表符分隔  
     - 第1列：样本名称  
     - 第2列：RNA 数据组装的 Trinity 序列的绝对路径
   - **RNAseq.tab**  
     - 三列，制表符分隔  
     - 第1列：样本名称  
     - 第2列：转录数据 Read 1 的绝对路径  
     - 第3列：转录数据 Read 2 的绝对路径  
     - 功能：提供完整的样本清单，明确每个样本及其对应的原始测序文件位置
   - **homolog.tab**  
     - 两列，制表符分隔  
     - 第1列：样本名称  
     - 第2列：同源注释所需的蛋白质序列数据的绝对路径
   - **genome.tab**  
     - 四列，制表符分隔  
     - 第1列、第2列：样本名称  
     - 第3列：基因组序列的绝对路径  
     - 第4列：基因组重复注释序列的绝对路径
   - **weights.txt**  
     - 三列，包括证据类别、类型和权重  
     - **证据类别**：`ABINITIO_PREDICTION`、`PROTEIN`、`TRANSCRIPT`  
     - **类型**：`GeneMark.hmm`、`GeneMark.hmmET`、`AUGUSTUS`、`GlimmerHMM`、`StringTie`、`transdecoder`  
     - **权重建议**：权重（pasa） >> 权重（同源注释） ≥ 权重（从头预测）

### 3.2 流程执行

软件使用时按照要求提供输入文件投递到节点上即可自动化完成基因结构注释.

示例命令：

```bash
python /pipline/annotation_gene/ann_gene.py \
  --Outputdir $PWD/ath1 \
  --Genome tab/genome.tab \
  --Homolog tab/homolog.tab \
  --RNAseq tab/RNAseq.tab \
  --EST tab/EST.tab \
  --config /pipline/annotation_gene/config/config.txt
```
软件输入参数解读：

--Outputdir 表示需要提供的结果输出目录，可以提供结果输出的绝对路径；

--Genome 表示需要提供的基因组物种和序列信息，可以按照要求存放在配置文件中；

--Homolog 表示需要提供的基因组注释使用的参考蛋白序列所在路径的文件；

--RNAseq 表示需要提供的基因组注释使用的转录组数据所在路径的文件；

--config 用于指定配置文件的路径。

程序运行完成后，可以获得不同证据支持的基因结构注释的gff文件和基因及蛋白序列文件。

### 3.2 HACGA1运行结果介绍

在软件运行的过程中，系统会根据设定的顺序自动创建多个目录，包括 01RNAdenovo、02trinity_pasa、03augustus、04genemarket、05genemarkep、06glimmerhmm、07evm和08evm_pasa，并在每个目录下自动执行相应的计算任务

![img](pic/pic2_1.jpg) ![img](pic/pic2_2.jpg) 

**图2 软件计算和输出结果目录**

在01RNAdenovo文件夹中将自动产生01_index、02_hisat2和03_stringtie文件夹，在这些文件中分别进行参考基因组索引准备、序列比对以及转录本组装和注释。StringTie组装产生转录本结构注释文件（merge.stringtie.assemble.gtf）以及序列文件（transcripts.fasta）。最终以获取到的RNADenovo.gff文件为最终的基于转录组数据组装得到的转录本注释文件（图3）。

![img](pic/pic3.jpg) 

**图3 RNADenovo.gff 基因组注释文件内容示例**

在02trinity_pasa文件中，这个目录包含了使用PASA流程进行转录本组装和注释的全套文件，主要分为两个功能模块，第一个是序列比对，主要包括gmap.spliced_alignments.gff3（GMAP比对结果）、minimap2.spliced_alignments.gff3（Minimap2比对结果）以及blat.spliced_alignments.gff3（BLAT比对结果）。第二个是组装的核心输出文件，主要包括pasa_cotton.assemblies.fasta（组装的转录本序列）、pasa_cotton.pasa_assemblies.gff3/bed/gtf（三种格式的注释文件）以及pasa_cotton.pasa_assemblies_described.txt（组装描述文件）以及TransDecoder预测得到的蛋白序列、CDS序列以及基因组坐标注释类文件。pasa.1.end.gff（图4）为最终利用转录本序列进行基因结构预测的注释文件。

![img](pic/pic4.jpg) 

**图4 pasa.1.end.gff 基因组注释文件内容示例**

在03augustus文件夹中，该文件夹包含Augustus基因预测工具的运行文件及训练数据,其核心脚本和输出文件，主要包括hints.gff（整合RNA-seq或同源蛋白等外部证据）、运行脚本如all.shell、hints.gff.sh、augustus.gff.sh以及输出结果：augustus.gff（最终基因预测注释）、augustus.out（原始输出）、stat.out（统计报告）。其中，augustus.gff（图5）为最终通过从头预测生成的基因结构文件。

![img](pic/pic5.jpg) 

**图5 augustus.gff 基因组注释文件内容示例**

04genemarket文件夹存储了GeneMark-ET/ES基因预测流程的完整运行文件和数据，主要用于对应物种蛋白质编码基因的预测与注释。主目录包含核心输入输出文件：intron.gff（内含子预测结果）、run.cfg（参数配置文件）、gmhmm.mod（训练完成的隐马尔可夫模型）以及最终预测结果genemark.gff（基因结构注释）。子目录output用于存储优化后的预测模型gmhmm.mod；info用于记录训练与预测过程的详细日志（如dna.gc.csv统计GC含量、training.trace跟踪模型迭代）；run目录包含分阶段训练的模型文件；data目录提供输入数据，包括基因组序列dna.fna、外部证据et.gff和训练集training.fna。最终，genemarket.gff（图6）作为最终的基于转录组数据的基因结构注释文件。

![img](pic/pic6.jpg) 

**图6 augustus.gff 基因组注释文件内容示例**

05genemarkep，该文件夹存储了GeneMark-ES/EP（基于自训练和外部证据的基因预测工具）的完整分析流程文件。核心文件包括data目录下的基因组序列dna.fna、训练集training.fna，以及外部证据文件（如prothint.gff、ep.gff、plus.gff）；模型训练：run目录包含分阶段训练的隐马尔可夫模型文件；预测结果主目录生成genemark_es.gtf（ES模式预测）和genemarkep.gff（EP模式预测），以及统计报告stat.out；prothint子目录存储同源蛋白比对证据（如diamond.out、spaln.gff），用于提升预测准确性。genemarkep.gff（图7）作为最终根据蛋白参考集生成的基因结构注释文件。

![img](pic/pic7.jpg) 

**图7 genemarkep.gff 基因组注释文件内容示例**

在06glimmerhmm中，该文件夹存储了GlimmerHMM基因预测工具的完整运行文件。主目录包含核心运行脚本和结果文件：glimmerhmm.gff.sh、glimmerhmm.out（原始输出）和stat.out（统计报告）。train子目录存放训练数据与模型文件：cds_00.fa和cds_01.fa为训练用编码序列，generate_glimmer_scripts.sh和run.sh控制训练流程；new_00目录存储训练生成的模型文件（如coding_0_40.model和noncoding_40_100.model）和关键参数文件（如exons.dat记录外显子特征，acc/don.errors统计剪接位点错误率，atg/stop.markov标记起始/终止密码子模型）。glimmerhmm.gff（图8）为最终从头预测的基因结构文件。

![img](pic/pic8.jpg) 

**图8 glimmerhmm.gff 基因组注释文件内容示例**

 在07evm文件夹中，存储了EVM基因预测整合流程的完整数据，用于整合多种证据（蛋白比对、转录本、基因预测）生成最终的棉花基因注释。主目录包含核心输入文件：protein_alignments.gff（同源蛋白比对证据）、transcript_alignments.gff（RNA-seq转录本支持）、gene_predictions.gff（Augustus/GeneMark等工具的预测结果），以及控制文件weights.txt（证据权重配置）和commands.sh。运行后生成的分区文件（如group3_1-10000000）按基因组坐标分割证据数据以提高效率，最终输出evm.out.gff3（整合的基因模型）。子目录（如group3/、ptg000162l/）对应不同染色体或scaffold的分区分析，均包含输入证据文件、参考基因组cottons.racon.2.adjust.fa及分区结果。evm.gff3（图9）为最终通过 EVM 软件整合不同预测方法的基因结构文件。

![img](pic/pic9.jpg) 

**图9 evm.gff3 基因组注释文件内容示例**

08evm_pasa文件夹存储了EVM与PASA联合分析的完整流程文件，主要用于整合多源证据生成棉花基因组的高置信度注释。核心文件包括：evm_pasa.end.sh（主控脚本）和alignAssembly.config/annotationCompare.config（PASA配置文件）用以指导证据整合与基因模型更新；输出结果包含PASA修正的基因结构、最优蛋白预测和统计文件protein.best.gff.stat.out/Genes_annotation.statistics.xls（注释质量评估）。表明该流程成功整合了蛋白同源性、转录本比对等证据，最终生成经过PASA验证的可靠基因组注释。protein.best.gff（图10）本流程最终整合的基因注释结果。

![img](pic/pic10.jpg)  

**图10 protein.best.gff基因组注释文件内容示例**

这些注释文件最终为基因结构的完善和验证提供了全面的数据支持。



## 常见问题与建议

- **软件版本冲突**：建议使用独立 Conda 环境。
- **路径问题**：所有配置路径需为绝对路径。
- **大数据加速**：推荐配置 Slurm 集群，提高运行效率。
