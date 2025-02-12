# annotation_gene
A software for structural annotation of the assembled genome.
This workflow employs a comprehensive genomic annotation strategy that integrates multiple data sources and software tools to automate gene structure annotation for eukaryotic genomes. Three strategies are used simultaneously for gene annotation: ab initio annotation, transcriptomic evidence annotation, and homologous protein sequence annotation. Prior to structural annotation, RepeatModeler2 is used to construct new repeat sequence models, and the genome undergoes repeat sequence masking to reduce false positive results in subsequent analyses.

Gene structure prediction combines ab initio prediction with homologous prediction. For RNA-seq data, Hisat2 and StringTie2 are used for alignment and annotation, generating preliminary transcript annotations. Then, GeneMark.ET is used for ab initio gene structure prediction, producing the first gene prediction result. Additionally, RNA-seq data is assembled using Trinity, followed by PASA analysis. If full-length transcript sequencing (ISO-seq) data is available, high-quality full-length transcripts can be generated directly using SMRTLink to validate and update the RNA-seq transcript results. Next, based on the integrated high-quality transcript data, PASA is used for initial open reading frame (ORF) prediction. Finally, based on the annotations from PASA, Augustus and GlimmerHMM are used to predict gene structures, generating the second and third prediction results.

For homologous prediction, this method selects protein sequence sets from the OrthoDB v10 database or from model species and closely related species of the target organism to obtain high-quality homologous protein information. Genemark.EP is used to align homologous protein sequences with the reference genome sequence, providing exon, intron, and other boundary information, which constitutes the fourth prediction result.

Finally, EvidenceModeler integrates all the evidence (including ab initio predictions, transcriptomic evidence, and homologous protein evidence) to generate the final gene model. PASA is further used to process and filter these gene models, performing untranslated region (UTR) annotation and alternative splicing annotation to ensure that the final integrated gene set is of high quality. The accuracy of the pipeline is validated by manually curating immune-related gene annotations in the pig genome, especially genes related to the major histocompatibility complex (MHC), such as the swine leukocyte antigen (SLA).

##  software version
hisat2-build  2.2.1
Python2  2.7.5
transdecoder  5.7.1
stringtie  2.2.3
gffread  0.12.7
pasa  2.5.2
augustus  3.5.0
perl  5.32.1
evidencemodeler  2.1.0

##
To ensure compatibility and stability of the project, the following environment configurations are recommended:

Recommended Python Version: Python 3.12.0
Recommended Operating System: CentOS 7
