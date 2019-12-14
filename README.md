---
title: 'DeepTE'
disqus: hackmd
---

DeepTE
===
DeepTE is aimed to classify transposons with unknown classification *via* Convolutional Neural Network.


# Introduction
Transposable elements (TEs) classification is an essential step decoding their roles in a genome. With reference genomes from non-model species available, it has begun to overstep efforts to annotate TEs, and more tools are needed to efficiently handle the emerged sequence information. We developed a novel tool, DeepTE, which classifies unknown TEs on basis of convolutional neural network. DeepTE utilized co-occurrence of k-mers towards TE sequences as input vector, and seven k-mer size was testified to be suitable for the classification. Eight models have been trained for different TE classification purposes. DeepTE applied domains from TEs to correct false classification. An additional model was also trained to distinguish between non-TEs and TEs targeting plant species. 

Given exclusive TEs of different species type, it can therefore classify seven orders, and 11-24 superframilies towards Plants, Metazoans, Fungi, and Others. This tool successfully leverages convolutional neural network to TE classification, assisting to precisely identify and annotate TEs in a new sequenced eukaryotic genome.

# Dependence and requirements
DeepTE is developed in Python with modules and external tools.

Before running this pipeline, a dependency check should be performed first to make sure every dependency is correctly installed.

For information about installing the dependencies, please see below. The version numbers listed below represents the version this pipeline is developed with, and using the newest version is recommended.

## Requirements
**Python** (v3.7 or more)  
Modules can be installed using pip: pip install -r requirements.txt or pip install [module_name]  
**Module version**  
biopython (1.72)  
keras (2.2.4)  
tensorflow (1.14.0)  
numpy (1.16.0)  

## Optional requirements
**HMMER** v3.1b1

**Model_dir**  
Download the model dir from link  
Plants:
https://drive.google.com/file/d/1voj86STKcQH8lAhvY6yl5E65nzaM6o0B/view?usp=sharing
Metazoans:
https://drive.google.com/file/d/1ExRwC3szJ4XMa3ikxM9Ccu31lY79rdw9/view?usp=sharing
Fungi:
https://drive.google.com/file/d/1uvnm99ypauIKtqCxoybdtT-mEMdoupip/view?usp=sharing
Others:
https://drive.google.com/file/d/1Q6HW1NhNs0a6Ykrw7jGEKKPWxawpWiuM/view?usp=sharing
UNS model:
https://drive.google.com/file/d/1uXTEtNQtJc2DO-JpT0s4Kv1k2ogUjCLr/view?usp=sharing

# Usage
```
usage:
**DeepTE**
DeepTE.py [-h] required: [-d working_dir][-o output_dir]
                         [-i ipt_seq][-sp sp_type]
                         ([-m model_name]|[-m_dir model_dir])
               optional: [-modify domain_file]
                         [-fam te_fam][-UNS yes]

arguments:
-h, --help        Show this help message and exit.

-d                Working directory to store intermediate files of each step. 
                  Default: ./.

-o                Output directory to store the output files. 
                  Default: ./.

-i                Input sequences that are unknown TE or DNA sequences.

-sp               P or M or F or O. P:Plants, M:Metazoans, F:Fungi, and O: Others.

-m                Provide one of model names: 
                  '-m P' or '-m M' or '-m F' or '-m O' or '-m U'.
                  This argument will directly download the model dir.
                  Users do not need to initiate '-m_dir'.
                  If users do not want to directly download model, please use '-m_dir', but users need to download model directory by themselves.

-m_dir            Provide model_dir that could be downloaded from website (optional requirements). 
                  If users set -UNS yes, please provide UNS_model directory that can be downlowed in the above link.

-fam              Provide TE family name for the input te sequence
                  Default: All
                  All: the input squence is unknown TEs
                  ClassI: the input sequence is ClassI TEs
                  ClassII: the input sequence is ClassII subclass1 TEs
                  LTR: the input sequence is LTR TEs
                  nLTR: the input sequence is nLTR TEs
                  LINE: the input sequence is LINE TEs
                  SINE: the input sequence is SINE TEs
                  Domain: the input sequence is Class II subclass1 TEs with specified super families

-modify           If set this argument, users need to provide domain file generated from another script: DeepTE_domain.py.

-UNS              If set this argument, users need change the -i to the the DNA sequences; 
                  This function will classify the sequences into TEs, CDS, or Intergenic sequences; -sp and -fam do not need to provide.
                  Note: this model is used for plants rather than metazoans and fungi.

**DeepTE_domain**
DeepTE_domain.py [-h] required: [-d working_dir][-o output_dir]
                                [-i ipt_seq][-s supfile_dir]
                                [--hmmscan hmmscan]
arguments:
-h, --help        Show this help message and exit.

-d                Working directory to store intermediate files of each step. 
                  Default: ./.

-o                Output directory to store the output files. 
                  Default: ./.

-i                Input sequences that are unknown TE sequences.

-s                Provide supplementary dir that contains required files.

--hmmscan         File path to hmmscan executable, Default: /usr/bin/hmmscan"


```

# Examples
**DeepTE.py**  
**Input data**  
Sequence data (fasta format)  

**Output data**  
Working directory  
a. opt_input_CNN_data.txt (input data that is transfered from user provided input data)  
b. store_temp_opt_dir (a directory contains prediction results for each TE group)
c. download_X_model_dir (store downloaded models. X represents P, M, F, O, or U)

Output directory  
a. opt_DeepTE.txt (a txt file with two columns. first column: original name; second column: predicted name with DeepTE)  
b. opt_DeepTE.fasta (a fasta file with new predicted TE name)  

**Command**
- [ ] Classify unknown TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P  
- [ ] Classify Class I TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam ClassI
- [ ] Classify Class II subclass1 TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam ClassII
- [ ] Classify LTR TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam LTR
- [ ] Classify nLTR TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam nLTR
- [ ] Classify LINE TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam LINE
- [ ] Classify SINE TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam SINE
- [ ] Classify TEs into MITEs and nMITEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam Domain  
- [ ] Classify Unknown sequences into TEs, Coding sequences, or Intergenic sequences  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -m U -UNS yes  

**DeepTE_domain.py**  
**Input data**  
Sequence data (fasta format)  
**Command**  
DeepTE_domain.py -d working_dir -o output_dir -i input_seq.fasta -s supfile_dir --hmmscan   Path/to/hmmscan  
**Output data**  
a. TE domain file. (first column: orignial name; second column: domain information for different open reading frames)






# Work flows
![](https://i.imgur.com/RlCZblM.png)


**Figure** A pipeline for classifying unknown TEs and sequences based on trained nine models. The unknown TEs go through eight models to be classified into different families. Two correction steps are conducted during classification. In Class model, TR domain exists in predicted Class I TEs that will be corrected to Class II_sub1 TEs, while RT domain exists in predicted Class II_sub1 TEs that will be corrected to Class I TEs. In ClassI model, EN domain exists in predicted LTR TEs will be correct to nLTR TEs. The unknown sequences go through UNS model to be classified into TEs, coding sequences (CDS), and intergenic sequences (INS).


## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::


