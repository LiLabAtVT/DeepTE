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
**Phython** (v3.7 or more)
Modules can be installed using pip: pip install -r requirements.txt or pip install [module_name]
**Model_dir**
Download the model dir from link

Plants:
https://drive.google.com/file/d/1E6ZPGjlFcszJjthJtKO3KOMh05xYzq_L/view?usp=sharing

Metazoans:
https://drive.google.com/file/d/18qCmKnGjZykoDMnKeDHWpKmuaOb_hjUg/view?usp=sharing

Fungi:
https://drive.google.com/file/d/1DqeZPhZA1J9QFQXE2cPsCB_0em52aYai/view?usp=sharing

Others:
https://drive.google.com/file/d/1Q6HW1NhNs0a6Ykrw7jGEKKPWxawpWiuM/view?usp=sharing

## Optional requirements
**HMMER** v3.1b1



# Usage
```
usage:
**DeepTE**
DeepTE.py [-h] required: [-d working_dir][-o output_dir]
                         [-i ipt_seq][-m model_dir]
                         [-sp sp_type][-fam te_fam]
               optional: [-modify domain_file]

arguments:
-h, --help        Show this help message and exit.

-d                Working directory to store intermediate files of each step. 
                  Default: ./.

-o                Output directory to store the output files. 
                  Default: ./.

-i                Input sequences that are unknown TEs or DNA sequences.

-m                Provide model_dir that could be downloaded from website.

-sp               P or M or F or O. P:Plants, M:Metazoans, F:Fungi, and O: Others.

-fam              Provide TE family name for the input te sequence
                  Default: All
                  All: the input squence is unknown TEs
                  ClassI: the input sequence is ClassI TEs
                  ClassII: the input sequence is ClassII TEs
                  LTR: the input sequence is LTR TEs
                  nLTR: the input sequence is nLTR TEs
                  LINE: the input sequence is LINE TEs
                  SINE: the input sequence is SINE TEs
                  Domain: the input sequence is Class II TEs with specified super families

-modify           If set this argument, users need to provide domain file generated from another script: DeepTE_domain.py


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

-i                Input sequences that are unknown TEs or DNA sequences.

-s                Provide supplementary dir that contains required files.

--hmmscan         File path to hmmscan executable, Default: /usr/bin/hmmscan"


```

# Examples
DeepTE.py  
**Input data**  
Sequence data (fasta format)  
**Command**  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -m model_dir -sp P -fam ClassI  
**Output data**  
a. name list (first column: original name; second column: predicted name with DeepTE)  
b. fasta file with new predicted TE name  

DeepTE_domain.py  
**Input data**  
Sequence data (fasta format)  
**Command**  
DeepTE_domain.py -d working_dir -o output_dir -i input_seq.fasta -s supfile_dir --hmmscan   Path/to/hmmscan  
**Output data**  
a. TE domain file. (first column: orignial name; second column: domain information for different open reading frames)






# Work flows
![](https://i.imgur.com/eqINPxf.png)

**Figure** A pipeline for classifying unknown TEs based on trained models. The unknown TEs go through eight models to be classified into different families. Two correction steps are conducted during classification. In Class model, TR domain exists in predicted Class I TEs that will be corrected to Class II TEs, while RT domain exists in predicted Class II TEs that will be corrected to Class I TEs. In ClassI model, EN domain exists in predicted LTR TEs will be correct to nLTR TEs.


## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::


