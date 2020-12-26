

DeepTE
===
DeepTE is aimed to classify transposons with unknown classification *via* Convolutional Neural Network.

# Updating
**12/25/2020**  
Add a proability threshold to annotate TEs. For example, a TE has a probability (0.6) to be ClassI, If users set 0.7 as the threshold, this TE will be labeled as 'unknown', Default: 0.6.

**08/15/2020**  
In the 'store_temp_opt_dir' of the working_dir, we add information about probability of each family the input TE belongs to.   

# Introduction
Transposable elements (TEs) classification is an essential step decoding their roles in a genome. With reference genomes from non-model species available, it has begun to overstep efforts to annotate TEs, and more tools are needed to efficiently handle the emerged sequence information. We developed a novel tool, DeepTE, which classifies unknown TEs on basis of convolutional neural network. DeepTE utilized co-occurrence of k-mers towards TE sequences as input vector, and seven k-mer size was testified to be suitable for the classification. Eight models have been trained for different TE classification purposes. DeepTE applied domains from TEs to correct false classification. An additional model was also trained to distinguish between non-TEs and TEs targeting plant species. 

Given exclusive TEs of different species type, it can therefore classify seven orders, and 11-24 superframilies towards Plants, Metazoans, Fungi, and Others. This tool successfully leverages convolutional neural network to TE classification, assisting to precisely identify and annotate TEs in a new sequenced eukaryotic genome.

# Dependence and requirements
DeepTE is developed in Python with modules and external tools.

Before running this pipeline, a dependency check should be performed first to make sure every dependency is correctly installed.

For information about installing the dependencies, please see below. The version numbers listed below represents the version this pipeline is developed with, and using the newest version is recommended.

## Requirements
### Use conda to install required packages (Recommend)
Install **conda**: https://www.anaconda.com/products/individual  
conda create -n py36 python=3.6  
conda activate py36  
conda install tensorflow-gpu=1.14.0  
conda install biopython  
conda install keras=2.2.4  
conda install numpy=1.16.0  
### Use pip to install required packages
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
Download the model dir from the cyVerse link  
Plants:  
https://de.cyverse.org/dl/d/89D2FE7A-41BA-4F64-80E2-B9C26D49E99F/Plants_model.tar.gz  
Metazoans:  
https://de.cyverse.org/dl/d/441459EF-6DDD-41A5-A9AB-1D5D13049F18/Metazoans_model.tar.gz  
Fungi:  
https://de.cyverse.org/dl/d/8B112733-063A-4DE9-89EC-22A062D8807B/Fungi_model.tar.gz  
Others:
https://de.cyverse.org/dl/d/34CF8ACB-0B1F-4210-8359-366A70539F01/Others_model.tar.gz  
UNS models:
https://de.cyverse.org/dl/d/3280369B-030A-4ADF-8B6F-EDD4EC21DC4A/UNS_model.tar.gz  

Download the model dir from the google link  
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
                         [-fam te_fam][-UNS yes][-prop_thr value]

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
                  ClassI: the input sequence is ClassI TEs
                  ClassII: the input sequence is ClassII subclass1 TEs
                  LTR: the input sequence is LTR TEs
                  nLTR: the input sequence is nLTR TEs
                  LINE: the input sequence is LINE TEs
                  SINE: the input sequence is SINE TEs
                  Domain: the input sequence is Class II subclass1 TEs with specified super families
                  If users do not initiate '-fam', DeepTE will regard your input sequences are unknown TEs.

-modify           If set this argument, users need to provide domain file generated from another script: DeepTE_domain.py.

-UNS              If set this argument, users need change the -i to the the DNA sequences; 
                  This function will classify the sequences into TEs, CDS, or Intergenic sequences; -sp and -fam do not need to provide.
                  Note: this model is used for plants rather than metazoans and fungi.

-prop_thr         Specify a probability threshold to annotate TE.
                  For example: a TE has a probability (0.6) to be ClassI.
                  If users set 0.7 as the threshold, 
                  this TE will be labeled as 'unknown', Default: 0.6.


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
b. store_temp_opt_dir (a directory contains prediction results for each TE group; b this directory also contains probability of each family the input TE belongs to)
c. download_X_model_dir (store downloaded models. X represents P, M, F, O, or U)

Output directory  
a. opt_DeepTE.txt (a txt file with two columns. first column: original name; second column: predicted name with DeepTE)  
b. opt_DeepTE.fasta (a fasta file with new predicted TE name)  

**Command**
- [ ] Classify unknown TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P  
**Or**  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m_dir Plants_model/

- [ ] Classify Class I TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam ClassI  
**Or**  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m_dir Plants_model/ -fam ClassI

- [ ] Classify Class II subclass1 TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam ClassII  
**Or**  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m_dir Plants_model/ -fam ClassII  

- [ ] Classify LTR TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam LTR  
**Or**  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m_dir Plants_model/ -fam LTR

- [ ] Classify nLTR TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam nLTR  
**Or**  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m_dir Plants_model/ -fam nLTR

- [ ] Classify LINE TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam LINE  
**Or**  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m_dir Plants_model/ -fam LINE

- [ ] Classify SINE TEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam SINE  
**Or**  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m_dir Plants_model/ -fam SINE

- [ ] Classify TEs into MITEs and nMITEs  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m P -fam Domain  
**Or**  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -sp P -m_dir Plants_model/ -fam Domain  

- [ ] Classify Unknown sequences into TEs, Coding sequences, or Intergenic sequences  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -m U -UNS yes  
**Or**  
DeepTE.py -d working_dir -o output_dir -i input_seq.fasta -m_dir UNS_model/ -UNS yes  

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

## Citation
Haidong Yan, Aureliano Bombarely, Song Li 2020 DeepTE: a computational method for de novo classification of transposons with convolutional neural network. Bioinformatics, Volume 36, Issue 15, 1 August 2020, Pages 4269–4275.

## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::


