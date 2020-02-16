#!/usr/bin/env python

##this script is to generate the kmer vector
import itertools


##word_seq generates eg. ['AA', 'AT', 'TC', 'CG', 'GT']
def word_seq(seq, k, stride=1):
    i = 0
    words_list = []
    while i <= len(seq) - k:
        words_list.append(seq[i: i + k])
        i += stride
    return (words_list)


##generate all the combinations of ATCG, we will input the k-mer number
def generate_kmer_dic (repeat_num):

    ##initiate a dic to store the kmer dic
    ##kmer_dic = {'ATC':0,'TTC':1,...}
    kmer_dic = {}

    bases = ['A','G','C','T']
    kmer_list = list(itertools.product(bases, repeat=int(repeat_num)))
    for eachitem in kmer_list:
        #print(eachitem)
        each_kmer = ''.join(eachitem)
        kmer_dic[each_kmer] = 0

    return (kmer_dic)

##add the number into the kmer in the kmer dic
##generate the vector from the kmer_dic
##this will generat_mat for only one sample
def generate_mat (words_list,kmer_dic):
    for eachword in words_list:
        kmer_dic[eachword] += 1

    num_list = []  ##this dic stores num_dic = [0,1,1,0,3,4,5,8,2...]
    for eachkmer in kmer_dic:
        num_list.append(kmer_dic[eachkmer])

    return (num_list)

##generate matrix for all samples
def generate_mats (seqs):

    seq_mats = []
    for eachseq in seqs:
        words_list = word_seq(eachseq, 7, stride=1)  ##change the k to 3
        kmer_dic = generate_kmer_dic(7)  ##this number should be the same as the window slide number
        num_list = generate_mat(words_list,kmer_dic)

        ##store the all the samples into seq_mats
        ##seq_mats = [[0,1,3,4],[3,4,5,6],...]
        seq_mats.append(num_list)

    return (seq_mats)


##generate the COV labels containing four kinds of labels
def conv_labels(labels,input_data_nm):
    converted = []
    for label in labels:



        if input_data_nm == 'All':
            if 'ClassI_' in label:
                converted.append(0)
            elif 'ClassII_' in label:
                converted.append(1)
            elif 'ClassIII_' in label:
                converted.append(2)


        if input_data_nm == 'ClassI':
            if '_LTR_' in label:
                converted.append(0)
            elif '_nLTR_' in label:
                converted.append(1)

        if input_data_nm =='LTR':
            if 'LTR_Copia' in label:
                converted.append(0)
            elif 'LTR_Gypsy' in label:
                converted.append(1)

        if input_data_nm == 'nLTR':
            if 'nLTR_LINE' in label:
                converted.append(0)
            elif 'nLTR_SINE' in label:
                converted.append(1)
            elif 'nLTR_DIRS' in label:
                converted.append(2)
            elif 'nLTR_PLE' in label:
                converted.append(3)

        if input_data_nm == 'LINE':
            if 'LINE_L1' in label:
                converted.append(0)
            elif 'LINE_I' in label:
                converted.append(1)

        if input_data_nm == 'SINE':
            if 'SINE_tRNA' in label:
                converted.append(0)
            elif 'SINE_7SL' in label:
                converted.append(1)

        if input_data_nm == 'ClassII':
            if label == 'ClassII_DNA_MITE_TcMar' or label == 'ClassII_DNA_nMITE_TcMar':
                converted.append(0)
            elif label == 'ClassII_DNA_MITE_hAT' or label == 'ClassII_DNA_nMITE_hAT':
                converted.append(1)
            elif label == 'ClassII_DNA_MITE_Mutator' or label == 'ClassII_DNA_nMITE_Mutator':
                converted.append(2)
            elif label == 'ClassII_DNA_MITE_P' or label == 'ClassII_DNA_nMITE_P':
                converted.append(3)
            elif label == 'ClassII_DNA_nMITE_Harbinger' or label == 'ClassII_DNA_MITE_Harbinger':
                converted.append(4)
            elif label == 'ClassII_DNA_MITE_CACTA' or label == 'ClassII_DNA_nMITE_CACTA':
                converted.append(5)

        if input_data_nm == 'Domain':
            if '_MITE_' in label:
                converted.append(0)
            elif '_nMITE_' in label:
                converted.append(1)


    return converted

