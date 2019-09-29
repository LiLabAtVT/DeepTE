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




