#!/usr/bin/env python

##import modules
from Bio import SeqIO


def change_format_for_ncc (input_ori_seq_file):
    ##initiate a dic to store output file
    final_format_dic = {}
    seq_count = 0
    for seq_record in SeqIO.parse(input_ori_seq_file,'fasta'):
        seq_count += 1
        label = seq_record.id
        final_format_dic[str(seq_count)] = {'label':label,'seq':str(seq_record.seq)}

    return (final_format_dic)

def generate_target_line (final_format_dic):

    final_format_line_dic = {}
    seq_count = 0
    for eachid in final_format_dic:
        seq_count += 1
        final_line = final_format_dic[eachid]['label'] + ',' + final_format_dic[eachid]['seq']
        final_format_line_dic[str(seq_count)] = final_line

    return (final_format_line_dic)



