#!/usr/bin/env python

##updation 9.26 generate UNS model output

##this script is extract the output from output_DeepTE and combine all the information
import subprocess
import glob
from Bio import SeqIO
import re

def store_infor (file):
    store_infor_dic = {}
    with open (file,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            store_infor_dic[col[0]] = col[1]
    return (store_infor_dic)

def store_mite_nmite_infor (file):
    store_infor_dic = {}
    with open (file,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            mt = re.match('.+_(.+)',col[1])
            mite_type = mt.group(1)
            store_infor_dic[col[0]] = mite_type
    return (store_infor_dic)


def extract_combine_infor (input_opt_DeepTE_dir,input_opt_combine_DeepTE_dir):

    opt_fl_list = glob.glob(input_opt_DeepTE_dir + '/*')
    if input_opt_DeepTE_dir + '/ClassII_results.txt' in opt_fl_list and input_opt_DeepTE_dir + '/Domain_results.txt' in opt_fl_list:

        ##combine ClassII_results and domain results
        ##store the name
        store_infor_ClassII_dic = store_infor (input_opt_DeepTE_dir + '/ClassII_results.txt')
        store_infor_Domain_dic = store_mite_nmite_infor(input_opt_DeepTE_dir + '/Domain_results.txt')

        store_new_infor_dic = {}
        for eachnm in store_infor_ClassII_dic:
            new_nm = store_infor_ClassII_dic[eachnm] + '_' + store_infor_Domain_dic[eachnm]
            store_new_infor_dic[eachnm] = new_nm

        ##write out combined information file
        with open (input_opt_DeepTE_dir + '/ClassII_domain_results.txt','w+') as opt:
            for eachnm in store_new_infor_dic:
                opt.write(eachnm + '\t' + store_new_infor_dic[eachnm] + '\n')

    ##combine the file
    final_fl_list = []
    opt_fl_list = glob.glob(input_opt_DeepTE_dir + '/*')
    for eachfl in opt_fl_list:
        if 'DIRS_results.txt' in eachfl:
            DIRS_results_fl = input_opt_DeepTE_dir + '/DIRS_results.txt'
            final_fl_list.append(DIRS_results_fl)
        if 'helitron_results.txt' in eachfl:
            helitron_results_fl = input_opt_DeepTE_dir + '/helitron_results.txt'
            final_fl_list.append(helitron_results_fl)
        if 'LINE_results.txt' in eachfl:
            line_results_fl = input_opt_DeepTE_dir + '/LINE_results.txt'
            final_fl_list.append(line_results_fl)
        if '/LTR_results.txt' in eachfl:
            ltr_results_fl = input_opt_DeepTE_dir + '/LTR_results.txt'
            final_fl_list.append(ltr_results_fl)
        if 'SINE_results.txt' in eachfl:
            sine_results_fl = input_opt_DeepTE_dir + '/SINE_results.txt'
            final_fl_list.append(sine_results_fl)
        if 'PLE_results.txt' in eachfl:
            ple_results_fl = input_opt_DeepTE_dir + '/PLE_results.txt'
            final_fl_list.append(ple_results_fl)
        if 'ClassII_domain_results.txt' in eachfl:
            classII_domain_results_fl = input_opt_DeepTE_dir + '/ClassII_domain_results.txt'
            final_fl_list.append(classII_domain_results_fl)

    if len(final_fl_list) != 0:
        fl_string = ''
        for eachfl in final_fl_list:
            fl_string = fl_string + ' ' + eachfl

        cmd = 'cat ' + fl_string + ' > ' + input_opt_combine_DeepTE_dir + '/opt_DeepTE.txt'
        subprocess.call(cmd,shell=True)

    ##if only domain exit
    if len(opt_fl_list) == 1:
        if input_opt_DeepTE_dir + '/Domain_results.txt' in opt_fl_list:
            cmd = 'cp ' + opt_fl_list[0] + ' ' + input_opt_combine_DeepTE_dir + '/opt_DeepTE.txt'
            subprocess.call(cmd, shell=True)

##define a function to generate the new fasta file with new name
##the new name will be added to the older name with -
def generate_fasta (input_fasta_file,opt_DeepTE_file,input_opt_combine_DeepTE_dir):

    ##initiate a dic store the old and new name
    name_dic = {}
    with open (opt_DeepTE_file,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            name_dic[col[0]] = col[1]

    orignial_seq_dic = {}
    for seq_record in SeqIO.parse(input_fasta_file,'fasta'):
        orignial_seq_dic[seq_record.id] = str(seq_record.seq)

    new_seq_dic = {}
    for eachid in orignial_seq_dic:
        new_nm = eachid + '__' + name_dic[eachid]
        new_seq_dic[new_nm] = orignial_seq_dic[eachid]

    ##write out fasta file
    with open (input_opt_combine_DeepTE_dir + '/opt_DeepTE.fasta','w+') as opt:
        for eachid in new_seq_dic:
            opt.write('>' + eachid + '\n' + new_seq_dic[eachid] + '\n')


##updation 9.26 generate UNS model output
def generate_fasta_UNS (input_fasta_file,input_opt_DeepTE_dir,input_opt_combine_DeepTE_dir):

    cmd = 'cp ' + input_opt_DeepTE_dir + '/UNS_results.txt ' + input_opt_combine_DeepTE_dir + '/opt_DeepTE.txt'
    subprocess.call(cmd, shell=True)

    ##initiate a dic store the old and new name
    name_dic = {}
    with open (input_opt_combine_DeepTE_dir + '/opt_DeepTE.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            name_dic[col[0]] = col[1]

    orignial_seq_dic = {}
    for seq_record in SeqIO.parse(input_fasta_file,'fasta'):
        orignial_seq_dic[seq_record.id] = str(seq_record.seq)

    new_seq_dic = {}
    for eachid in orignial_seq_dic:
        new_nm = eachid + '__' + name_dic[eachid]
        new_seq_dic[new_nm] = orignial_seq_dic[eachid]

    ##write out fasta file
    with open(input_opt_combine_DeepTE_dir + '/opt_DeepTE.fasta', 'w+') as opt:
        for eachid in new_seq_dic:
            opt.write('>' + eachid + '\n' + new_seq_dic[eachid] + '\n')






