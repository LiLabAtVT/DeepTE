#!/usr/bin/env python

##updation 122520 consider the prop filtration
##updation 12.14 revise script to generate O output add model_nm in the function
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


##updation 122520
def store_unknown (file,type):
    store_infor_dic = {}
    with open (file,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if col[1] == 'unknown':

                if type == 'All':
                    store_infor_dic[col[0]] = 'unknown'
                if type == 'ClassI':
                    store_infor_dic[col[0]] = 'ClassI'
                #if type == 'ClassII':
                #    store_infor_dic[col[0]] = 'ClassII'
                if type == 'nLTR':
                    store_infor_dic[col[0]] = 'ClassI_nLTR'
    return (store_infor_dic)


def store_mite_nmite_infor (file):
    store_infor_dic = {}
    with open (file,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            ##updation 122520
            if re.match('.+_(.+)',col[1]):
                mt = re.match('.+_(.+)',col[1])
                mite_type = mt.group(1)
            else:
                ##mite_type == 'unknown'
                mite_type = col[1]

            store_infor_dic[col[0]] = mite_type
    return (store_infor_dic)

##updation 122520
def modi_results (eachfl,new_annot_nm,input_opt_DeepTE_dir,fam_nm):

    store_modi_line_list = []
    with open(eachfl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if col[1] == 'unknown':
                annot_nm = new_annot_nm
            else:
                annot_nm = col[1]
            final_modi_line = col[0] + '\t' + annot_nm
            store_modi_line_list.append(final_modi_line)
    with open(input_opt_DeepTE_dir + '/' + fam_nm +'_modi_results.txt', 'w+') as opt:
        for eachline in store_modi_line_list:
            opt.write(eachline + '\n')

def extract_combine_infor (input_opt_DeepTE_dir,input_opt_combine_DeepTE_dir,sp_type):

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

    ##updation 12.14
    if sp_type != 'O':
        for eachfl in opt_fl_list:

            ##updation 122520
            ##modify the output to allow the unknown to be the relative name based on the file name

            if 'DIRS_results.txt' in eachfl:
                DIRS_results_fl = input_opt_DeepTE_dir + '/DIRS_results.txt'
                final_fl_list.append(DIRS_results_fl)

            if 'helitron_results.txt' in eachfl:
                helitron_results_fl = input_opt_DeepTE_dir + '/helitron_results.txt'
                final_fl_list.append(helitron_results_fl)

            if 'LINE_results.txt' in eachfl:
                new_annot_nm = 'ClassI_nLTR_LINE'
                fam_nm = 'LINE'
                modi_results(eachfl, new_annot_nm, input_opt_DeepTE_dir, fam_nm)
                line_results_fl = input_opt_DeepTE_dir + '/' + fam_nm + '_modi_results.txt'
                final_fl_list.append(line_results_fl)

            if '/LTR_results.txt' in eachfl:
                new_annot_nm = 'ClassI_LTR'
                fam_nm = 'LTR'
                modi_results(eachfl, new_annot_nm, input_opt_DeepTE_dir, fam_nm)
                ltr_results_fl = input_opt_DeepTE_dir + '/' + fam_nm + '_modi_results.txt'
                final_fl_list.append(ltr_results_fl)

            if 'SINE_results.txt' in eachfl:
                new_annot_nm = 'ClassI_nLTR_SINE'
                fam_nm = 'SINE'
                modi_results(eachfl, new_annot_nm, input_opt_DeepTE_dir, fam_nm)
                ltr_results_fl = input_opt_DeepTE_dir + '/' + fam_nm + '_modi_results.txt'
                final_fl_list.append(ltr_results_fl)

            if 'PLE_results.txt' in eachfl:
                ple_results_fl = input_opt_DeepTE_dir + '/PLE_results.txt'
                final_fl_list.append(ple_results_fl)
            if 'ClassII_domain_results.txt' in eachfl:
                classII_domain_results_fl = input_opt_DeepTE_dir + '/ClassII_domain_results.txt'
                final_fl_list.append(classII_domain_results_fl)


            ##updation 122520
            ##once the unknown occurs, the TE will not go to the next level, we need to consider this case
            if 'All_results.txt' in eachfl:
                All_results_fl = input_opt_DeepTE_dir + '/All_results.txt'
                store_infor_dic = store_unknown(All_results_fl, 'All')
                with open(input_opt_DeepTE_dir + '/All_unknown_results.txt', 'w+') as opt:
                    for eachnm in store_infor_dic:
                        opt.write(eachnm + '\t' + store_infor_dic[eachnm] + '\n')
                final_fl_list.append(input_opt_DeepTE_dir + '/All_unknown_results.txt')
            if 'ClassI_results.txt' in eachfl:
                ClassI_results_fl = input_opt_DeepTE_dir + '/ClassI_results.txt'
                store_infor_dic = store_unknown(ClassI_results_fl, 'ClassI')
                with open(input_opt_DeepTE_dir + '/ClassI_unknown_results.txt', 'w+') as opt:
                    for eachnm in store_infor_dic:
                        opt.write(eachnm + '\t' + store_infor_dic[eachnm] + '\n')
                final_fl_list.append(input_opt_DeepTE_dir + '/ClassI_unknown_results.txt')
            if 'nLTR_results.txt' in eachfl:
                nLTR_results_fl = input_opt_DeepTE_dir + '/nLTR_results.txt'
                store_infor_dic = store_unknown(nLTR_results_fl, 'nLTR')
                with open(input_opt_DeepTE_dir + '/nLTR_unknown_results.txt', 'w+') as opt:
                    for eachnm in store_infor_dic:
                        opt.write(eachnm + '\t' + store_infor_dic[eachnm] + '\n')
                final_fl_list.append(input_opt_DeepTE_dir + '/nLTR_unknown_results.txt')
            ##no need to store the ClassII since the domain results have helped to solve this problem.


    else:
        for eachfl in opt_fl_list:
            if 'DIRS_results.txt' in eachfl:
                DIRS_results_fl = input_opt_DeepTE_dir + '/DIRS_results.txt'
                final_fl_list.append(DIRS_results_fl)
            if 'helitron_results.txt' in eachfl:
                helitron_results_fl = input_opt_DeepTE_dir + '/helitron_results.txt'
                final_fl_list.append(helitron_results_fl)

            if '/nLTR_results.txt' in eachfl:
                new_annot_nm = 'ClassI_nLTR'
                fam_nm = 'nLTR'
                modi_results(eachfl, new_annot_nm, input_opt_DeepTE_dir, fam_nm)
                ltr_results_fl = input_opt_DeepTE_dir + '/' + fam_nm + '_modi_results.txt'
                final_fl_list.append(ltr_results_fl)

            if '/LTR_results.txt' in eachfl:
                new_annot_nm = 'ClassI_LTR'
                fam_nm = 'LTR'
                modi_results(eachfl, new_annot_nm, input_opt_DeepTE_dir, fam_nm)
                ltr_results_fl = input_opt_DeepTE_dir + '/' + fam_nm + '_modi_results.txt'
                final_fl_list.append(ltr_results_fl)

            if 'PLE_results.txt' in eachfl:
                ple_results_fl = input_opt_DeepTE_dir + '/PLE_results.txt'
                final_fl_list.append(ple_results_fl)
            if 'ClassII_domain_results.txt' in eachfl:
                classII_domain_results_fl = input_opt_DeepTE_dir + '/ClassII_domain_results.txt'
                final_fl_list.append(classII_domain_results_fl)

            ##updation 122520
            ##once the unknown occurs, the TE will not go to the next level, we need to consider this case
            if 'All_results.txt' in eachfl:
                All_results_fl = input_opt_DeepTE_dir + '/All_results.txt'
                store_infor_dic = store_unknown(All_results_fl, 'All')
                with open(input_opt_DeepTE_dir + '/All_unknown_results.txt', 'w+') as opt:
                    for eachnm in store_infor_dic:
                        opt.write(eachnm + '\t' + store_infor_dic[eachnm] + '\n')
                final_fl_list.append(input_opt_DeepTE_dir + '/All_unknown_results.txt')
            if 'ClassI_results.txt' in eachfl:
                ClassI_results_fl = input_opt_DeepTE_dir + '/ClassI_results.txt'
                store_infor_dic = store_unknown(ClassI_results_fl, 'ClassI')
                with open(input_opt_DeepTE_dir + '/ClassI_unknown_results.txt', 'w+') as opt:
                    for eachnm in store_infor_dic:
                        opt.write(eachnm + '\t' + store_infor_dic[eachnm] + '\n')
                final_fl_list.append(input_opt_DeepTE_dir + '/ClassI_unknown_results.txt')


    if len(final_fl_list) != 0:
        fl_string = ''
        for eachfl in final_fl_list:
            fl_string = fl_string + ' ' + eachfl

        cmd = 'cat ' + fl_string + ' > ' + input_opt_DeepTE_dir + '/opt_DeepTE.txt'
        subprocess.call(cmd,shell=True)

        ##updation 122520
        ##modify the unknown for the mite
        store_final_line_list = []
        with open (input_opt_DeepTE_dir + '/opt_DeepTE.txt','r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                annot = col[1]

                if annot == 'unknown_unknown':
                    new_annot = 'ClassII'
                else:
                    if re.match('unknown_.+',annot):
                        mt = re.match('unknown_(.+)',annot)
                        mite_st = mt.group(1)
                        new_annot = 'ClassII_' + mite_st
                    else:
                        new_annot = annot

                final_line = col[0] + '\t' + new_annot
                store_final_line_list.append(final_line)

        with open (input_opt_combine_DeepTE_dir + '/opt_DeepTE.txt','w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')

        ##debug
        #store_dup_te_dic = {}
        #for eachline in store_final_line_list:
        #    col = eachline.strip().split()
        #    if col[0] in store_dup_te_dic:
        #        store_dup_te_dic[col[0]] += 1
        #    else:
        #        store_dup_te_dic[col[0]] = 1
        #for eachte in store_dup_te_dic:
        #    if store_dup_te_dic[eachte] > 1:
        #        print(eachte)

    ##if only domain exit
    if len(opt_fl_list) == 3:
        if input_opt_DeepTE_dir + '/Domain_results.txt' in opt_fl_list:
            cmd = 'cp ' + input_opt_DeepTE_dir + '/Domain_results.txt' + ' ' + input_opt_combine_DeepTE_dir + '/opt_DeepTE.txt'
            print(cmd)
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






