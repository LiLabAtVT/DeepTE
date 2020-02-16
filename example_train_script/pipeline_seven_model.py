#!/usr/bin/env python

##this script will generate another pipeline that will be used to class different orders using eight models

import re
import glob
import os
import sys
import subprocess
import random
from classify_TE_keras_model_predict_kmer import train_model


input_CNN_train_data_file = sys.argv[1]
input_CNN_test_data_file = sys.argv[2]
input_CNN_train_data_dir = sys.argv[3]
input_CNN_test_data_dir = sys.argv[4]

#input_CNN_data_dir = sys.argv[2] ##this stores CNN data from each data set
#input_CNN_data_shuff_dir = sys.argv[3] ##this stores shuffle CNN data from each data set
input_store_class_report_dir = sys.argv[5] ##this stores the reports for each model


##the shuffle line file will not be used
##define a function to shuffle all lines in a file
def shuffle_line_file (input_file,store_dir):

    line_dic = {}

    count = 0
    with open (input_file,'r') as ipt:
        for eachline in ipt:
            count += 1
            eachline = eachline.strip('\n')
            line_dic[str(count)] = eachline

    line_id_list = list(line_dic.keys())

    random.shuffle(line_id_list)

    ##get the input_file name
    input_file_nm = ''
    if '/' in input_file:
        mt = re.match('.+\/(.+)',input_file)
        input_file_nm = mt.group(1)

    with open (store_dir + '/shuffle_' + input_file_nm,'w+') as opt:
        for eachid in line_id_list:
            opt.write(line_dic[eachid] + '\n')

    #new_file_dataset = store_dir + '/shuffle_' + input_file_nm
    #return (new_file_dataset)


##define a function to classify dataset into eight models that will be trained.
def divide_dataset (input_CNN_data_file,input_CNN_data_dir):

    list_divide_dataset_list = ['All','ClassI','LTR','nLTR','LINE','SINE','ClassII','Domain']

    for eachitem in list_divide_dataset_list:

        ##initiate a dic to store the final output line of the dataset
        divide_CNN_data_dic = {}

        line_count = 0
        with open(input_CNN_data_file, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split(',')
                fam_nm = col[0]

                line_count += 1

                if eachitem == 'All': ##store all the lines
                    divide_CNN_data_dic[str(line_count)] = eachline

                if eachitem == 'ClassI': ##store classI
                    if 'ClassI_' in fam_nm:
                        divide_CNN_data_dic[str(line_count)] = eachline

                if eachitem == 'LTR': ##store LTR
                    if '_LTR_' in fam_nm:
                        divide_CNN_data_dic[str(line_count)] = eachline

                if eachitem == 'nLTR': ##store noLTR
                    if '_nLTR_' in fam_nm:
                        divide_CNN_data_dic[str(line_count)] = eachline

                if eachitem == 'LINE': ##store noLTR
                    if 'nLTR_LINE' in fam_nm:
                        divide_CNN_data_dic[str(line_count)] = eachline

                if eachitem == 'SINE': ##store noLTR
                    if 'nLTR_SINE' in fam_nm:
                        divide_CNN_data_dic[str(line_count)] = eachline

                if eachitem == 'ClassII':
                    if 'ClassII_' in fam_nm:
                        divide_CNN_data_dic[str(line_count)] = eachline

                if eachitem == 'Domain':
                    if 'MITE' in fam_nm:
                        divide_CNN_data_dic[str(line_count)] = eachline


        ##write to a file
        with open (input_CNN_data_dir + '/' + eachitem + '_CNN_data.txt','w+') as opt:
            for eachid in divide_CNN_data_dic:
                opt.write(divide_CNN_data_dic[eachid] + '\n')

    return (list_divide_dataset_list)


def shuff_all_dataset (list_divide_train_dataset_list,input_CNN_train_data_dir,input_CNN_test_data_dir):

    ##mkdir for the shuffle dir
    input_CNN_shuffle_train_data_dir = 'input_CNN_shuffle_train_data_dir'
    if not os.path.exists(input_CNN_shuffle_train_data_dir):
        os.makedirs(input_CNN_shuffle_train_data_dir)

    input_CNN_shuffle_test_data_dir = 'input_CNN_shuffle_test_data_dir'
    if not os.path.exists(input_CNN_shuffle_test_data_dir):
        os.makedirs(input_CNN_shuffle_test_data_dir)

    ##shuffle file
    for eachnm in list_divide_train_dataset_list:
        shuffle_line_file(input_CNN_train_data_dir + '/' + eachnm + '_CNN_data.txt', input_CNN_shuffle_train_data_dir)
        shuffle_line_file(input_CNN_test_data_dir + '/' + eachnm + '_CNN_data.txt', input_CNN_shuffle_test_data_dir)

    return (input_CNN_shuffle_train_data_dir,input_CNN_shuffle_test_data_dir)



##def a function to train each model
def train_model_pipeline (list_divide_train_dataset_list,input_CNN_shuffle_train_data_dir,input_CNN_shuffle_test_data_dir,input_store_class_report_dir):

    ##initiate a dic to store all the score results for each model
    store_all_score_dic = {}

    for eachnm in list_divide_train_dataset_list:

        ##shuffle the training and test data
        train_dataset = input_CNN_shuffle_train_data_dir + '/shuffle_' + eachnm + '_CNN_data.txt'

        test_dataset = input_CNN_shuffle_test_data_dir + '/shuffle_' + eachnm + '_CNN_data.txt'

        train_model(train_dataset, test_dataset, input_store_class_report_dir, store_all_score_dic)

    return (store_all_score_dic)


list_divide_train_dataset_list = divide_dataset (input_CNN_train_data_file,input_CNN_train_data_dir)
list_divide_test_dataset_list = divide_dataset (input_CNN_test_data_file,input_CNN_test_data_dir)
input_CNN_shuffle_train_data_dir,input_CNN_shuffle_test_data_dir = shuff_all_dataset (list_divide_train_dataset_list,input_CNN_train_data_dir,input_CNN_test_data_dir)


store_all_score_dic = train_model_pipeline (list_divide_train_dataset_list,input_CNN_shuffle_train_data_dir,input_CNN_shuffle_test_data_dir
                                            ,input_store_class_report_dir)

##write out score dic
with open ('opt_score_all.txt','w+') as opt:
    for eachnm in store_all_score_dic:
        opt.write(eachnm + '\t' + store_all_score_dic[eachnm] + '\n')

