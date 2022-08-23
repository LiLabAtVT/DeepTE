#!/usr/bin/env python

##this script will generate another pipeline that will be used to class different orders using eight models

import re
import glob
import os
import sys
import subprocess
import random
from classify_TE_keras_model_predict_kmer import train_model




input_CNN_data_file = sys.argv[1]
input_store_class_report_dir = sys.argv[2] ##this stores the reports for each model


##def a function to train each model
def train_model_pipeline (input_CNN_data_file,input_store_class_report_dir):

    ##initiate a dic to store all the score results for each model
    store_all_score_dic = {}

    train_model(input_CNN_data_file, input_store_class_report_dir, store_all_score_dic)

    return (store_all_score_dic)


store_all_score_dic = train_model_pipeline (input_CNN_data_file,input_store_class_report_dir)

##write out score dic
with open ('opt_score_all.txt','w+') as opt:
    for eachnm in store_all_score_dic:
        opt.write(eachnm + '\t' + store_all_score_dic[eachnm] + '\n')

