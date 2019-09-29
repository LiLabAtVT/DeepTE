#!/usr/bin/env python

##this script is to wrap all the seven models
##this script will directly test all the input sequence
##do not classify the sequence into right and wrong
##later add domain detection for the sequence


import re
import glob
import numpy as np
from keras.models import load_model

#from seq_reader_kmer import load_data        # parsing data file
#from one_hot_rep_kmer import generate_mats   # converting to correct format

from scripts import DeepTE_seq_reader_kmer
from scripts import DeepTE_one_hot_rep_kmer


def generate_input_data_without_load_data (x):
    ##for the training dataset
    #X, y = load_data(input_dataset)  # sequences, labels
    # print(X)            ##Note: X is the ['AA','BB','CC]
    # print(y)            ##Note: y is the ['C1','C2','C3']
    X = DeepTE_one_hot_rep_kmer.generate_mats(x)  # convert to array of representation matrices
    # convert to integer labels
    #y = conv_labels(y, input_data_nm,name_number_dic,predic_list)
    # work with np arrays
    X = np.asarray(X)
    #Y = np.asarray(y)

    #return (X,Y)
    return (X)


##generate name_number_dic
##this will be helpful for the predicted classes identification
##different species types have different model types
def generate_name_number_dic (model_nm):

    name_number_dic = {}
    name_number_dic[model_nm] = {}

    #if input_spe_type == 'M':
    if model_nm == 'UNS':
        name_number_dic[model_nm][str(0)] = 'TE'
        name_number_dic[model_nm][str(1)] = 'CDS'
        name_number_dic[model_nm][str(2)] = 'INT'

    return (name_number_dic)


def predict_te (model,model_nm,x_test_list,y_test_nm_list):

    ##x_test_list: eg. ['TT','AA','CC']
    ##y_test_nm_list contains the name information
    ##y_test_nm_list: eg. ['C1','C2','C3']

    ##do not consider the right and wrong
    #x_right_list = []
    #y_right_list = []
    #y_right_nm_list = [] ##create this name list is to be classified based on other name

    x_new_list = []
    y_new_nm_list = []

    ##initiate a new dic to store the results for this prediction that could be used to calculate the parameters for each model
    #store_right_results_dic = {}
    ##initiate a new dic to store the wrong results
    #store_wrong_results_dic = {}

    store_results_dic = {}


    #################################
    ##step 1: generate the input data
    model = load_model(model)
    ##X_test is the input data for the prediction of the model
    X_test = generate_input_data_without_load_data(x_test_list)  ##different model_nm generates different class code


    #print('the name_number_dic is ')
    #print(name_number_dic)


    ##change X_test vector
    X_test = X_test.reshape(X_test.shape[0], 1, 16384, 1)  ##kmer == 3 so it would be 64
    X_test = X_test.astype('float64')


    ####################################
    ##step 2: generate the predict class
    predicted_classes = model.predict(X_test)
    predicted_classes = np.argmax(np.round(predicted_classes), axis=1)

    #######################################################################################
    ##step 3: extract rigth predicted order that will be used for the next prediction round
    predicted_classes_list = predicted_classes.tolist()

    #y_test_list = Y_test.tolist()  ##y_test_list contains [0,1,2,1,2,3]

    name_number_dic = generate_name_number_dic(model_nm)

    for i in range(0, len(predicted_classes_list)):
        x_new_list.append(x_test_list[i])
        y_new_nm_list.append(y_test_nm_list[i])
        store_results_dic[str(i)] = str(y_test_nm_list[i])  + '\t' + name_number_dic[model_nm][str(predicted_classes_list[i])]

    return (x_new_list,y_new_nm_list,store_results_dic,predicted_classes_list)



def classify_pipeline (input_model_dir,input_dataset,input_store_predict_dir):

    ##Note:
    #list_model_name = ['All','ClassI','LTR','noLTR','DNA','MITE','noMITE']
    ##rihgt_nm_list for y is the ['C1','C2','C3']
    ##right_list for y is the ['0','1','2']

    ##initate a dic to store all the information including wrong and right results
    #store_all_information_dic = {}

    ##store the model direction
    model_file_dic = {}
    model_fl_list = glob.glob(input_model_dir + '/*')
    for eachmodel_path in model_fl_list:
        ##get model name
        mt = re.match('.+/(.+)_model.h5', eachmodel_path)
        model_nm = mt.group(1)
        model_file_dic[model_nm] = eachmodel_path

    x_all_test_list, y_all_test_nm_list = DeepTE_seq_reader_kmer.load_data(input_dataset)


    ######################
    ##detect TE in the All
    model_name = 'UNS'
    x_all_right_list, y_all_right_nm_list, store_all_results_dic,predicted_classes_list = \
        predict_te(model_file_dic[model_name], model_name, x_all_test_list, y_all_test_nm_list)

    with open(input_store_predict_dir + '/' + model_name + '_results.txt', 'w+') as opt:
        for eachid in store_all_results_dic:
            # store_all_information_dic[new_id_name] = store_all_results_dic[eachid]
            opt.write(store_all_results_dic[eachid] + '\n')









