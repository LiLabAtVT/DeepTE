#!/usr/bin/env python

##updation 122520 provide prop selection for the results

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
def generate_name_number_dic (model_nm,input_spe_type):

    name_number_dic = {}
    name_number_dic[model_nm] = {}


    if input_spe_type == 'M':
        if model_nm == 'All':
            name_number_dic[model_nm][str(0)] = 'ClassI'
            name_number_dic[model_nm][str(1)] = 'ClassII'
            name_number_dic[model_nm][str(2)] = 'ClassIII'
        if model_nm == 'ClassI':
            name_number_dic[model_nm][str(0)] = 'ClassI_LTR'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR'
        if model_nm == 'LTR':
            name_number_dic[model_nm][str(0)] = 'ClassI_LTR_Copia'
            name_number_dic[model_nm][str(1)] = 'ClassI_LTR_Gypsy'
            name_number_dic[model_nm][str(2)] = 'ClassI_LTR_ERV'
            name_number_dic[model_nm][str(3)] = 'ClassI_LTR_BEL'
        if model_nm == 'nLTR':
            name_number_dic[model_nm][str(0)] = 'ClassI_nLTR_LINE'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR_SINE'
            name_number_dic[model_nm][str(2)] = 'ClassI_nLTR_DIRS'
            name_number_dic[model_nm][str(3)] = 'ClassI_nLTR_PLE'
        if model_nm == 'LINE':
            name_number_dic[model_nm][str(0)] = 'ClassI_nLTR_LINE_R2'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR_LINE_RTE'
            name_number_dic[model_nm][str(2)] = 'ClassI_nLTR_LINE_Jockey'
            name_number_dic[model_nm][str(3)] = 'ClassI_nLTR_LINE_L1'
            name_number_dic[model_nm][str(4)] = 'ClassI_nLTR_LINE_I'
        if model_nm == 'SINE':
            name_number_dic[model_nm][str(0)] = 'ClassI_nLTR_SINE_tRNA'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR_SINE_7SL'
            name_number_dic[model_nm][str(2)] = 'ClassI_nLTR_SINE_5S'
        if model_nm == 'ClassII':
            name_number_dic[model_nm][str(0)] = 'ClassII_DNA_TcMar'
            name_number_dic[model_nm][str(1)] = 'ClassII_DNA_hAT'
            name_number_dic[model_nm][str(2)] = 'ClassII_DNA_Mutator'
            name_number_dic[model_nm][str(3)] = 'ClassII_DNA_Merlin'
            name_number_dic[model_nm][str(4)] = 'ClassII_DNA_Transib'
            name_number_dic[model_nm][str(5)] = 'ClassII_DNA_P'
            name_number_dic[model_nm][str(6)] = 'ClassII_DNA_PiggyBac'
            name_number_dic[model_nm][str(7)] = 'ClassII_DNA_Harbinger'
            name_number_dic[model_nm][str(8)] = 'ClassII_DNA_CACTA'
        if model_nm == 'Domain':
            name_number_dic[model_nm][str(0)] = 'ClassII_DNA_MITE'
            name_number_dic[model_nm][str(1)] = 'ClassII_DNA_nMITE'

    if input_spe_type == 'P':
        if model_nm == 'All':
            name_number_dic[model_nm][str(0)] = 'ClassI'
            name_number_dic[model_nm][str(1)] = 'ClassII'
            name_number_dic[model_nm][str(2)] = 'ClassIII'
        if model_nm == 'ClassI':
            name_number_dic[model_nm][str(0)] = 'ClassI_LTR'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR'
        if model_nm == 'LTR':
            name_number_dic[model_nm][str(0)] = 'ClassI_LTR_Copia'
            name_number_dic[model_nm][str(1)] = 'ClassI_LTR_Gypsy'
        if model_nm == 'nLTR':
            name_number_dic[model_nm][str(0)] = 'ClassI_nLTR_LINE'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR_SINE'
            name_number_dic[model_nm][str(2)] = 'ClassI_nLTR_DIRS'
            name_number_dic[model_nm][str(3)] = 'ClassI_nLTR_PLE'
        if model_nm == 'LINE':
            name_number_dic[model_nm][str(0)] = 'ClassI_nLTR_LINE_L1'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR_LINE_I'
        if model_nm == 'SINE':
            name_number_dic[model_nm][str(0)] = 'ClassI_nLTR_SINE_tRNA'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR_SINE_7SL'
        if model_nm == 'ClassII':
            name_number_dic[model_nm][str(0)] = 'ClassII_DNA_TcMar'
            name_number_dic[model_nm][str(1)] = 'ClassII_DNA_hAT'
            name_number_dic[model_nm][str(2)] = 'ClassII_DNA_Mutator'
            name_number_dic[model_nm][str(3)] = 'ClassII_DNA_P'
            name_number_dic[model_nm][str(4)] = 'ClassII_DNA_Harbinger'
            name_number_dic[model_nm][str(5)] = 'ClassII_DNA_CACTA'
        if model_nm == 'Domain':
            name_number_dic[model_nm][str(0)] = 'ClassII_DNA_MITE'
            name_number_dic[model_nm][str(1)] = 'ClassII_DNA_nMITE'

    if input_spe_type == 'F':
        if model_nm == 'All':
            name_number_dic[model_nm][str(0)] = 'ClassI'
            name_number_dic[model_nm][str(1)] = 'ClassII'
            name_number_dic[model_nm][str(2)] = 'ClassIII'
        if model_nm == 'ClassI':
            name_number_dic[model_nm][str(0)] = 'ClassI_LTR'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR'
        if model_nm == 'LTR':
            name_number_dic[model_nm][str(0)] = 'ClassI_LTR_Copia'
            name_number_dic[model_nm][str(1)] = 'ClassI_LTR_Gypsy'
        if model_nm == 'nLTR':
            name_number_dic[model_nm][str(0)] = 'ClassI_nLTR_LINE'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR_SINE'
            name_number_dic[model_nm][str(2)] = 'ClassI_nLTR_DIRS'
            name_number_dic[model_nm][str(3)] = 'ClassI_nLTR_PLE'
        if model_nm == 'LINE':
            name_number_dic[model_nm][str(0)] = 'ClassI_nLTR_LINE_L1'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR_LINE_I'
        if model_nm == 'SINE':
            name_number_dic[model_nm][str(0)] = 'ClassI_nLTR_SINE_tRNA'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR_SINE_7SL'
        if model_nm == 'ClassII':
            name_number_dic[model_nm][str(0)] = 'ClassII_DNA_TcMar'
            name_number_dic[model_nm][str(1)] = 'ClassII_DNA_hAT'
            name_number_dic[model_nm][str(2)] = 'ClassII_DNA_Mutator'
            name_number_dic[model_nm][str(3)] = 'ClassII_DNA_Transib'
            name_number_dic[model_nm][str(4)] = 'ClassII_DNA_Harbinger'
            name_number_dic[model_nm][str(5)] = 'ClassII_DNA_CACTA'
            name_number_dic[model_nm][str(6)] = 'ClassII_DNA_Crypton'
        if model_nm == 'Domain':
            name_number_dic[model_nm][str(0)] = 'ClassII_DNA_MITE'
            name_number_dic[model_nm][str(1)] = 'ClassII_DNA_nMITE'

    if input_spe_type == 'O':
        if model_nm == 'All':
            name_number_dic[model_nm][str(0)] = 'ClassI'
            name_number_dic[model_nm][str(1)] = 'ClassII'
        if model_nm == 'ClassI':
            name_number_dic[model_nm][str(0)] = 'ClassI_LTR'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR'
        if model_nm == 'LTR':
            name_number_dic[model_nm][str(0)] = 'ClassI_LTR_Copia'
            name_number_dic[model_nm][str(1)] = 'ClassI_LTR_Gypsy'
        if model_nm == 'nLTR':
            name_number_dic[model_nm][str(0)] = 'ClassI_nLTR_LINE'
            name_number_dic[model_nm][str(1)] = 'ClassI_nLTR_SINE'
            name_number_dic[model_nm][str(2)] = 'ClassI_nLTR_DIRS'
            name_number_dic[model_nm][str(3)] = 'ClassI_nLTR_PLE'
        if model_nm == 'ClassII':
            name_number_dic[model_nm][str(0)] = 'ClassII_DNA_TcMar'
            name_number_dic[model_nm][str(1)] = 'ClassII_DNA_hAT'
            name_number_dic[model_nm][str(2)] = 'ClassII_DNA_Mutator'
            name_number_dic[model_nm][str(3)] = 'ClassII_DNA_Merlin'
            name_number_dic[model_nm][str(4)] = 'ClassII_DNA_PiggyBac'
            name_number_dic[model_nm][str(5)] = 'ClassII_DNA_Harbinger'
        if model_nm == 'Domain':
            name_number_dic[model_nm][str(0)] = 'ClassII_DNA_MITE'
            name_number_dic[model_nm][str(1)] = 'ClassII_DNA_nMITE'

    return (name_number_dic)


def predict_te (model,model_nm,x_test_list,y_test_nm_list,input_spe_type,y_all_test_nm_list,prop_thr):

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

    ##change X_test vector
    X_test = X_test.reshape(X_test.shape[0], 1, 16384, 1)  ##kmer == 3 so it would be 64
    X_test = X_test.astype('float64')


    ####################################
    ##step 2: generate the predict class
    Y_pred_keras = model.predict(X_test)
    predicted_classes = np.argmax(np.round(Y_pred_keras), axis=1)

    #######################################################################################
    ##step 3: extract right predicted order that will be used for the next prediction round
    predicted_classes_list = predicted_classes.tolist()
    # y_test_list = Y_test.tolist()  ##y_test_list contains [0,1,2,1,2,3]

    ##updation 122520 transfer the prop less than a threshold to be unknown for a class
    max_value_predicted_classes = np.amax(Y_pred_keras, axis=1)
    order = -1
    ls_thr_order_list = []
    for i in range(len(max_value_predicted_classes)):
        order += 1
        if max_value_predicted_classes[i] < float(prop_thr):
            ls_thr_order_list.append(order)

    new_predicted_classes_list = []
    order = -1
    for i in range(len(predicted_classes)):
        order += 1
        if order in ls_thr_order_list:
            new_class = 'unknown'
        else:
            new_class = predicted_classes[i]
        new_predicted_classes_list.append(new_class)


    name_number_dic = generate_name_number_dic(model_nm,input_spe_type)

    ##updation 122520
    for i in range(0, len(new_predicted_classes_list)):
        x_new_list.append(x_test_list[i])
        y_new_nm_list.append(y_test_nm_list[i])

        predicted_class = new_predicted_classes_list[i]
        if predicted_class != 'unknown':
            store_results_dic[str(i)] = str(y_test_nm_list[i])  + '\t' + name_number_dic[model_nm][str(new_predicted_classes_list[i])]
        else:
            store_results_dic[str(i)] = str(y_test_nm_list[i]) + '\t' + 'unknown'


    #for i in range(0, len(predicted_classes_list)):
    #    x_new_list.append(x_test_list[i])
    #    y_new_nm_list.append(y_test_nm_list[i])
    #    store_results_dic[str(i)] = str(y_test_nm_list[i])  + '\t' + name_number_dic[model_nm][str(predicted_classes_list[i])]


    ########################################################
    ##step 4: generate probability for prediction of each TE
    ##load the TE name
    store_prob_line_list = []
    first_line = 'TE_name'
    for i in range(len(list(name_number_dic[model_nm].keys()))):
        first_line = first_line + '\t' + name_number_dic[model_nm][str(i)]
    store_prob_line_list.append(first_line)

    for i in range(Y_pred_keras.shape[0]):
        prob_line = y_all_test_nm_list[i]
        for j in range(len(list(name_number_dic[model_nm].keys()))):
            prob_line = prob_line + '\t' + str(Y_pred_keras[i, j])
        store_prob_line_list.append(prob_line)

    ##updation 122520
    ##replace the predicted_classes_list to the new_predicted_classes_list
    return (x_new_list,y_new_nm_list,store_results_dic,new_predicted_classes_list,store_prob_line_list)




def classify_pipeline (input_model_dir,input_dataset,input_store_predict_dir,input_spe_type,fam_nm,prop_thr):

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

    ##select out classI to analyze
    x_classI_test_list = []
    y_classI_test_nm_list = []

    x_LTR_test_list = []
    y_LTR_test_nm_list = []

    x_nLTR_test_list = []
    y_nLTR_test_nm_list = []

    x_LINE_test_list = []
    y_LINE_test_nm_list = []

    x_SINE_test_list = []
    y_SINE_test_nm_list = []

    x_classII_test_list = []
    y_classII_test_nm_list = []

    store_helitron_results_dic = {}
    store_DIRS_results_dic = {}
    store_PLE_results_dic = {}

    ######################
    ##detect TE in the All
    model_name = 'All'
    x_all_right_list, y_all_right_nm_list, store_all_results_dic,predicted_classes_list,store_prob_line_list = \
        predict_te(model_file_dic[model_name], model_name, x_all_test_list, y_all_test_nm_list,input_spe_type,y_all_test_nm_list,prop_thr)

    ##updation 081520
    with open(input_store_predict_dir + '/' + model_name + '_probability_results.txt', 'w+') as opt:
        for eachline in store_prob_line_list:
            opt.write(eachline + '\n')

    ##y_all_right_list = [0,1,2,2,1,0,1]  has three items
    #print('the y_all_right_list len is ')
    #print(str(len(y_all_right_list)))

    count_classI = 0

    for i in range(len(y_all_right_nm_list)):
        if predicted_classes_list[i] == 0:
            count_classI += 1
            x_classI_test_list.append(x_all_right_list[i])
            y_classI_test_nm_list.append(y_all_right_nm_list[i])

        if predicted_classes_list[i] == 1:
            x_classII_test_list.append(x_all_right_list[i])
            y_classII_test_nm_list.append(y_all_right_nm_list[i])

        ##since O has no information for Helitron classification
        if input_spe_type != 'O':
            if predicted_classes_list[i] == 2:
                store_helitron_results_dic[str(i)] = y_all_right_nm_list[i] + '\t' + 'ClassIII_Helitron'

    ##since O has no information for Helitron classification
    if input_spe_type != 'O':
        ##write new results
        if fam_nm == 'All':
            with open(input_store_predict_dir + '/helitron_results.txt', 'w+') as opt:
                for eachid in store_helitron_results_dic:
                    #store_all_information_dic[new_id_name] = store_helitron_results_dic[eachid]
                    opt.write(store_helitron_results_dic[eachid] + '\n')

    ##write out All results
    if fam_nm == 'All':
        with open (input_store_predict_dir + '/' + model_name + '_results.txt','w+') as opt:
            for eachid in store_all_results_dic:
                #store_all_information_dic[new_id_name] = store_all_results_dic[eachid]
                opt.write(store_all_results_dic[eachid] + '\n')



    ##########################
    ##detect TE in the class I
    model_name = 'ClassI'

    if fam_nm == 'All' or fam_nm == 'ClassI':

        if fam_nm == 'ClassI':
            x_classI_ipt_test_list, y_classI_ipt_test_nm_list = DeepTE_seq_reader_kmer.load_data(input_dataset)
        else:
            x_classI_ipt_test_list = x_classI_test_list
            y_classI_ipt_test_nm_list = y_classI_test_nm_list

        ##updation
        if x_classI_ipt_test_list != []:

            x_classI_right_list, y_classI_right_nm_list, store_classI_results_dic,predicted_classes_list,store_prob_line_list = \
                predict_te(model_file_dic[model_name], model_name, x_classI_ipt_test_list, y_classI_ipt_test_nm_list,input_spe_type,y_classI_ipt_test_nm_list,prop_thr)

            for i in range(len(predicted_classes_list)):
                if predicted_classes_list[i] == 0:
                    x_LTR_test_list.append(x_classI_right_list[i])
                    y_LTR_test_nm_list.append(y_classI_right_nm_list[i])

                if predicted_classes_list[i] == 1:
                    x_nLTR_test_list.append(x_classI_right_list[i])
                    y_nLTR_test_nm_list.append(y_classI_right_nm_list[i])

            ##write out wrong results
            with open (input_store_predict_dir + '/' + model_name + '_results.txt','w+') as opt:
                for eachid in store_classI_results_dic:
                    #store_all_information_dic[new_id_name] = store_classI_results_dic[eachid]
                    opt.write(store_classI_results_dic[eachid] + '\n')

            ##udpation 081520
            with open(input_store_predict_dir + '/' + model_name + '_probability_results.txt', 'w+') as opt:
                for eachline in store_prob_line_list:
                    opt.write(eachline + '\n')



    ##################
    ##detect TE in LTR
    model_name = 'LTR'

    if fam_nm == 'All' or fam_nm == 'ClassI' or fam_nm == 'LTR':

        if fam_nm == 'LTR':
            x_LTR_ipt_test_list, y_LTR_ipt_test_nm_list = DeepTE_seq_reader_kmer.load_data(input_dataset)
        else:
            x_LTR_ipt_test_list = x_LTR_test_list
            y_LTR_ipt_test_nm_list = y_LTR_test_nm_list

        ##updation
        if x_LTR_ipt_test_list != []:

            x_LTR_right_list, y_LTR_right_nm_list, store_LTR_results_dic,predicted_classes_list,store_prob_line_list = \
                predict_te(model_file_dic[model_name], model_name, x_LTR_ipt_test_list, y_LTR_ipt_test_nm_list,input_spe_type,y_LTR_ipt_test_nm_list,prop_thr)

            ##write right results
            with open(input_store_predict_dir + '/' + model_name + '_results.txt', 'w+') as opt:
                for eachid in store_LTR_results_dic:
                    #store_all_information_dic[new_id_name] = store_LTR_results_dic[eachid]
                    opt.write(store_LTR_results_dic[eachid] + '\n')

            ##udpation 081520
            with open(input_store_predict_dir + '/' + model_name + '_probability_results.txt', 'w+') as opt:
                for eachline in store_prob_line_list:
                    opt.write(eachline + '\n')


    ####################
    ##detect TE in noLTR
    model_name = 'nLTR'

    if fam_nm == 'All' or fam_nm == 'ClassI' or fam_nm == 'nLTR':

        if fam_nm == 'nLTR':
            x_nLTR_ipt_test_list,y_nLTR_ipt_test_nm_list =  DeepTE_seq_reader_kmer.load_data(input_dataset)
        else:
            x_nLTR_ipt_test_list = x_nLTR_test_list
            y_nLTR_ipt_test_nm_list = y_nLTR_test_nm_list

        ##updation
        if x_nLTR_ipt_test_list != []:

            x_nLTR_right_list, y_nLTR_right_nm_list, store_nLTR_results_dic,predicted_classes_list,store_prob_line_list = \
                predict_te(model_file_dic[model_name], model_name, x_nLTR_ipt_test_list, y_nLTR_ipt_test_nm_list,input_spe_type,y_nLTR_ipt_test_nm_list,prop_thr)

            for i in range(len(predicted_classes_list)):
                if predicted_classes_list[i] == 0:
                    x_LINE_test_list.append(x_nLTR_right_list[i])
                    y_LINE_test_nm_list.append(y_nLTR_right_nm_list[i])

                if predicted_classes_list[i] == 1:
                    x_SINE_test_list.append(x_nLTR_right_list[i])
                    y_SINE_test_nm_list.append(y_nLTR_right_nm_list[i])

                if predicted_classes_list[i] == 2:
                    store_DIRS_results_dic[str(i)] = y_nLTR_right_nm_list[i] + '\t' + 'ClassI_nLTR_DIRS'

                if predicted_classes_list[i] == 3:
                    store_PLE_results_dic[str(i)] = y_nLTR_right_nm_list[i] + '\t' + 'ClassI_nLTR_PLE'


            ##write right results
            with open(input_store_predict_dir + '/' + model_name + '_results.txt', 'w+') as opt:
                for eachid in store_nLTR_results_dic:
                    #store_all_information_dic[new_id_name] = store_noLTR_results_dic[eachid]
                    opt.write(store_nLTR_results_dic[eachid] + '\n')

            ##since O has no information for Helitron classification
            if input_spe_type != 'O':
                with open(input_store_predict_dir + '/DIRS_results.txt', 'w+') as opt:
                    for eachid in store_DIRS_results_dic:
                        #store_all_information_dic[new_id_name] = store_helitron_results_dic[eachid]
                        opt.write(store_DIRS_results_dic[eachid] + '\n')

                with open(input_store_predict_dir + '/PLE_results.txt', 'w+') as opt:
                    for eachid in store_PLE_results_dic:
                        #store_all_information_dic[new_id_name] = store_helitron_results_dic[eachid]
                        opt.write(store_PLE_results_dic[eachid] + '\n')

            ##udpation 081520
            with open(input_store_predict_dir + '/' + model_name + '_probability_results.txt', 'w+') as opt:
                for eachline in store_prob_line_list:
                    opt.write(eachline + '\n')


    ###################
    ##detect TE in LINE
    ##since O has no information for Helitron classification
    if input_spe_type != 'O':
        model_name = 'LINE'

        if fam_nm == 'All' or fam_nm == 'ClassI' or fam_nm == 'nLTR' or fam_nm == 'LINE':

            if fam_nm == 'LINE':
                x_LINE_ipt_test_list,y_LINE_ipt_test_nm_list = DeepTE_seq_reader_kmer.load_data(input_dataset)
            else:
                x_LINE_ipt_test_list = x_LINE_test_list
                y_LINE_ipt_test_nm_list = y_LINE_test_nm_list

            ##updation
            if x_LINE_ipt_test_list != []:

                x_LINE_right_list, y_LINE_right_nm_list, store_LINE_results_dic,predicted_classes_list,store_prob_line_list = \
                    predict_te(model_file_dic[model_name], model_name, x_LINE_ipt_test_list, y_LINE_ipt_test_nm_list,input_spe_type,y_LINE_ipt_test_nm_list,prop_thr)

                ##write right results
                with open(input_store_predict_dir + '/' + model_name + '_results.txt', 'w+') as opt:
                    for eachid in store_LINE_results_dic:
                        #store_all_information_dic[new_id_name] = store_LTR_results_dic[eachid]
                        opt.write(store_LINE_results_dic[eachid] + '\n')

                ##udpation 081520
                with open(input_store_predict_dir + '/' + model_name + '_probability_results.txt', 'w+') as opt:
                    for eachline in store_prob_line_list:
                        opt.write(eachline + '\n')

        ###################
        ##detect TE in LINE
        model_name = 'SINE'

        if fam_nm == 'All' or fam_nm == 'ClassI' or fam_nm == 'nLTR' or fam_nm == 'SINE':

            if fam_nm == 'SINE':
                x_SINE_ipt_test_list,y_SINE_ipt_test_nm_list = DeepTE_seq_reader_kmer.load_data(input_dataset)
            else:
                x_SINE_ipt_test_list = x_SINE_test_list
                y_SINE_ipt_test_nm_list = y_SINE_test_nm_list

            ##updation
            if x_SINE_ipt_test_list != []:

                x_SINE_right_list, y_SINE_right_nm_list, store_SINE_results_dic,predicted_classes_list,store_prob_line_list = \
                    predict_te(model_file_dic[model_name], model_name, x_SINE_ipt_test_list, y_SINE_ipt_test_nm_list,input_spe_type,y_SINE_ipt_test_nm_list,prop_thr)

                ##write right results
                with open(input_store_predict_dir + '/' + model_name + '_results.txt', 'w+') as opt:
                    for eachid in store_SINE_results_dic:
                        #store_all_information_dic[new_id_name] = store_LTR_results_dic[eachid]
                        opt.write(store_SINE_results_dic[eachid] + '\n')

                ##udpation 081520
                with open(input_store_predict_dir + '/' + model_name + '_probability_results.txt', 'w+') as opt:
                    for eachline in store_prob_line_list:
                        opt.write(eachline + '\n')


    ###########################
    ##detect TE in the class II
    model_name = 'ClassII'

    if fam_nm == 'All' or fam_nm == 'ClassII':

        if fam_nm == 'ClassII':
            x_classII_ipt_test_list, y_classII_ipt_test_nm_list = DeepTE_seq_reader_kmer.load_data(input_dataset)
        else:
            x_classII_ipt_test_list = x_classII_test_list
            y_classII_ipt_test_nm_list = y_classII_test_nm_list

        ##updation
        if x_classII_ipt_test_list != []:

            x_classII_right_list, y_classII_right_nm_list, store_classII_results_dic,predicted_classes_list,store_prob_line_list = \
                predict_te(model_file_dic[model_name], model_name, x_classII_ipt_test_list, y_classII_ipt_test_nm_list,input_spe_type,y_classII_ipt_test_nm_list,prop_thr)


            ##write out wrong results
            with open(input_store_predict_dir + '/' + model_name + '_results.txt', 'w+') as opt:
                for eachid in store_classII_results_dic:
                    #store_all_information_dic[new_id_name] = store_DNA_results_dic[eachid]
                    opt.write(store_classII_results_dic[eachid] + '\n')

            ##udpation 081520
            with open(input_store_predict_dir + '/' + model_name + '_probability_results.txt', 'w+') as opt:
                for eachline in store_prob_line_list:
                    opt.write(eachline + '\n')



    #########################
    ##classify MITE and nMITE
    model_name = 'Domain'

    if fam_nm == 'All' or fam_nm == 'Domain' or fam_nm == 'ClassII':

        x_classII_ipt_test_list = ''
        y_classII_ipt_test_nm_list = ''
        if fam_nm == 'All':
            x_classII_ipt_test_list = x_classII_test_list
            y_classII_ipt_test_nm_list = y_classII_test_nm_list
        if fam_nm == 'Domain' or fam_nm == 'ClassII':
            x_classII_ipt_test_list, y_classII_ipt_test_nm_list = DeepTE_seq_reader_kmer.load_data(input_dataset)

        ##updation
        if x_classII_ipt_test_list != []:

            x_domain_right_list, y_domain_right_nm_list, store_all_results_dic, predicted_classes_list,store_prob_line_list = \
                predict_te(model_file_dic[model_name], model_name, x_classII_ipt_test_list, y_classII_ipt_test_nm_list, input_spe_type,y_classII_ipt_test_nm_list,prop_thr)

            with open (input_store_predict_dir + '/' + model_name + '_results.txt','w+') as opt:
                for eachid in store_all_results_dic:
                    opt.write(store_all_results_dic[eachid] + '\n')

            ##udpation 081520
            with open(input_store_predict_dir + '/' + model_name + '_probability_results.txt', 'w+') as opt:
                for eachline in store_prob_line_list:
                    opt.write(eachline + '\n')







