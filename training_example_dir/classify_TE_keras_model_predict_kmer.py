#!/usr/bin/env python


##change the parameters

################################################
##this script is to train the data and test data

from seq_reader_kmer import load_data        # parsing data file
from one_hot_rep_kmer import generate_mats, conv_labels   # converting to correct format


from sklearn.model_selection import StratifiedKFold     # cross validation

import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras.utils import np_utils
from keras.datasets import mnist
from keras.models import load_model

import re

##import classification_report helps us to detect the accuracy for each specific class
from sklearn.metrics import classification_report



def train_model (input_dataset,input_store_class_report_dir,store_all_score_dic):


    ##get the name of input data set
    print('the input dataset is ' + input_dataset)
    mt = re.match('.+/ipt_shuffle_(.+)_CNN_data.txt',input_dataset)
    input_data_nm = mt.group(1)

    ##change the class_num according to the input_data_nm

    class_num = 0 ##input the class number which is a binary classification
    #['All', 'ClassI', 'ClassII', 'LTR', 'noLTR', 'DNA', 'MITE', 'noMITE']
    if input_data_nm == 'All':
        class_num = 18
    if input_data_nm == 'ClassI':
        class_num = 2
    if input_data_nm == 'ClassII':
        class_num = 5
    if input_data_nm == 'LTR':
        class_num = 4
    if input_data_nm == 'noLTR':
        class_num = 3
    if input_data_nm == 'DNA':
        class_num = 5
    if input_data_nm == 'MITE':
        class_num = 5
    if input_data_nm == 'noMITE':
        class_num = 5


    #######################################
    # 1. Load data into train and test sets
    X, y = load_data(input_dataset) # sequences, labels
    #print(X)            ##Note: X is the ['AA','BB','CC]
    #print(y)            ##Note: y is the ['C1','C2','C3']


    X = generate_mats(X)     # convert to array of representation matrices


    # convert to integer labels
    y = conv_labels(y,input_data_nm)
    # work with np arrays
    X = np.asarray(X)
    #print(X)
    #print(np.shape(X))
    Y = np.asarray(y)
    #print(Y)

    ##calculate the dateset number and decide 10% as the test data
    seq_count = 0
    with open (input_dataset,'r') as ipt:
        for eachline in ipt:
            seq_count += 1

    divide_data_most_part = int(seq_count * 0.9)


    ##generate 13000 as train a total of 17036
    X_train = X[0:divide_data_most_part]  ##change the sample
    X_test = X[divide_data_most_part:]
    Y_train = y[0:divide_data_most_part]
    Y_test = y[divide_data_most_part:]

    #print('not change shape X_train is ')
    #print(X_train)
    #print(np.shape(X_train))
    #print('\n')

    ##########################
    # 2. Preprocess input data
    X_train = X_train.reshape(X_train.shape[0], 1, 16384, 1)  ##shape[0] indicates sample number
    #print('X_train is ')
    #print(X_train)

    X_test = X_test.reshape(X_test.shape[0], 1, 16384, 1)     ##kmer == 3 so it would be 64
    X_train = X_train.astype('float64')
    X_test = X_test.astype('float64')


    ###########################
    # 3. Preprocess class labels; i.e. convert 1-dimensional class arrays to 3-dimensional class matrices
    Y_train_one_hot = np_utils.to_categorical(Y_train, int(class_num))  # four labels
    Y_test_one_hot = np_utils.to_categorical(Y_test, int(class_num)) #four labels

    ###########################
    # 4. Define model architecture
    model = Sequential()

    model.add(Conv2D(100, (1, 3), activation='relu', input_shape=(1, 16384, 1)))
    model.add(MaxPooling2D(pool_size=(1, 2)))
    model.add(Conv2D(150, (1, 3), activation='relu'))
    model.add(MaxPooling2D(pool_size=(1, 2)))
    model.add(Conv2D(225, (1, 3), activation='relu'))
    model.add(MaxPooling2D(pool_size=(1, 2)))
    model.add(Dropout(0.5))

    model.add(Flatten())
    model.add(Dense(128, activation='relu'))
    model.add(Dropout(0.5))
    ##You can add a dropout layer to overcome the problem of overfitting to some extent. Dropout randomly turns off
    # a fraction of neurons during the training process, reducing the dependency on the training set by some amount.
    # How many fractions of neurons you want to turn off is decided by a hyperparameter, which can be tuned accordingly.
    # This way, turning off some neurons will not allow the network to memorize the training data since not all the neurons
    # will be active at the same time and the inactive neurons will not be able to learn anything.
    #This way, turning off some neurons will not allow the network to memorize the training data
    # since not all the neurons will be active at the same time and the inactive neurons will not be able to learn anything.


    model.add(Dense(int(class_num), activation='softmax'))
    # since 4 classes ##the output have four unit
    ##Your output's are integers for class labels. Sigmoid logistic function outputs values in range (0,1).
    # The output of the softmax is also in range (0,1), but the softmax function adds another constraint on outputs:-
    # the sum of outputs must be 1. Therefore the output of softmax can be interpreted as probability of the input
    # for each class.

    model.compile(loss='categorical_crossentropy',optimizer='adam',metrics=['accuracy'])

    ###########################
    # 6. Fit model on training data
    model.fit(X_train, Y_train_one_hot, validation_data=(X_test, Y_test_one_hot),
              batch_size=32, epochs=10, verbose=1)   ##epoch An epoch is an iteration over the entire x and y data provided
                                                    ##batch size if we have 1000 samples, we set batch size to 100, so we will
                                                    ##run 100 first and then the second 100, so this will help us to reduce the
                                                    ##the memory we use
    ###########################
    # 7. Evaluate model on test data
    score = model.evaluate(X_test, Y_test_one_hot, verbose=1)
    print ("\nscore = " + str(score))

    store_all_score_dic[input_data_nm] = str(score)

    #with open (input_store_class_report_dir + '/opt_score_results.txt','w+') as opt:
    #    opt.write("\nscore = " + str(score) + '\n')


    ###########################
    # 7.5.  save the model
    model.save(input_store_class_report_dir + '/' + input_data_nm + '_model.h5')

    ###########################
    # 8. Classification report

    ##predict label comparing with true data
    predicted_classes = model.predict(X_test)
    predicted_classes = np.argmax(np.round(predicted_classes),axis=1)

    #print('predicted_classes is ')
    #print(predicted_classes)
    #print('test_Y is ')
    #print(Y_test)

    ##change the Y_test to array
    ##Y_test is a list
    Y_test = np.asarray(Y_test)


    correct = np.where(predicted_classes==Y_test)[0]

    print ("Found %d correct labels" % len(correct))
    for i, correct in enumerate(correct[:int(class_num)]):
        #plt.subplot(3,3,i+1)
        #plt.imshow(test_X[correct].reshape(28,28), cmap='gray', interpolation='none')
        print("Predicted {}, Class {}".format(predicted_classes[correct], Y_test[correct]))


    target_names = ["Class {}".format(i) for i in range(int(class_num))] ##there are four classes

    with open (input_store_class_report_dir + '/' + input_data_nm + '_class_report.txt','w+') as opt:
        opt.write(classification_report(Y_test, predicted_classes, target_names=target_names) + '\n')