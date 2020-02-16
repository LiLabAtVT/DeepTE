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
import matplotlib.pyplot as plt

import re

##import classification_report helps us to detect the accuracy for each specific class
from sklearn.metrics import classification_report


def write_acc_epoch (input_eval_list,input_store_class_report_dir,input_data_nm,eval_name):

    ##store the accuracy, and loss to the opt
    epoch_dic = {}
    epoch_acc_count = 0
    for eacheval in input_eval_list:
        epoch_acc_count += 1
        epoch_dic[str(epoch_acc_count)] = str(eacheval)
    ##write out results
    with open(input_store_class_report_dir + '/' + input_data_nm + '_epoch_' + eval_name + '.txt', 'w+') as opt:
        for eachepoch in epoch_dic:
            opt.write(eachepoch + '\t' + epoch_dic[eachepoch] + '\n')



def generate_input_data (input_dataset,input_data_nm):

    ##for the training dataset
    X, y = load_data(input_dataset)  # sequences, labels
    # print(X)            ##Note: X is the ['AA','BB','CC]
    # print(y)            ##Note: y is the ['C1','C2','C3']
    X = generate_mats(X)  # convert to array of representation matrices
    # convert to integer labels
    y = conv_labels(y, input_data_nm)
    # work with np arrays
    X = np.asarray(X)
    # print(X)
    # print(np.shape(X))
    Y = np.asarray(y)
    # print(Y)

    return (X,y)


def train_model (input_train_dataset,input_test_dataset,input_store_class_report_dir,store_all_score_dic):


    ##get the name of input train data set
    mt = re.match('.+\/shuffle_(.+)_CNN_data.txt',input_train_dataset)
    input_data_nm = mt.group(1)

    ##change the class_num according to the input_data_nm

    class_num = 0 ##input the class number which is a binary classification
    #['All', 'ClassI', 'ClassII', 'LTR', 'noLTR', 'DNA', 'MITE', 'noMITE']
    if input_data_nm == 'All':
        class_num = 3
    if input_data_nm == 'ClassI':
        class_num = 2
    if input_data_nm == 'LTR':
        class_num = 2
    if input_data_nm == 'nLTR':
        class_num = 4
    if input_data_nm == 'LINE':
        class_num = 2
    if input_data_nm == 'SINE':
        class_num = 2
    if input_data_nm == 'ClassII':
        class_num = 6
    if input_data_nm == 'Domain':
        class_num = 2


    X_train, Y_train = generate_input_data(input_train_dataset, input_data_nm)
    X_test, Y_test = generate_input_data(input_test_dataset, input_data_nm)


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
    te_train = model.fit(X_train, Y_train_one_hot, validation_data=(X_test, Y_test_one_hot),
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

    print('the corret lable is')
    print(correct)

    wrong = np.where(predicted_classes != Y_test)[0]
    print('the wrong label is')
    print(wrong)


    ##generate a table showing with two columns showing the predicted class and the original class
    ##initiate a dic to store the correct and wrong lists
    #predicted_classes_list = predicted_classes.tolist()  # list
    #Y_test_list = Y_test.tolist()
    #len(Y_test_list)

    #for i in range(0,len(Y_test_list)):
    #    print(i)


    print('the predicted_classes is ')
    print(predicted_classes)

    print('the test_classes is ')
    print(Y_test)





    print ("Found %d correct labels" % len(correct))
    for i, correct in enumerate(correct[:int(class_num)]):
        #plt.subplot(3,3,i+1)
        #plt.imshow(test_X[correct].reshape(28,28), cmap='gray', interpolation='none')
        print("Predicted {}, Class {}".format(predicted_classes[correct], Y_test[correct]))


    target_names = ["Class {}".format(i) for i in range(int(class_num))] ##there are four classes

    with open (input_store_class_report_dir + '/' + input_data_nm + '_class_report.txt','w+') as opt:
        opt.write(classification_report(Y_test, predicted_classes, target_names=target_names) + '\n')


    ##generate train loss figure for each fam order
    accuracy = te_train.history['acc']
    val_accuracy = te_train.history['val_acc']
    loss = te_train.history['loss']
    val_loss = te_train.history['val_loss']
    epochs = range(len(accuracy))



    ##############
    ##generate the
    write_acc_epoch(accuracy, input_store_class_report_dir, input_data_nm, 'accuracy')
    write_acc_epoch(val_accuracy, input_store_class_report_dir, input_data_nm, 'val_accuracy')
    write_acc_epoch(loss, input_store_class_report_dir, input_data_nm, 'loss')
    write_acc_epoch(val_loss, input_store_class_report_dir, input_data_nm, 'val_loss')


    ##create one figure
    fig = plt.figure()
    plt.plot(epochs, accuracy, 'bo', label='Training accuracy')
    plt.plot(epochs, val_accuracy, 'b', label='Validation accuracy')
    plt.title('Training and validation accuracy')
    plt.legend()
    fig.savefig('plot_acc_' + input_data_nm + '.png')

    fig = plt.figure()
    plt.plot(epochs, loss, 'bo', label='Training loss')
    plt.plot(epochs, val_loss, 'b', label='Validation loss')
    plt.title('Training and validation loss')
    plt.legend()
    fig.savefig('plot_loss_' + input_data_nm + '.png')
