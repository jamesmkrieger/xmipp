#!/usr/bin/env python2

from keras.callbacks import TensorBoard, ModelCheckpoint
from keras.models import Model, Sequential
from keras.layers import Input, Conv2D, MaxPooling2D, BatchNormalization, Dropout, Flatten, Dense, concatenate, Subtract, SeparableConv2D, GlobalAveragePooling2D, Lambda, merge
from keras.optimizers import *
import tensorflow as tf
import cv2
import math
import numpy as np
import os
import random
import string
import sys
import xmippLib
import time
from keras.models import load_model
from time import time
import keras
from keras import callbacks
from keras.callbacks import Callback
from keras import regularizers

#batch_size = 516 # Number of boxes per batch


class EarlyStoppingByLossVal(Callback):
    def __init__(self, monitor='val_loss', value=0.30, verbose=0):
        super(Callback, self).__init__()
        self.monitor = monitor
        self.value = value
        self.verbose = verbose

    def on_epoch_end(self, epoch, logs={}):
        current = logs.get(self.monitor)
        if current is None:
            warnings.warn("Early stopping requires %s available!" % self.monitor, RuntimeWarning)

        if current < self.value:
            if self.verbose > 0:
                print("Epoch %05d: early stopping THR" % epoch)
            self.model.stop_training = True


class DataGenerator(keras.utils.Sequence):
    'Generates data for Keras'
    def __init__(self, list_IDs_zeros, list_IDs_ones, labels, batch_size, dim, 
                 shuffle, pathsTrain1, pathsTrain2):
        'Initialization'
        self.dim = dim
        self.batch_size = batch_size
        self.labels = labels
        self.list_IDs_zeros = list_IDs_zeros
        self.list_IDs_ones = list_IDs_ones
        self.shuffle = shuffle
	self.pathsTrain1 = pathsTrain1
	self.pathsTrain2 = pathsTrain2
        self.on_epoch_end()

    def __len__(self):
        'Denotes the number of batches per epoch'
        return int(np.floor((len(self.labels)) / self.batch_size))

    def __getitem__(self, index):
        'Generate one batch of data'
        # Generate indexes of the batch
        indexes_zeros = self.indexes_zeros[index*(self.batch_size//2):(index+1)*(self.batch_size//2)]
	indexes_ones = self.indexes_ones[index*(self.batch_size//2):(index+1)*(self.batch_size//2)]

        # Find list of IDs
	#list_IDs_temp = [self.list_IDs[k] for k in indexes]
	list_IDs_temp = []
	for i in range(int(self.batch_size//2)):
            list_IDs_temp.append(indexes_zeros[i])
	for i in range(int(self.batch_size//2)):
            list_IDs_temp.append(indexes_ones[i])

        # Generate data
        Xexp, y = self.__data_generation(list_IDs_temp)

        return Xexp, y

    def on_epoch_end(self):
        'Updates indexes after each epoch'
        self.indexes_zeros = self.list_IDs_zeros
	self.indexes_ones = self.list_IDs_ones
        if self.shuffle == True:
            np.random.shuffle(self.indexes_zeros)
            np.random.shuffle(self.indexes_ones)

    def __data_generation(self, list_IDs_temp):
        'Generates data containing batch_size samples' # X : (n_samples, *dim, n_channels)
        # Initialization
	Xexp = [np.zeros((self.batch_size,self.dim,self.dim,1),dtype=np.float64) for i in range(2)]
        y = np.empty((self.batch_size), dtype=np.int64)

        # Generate data
        for i, ID in enumerate(list_IDs_temp):
            # Store sample
            Iexp = np.reshape(xmippLib.Image(self.pathsTrain1[ID]).getData(),(self.dim,self.dim,1))
	    Iexp = (Iexp-np.mean(Iexp))/np.std(Iexp)
	    Xexp[0][i,:,:,:] = Iexp
            Iexp2 = np.reshape(xmippLib.Image(self.pathsTrain2[ID]).getData(),(self.dim,self.dim,1))
	    Iexp2 = (Iexp2-np.mean(Iexp2))/np.std(Iexp2)
	    Xexp[1][i,:,:,:] = Iexp2

            # Store class
            y[i] = self.labels[ID]

        return Xexp, y


def createValidationData(pathsTrain1, pathsTrain2, labels_vector, Xdim, percent=0.1):
    sizeValData = int(round(len(pathsTrain1)*percent))
    print("sizeValData",sizeValData)
    sys.stdout.flush()
    val_img_exp = []
    val_labels = []
    labels = np.array(labels_vector)	
    vectorOnes = np.where(labels==1)
    numberOnes = int(len(vectorOnes[0])*percent)
    vectorZeros = np.where(labels==0)
    numberZeros = int(len(vectorZeros[0])*percent)
    if numberZeros>numberOnes:
        numberZeros= numberOnes
    elif numberOnes>numberZeros:
        numberOnes=numberZeros
    numberZeros=100 #AJ just to test
    numberOnes=100 #AJ just to test
    pairs=[np.zeros((numberOnes+numberZeros, Xdim, Xdim, 1)) for i in range(2)]
    targets=np.zeros((numberOnes+numberZeros,))
    for i in range(numberOnes):
        labels = np.array(labels_vector)
	vectorOnes = np.where(labels==1)
	vectorOnes = vectorOnes[0]
	k = np.random.randint(0,len(vectorOnes))
	k = vectorOnes[k]
	Iexp = xmippLib.Image(pathsTrain1[k])
	Iexp = np.reshape(Iexp.getData(),(Xdim,Xdim,1))
	Iexp = (Iexp-np.mean(Iexp))/np.std(Iexp)
	pairs[0][i,:,:,:] = Iexp
        targets[i]=labels_vector[k]
	Iexp2 = xmippLib.Image(pathsTrain2[k])
	Iexp2 = np.reshape(Iexp2.getData(),(Xdim,Xdim,1))
	Iexp2 = (Iexp2-np.mean(Iexp2))/np.std(Iexp2)
	pairs[1][i,:,:,:] = Iexp2
	del pathsTrain1[k]
	del labels_vector[k]
	del pathsTrain2[k]

    for i in range(numberZeros):
        labels = np.array(labels_vector)
	vectorZeros = np.where(labels==0)
	vectorZeros = vectorZeros[0]
	k = np.random.randint(0,len(vectorZeros))
	k = vectorZeros[k]
	Iexp = xmippLib.Image(pathsTrain1[k])
	Iexp = np.reshape(Iexp.getData(),(Xdim,Xdim,1))
	Iexp = (Iexp-np.mean(Iexp))/np.std(Iexp)
	pairs[0][i,:,:,:] = Iexp
        targets[i]=labels_vector[k]
	Iexp2 = xmippLib.Image(pathsTrain2[k])
	Iexp2 = np.reshape(Iexp2.getData(),(Xdim,Xdim,1))
	Iexp2 = (Iexp2-np.mean(Iexp2))/np.std(Iexp2)
	pairs[1][i,:,:,:] = Iexp2
	del pathsTrain1[k]
	del labels_vector[k]
	del pathsTrain2[k]

    return np.asarray(pairs).astype('float64'), np.asarray(targets).astype('int64')



def constructModel(Xdim):

    input_shape = (Xdim,Xdim, 1)
    left_input = Input(input_shape)
    right_input = Input(input_shape)
    #build convnet to use in each siamese 'leg'
    convnet = Sequential()
    convnet.add(Conv2D(64, (10,10), activation='relu', input_shape=input_shape))
    convnet.add(MaxPooling2D())
    convnet.add(Conv2D(128, (7,7), activation='relu'))
    convnet.add(MaxPooling2D())
    convnet.add(Conv2D(128, (4,4), activation='relu'))
    convnet.add(MaxPooling2D())
    convnet.add(Conv2D(256, (4,4), activation='relu'))
    convnet.add(Flatten())
    convnet.add(Dense(4096, activation="sigmoid"))
    #encode each of the two inputs into a vector with the convnet
    encoded_l = convnet(left_input)
    encoded_r = convnet(right_input)
    #merge two encoded inputs with the l1 distance between them
    L1_distance = lambda x: K.abs(x[0]-x[1])
    both = merge([encoded_l, encoded_r], mode = L1_distance, output_shape=lambda x: x[0])
    prediction = Dense(1, activation='sigmoid')(both)
    siamese_net = Model(input=[left_input, right_input], output=prediction)
    #optimizer = SGD(0.0004,momentum=0.6,nesterov=True,decay=0.0003)

    return siamese_net

if __name__=="__main__":
    fnTrain = sys.argv[1]
    fnLabels = sys.argv[2]
    fnODir = sys.argv[3]
    modelFn = sys.argv[4]
    numEpochs = int(sys.argv[5])
    Xdim = int(sys.argv[6])
    batch_size = int(sys.argv[7])

    trainSet = open(fnTrain, "r")
    labels = open(fnLabels, "r")

    #AJ new code to generate data for validation set
    pathsTrain1 = []
    pathsTrain2 = []
    labels_vector = []
    cont=0
    linesTrain = trainSet.readlines()
    linesLabels = labels.readlines()
    for line, lineLabel in zip(linesTrain, linesLabels):
        pathsTrain1.append(line.split()[0])
        pathsTrain2.append(line.split()[1])
        labels_vector.append(int(lineLabel))
        cont+=1
    print(len(labels_vector))
    sys.stdout.flush()
    start_time = time()
    x_val, y_val = createValidationData(pathsTrain1, pathsTrain2, labels_vector, Xdim, 0.2)
    elapsed_time = time() - start_time
    print("Time in createValidationData: %0.10f seconds." % elapsed_time)
    np.savetxt(os.path.join(fnODir,'pruebaYval.txt'), y_val)
    print(len(labels_vector))
    sys.stdout.flush()
    #END AJ

    # Parameters
    params = {'dim': Xdim,
          'batch_size': batch_size,
          'shuffle': True,
	  'pathsTrain1': pathsTrain1,
	  'pathsTrain2': pathsTrain2}
    # Datasets
    list_IDs_zeros = np.where(np.array(labels_vector)==0)
    list_IDs_ones = np.where(np.array(labels_vector)==1)
    list_IDs_zeros = list_IDs_zeros[0]
    list_IDs_ones = list_IDs_ones[0]
    list_IDs_zeros_orig = list_IDs_zeros
    list_IDs_ones_orig = list_IDs_ones
    lenTotal = len(list_IDs_zeros)+len(list_IDs_ones)
    if len(list_IDs_zeros)<lenTotal:
	for i in range((lenTotal//len(list_IDs_zeros))-1):
	    list_IDs_zeros = np.append(list_IDs_zeros, list_IDs_zeros_orig)
	list_IDs_zeros = np.append(list_IDs_zeros, list_IDs_zeros[0:(lenTotal%len(list_IDs_zeros))])
    if len(list_IDs_ones)<lenTotal:
	for i in range((lenTotal//len(list_IDs_ones))-1):
	    list_IDs_ones = np.append(list_IDs_ones, list_IDs_ones_orig)
	list_IDs_ones = np.append(list_IDs_ones, list_IDs_ones[0:(lenTotal%len(list_IDs_ones))])
    print(len(list_IDs_zeros), len(list_IDs_ones))
    labels = labels_vector
    # Generator
    training_generator = DataGenerator(list_IDs_zeros, list_IDs_ones, labels, **params)

    print('Train mode')
    start_time = time()
    model = constructModel(Xdim)

    name_model = os.path.join(fnODir, modelFn+'.h5')
    
    #callbacks_list = [callbacks.ModelCheckpoint(filepath=name_model, monitor='val_loss', save_best_only=True),
    #		      EarlyStoppingByLossVal(monitor='val_loss', value=0.30)]

    callbacks_list = [callbacks.ModelCheckpoint(filepath=name_model, monitor='val_loss', save_best_only=True)]


    model.summary()
    adam_opt = Adam(lr=0.001)
    #model.compile(loss='binary_crossentropy', optimizer='Adam')
    model.compile(loss='mean_absolute_error', optimizer=adam_opt)

    steps = round(len(pathsTrain1)/batch_size)
    history = model.fit_generator(generator = training_generator, steps_per_epoch = steps, epochs=numEpochs, verbose=1, validation_data = ([x_val[0], x_val[1]], y_val), callbacks=callbacks_list, workers=4, use_multiprocessing=True)
    elapsed_time = time() - start_time
    print("Time in training model: %0.10f seconds." % elapsed_time)

    model = load_model(name_model)
    Ypred = model.predict(x_val)
    np.savetxt(os.path.join(fnODir,'pruebaYpred.txt'), Ypred)
    #mae= np.mean(np.absolute(Ypred-y_val))
    from sklearn.metrics import mean_absolute_error
    mae = mean_absolute_error(y_val, Ypred)
    print("Final model mean absolute error val_loss", mae)



