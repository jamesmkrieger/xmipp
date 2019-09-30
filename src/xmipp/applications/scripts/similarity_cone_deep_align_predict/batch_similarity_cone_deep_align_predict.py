#!/usr/bin/env python2

from keras.callbacks import TensorBoard, ModelCheckpoint
from keras.models import Model
from keras.layers import Input, Conv2D, MaxPooling2D, BatchNormalization, Dropout, Flatten, Dense
from keras.optimizers import Adam
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
#import pyworkflow.em.metadata as md
from shutil import copy

def loadData(mdIn, mdExp):
    XProj=None
    XExp=None
    Nproj = mdIn.size()
    Nexp = mdExp.size()
    idx = 0
    for objId in mdIn:
        fnImg = mdIn.getValue(xmippLib.MDL_IMAGE,objId)
        I = xmippLib.Image(fnImg)
        if XProj is None:
            Xdim, Ydim, _, _ = I.getDimensions()
            XProj = np.zeros((Nproj,Xdim,Ydim,1),dtype=np.float64)
	XProj[idx,:,:,0] = I.getData()
        idx+=1

    idx = 0
    for objId in mdExp:
        fnImg = mdExp.getValue(xmippLib.MDL_IMAGE,objId)
        I = xmippLib.Image(fnImg)
        if XExp is None:
            Xdim, Ydim, _, _ = I.getDimensions()
            XExp = np.zeros((Nexp,Xdim,Ydim,1),dtype=np.float64)
	XExp[idx,:,:,0] = I.getData()
        idx+=1

    return XProj, XExp, Xdim, Ydim, Nproj, Nexp


if __name__=="__main__":
    fnExp = sys.argv[1]
    fnODir = sys.argv[2]
    Xdim = int(sys.argv[3])
    numClassif = int(sys.argv[4])
    numMax=int(sys.argv[5])
    Nexp = int(sys.argv[6])

    print('Predict mode')
    newImage = xmippLib.Image()
    testSet = open(fnExp, "r")

    sizeBatch=1000
    maxBatchs=np.ceil(float(Nexp)/float(sizeBatch))
    Ypred = np.zeros((Nexp),dtype=np.float64)   
    refPred = np.zeros((Nexp,(numMax*2)+1),dtype=np.float64)    
    models=[]
    for i in range(numClassif):
	if os.path.exists(os.path.join(fnODir,'modelCone%d.h5'%(i+1))):
	    models.append(load_model(os.path.join(fnODir,'modelCone%d.h5'%(i+1))))
    if Nexp>sizeBatch:
        oneXExp=[np.zeros((sizeBatch, Xdim, Xdim, 1),dtype=np.float64) for i in range(2)]
        YpredAux = np.zeros((sizeBatch,numClassif),dtype=np.float64)

    idxExp = 0
    countBatch=0
    numBatch = 0
    done = 0
    count=0
    lines = trainSet.readlines()
    for line in lines:
	if numBatch==(maxBatchs-1) and done==0:
	    oneXExp=[np.zeros((Nexp-idxExp, Xdim, Xdim, 1),dtype=np.float64) for i in range(2)]
            YpredAux = np.zeros((Nexp-idxExp,numClassif),dtype=np.float64)
	    done=1
	fnExp = line.split()[0]
	Iexp = xmippLib.Image(fnExp)
        oneXExp[0][countBatch,:,:,0] = Iexp.getData()
	oneXExp[0][countBatch,:,:,0] = (oneXExp[0][countBatch,:,:,0]-np.mean(oneXExp[0][countBatch,:,:,0]))/np.std(oneXExp[0][countBatch,:,:,0])
	fnProj = line.split()[1]
	Iproj = xmippLib.Image(fnProj)
        oneXExp[1][countBatch,:,:,0] = Iproj.getData()
	oneXExp[1][countBatch,:,:,0] = (oneXExp[1][countBatch,:,:,0]-np.mean(oneXExp[1][countBatch,:,:,0]))/np.std(oneXExp[1][countBatch,:,:,0])
	countBatch+=1
	idxExp+=1
	refPred[idxExp-1,0] = count
        count=count+1
	if ((idxExp%sizeBatch)==0 or idxExp==Nexp):
	    countBatch = 0
            for i in range(numClassif):
	        model = models[i]
                out = model.predict([oneXExp])
		YpredAux[:,i] = out[:,0]
	    if numBatch==(maxBatchs-1):
		for n in range(numMax):
                    Ypred[numBatch*sizeBatch:Nexp] = np.max(YpredAux, axis=1)
	            auxPos = np.argmax(YpredAux, axis=1)
	            refPred[numBatch*sizeBatch:Nexp, (n*2)+1] = np.argmax(YpredAux, axis=1)+1
	            refPred[numBatch*sizeBatch:Nexp, (n*2)+2] = Ypred[numBatch*sizeBatch:Nexp]
		    for i,pos in enumerate(auxPos):
		        YpredAux[i,pos]=0.0
	    else:
		for n in range(numMax):
                    Ypred[idxExp-sizeBatch:idxExp] = np.max(YpredAux, axis=1)
	            auxPos = np.argmax(YpredAux, axis=1)
   		    #print(n, YpredAux, len(auxPos))
	            refPred[idxExp-sizeBatch:idxExp, (n*2)+1] = np.argmax(YpredAux, axis=1)+1
	            refPred[idxExp-sizeBatch:idxExp, (n*2)+2] = Ypred[idxExp-sizeBatch:idxExp]
		    for i,pos in enumerate(auxPos):
		        YpredAux[i,pos]=0.0
	    numBatch+=1

    np.savetxt(os.path.join(fnODir,'conePrediction.txt'), refPred)








