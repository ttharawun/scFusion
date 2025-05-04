# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 21:38:29 2019

@author: BioMed-X
"""


from keras.models import Sequential,Model
from keras.layers import Embedding,Dropout,Bidirectional,Flatten,Dense,LSTM,TimeDistributed, Activation,Input,merge,concatenate
from keras.callbacks import ModelCheckpoint,CSVLogger
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D
from tensorflow.keras.optimizers import Adam
import numpy as np
from tensorflow.keras.utils import to_categorical
import os
import sys
import random
import tensorflow as tf

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

def Cla_LSTM():
    input_layer = Input(shape=(61,))
    embedded = Embedding(input_dim=5, output_dim=5)(input_layer)
    lstm = Bidirectional(LSTM(32))(embedded)
    dense = Dense(32, activation='relu')(lstm)
    output = Dense(2, activation='softmax')(dense)
    model = Model(inputs=input_layer, outputs=output)
    return model

# In your __main__ block (replace fit block):
model = Cla_LSTM()
model.compile(loss='binary_crossentropy', optimizer=Adam(learning_rate=0.001), metrics=['accuracy'])

print("Training on 100 samples only (debug mode)")
model.fit(
    x=Tra_x[:100], y=Tra_y[:100],
    batch_size=16,
    epochs=2,
    validation_data=(Tst_x[:50], Tst_y[:50]),
    verbose=1,
    shuffle=True
)

if __name__ == '__main__':
    np.random.seed(1122)
    random.seed(1122)
    tf.random.set_seed(1122)
    npydir = sys.argv[1]
    weightfile = sys.argv[2]
    epochoutdir = sys.argv[3]
    itere = int(sys.argv[4])
    Good_for_Tra = np.load(npydir + '/Good_for_Tra.npy')
    Simu_for_Tra = np.load(npydir + '/Simu_for_Tra.npy')
    Good_for_Tst = np.load(npydir + '/Good_for_Tst.npy')
    Simu_for_Tst = np.load(npydir + '/Simu_for_Tst.npy')
    Tra_x = np.squeeze(np.concatenate((Good_for_Tra,Simu_for_Tra),axis=0))
    Tra_y = np.concatenate( (np.zeros((Good_for_Tra.shape[0],1)),np.ones((Simu_for_Tra.shape[0],1))), axis=0)
    Tst_x = np.squeeze(np.concatenate((Good_for_Tst,Simu_for_Tst),axis=0))
    Tst_y = np.concatenate((np.zeros((Good_for_Tst.shape[0], 1)), np.ones((Simu_for_Tst.shape[0], 1))), axis=0)
    Tra_y = to_categorical(Tra_y)
    Tst_y = to_categorical(Tst_y)

    LIST = list(range(Tra_x.shape[0]))
    np.random.shuffle(LIST)

    Tra_x = Tra_x[LIST,:]
    Tra_y = Tra_y[LIST,:]
    # Tra_x_input1 = Tra_x[..., 0]
    # Tra_x_input2 = Tra_x[..., 1][...,np.newaxis]

    LIST = list(range(Tst_x.shape[0]))
    np.random.shuffle(LIST)
    Tst_x = Tst_x[LIST,:]
    Tst_y = Tst_y[LIST,:]

    model = Cla_LSTM()
    model.load_weights(weightfile)

    ADAM = Adam(learning_rate=0.0001)
    model_checkpoint = ModelCheckpoint(filepath=epochoutdir + '/RetrainWeight-{epoch:03d}.hdf5', verbose=1, monitor='val_loss', save_best_only=True)
    model_checkpoint2 = ModelCheckpoint(filepath=epochoutdir + '/RetrainWeight.hdf5', verbose=1, monitor='val_loss', save_best_only=True)
    model.compile(loss='binary_crossentropy', optimizer=ADAM, metrics=['accuracy'])
    csv_loger=CSVLogger('log.csv',append=True,separator=';')

    # 训练模型
    batch_size = 500
    epochs = itere
    np.random.seed(1122)
    # model.fit(x=[Tra_x_input1,Tra_x_input2], y=Tra_y,batch_size=batch_size,epochs=epochs, verbose=1 ,callbacks=[model_checkpoint,csv_loger], validation_split=0.25, shuffle=True)
    model.fit(x=Tra_x, y=Tra_y,batch_size=batch_size,epochs=epochs,validation_data=(Tst_x, Tst_y),verbose=1 ,callbacks=[model_checkpoint,model_checkpoint2, csv_loger], shuffle=False)


# 

'''
model.fit(x_train, y_train, epochs=20, batch_size=128)
score = model.evaluate(x_test, y_test, batch_size=128)


model.add(Conv1D(filters=2, kernel_size=1 , activation='relu', name='conv1_1'))

model.add(Conv1D(filters=64, kernel_size=5 , activation='relu', name='conv1_1'))
model.add(MaxPooling1D(2))

model.add(Conv1D(filters=64, kernel_size=5 , activation='relu', name='conv1_1'))

model.add(Conv1D(filters=2, kernel_size=2 , activation='relu', name='conv1_1'))

model.add(Bidirectional(LSTM(20,return_sequences=True),merge_mode='concat'))

model.summary()

model.add(Activation('softmax'));
'''

    
