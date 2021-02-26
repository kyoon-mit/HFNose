import sys
import argparse
from types import SimpleNamespace
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.callbacks import EarlyStopping
# Taken from https://github.com/leimao/Frozen_Graph_TensorFlow/blob/master/TensorFlow_v2/train.py
from tensorflow.python.framework.convert_to_constants import convert_variables_to_constants_v2
import os
import numpy as np
import pandas as pd
from glob import glob
from math import sqrt
import json    
from utils import ParserWrapper

#Read https://stackoverflow.com/questions/58255821/how-to-use-k-get-session-in-tensorflow-2-0-or-how-to-migrate-it
from tensorflow.python.keras import backend as K

from preprocessing_pb import checkDirAndCreate

def training(args):
    if not isinstance(args, argparse.Namespace):
        args = SimpleNamespace(**args)

    # If you want to reserve only the 30% of the GPU memory uncomment the next 4 lines
    # from keras.backend.tensorflow_backend import set_session
    # config = tf.ConfigProto()
    # config.gpu_options.per_process_gpu_memory_fraction = 0.3
    # set_session(tf.Session(config=config))
    epochs = 15

    pwrap = ParserWrapper(norm_file=None, task=args.task)
    meanen, stden = pwrap.get_norm_params()
        
    # Taken from https://leimao.github.io/blog/Save-Load-Inference-From-TF2-Frozen-Graph/
    # See also https://github.com/leimao/Frozen_Graph_TensorFlow/blob/master/TensorFlow_v2/example_2.py#L86-L87 for handling multiple inputs
    def _saveModel(model, args):
        # Convert Keras model to ConcreteFunction
        full_model = tf.function(lambda x: model(x))
        inputs=None
        if args.sigmas:
            inputs = ( tf.TensorSpec(model.inputs[0].shape, model.inputs[0].dtype, name="input"),
                       tf.TensorSpec(model.inputs[1].shape, model.inputs[1].dtype, name="input_trackster") )
        else:
            inputs = ( tf.TensorSpec(model.inputs[0].shape, model.inputs[0].dtype, name="input"), )
        full_model = full_model.get_concrete_function(x=inputs)

        # Get frozen ConcreteFunction
        frozen_func = convert_variables_to_constants_v2(full_model)
        frozen_func.graph.as_graph_def()

        if args.debug:
            layers_debug = [op.name for op in frozen_func.graph.get_operations()]
            print("-" * 50)
            print("Frozen model layers: ")
            for layer in layers_debug:
                print(layer)

            print("-" * 50)
            print("Frozen model inputs: ")
            print(frozen_func.inputs)
            print("Frozen model outputs: ")
            print(frozen_func.outputs)

        frozenFolder = os.path.join(args.outputDir, args.frozenModelFolder)
        checkDirAndCreate( frozenFolder )
        graph=frozen_func.graph
        with graph.as_default():
            if args.task != 'er':
                id_op = graph.get_operation_by_name("model/pid_output/Softmax")
                assert(id_op)
                tf.identity(id_op.outputs[0], name="output/id_probabilities")
            if args.task != 'pid':
                er_op = graph.get_operation_by_name("model/regressed_energy/add")
                assert(er_op)
                tf.identity(er_op.outputs[0], name="output/regressed_energy")
        # Save frozen graph from frozen ConcreteFunction to hard drive
        tf.io.write_graph(graph_or_graph_def=graph,
                logdir = frozenFolder,
                name = args.frozenModelName,
                as_text=False)
        tf.io.write_graph(graph_or_graph_def=graph,
                logdir = frozenFolder,
                name = args.frozenModelName.replace('.pb','.txt'),
                as_text=True)

        name = os.path.join(frozenFolder, args.frozenModelName)
        with open(args.targetsFile, 'a') as f:
            f.write(name + '\n')
        name = name.replace('.pb','.txt')
        with open(args.targetsFile, 'a') as f:
            f.write(name + '\n')

    # Start the real event-processing and training
    batch_size = 512
    counts = json.loads(open("{}/counts.json".format(args.inputDir),"r").read())
    steps_per_epoch =  counts["train"] // batch_size
    eval_steps_per_epoch = counts["eval"] // batch_size

    # Create datasets from TFRecord files.
    files = []
    files.extend(glob("{}/train-*".format(args.inputDir)))

    if args.debug:
        print(files)

    train_ds = tf.data.TFRecordDataset(files)
    eval_ds = tf.data.TFRecordDataset(glob('{}/eval-*'.format(args.inputDir)))
    if args.algorithm == 'nn':
        train_ds = train_ds.map(pwrap.parse_example_training)
        train_ds = train_ds.batch(batch_size).repeat()
        eval_ds = eval_ds.map(pwrap.parse_example_training)
        eval_ds = eval_ds.batch(batch_size)

    if args.debug:
        # Read a single batch of examples from the training set and display shapes.
        for features_in, features_out in train_ds:
          break
        print('img.shape (batch_size, image_height, image_width, channels) =',
              features_in['input'].shape)
        if args.sigmas:
            print('trackster.shape (batch_size, energy, sigma1, sigma2, sigma3) =',
                  features_in['input_trackster'].shape)
            print("RawEnergy:", features_in['input_trackster'][0])
            print("Sigmas:", features_in['input_trackster'][1:])
        print("Energy:", features_in['input'][0,:,:,0])
        print("Eta:", features_in['input'][0,:,:,1])
        print("Phi:", features_in['input'][0,:,:,2])
        if args.task != 'er':
            print(features_out['pid_output'])
        if args.task != 'pid':
            print(features_out['regressed_energy'])

    if args.algorithm=='nn':
        model = deepLearningModel_(energyNorm=(meanen, stden), nclasses=pwrap.nclasses,
                                   epochs=epochs, pcaInput=args.sigmas, shower_window=args.shower_window,
                                   task=args.task)
        model.summary()

        early_stopping = EarlyStopping(monitor='loss', patience=5)
        history = model.fit(train_ds,
                            validation_data=eval_ds,
                            steps_per_epoch=steps_per_epoch,
                            validation_steps=eval_steps_per_epoch,
                            epochs=epochs,
                            callbacks=[early_stopping],
                            verbose=True) #, validation_split=0.1,  shuffle=True, verbose=1)

        name = os.path.join(args.modelFolder, args.modelHistory)
        checkDirAndCreate( args.modelFolder )

    elif args.algorithm == 'svr':
        train_ds = train_ds.cache()
        train_ds = list(train_ds.as_numpy_iterator())
        print(train_ds[0])
        #decode TFRecordDataset: https://stackoverflow.com/questions/47878207/how-to-read-tfrecord-data-into-tensors-numpy-arrays
        quit()
        model = SVRModel_(x, y)
        model.fit()
        
    history_save = pd.DataFrame(history.history).to_hdf(name, "history", append=False)
    with open(args.targetsFile, 'a') as f:
        f.write(name + '\n')
    name = os.path.join(args.modelFolder, args.modelName)
    model.save(name)
    with open(args.targetsFile, 'a') as f:
        f.write(name + '\n')

    _saveModel(model, args)

def SVRModel_(x, y):
    from sklearn.svm import SVR
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    return make_pipeline(StandardScaler(), SVR(C=1.0, epsilon=0.2))

def deepLearningModel_(energyNorm, nclasses, epochs, pcaInput=False, shower_window=False, task='all'):
    width, height, channels, epochs = 50, 10, 3, epochs

    conv_base = dict(activation='relu', padding='same', kernel_initializer='glorot_uniform')
    dens_base = dict(activation='relu', kernel_initializer='glorot_uniform')

    def rescaleEnergy_(x, mean=energyNorm[0], std=energyNorm[1]):
        return x * std + mean

    def batchNorm_(l, axis):
        l = tf.keras.activations.relu(l)
        l = layers.BatchNormalization(axis=axis)(l)
        return tf.keras.layers.Dropout(0.5)(l)

    input_img = keras.Input(shape=(width, height, channels), name='input')
    if pcaInput:
        input_trackster = keras.Input(shape=(4), name='input_trackster')

    enTensor, etaTensor, phiTensor = input_img[:,:,:,0], input_img[:,:,:,1], input_img[:,:,:,2]
    enTensor_bn, etaTensor_bn, phiTensor_bn = batchNorm_(enTensor, axis=[1,2]), batchNorm_(etaTensor, axis=[1,2]), batchNorm_(phiTensor, axis=[1,2])
    input_img_bn = tf.stack([enTensor_bn, etaTensor_bn, phiTensor_bn], axis=3)#batchNorm_(input_img, axis=3)
        
    """
    enTensor, etaTensor, phiTensor = input_img[:,:,:,0], input_img[:,:,:,1], input_img[:,:,:,2]
    enSum = tf.reduce_sum(enTensor, axis=2)

    #energy layer fractions
    enSumEE = tf.reduce_sum( tf.slice(enSum, begin=[0,0], size=[tf.shape(enSum)[0],28]), keepdims=True, axis=1 )
    enSumEEFraction = tf.math.divide_no_nan( enSumEE, tf.reduce_sum(enSum, keepdims=True, axis=1) )
    enSumFirstLayer = tf.reduce_sum( tf.slice(enSum, begin=[0,0], size=[tf.shape(enSum)[0],1]), keepdims=True, axis=1 )
    enSumFirstLayerFraction = tf.math.divide_no_nan( enSumFirstLayer, tf.reduce_sum(enSum, keepdims=True, axis=1) )
    
    enMaxima = tf.reduce_max(enTensor, axis=2)
    enMaximaFraction = tf.math.divide_no_nan(enMaxima, enSum) #fraction of the energy of the largest cluster per layer
    
    if shower_window: #throws away data far from the shower maximum
        dist = (12,8)
        enShape = tf.shape(enTensor)
        showerMaxIdx = tf.clip_by_value( tf.cast(tf.argmax(enSum, axis=1), dtype=tf.int32), dist[0], width-dist[1] )
        distSum = dist[0] + dist[1]
        distBelow = showerMaxIdx - tf.expand_dims(tf.convert_to_tensor(dist[0]), axis=0)
        distAbove = showerMaxIdx + tf.expand_dims(tf.convert_to_tensor(dist[1]), axis=0)

        batch_pile = []
        for i in range(distSum):
            batch_pile.append( tf.range(start=0, limit=enShape[0]) )
        batch_pile = tf.stack(batch_pile, axis=1)

        layer_pile = []
        for i in range(distSum):
            layer_pile.append( distBelow + i*tf.ones_like(distBelow) )
        layer_pile = tf.stack(layer_pile, axis=1)

        indices = tf.concat([tf.expand_dims(batch_pile, axis=2), tf.expand_dims(layer_pile, axis=2)], axis=2)

        #filter inputs around the shower maximum
        enTensor = tf.gather_nd(enTensor, indices, batch_dims=0)
        etaTensor = tf.gather_nd(etaTensor, indices, batch_dims=0)
        phiTensor = tf.gather_nd(phiTensor, indices, batch_dims=0)
        enMaximaFraction = tf.gather_nd(enMaximaFraction, indices, batch_dims=0)

        enSumEEFraction = tf.tile(enSumEEFraction, multiples=[1,distSum])
        enSumFirstLayerFraction = tf.tile(enSumFirstLayerFraction, multiples=[1,distSum])
        enSum = tf.reduce_sum(enTensor, axis=2)
    else:
        enSumEEFraction = tf.tile(enSumEEFraction, multiples=[1,width])
        enSumFirstLayerFraction = tf.tile(enSumFirstLayerFraction, multiples=[1,width])

    #layer cluster multiplicity
    lcMult = tf.math.count_nonzero(enTensor, axis=2)

    #weighted mean and variance for eta and phi values
    etaMean = tf.math.divide_no_nan( tf.reduce_sum( enTensor*etaTensor, axis=2 ), enSum )
    phiMean = tf.math.divide_no_nan( tf.reduce_sum( enTensor*phiTensor, axis=2 ), enSum )
    enSumSquared = tf.reduce_sum(tf.math.square(enTensor), axis=2)
    varianceDenominator = enSum - tf.math.divide_no_nan(enSumSquared, enSum)
    etaVariance = tf.math.divide_no_nan( tf.reduce_sum( enTensor * ( etaTensor - tf.expand_dims(etaMean, axis=2) )**2, axis=2),
                                         varianceDenominator )
    phiVariance = tf.math.divide_no_nan( tf.reduce_sum( enTensor * ( phiTensor - tf.expand_dims(phiMean, axis=2) )**2, axis=2),
                                         varianceDenominator )

    enSum = (enSum-energyNorm[0])/energyNorm[1]
    inputTensors = [enSum, tf.cast(lcMult, dtype=tf.float32), enMaximaFraction, etaMean, etaVariance, phiMean, phiVariance, enSumEEFraction, enSumFirstLayerFraction]
    if shower_window:
        inputTensors.append(tf.cast(layer_pile, dtype=tf.float32))
    inputTensor = tf.stack(inputTensors, axis=2)
    """
    #conv = layers.Conv1D(filters=32, kernel_size=5, name='conv1', **conv_base)(inputTensor)
    #conv = layers.Conv1D(filters=24, kernel_size=3, data_format='channels_last', name='conv2', **conv_base)(conv)
    #conv = layers.Conv1D(filters=16, kernel_size=3, data_format='channels_last', name='conv3', **conv_base)(conv)

    conv = layers.Conv2D(filters=16, kernel_size=(5,1), name='conv1', **conv_base)(input_img_bn)
    conv = layers.Conv2D(filters=16, kernel_size=(3,3), data_format='channels_last', name='conv2', **conv_base)(conv)
    conv = layers.Conv2D(filters=16, kernel_size=(3,3), data_format='channels_last', name='conv3', **conv_base)(conv)

    flat = layers.Flatten()(conv)
    if pcaInput:
        flat = layers.Concatenate()([flat, input_trackster])

    dense = layers.Dense(512, name='dense1', **dens_base)(flat)
    dense = layers.Dense(128, name='dense2', **dens_base)(dense)
    #dense = tf.concat([dense, enSumEEFraction, enSumFirstLayerFraction], axis=1)

    if task != 'er':
        dense_id = layers.Dense(64, name='dense_id1', **dens_base)(dense)
        dense_id = layers.Dense(16, name='dense_id2', **dens_base)(dense_id)
        pid = layers.Dense(nclasses, activation='softmax', name='pid_output')(dense_id)

    if task != 'pid':
        """
        pid_argmax = layers.Lambda(lambda x: tf.expand_dims( tf.cast(tf.argmax(x, axis=1), tf.float32), axis=1 ),
        name='pid_output_argmax', trainable=False)(pid)
        remix = layers.Concatenate()([dense, pid_argmax])
        """
        dense_er = layers.Dense(64, name='dense_er1', **dens_base)(dense)
        dense_er = layers.Dense(16, name='dense_er2', **dens_base)(dense_er)
        enreg = layers.Dense(1, name='regressed_energy_norm')(dense_er)
        enreg_rescaled = layers.Lambda(rescaleEnergy_, name='regressed_energy', trainable=False)(enreg)
        
        if energyNorm[0]==0 and energyNorm[1]==1:
            metric = tf.keras.losses.MeanAbsolutePercentageError()
        else:
            metric = tf.keras.losses.MeanAbsoluteError()

    if task == 'pid':
        model_outputs = [pid]
        loss = {'pid_output': 'categorical_crossentropy'}
        loss_weights = {'pid_output': 1.}
        metrics = {'pid_output': 'accuracy'}
    elif task == 'er':
        model_outputs = [enreg, enreg_rescaled]
        loss = {'regressed_energy_norm': metric}
        loss_weights = {'regressed_energy_norm': 1.}
        metrics = {'regressed_energy_norm': metric}
    else:
        model_outputs = [pid, enreg, enreg_rescaled]
        loss = {'pid_output': 'categorical_crossentropy', 'regressed_energy_norm': metric}
        loss_weights = {'pid_output': 1., 'regressed_energy_norm': 1.}
        metrics = {'pid_output': 'accuracy', 'regressed_energy_norm': metric}
            
    if pcaInput:
        model = keras.Model(inputs=[input_img, input_trackster], outputs=model_outputs)
    else:
        model = keras.Model(inputs=[input_img], outputs=model_outputs)

    model.compile(loss=loss, loss_weights=loss_weights, metrics=metrics,
                  optimizer=keras.optimizers.Adam(learning_rate=0.001))

    return model


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train PID and Energy regression on tf.Example/ProtocolBuffer files")
    parser.add_argument("-i", "--inputDir", help="Input directory where to search for input pb files.", required=True)
    parser.add_argument("-o", "--outputDir", help="Output directory", required=True)
    parser.add_argument("--normFile", help="File containing mean and std of the energy for normalization purposes (full path).", required=True)
    parser.add_argument("--modelName", help="name of the output model [.h5]", required=True)
    parser.add_argument("--modelFolder", help="path to the output model's folder [.h5]", required=True)
    parser.add_argument("--modelHistory", help="path to the output model's history [.h5]", required=True)
    parser.add_argument("--frozenModelName", help="name of the output model [.pb]", required=True)
    parser.add_argument("--frozenModelFolder", help="path to the output model's folder [.pb]", required=True)
    parser.add_argument("-d", "--debug", help="Activate debug mode", required=True)
    parser.add_argument("-s", "--sigmas", help="Add trackster's energy and PCAs sigmas to the input variables for the model.", required=True)
    parser.add_argument("--algorithm", help="Algorithm to perform the energy regression.", required=True)
    parser.add_argument("--shower_window", help="Drop all the training information from layers far away from the shower maximum.", required=True)
    parser.add_argument("--task", help='Choose the task which will be trained.', choices=['all', 'pid', 'er'], required=True)
    args = parser.parse_args()

    training(args) 
