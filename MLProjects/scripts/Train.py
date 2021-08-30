#!/usr/bin/env python3

# Original Repo: https://gitlab.cern.ch/gouskos/hgcal_reco_analysis/HGCMLAnalyzer/tools/pid/training.py
# TODO: more work needs to be done

import argparse
import os

import tensorflow as tf
import uproot
import numpy as np
import pandas as pd
from glob import glob
from tqdm import tqdm
  
def create_input_images (given_df, max_layers=50, max_layerclusters=10):
    """
    Create input images to feed into NN, given a dataframe containing events.    

    """
    
    # Get number of events
    unique_entries = given_df.index.get_level_values('entry').unique().to_list()
    nevents = len(unique_entries)
    
    # Keep only relevant variables
    trimmed_df = pd.DataFrame({'layer': given_df.lc_layer,
                              'energy': given_df.lc_energy,
                              'eta':    given_df.lc_eta,
                              'phi':    given_df.lc_phi})
    trimmed_df = trimmed_df.sort_values(['layer', 'energy'],
                                       ascending=[True, False]).reset_index(drop=True)

    # Procedure to nicely arrange indices
    idx = trimmed_df.groupby('layer').indices.values()
    idx_tmp = []
    for sublist in idx:
        for item in sublist[:min(len(sublist), max_layerclusters)]:
            idx_tmp.append(item)
    idx = np.array(idx_tmp)
    trimmed_df = trimmed_df.iloc[idx]
    global_idx = []
    for l, lc in trimmed_df.groupby("layer").indices.items():
        global_idx.extend([j + (l-1)*max_layerclusters for j in range(len(lc))])
    trimmed_df = trimmed_df.set_index(np.array(global_idx).astype(int))
    
    # Set layer number correctly
    l = [layer for layer in range(1, max_layers+1) for layercluster in range(max_layerclusters)]
    dummy = [0] * max_layers * max_layerclusters
    padded = pd.DataFrame({'layer':  l,
                           'energy': dummy,
                           'eta':    dummy,
                           'phi':    dummy})
    padded.iloc[trimmed_df.index] = trimmed_df
    
    training_img = np.array(padded.drop(columns='layer')).reshape(50,10,3)
    energy_label = given_df.cp_energy[0]
    pdgid_label  = abs(given_df.cp_pdgid[0])
    
    return training_img, energy_label, pdgid_label


def preprocess_rootfile (filepath, max_missingEnergyFraction=0.45, max_layers=50, max_layerclusters=10):
    """
    Preprocess single or multiple rootfile(s).
    Returns: numpy array
    """
    # Variables to keep in dataframe
    # keepvars = ['ts_energy', 'ts_sigma1', 'ts_sigma2', 'ts_sigma3',
    #             'cp_missingEnergyFraction', 'cp_energy', 'cp_pdgid',
    #             'lc_energy', 'lc_eta', 'lc_phi', 'lc_layer']
    keepvars = ['cp_missingEnergyFraction', 'cp_energy', 'cp_pdgid',
                'lc_energy', 'lc_eta', 'lc_phi', 'lc_layer']
    
    rootfiles = glob(filepath)
    
    df_list = []
    for rootfile in tqdm(rootfiles):
        filename = os.path.basename(rootfile).replace(".root", "")
        try:
            # if file size too large, read in chunks
            with uproot.open(filename) as infile:
                tree = infile['HFNoseNtuple/pidtree']
                
                # Convert root file to dataframe
                df = tree.arrays(filter_name=keepvars, library='pd')

                # Keep events with (missing rechit energy)/(caloparticle energy)
                # below given threshold
                df = df[df['cp_missingEnergyFraction'] < max_missingEnergyFraction]
                df.drop(['cp_missingEnergyFraction'], axis=1)   
        except:
            print('File {} had a problem.'.format(rootfile))
            raise

    output_df = pd.concat(df_list)
    
    # TODO: Need to shuffle and batch normalize
    
    return output_df
    

def create_model ():
    """
    Function to create model.
    """
    conv_base = dict(activation='relu', padding='same', kernel_initializer='glorot_uniform')
    dens_base = dict(activation='relu', kernel_initializer='glorot_uniform')

    # CNN part of the model
    model = models.Sequential()
    model.add(layers.Conv2D(16, (5, 1), input_shape=(50, 10, 3), name='conv1', **conv_base))
    model.add(layers.Conv2D(16, (3, 3), data_format='channels_last', name='conv2', **conv_base))
    model.add(layers.Conv2D(16, (3, 3), data_format='channels_last', name='conv3', **conv_base))

    model.add(layers.Flatten())
    model.add(layers.Dense(512, name='dense1', **dens_base))
    model.add(layers.Dense(128, name='dense2', **dens_base))

    # Create energy regression part
    ereg = model
    ereg.add(layers.Dense(64, name='dense_er1', **dens_base))
    ereg.add(layers.Dense(16, name='dense_er2', **dens_base))
    ereg.add(layers.Dense(1, name='regressed_energy_norm'))
    
    # TODO: add PID part
    
    ereg.compile(loss=tf.keras.losses.CategoricalCrossentropy(),
                 loss_weights=[1.],
                 metrics=['accuracy'],
                 optimizer=keras.optimizers.Adam(learning_rate=0.001))
                 
    return ereg


# TODO: define training and testing functions
def train ():
    return
    

def test ():
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # TODO: add stuff to parser
    args = parser.parse_args()
    preprocess(args)
