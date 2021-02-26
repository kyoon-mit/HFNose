import argparse
import sys
import uproot
import os, argparse
from types import SimpleNamespace
import numpy as np
import pandas as pd
from glob import glob
from tqdm import tqdm
from utils import make_example_pid_and_pad
from time import time
import tensorflow as tf

def checkDirAndCreate(d):
    if not os.path.exists(d):
        os.mkdir(d)

def preprocess(args):
    """ Preprocess the flat ROOT files in input

    The input ROOT files are organized as follows:

    - run
    - event
    - trackster_id
    - [layers   of the layerClusters belonging to this trackster]
    - [etas     of the layerClusters belonging to this trackster]
    - [phis     of the layerClusters belonging to this trackster]
    - [energies of the layerClusters belonging to this trackster]
    - energy of the caloParticle linked to this trackster
    - pdgId of the caloParticle linked to this trackster
    - raw_energy of the trackster

    FLAT PROTOCOL BUFFER PADDED FORMAT

    This is the format that will be used in input to the training. It should be structured as follows:

    - trackster number
    - energy of the caloParticle
    - pdgId of the caloParticle
    - raw_energy of the trackster
    - sigma1 of the trackster (along the "major" PCA component)
    - sigma2 of the trackster
    - sigma3 of the trackster
    - layer of a layerCluster
    - energy of a layerCluster
    - eta of a layerCluster
    - phi of a layerCluster

    This structure is padded in the sense that, for every trackster, there will be 10 layerClusters for each layer, for all 50 layers of HGCAL. The target final size of the DataFrame is:

    rows of target DataFrame = number of tracksters in input ROOT files * 50 * 10

    The number of input trackster in the input ROOT files has to be deduced by the ROOT file itself. Something along the lines:

    df.shape[0]

    done on the input ROOT file.
    """
    if not isinstance(args, argparse.Namespace):
        args = SimpleNamespace(**args)
    
    if args.debug:
        pd.set_option('display.max_rows', 200)

    rootfiles_path = os.path.join(args.inputDir, args.suffix + '*.root')
    rootfiles = glob(rootfiles_path)
    if len(rootfiles) == 0:
        print('Input directory: {} Input suffix: {}'.format(args.inputDir, args.suffix))
        raise ValueError('[preprocessing_pb.py]: No input files found in {}.'.format(rootfiles, rootfiles_path))

    max_perlayer = 10
    number_layers = 50

    keepVars = ['ts_energy', 'ts_sigma1', 'ts_sigma2', 'ts_sigma3', 'cp_missingEnergyFraction']
    keepVars += ['cp_energy', 'cp_pdgid', 'lc_energy', 'lc_eta', 'lc_phi', 'lc_layer']
    
    for rootfile in tqdm(rootfiles):
        start_time = time()
        filename = os.path.basename(rootfile).replace(".root", "")
        if args.debug:
            print('File: ', rootfile, filename)
        try:
            with uproot.open(rootfile) as open_file:
                directory = open_file[args.dir]
                tree = directory[args.tree]
                df = tree.arrays(filter_name=keepVars, library='pd')
                df = df[df['cp_missingEnergyFraction'] < args.maxMissingEnergyFraction]
                df.drop(['cp_missingEnergyFraction'], axis=1)
        except:
            print('File {} had a problem.'.format(rootfile))
            raise

        unique_entries = df.index.get_level_values('entry').unique().to_list()
        events = len(unique_entries)
        if args.debug:
            print(df)

        checkDirAndCreate( args.outputDir )
        checkDirAndCreate( os.path.join(args.outputDir, 'padded') )

        name = os.path.join(args.outputDir, 'padded', filename + '_padded.pb')
        writer = tf.io.TFRecordWriter(name)
        with open(args.targetsFile, 'a') as f:
            f.write(name + '\n')

        for entry in unique_entries:
            example = make_example_pid_and_pad(df.iloc[[entry]], number_layers, max_perlayer, args.debug)
            writer.write(example.SerializeToString())
            if args.debug:
                print(example)

        end_time = time()
        name = name.replace('.pb', '.log')
        with open(name, 'w') as logfile:
            logfile.write('Events: {} Time: {} Rate: {}'.format(events, (end_time-start_time), events/float(end_time-start_time)))
        with open(args.targetsFile, 'a') as f:
            f.write(name + '\n')

    print("Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PreProcessing Flat ROOT Ntuple in H5 format for ML")
    parser.add_argument("-i", "--inputDir", help="Input directory where to search for input ROOT files", required=True)
    parser.add_argument("-o", "--outputDir", help="Output directory", required=True)
    parser.add_argument("--targetsFile", help="Directory for luigi target files (outputs)", required=True)
    parser.add_argument("-s", "--suffix", help="Suffix to identify a sample", required=True)
    parser.add_argument("-d", "--dir", help="Directory to look for tree in input ROOT file", default="pid")
    parser.add_argument("-t", "--tree", help="Name of the tree to look for in input ROOT file", default="tree")
    parser.add_argument("-g", "--debug", help="Activate debug mode", action='store_true')
    parser.add_argument("-m", "--maxMissingEnergyFraction", help="Maximum value for the missing energy of a CaloParticle wrt the energy of the associated Trackster.", required=True, type=float, default=1.)
    args = parser.parse_args()
    preprocess(args)
