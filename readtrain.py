from readtrain import *
import gzip
import sklearn
import argparse
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow.keras.layers import Dense
from tensorflow.keras.models import Sequential
from keras.layers.convolutional import MaxPooling1D
from tensorflow.keras.callbacks import EarlyStopping
from keras.layers import Conv1D, BatchNormalization,Reshape
from keras.layers.core import Dense, Dropout, Activation, Flatten
"""
python readtrain.py --MethFeature Me.Kmat.txt.gz \
                    --UnMethFeature UnMe.Kmat.txt.gz \
                    --Tempdir /temp/folder \
                    --OutModel Trained.zmwmeth.model
"""
def get_args():
    parser = argparse.ArgumentParser()


# top
Top m64284e_230525_074225_8 176     12      11      GTTTTCGGGGCCCGACCCTTA   0.54,-0.76,2.48,0.15,-0.63,-0.76,-0.37,0.41,2.22,0.8,0.28,0.15,0.54,-0.63,-0.89,0.9>

chr1    12265   m64284e_230525_074225_108202546 15      16      CGACACTATGCACCGGCCGGG   -0.62,-0.5,-0.74,-0.74,-0.46,-0.46,-0.54,-0.05,0.07,-0.01,-0.54,0.0>


# step 1 :
