from sitetrain import *
import argparse
import numpy as np
import tensorflow as tf
from keras import models
from keras import layers
from keras.models import model_from_json
from tensorflow.keras.callbacks import EarlyStopping
#
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-npzl", "--npzlist", help="HSme npz file list")
    parser.add_argument("-r","--windowRadius", help="how many CpG sites in a sliding window? (2 * r + 1)")
    parser.add_argument("-bz", "--batchsize", help="batch size of the input training, say how many npz files will be fed to the model at once.")
    parser.add_argument("-c","--cpusnum",help="how many workers/cpus for this program")
    parser.add_argument("-o", "--outputfile")
    return parser.parse_args()
#
if __name__ == "__main__":
    """
    usage:
    python sitetrain --npzlist Sample.npz.list \
                     -r 5 \
                     -bz 5 \
                     -c 10 \
                     -o out.sitelevel.model
    """

    # -------------------------------------------------------------------------
    args = get_args()
    npzlistfile=args.npzlist
    R=int(args.windowRadius)
    Batch_Size = int(args.batchsize)
    cpus=int(args.cpusnum)
    outputfile = args.outputfile
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # 1. read in all the npz files
    with open(npzlistfile,"rt") as f:
        file_list = f.read().splitlines()
    # 2. split into training and validation
    total_file_line_num = len(file_list)
    train_file_list = file_list[0:int(total_file_line_num*0.9)]
    validation_file_list = file_list[(int(total_file_line_num*0.9)+1):]
    # 3. create train_dataset and test_dataset
    train_dataset = Generator(train_file_list, Batch_Size)
    test_dataset = Generator(validation_file_list, Batch_Size)
    #
    # -------------------------------------------------------------------------
    # 4. build model
    CpGModel = F_CpGModel(winSize=2*R+1, FeatureSize=32)
    early_stopping = EarlyStopping(patience=10, monitor='val_loss', restore_best_weights=True)
    # 5. feed data
    # fit_generator
    CpGModel.fit_generator(train_dataset,
                    steps_per_epoch=len(train_file_list)//Batch_Size,
                    epochs=100,
                    verbose=2,
                    validation_data=test_dataset,
                    validation_steps=len(validation_file_list)//Batch_Size,
                    callbacks=[early_stopping, CheckNanCallback(train_dataset)],use_multiprocessing=True, workers=cpus,max_queue_size=10000)
    # 6. saved as tf.model
    out_model_file_name = outdir + "/JointCpG.model"
    tf.saved_model.save(CpGModel, out_model_file_name)
else:
    pass

#
class Generator(tf.keras.utils.Sequence):
    # Class is a dataset wrapper for better training performance
    def __init__(self, file_list, batch_size):
        self.file_list = file_list
        self.batch_size = batch_size
        self.indices = list(range(len(file_list)))

    def __len__(self):
        return len(self.file_list) // self.batch_size

    def __getitem__(self, idx):

        inds = self.indices[idx * self.batch_size:(idx + 1) * self.batch_size]
        # print('Loading batch', idx, 'of', len(self))
        # print('Loading file', self.file_list[inds[0]])
        data = np.load(self.file_list[inds[0]], allow_pickle=True)
        batch_x = data['X']
        batch_y = data['Y']
        # remove rows with NaN values, by chatGPT
        nan_mask_y = ~np.isnan(batch_y)
        nan_mask_x = ~np.isnan(batch_x).any(axis=(1,2))
        nan_mask = np.logical_and(nan_mask_x, nan_mask_y)
        batch_x = batch_x[nan_mask]
        batch_y = batch_y[nan_mask]
        #
        for i in range(1,len(inds)):
            data = np.load(self.file_list[inds[i]], allow_pickle=True)
            batch_x = np.vstack((batch_x, data['X']))
            batch_y = np.concatenate((batch_y, data['Y']))
        return batch_x, batch_y

    def on_epoch_end(self):
        np.random.shuffle(self.indices)

class CheckNanCallback(tf.keras.callbacks.Callback):
    def __init__(self, generator):
        super().__init__()
        self.generator = generator

    def on_train_batch_end(self, batch, logs=None):
        if logs is not None and np.isnan(logs.get('loss')):
            current_file = self.generator.file_list[self.generator.indices[batch*self.generator.batch_size]]
            print(f'NaN encountered in {current_file} batch {batch}')
#
#
def F_CpGModel(winSize=11, FeatureSize=32):
    # Trial ID: 15
    # Hyperparameters: {'num_conv_layers': 2,
    # 'conv_filters': 96,
    # 'kernel_size': 4,
    # 'dropout_rate': 0.3,
    # 'num_dense_layers': 3,
    # 'dense_0_units': 30,
    # 'dense_1_units': 10,
    # 'dense_2_units': 10}
    # Score: 0.07388371229171753
    #
    num_conv_layers = 2
    conv_filters = 96
    kernel_size = 4
    dropout_rate = 0.3
    num_dense_layers = 3
    #
    seq_len = winSize
    input_features = FeatureSize
    #
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Input(shape=(seq_len, input_features)))
    #
    for i in range(num_conv_layers):
        if i == 0:
            model.add(tf.keras.layers.Conv1D(filters=conv_filters, kernel_size=kernel_size, padding='same', activation='relu'))
        else:
            model.add(tf.keras.layers.Conv1D(filters=conv_filters, kernel_size=kernel_size, padding='same', activation='relu'))
        model.add(tf.keras.layers.BatchNormalization())
    #
    model.add(tf.keras.layers.Dropout(dropout_rate))
    model.add(tf.keras.layers.MaxPooling1D(pool_size=2))
    model.add(tf.keras.layers.Flatten())
    #
    model.add(tf.keras.layers.Dense(30, activation='relu'))
    model.add(tf.keras.layers.Dense(10, activation='relu'))
    model.add(tf.keras.layers.Dense(10, activation='relu'))
    #
    model.add(tf.keras.layers.Dense(1, activation='linear'))
    model.compile(optimizer='sgd', loss=tf.keras.losses.MeanAbsoluteError(),metrics=[tf.keras.metrics.MeanAbsoluteError()])
    return model
