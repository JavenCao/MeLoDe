from readtrain import *
import numpy as np
import gzip
import os
import subprocess
import tensorflow as tf
import argparse
import shutil
from tensorflow.keras.layers import Dense
from tensorflow.keras.models import Sequential
from keras.layers.convolutional import MaxPooling1D
from tensorflow.keras.callbacks import EarlyStopping
from keras.layers import Conv1D, BatchNormalization,Reshape
from keras.layers.core import Dense, Dropout, Activation, Flatten

#
if __name__ == "__main__":
    """
    usage:
    python readtrain.py --MethFeature Me.Kmat.txt.gz \
                        --UnMethFeature UnMe.Kmat.txt.gz \
                        --OutDir /output/folder/
    """
    # mfile="/storage/yangjianLab/caoyujie/project/meth/6_Freeze/1_Single_molecule_level/1_K_Feature/meHG01109_hifi/test.meHG01109_hifi.Kmat.txt.gz"
    # umfile="/storage/yangjianLab/caoyujie/project/meth/6_Freeze/1_Single_molecule_level/1_K_Feature/wgaHG01109_hifi/test.wgaHG01109_hifi.Kmat.txt.gz"
    # outdir="/storage/yangjianLab/caoyujie/project/meth/6_Freeze/1_Single_molecule_level/1_K_Feature/test"
    args = get_args()
    mfile = args.MethFeature
    umfile = args.UnMethFeature
    outdir = args.OutDir
    outdir_tmp = outdir + "/tmp"# npz files will be writen into this folder
    os.makedirs(outdir_tmp, exist_ok=True)
    strand_flag = args.StrandResolution
    cpus = int(args.cpus)
    #
    # - - - - - - - - - - - - - -
    # step1: create npz files
    winsize = Create_Read_npz_file(FetureFile = mfile, status = 1, maxZMW = 1000, outdir = outdir_tmp)
    winsize = Create_Read_npz_file(FetureFile = umfile, status = 0, maxZMW = 1000, outdir = outdir_tmp)
    #
    # - - - - - - - - - - - - - -
    # step2: create npz file list and shuffle it
    command = f"find {outdir_tmp} -name '*.npz' -printf '%p\n' | shuf > {outdir}/train.npz.file.list"
    subprocess.run(command, shell=True)
    #
    # - - - - - - - - - - - - - -
    # step3: read in train.npz.file.list and do the training process
    Train_npz_file=f"{outdir}/train.npz.file.list"
    FileBatch_Size = 3 # each time use 3*1000=3000 ZMWs to train
    zmwmeth_model = ZMWmeth_model_training(Train_npz_file, winsize, FileBatch_Size, strand_flag, cpus)
    #
    # - - - - - - - - - - - - - -
    # step4: saved as tf.model
    out_model_file_name = outdir + "/ZMWmeth.model"
    tf.saved_model.save(zmwmeth_model, out_model_file_name)
    #
    # - - - - - - - - - - - - - -
    # step5: remove outdir_tmp and train_npz_file
    if os.path.exists(outdir_tmp):
        shutil.rmtree(outdir_tmp)
    if os.path.exists(train_npz_file):
        os.remove(train_npz_file)
else:
    pass
#
#
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-Meth", "--MethFeature", help="Fully-Methylated sample")
    parser.add_argument("-UnMeth", "--UnMethFeature", help="Fully-unMethylated sample")
    parser.add_argument("-strand","--StrandResolution", action='store_true', default=False)
    parser.add_argument("-c","--cpus",default=1,help="cpus")
    parser.add_argument("-out", "--OutDir", help="Output directory")
    return parser.parse_args()
#
#
def Create_Read_npz_file(FetureFile, status, maxZMW, outdir):
    """
    2024-01-11
    Create npz file for read-level training
    """
    #
    i_c = -1 # Count of ZMW
    i_c_file = 0 # File index for the output ZMWs
    #
    # with gzip.open(FetureFile,"rt") as fp:
    with open(FetureFile,"rt") as fp:
        # sacrifice the 1st line to get the window size for initiate
        example_line = next(fp)
        winSize =  len(example_line.rstrip().split("\t")[5])# 21
        x = np.ndarray((maxZMW, 13, winSize), dtype = np.float16)
        y = [status] * maxZMW
        #
        for line in fp:
            i_c = i_c + 1
            _,_,_,fn,rn,F_template_seq, F_ipdlist, F_pwlist, R_template_seq, R_ipdlist, R_pwlist = line.rstrip().split("\t")
            npasses = int(fn) + int(rn)
            mat = Format_closedZMW(F_template_seq, F_ipdlist.split(','), F_pwlist.split(','), R_template_seq, R_ipdlist.split(','), R_pwlist.split(','), npasses)
            x[i_c] = mat
            if i_c >= maxZMW-1:
                outputfile_name = f"{outdir}/status_{status}.Feature.ZMW.chunk{i_c_file}.npz"
                np.savez_compressed(outputfile_name, X=x, y=np.array(y))
                i_c = 0
                i_c_file = i_c_file + 1
                x = np.zeros((maxZMW, 13, winSize), dtype=np.float16)
                y = [status] * maxZMW
        if i_c > 0:
            outputfile_name = f"{outdir}/status_{status}.Feature.ZMW.chunk{i_c_file}.npz"
            np.savez_compressed(outputfile_name, X=x[:i_c], y=np.array(y))
    return winSize
#
#
def Format_closedZMW(F_template_seq, F_ipdlist, F_pwlist, R_template_seq, R_ipdlist, R_pwlist, npasses):
    """Freeze,2024-01-04"""
    F_seq_context = list(F_template_seq)
    F_ipdlist = np.asarray(F_ipdlist, dtype=np.float16)
    F_pwlist = np.asarray(F_pwlist, dtype=np.float16)
    #
    R_seq_context = list(R_template_seq)
    R_ipdlist = np.asarray(R_ipdlist, dtype=np.float16)
    R_pwlist = np.asarray(R_pwlist, dtype=np.float16)
    #
    total_win_size = len(F_ipdlist)
    npasses_list = np.full(total_win_size, npasses)
    #
    F_seq_mat_array = np.zeros((4, total_win_size))
    R_seq_mat_array = np.zeros((4, total_win_size))
    #
    seq_dict = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    for index in range(total_win_size):
        F_seq_mat_array[seq_dict[F_seq_context[index]], index] = 1
        R_seq_mat_array[seq_dict[R_seq_context[index]], index] = 1
    mat = np.vstack([F_seq_mat_array, F_ipdlist, F_pwlist, R_seq_mat_array, R_ipdlist, R_pwlist, npasses_list])
    return mat
#
#
def ZMWmeth_model_training(Train_npz_file, winsize, FileBatch_Size, strand_flag, cpus):
    # 1. read in all the npz files
    with open(Train_npz_file,"rt") as f:
        file_list = f.read().splitlines()
    #
    # 2. split into training and validation
    total_file_line_num = len(file_list)
    train_file_list = file_list[0:int(total_file_line_num*0.9)]
    validation_file_list = file_list[(int(total_file_line_num*0.9)+1):]
    #
    # 3. create train_dataset and test_dataset, which are generators
    train_dataset = Generator(train_file_list, FileBatch_Size, strand_flag)
    test_dataset = Generator(validation_file_list, FileBatch_Size, strand_flag)
    #
    if not strand_flag:
        model = ZMW_model(winSize=winsize, FeatureSize=13)
    if strand_flag:
        model = ZMW_model(winSize=winsize, FeatureSize=7)
    early_stopping = EarlyStopping(patience=20, monitor='val_loss', restore_best_weights=True)
        # Train the model
    model.fit(train_dataset,
            steps_per_epoch=len(train_file_list)//FileBatch_Size,
            epochs=70,
            validation_data=test_dataset,
            validation_steps=len(validation_file_list)//FileBatch_Size,
            verbose=2,
            callbacks=[early_stopping, CheckNanCallback(train_dataset)],
            use_multiprocessing=True, workers=cpus,
            max_queue_size=10000)
    return model
#
#
class Generator(tf.keras.utils.Sequence):
    # Class is a dataset wrapper for better training performance
    def __init__(self, file_list, batch_size, strand_flag):
        self.file_list = file_list
        self.batch_size = batch_size
        self.indices = list(range(len(file_list)))
        self.strand = strand_flag
    #
    def __len__(self):
        return len(self.file_list) // self.batch_size
    #
    def on_epoch_end(self):
        np.random.shuffle(self.indices)
    #
    def __getitem__(self, idx):
        inds = self.indices[idx * self.batch_size:(idx + 1) * self.batch_size]
        data = np.load(self.file_list[inds[0]], allow_pickle=True)
        batch_x = data['X']
        batch_y = data['y']
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
            batch_y = np.concatenate((batch_y, data['y']))
        #
        if not self.strand:
            # final return is in the format of (N,21,13)
            return np.transpose(batch_x,(0,2,1)), batch_y
        else:
            X = batch_x#(10000,13,21)
            y = batch_y#(N,)
            N = X.shape[0]# 10000
            N_row = X.shape[1]# 13 for closedZMW
            N_col = X.shape[2]# 21 winSize
            X_strand = np.ndarray(shape=(2*N, 7, N_col))# X_strand.shape = (200000, 7, 21)
            y_strand = np.ndarray(shape=(2*N,))
            for index in range(0,N):
                Forward_strand_mat = np.vstack([X[index][0:6],X[index][N_row-1]])
                Reverse_strand_mat = np.vstack([X[index][6:12],X[index][N_row-1]])
                X_strand[2*index] = Forward_strand_mat
                X_strand[2*index+1] = Reverse_strand_mat
                y_strand[2*index] = y[index]
                y_strand[2*index+1] = y[index]
            # final return is in the format of (2N,21,7)
            return np.transpose(X_strand,(0,2,1)),y_strand
#
#
class CheckNanCallback(tf.keras.callbacks.Callback):
    def __init__(self, generator):
        super().__init__()
        self.generator = generator
    def on_train_batch_end(self, batch, logs=None):
        if logs is not None and np.isnan(logs.get('loss')):
            current_file = self.generator.file_list[self.generator.indices[batch*self.generator.batch_size]]
            print(f'NaN encountered in {current_file} batch {batch}')


# 2023-05-09, after 20 trails
# Trial 18 summary
# Hyperparameters:
# num_conv_layers: 5
# conv_filters: 256
# kernel_size: 5
# dropout_rate: 0.2
# num_dense_layers: 2
# dense_0_units: 128
# dense_1_units: 32
# Score: 0.8616424202919

def ZMW_model(winSize=21, FeatureSize=13):
    #
    num_conv_layers = 5
    conv_filters = 256
    kernel_size = 5
    dropout_rate = 0.2
    #
    # Define input shape
    seq_len = winSize
    input_features = FeatureSize
    origin = tf.keras.layers.Input((seq_len, input_features))
    #
    conv_input = origin[:,:,:(input_features-1)]
    for i in range(num_conv_layers):
        if i == 0:
            conv = tf.keras.layers.Conv1D(filters=conv_filters, kernel_size=kernel_size, input_shape=(seq_len, input_features-1), padding='same', activation='relu')(conv_input)
        else:
            conv = tf.keras.layers.Conv1D(filters=conv_filters, kernel_size=kernel_size, padding='same', activation='relu')(conv)
        conv = tf.keras.layers.BatchNormalization()(conv)
    #
    drop1 = tf.keras.layers.Dropout(dropout_rate)(conv)
    pool2 = tf.keras.layers.MaxPooling1D(pool_size=2)(drop1)
    #
    flat = tf.keras.layers.Flatten()(pool2)
    depth = origin[:, 0, input_features-1]
    depth2 = tf.keras.layers.Flatten()(depth)
    concat = tf.keras.layers.concatenate([flat, depth2], axis=1)
    #
    dense = concat
    dense = tf.keras.layers.Dense(128, activation='relu')(dense)
    dense = tf.keras.layers.Dense(32, activation='relu')(dense)
    #
    drop2 = tf.keras.layers.Dropout(dropout_rate)(dense)
    dense3 = tf.keras.layers.Dense(1, activation='sigmoid')(drop2)
    # Compile Model
    model = tf.keras.Model(origin, dense3)
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    return model
