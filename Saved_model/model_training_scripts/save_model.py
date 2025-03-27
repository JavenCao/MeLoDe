modeljsonFile="/storage/yangjianLab/caoyujie/project/meth/5_Replication/8_CKCG_high_2Dmodel/model_training/H_model_training/H_training_26M.V1C1.npz.file.list.model.HD.CpGModel.json"
modelweightsFile=modeljsonFile + ".weights.h5"
from keras.models import model_from_json
json_file = open(modeljsonFile, 'r')
loaded_model_json = json_file.read()
json_file.close()
loaded_model = model_from_json(loaded_model_json)
loaded_model.load_weights(modelweightsFile)


import tensorflow as tf
saved_model="/storage/yangjianLab/caoyujie/project/meth/5_Replication/8_CKCG_high_2Dmodel/model_training/H_model_training/H_model"
# Convert the Keras model to the TensorFlow SavedModel format
tf.saved_model.save(loaded_model, saved_model)

