modeljsonFile="/storage/yangjianLab/caoyujie/project/meth/5_Replication/3_Training_ground/zmw_model/factors/6_strand_resolution/2_models/HG005_strand.json"
modelweightsFile=modeljsonFile + ".weights.h5"
from keras.models import model_from_json
json_file = open(modeljsonFile, 'r')
loaded_model_json = json_file.read()
json_file.close()
loaded_model = model_from_json(loaded_model_json)
loaded_model.load_weights(modelweightsFile)


import tensorflow as tf
saved_model="/storage/yangjianLab/caoyujie/project/meth/5_Replication/3_Training_ground/zmw_model/factors/6_strand_resolution/2_models/hifi_ds"
# Convert the Keras model to the TensorFlow SavedModel format
tf.saved_model.save(loaded_model, saved_model)

