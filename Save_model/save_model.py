modelweightsFile="/storage/yangjianLab/caoyujie/project/meth/5_Replication/3_Training_ground/zmw_model/factors/4_ancesty/7_train_by_all/2_ZMW_model/ZMW_model_5_inds.json.weights.h5"
modeljsonFile="/storage/yangjianLab/caoyujie/project/meth/5_Replication/3_Training_ground/zmw_model/factors/4_ancesty/7_train_by_all/2_ZMW_model/ZMW_model_5_inds.json"
from keras.models import model_from_json
json_file = open(modeljsonFile, 'r')
loaded_model_json = json_file.read()
json_file.close()
loaded_model = model_from_json(loaded_model_json)
loaded_model.load_weights(modelweightsFile)


import tensorflow as tf
saved_model="/storage/yangjianLab/caoyujie/project/meth/5_Replication/3_Training_ground/zmw_model/factors/4_ancesty/7_train_by_all/2_ZMW_model/Saved_model/Saved_model"
# Convert the Keras model to the TensorFlow SavedModel format
tf.saved_model.save(loaded_model, saved_model)

