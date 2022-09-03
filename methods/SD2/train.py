import os
import sys
import time
import numpy as np
import pickle as pkl
import tensorflow as tf
from utils import *
from models import DSTG
import pandas as pd
import warnings
import argparse
parser = argparse.ArgumentParser()
warnings.filterwarnings("ignore")

# Set random seed
seed = 123
np.random.seed(seed)
tf.compat.v1.set_random_seed(seed)
tf.set_random_seed(seed)
#parser.add_argument("-m", help="this is the mode of position", type=str, required=True)
#args = parser.parse_args()
mode = 'pseudo'
modes = ['anterior1','anterior2','posterior1','posterior2','kidney',
         'pseudo','PDAC_A_GSM3036911','PDAC_A_GSM4100721','PDAC_A_GSM4100722','PDAC_B_GSM3405534','PDAC_B_GSM4100723','PDAC_B_GSM4100724']
if mode in modes:
    print('mode is right!')
# Settings
flags = tf.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_string('dataset', 'Datadir', 'Input data')
flags.DEFINE_string('result', 'DSTG_Result', 'Output result')
flags.DEFINE_string('model', 'DSTG', 'Model string.')
flags.DEFINE_float('learning_rate', 0.1, 'Initial learning rate.')
flags.DEFINE_integer('epochs', 100, 'Number of epochs to train.')
flags.DEFINE_integer('hidden1', 10, 'Number of units in hidden layer 1.')
flags.DEFINE_float('dropout', 0, 'Dropout rate (1 - keep probability).')
flags.DEFINE_float('weight_decay', 0,
                   'Weight for L2 loss on embedding matrix.')
flags.DEFINE_integer('early_stopping', 50,
                     'Tolerance for early stopping (# of epochs).')
# Load data


#real_ST indicates whether this experiment uses the real ST data or not
#mode indicates the position experiment uses

adj, features, labels_binary_train, labels_binary_val, \
labels_binary_test, train_mask, pred_mask, val_mask, test_mask, new_label, true_label,celltypes = load_data(
    FLAGS.dataset,mode)

support = [preprocess_adj(adj)]
num_supports = 1
model_func = DSTG

# Define placeholders
placeholders = {
    'support':
    [tf.sparse_placeholder(tf.float32) for _ in range(num_supports)],
    'features':
    tf.sparse_placeholder(tf.float32,
                          shape=tf.constant(features[2], dtype=tf.int64)),
    'labels':
    tf.placeholder(tf.float32, shape=(None, labels_binary_train.shape[1])),
    'labels_mask':
    tf.placeholder(tf.int32),
    'dropout':
    tf.placeholder_with_default(0., shape=()),
    'num_features_nonzero':
    tf.placeholder(tf.int32)  # helper variable for sparse dropout
}

# Create model
model = model_func(placeholders, input_dim=features[2][1], logging=True)


# Define model evaluation function
def evaluate(features, support, labels, mask, placeholders):
    t_test = time.time()
    feed_dict_val = construct_feed_dict(features, support, labels, mask,
                                        placeholders)
    outs_val = sess.run([model.loss, model.accuracy], feed_dict=feed_dict_val)
    return outs_val[0], outs_val[1], (time.time() - t_test)


# Initialize session
sess = tf.Session()
# Init variables
sess.run(tf.global_variables_initializer())

train_accuracy = []
train_loss = []
val_accuracy = []
val_loss = []
test_accuracy = []
test_loss = []

# Train model
for epoch in range(FLAGS.epochs):
    t = time.time()
    # Construct feed dictionary
    feed_dict = construct_feed_dict(features, support, labels_binary_train,
                                    train_mask, placeholders)
    feed_dict.update({placeholders['dropout']: FLAGS.dropout})
    # Training step
    outs = sess.run([model.opt_op, model.loss, model.accuracy],
                    feed_dict=feed_dict)
    train_accuracy.append(outs[2])
    train_loss.append(outs[1])
    # Validation
    cost, acc, duration = evaluate(features, support, labels_binary_val,
                                   val_mask, placeholders)
    val_loss.append(cost)
    val_accuracy.append(acc)
    # test_cost, test_acc, test_duration = evaluate(features, support,
    #                                               labels_binary_test,
    #                                               test_mask, placeholders)
    # test_accuracy.append(test_acc)
    # test_loss.append(test_cost)
    print("Epoch:", '%04d' % (epoch + 1), "train_loss=",
          "{:.5f}".format(outs[1]), "train_acc=", "{:.5f}".format(outs[2]),
          "val_loss=", "{:.5f}".format(cost), "val_acc=", "{:.5f}".format(acc),
          "time=", "{:.5f}".format(time.time() - t))
    if epoch > FLAGS.early_stopping and val_loss[-1] > np.mean(
            val_loss[-(FLAGS.early_stopping + 1):-1]):
        print("Early stopping...")
        break

print("Finished Training....")

#' --------------- --------------- ---------------
#'  all outputs : prediction and activation
#' --------------- --------------- ---------------

all_mask = np.array([True] * len(train_mask))
labels_binary_all = new_label

feed_dict_all = construct_feed_dict(features, support, labels_binary_all,
                                    all_mask, placeholders)
feed_dict_all.update({placeholders['dropout']: FLAGS.dropout})

activation_output = sess.run(model.activations, feed_dict=feed_dict_all)[1]
predict_output = sess.run(model.outputs, feed_dict=feed_dict_all)

#' ------- accuracy on prediction masks ---------
ab = sess.run(tf.nn.softmax(predict_output))

true_label1 = np.array(true_label)

result_file1 = '{}/predict_output.csv'.format(FLAGS.result)
result_file2 = '{}/true_output.csv'.format(FLAGS.result)
final_pred_output = pd.DataFrame(ab[pred_mask],columns=celltypes)
final_pred_output.to_csv(result_file1,header=True)
final_true_output = pd.DataFrame(true_label1[pred_mask],columns=celltypes)
final_true_output.to_csv(result_file2,header=True)