# a crude implementation of rnn based generator

# import stuff
import pandas as pd
import sys
import tensorflow as tf
import numpy as np
import os
import json
from find_files import find_files as fifi

# set df display to max
pd.set_option('max_column', None)


def vocab_dict_maker(infile, charsize):
    '''
    makes tags, vocabs, char2dict etc
    :param infile:
    :return:
    '''
    tag = '_'.join(infile.split('/')[-1].split('_')[:-1])
    tag = tag + '_%s' % charsize
    charsize = int(charsize)
    #seqs = open(infile).read()[:charsize]
    contents = open(infile).read().split(' ')
    seqs = ' '.join(contents[:charsize])
    vocab = sorted(set(seqs))
    char2idx = {u:i for i,u in enumerate(vocab)}
    idx2char = np.array(vocab)
    text_as_int = np.array([char2idx[c] for c in seqs])
    return tag, seqs, vocab, char2idx,idx2char, text_as_int



# split input target
def split_input_target(seq):
    '''
    split input output, return input-output pairs
    :param seq:
    :return:
    '''
    input_seq = seq[:-1]
    target_seq = seq[1:]
    return input_seq, target_seq

def create_model(vocab_size, embedding_dim, rnn_units, batch_size):
    '''
    create model according to params
    :param vocab_size:
    :param embedding_dim:
    :param rnn_units:
    :param batch_size:
    :return:
    '''
    model = tf.keras.Sequential([
        tf.keras.layers.Embedding(vocab_size, embedding_dim, batch_input_shape = [batch_size,None]),
        tf.keras.layers.LSTM(rnn_units,
                             return_sequences = True,
                             stateful =True,
                             recurrent_initializer = 'glorot_uniform'),
        tf.keras.layers.Dense(vocab_size)
    ])
    return model

def loss(labels, logits):
    '''
    loss function for the net output
    :param labels:
    :param logits:
    :return:
    '''
    loss = tf.keras.losses.sparse_categorical_crossentropy(labels, logits, from_logits = True)
    return loss




# run stuff
def train_runner(infile_path):
    '''
    runs trainer
    :return:
    '''
    number_of_seqs = int(infile_path.split('_')[-1].split('.')[0])
    print(number_of_seqs)
    # indices = [number_of_seqs/(10**2), number_of_seqs/(10**1), number_of_seqs]
    # indices = [int(item) for item in indices]
    indices = [number_of_seqs]
    for number_of_seq in indices[:]:
        tag, seqs, vocab, char2idx, idx2char, text_as_int  =  vocab_dict_maker(infile=infile_path, charsize=number_of_seq)
        print(tag)
        vocab_size = len(vocab)
        embedding_dim = 256
        rnn_units = 1024
        seq_length = 42 # window size (w)
        examples_per_epoch = len(seqs)//seq_length+1
        char_dataset = tf.data.Dataset.from_tensor_slices(text_as_int)
        sequences = char_dataset.batch(seq_length+1, drop_remainder=True)
        dataset = sequences.map(split_input_target)
        batch_size = 64
        buffer_size = 1000
        dataset = dataset.shuffle(buffer_size).batch(batch_size, drop_remainder=True)
        # get dataset size since there is no dim/shape attributes in tf.dataset
        dataset_size = 0
        for item in dataset:
            dataset_size += 1
        print(dataset_size)
        scaller = 1
        # split train, val, test
        train_size = int(0.7/scaller * dataset_size)
        val_size = int(0.15/scaller * dataset_size)
        test_size = int(0.15/scaller * dataset_size)
        print('Trains batches {}, val batches {}, test batches {}'.format(train_size, val_size, test_size))
        train_dataset = dataset.take(train_size)
        test_dataset = dataset.skip(train_size)
        val_dataset = test_dataset.skip(val_size)
        test_dataset = test_dataset.take(test_size)
        #####
        model = create_model(vocab_size, embedding_dim,rnn_units, batch_size)
        model.summary()
        # sanity checks
        for input_example_batch, output_example_batch in dataset.take(1):
            example_batch_predictions = model(input_example_batch)
            print(example_batch_predictions.shape)
            sampled_indices = tf.random.categorical(example_batch_predictions[0], num_samples=1)
            print(sampled_indices)
            sampled_indices = tf.squeeze(sampled_indices, axis=-1).numpy()
            print(''.join(idx2char[sampled_indices]))
            print(''.join(idx2char[input_example_batch[0].numpy()]))
        ####
        example_loss = loss(output_example_batch, example_batch_predictions)
        print(example_loss.numpy().mean())
        model.compile(optimizer='adam', loss=loss)
        checkpoint_dir = '%s_training_checkpoints' % tag
        checkpoint_prefix =  os.path.join(checkpoint_dir, 'ckpt_{epoch}')
        print(checkpoint_prefix)
        checkpoint_callback = tf.keras.callbacks.ModelCheckpoint(
            filepath = checkpoint_prefix,
            save_weights_only = True,
            verbose = 1,
            save_best_only=False
        )



        nb_epoch = 20

        # catch and write history

        history = model.fit(train_dataset, epochs = nb_epoch,
                            callbacks = [checkpoint_callback],
                            validation_data=val_dataset)
        history_outfile = 'outfiles/%s_history.txt' % tag
        history_contents =[]
        for key in history.history:
            for i,val in enumerate(history.history[key]):
                history_content = [key, i+1, val]
                history_contents.append(history_content)
        historydf = pd.DataFrame(history_contents, columns=['loss_cat', 'epoch', 'value'])
        historydf.to_csv(history_outfile, index=False)


        # write model params

        model_params = {
        'nb_epoch': nb_epoch,
        'batch_size':batch_size,
        'buffer_size':buffer_size,
        'vocab_size':vocab_size,
        'embedding_dim':embedding_dim,
        'rnn_units':rnn_units
        }

        model_params_outname = 'outfiles/%s_model_params.json' % tag
        char2idx_outname = 'outfiles/%s_char2idx.json' % tag
        json.dump(model_params, open(model_params_outname, 'w'))
        json.dump(char2idx, open(char2idx_outname, 'w'))


def train_runner_corpus(infile_path):
    '''
    runs trainer
    :return:
    '''
    number_of_seq = int(infile_path.split('_')[-2])
    tag, seqs, vocab, char2idx, idx2char, text_as_int  =  vocab_dict_maker(infile=infile_path, charsize=number_of_seq)
    print(tag)
    vocab_size = len(vocab)
    embedding_dim = 256
    rnn_units = 1024
    seq_length = 42 # window size (w)
    examples_per_epoch = len(seqs)//seq_length+1
    char_dataset = tf.data.Dataset.from_tensor_slices(text_as_int)
    sequences = char_dataset.batch(seq_length+1, drop_remainder=True)
    dataset = sequences.map(split_input_target)
    batch_size = 64
    buffer_size = 1000
    dataset = dataset.shuffle(buffer_size).batch(batch_size, drop_remainder=True)
    # get dataset size since there is no dim/shape attributes in tf.dataset
    dataset_size = 0
    for item in dataset:
        dataset_size += 1
    print(dataset_size)
    scaller = 1
    # split train, val, test
    train_size = int(0.7/scaller * dataset_size)
    val_size = int(0.15/scaller * dataset_size)
    test_size = int(0.15/scaller * dataset_size)
    print('Trains batches {}, val batches {}, test batches {}'.format(train_size, val_size, test_size))
    train_dataset = dataset.take(train_size)
    test_dataset = dataset.skip(train_size)
    val_dataset = test_dataset.skip(val_size)
    test_dataset = test_dataset.take(test_size)
    #####
    model = create_model(vocab_size, embedding_dim,rnn_units, batch_size)
    model.summary()
    # sanity checks
    for input_example_batch, output_example_batch in dataset.take(1):
        example_batch_predictions = model(input_example_batch)
        print(example_batch_predictions.shape)
        sampled_indices = tf.random.categorical(example_batch_predictions[0], num_samples=1)
        print(sampled_indices)
        sampled_indices = tf.squeeze(sampled_indices, axis=-1).numpy()
        print(''.join(idx2char[sampled_indices]))
        print(''.join(idx2char[input_example_batch[0].numpy()]))
    ####
    example_loss = loss(output_example_batch, example_batch_predictions)
    print(example_loss.numpy().mean())
    model.compile(optimizer='adam', loss=loss)
    checkpoint_dir = '%s_corpus_training_checkpoints' % tag
    checkpoint_prefix =  os.path.join(checkpoint_dir, 'ckpt_{epoch}')
    print(checkpoint_prefix)
    checkpoint_callback = tf.keras.callbacks.ModelCheckpoint(
        filepath = checkpoint_prefix,
        save_weights_only = True,
        verbose = 1,
        save_best_only=False
    )



    nb_epoch = 20

    # catch and write history

    history = model.fit(train_dataset, epochs = nb_epoch,
                        callbacks = [checkpoint_callback],
                        validation_data=val_dataset)
    history_outfile = 'outfiles/%s_history.txt' % tag
    history_contents =[]
    for key in history.history:
        for i,val in enumerate(history.history[key]):
            history_content = [key, i+1, val]
            history_contents.append(history_content)
    historydf = pd.DataFrame(history_contents, columns=['loss_cat', 'epoch', 'value'])
    historydf.to_csv(history_outfile, index=False)


    # write model params

    model_params = {
    'nb_epoch': nb_epoch,
    'batch_size':batch_size,
    'buffer_size':buffer_size,
    'vocab_size':vocab_size,
    'embedding_dim':embedding_dim,
    'rnn_units':rnn_units
    }

    model_params_outname = 'outfiles/%s_model_params.json' % tag
    char2idx_outname = 'outfiles/%s_char2idx.json' % tag
    json.dump(model_params, open(model_params_outname, 'w'))
    json.dump(char2idx, open(char2idx_outname, 'w'))



# run stuff
#train_runner('outfiles/1OB1_5pctsData_top_70334.txt')
#train_runner('outfiles/1OB1_5pctsData_bot_70334.txt')
# antigens set runs
#infiles = fifi('outfiles', 'corpus.txt')
infiles = fifi('datasets/eleven/nsamples2', 'corpus.txt')
infiles= [item for item in infiles if 'sample' in item]
infiles= [item for item in infiles if '2YPV' in item]
print(infiles)
print(len(infiles))
for infile in infiles[:]:
    train_runner_corpus(infile)


