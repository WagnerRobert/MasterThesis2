import os

__author__ = 'delur'
import cPickle as pickle

def write_picklefile(pickle_object, filename, paths):
    path = os.path.join(paths["working_dir"] , "pickles")
    if not os.path.exists(path):
        os.makedirs(path)


    output = open(os.path.join(path, filename + ".pkl"), 'wb')
    pickle.dump(pickle_object, output)
    output.close()