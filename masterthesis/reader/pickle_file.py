import os
import cPickle as pickle

__author__ = 'delur'


def read_picklefile(filename, paths):
    path = os.path.join(paths["working_dir"] , "pickles")
    infile = open(os.path.join(path, filename + ".pkl"), 'r')
    data = pickle.load(infile)
    infile.close()
    return data

