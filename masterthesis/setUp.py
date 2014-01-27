__author__ = 'delur'

import os
def setUp(path):
    if not os.path.exists(os.path.join(path, "kmers")):
        os.makedirs(os.path.join(path, "kmers"))