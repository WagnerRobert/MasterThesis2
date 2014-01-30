__author__ = 'delur'

import os
def setUp(constants):
    if not os.path.exists(os.path.join(constants["working_dir"], "kmers")):
        os.makedirs(os.path.join(constants["working_dir"], "kmers"))
    if not os.path.exists(constants["pdf"]):
        os.makedirs(constants["pdf"])