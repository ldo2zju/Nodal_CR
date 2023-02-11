#!/usr/bin/env python
# coding=utf-8
import click
import fileinput
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import shutil
import six
import string
import subprocess as sp
import sys
import time
import warnings
from datetime import datetime
import glob
from matplotlib.cbook import MatplotlibDeprecationWarning
from xml.etree import cElementTree as ET
from multiprocessing import Process, Pool
import threading
import os.path

def get_featurecount(arg):
    cmd = " ".join(map(str,arg))
    sp.check_call(cmd,shell=True)
    get_file = True
    return get_file

def get_cut_featurecount(arg2):
    cmd = " ".join(map(str,arg2))
    sp.check_call(cmd,shell=True)
    return True


@click.command()

@click.option("--sorted_bam_file_path",
              help="the sorted_bam_file_path")

@click.option("--threads",
              default=10,
              type=int,
              show_default=True,
              help="number of threads used for featurecounts")
def get_featurecounts_file(sorted_bam_file_path,threads):
    gtf_file_path = str("path/Danio_rerio.GRCz11.92.chr.gtf")
    file_path_str = str(sorted_bam_file_path)
    file_list = os.listdir(sorted_bam_file_path)
    file_path_list = []
    for f in file_list:
        path = os.path.join(file_path_str,f)
        file_path_list.append(path)
    for each in file_path_list:
        for each_file in glob.iglob(os.path.join(each,"*sort.bam")):
            print(each_file)
            '''sort_bam_file_path = os.path.join(each,each_file)'''
            feature_count_file_name = str((os.path.basename(each)) + str("_featureCounts.txt"))
            feature_count_cut_file_name = str((os.path.basename(each)) + str("_cut_featureCounts.txt"))
            feature_count_file_path = os.path.join(each,feature_count_file_name)
            feature_count_cut_file_path = os.path.join(each,feature_count_cut_file_name)
            arg = ["featureCounts ","-T",threads,"-p","-t","exon","-g","gene_id ","-a",gtf_file_path,"-o",feature_count_file_path,each_file]
            get_file = get_featurecount(arg)
        for each_file in glob.iglob(os.path.join(each,"*summary")):
            print(each_file)
            if os.path.exists(each_file):
                arg2 = ["cut ","-f","1,7",feature_count_file_path,"|grep","-v","'^#'",">",feature_count_cut_file_path]
                get_cut_featurecount(arg2)      
            else:
                pass

if __name__ == '__main__':
    get_featurecounts_file()
                                          














        
