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

def get_bam(arg1):
    cmd = " ".join(map(str,arg1))
    sp.check_call(cmd,shell=True)
    bam_file_exist = True
    return bam_file_exist

def sort_bam(arg2):
    cmd = " ".join(map(str,arg2))
    sp.check_call(cmd,shell=True)
    sort_file_exist = True
    return sort_file_exist
def get_sort_bai(arg4):
    cmd = " ".join(map(str,arg4))
    sp.check_call(cmd,shell=True)
    return True

def get_counts(arg3):
    cmd = " ".join(map(str,arg3))
    sp.check_call(cmd,shell=True)
    return True

@click.command()

@click.option("--file_path",
              help="the star_out_file path")

def get_reads_counts(file_path):
    file_path_list = []
    file_list = os.listdir(file_path)
    gff_file_path = str("path/Danio_rerio.GRCz11.92.chr.gtf")
    file_path_str = str(file_path)
    for f in file_list:
        path = os.path.join(file_path_str,f)
        file_path_list.append(path)
    print(file_path_list)
    
    for each in file_path_list:
        for each_file in glob.iglob(os.path.join(each,"*.sam")):
            bam_file_name = str(".".join(os.path.basename(each_file).split(".")[0:-1])) + str(".bam")
            bam_file_path = os.path.join(each,bam_file_name)
            arg1 = ["samtools","view","-bS",each_file,">",bam_file_path]
            bam_file_exist = get_bam(arg1)            
            if bam_file_exist:
                bam_file_path = os.path.abspath(os.path.join(''+each+'',''+bam_file_name+''))
                sort_file = str(".".join(os.path.basename(bam_file_path).split(".")[0:-1])) + str(".sort.bam")
                sort_file_name = str(".".join(os.path.basename(bam_file_path).split(".")[0:-1])) + str(".sort") + str(".bam")
                sort_file_path = os.path.join(each,sort_file)
                sort_file_name_path = os.path.join(each,sort_file_name)
                arg2 = [" samtools","sort","-@","20",bam_file_path,"-o",sort_file_path]
                sort_file_exist = sort_bam(arg2)
'''                if sort_file_exist:
                    count_file = str(".".join(os.path.basename(bam_file_path).split(".")[0:-1])) + str(".count")
                    count_file_path = os.path.join(each,count_file)
                    arg3 = ["htseq-count","-f","bam",sort_file_name_path,gff_file_path,">",count_file_path]
                    get_counts(arg3)
'''

    
if __name__ == '__main__':
    get_reads_counts()







    
    
