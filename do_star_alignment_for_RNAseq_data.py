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


def star(arg):
    cmd = " ".join(map(str,arg))
    sp.check_call(cmd,shell=True)
    return True


@click.command()
@click.option("--threads",
              default=10,
              type=int,
              show_default=True,
              help="number of threads used for trmmomatic")

@click.option("--trim_fastq_file_path",
              help="the data_file path")



def star_alignment(trim_fastq_file_path,threads):
    i = 0
    genomedir_path = str("path/zebrafish_genome_index/")
    data_file_list = []
    trim_fastq_file_path_str = str(trim_fastq_file_path)
    for f in glob.iglob(os.path.join(trim_fastq_file_path_str,"*fq")):
        data_file_list.append(str(os.path.basename(f).split("-")[0]))
    data_file_num = len(data_file_list)
    data_file_list_nore = list(set(data_file_list))
    print(data_file_list_nore)
    for file in data_file_list_nore:
        for data_file in glob.iglob(os.path.join(trim_fastq_file_path_str,"*fq")):
            if str(file)==str(os.path.basename(data_file).split("-")[0]):
                if str(os.path.basename(data_file).split("_")[-1]) == str("1.clean.fq"):
                    data_file_r1 = data_file
                    i +=1
                elif str(os.path.basename(data_file).split("_")[-1]) == str("2.clean.fq"):
                    data_file_r2 = data_file
                    i +=1
                if i <= data_file_num and (i%2==0):
                    outfile_name = str("_".join(os.path.basename(data_file).split("_")[0:-2]))
                    arg = ["STAR","--runThreadN",threads,"--genomeDir",genomedir_path,"--readFilesIn",data_file_r1,data_file_r2,"--outFileNamePrefix",outfile_name]
                    star(arg)
                    outfile_path = os.path.abspath(os.path.join(trim_fastq_file_path, ''+outfile_name+''))
                    os.mkdir(outfile_path)
                    outfile_name_remove = str(''+outfile_name+'')+str("*")
                    for f in glob.iglob(os.path.join(trim_fastq_file_path, ''+outfile_name_remove+'')):
                        shutil.move(f,outfile_path)
                elif i <= data_file_num and (not(i%2==0)):
                    pass
if __name__ == '__main__':
    star_alignment()
                
        





