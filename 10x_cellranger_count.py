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
    cmd = "".join(map(str,arg))
    sp.check_call(cmd,shell=True)
    return True


@click.command()

###the file path include fastq files
@click.option("--fastqs",
              help="the data_file path")

@click.option("--transcriptome",
              help="the refference index file path")


def do_10x_count(fastqs,transcriptome):
    i = 0
    data_file_list = []
    fastqs_path_str = str(fastqs)
    for data_file_path in glob.iglob(os.path.join(fastqs_path_str,"*")):
        print(data_file_path)
        sample_id = str(os.path.basename(data_file_path))
        arg = ["path/cellranger-3.0.2/cellranger count"," --id=",sample_id," --transcriptome=",transcriptome," --fastqs=",fastqs," --sample=",sample_id]
        star(arg)
if __name__ == '__main__':
    do_10x_count()
                
        





