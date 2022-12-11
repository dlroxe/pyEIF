import community 
import fastcluster
import graph_tool
import igraph
import lifelines
import markov_clustering 
import matplotlib
from matplotlib_venn import venn3, venn3_circles
import mygene
import networkx
import numpy
import os
import pandas
import powerlaw  # Power laws are probability distributions with the form:p(x)∝x−α
import pygraphviz 
import pylab
import random
import seaborn 
import scipy
import scipy.spatial.distance as ssd 
import sklearn
import sklearn.cluster 
import sklearn.decomposition 
import subprocess


# critical parameter setting!
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams['figure.figsize'] = 15, 15
pandas.options.mode.chained_assignment = None  # default='warn'

# set current direcotry
data_file_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_data"
output_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_output"
data_output_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_output/ProcessedData"
