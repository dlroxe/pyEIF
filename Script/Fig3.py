#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 14:48:10 2020

@author: suwu
"""
import community 
import fastcluster
import graph_tool
import igraph
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


# critical parameter setting!
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams['figure.figsize'] = 15, 15
pandas.options.mode.chained_assignment = None  # default='warn'



# set current direcotry
data_file_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_data"
output_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_output/Fig3"
data_output_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_output/ProcessedData"

# retrieve published proteomics data ######
def ccle_pro():
    # download the proteomics data
    # https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz
    CCLE_PRO = pandas.read_csv(os.path.join(data_file_directory, 
                                            "protein_quant_current_normalized.csv"))
    # concatenate the Gene_Symbol (with duplicate names) and the Uniprot_Acc
    CCLE_PRO["Gene_Symbol"] = (CCLE_PRO["Gene_Symbol"].fillna("") + " " + CCLE_PRO["Uniprot_Acc"])
    CCLE_PRO = CCLE_PRO.set_index("Gene_Symbol")
    CCLE_PRO_subset = CCLE_PRO[CCLE_PRO.columns.drop(list(CCLE_PRO.filter(regex="Peptides")))]
    CCLE_PRO_subset = CCLE_PRO_subset.drop(columns=["Protein_Id", 
                                                    "Description", 
                                                    "Group_ID", 
                                                    "Uniprot", 
                                                    "Uniprot_Acc"]).T
    CCLE_PRO_subset.index.names = ["ccle_name"]
    return CCLE_PRO_subset

CCLE_PRO_subset = ccle_pro()
# CCLE_PRO_subset[CCLE_PRO_subset.duplicated()]

# find eif containing columns
EIF_cols = [col for col in CCLE_PRO_subset.columns if "EIF3" in col]
print(EIF_cols)

EIF2 = pandas.Series(["EIF4G1 Q04637-9", "EIF4A1 P60842", "EIF4E P06730-2", 
                      "EIF4EBP1 Q13541", 'EIF3E P60228', 'EIF3H O15372'])


# retrieve EIF4F proteomics data for correlation scatter plot ######
def reg_coef(x, y, label=None, color=None, **kwargs):
    ax = matplotlib.pyplot.gca()
    r, p = scipy.stats.pearsonr(x, y)
    ax.annotate(  #'data = (%.2f, %.2f)'%(r, p),
        "r = {:.3f}".format(r), xy=(0.5, 0.5), xycoords="axes fraction", ha="center"
    )
    ax.annotate(  #'data = (%.2f, %.2f)'%(r, p),
        "p = {:.3f}".format(p), xy=(0.5, 0.4), xycoords="axes fraction", ha="center"
    )
    ax.set_axis_off()


def eif_ccle_scatter(data, proteins):
    CCLE_EIF_PRO = data[proteins]
    CCLE_EIF_PRO.index.name = 'ccle'
    CCLE_EIF_PRO.reset_index(inplace=True)
    CCLE_EIF_PRO["ccle"] = CCLE_EIF_PRO["ccle"].str.split("_").str[1]
    CCLE_EIF_PRO.dropna(inplace=True)
    matplotlib.rcParams.update({'font.size': 12})
    g = seaborn.PairGrid(CCLE_EIF_PRO)
    # g.map(sbn.scatterplot)
    g.map_diag(seaborn.histplot)
    # g.map_lower(sbn.scatterplot)
    g.map_lower(seaborn.regplot)
    g.map_upper(reg_coef)
    g.tight_layout
    #g.add_legend()
    seaborn.set_style("ticks")
    matplotlib.pyplot.show()
    matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), 
    "CCLE_scatter.pdf"), dpi=300)
    # scatter plot of two protein expression across cell lines with color
    g = seaborn.PairGrid(CCLE_EIF_PRO, hue="ccle", diag_sharey=False, corner = True)
    g.map_lower(seaborn.scatterplot)
    g.map_diag(seaborn.histplot)
    g.map_upper(reg_coef)
    g.tight_layout
    g.add_legend()
    matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), 
    "CCLE_color_scatter.pdf"), dpi=300)
    matplotlib.pyplot.show()

eif_ccle_scatter(data = CCLE_PRO_subset, proteins = EIF2)


# identify CORs for eIF4F ######
def eif_corr_sum(df, gene_list):
    y = pandas.DataFrame()
    for gene in gene_list:
        print(gene)
        x = df.corrwith(df[gene])
        x.name = gene
        y = pandas.concat([y, x], axis=1)
    y.dropna(inplace=True)
    y.index.names = ["Gene_Symbol"]
    # pandas.DataFrame.ge: get Greater than or equal to of dataframe and other
    # pandas.DataFrame.le: get Less than or equal to of dataframe and other
    z = y[y.ge(0.5).any(axis=1)|y.le(-0.5).any(axis=1)]
    return (y, z)

EIF_COR_PRO, EIF_COR_PRO_sig = eif_corr_sum(df = CCLE_PRO_subset, gene_list = EIF2)
EIF_COR_PRO_sig.to_csv(os.path.join(os.path.expanduser(data_output_directory), "EIF_COR_PRO_sig.csv"))


# identify posCORs for eIF4F ######
# use pearsonr to identify r and p-value, select r > 0.5 and p < 0.05 #####
def eif_corr_sig_p(eif):
    y = pandas.DataFrame(columns=["r", "p"])
    r1 = scipy.stats.pearsonr(CCLE_PRO_subset[eif].dropna(), CCLE_PRO_subset[eif].dropna())
    r1 = pandas.Series(r1, index=["r", "p"])
    r1.name = eif
    y = y.append(r1)
    x = CCLE_PRO_subset.columns.to_list()
    x1 = x.remove(eif)
    for gene in x:
        print(gene)
        CCLE_PRO_two = CCLE_PRO_subset[[gene, eif]].dropna()
        if CCLE_PRO_two.empty == True:
            print(gene, "too many NAs")
        else:
            r = scipy.stats.pearsonr(CCLE_PRO_two[gene], CCLE_PRO_two[eif])
            r = pandas.Series(r, index=["r", "p"])
            r.name = gene
            y = y.append(r)
    COR = y.loc[(y["p"] <= 0.05)]
    posCOR = COR.loc[(COR["r"] >= 0.5)]
    negCOR = COR.loc[(COR["r"] <= -0.5)]
    return (posCOR, negCOR)

eIF4G1_pos, eIF4G1_neg = eif_corr_sig_p("EIF4G1 Q04637-9")
eIF4A1_pos, eIF4A1_neg = eif_corr_sig_p("EIF4A1 P60842")
eIF4E_pos, eIF4E_neg = eif_corr_sig_p("EIF4E P06730-2")
4EBP1_pos, 4EBP1_neg = eif_corr_sig_p("EIF4EBP1 Q13541")


def plot_EIF_Venn():
    fig = matplotlib.pyplot.figure(figsize=(12, 10))
    # eIF4G1_pos, eIF4G1_neg = eif_corr_sig_p('EIF4G1 Q04637-9')
    # eIF4A1_pos, eIF4A1_neg = eif_corr_sig_p('EIF4A1 P60842')
    # eIF4E_pos, eIF4E_neg = eif_corr_sig_p('EIF4E P06730-2')
    # EBP1_pos, EBP1_neg = eif_corr_sig_p('EIF4EBP1 Q13541')
    v = venn3(
        [set(eIF4G1_pos.index), set(eIF4A1_pos.index), set(eIF4E_pos.index)],
        ("posCOR_eIF4G1", "posCOR_eIF4A1", "posCOR_eIF4E"),
    )
    matplotlib.pyplot.show()
    matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), 
    "CCLE_pos_Venn.pdf"), dpi=300)
    fig = matplotlib.pyplot.figure(figsize=(12, 10))
    v = venn3(
        [set(eIF4G1_neg.index), set(eIF4A1_neg.index), set(eIF4E_neg.index)],
        ("negCOR_eIF4G1", "negCOR_eIF4A1", "negCOR_eIF4E"),
    )
    matplotlib.pyplot.show()
    matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), 
    "CCLE_neg_Venn.pdf"), dpi=300)

plot_EIF_Venn()

# Displaying dataframe as an heatmap
# with diverging colourmap as coolwarm
h = seaborn.clustermap(
    EIF_COR_PRO_sig,
    figsize=(10, 10),
    #method="centroid",
    #metric="euclidean",
    tree_kws=dict(linewidths=0.5, colors=(0.2, 0.2, 0.4)),
    cmap="coolwarm",
    #cbar_pos=(.2, .2, .03, .4),
    yticklabels=False
)
matplotlib.pyplot.show()
matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), 
"COR_heatmap.pdf"), 
dpi=300)


##############################
###### network analysis ######
##############################
# prepare the reference network data
def protein_interaction_reference():
    # HuRI = pd.read_csv("~/Downloads/HI-union.tsv", sep='\t', header = None)
    # download the reference file
    # https://stringdb-static.org/download/protein.physical.links.detailed.v11.0/9606.protein.physical.links.detailed.v11.5.txt.gz
    HuRI = pandas.read_csv(os.path.join(data_file_directory,"9606.protein.physical.links.detailed.v11.5.txt"), sep=" ")
    HuRI.head
    HuRI['protein1'] = HuRI['protein1'].str.split('\.').str[-1].str.strip()
    HuRI['protein2'] = HuRI['protein2'].str.split('\.').str[-1].str.strip()
    codes1, uniques1 = pandas.factorize(HuRI['protein1'])
    codes2, uniques2 = pandas.factorize(HuRI['protein2'])
    HuRI = HuRI[['protein1', 'protein2', 'experimental', 'database']]
    # Mapping ensembl gene ids to gene symbols¶
    mg = mygene.MyGeneInfo()
    node1 = mg.querymany(uniques1,
                         scopes='ensembl.protein',
                         fields='symbol',
                         species='human',
                         as_dataframe=True)
    node2 = mg.querymany(uniques2,
                         scopes='ensembl.protein',
                         fields='symbol',
                         species='human',
                         as_dataframe=True)
    dict1 = pandas.Series(node1.symbol.values,
                      index=node1.index).to_dict()
    dict2 = pandas.Series(node2.symbol.values,
                      index=node2.index).to_dict()
    HuRI['protein1'] = HuRI['protein1'].map(dict1)
    HuRI['protein2'] = HuRI['protein2'].map(dict2)
    HuRI = HuRI.dropna()  # some na in 'protein2' column, remove them
    # chose experimental and database values over 1
    HuRI = HuRI[(HuRI['experimental'] != 0) & (HuRI['database'] != 0)]
    CCLE_pro = pandas.Series(CCLE_PRO_subset.columns).apply(
        lambda x: x.split(' ')[0])
    HuRI_CCLE = HuRI[HuRI['protein1'].isin(
        CCLE_pro) & HuRI['protein2'].isin(CCLE_pro)]
    HuRI_CCLE = HuRI_CCLE.dropna()
    #HuRI_CCLE[HuRI_CCLE['protein1'] != HuRI_CCLE['protein2']]
    return (HuRI, HuRI_CCLE)

HuRI, HuRI_CCLE = protein_interaction_reference()

# retrieve the clustering genes
EIF_COR_sig = pandas.DataFrame(EIF_COR_PRO_sig.index.str.split(" ").str[0]).dropna()
EIF_COR_sig.columns = ['gene']

cluster1 = pandas.read_csv(os.path.join(os.path.expanduser(data_output_directory), "cluster1.csv")) 
cluster2 = pandas.read_csv(os.path.join(os.path.expanduser(data_output_directory), "cluster2.csv"))
cluster3 = pandas.read_csv(os.path.join(os.path.expanduser(data_output_directory), "cluster3.csv")) 
cluster4 = pandas.read_csv(os.path.join(os.path.expanduser(data_output_directory), "cluster4.csv")) 

cluster_Net = HuRI_CCLE[HuRI_CCLE['protein1'].isin(EIF_COR_sig['gene']) & HuRI_CCLE['protein2'].isin(EIF_COR_sig['gene'])]
cluster1_Net = HuRI_CCLE[HuRI_CCLE['protein1'].isin(cluster1['gene']) & HuRI_CCLE['protein2'].isin(cluster1['gene'])]
cluster2_Net = HuRI_CCLE[HuRI_CCLE['protein1'].isin(cluster2['gene']) & HuRI_CCLE['protein2'].isin(cluster2['gene'])]
cluster3_Net = HuRI_CCLE[HuRI_CCLE['protein1'].isin(cluster3['gene']) & HuRI_CCLE['protein2'].isin(cluster3['gene'])]
cluster4_Net = HuRI_CCLE[HuRI_CCLE['protein1'].isin(cluster4['gene']) & HuRI_CCLE['protein2'].isin(cluster4['gene'])]

codes, uniques = pandas.factorize(cluster3_Net['protein1'])
cluster3_Net['protein2'].value_counts()


## plot reference interaction map
    G = networkx.from_pandas_edgelist(HuRI_CCLE,
                                'protein1',
                                'protein2',
                                # edge_attr = ['experimental','database']
                                )
    C = networkx.from_pandas_edgelist(cluster_Net,
                                'protein1',
                                'protein2',
                                # edge_attr = ['experimental','database']
                                )
    C1 = networkx.from_pandas_edgelist(cluster1_Net,
                                'protein1',
                                'protein2',
                                # edge_attr = ['experimental','database']
                                )     
    C2 = networkx.from_pandas_edgelist(cluster2_Net,
                                'protein1',
                                'protein2',
                                # edge_attr = ['experimental','database']
                                )                                       
    C3 = networkx.from_pandas_edgelist(cluster3_Net,
                                'protein1',
                                'protein2',
                                # edge_attr = ['experimental','database']
                                )     
    C4 = networkx.from_pandas_edgelist(cluster4_Net,
                                'protein1',
                                'protein2',
                                # edge_attr = ['experimental','database']
                                )                                     
###test
    C = networkx.from_pandas_edgelist(cluster_Net,source='protein1',
                                   target='protein2',#edge_attr=True,
                                   create_using=networkx.DiGraph())
    C1 = networkx.from_pandas_edgelist(cluster1_Net,source='protein1',
                                   target='protein2',#edge_attr=True,
                                   create_using=networkx.DiGraph())   
    C2 = networkx.from_pandas_edgelist(cluster2_Net,source='protein1',
                                   target='protein2',#edge_attr=True,
                                   create_using=networkx.DiGraph())   
    C3 = networkx.from_pandas_edgelist(cluster3_Net,source='protein1',
                                   target='protein2',#edge_attr=True,
                                   create_using=networkx.DiGraph())    
    C4 = networkx.from_pandas_edgelist(cluster4_Net,source='protein1',
                                   target='protein2',#edge_attr=True,
                                   create_using=networkx.DiGraph()) 
###test                                   

    pos = networkx.kamada_kawai_layout(G)
    posC = networkx.kamada_kawai_layout(C)

    matplotlib.pyplot.figure(figsize=(10, 10))
    matplotlib.pyplot.axis('off')
    networkx.draw_networkx(C,
                           posC,
                           width= 0.5, 
                           node_color ='#cccccc', 
                           edge_color="#cccccc",
                           node_size = 5, 
                           with_labels = False,
                           connectionstyle="arc3,rad=0.5")   
    networkx.draw_networkx(C1,
                           posC, 
                           width= 0.5, 
                           node_color ='green', 
                           edge_color="lightgreen",
                           node_size = 5, 
                           label='cluster 1',
                           with_labels = False,
                           connectionstyle="arc3,rad=0.5") 
    networkx.draw_networkx(C2,
                           posC, 
                           width= 0.5, 
                           node_color ='orange', 
                           edge_color="gold",
                           node_size = 5, 
                           label='cluster 2',
                           with_labels = False,
                           connectionstyle="arc3,rad=0.5")            
    networkx.draw_networkx(C3,
                           posC, 
                           width= 0.5, 
                           node_color ='blue', 
                           edge_color="skyblue",
                           node_size = 5, 
                           label='cluster 3',
                           with_labels = False,
                           connectionstyle="arc3,rad=0.5")            
    networkx.draw_networkx(C4,
                           posC, 
                           width= 0.5, 
                           node_color ='red', 
                           edge_color="pink",
                           node_size = 5, 
                           label='cluster 4',
                           with_labels = False,
                           connectionstyle="arc3,rad=0.5")                                  
    matplotlib.pyplot.legend()
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), 'CORs_network_kk.pdf'),
                    dpi=300,
                    edgecolor='w',
                    #transparent=True,
                    bbox_inches='tight')   
    matplotlib.pyplot.show()
    
def plot_cluster_net(cluster, pos, label, edge_color, node_color):
  # Degree Centrality: the number of edges a node has
  # Important nodes have many connections
  degree = networkx.degree_centrality(cluster)
  degree_names, degree_ranks = zip(*degree.items())
  degree_df = pandas.DataFrame(data={'protein': degree_names, 'rank': degree_ranks})
  #degree_df['rank'].values*100
  # PageRank assigns a score of importance to each node. 
  # Important nodes are those with many in-links from important pages.
  pagerank = networkx.pagerank(cluster)
  # save the names and their respective scores separately
  pagerank_names, pagerank_ranks = zip(*pagerank.items())
  pagerank_df = pandas.DataFrame(data={'protein': pagerank_names, 'rank': pagerank_ranks})
  #pagerank_df['rank'].values*10000
    
  # Closeness centrality identifies nodes that are, on average, closest to other nodes.
  # important nodes are close to other nodes
  closeness = networkx.closeness_centrality(cluster)
  closeness_names, closeness_ranks = zip(*closeness.items())
  closeness_df = pandas.DataFrame(data={'protein': closeness_names, 'rank': closeness_ranks})
  #closeness_df['rank'].values*100
    
  # Betweenness centrality identifies bridges and brokers:
  # important nodes are those that connect other nodes.
  betweenness = networkx.betweenness_centrality(cluster, normalized=False)
  betweenness_names, betweenness_ranks = zip(*betweenness.items())
  betweenness_df = pandas.DataFrame(data={'protein': betweenness_names, 'rank': betweenness_ranks})
  #betweenness_df['rank'].values/100

  # eigenvector centrality measures the transitive influence of nodes
  # Eigenvector centrality identifies nodes that are connected to other well-connected nodes
  eigenvector = networkx.eigenvector_centrality(cluster)
  eigenvector_names, eigenvector_ranks = zip(*eigenvector.items())
  eigenvector_df = pandas.DataFrame(data={'protein': eigenvector_names, 'rank': eigenvector_ranks})
  #eigenvector_df['rank'].values*100

  pos_higher = {}
  for k, v in pos.items():
    if(v[1]>0):
        pos_higher[k] = (v[0]-0.015, v[1]+0.015)
    else:
        pos_higher[k] = (v[0]-0.015, v[1]-0.015)

  matplotlib.pyplot.figure(figsize=(10, 10))
  matplotlib.pyplot.axis('off')
  networkx.draw_networkx_edges(cluster,
                           pos,
                           width= 0.5,
                           edge_color= edge_color)    
  networkx.draw_networkx_nodes(cluster, 
                           pos,
                           node_color =node_color, 
                           label=label,
                           node_size = degree_df['rank'].values*100) 
  networkx.draw_networkx_labels(cluster, pos_higher, font_size=3.5)
  matplotlib.pyplot.tight_layout()
  matplotlib.pyplot.savefig(
      os.path.join(os.path.expanduser(output_directory), 
                   (label + "_network_kk.pdf")),
      dpi=300,
      edgecolor='w',
      bbox_inches='tight')
  matplotlib.pyplot.show()
  
  #community detection
  partition = community.best_partition(cluster)
  # color the nodes according to their partition
  matplotlib.pyplot.figure(figsize=(10, 10))
  matplotlib.pyplot.axis('off')
  cmap = matplotlib.cm.get_cmap('Set2', max(partition.values()) + 1)
  networkx.draw_networkx_nodes(cluster, 
                              pos, nodelist= partition.keys(), 
                              node_size=degree_df['rank'].values*100, 
                              cmap=cmap, node_color = list(partition.values()))
  networkx.draw_networkx_edges(cluster, 
                              pos, 
                              alpha=0.5, 
                              width= 0.5, 
                              edge_color="#cccccc")
  networkx.draw_networkx_labels(cluster, pos_higher, font_size=3.5)
  matplotlib.pyplot.tight_layout()
  matplotlib.pyplot.savefig(os.path.join(
      os.path.expanduser(output_directory), 
      (label + "_network_community.pdf")),
      dpi=300,
      edgecolor='w',
      bbox_inches='tight')    
  matplotlib.pyplot.show()
  
plot_cluster_net(G, pos, "all CORS", "grey", "black")    
plot_cluster_net(C1, posC, "cluster 1", "lightgreen", "green")    
plot_cluster_net(C2, posC, "cluster 2", "gold", "orange")  
plot_cluster_net(C3, posC, "cluster 3", "skyblue", "blue")  
plot_cluster_net(C4, posC, "cluster 4", "pink", "red")  

def plot_degree_distribution(cluster,pos,label,edge_color,node_color):
  # PageRank assigns a score of importance to each node. 
  # Important nodes are those with many in-links from important pages.
  pagerank = networkx.pagerank(cluster)
  # save the names and their respective scores separately
  pagerank_names, pagerank_ranks = zip(*pagerank.items())
  pagerank_df = pandas.DataFrame(data={'protein': pagerank_names, 'rank': pagerank_ranks})
  #pagerank_df['rank'].values*10000
  # used for degree distribution and powerlaw test
  degree_freq = networkx.degree_histogram(cluster)
  degrees = range(len(degree_freq))
  fig, ax = matplotlib.pyplot.subplots(figsize=(15, 15))
  matplotlib.pyplot.loglog(degrees, degree_freq, 'go-')
  matplotlib.pyplot.title('Degree distribution plots')
  matplotlib.pyplot.xlabel('log of the degree')
  matplotlib.pyplot.ylabel('Log of number of nodes of that degree')
    # draw graph in inset
  matplotlib.pyplot.axes([0.25,0.25,0.65,0.65])
  matplotlib.pyplot.axis('off')
  networkx.draw_networkx_edges(cluster,
                           pos,  
                           width= 0.5, 
                           edge_color=edge_color, 
                           connectionstyle="arc3,rad=0.5")                           
  networkx.draw_networkx_nodes(cluster, 
  pos, 
  node_color=node_color, 
  label=label, 
  node_size= pagerank_df['rank'].values*10000)    
  matplotlib.pyplot.show()

plot_degree_distribution(cluster = G, pos= pos, label = "all proteins", edge_color = "grey", node_color = "black")    
plot_degree_distribution(cluster = C1, pos= posC, label = "cluster 1", edge_color = "lightgreen", node_color = "green")    
plot_degree_distribution(cluster = C2, pos= posC, label = "cluster 2", edge_color = "gold", node_color = "orange")  
plot_degree_distribution(cluster = C3, pos= posC, label = "cluster 3", edge_color = "skyblue", node_color = "blue")  
plot_degree_distribution(cluster = C4, pos= posC, label = "cluster 4", edge_color = "pink", node_color = "red")  



##
def plot_degree_histogram(cluster, pos, label, edge_color, node_color):
  degree_sequence = sorted([d for n, d in cluster.degree()], reverse=True)  # degree sequence
  # Degree Centrality: the number of edges a node has
  # Important nodes have many connections
  degree = networkx.degree_centrality(cluster)
  degree_names, degree_ranks = zip(*degree.items())
  degree_df = pandas.DataFrame(data={'protein': degree_names, 'rank': degree_ranks})
  fig, ax = matplotlib.pyplot.subplots(figsize=(10, 10))
  matplotlib.pyplot.bar(*numpy.unique(degree_sequence, return_counts=True))
  matplotlib.pyplot.title("Degree histogram " + label)
  matplotlib.pyplot.ylabel("Node count")
  matplotlib.pyplot.xlabel("Degree")
  matplotlib.pyplot.xlim([0, 80])
  matplotlib.pyplot.ylim([0, 60])
  matplotlib.pyplot.axes([0.25,0.25,0.65,0.65])
  matplotlib.pyplot.axis('off')
  networkx.draw_networkx_edges(cluster,
                             pos=pos,
                             width= 0.5,
                             edge_color= edge_color)    
  networkx.draw_networkx_nodes(cluster, pos=pos,
                             node_color =node_color, 
                             label=label,
                             node_size = degree_df['rank'].values*100) 
  matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), (label + "_Degree histogram.pdf")),
                      dpi=300,
                      edgecolor='w',
                      #transparent=True,
                      bbox_inches='tight') 
  matplotlib.pyplot.show()
  
plot_degree_histogram(cluster = G, pos= pos, label = "all proteins", edge_color = "grey", node_color = "black")    
plot_degree_histogram(cluster = C1, pos= posC, label = "cluster 1", edge_color = "lightgreen", node_color = "green")    
plot_degree_histogram(cluster = C2, pos= posC, label = "cluster 2", edge_color = "gold", node_color = "orange")  
plot_degree_histogram(cluster = C3, pos= posC, label = "cluster 3", edge_color = "skyblue", node_color = "blue")  
plot_degree_histogram(cluster = C4, pos= posC, label = "cluster 4", edge_color = "pink", node_color = "red")  


matplotlib.pyplot.clf()
matplotlib.pyplot.cla()
matplotlib.pyplot.close()




    # calculates the betweenness centrality and shows the 10 individuals with the highest value
    # Betweenness centrality identifies bridges and brokers:
    # edges and nodes that connect otherwise poorly connected parts of a network
    betweenness = networkx.betweenness_centrality(C3, normalized=False)
    top_betweenness = sorted(betweenness.items(),
                             key=lambda x: x[1],
                             reverse=True)[0:20]
    protein, betweenness = zip(*top_betweenness)
    x_pos = numpy.arange(len(protein))
    matplotlib.pyplot.figure(figsize=(10, 4))
    matplotlib.pyplot.bar(x_pos, betweenness, align='center')
    matplotlib.pyplot.xticks(x_pos, protein, rotation=45, ha="right", fontsize=8)
    matplotlib.pyplot.ylabel('top betweenness')
    matplotlib.pyplot.show()


# Closeness centrality identifies nodes that are, on average, closest to other nodes.
    closeness = networkx.closeness_centrality(C3)
    top_closeness = sorted(closeness.items(),
                           key=lambda x: x[1],
                           reverse=True)[0:20]
    protein, closeness = zip(*top_closeness)
    x_pos = numpy.arange(len(protein))
    matplotlib.pyplot.figure(figsize=(10, 4))
    matplotlib.pyplot.bar(x_pos, closeness, align='center')
    matplotlib.pyplot.xticks(x_pos, protein, rotation=45, ha="right", fontsize=8)
    matplotlib.pyplot.ylabel('top closeness')
    matplotlib.pyplot.show()

# highly-connected hubs by a measure called eigenvector centrality
# Eigenvector centrality identifies nodes that are connected to other well-connected nodes
    eigenvector = networkx.eigenvector_centrality(C)
    # sort eigenvector centrality from highest to lowest
    top_eigenvector = sorted(eigenvector.items(),
                             key=lambda x: x[1],
                             reverse=True)[0:20]
    # save the names and their respective scores separately
    # reverse the tuples to go from most frequent to least frequent
    protein, eigenvector = zip(*top_eigenvector)
    x_pos = numpy.arange(len(protein))
    matplotlib.pyplot.figure(figsize=(10, 4))
    matplotlib.pyplot.bar(x_pos, eigenvector, align='center')
    matplotlib.pyplot.xticks(x_pos, protein, rotation=45, ha="right", fontsize=8)
    matplotlib.pyplot.ylabel('top eigenvector')
    matplotlib.pyplot.show()


    def centrality_histogram(x, title):
        matplotlib.pyplot.figure(figsize=(10, 10))
        #matplotlib.pyplot.axis('off')
        matplotlib.pyplot.hist(x, 50, density=True, edgecolor='black')
        matplotlib.pyplot.title(title)
        matplotlib.pyplot.xlabel("Centrality")
        matplotlib.pyplot.ylabel("Density")
        matplotlib.pyplot.show()

    # how much centrality is concentrated in one or a few nodes.
    centrality_histogram(networkx.eigenvector_centrality(C3).values(),
        title="eigenvector centrality histogram (cluster 3)")
    centrality_histogram(networkx.betweenness_centrality(C3).values(),
        title="betweenness centrality histogram (cluster 3)")
    centrality_histogram(networkx.closeness_centrality(C3).values(),
        title="closeness centrality histogram (cluster 3)")

    #  the number of triangles between a node and its neighbors
    triangles = networkx.triangles(C3)
    sorted(triangles.items(), key=lambda x: x[1], reverse=True)[0:10]

    # draws a histogram of all shortest path lengths within a network
    def path_length_histogram(x, title=None):
        # Find path lengths
        length_source_target = dict(networkx.shortest_path_length(x))
        # Convert dict of dicts to flat list
        all_shortest = sum([
            list(length_target.values())
            for length_target
            in length_source_target.values()],
        [])
        # Calculate integer bins
        high = max(all_shortest)
        bins = [-0.5 + i for i in range(high + 2)]
        # Plot histogram
        matplotlib.pyplot.figure(figsize=(10, 10))
        matplotlib.pyplot.hist(all_shortest, bins=bins, rwidth=0.8)
        matplotlib.pyplot.title(title)
        matplotlib.pyplot.xlabel("Distance")
        matplotlib.pyplot.ylabel("Count")
        matplotlib.pyplot.show()

path_length_histogram(C1, title="Histogram of all shortest path (cluster 1)")
path_length_histogram(C2, title="Histogram of all shortest path (cluster 2)")
path_length_histogram(C3, title="Histogram of all shortest path (cluster 3)")
path_length_histogram(C4, title="Histogram of all shortest path (cluster 4)")

    



















def EIF4F_interaction_reference(x):
    # HuRI = protein_interaction_reference()
    # construct network for all proteins interacting eIF4G1
    # find all eIF4G1 interacting proteins from reference data
    # x = ["EIF4A1"]
    EIF4F_RI = HuRI_CCLE.loc[HuRI_CCLE['protein1'].isin(x)]
    # construct a list of proteins with interaction to eIF4G1
    EIF4F_RI_List = EIF4F_RI['protein2'].append(pd.Series(x))
    # construct a network file recording all protein interactions within EIF4G_RI_list
    EIF4F_InterPro_Net = HuRI_CCLE[HuRI_CCLE['protein1'].isin(
        EIF4F_RI_List) & HuRI_CCLE['protein2'].isin(EIF4F_RI_List)]
    EIF4F_InterPro_Net.reset_index(inplace=True)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net.drop('index', 1)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net[['protein1',
                                             'protein2',
                                             'experimental',
                                             'database']]
    # Remove reverse duplicates from dataframe
    cols = ['protein1', 'protein2']
    EIF4F_InterPro_Net[cols] = np.sort(EIF4F_InterPro_Net[cols].values, axis=1)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net.drop_duplicates()
    EIF4F_InterPro_Net = EIF4F_InterPro_Net[(EIF4F_InterPro_Net['experimental'] != 0) |
                                            (EIF4F_InterPro_Net['database'] != 0)]
    melted_df = pd.melt(EIF4F_InterPro_Net,
                        id_vars=['protein1', 'protein2'],
                        value_vars=['experimental', 'database'],
                        var_name='sources',
                        value_name='scores')
    melted_df.loc[melted_df.sources == 'experimental', 'color'] = "r"
    melted_df.loc[melted_df.sources == 'database', 'color'] = "b"
    melted_df_r = melted_df[melted_df['color'] == "r"]
    melted_df_r = melted_df_r[(melted_df_r['scores'] != 0)]
    melted_df_b = melted_df[melted_df['color'] == "b"]
    melted_df_b = melted_df_b[(melted_df_b['scores'] != 0)]
    return (EIF4F_RI_List, EIF4F_InterPro_Net, melted_df_r, melted_df_b)
# uniques1, uniques2, EIF4F_InterPro_Net = EIF4F_interaction_reference (["EIF4G1"])


def plot_EIF4F_ref_network(x):
    # uniques1, uniques2, EIF4F_InterPro_Net = EIF4F_interaction_reference (x)
    ## plot nodes interacting eIF4G1 by experiments and database
    ## plot edges for experimental interactions
    EIF4F_InterPro_Net, melted_df_r, melted_df_b = EIF4F_interaction_reference(["EIF4A1"])
    G = nx.from_pandas_edgelist(
        EIF4F_InterPro_Net,
        "protein1",
        "protein2",
        edge_attr=["experimental", "database"],
    )
    G.add_nodes_from(nodes_for_adding=EIF4F_InterPro_Net.protein1.tolist())

    labels = [i for i in dict(G.nodes).keys()]
    labels = {i: i for i in dict(G.nodes).keys()}

    fig, ax = plt.subplots(figsize=(20, 20))
    ax.set_title(x[0] + " ref network")
    pos = nx.kamada_kawai_layout(G)
    # pos = nx.spring_layout(G, seed = 50)
    nx.draw_networkx_nodes(
        G, pos, ax=ax, label=True, node_color="#cccccc", node_size=1000
    )
    nx.draw_networkx_edges(G, pos, ax=ax)
    _ = nx.draw_networkx_labels(G, pos, labels, ax=ax)
    # plt.savefig(x[0] + '_ref_network.pdf')
    ax.set_facecolor('white')
    fig.set_facecolor('white')
    plt.show()

    ## plot Draw POPULAR protein
    protein = list(EIF4F_InterPro_Net.protein2.unique())
    fig, ax = plt.subplots(figsize=(40, 40))
    # fig, ax = plt.figure(figsize = (12, 10))
    pos = nx.kamada_kawai_layout(G)
    # pos = nx.spring_layout(G, seed = 50)
    # Draw every protein
    nx.draw_networkx_nodes(
        G, pos, ax=ax, label=True, node_color="#cccccc", node_size=4000
    )
    # Draw POPULAR protein
    popular_protein = [item for item in protein if G.degree(item) > 20]
    nx.draw_networkx_nodes(
        G, pos, nodelist=popular_protein, node_color="#00b4d9", node_size=4000
    )
    nx.draw_networkx_edges(G, pos, ax=ax, width=1)
    _ = nx.draw_networkx_labels(G, pos, labels, ax=ax)
    #plt.savefig(x[0] + "_ref_network.pdf")
    ax.set_facecolor('white')
    fig.set_facecolor('white')
    plt.show()

    # plot network with experimental interaction data only
    G = nx.from_pandas_edgelist(
        melted_df_r,
        "protein1",
        "protein2",
        edge_attr=["sources", "color"],
        create_using=nx.MultiGraph(),
    )
    G.add_nodes_from(nodes_for_adding=EIF4F_InterPro_Net.protein1.tolist())

    # weights = nx.get_edge_attributes(G,'weight').values()
    labels = [i for i in dict(G.nodes).keys()]
    labels = {i: i for i in dict(G.nodes).keys()}
    pos = nx.fruchterman_reingold_layout(G)
    # pos = nx.spring_layout(G, seed = 100)
    # pos= graphviz_layout(G, prog = "neato")
    # pos = nx.kamada_kawai_layout(G)
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.set_title(x[0] + " ref network")
    nx.draw_networkx_nodes(
        G, pos, ax=ax, label=True, node_color="#cccccc", node_size=1000
    )
    nx.draw_networkx_edges(G, pos, edge_color="r", ax=ax)
    _ = nx.draw_networkx_labels(G, pos, labels, ax=ax)
    # save_name = ["EIF4A1"][0] + '_ref_network.pdf'
    # plt.savefig(x[0] + '_ref_network.pdf')
    ax.set_facecolor('white')
    fig.set_facecolor('white')
    plt.show()


# plot_EIF4F_ref_network (["EIF4A1","EIF4G1"])
plot_EIF4F_ref_network(x = ["EIF4A1"])
plot_EIF4F_ref_network(["EIF4G1"])
plot_EIF4F_ref_network(["EIF4E"])
plot_EIF4F_ref_network(["EIF4BP1"])
plot_EIF4F_ref_network(["PIK3C3"])


def EIF4F_inter_coexp(x):
    ## construct a network file recording all protein-protein interactions within cluster0
    EIF4F_InterPro_Net = HuRI[HuRI["protein1"].isin(x) & HuRI["protein2"].isin(x)]
    EIF4F_InterPro_Net.reset_index(inplace=True)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net.drop("index", 1)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net[
        ["protein1", "protein2", "experimental", "database"]
    ]
    # Remove reverse duplicates from dataframe
    cols = ["protein1", "protein2"]
    EIF4F_InterPro_Net[cols] = np.sort(EIF4F_InterPro_Net[cols].values, axis=1)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net.drop_duplicates()
    EIF4F_InterPro_Net = EIF4F_InterPro_Net[
        (EIF4F_InterPro_Net["experimental"] != 0)
        & (EIF4F_InterPro_Net["database"] != 0)
    ]
    melted_df = pd.melt(
        EIF4F_InterPro_Net,
        id_vars=["protein1", "protein2"],
        value_vars=["experimental", "database"],
        var_name="sources",
        value_name="scores",
    )
    melted_df.loc[melted_df.sources == "experimental", "color"] = "r"
    melted_df.loc[melted_df.sources == "database", "color"] = "b"
    melted_df_r = melted_df[melted_df["color"] == "r"]
    melted_df_r = melted_df_r[(melted_df_r["scores"] != 0)]
    melted_df_b = melted_df[melted_df["color"] == "b"]
    melted_df_b = melted_df_b[(melted_df_b["scores"] != 0)]
    return (EIF4F_InterPro_Net, melted_df_r, melted_df_b)


# uniques1, uniques2, EIF4F_InterPro_Net = EIF4F_inter_coexp (cluster2)


def plot_EIF4F_coexp_network(x, y):
    x = pd.Series(eIF4G1_pos.index).apply(lambda x: x.split(" ")[0])
    EIF4F_Coexp_Net, melted_df_r, melted_df_b = EIF4F_inter_coexp(x)
    ## plot nodes interacting eIF4G1 by experiments and database
    ## plot edges for experimental interactions
    G = nx.from_pandas_edgelist(
        EIF4F_Coexp_Net, "protein1", "protein2", edge_attr=["experimental", "database"]
    )
    G.add_nodes_from(nodes_for_adding=EIF4F_Coexp_Net.protein1.tolist())

    protein = list(EIF4F_Coexp_Net.protein2.unique())

    labels = [i for i in dict(G.nodes).keys()]
    labels = {i: i for i in dict(G.nodes).keys()}

    fig, ax = plt.subplots(figsize=(80, 80))
    # fig, ax = plt.figure(figsize = (12, 10))
    pos = nx.kamada_kawai_layout(G)
    # pos = nx.spring_layout(G, seed = 50)
    # Draw every protein
    nx.draw_networkx_nodes(
        G, pos, ax=ax, label=True, node_color="#cccccc", node_size=4000
    )
    # Draw POPULAR protein
    popular_protein = [item for item in protein if G.degree(item) > 20]
    nx.draw_networkx_nodes(
        G, pos, nodelist=popular_protein, node_color="orange", node_size=4000
    )
    nx.draw_networkx_edges(G, pos, ax=ax, width=1, edge_color="#cccccc")
    _ = nx.draw_networkx_labels(G, pos, labels, font_size=16, ax=ax)
    plt.savefig(y + "_corr_network_graph.pdf")
    plt.show()

    ## plot nodes interacting eIF4G1 by experiment
    ## plot edges of experimental interaction.
    G = nx.from_pandas_edgelist(
        melted_df_r,
        "protein1",
        "protein2",
        edge_attr=["sources", "color"],
        create_using=nx.MultiGraph(),
    )
    G.add_nodes_from(nodes_for_adding=melted_df_r.protein1.tolist())

    # weights = nx.get_edge_attributes(G,'weight').values()
    labels = [i for i in dict(G.nodes).keys()]
    labels = {i: i for i in dict(G.nodes).keys()}
    # pos = nx.kamada_kawai_layout(G)
    # pos = nx.spring_layout(G, seed = 50)
    pos = graphviz_layout(G, prog="neato")
    fig, ax = plt.subplots(figsize=(100, 100))
    nx.draw_networkx_nodes(
        G,
        pos,
        ax=ax,
        label=True,
        node_color="lightgreen",
        node_size=5000,
        # cmap = plt.cm.Blues
    )
    nx.draw_networkx_edges(
        G,
        pos,
        # edge_color = "r",
        ax=ax,
    )
    _ = nx.draw_networkx_labels(G, pos, labels, font_size=18, font_weight="bold", ax=ax)
    plt.savefig(y + "_corr_graphviz.pdf")
    plt.show()


plot_EIF4F_coexp_network(cluster2, "cluster2")
list(map(plot_EIF4F_coexp_network, [cluster0, cluster1, cluster2, cluster3]))
plot_EIF4F_coexp_network(pd.Series(eIF4G1_pos.index).apply(lambda x: x.split(" ")[0]))
plot_EIF4F_coexp_network(pd.Series(eIF4A1_pos.index).apply(lambda x: x.split(" ")[0]))
plot_EIF4F_coexp_network(pd.Series(eIF4E_pos.index).apply(lambda x: x.split(" ")[0]))
plot_EIF4F_coexp_network(pd.Series(EBP1_pos.index).apply(lambda x: x.split(" ")[0]))


def plot_protein_corr_network(x):
    protein_pos, protein_neg = eif_corr_sig_p(x)
    protein_pos = pd.Series(protein_pos.index).apply(lambda x: x.split(" ")[0])
    plot_EIF4F_coexp_network(protein_pos, x)


plot_protein_corr_network("EIF4G1 Q04637-9")
plot_protein_corr_network("EIF4A1 P60842")
plot_protein_corr_network("EIF4E P06730-2")


def plot_combined_corr_network():
    EIF2 = pd.Series(["EIF4G1 Q04637-9", "EIF4A1 P60842", "EIF4E P06730-2"])
    s1 = pd.Series()
    for gene in EIF1:
        protein_pos, protein_neg = eif_corr_sig_p(gene)
        protein_pos = pd.Series(protein_pos.index).apply(lambda x: x.split(" ")[0])
        s1 = s1.append(protein_pos)
    s2 = s1.drop_duplicates()
    EIF4F_Coexp_Net, melted_df_r, melted_df_b = EIF4F_inter_coexp(s2)
    s3 = EIF4F_Coexp_Net.protein1.drop_duplicates()
    # extract only eIF4g1 interacting
    EIF4G1_protein_pos, protein_neg = eif_corr_sig_p("EIF4G1 Q04637-9")
    EIF4G1_protein_pos = pd.Series(EIF4G1_protein_pos.index).apply(
        lambda x: x.split(" ")[0]
    )
    EIF4G1_Coexp_Net = s3[s3.isin(EIF4G1_protein_pos)]

    # extract only eIF4a1 interacting
    EIF4A1_protein_pos, protein_neg = eif_corr_sig_p("EIF4A1 P60842")
    EIF4A1_protein_pos = pd.Series(EIF4A1_protein_pos.index).apply(
        lambda x: x.split(" ")[0]
    )
    EIF4A1_Coexp_Net = s3[s3.isin(EIF4A1_protein_pos)]
    # extract only eIF4a1 interacting
    EIF4E_protein_pos, protein_neg = eif_corr_sig_p("EIF4E P06730-2")
    EIF4E_protein_pos = pd.Series(EIF4E_protein_pos.index).apply(
        lambda x: x.split(" ")[0]
    )
    EIF4E_Coexp_Net = s3[s3.isin(EIF4E_protein_pos)]
    # extract only eIF4a1 interacting
    EIF4EBP1_protein_pos, protein_neg = eif_corr_sig_p("EIF4EBP1 Q13541")
    EIF4EBP1_protein_pos = pd.Series(EIF4EBP1_protein_pos.index).apply(
        lambda x: x.split(" ")[0]
    )
    EIF4EBP1_Coexp_Net = s3[s3.isin(EIF4EBP1_protein_pos)]
    ## plot nodes interacting eIF4G1 by experiments and database
    ## plot edges for experimental interactions
    G = nx.from_pandas_edgelist(
        EIF4F_Coexp_Net, "protein1", "protein2", edge_attr=["experimental", "database"]
    )
    G.add_nodes_from(nodes_for_adding=s3)

    protein = list(EIF4F_Coexp_Net.protein2.unique())

    labels = [i for i in dict(G.nodes).keys()]
    labels = {i: i for i in dict(G.nodes).keys()}

    fig, ax = plt.subplots(figsize=(100, 100))
    # fig, ax = plt.figure(figsize = (12, 10))
    pos = nx.kamada_kawai_layout(G)
    # pos = nx.spring_layout(G, seed = 50)

    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=EIF4A1_Coexp_Net,
        node_color="cyan",  # alpha = 0.5,
        node_size=4000,
    )
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=EIF4G1_Coexp_Net,
        node_color="yellow",
        alpha=0.5,
        node_size=4000,
    )
    nx.draw_networkx_nodes(
        G, pos, nodelist=EIF4E_Coexp_Net, node_color="coral", alpha=0.5, node_size=4000
    )
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=EIF4EBP1_Coexp_Net,
        node_color="violet",
        alpha=0.5,
        node_size=4000,
    )
    nx.draw_networkx_edges(G, pos, ax=ax, width=1, edge_color="#cccccc")
    _ = nx.draw_networkx_labels(G, pos, labels, font_weight="bold", font_size=16, ax=ax)
    plt.savefig("all_corr_network_graph.pdf")
    plt.show()

    def rescale(l, newmin, newmax):
        arr = list(l)
        return [
            (x - min(arr)) / (max(arr) - min(arr)) * (newmax - newmin) + newmin
            for x in arr
        ]

    bc = nx.betweenness_centrality(G)  # betweeness centrality
    s = rescale([v for v in bc.values()], 1500, 7000)

    components = nx.connected_components(G)
    # compare among components and find the one having maximun length(LC)
    largest_component = max(components, key=len)
    # largest_component
    # Q1.draw LC
    subgraph = G.subgraph(largest_component)
    # pos = nx.spring_layout(subgraph) # force nodes to separte
    # pos= graphviz_layout(G, prog = "neato")
    pos = nx.kamada_kawai_layout(G)
    betCent = nx.betweenness_centrality(subgraph, normalized=True, endpoints=True)
    s = [v * 10000 for v in betCent.values()]


def EIF4F_inter_coexp_sub(x, y):
    ## construct a network file recording all protein-protein interactions within cluster0
    EIF4F_InterPro = HuRI[HuRI["protein1"].isin(x) & HuRI["protein2"].isin(y)]
    EIF4F_InterPro = EIF4F_InterPro.dropna()
    ## construct a list of proteins with interaction to eIF4G1
    EIF4F_InterPro_List = EIF4F_InterPro["protein2"].append(pd.Series(x))
    ## construct a network file recording all protein interactions within EIF4G_RI_list
    EIF4F_InterPro_Net = HuRI[
        HuRI["protein1"].isin(EIF4F_InterPro_List)
        & HuRI["protein2"].isin(EIF4F_InterPro_List)
    ]

    EIF4F_InterPro_Net.reset_index(inplace=True)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net.drop("index", 1)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net[
        ["protein1", "protein2", "experimental", "database"]
    ]
    # Remove reverse duplicates from dataframe
    cols = ["protein1", "protein2"]
    EIF4F_InterPro_Net[cols] = np.sort(EIF4F_InterPro_Net[cols].values, axis=1)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net.drop_duplicates()
    EIF4F_InterPro_Net = EIF4F_InterPro_Net[
        (EIF4F_InterPro_Net["experimental"] != 0)
        & (EIF4F_InterPro_Net["database"] != 0)
    ]

    melted_df = pd.melt(
        EIF4F_InterPro_Net,
        id_vars=["protein1", "protein2"],
        value_vars=["experimental", "database"],
        var_name="sources",
        value_name="scores",
    )
    melted_df.loc[melted_df.sources == "experimental", "color"] = "r"
    melted_df.loc[melted_df.sources == "database", "color"] = "b"
    melted_df_r = melted_df[melted_df["color"] == "r"]
    melted_df_r = melted_df_r[(melted_df_r["scores"] != 0)]
    melted_df_b = melted_df[melted_df["color"] == "b"]
    melted_df_b = melted_df_b[(melted_df_b["scores"] != 0)]
    return (EIF4F_InterPro_Net, melted_df_r, melted_df_b)


def plot_EIF4F_coexp_sub_network(x, y):
    # EIF4F_InterPro_Net, melted_df_r, melted_df_b = EIF4F_inter_coexp (cluster2)
    EIF4F_InterPro_Net, melted_df_r, melted_df_b = EIF4F_inter_coexp_sub(x, y)
    ## plot nodes interacting eIF4G1 by experiments and database
    ## plot edges for experimental interactions
    G = nx.from_pandas_edgelist(
        EIF4F_InterPro_Net,
        "protein1",
        "protein2",
        edge_attr=["experimental", "database"],
    )
    G.add_nodes_from(nodes_for_adding=EIF4F_InterPro_Net.protein1.tolist())
    protein = list(EIF4F_InterPro_Net.protein1.unique())

    labels = [i for i in dict(G.nodes).keys()]
    labels = {i: i for i in dict(G.nodes).keys()}

    fig, ax = plt.subplots(figsize=(20, 20))
    # pos = nx.spring_layout(G,iterations=50)
    pos = nx.kamada_kawai_layout(G)
    # pos = nx.spring_layout(G, seed = 50)
    # Draw every protein

    nx.draw_networkx_nodes(
        G, pos, ax=ax, label=True, node_color="#cccccc", node_size=100
    )
    nx.draw_networkx_nodes(G, pos, nodelist=x, node_color="orange", node_size=100)
    # Draw POPULAR protein
    popular_protein = [item for item in protein if G.degree(item) > 10]
    nx.draw_networkx_nodes(
        G, pos, nodelist=popular_protein, node_color="orange", node_size=100
    )
    nx.draw_networkx_edges(G, pos, ax=ax, width=1, edge_color="#cccccc")
    _ = nx.draw_networkx_labels(G, pos, labels, ax=ax)

    # Compute largest connected component of the network  (LC)
    # lsit the components in network (g)
    components = nx.connected_components(G)
    # compare among components and find the one having maximun length(LC)
    largest_component = max(components, key=len)
    # largest_component
    # Q1.draw LC
    subgraph = G.subgraph(largest_component)
    # pos = nx.spring_layout(subgraph) # force nodes to separte
    # pos= graphviz_layout(G, prog = "neato")
    pos = nx.kamada_kawai_layout(G)
    betCent = nx.betweenness_centrality(subgraph, normalized=True, endpoints=True)
    node_color = [20000.0 * G.degree(v) for v in subgraph]
    node_size = [v * 10000 for v in betCent.values()]
    plt.figure(figsize=(20, 15))
    nx.draw_networkx(
        subgraph, pos=pos, with_labels=False, node_color=node_color, node_size=node_size
    )
    plt.axis("off")

    ## plot nodes interacting eIF4G1 by experiment and database
    ## plot edges of experimental interaction.
    G = nx.from_pandas_edgelist(
        melted_df_r,
        "protein1",
        "protein2",
        edge_attr=["sources", "color"],
        create_using=nx.MultiGraph(),
    )
    G.add_nodes_from(nodes_for_adding=melted_df_r.protein1.tolist())
    labels = [i for i in dict(G.nodes).keys()]
    labels = {i: i for i in dict(G.nodes).keys()}
    pos = nx.kamada_kawai_layout(G)
    # pos = nx.spring_layout(G, seed = 50)
    # pos= graphviz_layout(G, prog = "neato")
    fig, ax = plt.subplots(figsize=(20, 20))
    nx.draw_networkx_nodes(G, pos, ax=ax, label=True, node_size=100, cmap=plt.cm.Blues)
    nx.draw_networkx_edges(G, pos, edge_color="r", ax=ax)
    _ = nx.draw_networkx_labels(G, pos, labels, ax=ax)


plot_EIF4F_coexp_sub_network(["EIF4G1"], cluster1)
plot_EIF4F_coexp_sub_network(["EIF4E"], cluster3)


### IRES ###
def plot_ires_heatmap():
    IRES_list = pd.read_csv("~/Downloads/human_IRES_info.csv")
    EIF_COR_PRO_sig["Gene_Symbol"] = EIF_COR_PRO_sig.index
    EIF_COR_PRO_sig["Gene_Symbol"] = EIF_COR_PRO_sig["Gene_Symbol"].apply(
        lambda x: x.split(" ")[0]
    )
    EIF_COR_PRO_sig.set_index("Gene_Symbol", inplace=True)

    ires_merge = pd.merge(EIF_COR_PRO_sig, IRES_list, how="inner", on=["Gene_Symbol"])
    ires_merge.set_index("Gene_Symbol", inplace=True)

    sbn.clustermap(
        ires_merge.iloc[:, 0:4],
        method="centroid",
        metric="euclidean",
        tree_kws=dict(linewidths=0.5, colors=(0.2, 0.2, 0.4)),
        cmap="coolwarm",
    )


plot_ires_heatmap()


def plot_tsne(data):
    tsne_obj = TSNE(n_components=2, random_state=0).fit_transform(data)
    # tsne_obj
    # tsne_em = TSNE(n_components = 2, perplexity   = 30.0, n_iter       = 1000, verbose      =1).fit_transform(data)
    tsne_df = pd.DataFrame({"X": tsne_obj[:, 0], "Y": tsne_obj[:, 1]})
    tsne_df.head()
    plt.figure(figsize=(16, 10))
    sbn.scatterplot(x="X", y="Y", data=tsne_df)


plot_tsne(EIF_COR_PRO_sig)


# Displaying dataframe as an heatmap
# with diverging colourmap as coolwarm
D = sch.distance.pdist(EIF_COR_PRO_sig, metric="euclidean")
L = sch.linkage(D, method="centroid")

sch.dendrogram(
    sch.linkage(D, method="centroid"),
    orientation="top",
    p=5,
    truncate_mode="level",
    color_threshold=0.665,
)

h = sbn.clustermap(
    EIF_COR_PRO_sig,
    figsize=(6, 10),
    method="centroid",
    metric="euclidean",
    tree_kws=dict(linewidths=0.5, colors=(0.2, 0.2, 0.4)),
    cmap="coolwarm",
)
plt.savefig("COR_heatmap.pdf")


def plot_heatmap(data):
    # clustering colums
    data_1D_X = ssd.pdist(data.T, "euclidean")
    X = sch.linkage(data_1D_X, method="centroid")
    # clustering rows
    data_1D_Y = ssd.pdist(data, "euclidean")
    Y = sch.linkage(data_1D_Y, method="centroid")
    # plot first dendrogram
    fig = plt.figure(figsize=(8, 8))

    # sch.set_link_color_palette(['m', 'c', 'y', 'k','g','r'])
    ax1 = fig.add_axes([0.09, 0.1, 0.2, 0.6])
    Z1 = sch.dendrogram(Y, orientation="left", color_threshold=0.65)
    ax1.set_xticks([])
    ax1.set_yticks([])

    # second dendrogram.
    ax2 = fig.add_axes([0.3, 0.71, 0.6, 0.1])
    Z2 = sch.dendrogram(X)
    ax2.set_xticks([])
    ax2.set_yticks([])

    # plot matrix
    axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.6])
    # sorts based of clustering
    idx1 = Z1["leaves"]
    idx2 = Z2["leaves"]
    D = data.values[idx1, :]
    D = D[:, idx2]
    im = axmatrix.matshow(D, aspect="auto", origin="lower", cmap="coolwarm")
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])


plot_heatmap(EIF_COR_PRO_sig)


###########################################
## load RNA-seq data for CCLE cell lines ##
###########################################
## https://portals.broadinstitute.org/ccle/data
def ccle_rna():
    ## https://depmap.org/portal/download/
    CCLE_RNA = pd.read_csv("~/Downloads/CCLE_expression_full.csv")
    ccle = CCLE_RNA.set_index("Unnamed: 0")
    ccle.index.names = ["DepMap_ID"]
    return ccle


CCLE_RNA = ccle_rna()


def eif_corr_rna(eif):
    y = pd.DataFrame()
    for gene in eif:
        print(gene)
        x = CCLE_RNA.corrwith(CCLE_RNA[gene])
        x.name = gene
        y = pd.concat([y, x], axis=1)
    return y


# gee = ['EIF4G1 (ENSG00000114867)']
# EIF_COR_RNA = eif_corr_rna(gee)
def plot_eif_corr_rna(gene):
    EIF_COR_RNA = eif_corr_rna(gene)
    EIF_COR_RNA_sig = EIF_COR_RNA.loc[
        (EIF_COR_RNA["EIF4G1 (ENSG00000114867)"] >= 0.3)
        | (EIF_COR_RNA["EIF4G1 (ENSG00000114867)"] <= -0.3)
        | (EIF_COR_RNA["EIF4A1 (ENSG00000161960)"] >= 0.3)
        | (EIF_COR_RNA["EIF4A1 (ENSG00000161960)"] <= -0.3)
        | (EIF_COR_RNA["EIF4E (ENSG00000151247)"] >= 0.3)
        | (EIF_COR_RNA["EIF4E (ENSG00000151247)"] <= -0.3)
        | (EIF_COR_RNA["EIF4EBP1 (ENSG00000187840)"] >= 0.3)
        | (EIF_COR_RNA["EIF4EBP1 (ENSG00000187840)"] <= -0.3)
    ]
    EIF_COR_RNA_sig.dropna(inplace=True)
    EIF_COR_RNA_sig.dtypes
    # EIF_COR_PRO['index1'] = EIF_COR_PRO.index
    # Displaying dataframe as an heatmap
    # with diverging colourmap as coolwarm
    sbn.clustermap(
        EIF_COR_RNA_sig,
        # metric='correlation',
        # standard_scale=1,
        cmap="coolwarm",
    )


EIF = pd.Index(
    [
        "EIF4G1 (ENSG00000114867)",
        "EIF4A1 (ENSG00000161960)",
        "EIF4E (ENSG00000151247)",
        "EIF4EBP1 (ENSG00000187840)",
    ]
)
plot_eif_corr_rna(EIF)

# pd.Series(CCLE_PRO["Gene_Symbol"]).is_unique
# pd.Series(CCLE_PRO_subset.index).is_unique

# ids = CCLE_PRO["Gene_Symbol"]
# du = CCLE_PRO[ids.isin(ids[ids.duplicated()])].sort_values("Gene_Symbol")


def gdsc_drugrespose():
    ## http://www.cancerrxgene.org/downloads
    GDSC_drugrespose = pd.read_csv(
        "~/Downloads/PANCANCER_IC_Mon Nov  9 16_46_39 2020.csv"
    )
    GDSC_drugrespose.dtypes  # The data type of each column.
    # convert "Drug name" from object to category
    # GDSC_drugrespose['Drug name'] = GDSC_drugrespose['Drug name','Cell line name'].astype('category')
    for col in ["Drug name", "Drug Id", "Cell line name", "Tissue"]:
        GDSC_drugrespose[col] = GDSC_drugrespose[col].astype("category")
    return GDSC_drugrespose


GDSC_drugrespose = gdsc_drugrespose()
drug_list = GDSC_drugrespose["Drug name"].cat.categories


def cpd_corr_pro(a_drug_list, protein):
    CCLE_PRO_EIF = CCLE_PRO_subset[EIF]
    CCLE_PRO_EIF.index = CCLE_PRO_EIF.index.str.split("_").str[0]
    CCLE_PRO_EIF.index.names = ["Cell line name"]
    y = pd.DataFrame()
    for drug in a_drug_list:
        print(drug)
        ## select one drug type
        GDSC_drugrespose_subset = GDSC_drugrespose.loc[
            GDSC_drugrespose["Drug name"] == drug
        ]
        GDSC_drugrespose_subset["Cell line name"] = GDSC_drugrespose_subset[
            "Cell line name"
        ].str.replace("-", "", regex=True)
        GDSC_drugrespose_subset.set_index("Cell line name", inplace=True)
        GDSC_drugrespose_subset = GDSC_drugrespose_subset.drop(
            GDSC_drugrespose_subset.columns[[0, 1, 2, 3, 4, 5]], axis=1
        )
        GDSC_drugrespose_subset = GDSC_drugrespose_subset.drop(
            ["Max conc", "RMSE", "Dataset version", "Z score", "AUC"], 1
        )
        # select one gene expression
        # CCLE_expression_subset = CCLE_expression.loc[CCLE_expression['Description'] == 'EIF4G1']
        result = pd.merge(GDSC_drugrespose_subset, CCLE_PRO_EIF, on="Cell line name")
        x = result.corrwith(result["IC50"])
        x.name = drug
        y = pd.concat([y, x], axis=1)
    return y


# cpd_pro_correlation = cpd_corr_pro(drug, CCLE_PRO_EIF)
CPD_COR_EIF = cpd_corr_pro(drug_list, CCLE_PRO_subset).T
# CPD_COR_EIF['indx'] = CPD_COR_EIF.index
CPD_COR_EIF = CPD_COR_EIF.drop(["IC50"], 1)
CPD_COR_EIF_sig = CPD_COR_EIF.loc[
    (CPD_COR_EIF["EIF4G1 Q04637-9"] >= 0.3)
    | (CPD_COR_EIF["EIF4G1 Q04637-9"] <= -0.3)
    # |(EIF_COR_PRO['EIF4G1 Q04637-8'] >= 0.3) |  (EIF_COR_PRO['EIF4G1 Q04637-8'] <= -0.3)
    | (CPD_COR_EIF["EIF4A1 P60842"] >= 0.3)
    | (CPD_COR_EIF["EIF4A1 P60842"] <= -0.3)
    | (CPD_COR_EIF["EIF4E P06730-2"] >= 0.3)
    | (CPD_COR_EIF["EIF4E P06730-2"] <= -0.3)
    | (CPD_COR_EIF["EIF4EBP1 Q13541"] >= 0.3)
    | (CPD_COR_EIF["EIF4EBP1 Q13541"] <= -0.3)
    # | (CPD_COR_EIF['MKNK2 Q9HBH9'] >= 0.3) |  (CPD_COR_EIF['MKNK2 Q9HBH9'] <= -0.3)
    # | (CPD_COR_EIF['MKNK1 Q9BUB5'] >= 0.3) |  (CPD_COR_EIF['MKNK1 Q9BUB5'] <= -0.3)
]
CPD_COR_EIF_sig.dropna(inplace=True)
sbn.clustermap(
    CPD_COR_EIF_sig,
    method="centroid",
    metric="euclidean",
    tree_kws=dict(linewidths=0.5, colors=(0.2, 0.2, 0.4)),
    cmap="coolwarm",
)


def ccle_rna_annotation():
    ## load annontation data for CCLE cell lines
    ## https://depmap.org/portal/download/
    CCLE_annotation = pd.read_csv("~/Downloads/sample_info.csv")
    CCLE_annotation_subset = CCLE_annotation[["DepMap_ID", "stripped_cell_line_name"]]
    CCLE_annotation_subset.set_index("DepMap_ID", inplace=True)
    CCLE_RNA_annotation = pd.merge(CCLE_RNA, CCLE_annotation_subset, on="DepMap_ID")
    CCLE_RNA_annotation.set_index("stripped_cell_line_name", inplace=True)
    CCLE_RNA_annotation.index.names = ["Cell line name"]
    return CCLE_RNA_annotation


CCLE_RNA_annotation = ccle_rna_annotation()
CCLE_RNA_EIF = CCLE_RNA_annotation[
    [
        "EIF4G1 (ENSG00000114867)",
        "EIF4A1 (ENSG00000161960)",
        "EIF4E (ENSG00000151247)",
        "EIF4EBP1 (ENSG00000187840)",
    ]
]


def cpd_corr_rna(a_drug_list, rna):
    y = pd.DataFrame()
    for drug in a_drug_list:
        print(drug)
        GDSC_drugrespose_subset = GDSC_drugrespose.loc[
            GDSC_drugrespose["Drug name"] == drug
        ]
        GDSC_drugrespose_subset["Cell line name"] = GDSC_drugrespose_subset[
            "Cell line name"
        ].str.replace("-", "", regex=True)
        GDSC_drugrespose_subset.set_index("Cell line name", inplace=True)
        GDSC_drugrespose_subset = GDSC_drugrespose_subset.drop(
            GDSC_drugrespose_subset.columns[[0, 1, 2, 3, 4, 5]], axis=1
        )
        GDSC_drugrespose_subset = GDSC_drugrespose_subset.drop(
            ["Max conc", "RMSE", "Dataset version"], 1
        )
        # select one gene expression
        # CCLE_expression_subset = CCLE_expression.loc[CCLE_expression['Description'] == 'EIF4G1']
        result = pd.merge(GDSC_drugrespose_subset, rna, on="Cell line name")
        x = result.corrwith(result["IC50"])
        x.name = drug
        y = pd.concat([y, x], axis=1)
    return y


cpd_rna_correlation = cpd_corr_rna(drug_list, CCLE_RNA_EIF).T
cpd_rna_correlation = cpd_corr_rna(["Camptothecin"], CCLE_RNA_EIF).T


EIF4F_corr = cpd_rna_correlation.loc[cpd_rna_correlation.index.isin(EIF4F["Name"])]
EIF4F_corr_T = EIF4F_corr.T
cpd_rna_correlation["index1"] = cpd_rna_correlation.index


# select one drug type
GDSC_drugrespose_subset = GDSC_drugrespose.loc[
    GDSC_drugrespose["Drug name"] == "Camptothecin"
]
GDSC_drugrespose_subset["Cell line name"] = GDSC_drugrespose_subset[
    "Cell line name"
].str.replace("-", "", regex=True)
GDSC_drugrespose_subset.set_index("Cell line name", inplace=True)
GDSC_drugrespose_subset = GDSC_drugrespose_subset.drop(
    GDSC_drugrespose_subset.columns[[0, 1, 2, 3, 4, 5]], axis=1
)
GDSC_drugrespose_subset = GDSC_drugrespose_subset.drop(
    ["Max conc", "RMSE", "Dataset version"], 1
)

# select one gene expression
# CCLE_expression_subset = CCLE_expression.loc[CCLE_expression['Description'] == 'EIF4G1']
CCLE_expression = CCLE_expression.set_index("Name")
CCLE_expression = CCLE_expression.drop(["Description"], 1)
CCLE_expression_t = CCLE_expression.T
CCLE_expression_t.dtypes
CCLE_expression_t.index = CCLE_expression_t.index.str.split("_").str[0]
CCLE_expression_t.index.names = ["Cell line name"]
result = pd.merge(GDSC_drugrespose_subset, CCLE_expression_t, on="Cell line name")
y = result.corrwith(result["IC50"])
y.name = "Camptothecin"
y = y.to_frame()


new_header = CCLE_expression_transposed.iloc[0]  # grab the first row for the header
CCLE_expression_transposed = CCLE_expression_transposed[
    1:
]  # take the data less the header row
CCLE_expression_transposed.columns = new_header  # set the header row as the df header
CCLE_expression_transposed_drop = CCLE_expression_transposed.drop(index="Description")
CCLE_expression_transposed_drop_reset = CCLE_expression_transposed_drop.reset_index()
CCLE_expression_transposed_drop_reset["Cell line name"] = (
    CCLE_expression_transposed_drop_reset["index"].str.split("_").str[0]
)
CCLE_expression_transposed_drop_reset_rename = CCLE_expression_transposed_drop_reset.drop(
    ["index"], axis=1
)
x = CCLE_expression_transposed_drop_reset_rename["Cell line name"]
y = list(CCLE_expression_transposed_drop.index.values)

result = pd.merge(
    GDSC_drugrespose_subset,
    CCLE_expression_transposed_drop_reset_rename,
    on="Cell line name",
)
result_1 = result.drop(result.columns[[0, 1, 2, 3, 4, 5, 6]], axis=1)
result_1 = result_1.drop(["Max conc", "RMSE", "Dataset version"], 1)

print(result_1[["ENSG00000210195.2", "ENSG00000210196.2"]].corr(result_1["IC50"]))

result_1.dtypes
