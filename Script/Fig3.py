import datatable
import os
import lifelines
import matplotlib
import mygene
import networkx
import numpy
import pandas
import seaborn 
import scipy
import subprocess

# critical parameter setting!
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams['figure.figsize'] = 15, 15
pandas.options.mode.chained_assignment = None  # default='warn'

# setup current direcotry
data_file_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_data"
output_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_output"

# retrieve published proteomics data ######
def ccle_pro():
    # download the proteomics data
    # https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz
    ccle_pro = pandas.read_csv(os.path.join(data_file_directory, 
                                            "protein_quant_current_normalized.csv"))
    # concatenate the Gene_Symbol (with duplicate names) and the Uniprot_Acc
    ccle_pro["Gene_Symbol"] = (ccle_pro["Gene_Symbol"].fillna("") + " " + ccle_pro["Uniprot_Acc"])
    ccle_pro = ccle_pro.set_index("Gene_Symbol")
    ccle_pro_subset = ccle_pro[ccle_pro.columns.drop(list(ccle_pro.filter(regex="Peptides")))]
    ccle_pro_subset = ccle_pro_subset.drop(columns=["Protein_Id", 
                                                    "Description", 
                                                    "Group_ID", 
                                                    "Uniprot", 
                                                    "Uniprot_Acc"]).T
    ccle_pro_subset.index.names = ["ccle_name"]
    return ccle_pro_subset

CCLE_PRO_subset = ccle_pro()

ccle_anno = pandas.read_csv(os.path.join(data_file_directory, "sample_info.csv"))
ccle_cnv = pandas.read_csv(os.path.join(data_file_directory, "CCLE_gene_cn.csv"))

# retrieve published depmap data ######
def dep_crispr():
  ccle_dep_crispr = pandas.read_csv(os.path.join(data_file_directory, "CRISPR_gene_effect.csv"))
  ccle_dep_crispr.dtypes
  ccle_dep_crispr_median = ccle_dep_crispr.median(numeric_only=True)
  ccle_dep_crispr_median = pandas.DataFrame(data=ccle_dep_crispr.median(numeric_only=True))
  ccle_dep_crispr_median.index = ccle_dep_crispr_median.index.str.split(" ").str[0]
  return(ccle_dep_crispr, ccle_dep_crispr_median)

ccle_dep_crispr, ccle_dep_crispr_median = dep_crispr()

def dep_rnai():
  ccle_dep_rnai = pandas.read_csv(os.path.join(data_file_directory, "D2_combined_gene_dep_scores.csv")).T
  new_header = ccle_dep_rnai.iloc[0] #grab the first row for the header
  ccle_dep_rnai = ccle_dep_rnai[1:] #take the data less the header row
  ccle_dep_rnai.columns = new_header #set the header row as the df header
  ccle_dep_rnai = ccle_dep_rnai.astype('float64')
  return(ccle_dep_rnai)

ccle_dep_rnai = dep_rnai()

##
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
    matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), "Fig3",
    "CCLE_scatter.pdf"), dpi=300)
    # scatter plot of two protein expression across cell lines with color
    g = seaborn.PairGrid(CCLE_EIF_PRO, hue="ccle", diag_sharey=False, corner = True)
    g.map_lower(seaborn.scatterplot)
    g.map_diag(seaborn.histplot)
    g.map_upper(reg_coef)
    g.tight_layout
    g.add_legend()
    matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), "Fig3",
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
EIF_COR_PRO_sig.to_csv(os.path.join(os.path.expanduser(output_directory),"ProcessedData","EIF_COR_PRO_sig.csv"))


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
matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), "Fig3",
"COR_heatmap.pdf"), 
dpi=300)


##### Call R.script to perform heatmap clustering and pathway analysis #####
subprocess.call("Rscript /home/suwu/github/pyEIF/Script/Fig3.R", shell=True)



##############################
###### network analysis ######
##############################
# prepare the reference network data
def protein_interaction_reference():
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
  dict1 = pandas.Series(node1.symbol.values,index=node1.index).to_dict()
  dict2 = pandas.Series(node2.symbol.values,index=node2.index).to_dict()
  HuRI['protein1'] = HuRI['protein1'].map(dict1)
  HuRI['protein2'] = HuRI['protein2'].map(dict2)
  HuRI = HuRI.dropna()  # some na in 'protein2' column, remove them
  # chose experimental and database values over 1
  HuRI = HuRI[(HuRI['experimental'] != 0) & (HuRI['database'] != 0)]
  CCLE_pro = pandas.Series(CCLE_PRO_subset.columns).apply(lambda x: x.split(' ')[0])
  HuRI_CCLE = HuRI[HuRI['protein1'].isin(CCLE_pro) & HuRI['protein2'].isin(CCLE_pro)]
  HuRI_CCLE = HuRI_CCLE.dropna()
  #HuRI_CCLE[HuRI_CCLE['protein1'] != HuRI_CCLE['protein2']]
  return(HuRI, HuRI_CCLE)

HuRI, HuRI_CCLE = protein_interaction_reference()

# retrieve the clustering genes
EIF_COR_sig = pandas.DataFrame(EIF_COR_PRO_sig.index.str.split(" ").str[0]).dropna()
EIF_COR_sig.columns = ['gene']

cluster1 = pandas.read_csv(os.path.join(os.path.expanduser(output_directory),"ProcessedData", "cluster1.csv")) 
cluster2 = pandas.read_csv(os.path.join(os.path.expanduser(output_directory),"ProcessedData", "cluster2.csv"))
cluster3 = pandas.read_csv(os.path.join(os.path.expanduser(output_directory),"ProcessedData", "cluster3.csv")) 
cluster4 = pandas.read_csv(os.path.join(os.path.expanduser(output_directory),"ProcessedData", "cluster4.csv")) 

cluster_Net = HuRI_CCLE[HuRI_CCLE['protein1'].isin(EIF_COR_sig['gene']) & HuRI_CCLE['protein2'].isin(EIF_COR_sig['gene'])]
cluster1_Net = HuRI_CCLE[HuRI_CCLE['protein1'].isin(cluster1['gene']) & HuRI_CCLE['protein2'].isin(cluster1['gene'])]
cluster2_Net = HuRI_CCLE[HuRI_CCLE['protein1'].isin(cluster2['gene']) & HuRI_CCLE['protein2'].isin(cluster2['gene'])]
cluster3_Net = HuRI_CCLE[HuRI_CCLE['protein1'].isin(cluster3['gene']) & HuRI_CCLE['protein2'].isin(cluster3['gene'])]
cluster4_Net = HuRI_CCLE[HuRI_CCLE['protein1'].isin(cluster4['gene']) & HuRI_CCLE['protein2'].isin(cluster4['gene'])]


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

#pos = networkx.kamada_kawai_layout(G)
posC = networkx.kamada_kawai_layout(C)

##
def plot_combined_cluster_net():
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
  matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), "Fig3",'CORs_network_kk.pdf'),
                    dpi=300,
                    edgecolor='w',
                    #transparent=True,
                    bbox_inches='tight')   
  matplotlib.pyplot.show()
plot_combined_cluster_net()

##
def plot_cluster_net_depscore(cluster, label, pos, edge_color, node_color):
  df1 = (pandas.DataFrame(list(cluster.degree), columns=['node','degree']).set_index('node'))
  cluster1_dep = pandas.merge(ccle_dep_crispr_median, 
                              df1, 
                              left_index=True, 
                              right_index=True)
  cluster1_dep.index.names = ['node']
  # Degree Centrality: the number of edges a node has
  # Important nodes have many connections
  degree = networkx.degree_centrality(cluster)
  degree_names, degree_ranks = zip(*degree.items())
  degree_df = pandas.DataFrame(data={'protein': degree_names, 'rank': degree_ranks})
  degree_df1 = degree_df[degree_df['protein'].isin(cluster1_dep.index)]
  cluster_degree_depscore = pandas.merge(cluster1_dep, 
                              degree_df1, 
                              left_index=True, 
                              right_on='protein')
  # position of labels above nodes
  pos_higher = {}
  for k, v in pos.items():
    if(v[1]>0):
      pos_higher[k] = (v[0]-0.015, v[1]+0.015)
    else:
      pos_higher[k] = (v[0]-0.015, v[1]-0.015)
  
  # color the plot by clusters
  matplotlib.pyplot.clf()  
  matplotlib.pyplot.figure(figsize=(10, 10))
  matplotlib.pyplot.axis('off')
  networkx.draw_networkx(cluster,
                             pos, 
                             width= 0.5, 
                             node_color=node_color, 
                             edge_color=edge_color,
                             node_size=(degree_df['rank'].values*100),
                             label=label,
                             with_labels = False) 
  networkx.draw_networkx_labels(cluster, pos_higher, font_size=3.5)
  matplotlib.pyplot.tight_layout()
  matplotlib.pyplot.savefig(os.path.join(
      os.path.expanduser(output_directory), "Fig3",
      (label + "_network_kk.pdf")),
      dpi=300,
      edgecolor='w',
      bbox_inches='tight')   
  matplotlib.pyplot.show()      

  # color the nodes according to their partition
  partition = community.best_partition(cluster)
  matplotlib.pyplot.clf()
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
      os.path.expanduser(output_directory), "Fig3",
      (label + "_network_community.pdf")),
      dpi=300,
      edgecolor='w',
      bbox_inches='tight')    
  matplotlib.pyplot.show()

  # color the node by depscores
  vmin = ccle_dep_crispr_median[0].min()
  vmax = ccle_dep_crispr_median[0].max()
  cmap = matplotlib.pyplot.cm.coolwarm_r
  matplotlib.pyplot.clf()
  matplotlib.pyplot.figure(figsize=(10, 10))
  matplotlib.pyplot.axis('off')
  networkx.draw_networkx(cluster, 
                         pos, 
                         width= 0.5, 
                         edge_color="#cccccc", 
                         with_labels=False, 
                         nodelist=cluster_degree_depscore['protein'],
                         node_size=cluster_degree_depscore['rank'].values*100,
                         node_color=cluster_degree_depscore[0],
                         cmap=cmap, 
                         vmin=vmin, 
                         vmax=vmax)
  # networkx.draw_networkx_labels(cluster, pos_higher, font_size=3.5)
  sm = matplotlib.pyplot.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax))
  sm.set_array([])
  sub_ax = matplotlib.pyplot.axes([0.90, 0.75, 0.02, 0.2]) 
  cbar = matplotlib.pyplot.colorbar(sm, cax=sub_ax)
  matplotlib.pyplot.tight_layout()
  matplotlib.pyplot.savefig(os.path.join(
      os.path.expanduser(output_directory), "Fig3",
      (label + "_depscore.pdf")),
      dpi=300,
      edgecolor='w',
      bbox_inches='tight')  
  matplotlib.pyplot.show()
      
##
plot_cluster_net_depscore(cluster=C1, label = "cluster 1", pos = posC, edge_color = "lightgreen", node_color = "green")
plot_cluster_net_depscore(cluster=C2, label = "cluster 2", pos = posC, edge_color = "gold", node_color = "orange")
plot_cluster_net_depscore(cluster=C3, label = "cluster 3", pos = posC, edge_color = "skyblue", node_color = "blue")
plot_cluster_net_depscore(cluster=C4, label = "cluster 4", pos = posC, edge_color = "pink", node_color = "red")


##
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

#plot_degree_distribution(cluster = G, pos= pos, label = "all proteins", edge_color = "grey", node_color = "black")    
plot_degree_distribution(cluster = C1, pos= posC, label = "cluster 1", edge_color = "lightgreen", node_color = "green")    
plot_degree_distribution(cluster = C2, pos= posC, label = "cluster 2", edge_color = "gold", node_color = "orange")  
plot_degree_distribution(cluster = C3, pos= posC, label = "cluster 3", edge_color = "skyblue", node_color = "blue")  
plot_degree_distribution(cluster = C4, pos= posC, label = "cluster 4", edge_color = "pink", node_color = "red")  


##
def plot_depscore_histogram(cluster, pos, label, edge_color, node_color):
  df1 = (pandas.DataFrame(list(cluster.degree), columns=['node','degree']).set_index('node'))
  cluster1_dep = pandas.merge(ccle_dep_crispr_median, 
                              df1, 
                              left_index=True, 
                              right_index=True)
  cluster1_dep.index.names = ['node']
  # Degree Centrality: the number of edges a node has
  # Important nodes have many connections
  degree_sequence = sorted([d for n, d in cluster.degree()], reverse=True)  # degree sequence
  degree = networkx.degree_centrality(cluster)
  degree_names, degree_ranks = zip(*degree.items())
  degree_df = pandas.DataFrame(data={'protein': degree_names, 'rank': degree_ranks})
  degree_df1 = degree_df[degree_df['protein'].isin(cluster1_dep.index)]
  cluster_degree_depscore = pandas.merge(cluster1_dep, 
                                         degree_df1, 
                                         left_index=True, 
                                         right_on='protein')
  cluster_degree_depscore.columns = ['depscore', 'degree', 'protein', 'rank']
  
  vmin = ccle_dep_crispr_median[0].min()
  vmax = ccle_dep_crispr_median[0].max()
  cmap = matplotlib.pyplot.cm.coolwarm_r
  matplotlib.pyplot.clf()
  fig, ax = matplotlib.pyplot.subplots(figsize=(10, 10))
  #seaborn.histplot(data=cluster_degree_depscore, x="depscore", fill=True, binwidth=0.035)
  cluster_degree_depscore['depscore'].hist(
    bins=80, color="royalblue",    
    grid = False,
    rwidth = 0.9)
  matplotlib.pyplot.title("Depscore histogram " + label)
  matplotlib.pyplot.ylabel("Node count")
  matplotlib.pyplot.xlabel("Depscore")
  matplotlib.pyplot.xlim([ccle_dep_crispr_median[0].max(), ccle_dep_crispr_median[0].min()])
  matplotlib.pyplot.ylim([0, 18])
  matplotlib.pyplot.axes([0.25,0.25,0.65,0.65])
  matplotlib.pyplot.axis('off')
  networkx.draw_networkx(cluster, 
                         pos, 
                         width= 0.5, 
                         edge_color="#cccccc", 
                         with_labels=False, 
                         nodelist=cluster_degree_depscore['protein'],
                         node_size=cluster_degree_depscore['rank'].values*100,
                         node_color=cluster_degree_depscore["depscore"],
                         cmap=cmap, 
                         vmin=ccle_dep_crispr_median[0].min(), 
                         vmax=ccle_dep_crispr_median[0].max())
  matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory),"Fig3", (label + "_Depscore histogram.pdf")),
                      dpi=300,
                      edgecolor='w',
                      #transparent=True,
                      bbox_inches='tight') 
  matplotlib.pyplot.show()

  matplotlib.pyplot.clf()
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
  matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), "Fig3",
  (label + "_Degree histogram.pdf")),
                      dpi=300,
                      edgecolor='w',
                      #transparent=True,
                      bbox_inches='tight') 
  matplotlib.pyplot.show()

#plot_depscore_histogram(cluster = G, pos= pos, label = "all proteins", edge_color = "grey", node_color = "black")    
plot_depscore_histogram(cluster=C1, pos=posC, label="cluster 1", edge_color="lightgreen", node_color="green")    
plot_depscore_histogram(cluster=C2, pos=posC, label="cluster 2", edge_color="gold", node_color="orange")  
plot_depscore_histogram(cluster=C3, pos=posC, label="cluster 3", edge_color="skyblue", node_color="blue")  
plot_depscore_histogram(cluster=C4, pos=posC, label="cluster 4", edge_color="pink", node_color="red")  
 

## Histogram and Density Curve for depscores
def depscore_hisplot(data, method):
  fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = matplotlib.pyplot.subplots(3, 2, figsize=(6, 10),  sharey=True)
  seaborn.histplot(data['EIF4G1 (1981)'], color="dodgerblue", fill=True, ax=ax1, kde=True)
  seaborn.rugplot(data['EIF4G1 (1981)'], color="royalblue", ax=ax1)
  seaborn.histplot(data['EIF4A1 (1973)'],  color="gold", fill=True, ax=ax2, kde=True)
  seaborn.rugplot(data['EIF4A1 (1973)'], color="royalblue", ax=ax2)
  seaborn.histplot(data['EIF4A2 (1974)'], color="green", fill=True, ax=ax3, kde=True)
  seaborn.rugplot(data['EIF4A2 (1974)'], color="royalblue", ax=ax3)
  seaborn.histplot(data['EIF4E (1977)'], color="grey", fill=True, ax=ax4, kde=True)
  seaborn.rugplot(data['EIF4E (1977)'], color="royalblue", ax=ax4)
  seaborn.histplot(data['EIF3E (3646)'],color="deeppink", fill=True, ax=ax5, kde=True)
  seaborn.rugplot(data['EIF3E (3646)'], color="royalblue", ax=ax5)
  seaborn.histplot(data['EIF3H (8667)'], color="cyan", fill=True, ax=ax6, kde=True)
  seaborn.rugplot(data['EIF3H (8667)'], color="royalblue", ax=ax6)
  # remove axis labels
  ax1.set(xlabel=None, ylabel=None)
  ax2.set(xlabel=None, ylabel=None)
  ax3.set(xlabel=None, ylabel=None)
  ax4.set(xlabel=None, ylabel=None)
  ax5.set(xlabel=None, ylabel=None)
  ax6.set(xlabel=None, ylabel=None)
  # Set common labels
  fig.text(0.5, 0.0, ('DepMap Score (' + method + ')'), ha='center')
  fig.text(0.0, 0.5, 'Count', va='center', rotation='vertical')
  ax1.legend(['EIF4G1'],frameon=False)
  ax2.legend(['EIF4A1'],frameon=False)
  ax3.legend(['EIF4A2'],frameon=False)
  ax4.legend(['EIF4E'],frameon=False)
  ax5.legend(['EIF3E'],frameon=False)
  ax6.legend(['EIF3H'],frameon=False)
  ax1.axvline(-1,color='red',ls='--')
  ax1.axvline(0,color='grey',ls='--')
  ax2.axvline(-1,color='red',ls='--')
  ax2.axvline(0,color='grey',ls='--')
  ax3.axvline(-1,color='red',ls='--')
  ax3.axvline(0,color='grey',ls='--')
  ax4.axvline(-1,color='red',ls='--')
  ax4.axvline(0,color='grey',ls='--')
  ax5.axvline(-1,color='red',ls='--')
  ax5.axvline(0,color='grey',ls='--')
  ax6.axvline(-1,color='red',ls='--')
  ax6.axvline(0,color='grey',ls='--')
  # Setting the values for all axes.
  matplotlib.pyplot.setp(((ax1, ax2, ax3), (ax4, ax5, ax6)),  xlim=(-3, 1))
  matplotlib.pyplot.tight_layout()
  matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory),"Fig3", (method + "_dep.pdf")),
                        dpi=300,
                        edgecolor='w',
                        #transparent=True,
                        bbox_inches='tight') 
  matplotlib.pyplot.show()

depscore_hisplot(data=ccle_dep_crispr, method="CRISPR")
depscore_hisplot(data=ccle_dep_rnai, method="RNAi")


matplotlib.pyplot.clf()
matplotlib.pyplot.cla()
matplotlib.pyplot.close()


### depmap sensitivity to RNAi 
def plot_dep_cnv(protein):
  result = pandas.merge(ccle_anno[['DepMap_ID', 'CCLE_Name']], 
                      ccle_cnv[['Unnamed: 0',protein]], 
                      left_on='DepMap_ID', 
                      right_on='Unnamed: 0',
                      suffixes=('_dep', '_cnv'))
  dep_cnv = pandas.merge(ccle_dep_rnai[[protein]], 
                       result[[protein, 'CCLE_Name']], 
                       left_index=True, 
                       right_on='CCLE_Name', 
                       suffixes=('_dep', '_cnv')).dropna()

  g = seaborn.JointGrid(data=dep_cnv, x=(protein+"_dep"), y=(protein+"_cnv"))
  g.plot_joint(seaborn.regplot, fit_reg = True, color="royalblue")
  g.plot_marginals(seaborn.histplot, kde=True, bins=80, color="royalblue")
  r, p = scipy.stats.pearsonr(dep_cnv[protein+"_dep"], dep_cnv[protein+"_cnv"])
  g.ax_joint.annotate(f'$\\rho = {r:.3f}, p = {p:.3f}$',
                      xy=(0.1, 0.9), xycoords='axes fraction',
                      ha='left', va='center')
  g.ax_joint.axvline(x=-1)
  matplotlib.pyplot.show()


plot_dep_cnv(protein='EIF4G1 (1981)')
plot_dep_cnv(protein='EIF3E (3646)')



