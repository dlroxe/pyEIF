ccle_dep_crispr = pandas.read_csv(os.path.join(data_file_directory, "CRISPR_gene_effect.csv"))
ccle_dep_crispr.dtypes
ccle_dep_crispr_median = ccle_dep_crispr.median(numeric_only=True)
ccle_dep_crispr_median = pandas.DataFrame(data=ccle_dep_crispr.median(numeric_only=True))
ccle_dep_crispr_median.index = ccle_dep_crispr_median.index.str.split(" ").str[0]

ccle_dep_rnai = pandas.read_csv(os.path.join(data_file_directory, "D2_combined_gene_dep_scores.csv")).T
new_header = ccle_dep_rnai.iloc[0] #grab the first row for the header
ccle_dep_rnai = ccle_dep_rnai[1:] #take the data less the header row
ccle_dep_rnai.columns = new_header #set the header row as the df header

## Histogram and Density Curve
#matplotlib.pyplot.figure(figsize=(6, 10))
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
  matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), (method + "_dep.pdf")),
                        dpi=300,
                        edgecolor='w',
                        #transparent=True,
                        bbox_inches='tight') 
  matplotlib.pyplot.show()

depscore_hisplot(data=ccle_dep_crispr, method="CRISPR")
depscore_hisplot(data=ccle_dep_rnai, method="RNAi")


EIF_cols = [col for col in ccle_dep_score.columns if "EIF4A2" in col]
print(EIF_cols)

matplotlib.pyplot.clf()
matplotlib.pyplot.cla()
matplotlib.pyplot.close()



###
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

  pos_higher = {}
  for k, v in pos.items():
    if(v[1]>0):
      pos_higher[k] = (v[0]-0.015, v[1]+0.015)
    else:
      pos_higher[k] = (v[0]-0.015, v[1]-0.015)
  
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
  # matplotlib.pyplot.tight_layout()
  matplotlib.pyplot.savefig(os.path.join(
      os.path.expanduser(output_directory), 
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
  # matplotlib.pyplot.tight_layout()
  matplotlib.pyplot.savefig(os.path.join(
      os.path.expanduser(output_directory), 
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
  sub_ax = matplotlib.pyplot.axes([0.90, 0.55, 0.02, 0.3]) 
  cbar = matplotlib.pyplot.colorbar(sm, cax=sub_ax)
  #matplotlib.pyplot.tight_layout()
  matplotlib.pyplot.savefig(os.path.join(
      os.path.expanduser(output_directory), 
      (label + "_depscore.pdf")),
      dpi=300,
      edgecolor='w',
      bbox_inches='tight')  
  matplotlib.pyplot.show()
      
##
plot_cluster_net_depscore(cluster=C1, label = "cluster 1", pos=posC, edge_color = "lightgreen", node_color = "green")  
plot_cluster_net_depscore(cluster=C2, label = "cluster 2", pos=posC, edge_color = "gold", node_color = "orange")  
plot_cluster_net_depscore(cluster=C3, label = "cluster 3", pos=posC, edge_color = "skyblue", node_color = "blue")  
plot_cluster_net_depscore(cluster=C4, label = "cluster 4", pos=posC, edge_color = "pink", node_color = "red")  
plot_cluster_net_depscore(cluster=G, pos = pos)  



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
  matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), (label + "_Depscore histogram.pdf")),
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
  matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), 
  (label + "_Degree histogram.pdf")),
                      dpi=300,
                      edgecolor='w',
                      #transparent=True,
                      bbox_inches='tight') 
  matplotlib.pyplot.show()

#plot_depscore_histogram(cluster = G, pos= pos, label = "all proteins", edge_color = "grey", node_color = "black")    
plot_depscore_histogram(cluster = C1, pos= posC, label = "cluster 1", edge_color = "lightgreen", node_color = "green")    
plot_depscore_histogram(cluster = C2, pos= posC, label = "cluster 2", edge_color = "gold", node_color = "orange")  
plot_depscore_histogram(cluster = C3, pos= posC, label = "cluster 3", edge_color = "skyblue", node_color = "blue")  
plot_depscore_histogram(cluster = C4, pos= posC, label = "cluster 4", edge_color = "pink", node_color = "red")  



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
  matplotlib.pyplot.clf()
  fig, ax = matplotlib.pyplot.subplots(figsize=(10, 10))
  g = seaborn.JointGrid(data=cluster_degree_depscore, x=('depscore'), y=('degree'))
  g.plot_joint(seaborn.regplot, fit_reg = True, color="royalblue")
  g.plot_marginals(seaborn.histplot, kde=True, bins=80, color="royalblue")
  ax = g.ax_joint
  ax.set_yscale('log')
  ax.set_xscale('log')
  r, p = scipy.stats.pearsonr(cluster_degree_depscore['depscore'], cluster_degree_depscore['degree'])
  g.ax_joint.annotate(f'$\\rho = {r:.3f}, p = {p:.3f}$',
                      xy=(0.1, 0.9), xycoords='axes fraction',
                      ha='left', va='center')
  g.ax_joint.axvline(x=-1)
  matplotlib.pyplot.show()




ccle_dep_score = pandas.read_csv(os.path.join(data_file_directory, "D2_combined_gene_dep_scores.csv")).T
## Convert Row to Column Header
new_header = ccle_dep_score.iloc[0] #grab the first row for the header
ccle_dep_score = ccle_dep_score[1:] #take the data less the header row
ccle_dep_score.columns = new_header #set the header row as the df header

EIF_cols = [col for col in ccle_dep_score.columns if "EIF4E" in col]
print(EIF_cols)


ccle_anno = pandas.read_csv(os.path.join(data_file_directory, "sample_info.csv"))
ccle_cnv = pandas.read_csv(os.path.join(data_file_directory, "CCLE_gene_cn.csv"))



def plot_dep_cnv(protein):
  result = pandas.merge(ccle_anno[['DepMap_ID', 'CCLE_Name']], 
                      ccle_cnv[['Unnamed: 0',protein]], 
                      left_on='DepMap_ID', 
                      right_on='Unnamed: 0',
                      suffixes=('_dep', '_cnv'))
  dep_cnv = pandas.merge(ccle_dep_score[[protein]], 
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

ccle_dep_score.dtypes


g = sns.jointplot(data=dep_cnv, x='EIF4G1 (1981)_dep', y='EIF4G1 (1981)_cnv', kind='reg', color='royalblue')
g.ax_joint.axvline(x=1)

# ax.annotate(stats.pearsonr)
g.ax_joint.annotate(f'$\\rho = {r:.3f}, p = {p:.3f}$',
                    xy=(0.1, 0.9), xycoords='axes fraction',
                    ha='left', va='center')
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.show()



