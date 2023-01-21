## EXPERIMENTAL
## use the following python code to perform pathway enrichment analysis
from typing import List
import gseapy


# Call this e.g. as follows:
# enr = generate_enriched_data(top_amp_genes['Gene'])
def generate_enriched_data(gene_list: List[str]):
  return gseapy.enrichr(
    gene_list=gene_list['Gene'],
    # or "./tests/data/gene_list.txt",
    gene_sets=['Reactome_2022'],
    organism='human',
    # don't forget to set organism to the one you desired! e.g. Yeast
    outdir=None,  # don't write to disk
  )


target_gene_list = []  # use a 'real' gene list for 'real results
enr = generate_enriched_data(target_gene_list)

enr.results.head(15)

ax = gseapy.dotplot(
  enr.results,
  column="Adjusted P-value",
  # x='Gene_set', # set x axis, so you could do a multi-sample/library comparsion
  size=10,
  top_term=8,
  figsize=(3, 5),
  cutoff=1,
  title="Reactome",
  xticklabels_rot=45,  # rotate xtick labels
  show_ring=True,  # set to False to revmove outer ring
  # marker='o',
)

ax = gseapy.barplot(
  enr.results,
  column="Adjusted P-value",
  group='Gene_set',
  # set group, so you could do a multi-sample/library comparsion
  size=8,
  top_term=10,
  figsize=(3, 5), cutoff=1,
  color=['darkred', 'darkblue']  # set colors for group
)
##
