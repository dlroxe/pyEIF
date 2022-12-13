import datatable
import os
import lifelines
import matplotlib
import pandas
import subprocess

# critical parameter setting!
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams['figure.figsize'] = 15, 15
pandas.options.mode.chained_assignment = None  # default='warn'

# setup current direcotry
data_file_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_data"
output_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_output"


##### Call R.script to perform heatmap clustering and pathway analysis #####
subprocess.call("Rscript /home/suwu/github/pyEIF/Script/Fig1.R", shell=True)

## KM survival analysis
TCGA_OS = pandas.read_table(os.path.join(data_file_directory, 
"Survival_SupplementalTable_S1_20171025_xena_sp")).set_index('sample')
TCGA_OS = TCGA_OS[["OS", "OS.time"]]

TCGA_CNV = datatable.fread(os.path.join(os.path.expanduser(output_directory),"ProcessedData","TCGA_CNV.csv"), header=True).to_pandas()
TCGA_CNV = TCGA_CNV.set_index(['C10'])
TCGA_CNV_EIF = TCGA_CNV[["EIF4G1", "EIF3E", "EIF3H"]]


TCGA_CNV_OS_EIF = pandas.merge(TCGA_OS, 
                               TCGA_CNV_EIF, 
                               left_index=True, 
                               right_index=True).dropna()

T = TCGA_CNV_OS_EIF["OS.time"]
E = TCGA_CNV_OS_EIF["OS"]

#This function creates a Kaplan Meier plot using the matplotlib and lifelines libraries in Python. 
# The Kaplan Meier method is used to estimate the survival function from censored data. 
# In this case, the data is censored by overall survival time after a diagnosis of cancer. 
# The function takes two gene names as input and creates a plot with four lines showing the estimated survival probabilities for the following groups:
    #Patients with amp and dup in gene 1
    #Patients with amp and dup  in gene 2
    #Patients with amp and dup  in both gene 1 and gene 2
    #Patients with diploid copies of both genes
# The function also performs a log-rank test on each of these groups and compares them to the diploid group to see if there is a significant difference in survival rates between the groups. 
# The results of these tests are displayed on the plot.
def combined_survival_plot(gene1, gene2):
  # Define the CNV status for each gene
  gene1_GAIN = (TCGA_CNV_OS_EIF[gene1] == "AMP")|(TCGA_CNV_OS_EIF[gene1] == "DUP")
  gene2_GAIN = (TCGA_CNV_OS_EIF[gene2] == "AMP")|(TCGA_CNV_OS_EIF[gene2] == "DUP")
  gene1_gene2_GAIN = gene1_GAIN & gene2_GAIN
  DIPLOID = (TCGA_CNV_OS_EIF[gene1] == "DIPLOID")&(TCGA_CNV_OS_EIF[gene2] == "DIPLOID")
 
  # Create a Kaplan Meier plot
  matplotlib.pyplot.clf()
  ax = matplotlib.pyplot.subplot()
  kmf = lifelines.KaplanMeierFitter()
  kmf.fit(T[gene1_GAIN], event_observed=E[gene1_GAIN], label=(gene1+"GAIN"))
  kmf.plot_survival_function(ax=ax)
  kmf.fit(T[gene2_GAIN], event_observed=E[gene2_GAIN], label=(gene2+"GAIN"))
  kmf.plot_survival_function(ax=ax)
  kmf.fit(T[gene1_gene2_GAIN], event_observed=E[gene1_gene2_GAIN], label=(gene1+gene2+"GAIN"))
  kmf.plot_survival_function(ax=ax)
  kmf.fit(T[DIPLOID], event_observed=E[DIPLOID], label="DIPLOID")
  kmf.plot_survival_function(ax=ax)
  ax.set(
      title='Kaplan Meier estimates by copy number status',
      xlabel='Overall survival (days after diagnosis)',
      ylabel='Estimated probability of survival'
  )  
  result1 = lifelines.statistics.logrank_test(T[gene1_GAIN], T[DIPLOID], E[gene1_GAIN], E[DIPLOID], alpha=.99)
  result2 = lifelines.statistics.logrank_test(T[gene2_GAIN], T[DIPLOID], E[gene2_GAIN], E[DIPLOID], alpha=.99)
  result3 = lifelines.statistics.logrank_test(T[gene1_gene2_GAIN], T[DIPLOID], E[gene1_gene2_GAIN], E[DIPLOID], alpha=.99)
  result4 = lifelines.statistics.logrank_test(T[gene1_GAIN], T[gene1_gene2_GAIN], E[gene1_GAIN], E[gene1_gene2_GAIN], alpha=.99)
  result5 = lifelines.statistics.logrank_test(T[gene2_GAIN], T[gene1_gene2_GAIN], E[gene2_GAIN], E[gene1_gene2_GAIN], alpha=.99)  
  ax.text(1, 0.85, (gene1 + ' GAIN vs DIPLOID '+ "p = %.3f" % result1.p_value))
  ax.text(1, 0.8, (gene2 + ' GAIN vs DIPLOID '+ "p = %.3f" % result2.p_value))
  ax.text(1, 0.75, (gene1 + gene2 + ' GAIN vs DIPLOID '+ "p = %.3f" %result3.p_value))
  ax.text(1, 0.7, (gene1 + ' GAIN vs '+ gene1 + gene2 + 'GAIN '+ "p = %.3f" %result4.p_value))
  ax.text(1, 0.65, (gene2 + ' GAIN vs '+ gene1 + gene2 + 'GAIN '+ "p = %.3f" %result5.p_value))  
  matplotlib.pyplot.show()
  matplotlib.pyplot.savefig(os.path.join(
      os.path.expanduser(output_directory), "Fig1",
      ("Gain_KM.pdf")),
      dpi=300,
      edgecolor='w',
      bbox_inches='tight') 
  
combined_survival_plot(gene1="EIF3E", gene2="EIF4G1")


#This function creates a Kaplan Meier plot using the matplotlib and lifelines libraries in Python. 
# The Kaplan Meier method is used to estimate the survival function from censored data. 
# In this case, the data is censored by overall survival time after a diagnosis of cancer. 
# The function takes two gene names as input and creates a plot with four lines showing the estimated survival probabilities for the following groups:
    #Patients with amp in gene 1
    #Patients with amp in gene 2
    #Patients with amp in both gene 1 and gene 2
    #Patients with diploid copies of both genes
# The function also performs a log-rank test on each of these groups and compares them to the diploid group to see if there is a significant difference in survival rates between the groups. 
# The results of these tests are displayed on the plot.
def combined_survival_plot2(gene1, gene2):
  gene1_GAIN = (TCGA_CNV_OS_EIF[gene1] == "AMP")
  gene2_GAIN = (TCGA_CNV_OS_EIF[gene2] == "AMP")
  gene1_gene2_GAIN = gene1_GAIN & gene2_GAIN
  DIPLOID = (TCGA_CNV_OS_EIF[gene1] == "DIPLOID")&(TCGA_CNV_OS_EIF[gene2] == "DIPLOID")

  matplotlib.pyplot.clf()
  ax = matplotlib.pyplot.subplot()
  kmf = lifelines.KaplanMeierFitter()
  kmf.fit(T[gene1_GAIN], event_observed=E[gene1_GAIN], label=(gene1+"GAIN"))
  kmf.plot_survival_function(ax=ax)
  kmf.fit(T[gene2_GAIN], event_observed=E[gene2_GAIN], label=(gene2+"GAIN"))
  kmf.plot_survival_function(ax=ax)
  kmf.fit(T[gene1_gene2_GAIN], event_observed=E[gene1_gene2_GAIN], label=(gene1+gene2+"GAIN"))
  kmf.plot_survival_function(ax=ax)
  kmf.fit(T[DIPLOID], event_observed=E[DIPLOID], label="DIPLOID")
  kmf.plot_survival_function(ax=ax)
  ax.set(
      title='Kaplan Meier estimates by copy number status',
      xlabel='Overall survival (days after diagnosis)',
      ylabel='Estimated probability of survival'
  )  
  result1 = lifelines.statistics.logrank_test(T[gene1_GAIN], T[DIPLOID], E[gene1_GAIN], E[DIPLOID], alpha=.99)
  result2 = lifelines.statistics.logrank_test(T[gene2_GAIN], T[DIPLOID], E[gene2_GAIN], E[DIPLOID], alpha=.99)
  result3 = lifelines.statistics.logrank_test(T[gene1_gene2_GAIN], T[DIPLOID], E[gene1_gene2_GAIN], E[DIPLOID], alpha=.99)
  result4 = lifelines.statistics.logrank_test(T[gene1_GAIN], T[gene1_gene2_GAIN], E[gene1_GAIN], E[gene1_gene2_GAIN], alpha=.99)
  result5 = lifelines.statistics.logrank_test(T[gene2_GAIN], T[gene1_gene2_GAIN], E[gene2_GAIN], E[gene1_gene2_GAIN], alpha=.99)  
  ax.text(1, 0.85, (gene1 + ' AMP vs DIPLOID '+ "p = %.3f" % result1.p_value))
  ax.text(1, 0.8, (gene2 + ' AMP vs DIPLOID '+ "p = %.3f" % result2.p_value))
  ax.text(1, 0.75, (gene1 + gene2 + ' AMP vs DIPLOID '+ "p = %.3f" %result3.p_value))
  ax.text(1, 0.7, (gene1 + ' AMP vs '+ gene1 + gene2 + 'AMP '+ "p = %.3f" %result4.p_value))
  ax.text(1, 0.65, (gene2 + ' AMP vs '+ gene1 + gene2 + 'AMP '+ "p = %.3f" %result5.p_value))  
  matplotlib.pyplot.show()
  matplotlib.pyplot.savefig(os.path.join(
      os.path.expanduser(output_directory), "Fig1",
      ("AMP_KM.pdf")),
      dpi=300,
      edgecolor='w',
      bbox_inches='tight') 
  
combined_survival_plot2(gene1="EIF3E", gene2="EIF4G1")

