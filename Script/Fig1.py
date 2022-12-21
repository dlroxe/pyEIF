import datatable
import lifelines
import matplotlib
import os
import pandas
import seaborn 
import subprocess

# critical parameter setting!
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams['figure.figsize'] = 15, 15
pandas.options.mode.chained_assignment = None  # default='warn'

# setup current direcotry
data_file_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_data"
output_directory = "~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_output"

# Call R.script to generate TCGA_CNV.csv data, perform heatmap clustering and pathway analysis #
subprocess.call("Rscript /home/suwu/github/pyEIF/Script/Fig1.R", shell=True)

TCGA_CNV = datatable.fread(os.path.join(data_file_directory, 
"Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"), header=True).to_pandas()
TCGA_CNV = TCGA_CNV.set_index('Sample')
TCGA_CNV = TCGA_CNV.T
TCGA_CNV = TCGA_CNV.replace([2], 'AMP')
TCGA_CNV = TCGA_CNV.replace([1], 'DUP')
TCGA_CNV = TCGA_CNV.replace([0], 'DIPLOID')
TCGA_CNV = TCGA_CNV.replace([-1], 'DEL')
TCGA_CNV = TCGA_CNV.replace([-2], 'HOMDEL')
TCGA_CNV.head()

count = TCGA_CNV.apply(pandas.value_counts, normalize=True).T

def cnv_combined_freq_plot(df, cnv):
    ax = matplotlib.pyplot.figure(figsize=(10, 10))
    # Iterate through the five cnvs
    for CNV in cnv:
        # Subset to the airline
        subset = df[CNV]
        # Draw the density plot
        ax = seaborn.kdeplot(subset, linewidth=2, fill=True, label = CNV) 
    # Plot formatting
    matplotlib.pyplot.legend(prop={'size': 16}, title = 'CNV')
    matplotlib.pyplot.title('Density plot of CNV status in TCGA samples')
    matplotlib.pyplot.xlabel('CNV frequency in all tumors')
    matplotlib.pyplot.ylabel("Probability (%)")
    matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), "Fig1", ("gain density plot.pdf")),
                              dpi=300,
                              edgecolor='w',
                              #transparent=True,
                              bbox_inches='tight') 
    matplotlib.pyplot.show()

cnv_combined_freq_plot(count, ['AMP', 'DUP'])
cnv_combined_freq_plot(count, ['AMP', 'DUP', 'DIPLOID', 'DEL', 'HOMDEL'])


##
def cnv_freq_plot(df, cnv, cutoff):
    ax = matplotlib.pyplot.figure(figsize=(10, 10))

    # Draw the density curve with it's area shaded
    seaborn.kdeplot(df[cnv], linewidth=2, fill=True)
    # Invoke kdeplot() again to get reference to axes 
    ax = seaborn.kdeplot(df[cnv])
    
    # Below code to shade partial region is from 
    # https://stackoverflow.com/a/49100655
    
    # Get all the lines used to draw the density curve 
    kde_lines = ax.get_lines()[-1]
    kde_x, kde_y = kde_lines.get_data()
    
    # Use Numpy mask to filter the lines for region 
    # reresenting height greater than 60 inches 
    mask = kde_x > cutoff
    filled_x, filled_y = kde_x[mask], kde_y[mask]
    
    # Shade the partial region 
    ax.fill_between(filled_x, y1=filled_y)
    
    # find the number of genes over cutoff threshhold
    number = len(df[df[cnv]>cutoff])
    
    ax.annotate(str(number)+' genes', xy=(0.06, 2), xytext=(0.07, 8),
            arrowprops=dict(facecolor='black', shrink=0.05))
    
    # vertical line at x = 60 for reference
    matplotlib.pyplot.axvline(x=cutoff, linewidth=2, linestyle='--')
    matplotlib.pyplot.title("Density plot of "+str(cnv)+" status in TCGA samples")
    matplotlib.pyplot.xlabel(str(cnv)+" frequency in all tumors")
    matplotlib.pyplot.ylabel("Probability (%)")
    matplotlib.pyplot.savefig(os.path.join(os.path.expanduser(output_directory), 
                                           "Fig1", 
                                           (cnv + "_density plot.pdf")),
                              dpi=300,
                              edgecolor='w',
                              #transparent=True,
                              bbox_inches='tight') 
    matplotlib.pyplot.show()

cnv_freq_plot(df=count, cnv="AMP", cutoff=0.05)




## KM survival analysis
#TCGA_OS = pandas.read_table(os.path.join(data_file_directory, 
#"Survival_SupplementalTable_S1_20171025_xena_sp")).set_index('sample')
TCGA_OS = datatable.fread(os.path.join(data_file_directory, 
"Survival_SupplementalTable_S1_20171025_xena_sp"), header=True).to_pandas()
TCGA_OS = TCGA_OS.set_index('sample')
TCGA_OS = TCGA_OS[["OS", "OS.time"]]


TCGA_CNV_EIF = TCGA_CNV[["EIF4G1", "EIF3E", "EIF3H"]]


#TCGA_CNV = datatable.fread(os.path.join(os.path.expanduser(output_directory),"ProcessedData","TCGA_CNV.csv"), header=True).to_pandas()
#TCGA_CNV = TCGA_CNV.set_index(['C10'])
#TCGA_CNV_EIF = TCGA_CNV[["EIF4G1", "EIF3E", "EIF3H"]]


TCGA_CNV_OS_EIF = pandas.merge(TCGA_OS, 
                               TCGA_CNV[["EIF4G1", "EIF3E", "EIF3H"]], 
                               left_index=True, 
                               right_index=True).dropna()

T = TCGA_CNV_OS_EIF["OS.time"]
E = TCGA_CNV_OS_EIF["OS"]

# This function creates a Kaplan Meier plot using the matplotlib and lifelines libraries in Python. 
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


# This function creates a Kaplan Meier plot using the matplotlib and lifelines libraries in Python. 
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

