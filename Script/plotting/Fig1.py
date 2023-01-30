import datatable
import enum
import lifelines
import matplotlib
import os
import pandas
import seaborn

from absl import app
from absl import flags
from absl import logging

# critical parameter setting!
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams['figure.figsize'] = 15, 15
pandas.options.mode.chained_assignment = 'raise'  # default='warn'

# TODO(dlroxe): Try to define FLAGs in a single place, and reference them from
#               everywhere that needs them.
FLAGS = flags.FLAGS
flags.DEFINE_string('data_directory',
                    # '~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_data',
                    os.path.join('~', 'Desktop', 'pyeif_data'),
                    'parent dir for data files')
flags.DEFINE_string('output_directory',
                    # '~/Documents/Bioinformatics_analysis/eIF4G-analysis/eIF4G_output',
                    os.path.join('~', 'Desktop', 'pyeif_output'),
                    'parent dir for output')


def cnv_combined_freq_plot(df, cnv):
  ax = matplotlib.pyplot.figure(figsize=(10, 10))
  # Iterate through the five cnvs
  for CNV in cnv:
    # Subset to the airline
    subset = df[CNV]
    # Draw the density plot
    ax = seaborn.kdeplot(subset, linewidth=2, fill=True, label=CNV)
    # Plot formatting
  matplotlib.pyplot.legend(prop={'size': 16}, title='CNV')
  matplotlib.pyplot.title('Density plot of CNV status in TCGA samples')
  matplotlib.pyplot.xlabel('CNV frequency in all tumors')
  matplotlib.pyplot.ylabel("Probability density (%)")
  matplotlib.pyplot.savefig(
    _abspath(
      os.path.join(FLAGS.output_directory, "Fig1", "gain-density-plot.pdf")),
    dpi=300,
    edgecolor='w',
    # transparent=True,
    bbox_inches='tight')
  matplotlib.pyplot.show()


##
def cnv_freq_plot(df, cnv, cutoff):
  # Draw the density curve with its area shaded
  seaborn.kdeplot(df[cnv], linewidth=2, fill=True)

  # Invoke kdeplot() again to get reference to axes
  ax = seaborn.kdeplot(df[cnv])

  # TODO(dlroxe): Why doesn't the following single line work, in lieu of the
  #               two lines above?
  #               ax = seaborn.kdeplot(df[cnv], linewidth=2, fill=True)

  # Below code to shade partial region is from
  # https://stackoverflow.com/a/49100655

  # Get all the lines used to draw the density curve
  kde_lines = ax.get_lines()[-1]
  kde_x, kde_y = kde_lines.get_data()

  # Use Numpy mask to filter the lines for region
  # representing height greater than 60 inches
  mask = kde_x > cutoff
  filled_x, filled_y = kde_x[mask], kde_y[mask]

  # Shade the partial region
  ax.fill_between(filled_x, y1=filled_y)

  # find the number of genes over cutoff threshold
  number = len(df[df[cnv] > cutoff])

  ax.annotate(str(number) + ' genes', xy=(0.06, 2), xytext=(0.07, 8),
              arrowprops=dict(facecolor='black', shrink=0.05))

  # vertical line at x = 60 for reference
  matplotlib.pyplot.axvline(x=cutoff, linewidth=2, linestyle='--')
  matplotlib.pyplot.title(
    "Density plot of " + str(cnv) + " status in TCGA samples")
  matplotlib.pyplot.xlabel(str(cnv) + " frequency in all tumors")
  matplotlib.pyplot.ylabel("Probability density (%)")
  matplotlib.pyplot.savefig(
    _abspath(
      os.path.join(FLAGS.output_directory, "Fig1", f'{cnv}-density plot.pdf')),
    dpi=300,
    edgecolor='w',
    # transparent=True,
    bbox_inches='tight')
  matplotlib.pyplot.show()


class SurvivalComparisonMode(enum.Enum):
  UNKNOWN = 0  # by convention, always make the zeroth enum value 'UNKNOWN'
  AMP = 1
  GAIN = 2  # AMP+DUP


# This function creates a Kaplan Meier plot using the matplotlib and lifelines
# libraries in Python.
#
# The Kaplan Meier method is used to estimate the survival function from
# censored data.
#
# In this case, the data is censored by overall survival time after a diagnosis
# of cancer.
#
# The function takes two gene names as input and creates a plot with four lines
# showing the estimated survival probabilities for the following groups:
#
# Patients with amp|gain in gene 1
# Patients with amp|gain in gene 2
# Patients with amp|gain in both gene 1 and gene 2
# Patients with diploid copies of both genes
#
# The function also performs a log-rank test on each of these groups and
# compares them to the diploid group to see if there is a significant difference
# in survival rates between the groups.
#
# The results of these tests are displayed on the plot.

def combined_survival_plot(
    gene1, gene2, tcga_cnv_os_eif, mode: SurvivalComparisonMode):
  if mode not in (SurvivalComparisonMode.AMP, SurvivalComparisonMode.GAIN):
    # TODO(dlroxe): find a more appropriate exception
    raise Exception("unknown comparison mode for survival plot")

  t = tcga_cnv_os_eif['OS.time']
  e = tcga_cnv_os_eif['OS']

  gene1_gain = tcga_cnv_os_eif[gene1] == "AMP"
  gene2_gain = tcga_cnv_os_eif[gene2] == "AMP"
  gene1_gene2_gain = gene1_gain & gene2_gain
  g1, g2 = tcga_cnv_os_eif[gene1], tcga_cnv_os_eif[gene2]
  diploid = (g1 == "DIPLOID") & (g2 == "DIPLOID")

  if mode == SurvivalComparisonMode.GAIN:
    gene1_gain = gene1_gain | (tcga_cnv_os_eif[gene1] == "DUP")
    gene2_gain = gene2_gain | (tcga_cnv_os_eif[gene2] == "DUP")
    gene1_gene2_gain = gene1_gain & gene2_gain

  # Create a Kaplan Meier plot
  matplotlib.pyplot.clf()
  ax = matplotlib.pyplot.subplot()
  kmf = lifelines.KaplanMeierFitter()
  kmf.fit(
    t[gene1_gain], event_observed=e[gene1_gain], label=(gene1 + mode.name))
  kmf.plot_survival_function(ax=ax)

  kmf.fit(
    t[gene2_gain], event_observed=e[gene2_gain], label=(gene2 + mode.name))
  kmf.plot_survival_function(ax=ax)

  kmf.fit(
    t[gene1_gene2_gain], event_observed=e[gene1_gene2_gain],
    label=(gene1 + gene2 + mode.name))
  kmf.plot_survival_function(ax=ax)

  kmf.fit(t[diploid], event_observed=e[diploid], label="DIPLOID")
  kmf.plot_survival_function(ax=ax)

  ax.set(
    title='Kaplan Meier estimates by copy number status',
    xlabel='Overall survival (days after diagnosis)',
    ylabel='Estimated probability of survival'
  )
  result1 = lifelines.statistics.logrank_test(
    t[gene1_gain], t[diploid], e[gene1_gain], e[diploid], alpha=.99)
  result2 = lifelines.statistics.logrank_test(
    t[gene2_gain], t[diploid], e[gene2_gain], e[diploid], alpha=.99)
  result3 = lifelines.statistics.logrank_test(
    t[gene1_gene2_gain], t[diploid], e[gene1_gene2_gain], e[diploid], alpha=.99)
  result4 = lifelines.statistics.logrank_test(
    t[gene1_gain], t[gene1_gene2_gain], e[gene1_gain], e[gene1_gene2_gain],
    alpha=.99)
  result5 = lifelines.statistics.logrank_test(
    t[gene2_gain], t[gene1_gene2_gain], e[gene2_gain], e[gene1_gene2_gain],
    alpha=.99)

  ax.text(1, 0.85,
          f'{gene1} {mode.name} vs DIPLOID p = {result1.p_value:.3f}')
  ax.text(1, 0.8,
          f'{gene2} {mode.name} vs DIPLOID p = {result2.p_value:.3f}')
  ax.text(1, 0.75,
          f'{gene1} {gene2} {mode.name} vs DIPLOID p = {result3.p_value:.3f}')
  ax.text(1, 0.7,
          f'{gene1} {mode.name} vs {gene1} {gene2} {mode.name} '
          f'p = {result4.p_value:.3f}')
  ax.text(1, 0.65,
          f'{gene2} {mode.name} vs {gene1} {gene2} {mode.name} '
          f'p = {result5.p_value:.3f}')

  matplotlib.pyplot.show()
  matplotlib.pyplot.savefig(
    _abspath(
      os.path.join(FLAGS.output_directory, 'Fig1', f'{mode.name}-KM.pdf')),
    dpi=300,
    edgecolor='w',
    bbox_inches='tight')


# TODO(dlroxe): This function is cropping up in several places.  Find a
#               One True Home for it.
def _abspath(path):
  # What an absurd incantation to resolve "~"; but OK. Thanks, StackOverflow.
  return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def main(unused_argv):
  # KM survival analysis
  logging.info('commencing KM survival analysis')

  supplemental_table_path = _abspath(
    os.path.join(
      FLAGS.data_directory, 'Survival_SupplementalTable_S1_20171025_xena_sp'))
  logging.info(f'reading: {supplemental_table_path}')

  tcga_os = datatable.fread(supplemental_table_path, header=True).to_pandas()
  tcga_os.set_index('sample', inplace=True)
  tcga_os = tcga_os[['OS', 'OS.time']]

  tcga_cnv_path = _abspath(
    os.path.join(
      FLAGS.output_directory, 'ProcessedData', 'TCGA-CNV-thresholds.csv'))
  logging.info(f'reading: {tcga_cnv_path}')
  tcga_cnv = datatable.fread(tcga_cnv_path, header=True).to_pandas()

  # TODO(dlroxe): adjust init_data.py, and coordinate here, so that the
  #               rename(index:sample) is no longer needed
  tcga_cnv = tcga_cnv.rename(columns={'index': 'sample'})
  tcga_cnv.set_index('sample', inplace=True)

  # TODO(dlroxe): 4g1, 3e, 3h should be flag-configurable, not hardcoded
  tcga_os_eif = pandas.merge(tcga_os,
                             tcga_cnv[['EIF4G1', 'EIF3E', 'EIF3H']],
                             left_index=True,
                             right_index=True).dropna()

  logging.info('generating plots')
  count = tcga_cnv.apply(pandas.value_counts, normalize=True).T
  cnv_combined_freq_plot(count, ['AMP', 'DUP'])
  cnv_combined_freq_plot(count, ['AMP', 'DUP', 'DIPLOID', 'DEL', 'HOMDEL'])
  cnv_freq_plot(df=count, cnv="AMP", cutoff=0.05)

  # TODO(dlroxe): 3e, 4g1 should be flag-configurable, not hardcoded
  for mode in (SurvivalComparisonMode.AMP, SurvivalComparisonMode.GAIN):
    combined_survival_plot(
      'EIF3E', 'EIF4G1',
      tcga_os_eif,
      mode)

  logging.info('KM survival analysis complete')


if __name__ == "__main__":
  app.run(main)
