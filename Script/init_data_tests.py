"""Tests for init_data.py."""

import pandas
from absl.testing import absltest

import init_data


class UnitTests(absltest.TestCase):
  def test_treshold_label_substitution(self):
    raw_threshold_data = pandas.DataFrame({
      'Sample': ['gene1', 'gene2', 'gene3', 'gene4', 'gene5'],
      'TCGA-A5-A0GI-01': [2.0, 0.0, 0.0, 0.0, 1.0],
      'TCGA-S9-A7J2-01': [-2.0, -1.0, -1.0, -1.0, -1.0],
      'TCGA-06-0150-01': [0.0, 0.0, 0.0, 0.0, 0.0],
    }).set_index('Sample')

    expected_threshold_data = pandas.DataFrame({
      'Sample': ['gene1', 'gene2', 'gene3', 'gene4', 'gene5'],
      'TCGA-A5-A0GI-01': ['AMP', 'DIPLOID', 'DIPLOID', 'DIPLOID', 'DUP'],
      'TCGA-S9-A7J2-01': ['HOMDEL', 'DEL', 'DEL', 'DEL', 'DEL'],
      'TCGA-06-0150-01': ['DIPLOID', 'DIPLOID', 'DIPLOID', 'DIPLOID', 'DIPLOID'],
    }).set_index('Sample')

    all_threshold_data = init_data.TcgaCnvParser.get_tcga_cnv(values_data_frame=raw_threshold_data)
    self.assertTrue(all_threshold_data.equals(expected_threshold_data),
                    all_threshold_data.compare(expected_threshold_data))  # the 'compare' diff is printed for errors
