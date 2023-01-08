"""Tests for init_data.py."""

from absl.testing import absltest

import entrez_lookup
import init_data
import pandas


class UnitTests(absltest.TestCase):
  def test_threshold_label_substitution(self):
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
      'TCGA-06-0150-01':
        ['DIPLOID', 'DIPLOID', 'DIPLOID', 'DIPLOID', 'DIPLOID'],
    }).set_index('Sample')

    all_threshold_data = init_data.TcgaCnvParser.get_tcga_cnv(
      values_data_frame=raw_threshold_data)
    self.assertTrue(all_threshold_data.equals(expected_threshold_data),
                    msg=all_threshold_data.compare(expected_threshold_data))

  def test_phenotype_join(self):
    threshold_data = pandas.DataFrame({
      'Sample': ['gene1', 'gene2', 'gene3', 'gene4', 'gene5'],
      'TCGA-A5-A0GI-01': ['AMP', 'DIPLOID', 'DIPLOID', 'DIPLOID', 'DUP'],
      'TCGA-S9-A7J2-01': ['HOMDEL', 'DEL', 'DEL', 'DEL', 'DEL'],
      'TCGA-06-0150-01':
        ['DIPLOID', 'DIPLOID', 'DIPLOID', 'DIPLOID', 'DIPLOID'],
    }).set_index('Sample').transpose()

    phenotype_data = pandas.DataFrame({
      'Sample': ['TCGA-A5-A0GI-01', 'TCGA-S9-A7J2-01', 'TCGA-06-0150-01'],
      'sample.type': ['type1', 'type1', 'type2'],
      'primary_disease': ['disease1', 'disease2', 'disease2'],
    }).set_index('Sample')

    expected_joined_data = pandas.DataFrame({
      'Sample':
        ['gene1', 'gene2', 'gene3', 'gene4', 'gene5', 'sample.type',
         'primary_disease'],
      'TCGA-A5-A0GI-01':
        ['AMP', 'DIPLOID', 'DIPLOID', 'DIPLOID', 'DUP', 'type1', 'disease1'],
      'TCGA-S9-A7J2-01':
        ['HOMDEL', 'DEL', 'DEL', 'DEL', 'DEL', 'type1', 'disease2'],
      'TCGA-06-0150-01':
        ['DIPLOID', 'DIPLOID', 'DIPLOID', 'DIPLOID', 'DIPLOID', 'type2',
         'disease2'],

    }).set_index('Sample').transpose()

    joined_data = init_data.TcgaCnvParser.merge_cnv_phenotypes(
      cnv_data=threshold_data, phenotype_data=phenotype_data)

    self.assertTrue(joined_data.equals(expected_joined_data),
                    msg='\n' + joined_data.compare(expected_joined_data))

  def test_top_genes(self):
    threshold_data = pandas.DataFrame({
      'Sample': ['gene1', 'gene2', 'gene3', 'gene4', 'gene5'],
      'TCGA-A5-A0GI-01': ['AMP', 'DIPLOID', 'DIPLOID', 'DIPLOID', 'DUP'],
      'TCGA-S9-A7J2-01': ['HOMDEL', 'DEL', 'DEL', 'DEL', 'DEL'],
      'TCGA-06-0150-01': ['DIPLOID', 'DIPLOID', 'DIPLOID', 'DIPLOID',
                          'DIPLOID'],
    }).set_index('Sample').transpose()
    entrez_handle = entrez_lookup.EntrezLookup(hs_file=None)
    top_genes = init_data.TcgaCnvParser.get_top_genes  # for convenience
    top_genes01 = top_genes(threshold_data, labels=['AMP'], percent=30,
                            genedb_handle=entrez_handle)
    top_genes02 = top_genes(threshold_data, labels=['AMP'], percent=40,
                            genedb_handle=entrez_handle)
    top_genes03 = top_genes(threshold_data, labels=['AMP', 'DUP'], percent=30,
                            genedb_handle=entrez_handle)
    top_genes04 = top_genes(threshold_data, labels=['AMP', 'DUP'], percent=40,
                            genedb_handle=entrez_handle)
    top_genes05 = top_genes(threshold_data, labels=['AMP', 'DIPLOID', 'DUP'],
                            percent=30, genedb_handle=entrez_handle)
    top_genes06 = top_genes(threshold_data, labels=['AMP', 'DIPLOID', 'DUP'],
                            percent=60, genedb_handle=entrez_handle)
    top_genes07 = top_genes(threshold_data, labels=['AMP', 'DIPLOID', 'DUP'],
                            percent=70, genedb_handle=entrez_handle)
    top_genes08 = top_genes(threshold_data, labels=['DIPLOID'], percent=30,
                            genedb_handle=entrez_handle)
    top_genes09 = top_genes(threshold_data, labels=['DIPLOID'], percent=60,
                            genedb_handle=entrez_handle)

    # TODO(dlroxe): These are simple assertions about counts.
    #               Add specific tests of gene names and percentages.
    self.assertEqual(1, len(top_genes01), msg='\n' + str(top_genes01))
    self.assertEqual(0, len(top_genes02), msg='\n' + str(top_genes02))
    self.assertEqual(2, len(top_genes03), msg='\n' + str(top_genes03))
    self.assertEqual(0, len(top_genes04), msg='\n' + str(top_genes04))
    self.assertEqual(5, len(top_genes05), msg='\n' + str(top_genes05))
    self.assertEqual(5, len(top_genes06), msg='\n' + str(top_genes06))
    self.assertEqual(0, len(top_genes07), msg='\n' + str(top_genes07))
    self.assertEqual(5, len(top_genes08), msg='\n' + str(top_genes08))
    self.assertEqual(3, len(top_genes09), msg='\n' + str(top_genes06))


if __name__ == '__main__':
  absltest.main()
