"""Tests for init_data.py."""

from absl.testing import absltest

import org_hs_eg_db_lookup
import tcga_cnv_parser
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
    }).set_index('Sample').astype('category')

    parser = tcga_cnv_parser.TcgaCnvParser('', '', '', '', '')
    all_threshold_data = parser.get_tcga_cnv(
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

    parser = tcga_cnv_parser.TcgaCnvParser('', '', '', '', '')
    joined_data = parser.merge_cnv_phenotypes(
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

    sample_count = len(threshold_data.index)
    melted_threshold_data = tcga_cnv_parser.TcgaCnvParser.melt_threshold_data(
      threshold_data)

    genedb_handle = org_hs_eg_db_lookup.OrgHsEgDbLookup(org_hs_eg_db_file=None)
    top_genes = tcga_cnv_parser.TcgaCnvParser.get_top_genes  # for convenience

    top_genes01 = top_genes(sample_count=sample_count, df=melted_threshold_data,
                            labels=['AMP'], percent=30,
                            genedb_handle=genedb_handle)
    top_genes02 = top_genes(sample_count=sample_count, df=melted_threshold_data,
                            labels=['AMP'], percent=40,
                            genedb_handle=genedb_handle)
    top_genes03 = top_genes(sample_count=sample_count, df=melted_threshold_data,
                            labels=['AMP', 'DUP'],
                            percent=30,
                            genedb_handle=genedb_handle)
    top_genes04 = top_genes(sample_count=sample_count, df=melted_threshold_data,
                            labels=['AMP', 'DUP'],
                            percent=40,
                            genedb_handle=genedb_handle)
    top_genes05 = top_genes(sample_count=sample_count, df=melted_threshold_data,
                            labels=['AMP', 'DIPLOID', 'DUP'],
                            percent=30, genedb_handle=genedb_handle)
    top_genes06 = top_genes(sample_count=sample_count, df=melted_threshold_data,
                            labels=['AMP', 'DIPLOID', 'DUP'],
                            percent=60, genedb_handle=genedb_handle)
    top_genes07 = top_genes(sample_count=sample_count, df=melted_threshold_data,
                            labels=['AMP', 'DIPLOID', 'DUP'],
                            percent=70, genedb_handle=genedb_handle)
    top_genes08 = top_genes(sample_count=sample_count, df=melted_threshold_data,
                            labels=['DIPLOID'],
                            percent=30,
                            genedb_handle=genedb_handle)
    top_genes09 = top_genes(sample_count=sample_count, df=melted_threshold_data,
                            labels=['DIPLOID'],
                            percent=60,
                            genedb_handle=genedb_handle)

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
