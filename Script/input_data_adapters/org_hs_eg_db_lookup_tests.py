"""Tests for hs_data_lookup.py."""

from absl.testing import absltest

import org_hs_eg_db_lookup
import pandas


class UnitTests(absltest.TestCase):

  def test_degenerate_case(self):
    handle = org_hs_eg_db_lookup.OrgHsEgDbLookup(org_hs_eg_db_file=None)
    self.assertEqual(pandas.NA,
                     handle.translate_gene_symbol_to_entrez_id('foo'))
