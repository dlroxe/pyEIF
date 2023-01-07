"""Tests for entrez_lookup.py."""

from absl import logging
from absl.testing import absltest

import os
import entrez_lookup


class UnitTests(absltest.TestCase):

  def test_degenerate_case(self):
    handle = entrez_lookup.EntrezLookup(hs_file=None)
    self.assertIsNone(handle.translate_gene_symbol_to_entrez_id('foo'))

  def test_first20(self):
    handle = entrez_lookup.EntrezLookup(
      hs_file=os.path.join('test_data', 'entrez_data', 'first_20_Hs.data'))
    self.assertEqual('10', handle.translate_gene_symbol_to_entrez_id('NAT2'))

  def test_first3_short(self):
    handle = entrez_lookup.EntrezLookup(
      hs_file=os.path.join('test_data', 'entrez_data', 'first_3_short_Hs.data'))
    self.assertEqual('10', handle.translate_gene_symbol_to_entrez_id('NAT2'))

  # TODO(dlroxe): This could be used to evaluate *every* Entrez ID.
  #               The short_Hs.data file has all the GENE_ID entries that
  #               the original Hs.data file has.
  def test_all_short(self):
    handle = entrez_lookup.EntrezLookup(
      hs_file=os.path.join('test_data', 'entrez_data', 'short_Hs.data'))
    self.assertEqual('10', handle.translate_gene_symbol_to_entrez_id('NAT2'))
