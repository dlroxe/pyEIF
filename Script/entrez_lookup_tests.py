"""Tests for entrez_lookup.py."""

from absl.testing import absltest

import entrez_lookup


class UnitTests(absltest.TestCase):
  def test_degenerate_case(self):
    entrez_handle = entrez_lookup.EntrezLookup(hs_file=None)
    self.assertTrue(
      entrez_handle.translate_gene_symbol_to_entrez_id('foo') is None)