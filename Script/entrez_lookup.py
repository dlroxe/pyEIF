"""
This file provides abstractions for "Entrez" lookups.

For example, one way of mapping gene symbols to Entrez IDs is to read a list of
entries from a local flat file such as 'Hs.data'.  Another way is to make a
symbol-specific call to a remote Entrez database.

The purpose of this module is to "hide" from the rest of the package the
particular implementation for Entrez lookups, so that it may be updated or
changed without impacting any other code.
"""
from Bio import UniGene
from typing import Dict, Optional
from absl import logging


class EntrezLookup:
  """
  This is a stateful container that implements all methods necessary for
  "Entrez" interactions.  For example, when instantiated, it reads a local
  'Hs.data' file and retains in memory a lookup dictionary for the data therein.
  """

  # TODO(dlroxe): Hmm, looks like 'Hs.data' started going stale around 2013.  Is
  #               there another source (that doesn't require an active network
  #               connection, preferably, so that data can be versioned and
  #               archived)?
  #               Hmm, perhaps 'AnnoKey' or something like it could cache?
  #               http://bjpop.github.io/annokey/
  #               In particular, perhaps this could be helpful:
  #               get_ncbi_gene_snapshot_xml.py

  def __init__(self, hs_file: Optional[str]):
    # It may be useful in tests to short-circuit lookups to immediate 'None'.
    self._hs_data_dict = self._init_hs_data(hs_file) if hs_file else None

  def translate_gene_symbol_to_entrez_id(self, gene_symbol):
    """Returns the Entrez ID for 'gene_symbol', or None if lookup fails."""
    if self._hs_data_dict is None:
      return None
    key = gene_symbol.split('|')[0]
    return self._hs_data_dict[
      key].gene_id if key in self._hs_data_dict else None

  @staticmethod
  def _init_hs_data(hs_file: str) -> Dict[str, UniGene.Record]:
    # Note, in the version available on Jan 1 2023, the Hs.data file included a
    # typo that made it unparseable.  In particular, "ID Hs.716839" included a
    # line containing "ACC == D22S272", which should have only one "=".
    # Furthermore, there seems to be a parse error with the last record
    # "no matter what".
    #
    # Of course, it should be assumed that such a large set of data will have at
    # least one small error *somewhere*.
    #
    # This function includes error-handling logic to make it easy to detect such
    # problems, because the built-in next() functionality for the generator
    # returned by UniGene.parse() is not sufficiently robust.  How nice it would
    # be, simply to write "for record in UniGene.parse()", but alas it's just
    # too brittle and we can't count on it.  Oh, well.
    logging.info('initializing HS data')
    hs_data_dict = {}
    hs_records = None
    with open(hs_file) as raw_hs_data:
      hs_records = UniGene.parse(raw_hs_data)
      # Available attributes in each record are documented here:
      # https://biopython.org/docs/1.76/api/Bio.UniGene.html#Bio.UniGene.Record
      try:
        count = 0
        last_id = ''
        while True:
          try:
            hs_record = next(hs_records)
            hs_data_dict[hs_record.symbol] = hs_record
            last_id = hs_record.ID
            count = count + 1
          except Exception as ex:
            logging.warning(
              'error parsing record %s (previous well-parsed ID: %s):\n%s',
              count, last_id, ex)
            break
      except StopIteration:
        logging.info('all HS records consumed')
    logging.info('HS data initialized')
    return hs_data_dict
