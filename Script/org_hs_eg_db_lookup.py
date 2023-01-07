"""
This file provides abstractions for org.Hs.eg.db lookups.

The BioConductor project maintains an 'org.Hs.eg.db' package for the R language,
which is versioned and periodically updated with a fresh snapshot of
publicly-available data.  At the core of that package is an SqlLite database.

This module assumes that the source code for that package has been downloaded
and installed, and provides lookup functions against the contents of the
SqlLite database it contains.
"""
from absl import logging
from typing import List

import sqlite3
import textwrap


class OrgHsEgDbLookup:
  """
  This is a stateful container that opens the org.Hs.eg.db database and
  organizes its contents (e.g. by joining 'alias' and 'genes' tables in a new
  in-memory table), then provides access to the data via specific interface
  functions.  Put another way: This is an abstraction barrier that 'hides' SQL
  statements from the rest of the package.
  """

  def __init__(self, org_hs_eg_db_file: str):
    self._memdb = sqlite3.connect(':memory:')
    self._init_memdb(org_hs_eg_db_file)

  def translate_gene_symbol_to_entrez_id(self, gene_symbol: str) -> str:
    """Returns the Entrez ID for 'gene_symbol', or None if lookup fails."""
    rows = self._execute(textwrap.dedent(f'''\
      select
        gene_id
      from
       alias2entrez
      where
       alias_symbol = '{gene_symbol}'
    '''), verbose=False)
    if len(rows) != 1:
      logging.warning(
        f'found {len(rows)} Entrez ID rows for {gene_symbol}: {rows}')
    return rows[0][0] if rows else None

  def _init_memdb(self, org_hs_eg_db_file: str) -> None:
    logging.info('accessing org_hs_eg_db: %s', org_hs_eg_db_file)
    # force read-only mode; prevent any possible modification of data on disk
    with sqlite3.connect(
        f'file:{org_hs_eg_db_file}?mode=ro', uri=True) as diskdb:
      logging.info('copying org_hs_eg_db to in-memory db.')
      self._memdb.executescript("".join(line for line in diskdb.iterdump()))
      self._memdb.commit()

      logging.info('creating alias-entrez lookup table')
      self._execute(cmd=textwrap.dedent('''
        create table alias2entrez as
          select
            alias_symbol,
            gene_id
          from
            alias cross join genes on alias._id = genes._id
        ;
        '''), verbose=True)
      self._memdb.commit()

  def _execute(self, cmd: str, verbose: bool) -> List[sqlite3.Row]:
    if verbose:
      logging.info('executing: %s', cmd)
    cur = self._memdb.execute(cmd)
    return cur.fetchall()
