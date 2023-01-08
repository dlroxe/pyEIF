"""
This file provides abstractions for reactome.db lookups.

The BioConductor project maintains a 'reactome.db' package for the R language,
which is versioned and periodically updated with a fresh snapshot of
publicly-available data.  At the core of that package is an SWLite database.

This module assumes that the source code for that package has been downloaded
and installed, and provides lookup functions against the contents of the
SQLite database it contains.
"""
from absl import logging
from typing import List

import sqlite3


# TODO(dlroxe): It might make sense to abandon this nascent wrapper before
#               completing its development.  The idea here is to permit
#               pathway enrichment using exactly the same data source that
#               corresponding R code uses (reactome.db), which offers
#               bug-for-bug compatibility.
#
#               A possible mapping at the SQLIte level might be:
#               GenomeEncodedEntity
#               [==>CatalystActivity]
#               ==> ReactionLikeEvent
#               ==> Pathway
#
#               However it's not so clear e.g. how the "p value cutoff"
#               from ReactomePA::enrichPathway() could be brought to bear.
#               Thus, this appraoch may not be easily feasible, so instead it
#               may make sense to use the "latest" data here (albeit with the
#               same problem for p value):
#
#               https://reactome.org/download/current/NCBI2Reactome.txt
#
#               Another way might be to use "rpy2" to invoke the underlying
#               R package just enough to make it print a CSV file with the
#               needed data frame.  That is probably more trouble than it's
#               worth; for example a shell script could invoke the R file to
#               do something similar.
#
#               Yet another approach would be to use the reactome.org
#               REST API to download all required data in JSON format,
#               and parse it e.g. with Protocol Buffers.  (This approach
#               probably has the least merit.)

class ReactomeLookup:
  """
  This is a stateful container that opens the reactome.db database and
  organizes its contents (e.g. by creating certain lookup tables),
  then provides access to the data via specific interface
  functions.  Put another way: This is an abstraction barrier that 'hides' SQL
  statements from the rest of the package.
  """

  def __init__(self, reactome_db_file: str):
    self._memdb = sqlite3.connect(':memory:')
    self._init_memdb(reactome_db_file)
    # TODO(dlroxe): remove test code below this line
    # res = self._execute('PRAGMA table_info(sqlite_schema);', verbose=True)
    # logging.info('schema:\n%s', '\n'.join([str(row) for row in res]))
    res = self._execute(
      'select name from sqlite_schema where type = "table" and name not like "sqlite_%";',
      verbose=True)
    logging.info('names:\n%s', '\n'.join([str(row) for row in res]))

  def _init_memdb(self, reactome_db_file: str) -> None:
    logging.info('accessing reactome_db: %s', reactome_db_file)
    # force read-only mode; prevent any possible modification of data on disk
    with sqlite3.connect(
        f'file:{reactome_db_file}?mode=ro', uri=True) as diskdb:
      # TODO(dlroxe): If this approach is retained, then find some way to
      #               trim down the data (e.g. just to Homo Sapiens)
      #               *before* copying it all.  The in-memory representation
      #               of the whole data base is ~7G in size, and clunky.
      logging.info('copying reactome_db to in-memory db.')
      self._memdb.executescript("".join(line for line in diskdb.iterdump()))
      self._memdb.commit()

      # ... other initialization here ...
      self._memdb.commit()

  def _execute(self, cmd: str, verbose: bool) -> List[sqlite3.Row]:
    if verbose:
      logging.info('executing: %s', cmd)
    cur = self._memdb.execute(cmd)
    return cur.fetchall()
