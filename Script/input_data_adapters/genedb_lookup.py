"""This is a simple interface class for gene database lookups."""


class GeneDBLookup:
  def translate_gene_symbol_to_entrez_id(self, gene_symbol: str) -> str:
    pass  # only implemented by subclasses
